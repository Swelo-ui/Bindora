"""
Comprehensive Phase 2 (Days 61-120) Scientific & Engineering Test Suite for Bindora Dock.

Tests:
1. Macrocyclic Ring Conformational Sampling (srETKDGv3 / ETKDGv3, MMFF94, RMSD pruning)
2. Covalent Docking Module (Warhead detection, catalytic nucleophile pairing, Bürgi-Dunitz / SN2 trajectories)
3. Deep Learning Rescoring & Multi-Engine Consensus Matrix (Vina, Vinardo, MM-GBSA, CNN, Strain penalty)
4. Automated PDBbind Core Set / CASF Affinity Correlation Engine (Pearson R, Spearman rho, RMSE, MAE)
5. MM-PBSA & Explicit Solvent OpenMM Simulation Script Generation
"""

import math
import tempfile
from pathlib import Path
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from backend.services.macrocycle import MacrocycleConformerEngine
from backend.services.covalent import CovalentDockingService
from backend.services.consensus import ConsensusScoringService
from backend.services.pdbbind_validation import PDBbindValidationEngine
from backend.services.md_export import OpenMMExportService
from backend.services.refinement import ComplexRefinementService
from backend.services.docking import DockingEngine


class TestMacrocycleSampling:
    """Test suite for macrocycle ring detection and conformer sampling."""

    def test_macrocycle_detection(self):
        # 14-membered aliphatic ring
        is_macro, sizes, _ = MacrocycleConformerEngine.is_macrocycle("C1CCCCCCCCCCCCC1")
        assert is_macro is True
        assert 14 in sizes

        # Non-macrocycle (benzene / cyclohexane / 8-membered ring)
        is_macro_benz, _, _ = MacrocycleConformerEngine.is_macrocycle("c1ccccc1")
        assert is_macro_benz is False
        is_macro_oct, _, _ = MacrocycleConformerEngine.is_macrocycle("C1CCCCCCC1")
        assert is_macro_oct is False

    def test_sample_macrocycle_conformers(self):
        # Sample 14-membered macrocycle
        res = MacrocycleConformerEngine.sample_macrocycle_conformers(
            "C1CCCCCCCCCCCCC1",
            num_confs=10,
            energy_window=15.0,
            rmsd_threshold=0.5,
            random_seed=42
        )
        assert res["is_macrocycle"] is True
        assert res["conformers_retained_after_rmsd_pruning"] > 0
        assert res["best_mol"] is not None
        assert res["best_mol"].GetNumConformers() == 1
        assert res["best_mol"].GetConformer().Is3D() is True
        assert len(res["relative_energies_kcal"]) > 0
        assert res["relative_energies_kcal"][0] == 0.0  # Lowest energy relative is 0.0

    def test_prepare_ligand_macrocycle_integration(self):
        # Verify prepare_ligand in DockingEngine auto-engages macrocycle engine
        prep = DockingEngine.prepare_ligand("C1CCCCCCCCCCCCC1")
        log = prep["prep_log"]
        assert log.get("is_macrocycle") is True
        assert 14 in log.get("macrocycle_ring_sizes", [])
        assert "Macrocycle" in log.get("conformer_algorithm", "")


class TestCovalentDocking:
    """Test suite for covalent warhead detection and reactive geometry analysis."""

    def test_warhead_detection(self):
        # 1. Michael acceptor (acrylamide)
        w_acry = CovalentDockingService.detect_warheads("C=CC(=O)Nc1ccccc1")
        assert len(w_acry) >= 1
        assert w_acry[0]["warhead_type"] == "michael_acceptor"

        # 2. Haloacetamide
        w_halo = CovalentDockingService.detect_warheads("ClCC(=O)Nc1ccccc1")
        assert len(w_halo) >= 1
        assert w_halo[0]["warhead_type"] == "haloacetamide"

        # 3. Sulfonyl fluoride
        w_sufex = CovalentDockingService.detect_warheads("c1ccccc1S(=O)(=O)F")
        assert len(w_sufex) >= 1
        assert w_sufex[0]["warhead_type"] == "sulfonyl_fluoride"

        # 4. Nitrile
        w_nitrile = CovalentDockingService.detect_warheads("N#Cc1ccccc1")
        assert len(w_nitrile) >= 1
        assert w_nitrile[0]["warhead_type"] == "nitrile"

        # 5. Non-covalent ligand (e.g. aspirin or paracetamol)
        w_non = CovalentDockingService.detect_warheads("CC(=O)Nc1ccc(O)cc1")
        assert len(w_non) == 0

    def test_evaluate_covalent_geometry_optimal(self):
        # Setup synthetic receptor with CYS 145 SG at (10, 10, 10)
        rec_pdb = "ATOM      1  SG  CYS A 145      10.000  10.000  10.000  1.00  0.00           S\n"
        # Setup docked pose with beta-carbon of acrylamide within 3.3 A
        lig_pdb = (
            "ATOM      1  C1  LIG L   1      10.000  10.000  13.300  1.00  0.00           C\n"
            "ATOM      2  C2  LIG L   1      10.000  11.200  14.000  1.00  0.00           C\n"
            "ATOM      3  C3  LIG L   1      10.000  11.200  15.500  1.00  0.00           C\n"
            "ATOM      4  O1  LIG L   1      10.000  12.200  16.200  1.00  0.00           O\n"
        )
        res = CovalentDockingService.evaluate_covalent_geometry(
            docked_pose_pdb_or_pdbqt=lig_pdb,
            receptor_pdb_or_pdbqt=rec_pdb,
            smiles="C=CC(=O)N"
        )
        assert res["is_covalent_candidate"] is True
        assert res["covalent_feasibility_score"] >= 0.70
        assert res["feasibility_assessment"] == "OPTIMAL_COVALENT_GEOMETRY"
        assert res["covalent_energy_bonus_kcal"] < 0.0

    def test_evaluate_covalent_geometry_distant(self):
        # Setup synthetic receptor with distant nucleophile (> 6 A away)
        rec_pdb = "ATOM      1  SG  CYS A 145      10.000  10.000  10.000  1.00  0.00           S\n"
        lig_pdb = (
            "ATOM      1  C1  LIG L   1      25.000  25.000  25.000  1.00  0.00           C\n"
            "ATOM      2  C2  LIG L   1      25.000  26.200  26.000  1.00  0.00           C\n"
        )
        res = CovalentDockingService.evaluate_covalent_geometry(
            docked_pose_pdb_or_pdbqt=lig_pdb,
            receptor_pdb_or_pdbqt=rec_pdb,
            smiles="C=CC(=O)N"
        )
        assert res["is_covalent_candidate"] is True
        assert res["covalent_feasibility_score"] == 0.0
        assert res["feasibility_assessment"] == "NO_REACTIVE_PAIR"


class TestConsensusScoring:
    """Test suite for multi-engine deep learning and consensus ranking."""

    def test_consensus_scoring_decoy_filtering(self):
        # Pose 1: Moderate Vina, excellent MM-GBSA, high CNN, low strain
        pose1 = {
            "mode": 1,
            "affinity_kcal": -9.0,
            "vinardo_affinity_kcal": -8.8,
            "mmgbsa_delta_g_kcal": -13.5,
            "ligand_strain_kcal": 1.2,
            "gnina": {"available": True, "cnn_score": 0.88}
        }
        # Pose 2: Deceptive false-positive Vina score, terrible MM-GBSA, high strain
        pose2 = {
            "mode": 2,
            "affinity_kcal": -10.5,
            "vinardo_affinity_kcal": -7.5,
            "mmgbsa_delta_g_kcal": -3.0,
            "ligand_strain_kcal": 8.5,  # High-strain decoy
            "gnina": {"available": True, "cnn_score": 0.25}
        }
        res = ConsensusScoringService.compute_pose_consensus([pose1, pose2])

        # Top consensus rank MUST be pose 1 due to high-strain decoy penalty on pose 2
        assert res[0]["mode"] == 1
        assert res[0]["consensus_rank"] == 1
        assert res[0]["consensus_confidence"] == "HIGH_CONFIDENCE"

        # Pose 2 must be flagged as high-strain decoy
        assert res[1]["mode"] == 2
        assert res[1]["consensus_confidence"] == "DECOY_HIGH_STRAIN"
        assert res[1]["scoring_breakdown"]["strain_penalty_kcal"] > 0.0

    def test_consensus_with_covalent_bonus(self):
        pose = {
            "mode": 1,
            "affinity_kcal": -8.0,
            "vinardo_affinity_kcal": -7.8,
            "mmgbsa_delta_g_kcal": -9.0,
            "ligand_strain_kcal": 1.5,
            "gnina": {"available": True, "cnn_score": 0.70}
        }
        res_standard = ConsensusScoringService.compute_pose_consensus([pose])
        score_std = res_standard[0]["consensus_score"]

        res_cov = ConsensusScoringService.compute_pose_consensus([pose], covalent_bonus_kcal=-3.0)
        score_cov = res_cov[0]["consensus_score"]

        # Covalent bonus should make composite score more favorable (more negative)
        assert score_cov < score_std


class TestPDBbindValidation:
    """Test suite for automated PDBbind core set and CASF correlation benchmarking."""

    def test_reference_dataset_availability(self):
        dataset = PDBbindValidationEngine.get_reference_dataset()
        assert len(dataset) >= 10
        pdb_ids = [d["pdb_id"] for d in dataset]
        assert "1HSG" in pdb_ids
        assert "1M17" in pdb_ids
        assert "1T46" in pdb_ids

    def test_calculate_correlation_metrics(self):
        y_true = [-13.6, -11.1, -12.1, -10.5, -9.7]
        y_pred = [-13.2, -10.8, -11.9, -10.2, -9.4]
        metrics = PDBbindValidationEngine.calculate_correlation_metrics(y_true, y_pred)

        assert metrics["pearson_r"] > 0.95
        assert metrics["spearman_rho"] > 0.95
        assert metrics["rmse_kcal"] < 1.0
        assert metrics["mae_kcal"] < 1.0

    def test_evaluate_benchmark_execution(self):
        predictions = [
            {"pdb_id": "1HSG", "predicted_delta_g_kcal": -13.0, "pose_rmsd_angstroms": 0.8},
            {"pdb_id": "1M17", "predicted_delta_g_kcal": -10.9, "pose_rmsd_angstroms": 1.1},
            {"pdb_id": "1T46", "predicted_delta_g_kcal": -12.0, "pose_rmsd_angstroms": 1.3},
            {"pdb_id": "2X00", "predicted_delta_g_kcal": -10.1, "pose_rmsd_angstroms": 0.7}
        ]
        res = PDBbindValidationEngine.evaluate_benchmark(predictions)
        assert res["status"] == "SUCCESS"
        assert res["complexes_evaluated"] == 4
        assert res["scoring_power_metrics"]["pearson_r"] > 0.90
        assert res["docking_power_metrics"]["docking_success_rate_percent"] == 100.0


class TestOpenMMExport:
    """Test suite for standalone OpenMM MD simulation and MM-PBSA export."""

    def test_generate_simulation_package(self):
        rec_pdb = (
            "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C\n"
            "ATOM      2  CB  ALA A   1      11.000  10.000  10.000  1.00  0.00           C\n"
        )
        lig_pdb = (
            "HETATM    3  C1  LIG L   1      15.000  10.000  10.000  1.00  0.00           C\n"
            "HETATM    4  C2  LIG L   1      16.000  10.000  10.000  1.00  0.00           C\n"
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            pkg = OpenMMExportService.generate_simulation_package(
                receptor_pdb=rec_pdb,
                docked_pose_pdb_or_pdbqt=lig_pdb,
                output_dir=tmp_dir,
                ligand_smiles="CC",
                job_name="test_sim"
            )
            assert pkg["status"] == "SUCCESS"
            dir_path = Path(tmp_dir)
            assert (dir_path / "receptor.pdb").exists()
            assert (dir_path / "ligand.sdf").exists() or (dir_path / "ligand.pdb").exists()
            assert (dir_path / "complex.pdb").exists()
            assert (dir_path / "run_openmm_md.py").exists()

            # Verify script contents
            script_text = (dir_path / "run_openmm_md.py").read_text(encoding="utf-8")
            assert "amber14-all.xml" in script_text
            assert "LangevinMiddleIntegrator" in script_text
            assert "MonteCarloBarostat" in script_text
            assert "PME" in script_text

    def test_complex_refinement_service_export_wrapper(self):
        rec_pdb = "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00  0.00           C\n"
        lig_pdb = "HETATM    2  C1  LIG L   1      15.000  10.000  10.000  1.00  0.00           C\n"
        with tempfile.TemporaryDirectory() as tmp_dir:
            pkg = ComplexRefinementService.export_openmm_md_package(
                receptor_pdb=rec_pdb,
                docked_pdb_or_pdbqt=lig_pdb,
                output_dir=tmp_dir,
                ligand_smiles="CC"
            )
            assert pkg["status"] == "SUCCESS"
            assert len(pkg["files_generated"]) >= 4
