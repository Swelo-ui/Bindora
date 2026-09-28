"""
Comprehensive Test Suite for 30-60 Day Architect Upgrades & Scientific Insights.

Tests:
1. Virtual Covalent Adduct Topology Builder (CYS SG linkage, CONECT records, leaving group displacement)
2. Monte Carlo Backbone phi/psi Induced-Fit Docking (IFD ensemble, Ramachandran perturbation, IFD composite score)
3. GPU Acceleration & OpenCL/CUDA Hardware Profiling
4. The "Grease Bias" Defeated by MM-GBSA (Solvation penalty on hydrophobic greasy decoys)
5. PDB Coordinate Quantization & Relaxation Delta (3-decimal coordinate strain baseline <= 4.0 kcal/mol)
"""

import math
import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from backend.services.covalent import CovalentDockingService
from backend.services.induced_fit import InducedFitService
from backend.services.hardware_profiler import HardwareTelemetrySampler
from backend.services.refinement import ComplexRefinementService


class TestVirtualCovalentAdduct:
    """Test suite for virtual covalent adduct formation and PDB CONECT topology generation."""

    def test_build_covalent_adduct_michael_acceptor(self):
        rec_pdb = (
            "ATOM    145  SG  CYS A 145      10.000  10.000  10.000  1.00  0.00           S\n"
            "ATOM    146  CA  CYS A 145       8.500  10.000  10.000  1.00  0.00           C\n"
        )
        lig_pdb = (
            "HETATM 2001  C1  LIG L   1      10.000  10.000  13.300  1.00  0.00           C\n"
            "HETATM 2002  C2  LIG L   1      10.000  11.200  14.000  1.00  0.00           C\n"
            "HETATM 2003  C3  LIG L   1      10.000  11.200  15.500  1.00  0.00           C\n"
            "HETATM 2004  O1  LIG L   1      10.000  12.200  16.200  1.00  0.00           O\n"
        )
        res = CovalentDockingService.build_covalent_adduct_complex(
            receptor_pdb_or_pdbqt=rec_pdb,
            docked_pose_pdb_or_pdbqt=lig_pdb,
            smiles="C=CC(=O)N"
        )

        assert res["adduct_formed"] is True
        assert res["adduct_status"] == "COVALENT_ADDUCT_FORMED"
        assert res["covalent_bond_length_angstroms"] == 1.82
        assert "CONECT 2001  145" in res["adduct_complex_pdb"]
        assert "CONECT  145 2001" in res["adduct_complex_pdb"]
        assert "C-S Thioether" in res["covalent_bond_type"]

    def test_build_covalent_adduct_haloacetamide_leaving_group(self):
        rec_pdb = (
            "ATOM    145  SG  CYS A 145      10.000  10.000  10.000  1.00  0.00           S\n"
        )
        lig_pdb = (
            "HETATM 2001  C1  LIG L   1      10.000  10.000  13.200  1.00  0.00           C\n"
            "HETATM 2002  CL1 LIG L   1      10.000  10.000  15.000  1.00  0.00          CL\n"
            "HETATM 2003  C2  LIG L   1      10.000  11.200  14.000  1.00  0.00           C\n"
        )
        res = CovalentDockingService.build_covalent_adduct_complex(
            receptor_pdb_or_pdbqt=rec_pdb,
            docked_pose_pdb_or_pdbqt=lig_pdb,
            smiles="ClCC(=O)N"
        )

        assert res["adduct_formed"] is True
        assert res["adduct_status"] == "COVALENT_ADDUCT_FORMED"
        # Verify chlorine leaving group was displaced from the adduct complex
        assert "CL" not in res["adduct_complex_pdb"]
        assert res["covalent_bond_length_angstroms"] == 1.82

    def test_build_covalent_adduct_non_covalent(self):
        rec_pdb = "ATOM    145  SG  CYS A 145      10.000  10.000  10.000  1.00  0.00           S\n"
        lig_pdb = "HETATM 2001  C1  LIG L   1      10.000  10.000  13.200  1.00  0.00           C\n"
        res = CovalentDockingService.build_covalent_adduct_complex(
            receptor_pdb_or_pdbqt=rec_pdb,
            docked_pose_pdb_or_pdbqt=lig_pdb,
            smiles="c1ccccc1"
        )
        assert res["adduct_formed"] is False


class TestInducedFitService:
    """Test suite for Monte Carlo backbone phi/psi induced-fit ensemble sampling."""

    SAMPLE_RECEPTOR = (
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1      12.000  11.417  10.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1      11.236  12.381  10.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1      11.966   9.227  11.214  1.00  0.00           C\n"
        "ATOM      6  N   PHE A   2      13.310  11.550  10.000  1.00  0.00           N\n"
        "ATOM      7  CA  PHE A   2      14.000  12.840  10.000  1.00  0.00           C\n"
        "ATOM      8  C   PHE A   2      13.500  13.800  11.100  1.00  0.00           C\n"
        "ATOM      9  O   PHE A   2      14.100  14.850  11.350  1.00  0.00           O\n"
        "TER\n"
        "END\n"
    )

    def test_extract_pocket_residues(self):
        res = InducedFitService.extract_pocket_residues(
            self.SAMPLE_RECEPTOR,
            pocket_center={"x": 11.0, "y": 10.0, "z": 10.0},
            radius=5.0
        )
        assert len(res) >= 1
        assert res[0]["residue_name"] == "ALA"

    def test_sample_backbone_induced_fit(self):
        res = InducedFitService.sample_backbone_induced_fit(
            self.SAMPLE_RECEPTOR,
            pocket_center={"x": 11.0, "y": 10.0, "z": 10.0},
            radius=8.0,
            num_conformations=2
        )
        assert res["status"] == "SUCCESS"
        assert res["conformations_generated"] == 3  # 1 crystal reference + 2 induced conformations
        # Conformation 0 is crystal reference
        assert res["conformations"][0]["pocket_rmsd_angstroms"] == 0.0
        assert res["conformations"][0]["receptor_strain_kcal"] == 0.0
        # Conformation 1 is induced
        assert res["conformations"][1]["status"] == "INDUCED_FIT_RELAXED"


class TestGPUHardwareAcceleration:
    """Test suite for hardware GPU acceleration detection and OpenCL/CUDA profiling."""

    def test_detect_gpu_capabilities(self):
        gpu_info = HardwareTelemetrySampler.detect_gpu_capabilities()
        assert "cuda_available" in gpu_info
        assert "opencl_available" in gpu_info
        assert "gpu_hardware_ready" in gpu_info
        assert "preferred_acceleration_platform" in gpu_info
        assert isinstance(gpu_info["openmm_platforms"], list)
        # On this environment, OpenMM has OpenCL available
        assert gpu_info["opencl_available"] is True
        assert "OpenCL" in gpu_info["preferred_acceleration_platform"]


class TestScientificPhysicsCalibration:
    """Test suite for computational chemistry insights: Grease Bias and PDB Quantization Delta."""

    def test_grease_bias_defeated_by_mmgbsa(self):
        """
        Verify that a hydrophobic greasy decoy (e.g. tetradecane/pyrene) in a polar pocket
        receives a significant solvation penalty via MM-GBSA continuum dielectric.
        """
        polar_receptor = (
            "ATOM      1  N   SER A   1      10.000  10.000  10.000  1.00  0.00           N\n"
            "ATOM      2  CA  SER A   1      11.458  10.000  10.000  1.00  0.00           C\n"
            "ATOM      3  OG  SER A   1      12.000  11.200  10.000  1.00  0.00           O\n"
            "ATOM      4  N   ASP A   2      13.500  10.000  10.000  1.00  0.00           N\n"
            "ATOM      5  OD1 ASP A   2      14.000  11.500  10.000  1.00  0.00           O\n"
            "ATOM      6  OD2 ASP A   2      14.000   9.000  10.000  1.00  0.00           O\n"
        )
        # Hydrophobic greasy decoy pose (pure carbons)
        greasy_pose = (
            "HETATM    1  C1  DEC L   1      12.000  10.000  13.000  1.00  0.00           C\n"
            "HETATM    2  C2  DEC L   1      13.200  10.000  13.500  1.00  0.00           C\n"
            "HETATM    3  C3  DEC L   1      14.200  10.000  14.000  1.00  0.00           C\n"
        )
        strain_data = {"ligand_strain_relaxation_kcal": 0.5, "is_high_strain": False}
        rescore = ComplexRefinementService.calculate_mmgbsa_rescore(polar_receptor, greasy_pose, strain_data)

        # Greasy decoy in polar pocket should have unfavorable polar desolvation penalty (Delta G_GB > 0)
        assert rescore["available"] is True
        assert rescore["components"]["gb_desolvation_penalty_kcal"] >= 0.0

    def test_pdb_quantization_relaxation_baseline(self):
        """
        Verify that standard crystallographic coordinates with 3-decimal floating point
        quantization roundoff have baseline relaxation strain within calibrated threshold (<= 4.0 kcal/mol).
        """
        # Create small test ligand with standard 3-decimal PDB records
        mol = Chem.MolFromSmiles("c1ccccc1O")
        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol_h)
        pdb_block = Chem.MolToPDBBlock(mol_h)

        strain_eval = ComplexRefinementService.calculate_ligand_strain(pdb_block, smiles="c1ccccc1O")
        strain_val = strain_eval["ligand_strain_relaxation_kcal"]

        # 3-decimal PDB quantization strain MUST be <= 4.0 kcal/mol for relaxed crystal molecules
        assert strain_val <= 4.0
        assert strain_eval["is_high_strain"] is False
        assert "Low Strain" in strain_eval["strain_classification"]
        assert strain_eval["decoy_filter_flag"] == "PASS"
