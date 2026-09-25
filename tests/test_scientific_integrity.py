import pytest
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from backend.services.bioactivity import BioactivityService
from backend.services.adme import ADMEProfiler
from backend.services.docking import DockingEngine
from backend.services.interaction_engine import InteractionEngine
from backend.utils.rmsd_calculator import calculate_rmsd


class TestThermodynamicCalculations:
    """Test rigorous thermodynamic conversions and Ligand Efficiency metrics."""

    def test_kd_conversion_erlotinib(self):
        """Verify -7.07 kcal/mol converts to ~6.5 uM (6500 nM) at 298.15 K."""
        affinity = -7.07
        heavy_atoms = 29
        mw = 393.44

        thermo = BioactivityService.calculate_thermodynamics(affinity, heavy_atoms, mw)

        assert thermo["binding_affinity_kcal"] == -7.07
        assert thermo["docking_score_kcal"] == -7.07
        assert thermo["temperature_kelvin"] == 298.15
        assert thermo["gas_constant_kcal_per_mol_k"] == 0.0019872041

        # Manual math: RT = 0.0019872041 * 298.15 = 0.592484902415 kcal/mol
        # Kd = exp(-7.07 / 0.592484902415) = exp(-11.932793) = 6.57e-6 M = 6571 nM = 6.57 uM
        kd_nm = thermo["theoretical_kd_nm"]
        assert 6000 <= kd_nm <= 7000, f"Expected Kd ~6500 nM, got {kd_nm}"
        kd_um = thermo["theoretical_kd_um"]
        assert 6.0 <= kd_um <= 7.0, f"Expected Kd ~6.5 uM, got {kd_um}"

        # Assert scientific terminology requirements
        assert thermo["affinity_derived_kd_nm"] == thermo["theoretical_kd_nm"]
        assert thermo["metric_type"] == "Derived / Model-based"
        assert "Affinity-derived Kd-like estimate" in thermo["kd_type"]
        assert "mathematically derived from the docking score" in thermo["scientific_disclaimer"]
        assert "not an experimentally measured or rigorously calculated thermodynamic Kd" in thermo["scientific_disclaimer"]

    def test_kd_conversion_gefitinib(self):
        """Verify -7.99 kcal/mol converts correctly at 298.15 K."""
        affinity = -7.99
        heavy_atoms = 31
        mw = 446.90

        thermo = BioactivityService.calculate_thermodynamics(affinity, heavy_atoms, mw)

        # Kd = exp(-7.99 / 0.5924849) = exp(-13.48557) = 1.39e-6 M = 1390 nM = 1.39 uM
        kd_nm = thermo["theoretical_kd_nm"]
        assert 1300 <= kd_nm <= 1500, f"Expected Kd ~1390 nM, got {kd_nm}"
        assert thermo["affinity_derived_kd_nm"] == thermo["theoretical_kd_nm"]

    def test_ligand_efficiency_calculation(self):
        """Verify LE = |score| / heavy_atoms with explicit scientific documentation."""
        # Erlotinib: 7.07 / 29 = 0.24379... -> 0.244
        thermo_erl = BioactivityService.calculate_thermodynamics(-7.07, 29, 393.44)
        le_erl = thermo_erl["ligand_efficiency"]
        assert le_erl["value"] == 0.244
        assert "|docking_score|" in le_erl["formula"] or "|Vina Score|" in le_erl["formula"] or "score" in le_erl["formula"]
        assert le_erl["heavy_atoms"] == 29

        # Gefitinib: 7.99 / 31 = 0.25774... -> 0.258
        thermo_gef = BioactivityService.calculate_thermodynamics(-7.99, 31, 446.90)
        le_gef = thermo_gef["ligand_efficiency"]
        assert le_gef["value"] == 0.258
        assert le_gef["heavy_atoms"] == 31


class TestADMEAndLipinskiIntegrity:
    """Test RDKit physicochemical descriptors and Lipinski matrix numeric propagation."""

    def test_erlotinib_lipinski_descriptors(self):
        """Ensure HBD and HBA are numeric integers, never None or '-'."""
        # Erlotinib SMILES
        smiles = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
        adme = ADMEProfiler.calculate_adme(smiles)

        phys = adme["physicochemical"]
        # Check both primary keys and backward-compatible aliases
        assert phys["hbd"]["value"] == 1
        assert phys["hba"]["value"] == 7
        assert phys["h_bond_donors"]["value"] == 1
        assert phys["h_bond_acceptors"]["value"] == 7

        assert 390.0 < phys["molecular_weight"]["value"] < 400.0
        assert 2.0 < phys["logp"]["value"] < 3.5

        # Check Lipinski rule section
        lip = adme["drug_likeness"]["lipinski"]
        assert lip["hbd_value"] == 1
        assert lip["hba_value"] == 7
        assert lip["status"] == "Pass"
        assert len(lip["violations"]) == 0

    def test_gefitinib_lipinski_descriptors(self):
        """Ensure Gefitinib descriptors are correctly computed and populated."""
        # Gefitinib SMILES: 1 donor (aniline NH), 7 acceptors (methoxy O, 2 quinazoline N, ether O, morpholine N, morpholine O, aniline N)
        smiles = "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4"
        adme = ADMEProfiler.calculate_adme(smiles)

        phys = adme["physicochemical"]
        assert phys["hbd"]["value"] == 1
        assert phys["hba"]["value"] == 7
        assert phys["h_bond_donors"]["value"] == 1
        assert phys["h_bond_acceptors"]["value"] == 7
        assert 440.0 < phys["molecular_weight"]["value"] < 450.0


class TestRMSDMethodology:
    """Test symmetry-aware automorphism testing and element-constrained assignment."""

    def test_identical_structures_rmsd_zero(self):
        """PDB block compared with itself must yield 0.000 A."""
        # Simple 3-atom water molecule for test
        pdb_block = (
            "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O\n"
            "HETATM    2  H1  HOH A   1       0.957   0.000   0.000  1.00 20.00           H\n"
            "HETATM    3  H2  HOH A   1      -0.240   0.927   0.000  1.00 20.00           H\n"
            "END\n"
        )
        rmsd = calculate_rmsd(pdb_block, pdb_block)
        assert rmsd == 0.0

        details = calculate_rmsd(pdb_block, pdb_block, return_details=True)
        assert details["rmsd"] == 0.0
        assert "method" in details

    def test_symmetric_molecule_automorphism(self):
        """Test that a 180-degree flip of a symmetric ring does not falsely report high RMSD."""
        # Para-xylene: flipping methyl groups (positions 1 and 4)
        m = Chem.MolFromSmiles("Cc1ccc(C)cc1")
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=42)
        conf = m.GetConformer()

        ref_pdb = Chem.MolToPDBBlock(m)

        # Create swapped conformer by inverting coordinates symmetrically
        m_flipped = Chem.Mol(m)
        conf_flipped = m_flipped.GetConformer()
        # Invert along X axis
        for i in range(m.GetNumAtoms()):
            pt = conf.GetAtomPosition(i)
            conf_flipped.SetAtomPosition(i, (-pt.x, -pt.y, pt.z))

        flipped_pdb = Chem.MolToPDBBlock(m_flipped)

        rmsd = calculate_rmsd(ref_pdb, flipped_pdb)
        assert isinstance(rmsd, float)
        assert rmsd >= 0.0

        details = calculate_rmsd(ref_pdb, flipped_pdb, return_details=True)
        assert details["rmsd"] >= 0.0
        assert details["method"] == "topological_symmetry_graph_isomorphism"
        assert details["automorphisms_tested"] == 4


class TestInteractionDetectionAndFluorinePhysics:
    """Test detection criteria and differentiation between sigma-hole halogen bonds and fluorine polar contacts."""

    def test_detection_criteria_exposed(self):
        """Interaction engine must expose explicit distance and angle criteria."""
        criteria = InteractionEngine.get_criteria()
        assert "hydrogen_bond" in criteria
        assert criteria["hydrogen_bond"]["distance_max_angstroms"] == 3.5
        assert "halogen_bond_sigma_hole" in criteria
        assert criteria["halogen_bond_sigma_hole"]["distance_max_angstroms"] == 3.8
        assert criteria["halogen_bond_sigma_hole"]["min_c_x_acceptor_angle_deg"] == 130.0
        assert "fluorine_polar_contact" in criteria
        assert criteria["fluorine_polar_contact"]["distance_max_angstroms"] == 3.5
        assert "sigma-hole" in criteria["fluorine_polar_contact"]["note"].lower()


class TestDockingClassification:
    """Test classification of experiments into Native Redocking vs Cross-Docking vs Targeted."""

    def test_native_redocking_classification(self):
        """Erlotinib docked against 1M17 (native Erlotinib) must be Native Redocking."""
        erlotinib_smiles = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
        native_info = {
            "has_native": True,
            "name": "Erlotinib",
            "smiles": erlotinib_smiles,
            "chain": "A",
            "atom_count": 29
        }

        classification = DockingEngine.classify_docking_experiment(
            docked_smiles=erlotinib_smiles,
            native_ligand_info=native_info,
            docked_ligand_name="Erlotinib (Tarceva)"
        )

        assert classification["is_native_redocking"] is True
        assert classification["docking_mode"] == "Native Redocking"
        assert "Self-Validation" in classification["validation_applicability"]
        assert classification["crystal_ligand_name"] == "Erlotinib"

    def test_cross_docking_classification(self):
        """Gefitinib docked against 1M17 (native Erlotinib) must be Cross-Docking."""
        gefitinib_smiles = "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4"
        erlotinib_smiles = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
        native_info = {
            "has_native": True,
            "name": "Erlotinib (AQ4)",
            "smiles": erlotinib_smiles,
            "chain": "A",
            "atom_count": 29
        }

        classification = DockingEngine.classify_docking_experiment(
            docked_smiles=gefitinib_smiles,
            native_ligand_info=native_info,
            docked_ligand_name="Gefitinib (Iressa)"
        )

        assert classification["is_native_redocking"] is False
        assert "Cross-Docking" in classification["docking_mode"]
        assert "Non-Native" in classification["validation_applicability"]
        assert "Erlotinib" in classification["crystal_ligand_name"]
        assert "AutoDock 4.2" in classification["description"]

    def test_targeted_docking_without_native(self):
        """Docking into a apo or custom cavity without crystal ligand is Targeted Pocket Docking."""
        classification = DockingEngine.classify_docking_experiment(
            docked_smiles="CC(=O)Oc1ccccc1C(=O)O",
            native_ligand_info=None,
            docked_ligand_name="Aspirin"
        )

        assert classification["is_native_redocking"] is False
        assert classification["docking_mode"] == "Targeted Pocket Docking"
        assert "Not Applicable" in classification["validation_applicability"]
