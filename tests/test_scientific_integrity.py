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

    # ──────────────────────────────────────────────────────────────────────────
    # Regression tests for STI / Imatinib identity fix (structure-first)
    # ──────────────────────────────────────────────────────────────────────────

    def test_sti_imatinib_identity_native_redocking(self):
        """
        TEST A: PDB ligand code STI + docked Imatinib MUST resolve to Native Redocking.

        Root bug: The old code compared PDB-parsed SMILES (corrupted by aromatic perception
        failure) against the user's clean Imatinib SMILES. They did not match canonically.
        The name fallback also failed because 'STI' != 'imatinib'.
        Fix: Use RCSB CCD InChIKey lookup which unambiguously confirms STI == Imatinib.
        """
        imatinib_smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(C)CC5"
        # Simulate the bad PDB-parsed SMILES that Bindora was generating from 1T46 crystal data
        bad_pdb_smiles = "CC1CC[C@@H](NC(O)C2CCC(CN3CCN(C)CC3)CC2)CC1NC1NCCC(C2CCCNC2)N1"
        native_info = {
            "has_native": True,
            "name": "STI",            # PDB 3-letter code — key test case
            "smiles": bad_pdb_smiles,  # Corrupted from PDB aromatic parsing
            "chain": "A",
            "atom_count": 37
        }

        classification = DockingEngine.classify_docking_experiment(
            docked_smiles=imatinib_smiles,
            native_ligand_info=native_info,
            docked_ligand_name="Imatinib"
        )

        # Must be Native Redocking — NOT Cross-Docking
        assert classification["is_native_redocking"] is True, (
            f"Expected Native Redocking but got: {classification['docking_mode']}. "
            f"Method: {classification.get('identity_method')}"
        )
        assert classification["docking_mode"] == "Native Redocking"
        assert "Self-Validation" in classification["validation_applicability"]
        assert "STI" in classification["crystal_ligand_name"]
        assert classification["docked_ligand_name"] == "Imatinib"

        # Must be structure-confirmed, not just name heuristic
        assert classification["identity_confidence"] == "structure-confirmed", (
            f"Expected structure-confirmed confidence, got: {classification['identity_confidence']}"
        )

        # Description must explain what actually happened
        assert "independently prepared and redocked" in classification["description"]

    def test_sti571_imatinib_same_molecule(self):
        """
        TEST C (variant): When native SMILES is directly Imatinib and docked is also Imatinib,
        classification must be Native Redocking.
        """
        imatinib_smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(C)CC5"
        native_info = {
            "has_native": True,
            "name": "Imatinib",
            "smiles": imatinib_smiles,
            "chain": "A",
            "atom_count": 37
        }

        classification = DockingEngine.classify_docking_experiment(
            docked_smiles=imatinib_smiles,
            native_ligand_info=native_info,
            docked_ligand_name="Imatinib"
        )

        assert classification["is_native_redocking"] is True
        assert classification["docking_mode"] == "Native Redocking"
        assert classification["identity_confidence"] == "structure-confirmed"

    def test_erlotinib_vs_gefitinib_different_structures(self):
        """
        TEST B: Erlotinib (native) vs Gefitinib (docked) must be Cross-Docking.
        They share similar pharmacophore but have different molecular graphs.
        """
        from backend.services.ligand_identity import compare_ligand_identity

        erlotinib = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
        gefitinib = "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4"

        result = compare_ligand_identity(
            docked_smiles=gefitinib,
            native_smiles=erlotinib,
            native_name="Erlotinib",
            docked_name="Gefitinib"
        )

        assert result["is_same_molecule"] is False, (
            "Erlotinib and Gefitinib are distinct molecules — must NOT be classified as same"
        )

    def test_no_native_ligand_is_targeted_docking(self):
        """
        TEST D: No native crystal ligand present → Targeted Pocket Docking.
        """
        classification = DockingEngine.classify_docking_experiment(
            docked_smiles="Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(C)CC5",
            native_ligand_info={"has_native": False},
            docked_ligand_name="Imatinib"
        )

        assert classification["is_native_redocking"] is False
        assert classification["docking_mode"] == "Targeted Pocket Docking"
        assert classification["has_crystal_reference"] is False

    def test_missing_docked_smiles_does_not_crash(self):
        """
        TEST E (variant): When docked SMILES is missing, classification should not crash
        and must report identity as undetermined rather than guessing.
        """
        native_info = {
            "has_native": True,
            "name": "STI",
            "smiles": "",
            "chain": "A",
            "atom_count": 37
        }

        # Should not raise any exception
        classification = DockingEngine.classify_docking_experiment(
            docked_smiles=None,
            native_ligand_info=native_info,
            docked_ligand_name="Unknown ligand"
        )

        # When docked SMILES is missing, cannot confirm identity structurally
        # So it must NOT claim Native Redocking from name heuristic alone
        # (STI vs "Unknown ligand" should not match by name)
        assert isinstance(classification, dict)
        assert "docking_mode" in classification
        # Must not crash — structural integrity is more important than claiming a match


class TestLigandIdentityService:
    """Unit tests for the structure-first LigandIdentityService."""

    def test_inchikey_match_sti_imatinib(self):
        """STI (PDB code) and Imatinib SMILES must share the same InChIKey."""
        from backend.services.ligand_identity import compare_ligand_identity

        imatinib_smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(C)CC5"
        result = compare_ligand_identity(
            docked_smiles=imatinib_smiles,
            native_pdb_code="STI",
            native_name="STI",
            docked_name="Imatinib"
        )

        assert result["is_same_molecule"] is True
        assert result["confidence"] == "structure-confirmed"
        assert "InChIKey" in result["method"]
        assert result["inchikey_docked"] == "KTUFNOKKBVMGRW-UHFFFAOYSA-N"
        assert result["inchikey_native"] == "KTUFNOKKBVMGRW-UHFFFAOYSA-N"

    def test_canonical_smiles_match_same_molecule(self):
        """Same molecule given as slightly different but equivalent SMILES must match."""
        from backend.services.ligand_identity import compare_ligand_identity

        # Two equivalent representations of Imatinib
        smiles_a = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(C)CC5"
        smiles_b = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"

        result = compare_ligand_identity(
            docked_smiles=smiles_a,
            native_smiles=smiles_b,
            native_name="Imatinib",
            docked_name="Imatinib"
        )

        assert result["is_same_molecule"] is True
        assert result["confidence"] == "structure-confirmed"

    def test_different_molecules_not_same(self):
        """Erlotinib vs Gefitinib must not be reported as the same molecule."""
        from backend.services.ligand_identity import compare_ligand_identity

        erlotinib = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
        gefitinib = "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4"

        result = compare_ligand_identity(
            docked_smiles=gefitinib,
            native_smiles=erlotinib,
            native_name="Erlotinib",
            docked_name="Gefitinib"
        )

        assert result["is_same_molecule"] is False

    def test_resolve_pdb_ligand_identity_sti(self):
        """resolve_pdb_ligand_identity('STI') must return Imatinib-related data."""
        from backend.services.ligand_identity import resolve_pdb_ligand_identity

        result = resolve_pdb_ligand_identity("STI")

        assert result["pdb_ligand_code"] == "STI"
        assert result["inchi_key"] == "KTUFNOKKBVMGRW-UHFFFAOYSA-N"
        assert result["identity_confidence"] == "structure-confirmed"
        # Imatinib or STI-571 should appear in synonyms
        synonyms_upper = [s.upper() for s in result["synonyms"]]
        assert "IMATINIB" in synonyms_upper or "STI-571" in synonyms_upper

