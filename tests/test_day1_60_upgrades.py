import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from backend.services.docking import DockingEngine
from backend.services.interaction_engine import InteractionEngine
from backend.services.refinement import ComplexRefinementService
from backend.services.bioactivity import BioactivityService
from backend.config import VINA_GPU_EXE, UNIDOCK_EXE, USE_GPU_DOCKING

def test_true_ligand_strain_and_coordinate_preservation():
    """Verify that calculate_ligand_strain evaluates true 3D strain without overwriting coordinates when SMILES is passed."""
    smi = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)

    # 1. When optimized, strain should be low (<= 4.0 kcal/mol)
    strain_res = ComplexRefinementService.calculate_ligand_strain(pdb_block, smiles=smi)
    assert strain_res["status"] == "Strain Calculated Successfully"
    assert strain_res["ligand_strain_relaxation_kcal"] is not None
    assert strain_res["ligand_strain_relaxation_kcal"] <= 4.0
    assert strain_res["is_high_strain"] is False
    assert strain_res["decoy_filter_flag"] == "PASS"

    # 2. When coordinates are strained via unfavorable dihedral angle, strain should be detected and flagged
    AllChem.SetDihedralDeg(mol.GetConformer(), 0, 1, 3, 4, 90.0)
    distorted_pdb = Chem.MolToPDBBlock(mol)

    strained_res = ComplexRefinementService.calculate_ligand_strain(distorted_pdb, smiles=smi)
    assert strained_res["ligand_strain_relaxation_kcal"] > 6.0
    assert strained_res["is_high_strain"] is True
    assert strained_res["decoy_filter_flag"] == "FLAG_HIGH_STRAIN_DECOY"
    assert "High Ligand Strain Alert" in strained_res["strain_warning"]

def test_mmgbsa_physics_rescoring():
    """Verify MM-GBSA implicit solvent rescoring computes vdW, electrostatics, desolvation, and SA terms."""
    rec_pdb = (
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  C   ALA A   1      12.000  11.362  10.000  1.00 20.00           C\n"
        "ATOM      4  O   ALA A   1      11.246  12.327  10.000  1.00 20.00           O\n"
        "END\n"
    )
    lig_mol = Chem.MolFromSmiles("CCO")
    lig_mol = Chem.AddHs(lig_mol)
    AllChem.EmbedMolecule(lig_mol, randomSeed=42)
    lig_pdb = Chem.MolToPDBBlock(lig_mol)

    mmgbsa = ComplexRefinementService.calculate_mmgbsa_rescore(rec_pdb, lig_pdb)
    assert mmgbsa["available"] is True
    assert "mmgbsa_delta_g_kcal" in mmgbsa
    assert "components" in mmgbsa
    assert "vdw_interaction_kcal" in mmgbsa["components"]
    assert "electrostatic_interaction_kcal" in mmgbsa["components"]
    assert "gb_desolvation_penalty_kcal" in mmgbsa["components"]
    assert "sa_nonpolar_hydrophobic_kcal" in mmgbsa["components"]

def test_catalytic_metal_coordination_detection():
    """Verify that interaction engine detects catalytic metal coordination (Zn2+, Mg2+) with ligand donors."""
    rec_with_zn = (
        "HETATM    1 ZN    ZN A 301      10.000  10.000  10.000  1.00 20.00          ZN\n"
        "ATOM      2  CA  HIS A 100      10.000  12.000  10.000  1.00 20.00           C\n"
        "END\n"
    )
    # Ligand oxygen 2.1 A from Zn
    lig_coord = (
        "HETATM    1  O1  LIG     1      10.000  10.000  12.100  1.00 20.00           O\n"
        "END\n"
    )
    contacts = InteractionEngine.analyze(rec_with_zn, lig_coord)
    assert contacts["total_metal_coordination_count"] == 1
    assert len(contacts["metal_coordinations"]) == 1
    coord = contacts["metal_coordinations"][0]
    assert coord["metal"] == "ZN"
    assert coord["distance"] == 2.1
    assert "Direct Coordination (Inner Sphere)" in coord["type"]
    assert "ZN 301:A" in contacts["interacting_residues"]

def test_automated_conserved_bridging_water_detection():
    """Verify that prepare_receptor detects and retains conserved bridging waters making >= 2 H-bonds to protein."""
    # Active site with two polar residues (donor and acceptor) and a bridging water
    pdb_content = (
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  C   ALA A   1      12.000  11.000  10.000  1.00 20.00           C\n"
        "ATOM      4  O   ALA A   1      10.000  12.800  10.000  1.00 20.00           O\n"
        "ATOM      5  N   GLY A   2      12.800  10.000  10.000  1.00 20.00           N\n"
        # Water at (11.0, 11.5, 10.0) making contacts to O (2.2 A) and N (2.3 A)
        "HETATM    6  O   HOH A 201      11.000  11.500  10.000  1.00 20.00           O\n"
        "END\n"
    )
    prep = DockingEngine.prepare_receptor(pdb_content, auto_conserved_waters=True)
    assert prep["prep_log"]["conserved_structural_waters_retained"] >= 1
    assert "OA" in prep["pdbqt_text"]

def test_bioactivity_kd_confidence_intervals():
    """Verify that BioactivityService outputs realistic 95% confidence intervals for Kd and pKd."""
    thermo = BioactivityService.calculate_thermodynamics(-10.5, 30, 450.0)
    assert "kd_confidence_interval_95" in thermo
    ci = thermo["kd_confidence_interval_95"]
    assert ci["kd_lower_nm"] < thermo["theoretical_kd_nm"] < ci["kd_upper_nm"]
    assert ci["pkd_lower"] < thermo["pkd"] < ci["pkd_upper"]
    assert "formatted_range" in ci
    assert "+/- 1.5" in ci["formatted_range"]
    assert "95%" in ci["confidence_level"]
