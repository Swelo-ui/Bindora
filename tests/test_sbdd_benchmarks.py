import pytest
from backend.services.refinement import ComplexRefinementService
from backend.services.consensus import ConsensusScoringService
from backend.services.docking import DockingEngine
from backend.services.fetcher import StructureFetcher

def test_decoy_discrimination_pentacene():
    """Verify that a purely hydrophobic flat polyaromatic decoy (Pentacene) with 0 polar contacts is flagged as a decoy."""
    pentacene_smiles = "c1ccc2cc3cc4ccccc4cc3cc2c1"
    
    # Minimal receptor pocket with polar residues (e.g. EGFR hinge or HIV-1 protease)
    sample_receptor_pdb = """ATOM      1  N   MET A 769      21.000   1.000  50.000  1.00 20.00           N
ATOM      2  O   THR A 766      22.500   2.500  51.500  1.00 20.00           O
ATOM      3  N   LYS A 721      20.500  -2.000  53.000  1.00 20.00           N
ATOM      4  O   ASP A 831      24.000  -1.500  52.000  1.00 20.00           O
END
"""
    # Pentacene conformer (pure carbons and hydrogens)
    pentacene_pdbqt = """ATOM      1  C   LIG     1      21.500   0.500  50.500  0.00  0.00    +0.000 C
ATOM      2  C   LIG     1      22.000   1.200  51.000  0.00  0.00    +0.000 C
ATOM      3  C   LIG     1      22.500   1.800  51.500  0.00  0.00    +0.000 C
ATOM      4  C   LIG     1      23.000   2.200  52.000  0.00  0.00    +0.000 C
ATOM      5  C   LIG     1      23.500   2.600  52.500  0.00  0.00    +0.000 C
ATOM      6  C   LIG     1      24.000   3.000  53.000  0.00  0.00    +0.000 C
ATOM      7  C   LIG     1      24.500   3.400  53.500  0.00  0.00    +0.000 C
ATOM      8  C   LIG     1      25.000   3.800  54.000  0.00  0.00    +0.000 C
END
"""
    mmgbsa = ComplexRefinementService.calculate_mmgbsa_rescore(
        receptor_pdb=sample_receptor_pdb,
        docked_pdb_or_pdbqt=pentacene_pdbqt,
        smiles=pentacene_smiles
    )
    
    assert mmgbsa["available"] is True
    assert mmgbsa["polar_contacts_count"] == 0
    assert mmgbsa["is_grease_ball_decoy"] is True
    assert mmgbsa["decoy_filter_verdict"] == "FLAGGED_GREASY_DECOY"
    assert mmgbsa["components"]["opportunistic_decoy_penalty_kcal"] > 0.0

def test_metal_coordination_physics():
    """Verify that a coordinating heteroatom within 1.8-2.6 A of a catalytic Zn2+ receives a coordination bonus."""
    receptor_zn_pdb = """HETATM 2472  ZN   ZN A 261      37.595  16.545 -14.676  1.00 20.00     2.000  Zn
ATOM      1  N   HIS A  94      35.800  15.500 -14.000  1.00 20.00           N
ATOM      2  N   HIS A  96      38.500  18.000 -14.200  1.00 20.00           N
END
"""
    # Coordinated sulfonamide oxygen at 2.1 A from Zn
    coordinated_ligand_pdbqt = """ATOM      1  OA  AZM     1      37.600  16.500 -12.576  0.00  0.00    -0.450 OA
ATOM      2  NA  AZM     1      36.500  15.500 -11.800  0.00  0.00    -0.350 NA
ATOM      3  SA  AZM     1      37.000  15.000 -10.500  0.00  0.00    +0.800 SA
ATOM      4  C   AZM     1      38.000  14.000  -9.500  0.00  0.00    +0.100 C
ATOM      5  C   AZM     1      38.500  13.000  -8.500  0.00  0.00    +0.000 C
ATOM      6  N   AZM     1      39.000  12.000  -7.500  0.00  0.00    -0.200 N
ATOM      7  O   AZM     1      39.500  11.000  -6.500  0.00  0.00    -0.300 O
ATOM      8  C   AZM     1      40.000  10.000  -5.500  0.00  0.00    +0.000 C
END
"""
    mmgbsa = ComplexRefinementService.calculate_mmgbsa_rescore(
        receptor_pdb=receptor_zn_pdb,
        docked_pdb_or_pdbqt=coordinated_ligand_pdbqt,
        smiles="CC(=O)Nc1nnc(s1)S(=O)(=O)N"
    )
    
    assert mmgbsa["available"] is True
    assert mmgbsa["components"]["metal_coordination_bonus_kcal"] <= -5.0

def test_torsional_entropy_penalty():
    """Verify that high-torsion flexible ligands (N_rot > 8) receive a torsional entropy penalty."""
    # Molecule with 12 rotatable bonds
    high_torsion_smiles = "CCCCCCCCCCCC(=O)O"
    
    rec_pdb = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  O   ALA A   1       2.000   2.000   2.000  1.00 20.00           O
END
"""
    lig_pdbqt = """ATOM      1  C   LIG     1       1.000   1.000   1.000  0.00  0.00    +0.000 C
ATOM      2  C   LIG     1       2.000   1.000   1.000  0.00  0.00    +0.000 C
ATOM      3  C   LIG     1       3.000   1.000   1.000  0.00  0.00    +0.000 C
ATOM      4  C   LIG     1       4.000   1.000   1.000  0.00  0.00    +0.000 C
ATOM      5  C   LIG     1       5.000   1.000   1.000  0.00  0.00    +0.000 C
ATOM      6  C   LIG     1       6.000   1.000   1.000  0.00  0.00    +0.000 C
ATOM      7  C   LIG     1       7.000   1.000   1.000  0.00  0.00    +0.000 C
ATOM      8  O   LIG     1       8.000   1.000   1.000  0.00  0.00    -0.300 OA
END
"""
    mmgbsa = ComplexRefinementService.calculate_mmgbsa_rescore(
        receptor_pdb=rec_pdb,
        docked_pdb_or_pdbqt=lig_pdbqt,
        smiles=high_torsion_smiles
    )
    
    assert mmgbsa["available"] is True
    assert mmgbsa["rotatable_bonds"] >= 10
    assert mmgbsa["components"]["torsional_entropy_penalty_kcal"] > 0.0

def test_consensus_matrix_demotes_decoys():
    """Verify that consensus scoring relegates decoys to the bottom of the ranking."""
    poses = [
        {
            "mode": 1,
            "affinity_kcal": -9.5,
            "decoy_filter_flag": "FLAGGED_GREASY_DECOY",
            "mmgbsa": {"is_grease_ball_decoy": True, "decoy_reason": "Zero polar contacts", "mmgbsa_delta_g_kcal": 2.5}
        },
        {
            "mode": 2,
            "affinity_kcal": -7.2,
            "decoy_filter_flag": "PASS",
            "mmgbsa": {"is_grease_ball_decoy": False, "mmgbsa_delta_g_kcal": -8.5}
        }
    ]
    ranked = ConsensusScoringService.compute_pose_consensus(poses)
    
    assert ranked[0]["mode"] == 2
    assert ranked[0]["consensus_rank"] == 1
    assert ranked[1]["mode"] == 1
    assert ranked[1]["consensus_rank"] == 2
    assert ranked[1]["consensus_confidence"] == "DECOY_GREASE_BALL"
