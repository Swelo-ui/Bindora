import pytest
import numpy as np
from rdkit import Chem

from backend.services.biophysical_prep import BiophysicalReceptorPreparer
from backend.services.interaction_engine import InteractionEngine
from backend.services.interaction_diagram import InteractionDiagramGenerator
from backend.services.pharmacophore import PharmacophoreService
from backend.services.batch_manager import BatchScreeningManager
from backend.services.docking import DockingEngine
from backend.services.refinement import ComplexRefinementService

def test_biophysical_receptor_polar_hydrogens_and_charges():
    # Synthetic dipeptide: SER - LYS
    protein_lines = [
        "ATOM      1  N   SER A   1      10.000  10.000  10.000  1.00 20.00           N",
        "ATOM      2  CA  SER A   1      11.458  10.000  10.000  1.00 20.00           C",
        "ATOM      3  C   SER A   1      12.000  11.362  10.000  1.00 20.00           C",
        "ATOM      4  O   SER A   1      11.246  12.327  10.000  1.00 20.00           O",
        "ATOM      5  CB  SER A   1      12.009   9.219   8.796  1.00 20.00           C",
        "ATOM      6  OG  SER A   1      11.500   7.890   8.796  1.00 20.00           O",
        "ATOM      7  N   LYS A   2      13.300  11.500  10.000  1.00 20.00           N",
        "ATOM      8  CA  LYS A   2      14.000  12.700  10.000  1.00 20.00           C",
        "ATOM      9  C   LYS A   2      15.500  12.500  10.000  1.00 20.00           C",
        "ATOM     10  O   LYS A   2      16.100  11.400  10.000  1.00 20.00           O",
        "ATOM     11  CB  LYS A   2      13.500  13.700   8.900  1.00 20.00           C",
        "ATOM     12  CG  LYS A   2      13.900  15.100   9.200  1.00 20.00           C",
        "ATOM     13  CD  LYS A   2      13.400  16.100   8.100  1.00 20.00           C",
        "ATOM     14  CE  LYS A   2      13.800  17.500   8.400  1.00 20.00           C",
        "ATOM     15  NZ  LYS A   2      13.300  18.500   7.300  1.00 20.00           N",
    ]

    pdb_lines, pdbqt_lines, stats = BiophysicalReceptorPreparer.prepare(protein_lines, pH=7.4)

    assert stats["polar_hydrogens_added"] > 0
    assert len(pdbqt_lines) > len(protein_lines)  # Polar hydrogens added

    # Verify polar hydrogens are typed as HD with positive partial charge
    hd_lines = [l for l in pdbqt_lines if l.strip().endswith("HD")]
    assert len(hd_lines) >= 4  # Serine OG-HG, Backbone NH, Lysine NZ-3Hs

    for hdl in hd_lines:
        charge = float(hdl[70:76])
        assert charge > 0.15  # Positive partial charge for donor hydrogen

    # Verify acceptor oxygens are typed as OA with negative partial charge
    oa_lines = [l for l in pdbqt_lines if l.strip().endswith("OA")]
    assert len(oa_lines) >= 3  # Serine O, OG, Lysine O
    for oal in oa_lines:
        charge = float(oal[70:76])
        assert charge < -0.40

def test_interaction_engine_plip_geometry_and_repulsion():
    # Receptor with an Aspartate (anionic) and a backbone NH donor
    rec_pdb = (
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  H   ALA A   1      10.000  10.000   9.000  1.00 20.00           H\n"
        "ATOM      3  CA  ALA A   1      11.458  10.000  10.000  1.00 20.00           C\n"
        "ATOM      4  O   ALA A   1      12.000  11.362  10.000  1.00 20.00           O\n"
        "ATOM      5  OD1 ASP A   2       5.000   5.000   5.000  1.00 20.00           O\n"
        "ATOM      6  OD2 ASP A   2       5.000   5.000   6.200  1.00 20.00           O\n"
    )

    # 1. Negative Control: Two pure carbonyl oxygens facing each other -> Repulsion, NO H-bond
    lig_repulsion = "HETATM    1  O1  LIG     1      12.000  13.200  10.000  1.00 20.00           O\n"
    res_rep = InteractionEngine.analyze(rec_pdb, lig_repulsion)
    assert len(res_rep["hydrogen_bonds"]) == 0

    # 2. Positive Control: Collinear Donor-H...Acceptor (N-H...O) -> Strong directional H-bond
    lig_hbond = "HETATM    1  O1  LIG     1      10.000  10.000   7.200  1.00 20.00           O\n"
    res_hb = InteractionEngine.analyze(rec_pdb, lig_hbond)
    assert len(res_hb["hydrogen_bonds"]) == 1
    hb = res_hb["hydrogen_bonds"][0]
    assert hb["angle_deg"] >= 120.0
    assert hb["distance"] <= 3.2
    assert hb["ligand_atom_idx"] == 0

    # 3. Salt Bridge: Cationic ligand nitrogen interacting with Aspartate carboxylate
    lig_cation = "HETATM    1  N1  LIG     1       5.000   5.000   2.500  1.00 20.00           N\n"
    res_sb = InteractionEngine.analyze(rec_pdb, lig_cation)
    assert len(res_sb["salt_bridges"]) == 1
    sb = res_sb["salt_bridges"][0]
    assert sb["type"] == "Salt Bridge"
    assert "ASP" in sb["res_name"]

def test_2d_interaction_diagram_exact_atom_mapping():
    # 2D diagram generator must map contact lines to exact interacting ligand atom index
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    interactions = {
        "hydrogen_bonds": [
            {
                "res_name": "SER",
                "res_num": 530,
                "chain": "A",
                "distance": 2.8,
                "ligand_atom_idx": 3  # Specific carbonyl oxygen
            }
        ],
        "hydrophobic_contacts": [
            {
                "res_name": "TYR",
                "res_num": 385,
                "chain": "A",
                "distance": 3.7,
                "ligand_atom_idx": 6  # Specific aromatic ring carbon
            }
        ]
    }

    svg = InteractionDiagramGenerator.generate_diagram_svg(smiles, interactions)
    assert "<svg" in svg
    assert "SER 530:A" in svg
    assert "TYR 385:A" in svg
    assert "LigPlot-style 2D Map" in svg

def test_3d_spatial_pharmacophore_model():
    actives = [
        {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "chembl_id": "CHEMBL25"},
        {"smiles": "OC(=O)c1ccccc1O", "chembl_id": "CHEMBL42"},
        {"smiles": "CC(=O)Nc1ccc(O)cc1", "chembl_id": "CHEMBL112"}
    ]

    profile = PharmacophoreService.build_consensus_profile(actives)
    assert "consensus_spheres" in profile
    assert len(profile["consensus_spheres"]) > 0
    assert "distance_constraints" in profile
    assert profile["model_type"].startswith("True 3D Spatial")

    # Match Aspirin against the profile
    cand_match = PharmacophoreService.match_candidate("CC(=O)Oc1ccccc1C(=O)O", profile)
    assert "match_score_pct" in cand_match
    assert cand_match["match_score_pct"] >= 70.0
    assert "spatial_sphere_alignment" in cand_match
    assert cand_match["spatial_rmsd_angstroms"] is not None

def test_batch_manager_async_lifecycle():
    manager = BatchScreeningManager.get_instance()
    job_id = manager.start_batch_job(
        receptor_pdbqt="ATOM      1  N   ALA A   1\n",
        receptor_pdb="ATOM      1  N   ALA A   1\n",
        pocket_center={"x": 0.0, "y": 0.0, "z": 0.0},
        pocket_size={"x": 20.0, "y": 20.0, "z": 20.0},
        ligand_list=[
            {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
            {"name": "Phenol", "smiles": "c1ccccc1O"}
        ],
        exhaustiveness=1
    )

    assert job_id.startswith("batch_")
    job = manager.get_job(job_id)
    assert job is not None
    assert job["total"] == 2
    assert job["status"] in ("queued", "running", "completed")

def test_complex_pose_refinement():
    rec_pdb = (
        "ATOM      1  N   SER A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  SER A   1      11.458  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  C   SER A   1      12.000  11.362  10.000  1.00 20.00           C\n"
        "ATOM      4  O   SER A   1      11.246  12.327  10.000  1.00 20.00           O\n"
    )
    docked_pose = (
        "HETATM    1  C1  LIG     1      10.000  10.000   5.000  1.00 20.00           C\n"
        "HETATM    2  O1  LIG     1      10.000  10.000   6.200  1.00 20.00           O\n"
    )

    refine_res = ComplexRefinementService.refine_pose(rec_pdb, docked_pose, smiles="CC(=O)O")
    assert "method" in refine_res
    assert refine_res["mmgbsa_dG_kcal"] is not None
    assert "status" in refine_res

def test_ptm_and_nonstandard_residue_preparation():
    # Tripeptide with MSE (Selenomethionine) and SEP (Phosphoserine)
    pdb_content = (
        "ATOM      1  N   MSE A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  MSE A   1      11.458  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  C   MSE A   1      12.000  11.362  10.000  1.00 20.00           C\n"
        "ATOM      4  O   MSE A   1      11.246  12.327  10.000  1.00 20.00           O\n"
        "ATOM      5  CB  MSE A   1      12.009   9.219   8.796  1.00 20.00           C\n"
        "ATOM      6  CG  MSE A   1      11.800   7.700   8.700  1.00 20.00           C\n"
        "ATOM      7  SE  MSE A   1      12.500   6.700   7.200  1.00 20.00          SE\n"
        "ATOM      8  CE  MSE A   1      11.500   5.100   7.400  1.00 20.00           C\n"
        "HETATM    9  N   SEP A   2      13.300  11.500  10.000  1.00 20.00           N\n"
        "HETATM   10  CA  SEP A   2      14.000  12.700  10.000  1.00 20.00           C\n"
        "HETATM   11  C   SEP A   2      15.500  12.500  10.000  1.00 20.00           C\n"
        "HETATM   12  O   SEP A   2      16.100  11.400  10.000  1.00 20.00           O\n"
        "HETATM   13  CB  SEP A   2      13.500  13.700   8.900  1.00 20.00           C\n"
        "HETATM   14  OG  SEP A   2      13.900  15.000   9.100  1.00 20.00           O\n"
        "HETATM   15  P   SEP A   2      13.200  16.300   8.400  1.00 20.00           P\n"
        "HETATM   16  O1P SEP A   2      11.700  16.100   8.400  1.00 20.00           O\n"
        "HETATM   17  O2P SEP A   2      13.800  17.500   9.200  1.00 20.00           O\n"
        "HETATM   18  O3P SEP A   2      13.600  16.400   6.900  1.00 20.00           O\n"
    )

    rec = DockingEngine.prepare_receptor(pdb_content)
    assert rec["prep_log"]["protein_atoms_retained"] >= 18
    pdbqt = rec["pdbqt_text"]

    # Verify MSE SE atom is typed as SA with negative partial charge
    se_lines = [l for l in pdbqt.splitlines() if " SE " in l or l.strip().endswith("SA")]
    assert len(se_lines) >= 1
    for sel in se_lines:
        assert float(sel[70:76]) < 0.0

    # Verify SEP phosphate P and oxygens
    p_lines = [l for l in pdbqt.splitlines() if " P " in l and "SEP" in l]
    assert len(p_lines) >= 1
    assert float(p_lines[0][70:76]) > 1.0  # Positive phosphorus

    o1p_lines = [l for l in pdbqt.splitlines() if "O1P" in l]
    assert len(o1p_lines) >= 1
    assert float(o1p_lines[0][70:76]) < -0.5  # Ionized phosphate oxygen

def test_catalytic_metal_retention_and_strict_columns():
    pdb_with_metals = (
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  C   ALA A   1      12.000  11.362  10.000  1.00 20.00           C\n"
        "ATOM      4  O   ALA A   1      11.246  12.327  10.000  1.00 20.00           O\n"
        # Catalytic Zinc near the pocket centroid (~ 2.0 A)
        "HETATM    5  ZN   ZN A 101      11.000  11.000  11.000  1.00 20.00          ZN\n"
        # Distant buffer magnesium far away (> 40 A)
        "HETATM    6  MG   MG A 102      55.000  65.000  75.000  1.00 20.00          MG\n"
    )

    rec = DockingEngine.prepare_receptor(pdb_with_metals)
    assert rec["prep_log"]["catalytic_metals_retained"] == 1
    assert rec["prep_log"]["ions_and_buffer_removed"] >= 1

    # Check strict fixed-width alignment for the retained Zinc
    zn_lines = [l for l in rec["pdbqt_text"].splitlines() if "ZN" in l]
    assert len(zn_lines) == 1
    zn_line = zn_lines[0]
    assert len(zn_line) >= 78
    assert zn_line[30:38].strip() == "11.000"  # X
    assert zn_line[38:46].strip() == "11.000"  # Y
    assert zn_line[46:54].strip() == "11.000"  # Z
    assert zn_line[70:76].strip() == "2.000"   # +2.000 charge
    assert zn_line[78:80].strip() == "Zn"      # AD4 Zn atom type

def test_advanced_3d_interaction_types():
    rec_pdb = (
        "ATOM      1  N   ARG A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  ARG A   1      11.458  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  CZ  ARG A   1      13.000  10.000  10.000  1.00 20.00           C\n"
        "ATOM      4  NH1 ARG A   1      13.500  11.000  10.000  1.00 20.00           N\n"
        "ATOM      5  NH2 ARG A   1      13.500   9.000  10.000  1.00 20.00           N\n"
    )
    # Ligand carboxylate oxygen for salt bridge
    lig_pdbqt = "HETATM    1  O1  LIG     1      13.500  10.000  12.500  1.00 20.00           O\n"

    contacts = InteractionEngine.analyze(rec_pdb, lig_pdbqt)
    assert len(contacts["salt_bridges"]) == 1
    sb = contacts["salt_bridges"][0]
    assert "start_coord" in sb and len(sb["start_coord"]) == 3
    assert "end_coord" in sb and len(sb["end_coord"]) == 3
    assert sb["distance"] <= 3.5
