import pytest
from backend.services.docking import DockingEngine
from backend.services.fetcher import StructureFetcher

def test_ligand_preparation():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    lig = DockingEngine.prepare_ligand(smiles)
    assert "pdbqt_text" in lig
    assert "ROOT" in lig["pdbqt_text"]
    assert "ENDROOT" in lig["pdbqt_text"]
    assert lig["heavy_atom_count"] == 13

def test_receptor_preparation_and_docking():
    # Fetch 1CX2 and prepare chain A
    rec_meta = StructureFetcher.fetch_rcsb_pdb("1CX2")
    assert rec_meta is not None
    rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"], target_chain="A")
    assert rec["atom_count"] > 1000
    assert "pdbqt_text" in rec

    # Prepare Aspirin
    lig = DockingEngine.prepare_ligand("CC(=O)Oc1ccccc1C(=O)O")

    # Run quick docking (exhaustiveness=1, 2 modes)
    pocket = rec["detected_pocket"]
    poses = DockingEngine.run_docking(
        rec["pdbqt_text"],
        lig["pdbqt_text"],
        pocket["center"],
        pocket["size"],
        exhaustiveness=1,
        num_modes=2
    )

    assert len(poses) >= 1
    top_pose = poses[0]
    assert top_pose["affinity_kcal"] < 0.0
    assert "pdb_block" in top_pose

    # Test interaction extraction
    contacts = DockingEngine.analyze_interactions(rec["cleaned_pdb"], top_pose["pdbqt_content"])
    assert "hydrogen_bonds" in contacts
    assert "hydrophobic_contacts" in contacts
