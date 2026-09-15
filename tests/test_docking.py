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

def test_preparation_transparency_logs_and_blind_box():
    rec_meta = StructureFetcher.fetch_rcsb_pdb("1CX2")
    rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"], target_chain="A")

    # Verify receptor prep_log
    assert "prep_log" in rec
    assert rec["prep_log"]["waters_removed"] >= 0
    assert "Standard physiological pH 7.4" in rec["prep_log"]["protonation_state"]

    # Verify blind docking bounding box
    assert "blind_docking_box" in rec
    bbox = rec["blind_docking_box"]
    assert bbox["size"]["x"] >= 26.0
    assert bbox["size"]["y"] >= 26.0

    # Verify native ligand extraction
    assert "native_ligand" in rec
    assert rec["native_ligand"]["has_native"] is True
    assert rec["native_ligand"]["name"] == "S58"

    # Verify ligand prep_log
    lig = DockingEngine.prepare_ligand("CC(=O)Oc1ccccc1C(=O)O")
    assert "prep_log" in lig
    assert lig["prep_log"]["input_format"] == "SMILES"
    assert "MMFF94" in lig["prep_log"]["energy_minimization"]

def test_redocking_validation():
    rec_meta = StructureFetcher.fetch_rcsb_pdb("1CX2")
    rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"], target_chain="A")
    native = rec["native_ligand"]
    assert native["has_native"] is True

    result = DockingEngine.run_redocking_validation(
        rec["pdbqt_text"],
        native["pdb_block"],
        native["center"],
        {"x": 20.0, "y": 20.0, "z": 20.0},
        exhaustiveness=2
    )

    assert "affinity_kcal" in result
    assert result["affinity_kcal"] < 0.0
    assert "rmsd_angstroms" in result
    assert result["rmsd_angstroms"] >= 0.0
    assert "validation_badge" in result
    assert "docked_pdb" in result

def test_replicate_sampling():
    rec_meta = StructureFetcher.fetch_rcsb_pdb("1CX2")
    rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"], target_chain="A")
    lig = DockingEngine.prepare_ligand("CC(=O)Oc1ccccc1C(=O)O")

    poses = DockingEngine.run_docking(
        rec["pdbqt_text"],
        lig["pdbqt_text"],
        rec["detected_pocket"]["center"],
        rec["detected_pocket"]["size"],
        exhaustiveness=1,
        num_modes=2,
        replicates=3
    )

    assert len(poses) >= 1
    top_pose = poses[0]
    assert "replicate_stats" in top_pose
    rs = top_pose["replicate_stats"]
    assert rs["replicates_count"] == 3
    assert len(rs["affinities_kcal"]) == 3
    assert rs["mean_affinity_kcal"] < 0.0
    assert rs["sd_affinity_kcal"] >= 0.0

def test_native_ligand_fractional_occupancy_7bv2():
    # 7BV2 has native ligand F86 (Remdesivir metabolite) with fractional occupancy (0.50)
    rec_meta = StructureFetcher.fetch_rcsb_pdb("7BV2")
    assert rec_meta is not None
    rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"])
    native = rec["native_ligand"]
    assert native["has_native"] is True
    assert native["name"] == "F86"

    # Test dedicated prepare_native_ligand pipeline (PDB -> SDF -> AddHs -> Meeko PDBQT)
    nat_prep = DockingEngine.prepare_native_ligand(native["pdb_block"])
    assert "ROOT" in nat_prep["pdbqt_text"]
    assert "ENDROOT" in nat_prep["pdbqt_text"]
    assert "TORSDOF" in nat_prep["pdbqt_text"]
    assert nat_prep["heavy_atom_count"] == 24

    # Run redocking validation (exhaustiveness=1 for fast test)
    result = DockingEngine.run_redocking_validation(
        rec["pdbqt_text"],
        native["pdb_block"],
        rec["detected_pocket"]["center"],
        rec["detected_pocket"]["size"],
        exhaustiveness=1
    )

    assert "affinity_kcal" in result
    assert result["affinity_kcal"] < 0.0
    assert "rmsd_angstroms" in result
    assert "validation_badge" in result
    assert "docked_pdb" in result


