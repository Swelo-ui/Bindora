import pytest
from backend.services.docking import DockingEngine
from backend.services.fetcher import StructureFetcher
from backend.services.pocket_detection import PocketDetectionService
from backend.services.ensemble import EnsembleDockingService
from backend.services.pharmacophore import PharmacophoreService

def test_similarity_search_aspirin():
    aspirin_smiles = 'CC(=O)Oc1ccccc1C(=O)O'
    results = StructureFetcher.search_similar_compounds(aspirin_smiles, threshold=80, max_records=4)
    assert isinstance(results, list)
    assert len(results) > 0
    first = results[0]
    assert 'cid' in first
    assert 'smiles' in first
    assert 'title' in first or 'name' in first

def test_pocket_detection_without_ligand():
    meta = StructureFetcher.fetch_rcsb_pdb('1CX2')
    assert meta is not None and 'pdb_content' in meta
    pdb_content = meta['pdb_content']
    pockets = PocketDetectionService.detect_pockets(pdb_content, max_pockets=3)
    assert len(pockets) > 0
    top = pockets[0]
    assert 'center' in top
    assert 'size' in top
    assert 'druggability_score' in top
    assert top['druggability_score'] > 0.0

def test_ensemble_structures_fetch():
    structures = EnsembleDockingService.fetch_ensemble_structures('P00533', limit=4)
    assert len(structures) > 0
    for s in structures:
        assert 'pdb_id' in s
        assert 'method' in s

def test_pharmacophore_feature_extraction_and_consensus():
    feats = PharmacophoreService.extract_features('CC(=O)Oc1ccccc1C(=O)O')
    assert isinstance(feats, dict)
    assert feats.get('Aromatic', 0) >= 1
    assert feats.get('Acceptor', 0) >= 2
    actives = [
        {'molecule_chembl_id': 'M1', 'canonical_smiles': 'CC(=O)Oc1ccccc1C(=O)O'},
        {'molecule_chembl_id': 'M2', 'canonical_smiles': 'Oc1ccccc1C(=O)O'},
        {'molecule_chembl_id': 'M3', 'canonical_smiles': 'Cc1ccc(cc1)C(=O)O'},
    ]
    profile = PharmacophoreService.build_consensus_profile(actives)
    assert 'core_requirements' in profile
    match = PharmacophoreService.match_candidate('CC(=O)Oc1ccccc1C(=O)O', profile)
    assert match['match_score_pct'] > 50.0

def test_docking_with_flexible_residues():
    meta = StructureFetcher.fetch_rcsb_pdb('1CX2')
    assert meta is not None and 'pdb_content' in meta
    pdb_content = meta['pdb_content']
    rec_prep = DockingEngine.prepare_receptor(pdb_content, target_chain='A')
    lig_prep = DockingEngine.prepare_ligand('c1ccccc1')
    center = rec_prep['detected_pocket']['center']
    size = {'x': 18.0, 'y': 18.0, 'z': 18.0}
    poses = DockingEngine.run_docking(
        receptor_pdbqt=rec_prep['pdbqt_text'],
        ligand_pdbqt=lig_prep['pdbqt_text'],
        center=center,
        size=size,
        exhaustiveness=1,
        num_modes=1,
        receptor_pdb=rec_prep['cleaned_pdb'],
        flexible_residues=['A:523']
    )
    assert len(poses) > 0
    top = poses[0]
    assert 'affinity_kcal' in top
    assert top['affinity_kcal'] < 0
