import pytest
import json
from backend.app import app

@pytest.fixture
def client():
    app.config["TESTING"] = True
    with app.test_client() as client:
        yield client

def test_api_health(client):
    res = client.get("/api/health")
    assert res.status_code == 200
    data = res.get_json()
    assert data["status"] == "healthy"
    assert data["vina_available"] is True

def test_api_benchmarks(client):
    res = client.get("/api/benchmarks")
    assert res.status_code == 200
    data = res.get_json()
    assert "benchmarks" in data
    assert len(data["benchmarks"]) >= 4

def test_api_pubchem_search(client):
    res = client.get("/api/search/pubchem?query=aspirin")
    assert res.status_code == 200
    data = res.get_json()
    assert data["cid"] == 2244
    assert "adme" in data

def test_api_adme_calculation(client):
    res = client.post("/api/pkpd/adme", json={"smiles": "CC(=O)Oc1ccccc1C(=O)O"})
    assert res.status_code == 200
    data = res.get_json()
    assert data["drug_likeness"]["lipinski"]["status"] == "Pass"

def test_api_reproducibility_versions(client):
    res = client.get("/api/reproducibility/versions")
    assert res.status_code == 200
    data = res.get_json()
    assert "autodock_vina" in data
    assert "rdkit" in data
    assert "gemmi" in data
    assert "utc_timestamp" in data

def test_api_batch_sanitization_and_consensus(client):
    from backend.services.fetcher import StructureFetcher
    from backend.services.docking import DockingEngine

    rec_meta = StructureFetcher.fetch_rcsb_pdb("1CX2")
    rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"], target_chain="A")

    payload = {
        "receptor_pdbqt": rec["pdbqt_text"],
        "receptor_pdb": rec["cleaned_pdb"],
        "center": rec["detected_pocket"]["center"],
        "size": {"x": 20.0, "y": 20.0, "z": 20.0},
        "ligands": [
            {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
            {"name": "MalformedDrug", "smiles": "INVALID((SMILES$$$"},
        ],
        "exhaustiveness": 1
    }

    res = client.post("/api/batch/run", json=payload)
    assert res.status_code == 200
    data = res.get_json()
    leaderboard = data.get("leaderboard", [])
    assert len(leaderboard) == 2

    # Aspirin should be valid with consensus score
    valid_row = [row for row in leaderboard if row["name"] == "Aspirin"][0]
    assert valid_row["valid"] is True
    assert valid_row["affinity_kcal"] < 0.0
    assert valid_row["consensus_score"] > 0.0

    # MalformedDrug should report explicit error rather than crashing or silent skip
    failed_row = [row for row in leaderboard if row["name"] == "MalformedDrug"][0]
    assert failed_row["valid"] is False
    assert "Failed" in failed_row["status"]
    assert failed_row["rank"] == "—"

