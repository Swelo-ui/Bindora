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
