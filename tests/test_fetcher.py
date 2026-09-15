import pytest
from backend.services.fetcher import StructureFetcher

def test_pubchem_search_aspirin():
    result = StructureFetcher.search_pubchem("aspirin")
    assert result is not None
    assert result["cid"] == 2244
    assert result["formula"] == "C9H8O4"
    assert "C(=O)O" in result["smiles"]
    assert result["weight"] > 170.0

def test_rcsb_pdb_fetch_1cx2():
    result = StructureFetcher.fetch_rcsb_pdb("1CX2")
    assert result is not None
    assert result["pdb_id"] == "1CX2"
    assert "CYCLOOXYGENASE" in result["title"].upper()
    assert "ATOM" in result["pdb_content"]

def test_uniprot_search():
    result = StructureFetcher.search_uniprot("BCR-ABL1")
    assert result is not None
    assert result["accession"] != ""
