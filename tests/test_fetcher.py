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

def test_uniprot_pdb_xref_1cx2():
    """
    Integration test: verifies that the UniProt xref:pdb-{pdb_id} query syntax
    documented in README Section 2 (Pathway Annotations row) is working.
    COX-2 / 1CX2 is used as a stable, known-good reference.
    This test backs the README claim — if it fails, the query syntax is broken.
    """
    result = StructureFetcher.fetch_uniprot_by_pdb_id("1CX2")
    assert result is not None, (
        "UniProt xref:pdb-1CX2 query returned None — "
        "check UNIPROT_BASE_URL and query syntax in fetcher.py"
    )
    assert result.get("accession"), "UniProt accession must be non-empty for 1CX2"
    assert result.get("function"), (
        "UniProt function field is empty for 1CX2 (COX-2). "
        "The xref query works but no functional annotation was returned."
    )
    # COX-2 is Prostaglandin H2 synthase — gene name PTGS2
    accession = result["accession"]
    assert len(accession) in (6, 10), f"Unexpected UniProt accession format: {accession}"
