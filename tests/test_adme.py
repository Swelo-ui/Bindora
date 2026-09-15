import pytest
from backend.services.adme import ADMEProfiler

def test_adme_aspirin():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    result = ADMEProfiler.calculate_adme(smiles)
    
    assert "error" not in result
    phys = result["physicochemical"]
    assert abs(phys["molecular_weight"]["value"] - 180.16) < 0.5
    assert phys["hbd"]["value"] == 1
    assert phys["hba"]["value"] == 3
    
    lipinski = result["drug_likeness"]["lipinski"]
    assert lipinski["status"] == "Pass"
    assert lipinski["violations_count"] == 0

    pk = result["pharmacokinetics"]
    assert pk["gi_absorption"]["level"] == "High"

    safety = result["medicinal_chemistry_safety"]
    assert safety["pains_alerts"]["count"] == 0

def test_adme_invalid_smiles():
    result = ADMEProfiler.calculate_adme("INVALID_NOT_A_SMILES_STRING")
    assert "error" in result
