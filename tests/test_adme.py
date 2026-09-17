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

def test_adme_sascore_and_expanded_alerts():
    # Aspirin
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    res = ADMEProfiler.calculate_adme(aspirin)
    assert "sascore" in res["physicochemical"]
    sas = res["physicochemical"]["sascore"]
    assert 1.0 <= sas["value"] <= 3.0
    assert sas["interpretation"] == "Easy synthetic accessibility"
    
    # Check 4 alert catalogs
    safety = res["medicinal_chemistry_safety"]
    assert "nih_alerts" in safety
    assert "zinc_alerts" in safety
    assert "total_alerts_count" in safety
    assert safety["total_alerts_count"] >= 1 # Brenk phenol ester

def test_boiled_egg_published_ellipse():
    # Atenolol: well-known non-BBB permeant, high GI absorption
    atenolol_smi = "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"
    res_at = ADMEProfiler.calculate_adme(atenolol_smi)
    assert res_at["pharmacokinetics"]["gi_absorption"]["level"] == "High"
    assert "Non-permeant" in res_at["pharmacokinetics"]["bbb_permeation"]["status"]

    # Diazepam: well-known BBB permeant and high GI absorption
    diazepam_smi = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"
    res_dz = ADMEProfiler.calculate_adme(diazepam_smi)
    assert res_dz["pharmacokinetics"]["gi_absorption"]["level"] == "High"
    assert "Permeant" in res_dz["pharmacokinetics"]["bbb_permeation"]["status"]

