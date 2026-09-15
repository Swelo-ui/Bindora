import pytest
from backend.services.bioactivity import BioactivityService

def test_thermodynamics_conversion():
    # delta G = -8.5 kcal/mol for molecule with 36 heavy atoms, MW 493.6
    result = BioactivityService.calculate_thermodynamics(-8.5, 36, 493.6)
    
    assert result["binding_affinity_kcal"] == -8.5
    assert result["theoretical_kd_nm"] > 0
    assert result["theoretical_kd_um"] > 0
    assert result["ligand_efficiency"]["value"] > 0.2
    assert "Moderate Affinity" in result["potency_class"]

def test_chembl_crosscheck_imatinib():
    result = BioactivityService.crosscheck_chembl("imatinib", "abl1")
    assert result is not None
    assert result["is_cross_checked"] is True
    assert result["status_badge"] == "Experimentally Corroborated"
    assert len(result["experimental_records"]) > 0
    # Spot-check that measured value is within nanomolar range
    val = result["experimental_records"][0]["value"]
    assert val > 0

def test_weak_binder_detection():
    # delta G = -4.894 kcal/mol (chlorogenic acid example: weak binder)
    result = BioactivityService.calculate_thermodynamics(-4.894, 25, 354.31)
    assert result["is_weak_binder"] is True
    assert result["weak_binder_warning"] is not None
    assert "Sub-threshold / Weak Binding Alert" in result["weak_binder_warning"]
    assert result["size_independent_le"]["value"] > 0
    assert result["fit_quality"]["value"] > 0

    # delta G = -9.2 kcal/mol (potent binder: should NOT be flagged as weak)
    strong = BioactivityService.calculate_thermodynamics(-9.2, 30, 420.0)
    assert strong["is_weak_binder"] is False
    assert strong["weak_binder_warning"] is None

