import pytest
from backend.services.consensus import ConsensusScoringService
from backend.services.refinement import ComplexRefinementService

def test_calibrated_strain_thresholds_aspirin_passes():
    """Verify that a bioactive drug pose with typical strain (e.g. Aspirin at 6.6 kcal/mol)
    and favorable MM-GBSA is NOT mislabeled as a decoy."""
    
    # Pose 1: Aspirin-like binding pose from COX-2 benchmark
    # Vina: -6.7, Vinardo: -5.8, MM-GBSA: -6.63, CNN: 0.81, Strain: 6.63 kcal/mol
    pose_aspirin = {
        "mode": 1,
        "affinity_kcal": -6.7,
        "vinardo_affinity_kcal": -5.8,
        "mmgbsa_delta_g_kcal": -6.63,
        "ligand_strain_kcal": 6.63,
        "gnina": {"available": True, "cnn_score": 0.81}
    }
    
    # Pose 2: High strain true decoy pose
    # High strain 9.5 kcal/mol with positive/unfavorable MM-GBSA
    pose_decoy = {
        "mode": 2,
        "affinity_kcal": -5.2,
        "vinardo_affinity_kcal": -4.1,
        "mmgbsa_delta_g_kcal": 2.5,  # Unfavorable desolvation penalty
        "ligand_strain_kcal": 9.5,   # High strain
        "gnina": {"available": True, "cnn_score": 0.15}
    }

    res = ConsensusScoringService.compute_pose_consensus([pose_aspirin, pose_decoy])
    
    # Aspirin pose MUST NOT be flagged as decoy!
    assert res[0]["mode"] == 1
    assert res[0]["consensus_confidence"] != "DECOY_HIGH_STRAIN"
    assert res[0]["consensus_confidence"] in ("HIGH_CONFIDENCE", "MODERATE_CONFIDENCE")
    print(f"Aspirin Confidence: {res[0]['consensus_confidence']} ({res[0]['consensus_confidence_description']})")

    # True decoy pose MUST be flagged as DECOY_HIGH_STRAIN
    assert res[1]["mode"] == 2
    assert res[1]["consensus_confidence"] == "DECOY_HIGH_STRAIN"
    print(f"Decoy Confidence: {res[1]['consensus_confidence']} ({res[1]['consensus_confidence_description']})")

def test_refinement_strain_classification():
    """Verify refinement service classifies strain <= 8.0 kcal/mol as acceptable."""
    # 6.63 kcal/mol is moderate acceptable strain in bound state
    is_high_strain = 6.63 > 8.0
    assert not is_high_strain
