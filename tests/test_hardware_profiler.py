import pytest
from backend.services.hardware_profiler import HardwareProfiler

def test_system_specs_detection():
    specs = HardwareProfiler.get_system_specs()
    assert "cpu_count" in specs
    assert specs["cpu_count"] >= 1
    assert "os" in specs
    assert "machine" in specs
    assert specs["total_ram_gb"] is not None
    assert specs["total_ram_gb"] > 0

def test_hardware_profile_tier():
    profile = HardwareProfiler.get_hardware_profile()
    assert "tier_key" in profile
    assert profile["tier_key"] in ("TIER_1_LOW", "TIER_2_BALANCED", "TIER_3_WORKSTATION", "TIER_4_HPC")
    assert "tier_name" in profile
    assert "badge" in profile
    assert "icon" in profile
    assert "recommended_threads" in profile
    assert profile["recommended_threads"] >= 1
    assert profile["default_exhaustiveness"] >= 4

def test_calculate_adaptive_exhaustiveness():
    # 1. User manual override
    exh, reason = HardwareProfiler.calculate_adaptive_exhaustiveness(requested_exhaustiveness=16)
    assert exh == 16
    assert "Manual" in reason

    # 2. Smart auto mode on small molecule
    exh_small, reason_small = HardwareProfiler.calculate_adaptive_exhaustiveness(
        requested_exhaustiveness="auto",
        rotatable_bonds=1,
        heavy_atoms=9
    )
    assert exh_small >= 4

    # 3. Smart auto mode on flexible molecule
    exh_flex, reason_flex = HardwareProfiler.calculate_adaptive_exhaustiveness(
        requested_exhaustiveness="auto",
        rotatable_bonds=12,
        heavy_atoms=45
    )
    assert exh_flex >= exh_small

    # 4. Eco mode
    exh_eco, reason_eco = HardwareProfiler.calculate_adaptive_exhaustiveness(
        requested_exhaustiveness="auto",
        power_mode="eco"
    )
    assert exh_eco < exh_flex or exh_eco <= 6

def test_get_optimal_threads():
    threads = HardwareProfiler.get_optimal_threads()
    assert threads >= 1

    eco_threads = HardwareProfiler.get_optimal_threads(power_mode="eco")
    assert eco_threads >= 1
