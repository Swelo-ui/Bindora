import os
import sys
import platform
import multiprocessing
from typing import Dict, Any, Optional, Tuple

class HardwareProfiler:
    """
    Intelligent Hardware Profiler for Bindora Dock.
    Detects CPU cores, system RAM, platform architecture, and dynamically
    allocates optimal threads, exhaustiveness, and batch workers based on
    the host hardware profile (Low-end laptop, standard, workstation, or cloud server).
    """

    @staticmethod
    def get_system_specs() -> Dict[str, Any]:
        """Detect CPU cores, total RAM, available RAM, and platform architecture."""
        cpu_count = os.cpu_count() or multiprocessing.cpu_count() or 1
        total_ram_gb = None
        avail_ram_gb = None

        if sys.platform == "win32":
            try:
                import ctypes
                class MEMORYSTATUSEX(ctypes.Structure):
                    _fields_ = [
                        ("dwLength", ctypes.c_ulong),
                        ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong),
                        ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                    ]
                stat = MEMORYSTATUSEX()
                stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
                ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
                total_ram_gb = round(stat.ullTotalPhys / (1024 ** 3), 2)
                avail_ram_gb = round(stat.ullAvailPhys / (1024 ** 3), 2)
            except Exception:
                pass
        elif sys.platform.startswith("linux"):
            try:
                with open("/proc/meminfo", "r", encoding="utf-8") as f:
                    mem_data = {}
                    for line in f:
                        parts = line.split(":")
                        if len(parts) == 2:
                            mem_data[parts[0].strip()] = parts[1].strip()
                    if "MemTotal" in mem_data:
                        kb = int(mem_data["MemTotal"].split()[0])
                        total_ram_gb = round(kb / (1024 ** 2), 2)
                    if "MemAvailable" in mem_data:
                        kb = int(mem_data["MemAvailable"].split()[0])
                        avail_ram_gb = round(kb / (1024 ** 2), 2)
            except Exception:
                pass
        elif sys.platform == "darwin":
            try:
                import subprocess
                out = subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True).strip()
                total_ram_gb = round(int(out) / (1024 ** 3), 2)
            except Exception:
                pass

        if total_ram_gb is None:
            total_ram_gb = 4.0 if cpu_count <= 2 else (8.0 if cpu_count <= 4 else 16.0)

        is_cloud = (
            bool(os.environ.get("CODESPACES")) or
            bool(os.environ.get("GITHUB_CODESPACES_PORT_FORWARDING_DOMAIN")) or
            bool(os.environ.get("COLAB_GPU")) or
            bool(os.environ.get("KUBERNETES_SERVICE_HOST")) or
            "github.dev" in os.environ.get("HOSTNAME", "")
        )

        return {
            "cpu_count": cpu_count,
            "total_ram_gb": total_ram_gb,
            "avail_ram_gb": avail_ram_gb,
            "os": platform.system(),
            "os_release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor() or platform.machine(),
            "is_cloud": is_cloud
        }

    @classmethod
    def get_hardware_profile(cls) -> Dict[str, Any]:
        """
        Classifies current host machine into an optimized computational tier
        and provides recommended execution parameters.
        """
        specs = cls.get_system_specs()
        cpu = specs["cpu_count"]
        ram = specs["total_ram_gb"] or 8.0
        is_cloud = specs["is_cloud"]

        if cpu <= 2 or ram <= 4.0:
            tier_key = "TIER_1_LOW"
            tier_name = "Constrained / Ultra-Portable Laptop"
            icon = "🌱"
            badge = f"{cpu} Cores (Low-Overhead Mode)"
            default_exhaustiveness = 4
            recommended_threads = cpu
            batch_workers = 1
            desc = "Lightweight profile: Keeps memory footprint low and prevents system lag on budget hardware."
        elif cpu <= 6 and ram <= 12.0:
            tier_key = "TIER_2_BALANCED"
            tier_name = "Standard Multi-Core Laptop" if not is_cloud else "Cloud Compute (Azure/Codespace)"
            icon = "⚡"
            badge = f"{cpu} Cores (Balanced Parallel)"
            default_exhaustiveness = 8
            recommended_threads = cpu
            batch_workers = 2
            desc = "Balanced profile: 100% multi-core acceleration with standard publication-grade exhaustiveness."
        elif cpu <= 16:
            tier_key = "TIER_3_WORKSTATION"
            tier_name = "High-Performance Research Workstation"
            icon = "🚀"
            badge = f"{cpu} Cores (High-Power Workstation)"
            default_exhaustiveness = 12
            recommended_threads = cpu - 1 if cpu >= 8 else cpu
            batch_workers = max(2, cpu // 3)
            desc = "Workstation profile: Deep conformational sampling with parallel multi-threading active."
        else:
            tier_key = "TIER_4_HPC"
            tier_name = "Enterprise HPC Server"
            icon = "🏎️"
            badge = f"{cpu} Cores (Enterprise HPC)"
            default_exhaustiveness = 16
            recommended_threads = cpu
            batch_workers = max(4, cpu // 4)
            desc = "Server-class profile: Ultra-deep Monte Carlo search with massive parallel throughput."

        return {
            "tier_key": tier_key,
            "tier_name": tier_name,
            "icon": icon,
            "badge": badge,
            "specs": specs,
            "recommended_threads": recommended_threads,
            "default_exhaustiveness": default_exhaustiveness,
            "batch_workers": batch_workers,
            "description": desc
        }

    @classmethod
    def calculate_adaptive_exhaustiveness(
        cls,
        requested_exhaustiveness: Optional[Any] = None,
        rotatable_bonds: int = 0,
        heavy_atoms: int = 0,
        power_mode: str = "smart"
    ) -> Tuple[int, str]:
        profile = cls.get_hardware_profile()
        tier_key = profile["tier_key"]
        tier_default = profile["default_exhaustiveness"]

        if requested_exhaustiveness is not None and str(requested_exhaustiveness).lower() not in ("auto", "0", ""):
            try:
                val = int(requested_exhaustiveness)
                if val > 0:
                    return max(1, min(64, val)), f"Manual override ({val})"
            except (ValueError, TypeError):
                pass

        mode = (power_mode or "smart").lower()
        if mode == "eco":
            return max(3, tier_default - 2), f"Eco/Battery saver mode ({tier_default - 2})"
        elif mode == "performance":
            return min(32, tier_default + 4), f"Maximum performance mode ({tier_default + 4})"

        base = tier_default
        cpu_n = profile["specs"]["cpu_count"]

        if rotatable_bonds <= 3 and heavy_atoms < 20:
            adaptive = max(4, base - 2)
            return adaptive, f"Smart Adaptive ({cpu_n} Cores): Rapid convergence for small/rigid molecule (exh={adaptive})"
        elif rotatable_bonds >= 10:
            if tier_key in ("TIER_3_WORKSTATION", "TIER_4_HPC"):
                adaptive = min(24, base + 4)
                return adaptive, f"Smart Adaptive ({cpu_n} Cores): High-power deep sampling for {rotatable_bonds} torsions (exh={adaptive})"
            else:
                return base, f"Smart Adaptive ({cpu_n} Cores): Balanced sampling for {rotatable_bonds} torsions (exh={base})"

        return base, f"Smart Adaptive ({cpu_n} Cores): Standard research-grade convergence (exh={base})"

    @classmethod
    def get_optimal_threads(cls, requested_cpu: Optional[Any] = None, power_mode: str = "smart") -> int:
        specs = cls.get_system_specs()
        max_cores = specs["cpu_count"]

        if requested_cpu is not None:
            try:
                c = int(requested_cpu)
                if 1 <= c <= max_cores:
                    return c
            except (ValueError, TypeError):
                pass

        mode = (power_mode or "smart").lower()
        if mode == "eco":
            return max(1, max_cores // 2)
        
        profile = cls.get_hardware_profile()
        return profile["recommended_threads"]
