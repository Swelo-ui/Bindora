"""
Bindora Hardware Profiler & Continuous Multi-Core Telemetry Sampler.
Provides smooth, real-time CPU and Memory diagnostics matching native host task manager metrics.
"""

import sys
import time
import threading
import platform
from typing import Dict, Any, List
import psutil

class HardwareTelemetrySampler:
    """
    Background daemon sampler that continuously monitors per-core and overall CPU load,
    as well as system RAM, eliminating 0%/100% burst jitter from instantaneous sampling.
    """
    _instance = None
    _init_lock = threading.Lock()

    def __init__(self):
        self.lock = threading.Lock()
        
        # 1. Host CPU Detection
        self.cpu_name = self._detect_cpu_name()
        self.platform = sys.platform
        self.architecture = platform.machine()
        
        self.logical_cores = psutil.cpu_count(logical=True) or 4
        self.physical_cores = psutil.cpu_count(logical=False) or self.logical_cores

        # 2. Telemetry State Buffers
        self.latest_per_core: List[float] = [12.0] * self.logical_cores
        self.latest_overall: float = 12.0
        self.total_ram_gb: float = 8.0
        self.used_ram_gb: float = 4.0
        self.available_ram_gb: float = 4.0
        self.ram_percent: float = 50.0

        # Prime initial values
        self._sample_once(is_priming=True)

        # 3. Launch background daemon worker
        self._running = True
        self.worker_thread = threading.Thread(
            target=self._telemetry_loop,
            daemon=True,
            name="BindoraHardwareTelemetryWorker"
        )
        self.worker_thread.start()

    def _detect_cpu_name(self) -> str:
        """Query native registry on Windows, /proc/cpuinfo on Linux, or sysctl on macOS for exact CPU string."""
        if sys.platform == "win32":
            try:
                import winreg
                key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
                val, _ = winreg.QueryValueEx(key, "ProcessorNameString")
                winreg.CloseKey(key)
                if val and str(val).strip():
                    return str(val).strip()
            except Exception:
                pass
        elif sys.platform.startswith("linux"):
            try:
                with open("/proc/cpuinfo", "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        if "model name" in line:
                            return line.split(":", 1)[1].strip()
            except Exception:
                pass
        elif sys.platform == "darwin":
            try:
                import subprocess
                val = subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string"], text=True).strip()
                if val:
                    return val
            except Exception:
                pass
        
        proc = platform.processor()
        return proc if proc else "Multi-Core System Processor"

    def _sample_once(self, is_priming: bool = False):
        try:
            # When interval is provided, psutil blocks for exactly that duration and measures delta ticks
            sample_interval = 0.1 if is_priming else 0.5
            raw_cores = psutil.cpu_percent(interval=sample_interval, percpu=True)
            vmem = psutil.virtual_memory()

            if raw_cores and len(raw_cores) == self.logical_cores:
                with self.lock:
                    if is_priming:
                        self.latest_per_core = [round(c, 1) for c in raw_cores]
                    else:
                        # Exponential moving average filter (alpha = 0.65) to smooth out sampling micro-jitter
                        smoothed = []
                        for i, curr_val in enumerate(raw_cores):
                            prev_val = self.latest_per_core[i] if i < len(self.latest_per_core) else curr_val
                            smooth_val = round((prev_val * 0.35) + (curr_val * 0.65), 1)
                            smoothed.append(smooth_val)
                        self.latest_per_core = smoothed

                    self.latest_overall = round(sum(self.latest_per_core) / len(self.latest_per_core), 1)
                    self.total_ram_gb = round(vmem.total / (1024 ** 3), 2)
                    self.used_ram_gb = round(vmem.used / (1024 ** 3), 2)
                    self.available_ram_gb = round(vmem.available / (1024 ** 3), 2)
                    self.ram_percent = round(vmem.percent, 1)
        except Exception:
            pass

    def _telemetry_loop(self):
        """Continuous background thread loop sampling every 500ms."""
        # Warm-up delay
        time.sleep(0.2)
        while self._running:
            try:
                self._sample_once(is_priming=False)
            except Exception:
                time.sleep(0.5)

    def get_telemetry(self, is_docking_active: bool = False) -> Dict[str, Any]:
        """Return non-blocking real-time hardware telemetry snapshot in microseconds."""
        with self.lock:
            per_core = list(self.latest_per_core)
            overall = self.latest_overall
            total_ram = self.total_ram_gb
            used_ram = self.used_ram_gb
            avail_ram = self.available_ram_gb
            ram_pct = self.ram_percent

        return {
            "cpu_name": self.cpu_name,
            "platform": self.platform,
            "architecture": self.architecture,
            "physical_cores": self.physical_cores,
            "logical_cores": self.logical_cores,
            "vina_allocated_threads": self.logical_cores,
            "is_docking_active": bool(is_docking_active),
            "docking_cores_in_use": self.logical_cores if is_docking_active else 0,
            "overall_cpu_percent": overall,
            "per_core_percent": per_core,
            "total_ram_gb": total_ram,
            "available_ram_gb": avail_ram,
            "used_ram_gb": used_ram,
            "ram_percent": ram_pct,
            "gpu": self.detect_gpu_capabilities()
        }

    @classmethod
    def detect_gpu_capabilities(cls) -> Dict[str, Any]:
        """Detect available GPU acceleration platforms (CUDA, OpenCL) and GPU docking engine binaries."""
        import shutil
        from pathlib import Path
        from backend.config import VINA_GPU_EXE, UNIDOCK_EXE, BIN_DIR

        cuda_available = False
        opencl_available = False
        platforms_found = []

        try:
            import openmm as mm
            num_p = mm.Platform.getNumPlatforms()
            for i in range(num_p):
                p_name = mm.Platform.getPlatform(i).getName()
                platforms_found.append(p_name)
                if p_name.upper() == "CUDA":
                    cuda_available = True
                elif p_name.upper() == "OPENCL":
                    opencl_available = True
        except Exception:
            pass

        # Check GPU docking binaries
        gpu_bin_found = None
        if VINA_GPU_EXE and Path(VINA_GPU_EXE).exists():
            gpu_bin_found = str(VINA_GPU_EXE)
        elif UNIDOCK_EXE and Path(UNIDOCK_EXE).exists():
            gpu_bin_found = str(UNIDOCK_EXE)
        else:
            for bname in ["vina-gpu", "AutoDock-Vina-GPU-2.1", "unidock"]:
                w = shutil.which(bname)
                if w:
                    gpu_bin_found = w
                    break

        preferred_accel = "CPU"
        if cuda_available:
            preferred_accel = "CUDA (NVIDIA Tensor/CUDA Cores)"
        elif opencl_available:
            preferred_accel = "OpenCL (Hardware GPU Acceleration)"

        return {
            "cuda_available": cuda_available,
            "opencl_available": opencl_available,
            "gpu_hardware_ready": cuda_available or opencl_available,
            "preferred_acceleration_platform": preferred_accel,
            "openmm_platforms": platforms_found,
            "gpu_docking_binary_present": bool(gpu_bin_found),
            "gpu_docking_binary_path": gpu_bin_found,
            "gpu_docking_ready": bool(gpu_bin_found and (cuda_available or opencl_available))
        }

    @classmethod
    def get_instance(cls) -> "HardwareTelemetrySampler":
        if cls._instance is None:
            with cls._init_lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance

