import os
import sys
import urllib.request
import subprocess
from pathlib import Path

import shutil

BASE_DIR = Path(__file__).resolve().parent.parent.parent
BIN_DIR = BASE_DIR / "bin"

_VINA_NOTIFIED = False

def get_platform_vina_info():
    if sys.platform == "win32":
        return BIN_DIR / "vina.exe", "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe"
    elif sys.platform == "darwin":
        return BIN_DIR / "vina", "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_x86_64"
    else:
        return BIN_DIR / "vina", "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64"

def ensure_vina(verbose: bool = False):
    """Ensure the AutoDock Vina binary exists and is functional across Windows, Linux, and macOS."""
    global _VINA_NOTIFIED
    BIN_DIR.mkdir(parents=True, exist_ok=True)
    
    # Check system PATH (common on Linux/Colab/Docker)
    sys_vina = shutil.which("vina")
    if sys_vina:
        if verbose and not _VINA_NOTIFIED:
            print(f"[VINA] Found system Vina on PATH: {sys_vina}")
            _VINA_NOTIFIED = True
        return str(sys_vina)

    vina_exe, vina_url = get_platform_vina_info()

    if vina_exe.exists() and vina_exe.stat().st_size > 100000:
        if verbose and not _VINA_NOTIFIED:
            print(f"[VINA] Found existing binary at: {vina_exe}")
            _VINA_NOTIFIED = True
        return str(vina_exe)
        
    print(f"[VINA] Downloading AutoDock Vina 1.2.7 ({sys.platform}) from:\n  {vina_url}")
    try:
        req = urllib.request.Request(vina_url, headers={"User-Agent": "Bindora/1.0"})
        with urllib.request.urlopen(req, timeout=60) as response, open(vina_exe, "wb") as out_file:
            data = response.read()
            out_file.write(data)
        if sys.platform != "win32":
            os.chmod(vina_exe, 0o755)
        print(f"[VINA] Download complete! Size: {vina_exe.stat().st_size / 1024 / 1024:.2f} MB")
    except Exception as e:
        print(f"[VINA ERROR] Failed to download Vina: {e}", file=sys.stderr)
        raise

    # Verify execution
    try:
        res = subprocess.run([str(vina_exe), "--help"], capture_output=True, text=True, timeout=10)
        if "AutoDock Vina" in res.stdout or "AutoDock Vina" in res.stderr:
            print("[VINA] Binary verified successfully!")
            return str(vina_exe)
        else:
            print("[VINA WARNING] Binary ran but output did not contain expected banner.")
    except Exception as e:
        print(f"[VINA WARNING] Verification run failed: {e}", file=sys.stderr)
        
    return str(vina_exe)

if __name__ == "__main__":
    path = ensure_vina()
    print("Vina ready at:", path)
