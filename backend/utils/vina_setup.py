import os
import sys
import urllib.request
import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent
BIN_DIR = BASE_DIR / "bin"
VINA_EXE = BIN_DIR / "vina.exe"
VINA_URL = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe"

def ensure_vina():
    """Ensure the AutoDock Vina binary exists and is functional."""
    BIN_DIR.mkdir(parents=True, exist_ok=True)
    
    if VINA_EXE.exists() and VINA_EXE.stat().st_size > 100000:
        print(f"[VINA] Found existing binary at: {VINA_EXE}")
        return str(VINA_EXE)
        
    print(f"[VINA] Downloading AutoDock Vina 1.2.7 for Windows from:\n  {VINA_URL}")
    try:
        req = urllib.request.Request(VINA_URL, headers={"User-Agent": "Bindora/1.0"})
        with urllib.request.urlopen(req, timeout=60) as response, open(VINA_EXE, "wb") as out_file:
            data = response.read()
            out_file.write(data)
        print(f"[VINA] Download complete! Size: {VINA_EXE.stat().st_size / 1024 / 1024:.2f} MB")
    except Exception as e:
        print(f"[VINA ERROR] Failed to download Vina: {e}", file=sys.stderr)
        raise

    # Verify execution
    try:
        res = subprocess.run([str(VINA_EXE), "--help"], capture_output=True, text=True, timeout=10)
        if "AutoDock Vina" in res.stdout or "AutoDock Vina" in res.stderr:
            print("[VINA] Binary verified successfully!")
            return str(VINA_EXE)
        else:
            print("[VINA WARNING] Binary ran but output did not contain expected banner.")
    except Exception as e:
        print(f"[VINA WARNING] Verification run failed: {e}", file=sys.stderr)
        
    return str(VINA_EXE)

if __name__ == "__main__":
    path = ensure_vina()
    print("Vina ready at:", path)
