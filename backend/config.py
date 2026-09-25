import os
from pathlib import Path

import sys
import shutil
import subprocess

# Ensure stdout and stderr are non-null in windowed (noconsole) mode
if sys.stdout is None or sys.stderr is None:
    appdata = os.environ.get("APPDATA")
    log_dir = Path(appdata) / "Bindora" if appdata else Path.home() / ".bindora"
    try:
        log_dir.mkdir(parents=True, exist_ok=True)
        backend_log = open(log_dir / "backend.log", "a", encoding="utf-8", buffering=1)
        if sys.stdout is None:
            sys.stdout = backend_log
        if sys.stderr is None:
            sys.stderr = backend_log
    except Exception:
        pass

def get_subprocess_kwargs():
    """
    Return kwargs to suppress console/terminal windows when spawning
    subprocesses (e.g. Vina, fpocket, WSL) on Windows desktop GUI apps.
    """
    kwargs = {}
    if sys.platform == "win32":
        kwargs["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)
        si = subprocess.STARTUPINFO()
        si.dwFlags |= getattr(subprocess, "STARTF_USESHOWWINDOW", 1)
        si.wShowWindow = getattr(subprocess, "SW_HIDE", 0)
        kwargs["startupinfo"] = si
    return kwargs

IS_FROZEN = getattr(sys, "frozen", False)
if IS_FROZEN:
    # PyInstaller unpacks data or runs from onedir folder (_internal)
    bundle_candidate = Path(getattr(sys, "_MEIPASS", Path(sys.executable).resolve().parent))
    if not (bundle_candidate / "frontend").exists() and (bundle_candidate / "_internal" / "frontend").exists():
        BASE_DIR = bundle_candidate / "_internal"
    elif not (bundle_candidate / "frontend").exists() and (bundle_candidate.parent / "_internal" / "frontend").exists():
        BASE_DIR = bundle_candidate.parent / "_internal"
    else:
        BASE_DIR = bundle_candidate

    BACKEND_DIR = BASE_DIR / "backend"
    FRONTEND_DIR = BASE_DIR / "frontend"
    BIN_DIR = BASE_DIR / "bin"
    if not BIN_DIR.exists() and (BASE_DIR.parent / "_internal" / "bin").exists():
        BIN_DIR = BASE_DIR.parent / "_internal" / "bin"
    
    # Store persistent writable user data in %APPDATA%/Bindora (avoiding read-only Program Files)
    appdata = os.environ.get("APPDATA")
    if appdata:
        DATA_DIR = Path(appdata) / "Bindora"
    else:
        DATA_DIR = Path.home() / ".bindora"
else:
    BASE_DIR = Path(__file__).resolve().parent.parent
    BACKEND_DIR = BASE_DIR / "backend"
    FRONTEND_DIR = BASE_DIR / "frontend"
    BIN_DIR = BASE_DIR / "bin"
    DATA_DIR = BASE_DIR / "data"

def _resolve_vina_path():
    env_vina = os.environ.get("VINA_EXE")
    if env_vina and os.path.exists(env_vina):
        return Path(env_vina)
    if sys.platform == "win32":
        candidates = [
            BIN_DIR / "vina.exe",
            BASE_DIR / "bin" / "vina.exe",
            BASE_DIR / "_internal" / "bin" / "vina.exe",
            DATA_DIR / "bin" / "vina.exe",
        ]
        for c in candidates:
            if c.exists():
                return c
    sys_vina = shutil.which("vina")
    if sys_vina:
        return Path(sys_vina)
    if sys.platform == "win32":
        return BIN_DIR / "vina.exe"
    local_vina = BIN_DIR / "vina"
    if local_vina.exists():
        return local_vina
    return BIN_DIR / "vina"

VINA_EXE = _resolve_vina_path()
GNINA_EXE = os.environ.get("GNINA_EXE", str(BIN_DIR / "gnina.exe"))

CACHE_DIR = DATA_DIR / "cache"
BENCHMARKS_DIR = DATA_DIR / "benchmarks"

# Ensure runtime directories exist
DATA_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)
BENCHMARKS_DIR.mkdir(parents=True, exist_ok=True)
if not IS_FROZEN:
    try:
        BIN_DIR.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass

# When frozen, seed and synchronize bundled benchmarks into persistent directory
if IS_FROZEN:
    bundled_benchmarks = BASE_DIR / "data" / "benchmarks"
    if not bundled_benchmarks.exists() and (BASE_DIR.parent / "_internal" / "data" / "benchmarks").exists():
        bundled_benchmarks = BASE_DIR.parent / "_internal" / "data" / "benchmarks"
    if bundled_benchmarks.exists():
        for b_file in bundled_benchmarks.glob("*.json"):
            dest_file = BENCHMARKS_DIR / b_file.name
            try:
                if not dest_file.exists() or dest_file.stat().st_size != b_file.stat().st_size or b_file.stat().st_mtime > dest_file.stat().st_mtime:
                    shutil.copy2(b_file, dest_file)
            except Exception:
                pass

try:
    from dotenv import load_dotenv
    load_dotenv(BASE_DIR / ".env")
except ImportError:
    pass

PROJECT_NAME = "Bindora"

# API Keys and Models
OPENROUTER_API_KEY = os.environ.get("OPENROUTER_API_KEY", "")
OPENROUTER_MODEL = os.environ.get("OPENROUTER_MODEL", "deepseek/deepseek-v4-flash-0731")
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY", "")

# Server configuration
HOST = os.environ.get("BINDORA_HOST", "127.0.0.1")
PORT = int(os.environ.get("BINDORA_PORT", "5000"))
DEBUG = os.environ.get("BINDORA_DEBUG", "False").lower() in ("true", "1", "yes")

# Security configuration
MAX_CONTENT_LENGTH = int(os.environ.get("BINDORA_MAX_CONTENT_LENGTH", str(32 * 1024 * 1024)))  # 32 MB default
CORS_ORIGINS = os.environ.get("BINDORA_CORS_ORIGINS", "http://localhost:5000,http://127.0.0.1:5000")

# Database configuration
_sqlite_path = (DATA_DIR / "bindora.db").resolve().as_posix()
DATABASE_URL = os.environ.get("DATABASE_URL", f"sqlite:///{_sqlite_path}")

# Public APIs
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
RCSB_FILE_URL = "https://files.rcsb.org/download"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
