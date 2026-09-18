import os
from pathlib import Path

import sys
import shutil

BASE_DIR = Path(__file__).resolve().parent.parent
BACKEND_DIR = BASE_DIR / "backend"
FRONTEND_DIR = BASE_DIR / "frontend"
BIN_DIR = BASE_DIR / "bin"

def _resolve_vina_path():
    env_vina = os.environ.get("VINA_EXE")
    if env_vina and os.path.exists(env_vina):
        return Path(env_vina)
    if sys.platform == "win32":
        win_bin = BIN_DIR / "vina.exe"
        if win_bin.exists():
            return win_bin
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

DATA_DIR = BASE_DIR / "data"
CACHE_DIR = DATA_DIR / "cache"
BENCHMARKS_DIR = DATA_DIR / "benchmarks"

# Ensure runtime directories exist
CACHE_DIR.mkdir(parents=True, exist_ok=True)
BENCHMARKS_DIR.mkdir(parents=True, exist_ok=True)
BIN_DIR.mkdir(parents=True, exist_ok=True)

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
HOST = os.environ.get("BINDORA_HOST", os.environ.get("ANUDOCK_HOST", "127.0.0.1"))
PORT = int(os.environ.get("BINDORA_PORT", os.environ.get("ANUDOCK_PORT", "5000")))
DEBUG = os.environ.get("BINDORA_DEBUG", os.environ.get("ANUDOCK_DEBUG", "False")).lower() in ("true", "1", "yes")

# Security configuration
MAX_CONTENT_LENGTH = int(os.environ.get("BINDORA_MAX_CONTENT_LENGTH", str(32 * 1024 * 1024)))  # 32 MB default
CORS_ORIGINS = os.environ.get("BINDORA_CORS_ORIGINS", "http://localhost:5000,http://127.0.0.1:5000")

# Database configuration
DATABASE_URL = os.environ.get("DATABASE_URL", f"sqlite:///{DATA_DIR / 'bindora.db'}")

# Public APIs
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
RCSB_FILE_URL = "https://files.rcsb.org/download"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
