import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
BACKEND_DIR = BASE_DIR / "backend"
FRONTEND_DIR = BASE_DIR / "frontend"
BIN_DIR = BASE_DIR / "bin"
VINA_EXE = BIN_DIR / "vina.exe"
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
DEBUG = os.environ.get("BINDORA_DEBUG", os.environ.get("ANUDOCK_DEBUG", "True")).lower() in ("true", "1", "yes")

# Public APIs
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
RCSB_FILE_URL = "https://files.rcsb.org/download"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
