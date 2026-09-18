# Design Document: Production Hardening & Architectural Upgrades

## High-Level Architecture

### System Overview

```
┌─────────────────────────────────────────────────────────────┐
│                     Docker Container                         │
│  ┌────────────────────────────────────────────────────────┐ │
│  │              Gunicorn WSGI Server (Port 5000)          │ │
│  │                    (2 workers, 120s timeout)            │ │
│  └────────────────────────────────────────────────────────┘ │
│                             │                                │
│  ┌────────────────────────────────────────────────────────┐ │
│  │                  Flask Application                      │ │
│  │  ┌──────────────┐  ┌──────────────┐  ┌─────────────┐  │ │
│  │  │   API Routes │  │  Validation  │  │   Services  │  │ │
│  │  │   (app.py)   │──│   Module     │──│  (docking,  │  │ │
│  │  │              │  │ (validators) │  │   adme, etc)│  │ │
│  │  └──────────────┘  └──────────────┘  └─────────────┘  │ │
│  │         │                                      │         │ │
│  │  ┌──────────────┐                    ┌─────────────┐   │ │
│  │  │  DB Layer    │                    │ AutoDock    │   │ │
│  │  │ (SQLAlchemy) │                    │ Vina Binary │   │ │
│  │  └──────────────┘                    └─────────────┘   │ │
│  └────────────────────────────────────────────────────────┘ │
│                             │                                │
│  ┌────────────────────────────────────────────────────────┐ │
│  │              Persistent Volumes                         │ │
│  │  ┌─────────────────┐       ┌─────────────────┐        │ │
│  │  │  bindora-data   │       │  bindora-cache  │        │ │
│  │  │  (bindora.db)   │       │  (PDB/SMILES)   │        │ │
│  │  └─────────────────┘       └─────────────────┘        │ │
│  └────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
                             │
                    ┌────────┴────────┐
                    │                 │
              ┌─────────┐      ┌──────────┐
              │ Browser │      │ CLI Tool │
              └─────────┘      └──────────┘
```

### Component Responsibilities

**1. Validation Module (`backend/utils/validators.py`)**
- Input sanitization and validation
- Chemical structure parsing
- Grid box parameter verification
- Error message formatting

**2. Database Layer (`backend/db/`)**
- SQLAlchemy ORM models
- Database session management
- CRUD operations for sessions, jobs, benchmarks
- Transaction handling

**3. API Routes (`backend/app.py` + `backend/routes/`)**
- HTTP request handling
- Input validation orchestration
- Response formatting
- Error handling

**4. Services Layer (existing)**
- Molecular docking (AutoDock Vina)
- ADME profiling (RDKit)
- Structure fetching (PubChem, RCSB)
- Interaction analysis

---

## Database Schema Design

### Entity-Relationship Diagram

```
┌─────────────────────────┐
│   docking_sessions      │
├─────────────────────────┤
│ id (UUID, PK)           │
│ timestamp (DATETIME)    │
│ pdb_id (VARCHAR)        │
│ ligand_name (VARCHAR)   │
│ ligand_smiles (TEXT)    │
│ affinity_kcal (FLOAT)   │
│ rmsd_lb (FLOAT)         │
│ rmsd_ub (FLOAT)         │
│ execution_duration_s    │
│ receptor_pdbqt_path     │
│ ligand_pdbqt_path       │
│ top_pose_pdbqt_path     │
└─────────────────────────┘

┌─────────────────────────┐
│ batch_screening_jobs    │
├─────────────────────────┤
│ id (UUID, PK)           │
│ job_id (VARCHAR UNIQUE) │
│ receptor_pdb_id (VARCHAR│
│ total_ligands (INTEGER) │
│ completed_ligands (INT) │
│ status (VARCHAR)        │
│ created_at (DATETIME)   │
│ completed_at (DATETIME) │
└─────────────────────────┘

┌─────────────────────────┐
│   benchmark_records     │
├─────────────────────────┤
│ id (INTEGER, PK)        │
│ benchmark_type (VARCHAR)│
│ pdb_id (VARCHAR)        │
│ ligand_name (VARCHAR)   │
│ rmsd_angstrom (FLOAT)   │
│ affinity_kcal (FLOAT)   │
│ timestamp (DATETIME)    │
└─────────────────────────┘
```

### Table Specifications

**docking_sessions**
- Primary Key: UUID for global uniqueness
- Indexed Fields: `pdb_id`, `timestamp` (for fast queries)
- Nullable Fields: `rmsd_lb`, `rmsd_ub` (redocking only)

**batch_screening_jobs**
- Primary Key: UUID
- Unique Constraint: `job_id` (client-generated identifier)
- Status Values: "queued", "running", "completed", "cancelled", "interrupted"

**benchmark_records**
- Primary Key: Auto-increment integer
- Composite Index: (`benchmark_type`, `timestamp`)
- Benchmark Types: "casf2016", "dude", "accuracy"

---

## Component Design

### 1. Validation Module

**File:** `backend/utils/validators.py`

**Classes:**
```python
class ValidationResult:
    """Standardized validation result container"""
    valid: bool
    error: Optional[str]
    data: Optional[Dict[str, Any]]
```

**Functions:**

```python
def validate_smiles(smiles: str) -> ValidationResult:
    """
    Validate SMILES string using RDKit parser
    
    Checks:
    - Length ≤ 5000 characters
    - RDKit parseable
    - Sanitizable (valence, aromaticity)
    - Heavy atom count ≥ 1
    
    Returns:
    - valid=True: canonical_smiles, mol_object, heavy_atom_count
    - valid=False: error message
    """

def validate_pdb_content(pdb_text: str) -> ValidationResult:
    """
    Validate PDB block structure
    
    Checks:
    - Length ≤ 10MB
    - Contains ATOM or HETATM records
    - Coordinate bounds within [-999, 999]
    - At least 3 atoms
    
    Returns:
    - valid=True: atom_count, residue_count
    - valid=False: error message
    """

def validate_grid_box(center: Dict, size: Dict) -> ValidationResult:
    """
    Validate docking grid box parameters
    
    Checks:
    - center.x, center.y, center.z are numeric
    - size.x, size.y, size.z are positive
    - 5.0 ≤ size ≤ 150.0 Angstroms
    - center coordinates within [-999, 999]
    
    Returns:
    - valid=True: validated_center, validated_size
    - valid=False: error message
    """
```

**Error Message Format:**
```json
{
  "valid": false,
  "error": "Human-readable description",
  "field": "smiles|pdb_content|center|size",
  "error_code": "INVALID_SMILES|INVALID_PDB|INVALID_GRID"
}
```

---

### 2. Database Layer

**File:** `backend/db/database.py`

```python
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

# Global session factory
engine = None
SessionLocal = None
Base = declarative_base()

def init_db(database_url: str):
    """Initialize database engine and create tables"""
    global engine, SessionLocal
    engine = create_engine(database_url, echo=False)
    SessionLocal = scoped_session(sessionmaker(bind=engine))
    Base.metadata.create_all(engine)

def get_db():
    """Get thread-safe database session"""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
```

**File:** `backend/db/models.py`

```python
from sqlalchemy import Column, String, Float, Integer, DateTime
from sqlalchemy.dialects.postgresql import UUID
import uuid
from datetime import datetime

class DockingSession(Base):
    __tablename__ = "docking_sessions"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    timestamp = Column(DateTime, default=datetime.utcnow, index=True)
    pdb_id = Column(String(10), index=True)
    ligand_name = Column(String(200))
    ligand_smiles = Column(String(5000))
    affinity_kcal = Column(Float)
    rmsd_lb = Column(Float, nullable=True)
    rmsd_ub = Column(Float, nullable=True)
    execution_duration_s = Column(Float)
    receptor_pdbqt_path = Column(String(500), nullable=True)
    ligand_pdbqt_path = Column(String(500), nullable=True)
    top_pose_pdbqt_path = Column(String(500), nullable=True)

class BatchScreeningJob(Base):
    __tablename__ = "batch_screening_jobs"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    job_id = Column(String(100), unique=True, index=True)
    receptor_pdb_id = Column(String(10))
    total_ligands = Column(Integer)
    completed_ligands = Column(Integer, default=0)
    status = Column(String(20), default="queued")
    created_at = Column(DateTime, default=datetime.utcnow)
    completed_at = Column(DateTime, nullable=True)

class BenchmarkRecord(Base):
    __tablename__ = "benchmark_records"
    
    id = Column(Integer, primary_key=True)
    benchmark_type = Column(String(50), index=True)
    pdb_id = Column(String(10))
    ligand_name = Column(String(200))
    rmsd_angstrom = Column(Float)
    affinity_kcal = Column(Float)
    timestamp = Column(DateTime, default=datetime.utcnow, index=True)
```

---

### 3. API Routes Enhancement

**File:** `backend/routes/session_routes.py`

```python
from flask import Blueprint, jsonify, request
from backend.db.database import get_db
from backend.db.models import DockingSession

session_bp = Blueprint('sessions', __name__, url_prefix='/api/sessions')

@session_bp.route('/list', methods=['GET'])
def list_sessions():
    """Paginated session history"""
    limit = int(request.args.get('limit', 50))
    offset = int(request.args.get('offset', 0))
    
    db = next(get_db())
    sessions = db.query(DockingSession)\
        .order_by(DockingSession.timestamp.desc())\
        .limit(limit)\
        .offset(offset)\
        .all()
    
    return jsonify([{
        'id': str(s.id),
        'timestamp': s.timestamp.isoformat(),
        'pdb_id': s.pdb_id,
        'ligand_name': s.ligand_name,
        'affinity_kcal': s.affinity_kcal
    } for s in sessions])

@session_bp.route('/<session_id>', methods=['GET'])
def get_session(session_id):
    """Get session details"""
    db = next(get_db())
    session = db.query(DockingSession).filter_by(id=session_id).first()
    
    if not session:
        return jsonify({'error': 'Session not found'}), 404
    
    return jsonify({
        'id': str(session.id),
        'timestamp': session.timestamp.isoformat(),
        'pdb_id': session.pdb_id,
        'ligand_name': session.ligand_name,
        'ligand_smiles': session.ligand_smiles,
        'affinity_kcal': session.affinity_kcal,
        'rmsd_lb': session.rmsd_lb,
        'rmsd_ub': session.rmsd_ub,
        'execution_duration_s': session.execution_duration_s
    })

@session_bp.route('/<session_id>', methods=['DELETE'])
def delete_session(session_id):
    """Delete session"""
    db = next(get_db())
    session = db.query(DockingSession).filter_by(id=session_id).first()
    
    if not session:
        return jsonify({'error': 'Session not found'}), 404
    
    db.delete(session)
    db.commit()
    
    return '', 204
```

---

## Low-Level Design

### Security Configuration Flow

```
Application Startup
    │
    ├─> Load environment variables (.env)
    │
    ├─> backend/config.py
    │   │
    │   ├─> DEBUG = os.getenv("BINDORA_DEBUG", "False") != "False"
    │   ├─> MAX_CONTENT_LENGTH = 32 * 1024 * 1024
    │   ├─> CORS_ORIGINS = os.getenv("CORS_ORIGINS", "http://localhost:5000")
    │   └─> DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///data/bindora.db")
    │
    ├─> backend/app.py
    │   │
    │   ├─> app.config["MAX_CONTENT_LENGTH"] = MAX_CONTENT_LENGTH
    │   ├─> CORS(app, origins=CORS_ORIGINS.split(","))
    │   └─> @app.errorhandler(413) -> "Payload too large"
    │
    └─> Ready to accept requests
```

### Validation Middleware Flow

```
HTTP Request → /api/structure/ligand
    │
    ├─> Extract JSON body
    │
    ├─> validators.validate_smiles(body["structure"])
    │   │
    │   ├─> Check length ≤ 5000
    │   ├─> RDKit parse: Chem.MolFromSmiles()
    │   ├─> Sanitize: Chem.SanitizeMol()
    │   ├─> Count heavy atoms
    │   │
    │   ├─> [VALID] → Return {valid: true, canonical_smiles: ...}
    │   └─> [INVALID] → Return {valid: false, error: "..."}
    │
    ├─> If valid=false
    │   └─> Return HTTP 400 with error message
    │
    └─> If valid=true
        └─> Proceed to DockingEngine.prepare_ligand()
```

### Database Transaction Pattern

```python
def record_docking_session(session_data: dict):
    """
    Record completed docking session with rollback safety
    """
    db = next(get_db())
    try:
        session = DockingSession(
            pdb_id=session_data.get('pdb_id'),
            ligand_name=session_data['ligand_name'],
            ligand_smiles=session_data['ligand_smiles'],
            affinity_kcal=session_data['affinity_kcal'],
            execution_duration_s=session_data['duration']
        )
        db.add(session)
        db.commit()
        return str(session.id)
    except Exception as e:
        db.rollback()
        raise e
    finally:
        db.close()
```

---

## Docker Architecture

### Multi-Stage Build Strategy

```dockerfile
# Stage 1: Base image with system dependencies
FROM python:3.11-slim AS base
RUN apt-get update && apt-get install -y \
    wget ca-certificates libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Stage 2: Vina binary installation
FROM base AS vina
WORKDIR /tmp
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    && mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina \
    && chmod +x /usr/local/bin/vina

# Stage 3: Python dependencies
FROM base AS dependencies
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Stage 4: Final application
FROM dependencies AS final
COPY --from=vina /usr/local/bin/vina /usr/local/bin/vina
RUN useradd -m -u 1000 bindora
WORKDIR /app
COPY --chown=bindora:bindora . /app
USER bindora
EXPOSE 5000
CMD ["gunicorn", "--config", "gunicorn_config.py", "backend.app:app"]
```

### Volume Mount Strategy

```yaml
volumes:
  bindora-data:
    # Stores bindora.db (SQLite database)
    # Persists across container restarts
    # Backed up via docker volume backup
  
  bindora-cache:
    # Stores cached PDB structures and API responses
    # Can be cleared without data loss
    # Improves performance for repeated queries
```

---

## Security Architecture

### Defense in Depth Layers

1. **Input Validation Layer**
   - SMILES length limits (5000 chars)
   - PDB size limits (10MB)
   - Grid coordinate bounds checking
   - Sanitization before processing

2. **Network Security Layer**
   - CORS origin whitelisting
   - Payload size limits (32MB)
   - Rate limiting (future enhancement)

3. **Application Security Layer**
   - DEBUG=False in production
   - No secrets in source code
   - Environment variable configuration
   - Error messages without stack traces

4. **Container Security Layer**
   - Non-root user execution (UID 1000)
   - Minimal base image (slim)
   - No unnecessary packages
   - Read-only filesystem (future enhancement)

### Configuration Security Matrix

| Config Variable | Default | Production Override | Security Impact |
|----------------|---------|---------------------|-----------------|
| DEBUG | False | False | Prevents stack trace leaks |
| CORS_ORIGINS | localhost:5000 | https://bindora.example.com | Prevents CSRF |
| MAX_CONTENT_LENGTH | 32MB | 32MB | Prevents DoS via large payloads |
| DATABASE_URL | sqlite:///data/bindora.db | sqlite:///data/bindora.db | Local file, no network exposure |

---

## Performance Considerations

### Database Indexing Strategy
- `docking_sessions.timestamp` (DESC) → Fast recent session queries
- `docking_sessions.pdb_id` → Fast filtering by target
- `batch_screening_jobs.job_id` (UNIQUE) → O(1) job lookup
- `benchmark_records.benchmark_type` → Fast benchmark filtering

### Connection Pooling
- SQLAlchemy scoped_session → Thread-safe session reuse
- Gunicorn worker count = 2 → Balance concurrency vs memory
- Worker timeout = 120s → Allow long docking calculations

### Caching Strategy (Existing)
- PDB structures cached in `data/cache/`
- PubChem API responses cached
- No database query caching (future enhancement)

---

## Deployment Architecture

### Production Stack

```
Internet
    │
    ├─> Reverse Proxy (Nginx/Traefik) [Future]
    │       │
    │       ├─> SSL/TLS Termination
    │       ├─> Rate Limiting
    │       └─> Load Balancing
    │
    └─> Docker Container (bindora-app)
            │
            ├─> Gunicorn (0.0.0.0:5000)
            │       │
            │       ├─> Worker 1 (Flask App)
            │       └─> Worker 2 (Flask App)
            │
            ├─> Volume: bindora-data (/app/data)
            │       └─> bindora.db (SQLite)
            │
            └─> Volume: bindora-cache (/app/data/cache)
                    └─> PDB files, API responses
```

### Environment Configuration Files

**.env (Production)**
```bash
BINDORA_HOST=0.0.0.0
BINDORA_PORT=5000
BINDORA_DEBUG=False
CORS_ORIGINS=https://bindora.example.com,https://www.bindora.example.com
DATABASE_URL=sqlite:////app/data/bindora.db
VINA_EXE=/usr/local/bin/vina
```

**gunicorn_config.py**
```python
import os

bind = f"{os.getenv('BINDORA_HOST', '0.0.0.0')}:{os.getenv('BINDORA_PORT', '5000')}"
workers = 2
worker_class = "sync"
timeout = 120
accesslog = "-"
errorlog = "-"
loglevel = "info"
```

---

## Backward Compatibility Strategy

### API Response Schema Preservation
- All existing endpoints return identical JSON structure
- New fields added as optional (not breaking)
- Database persistence invisible to API consumers

### CLI Compatibility
- All existing commands work unchanged
- New --db-session-id flag added (optional)
- Output format identical

### File Format Compatibility
- PDBQT files unchanged
- JSON report format unchanged
- Markdown benchmark reports unchanged

---

## Testing Strategy

### Unit Tests
- `test_validators.py` → All validation functions
- `test_database.py` → ORM models and queries
- `test_api_routes.py` → Session management endpoints

### Integration Tests
- `test_docking_with_db.py` → Docking + automatic session recording
- `test_batch_with_db.py` → Batch screening + job persistence

### System Tests
- `test_docker_build.py` → Container build verification
- `test_docker_compose.py` → Full stack startup
- `test_production_config.py` → Security settings verification

---

## Migration Path

### Phase 1: Code Changes (No Breaking)
1. Add validators.py
2. Add db/ directory with models
3. Modify app.py to add validation
4. Add session recording hooks

### Phase 2: Database Initialization
1. Run `python -m backend.db.database` to create tables
2. Verify bindora.db created
3. Test session recording

### Phase 3: Docker Deployment
1. Build Docker image
2. Test locally with docker-compose
3. Deploy to production

### Phase 4: Monitoring & Validation
1. Verify logs for errors
2. Check database growth
3. Monitor response times
4. Validate security settings

---

## Appendix: File Changes Summary

| File | Action | Description |
|------|--------|-------------|
| requirements.txt | MODIFY | Add gemmi, scipy, sqlalchemy, gunicorn |
| backend/config.py | MODIFY | Add security configs, database URL |
| backend/utils/validators.py | NEW | Input validation module |
| backend/db/database.py | NEW | SQLAlchemy setup |
| backend/db/models.py | NEW | ORM models |
| backend/routes/session_routes.py | NEW | Session API endpoints |
| backend/app.py | MODIFY | Add validation, register blueprints |
| Dockerfile | NEW | Multi-stage container build |
| docker-compose.yml | NEW | Orchestration configuration |
| .dockerignore | NEW | Exclude unnecessary files |
| gunicorn_config.py | NEW | Production WSGI config |
| tests/test_hardening.py | NEW | Validation tests |
| tests/test_database.py | NEW | Database tests |
| tests/test_docker.py | NEW | Container tests |

