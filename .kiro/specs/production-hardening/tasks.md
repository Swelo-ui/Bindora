# Implementation Tasks: Production Hardening & Architectural Upgrades

## Task Breakdown

### Phase 1: Dependencies & Security (Tasks 1-6)

#### Task 1: Complete requirements.txt with missing dependencies
**Status:** pending  
**Dependencies:** None  
**Files:**
- `requirements.txt`

**Acceptance Criteria:**
- Add gemmi>=0.6.5
- Add scipy>=1.11.0
- Add python-dotenv>=1.0.0
- Add colorama>=0.4.6
- Add sqlalchemy>=2.0.0
- Add gunicorn>=21.2.0

---

#### Task 2: Update backend/config.py with security hardening
**Status:** pending  
**Dependencies:** Task 1  
**Files:**
- `backend/config.py`

**Acceptance Criteria:**
- Change DEBUG default to False
- Add MAX_CONTENT_LENGTH = 32 * 1024 * 1024
- Add CORS_ORIGINS from environment with default "http://localhost:5000,http://127.0.0.1:5000"
- Add DATABASE_URL with default "sqlite:///data/bindora.db"

---

#### Task 3: Create input validation module
**Status:** pending  
**Dependencies:** Task 1  
**Files:**
- `backend/utils/validators.py` (NEW)

**Acceptance Criteria:**
- Implement ValidationResult class
- Implement validate_smiles() function
- Implement validate_pdb_content() function
- Implement validate_grid_box() function
- All functions return ValidationResult with proper error messages

---

#### Task 4: Integrate validation into API endpoints
**Status:** pending  
**Dependencies:** Task 2, Task 3  
**Files:**
- `backend/app.py`

**Acceptance Criteria:**
- Apply MAX_CONTENT_LENGTH to app.config
- Add @app.errorhandler(413) for payload too large
- Add validate_smiles to /api/structure/ligand endpoint
- Add validate_pdb_content to /api/structure/receptor endpoint
- Add validate_grid_box to /api/docking/run endpoint
- Add validate_smiles to /api/pkpd/adme endpoint
- All invalid inputs return HTTP 400 with error message

---

#### Task 5: Configure CORS restrictions
**Status:** pending  
**Dependencies:** Task 2, Task 4  
**Files:**
- `backend/app.py`

**Acceptance Criteria:**
- Replace CORS(app) with CORS(app, origins=CORS_ORIGINS.split(","))
- Verify CORS headers in response
- Test cross-origin request rejection

---

#### Task 6: Create validation test suite
**Status:** pending  
**Dependencies:** Task 3, Task 4, Task 5  
**Files:**
- `tests/test_hardening.py` (NEW)

**Acceptance Criteria:**
- Test validate_smiles with valid/invalid inputs
- Test validate_pdb_content with valid/invalid inputs
- Test validate_grid_box with valid/invalid inputs
- Test API endpoints reject invalid inputs with HTTP 400
- Test MAX_CONTENT_LENGTH triggers HTTP 413
- Test CORS origin restrictions
- Test DEBUG defaults to False

---

### Phase 2: Database Architecture (Tasks 7-12)

#### Task 7: Create database initialization module
**Status:** pending  
**Dependencies:** Task 1  
**Files:**
- `backend/db/__init__.py` (NEW)
- `backend/db/database.py` (NEW)

**Acceptance Criteria:**
- Implement init_db() function
- Implement get_db() context manager
- Use scoped_session for thread safety
- Create data/bindora.db on first call

---

#### Task 8: Define ORM models
**Status:** pending  
**Dependencies:** Task 7  
**Files:**
- `backend/db/models.py` (NEW)

**Acceptance Criteria:**
- Define DockingSession model with all required columns
- Define BatchScreeningJob model with all required columns
- Define BenchmarkRecord model with all required columns
- Add proper indexes (timestamp, pdb_id, job_id)
- Add UUID primary keys

---

#### Task 9: Initialize database on application startup
**Status:** pending  
**Dependencies:** Task 2, Task 7, Task 8  
**Files:**
- `backend/app.py`

**Acceptance Criteria:**
- Import init_db from backend.db.database
- Call init_db(DATABASE_URL) before routes
- Verify bindora.db created in data/ directory
- Verify tables created with correct schema

---

#### Task 10: Implement session recording in docking endpoint
**Status:** pending  
**Dependencies:** Task 8, Task 9  
**Files:**
- `backend/app.py` (modify /api/docking/run)

**Acceptance Criteria:**
- After successful docking, create DockingSession record
- Save pdb_id, ligand_name, ligand_smiles, affinity, duration
- Handle database errors with rollback
- Return session_id in response

---

#### Task 11: Create session management API routes
**Status:** pending  
**Dependencies:** Task 8, Task 9  
**Files:**
- `backend/routes/session_routes.py` (NEW)
- `backend/app.py` (register blueprint)

**Acceptance Criteria:**
- Implement GET /api/sessions/list with pagination
- Implement GET /api/sessions/<session_id>
- Implement DELETE /api/sessions/<session_id>
- Register session_bp blueprint in app.py
- All endpoints return proper JSON responses

---

#### Task 12: Create database test suite
**Status:** pending  
**Dependencies:** Task 7, Task 8, Task 10, Task 11  
**Files:**
- `tests/test_database.py` (NEW)

**Acceptance Criteria:**
- Test database initialization
- Test DockingSession CRUD operations
- Test BatchScreeningJob CRUD operations
- Test BenchmarkRecord CRUD operations
- Test session list pagination
- Test session retrieval by ID
- Test session deletion
- Test transaction rollback on errors

---

### Phase 3: Production Deployment (Tasks 13-18)

#### Task 13: Create Gunicorn production configuration
**Status:** pending  
**Dependencies:** Task 1  
**Files:**
- `gunicorn_config.py` (NEW)

**Acceptance Criteria:**
- Set bind from BINDORA_HOST and BINDORA_PORT
- Configure 2 workers
- Set timeout to 120 seconds
- Configure logging to stdout/stderr
- Set worker_class to "sync"

---

#### Task 14: Create multi-stage Dockerfile
**Status:** pending  
**Dependencies:** Task 1, Task 13  
**Files:**
- `Dockerfile` (NEW)

**Acceptance Criteria:**
- Use python:3.11-slim base image
- Install system dependencies (wget, ca-certificates, libgomp1)
- Download AutoDock Vina 1.2.5 Linux binary
- Install to /usr/local/bin/vina with +x permissions
- Create non-root user "bindora" (UID 1000)
- Copy requirements.txt and install dependencies
- Copy application code
- Set workdir to /app
- Switch to bindora user
- Expose port 5000
- Set CMD to gunicorn with config file

---

#### Task 15: Create Docker Compose orchestration
**Status:** pending  
**Dependencies:** Task 14  
**Files:**
- `docker-compose.yml` (NEW)

**Acceptance Criteria:**
- Define bindora-app service
- Map port 5000:5000
- Define bindora-data named volume mounted to /app/data
- Define bindora-cache named volume mounted to /app/data/cache
- Pass BINDORA_HOST, BINDORA_PORT, BINDORA_DEBUG env vars
- Load .env file if present

---

#### Task 16: Create .dockerignore file
**Status:** pending  
**Dependencies:** None  
**Files:**
- `.dockerignore` (NEW)

**Acceptance Criteria:**
- Exclude .git, __pycache__, *.pyc
- Exclude tests/, .pytest_cache
- Exclude bin/vina.exe (Windows binary)
- Exclude data/cache/* (will be volume)
- Exclude .env (use docker-compose env)

---

#### Task 17: Add container health check
**Status:** pending  
**Dependencies:** Task 14  
**Files:**
- `Dockerfile`

**Acceptance Criteria:**
- Add HEALTHCHECK instruction
- Use wget to check http://localhost:5000/api/health
- Set interval to 30 seconds
- Set timeout to 5 seconds
- Set start-period to 10 seconds
- Set retries to 3

---

#### Task 18: Create Docker deployment test suite
**Status:** pending  
**Dependencies:** Task 14, Task 15, Task 17  
**Files:**
- `tests/test_docker.py` (NEW)

**Acceptance Criteria:**
- Test Dockerfile builds successfully
- Test container starts and responds to health check
- Test volume mounts persist data
- Test environment variable configuration
- Test non-root user execution

---

### Phase 4: Integration & Validation (Tasks 19-21)

#### Task 19: Update .env.example with new variables
**Status:** pending  
**Dependencies:** Task 2  
**Files:**
- `.env.example`

**Acceptance Criteria:**
- Add BINDORA_DEBUG=False
- Add CORS_ORIGINS example
- Add DATABASE_URL example
- Add MAX_CONTENT_LENGTH example
- Add comments explaining each variable

---

#### Task 20: Run full test suite and fix regressions
**Status:** pending  
**Dependencies:** Task 6, Task 12, Task 18  
**Files:**
- All test files

**Acceptance Criteria:**
- Run pytest tests/ -v
- All existing tests pass
- All new hardening tests pass
- All database tests pass
- All Docker tests pass (if Docker available)
- No import errors
- No deprecation warnings

---

#### Task 21: Update documentation
**Status:** pending  
**Dependencies:** Task 20  
**Files:**
- `README.md`
- `CONTRIBUTING.md`

**Acceptance Criteria:**
- Update requirements section with new dependencies
- Add database configuration section
- Add Docker deployment instructions
- Add production configuration guide
- Add security best practices section
- Update quickstart to include Docker option

---

## Task Dependencies Graph

```
Task 1 (requirements.txt)
    ├─> Task 2 (config.py security)
    │   ├─> Task 4 (API validation)
    │   │   ├─> Task 5 (CORS)
    │   │   │   └─> Task 6 (validation tests)
    │   │   └─> Task 10 (session recording)
    │   ├─> Task 9 (db init)
    │   └─> Task 19 (.env.example)
    │
    ├─> Task 3 (validators.py)
    │   └─> Task 4 (API validation)
    │
    ├─> Task 7 (database.py)
    │   ├─> Task 8 (models.py)
    │   │   ├─> Task 9 (db init)
    │   │   ├─> Task 10 (session recording)
    │   │   ├─> Task 11 (session routes)
    │   │   │   └─> Task 12 (db tests)
    │   │   └─> Task 12 (db tests)
    │
    └─> Task 13 (gunicorn config)
        └─> Task 14 (Dockerfile)
            ├─> Task 15 (docker-compose)
            │   └─> Task 18 (docker tests)
            └─> Task 17 (healthcheck)
                └─> Task 18 (docker tests)

Task 6, Task 12, Task 18 ──> Task 20 (full test suite)
Task 20 ──> Task 21 (documentation)
```

## Execution Order

### Parallel Batch 1 (No dependencies)
- Task 1: requirements.txt
- Task 16: .dockerignore

### Parallel Batch 2 (Depends on Task 1)
- Task 2: config.py
- Task 3: validators.py
- Task 7: database.py
- Task 13: gunicorn config

### Sequential Batch 3
- Task 8: models.py (depends on Task 7)
- Task 9: db init (depends on Task 2, 7, 8)

### Parallel Batch 4
- Task 4: API validation (depends on Task 2, 3)
- Task 10: session recording (depends on Task 8, 9)
- Task 11: session routes (depends on Task 8, 9)
- Task 14: Dockerfile (depends on Task 1, 13)

### Parallel Batch 5
- Task 5: CORS (depends on Task 2, 4)
- Task 15: docker-compose (depends on Task 14)
- Task 17: healthcheck (depends on Task 14)

### Parallel Batch 6 (Testing)
- Task 6: validation tests (depends on Task 3, 4, 5)
- Task 12: database tests (depends on Task 7, 8, 10, 11)
- Task 18: docker tests (depends on Task 14, 15, 17)

### Sequential Batch 7
- Task 19: .env.example (depends on Task 2)
- Task 20: full test suite (depends on Task 6, 12, 18)
- Task 21: documentation (depends on Task 20)

## Estimated Effort

| Phase | Tasks | Effort | Priority |
|-------|-------|--------|----------|
| Phase 1: Security | 1-6 | 4-6 hours | Critical |
| Phase 2: Database | 7-12 | 6-8 hours | High |
| Phase 3: Docker | 13-18 | 4-5 hours | Medium |
| Phase 4: Integration | 19-21 | 2-3 hours | Low |
| **Total** | **21 tasks** | **16-22 hours** | - |

## Risk Mitigation

### High-Risk Tasks
- Task 4: API validation integration (breaking changes possible)
- Task 10: Session recording (database transaction errors)
- Task 14: Dockerfile (platform-specific issues)

### Mitigation Strategies
- Run tests after each task
- Use git branches for each phase
- Test database migrations on copy of production data
- Test Docker build on Linux, macOS, Windows if possible

