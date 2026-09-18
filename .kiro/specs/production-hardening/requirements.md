# Requirements Document

## Introduction

This document specifies requirements for transforming Bindora Dock from a research-grade prototype into a hardened, production-ready, containerized Computer-Aided Drug Design (CADD) platform. The transformation addresses critical gaps in dependency management, security posture, data persistence, and deployment readiness while maintaining full backward compatibility with existing functionality.

## Glossary

- **Application**: The Bindora Dock system consisting of Flask backend, frontend assets, and AutoDock Vina integration
- **Container**: A Docker container encapsulating the Application and its runtime dependencies
- **Database**: SQLite relational database storing session history and computational results
- **Input_Validator**: Centralized validation module for chemical structures and computational parameters
- **ORM**: SQLAlchemy Object-Relational Mapper providing database abstraction
- **Request**: HTTP request received by the Application API endpoints
- **Session**: A DockingSession database record representing a completed molecular docking simulation
- **WSGI_Server**: Gunicorn production WSGI HTTP server
- **User**: A researcher, scientist, or developer interacting with the Application

---

## Requirements

### Requirement 1: Dependency Completeness and Environment Configuration

**User Story:** As a developer, I want all runtime dependencies explicitly declared and installable, so that the Application runs without import errors or missing module failures.

#### Acceptance Criteria

1. THE Application SHALL include gemmi in requirements.txt with minimum version 0.6.0
2. THE Application SHALL include scipy in requirements.txt with minimum version 1.10.0
3. THE Application SHALL include python-dotenv in requirements.txt with minimum version 1.0.0
4. THE Application SHALL include colorama in requirements.txt with minimum version 0.4.0
5. THE Application SHALL include sqlalchemy in requirements.txt with minimum version 2.0.0
6. THE Application SHALL include gunicorn in requirements.txt with minimum version 21.0.0
7. WHEN the Application starts, THE Application SHALL successfully import all dependencies without ImportError
8. WHEN requirements.txt is installed via pip, THE Application SHALL start without dependency-related errors

---

### Requirement 2: Security Hardening and Default Configuration

**User Story:** As a security engineer, I want the Application to follow security best practices by default, so that deployments are not vulnerable to common web application attacks.

#### Acceptance Criteria

1. THE Application SHALL set DEBUG to False by default in backend/config.py
2. WHERE the BINDORA_DEBUG environment variable equals "true" or "1", THE Application SHALL enable DEBUG mode
3. THE Application SHALL define MAX_CONTENT_LENGTH configuration limiting request payloads to 33554432 bytes (32 MB)
4. WHEN a Request exceeds MAX_CONTENT_LENGTH, THE Application SHALL return HTTP 413 status code
5. THE Application SHALL define CORS_ORIGINS configuration variable for Cross-Origin Resource Sharing control
6. WHERE CORS_ORIGINS environment variable is unset, THE Application SHALL default to "http://localhost:5000,http://127.0.0.1:5000"
7. THE Application SHALL configure Flask-CORS to restrict origins to CORS_ORIGINS list
8. WHEN a Request originates from an origin not in CORS_ORIGINS, THE Application SHALL reject the Request with appropriate CORS headers

---

### Requirement 3: Input Validation Infrastructure

**User Story:** As a developer, I want centralized validation functions for chemical and computational inputs, so that invalid data is rejected before processing and error messages are consistent.

#### Acceptance Criteria

1. THE Application SHALL provide a validate_smiles function accepting a SMILES string parameter
2. WHEN validate_smiles receives a valid SMILES string, THE validate_smiles function SHALL return a dictionary with "valid" key set to True and "canonical_smiles" key containing the canonicalized SMILES
3. WHEN validate_smiles receives an invalid SMILES string, THE validate_smiles function SHALL return a dictionary with "valid" key set to False and "error" key containing a descriptive error message
4. THE Application SHALL provide a validate_pdb_content function accepting a PDB block string parameter
5. WHEN validate_pdb_content receives a string containing valid ATOM or HETATM records, THE validate_pdb_content function SHALL return a dictionary with "valid" key set to True
6. WHEN validate_pdb_content receives a string without ATOM or HETATM records, THE validate_pdb_content function SHALL return a dictionary with "valid" key set to False and "error" key containing a descriptive error message
7. THE Application SHALL provide a validate_grid_box function accepting center coordinates and size dimensions parameters
8. WHEN validate_grid_box receives center coordinates as a list of three numeric values and size dimensions as a list of three positive numeric values, THE validate_grid_box function SHALL return a dictionary with "valid" key set to True
9. WHEN validate_grid_box receives center coordinates with non-numeric values or size dimensions with non-positive values, THE validate_grid_box function SHALL return a dictionary with "valid" key set to False and "error" key containing a descriptive error message

---

### Requirement 4: API Endpoint Input Validation

**User Story:** As a system administrator, I want all API endpoints to validate inputs before processing, so that malformed requests are rejected early and do not cause server errors.

#### Acceptance Criteria

1. WHEN the /api/structure/ligand endpoint receives a Request, THE Application SHALL validate the "structure" field using validate_smiles before ligand preparation
2. WHEN the /api/structure/ligand endpoint receives invalid SMILES, THE Application SHALL return HTTP 400 status code with validation error message
3. WHEN the /api/structure/receptor endpoint receives a Request, THE Application SHALL validate the "pdb_content" field using validate_pdb_content before receptor preparation
4. WHEN the /api/structure/receptor endpoint receives invalid PDB content, THE Application SHALL return HTTP 400 status code with validation error message
5. WHEN the /api/docking/run endpoint receives a Request, THE Application SHALL validate the "center" and "size" fields using validate_grid_box before docking execution
6. WHEN the /api/docking/run endpoint receives invalid grid box parameters, THE Application SHALL return HTTP 400 status code with validation error message
7. WHEN the /api/pkpd/adme endpoint receives a Request, THE Application SHALL validate the "smiles" field using validate_smiles before ADME calculation
8. WHEN the /api/pkpd/adme endpoint receives invalid SMILES, THE Application SHALL return HTTP 400 status code with validation error message

---

### Requirement 5: Database Architecture and ORM Layer

**User Story:** As a developer, I want a relational database with SQLAlchemy ORM, so that session history and computational results are persisted across Application restarts.

#### Acceptance Criteria

1. THE Application SHALL use SQLAlchemy version 2.0 or higher as the ORM layer
2. THE Application SHALL create a SQLite database file at data/bindora.db on first startup
3. THE Database SHALL define a docking_sessions table with columns: id (primary key), timestamp (datetime), pdb_id (string), ligand_name (string), ligand_smiles (string), affinity_kcal (float), rmsd_lb (float), rmsd_ub (float), execution_duration_s (float), receptor_pdbqt_path (string), ligand_pdbqt_path (string), top_pose_pdbqt_path (string)
4. THE Database SHALL define a batch_screening_jobs table with columns: id (primary key), job_id (string unique), receptor_pdb_id (string), total_ligands (integer), completed_ligands (integer), status (string), created_at (datetime), completed_at (datetime nullable)
5. THE Database SHALL define a benchmark_records table with columns: id (primary key), benchmark_type (string), pdb_id (string), ligand_name (string), rmsd_angstrom (float), affinity_kcal (float), timestamp (datetime)
6. THE ORM SHALL use scoped_session for thread-safe database access
7. WHEN the Application starts, THE ORM SHALL automatically create all tables if they do not exist
8. WHEN a database query fails, THE ORM SHALL rollback the transaction and not corrupt the Database

---

### Requirement 6: Session History Recording and Management

**User Story:** As a researcher, I want completed docking simulations automatically saved to the database, so that I can review past experiments without re-running computations.

#### Acceptance Criteria

1. WHEN the /api/docking/run endpoint completes successfully, THE Application SHALL create a Session record in the docking_sessions table
2. THE Session record SHALL include timestamp, ligand SMILES, affinity in kcal/mol, RMSD bounds, and execution duration
3. WHERE the receptor was fetched from RCSB PDB, THE Session record SHALL include the pdb_id field
4. THE Application SHALL provide a /api/sessions/list endpoint accepting optional limit and offset query parameters
5. WHEN the /api/sessions/list endpoint receives a Request, THE Application SHALL return a JSON array of Session records ordered by timestamp descending
6. THE Application SHALL provide a /api/sessions/<session_id> endpoint accepting a session ID path parameter
7. WHEN the /api/sessions/<session_id> endpoint receives a Request with a valid session_id, THE Application SHALL return the complete Session record as JSON
8. WHEN the /api/sessions/<session_id> endpoint receives a Request with an invalid session_id, THE Application SHALL return HTTP 404 status code
9. THE Application SHALL provide a /api/sessions/<session_id> DELETE endpoint accepting a session ID path parameter
10. WHEN the /api/sessions/<session_id> DELETE endpoint receives a Request with a valid session_id, THE Application SHALL delete the Session record and return HTTP 204 status code

---

### Requirement 7: Batch Screening Job Persistence

**User Story:** As a researcher, I want batch screening jobs tracked in the database, so that job status persists across Application restarts and I can resume interrupted screens.

#### Acceptance Criteria

1. WHEN the /api/batch/start endpoint creates a batch job, THE Application SHALL create a batch_screening_jobs record with status "queued"
2. WHEN a batch screening job starts processing, THE Application SHALL update the batch_screening_jobs record status to "running"
3. WHEN a batch screening job completes all ligands, THE Application SHALL update the batch_screening_jobs record status to "completed" and set completed_at timestamp
4. WHEN a batch screening job is cancelled, THE Application SHALL update the batch_screening_jobs record status to "cancelled"
5. WHEN the /api/batch/status/<job_id> endpoint receives a Request, THE Application SHALL query the batch_screening_jobs table and return the persisted job status
6. WHERE the Application restarts WHILE a batch job status is "running", THE Application SHALL change the status to "interrupted" on startup

---

### Requirement 8: Benchmark Record Persistence

**User Story:** As a validation engineer, I want benchmark results automatically recorded in the database, so that validation history is preserved and can be queried for reproducibility verification.

#### Acceptance Criteria

1. WHEN a CASF-2016 benchmark run completes a complex redocking, THE Application SHALL create a benchmark_records entry with benchmark_type "casf2016"
2. WHEN a DUD-E virtual screening benchmark completes, THE Application SHALL create benchmark_records entries with benchmark_type "dude"
3. THE benchmark_records entry SHALL include pdb_id, ligand_name, rmsd_angstrom, affinity_kcal, and timestamp fields
4. THE Application SHALL provide a /api/benchmarks/history endpoint accepting optional benchmark_type query parameter
5. WHEN the /api/benchmarks/history endpoint receives a Request, THE Application SHALL return a JSON array of benchmark_records ordered by timestamp descending
6. WHERE the benchmark_type query parameter is provided, THE /api/benchmarks/history endpoint SHALL filter results to matching benchmark_type

---

### Requirement 9: Production WSGI Server Configuration

**User Story:** As a DevOps engineer, I want the Application to run under a production WSGI server with appropriate worker configuration, so that the deployment can handle concurrent requests reliably.

#### Acceptance Criteria

1. THE Application SHALL include a gunicorn_config.py configuration file at the project root
2. THE gunicorn_config.py file SHALL configure 2 worker processes
3. THE gunicorn_config.py file SHALL configure worker_class as "sync"
4. THE gunicorn_config.py file SHALL configure timeout to 120 seconds
5. THE gunicorn_config.py file SHALL configure bind address from HOST and PORT environment variables with defaults "0.0.0.0:5000"
6. THE gunicorn_config.py file SHALL configure accesslog to "-" for stdout logging
7. THE gunicorn_config.py file SHALL configure errorlog to "-" for stderr logging
8. THE Application SHALL include a start_production.sh script that launches gunicorn with gunicorn_config.py
9. WHEN start_production.sh is executed, THE WSGI_Server SHALL start and serve the Application on the configured bind address

---

### Requirement 10: Multi-Stage Docker Container

**User Story:** As a DevOps engineer, I want a multi-stage Dockerfile that produces a minimal production container, so that the Container image size is optimized and attack surface is reduced.

#### Acceptance Criteria

1. THE Container SHALL use python:3.11-slim as the base image
2. THE Container SHALL install system dependencies: wget, ca-certificates
3. THE Container SHALL download AutoDock Vina version 1.2.7 Linux x86_64 binary from official GitHub releases
4. THE Container SHALL install AutoDock Vina binary to /usr/local/bin/vina with executable permissions
5. THE Container SHALL create a non-root user named "bindora" with UID 1000
6. THE Container SHALL set working directory to /app
7. THE Container SHALL copy requirements.txt and install Python dependencies
8. THE Container SHALL copy Application source code to /app
9. THE Container SHALL change ownership of /app to bindora user
10. THE Container SHALL switch to bindora user before starting the Application
11. THE Container SHALL expose port 5000
12. THE Container SHALL set CMD to execute gunicorn using gunicorn_config.py configuration

---

### Requirement 11: Docker Compose Orchestration

**User Story:** As a developer, I want a docker-compose.yml configuration for single-command deployment, so that I can start the entire Application stack without manual container management.

#### Acceptance Criteria

1. THE Application SHALL provide a docker-compose.yml file at the project root
2. THE docker-compose.yml file SHALL define a service named "bindora-app"
3. THE bindora-app service SHALL build from the Dockerfile in the current directory
4. THE bindora-app service SHALL map host port 5000 to Container port 5000
5. THE bindora-app service SHALL define a named volume "bindora-data" mounted to /app/data
6. THE bindora-app service SHALL define a named volume "bindora-cache" mounted to /app/data/cache
7. THE bindora-app service SHALL pass BINDORA_HOST environment variable with default "0.0.0.0"
8. THE bindora-app service SHALL pass BINDORA_PORT environment variable with default "5000"
9. THE bindora-app service SHALL pass BINDORA_DEBUG environment variable with default "False"
10. WHERE a .env file exists in the project root, THE docker-compose.yml SHALL load environment variables from .env file
11. WHEN "docker-compose up" is executed, THE Container SHALL start and the Application SHALL be accessible on http://localhost:5000

---

### Requirement 12: Persistent Volume Data Management

**User Story:** As a system administrator, I want application data stored in Docker named volumes, so that data persists across Container restarts and upgrades.

#### Acceptance Criteria

1. WHEN the Container writes to /app/data/bindora.db, THE Database file SHALL persist in the bindora-data volume
2. WHEN the Container writes to /app/data/cache, THE cache files SHALL persist in the bindora-cache volume
3. WHEN the Container is stopped and started, THE Database SHALL retain all Session records
4. WHEN the Container is removed and recreated, THE named volumes SHALL preserve existing data
5. WHEN "docker-compose down -v" is executed, THE named volumes SHALL be deleted and data cleared

---

### Requirement 13: Container Health Monitoring

**User Story:** As a DevOps engineer, I want the Container to expose a health check endpoint, so that orchestration systems can monitor Application availability.

#### Acceptance Criteria

1. THE Dockerfile SHALL define a HEALTHCHECK instruction
2. THE HEALTHCHECK instruction SHALL execute a curl or wget request to http://localhost:5000/api/health every 30 seconds
3. THE HEALTHCHECK instruction SHALL have a timeout of 5 seconds
4. THE HEALTHCHECK instruction SHALL have a start period of 10 seconds
5. WHEN the Application is running and healthy, THE HEALTHCHECK SHALL return exit code 0
6. WHEN the Application is not responding, THE HEALTHCHECK SHALL return non-zero exit code after 3 consecutive failures

---

### Requirement 14: Environment Variable Configuration Precedence

**User Story:** As a DevOps engineer, I want environment variables to override default configuration, so that I can customize deployment settings without modifying source code.

#### Acceptance Criteria

1. WHERE the BINDORA_HOST environment variable is set, THE Application SHALL bind to the specified host address instead of the default
2. WHERE the BINDORA_PORT environment variable is set, THE Application SHALL bind to the specified port instead of the default
3. WHERE the BINDORA_DEBUG environment variable is set, THE Application SHALL enable or disable DEBUG mode based on the value
4. WHERE the CORS_ORIGINS environment variable is set, THE Application SHALL configure CORS to accept the specified comma-separated origin list
5. WHERE the VINA_EXE environment variable is set, THE Application SHALL use the specified Vina binary path
6. WHERE the DATABASE_URL environment variable is set, THE Application SHALL connect to the specified database URL instead of the default SQLite path

---

### Requirement 15: Validation Module Parser Pretty Printer

**User Story:** As a developer, I want the Input_Validator module to provide parsers and pretty printers for chemical structures, so that validation results include both parsed structures and formatted output for debugging.

#### Acceptance Criteria

1. THE Input_Validator SHALL provide a parse_smiles function accepting a SMILES string parameter
2. WHEN parse_smiles receives a valid SMILES string, THE parse_smiles function SHALL return an RDKit Mol object
3. WHEN parse_smiles receives an invalid SMILES string, THE parse_smiles function SHALL return None
4. THE Input_Validator SHALL provide a format_smiles function accepting an RDKit Mol object parameter
5. WHEN format_smiles receives a valid Mol object, THE format_smiles function SHALL return the canonical SMILES string
6. THE Input_Validator SHALL provide a parse_pdb function accepting a PDB block string parameter
7. WHEN parse_pdb receives a valid PDB block, THE parse_pdb function SHALL return a gemmi Structure object
8. THE Input_Validator SHALL provide a format_pdb function accepting a gemmi Structure object parameter
9. WHEN format_pdb receives a valid Structure object, THE format_pdb function SHALL return a formatted PDB block string
10. FOR ALL valid SMILES strings, parsing then formatting then parsing SHALL produce an equivalent Mol object (round-trip property)
11. FOR ALL valid PDB blocks with ATOM records, parsing then formatting then parsing SHALL preserve all ATOM records (round-trip property)

---

### Requirement 16: Backward Compatibility Preservation

**User Story:** As a researcher, I want all existing API endpoints and CLI functionality to work after production hardening, so that my scripts and workflows continue functioning without modification.

#### Acceptance Criteria

1. WHEN the production-hardened Application receives requests to existing API endpoints, THE Application SHALL return responses with identical JSON schema as the research prototype
2. THE Application SHALL preserve all existing CLI commands and output formats
3. THE Application SHALL maintain compatibility with existing docked pose file formats (PDBQT)
4. THE Application SHALL preserve existing benchmark report file formats (JSON and Markdown)
5. WHEN existing integration tests are executed against the production-hardened Application, THE tests SHALL pass without modification

---

### Requirement 17: Hardening Validation Test Suite

**User Story:** As a QA engineer, I want automated tests verifying all production hardening features, so that regressions are detected before deployment.

#### Acceptance Criteria

1. THE Application SHALL include a tests/test_hardening.py file testing dependency imports
2. THE test_hardening.py file SHALL include a test verifying DEBUG defaults to False
3. THE test_hardening.py file SHALL include a test verifying MAX_CONTENT_LENGTH rejects oversized payloads
4. THE test_hardening.py file SHALL include a test verifying CORS origin restrictions
5. THE test_hardening.py file SHALL include tests for validate_smiles function with valid and invalid inputs
6. THE test_hardening.py file SHALL include tests for validate_pdb_content function with valid and invalid inputs
7. THE test_hardening.py file SHALL include tests for validate_grid_box function with valid and invalid inputs
8. THE Application SHALL include a tests/test_database.py file testing SQLAlchemy ORM functionality
9. THE test_database.py file SHALL include tests for Session creation and retrieval
10. THE test_database.py file SHALL include tests for batch job persistence
11. THE Application SHALL include a tests/test_docker.py file testing Container build and startup
12. THE test_docker.py file SHALL include a test verifying the Container exposes port 5000 and responds to health checks

---

## Iteration and Feedback

This requirements document is subject to review and refinement. All stakeholders are invited to provide feedback on:

- Completeness of acceptance criteria
- Clarity of requirement specifications
- Technical feasibility and implementation approach
- Security and performance considerations
- Additional requirements for production readiness

All modifications will be tracked and incorporated through iterative review cycles before proceeding to design and implementation phases.
