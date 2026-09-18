"""
Test suite for database functionality
Tests SQLAlchemy ORM, session management, and CRUD operations
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest
import tempfile
import os
from datetime import datetime
from backend.db.database import init_db, get_db, Base
from backend.db.models import DockingSession, BatchScreeningJob, BenchmarkRecord


@pytest.fixture
def test_db():
    """Create a temporary test database"""
    # Create temporary database file
    db_fd, db_path = tempfile.mkstemp(suffix='.db')
    database_url = f"sqlite:///{db_path}"
    
    # Initialize database
    init_db(database_url)
    
    yield database_url
    
    # Cleanup
    from backend.db.database import close_db
    close_db()
    os.close(db_fd)
    try:
        os.unlink(db_path)
    except Exception:
        pass


class TestDatabaseInitialization:
    """Test database setup and initialization"""
    
    def test_init_db_creates_tables(self, test_db):
        """Test that init_db creates all tables"""
        with get_db() as db:
            # Check that tables exist by querying them
            result = db.query(DockingSession).count()
            assert result == 0  # Empty table
    
    def test_get_db_context_manager(self, test_db):
        """Test get_db context manager"""
        with get_db() as db:
            assert db is not None
            # DB should auto-commit on context exit


class TestDockingSessionModel:
    """Test DockingSession ORM model"""
    
    def test_create_docking_session(self, test_db):
        """Test creating a docking session record"""
        with get_db() as db:
            session = DockingSession(
                pdb_id="1CX2",
                ligand_name="Aspirin",
                ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
                affinity_kcal=-7.5,
                execution_duration_s=12.5,
                exhaustiveness=8,
                num_modes=9
            )
            db.add(session)
            db.commit()
            
            assert session.id is not None
            assert session.timestamp is not None
    
    def test_query_docking_sessions(self, test_db):
        """Test querying docking sessions"""
        with get_db() as db:
            # Create test sessions
            session1 = DockingSession(
                pdb_id="1CX2",
                ligand_name="Aspirin",
                ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
                affinity_kcal=-7.5,
                execution_duration_s=12.5
            )
            session2 = DockingSession(
                pdb_id="1HSG",
                ligand_name="Indinavir",
                ligand_smiles="CC(C)(C)NC(=O)C1CC2CCCCC2CN1CC(O)C(Cc1ccccc1)NC(=O)C(O)C(N)Cc1ccccc1",
                affinity_kcal=-9.2,
                execution_duration_s=15.3
            )
            db.add_all([session1, session2])
            db.commit()
            
            # Query all sessions
            sessions = db.query(DockingSession).all()
            assert len(sessions) == 2
            
            # Query by PDB ID
            cx2_sessions = db.query(DockingSession).filter_by(pdb_id="1CX2").all()
            assert len(cx2_sessions) == 1
            assert cx2_sessions[0].ligand_name == "Aspirin"
    
    def test_docking_session_to_dict(self, test_db):
        """Test DockingSession to_dict method"""
        with get_db() as db:
            session = DockingSession(
                pdb_id="1CX2",
                ligand_name="Aspirin",
                ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
                affinity_kcal=-7.5,
                execution_duration_s=12.5
            )
            db.add(session)
            db.commit()
            
            session_dict = session.to_dict()
            assert isinstance(session_dict, dict)
            assert session_dict['pdb_id'] == "1CX2"
            assert session_dict['ligand_name'] == "Aspirin"
            assert session_dict['affinity_kcal'] == -7.5
    
    def test_delete_docking_session(self, test_db):
        """Test deleting a docking session"""
        with get_db() as db:
            session = DockingSession(
                pdb_id="1CX2",
                ligand_name="Aspirin",
                ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
                affinity_kcal=-7.5,
                execution_duration_s=12.5
            )
            db.add(session)
            db.commit()
            session_id = session.id
            
            # Delete session
            db.delete(session)
            db.commit()
            
            # Verify deletion
            deleted_session = db.query(DockingSession).filter_by(id=session_id).first()
            assert deleted_session is None


class TestBatchScreeningJobModel:
    """Test BatchScreeningJob ORM model"""
    
    def test_create_batch_job(self, test_db):
        """Test creating a batch screening job"""
        with get_db() as db:
            job = BatchScreeningJob(
                job_id="batch_001",
                receptor_pdb_id="1CX2",
                total_ligands=100,
                completed_ligands=0,
                status="queued"
            )
            db.add(job)
            db.commit()
            
            assert job.id is not None
            assert job.created_at is not None
    
    def test_update_batch_job_status(self, test_db):
        """Test updating batch job status"""
        with get_db() as db:
            job = BatchScreeningJob(
                job_id="batch_002",
                receptor_pdb_id="1HSG",
                total_ligands=50,
                status="queued"
            )
            db.add(job)
            db.commit()
            
            # Update status
            job.status = "running"
            job.completed_ligands = 25
            db.commit()
            
            # Verify update
            updated_job = db.query(BatchScreeningJob).filter_by(job_id="batch_002").first()
            assert updated_job.status == "running"
            assert updated_job.completed_ligands == 25


class TestBenchmarkRecordModel:
    """Test BenchmarkRecord ORM model"""
    
    def test_create_benchmark_record(self, test_db):
        """Test creating a benchmark record"""
        with get_db() as db:
            record = BenchmarkRecord(
                benchmark_type="casf2016",
                pdb_id="1A1E",
                ligand_name="Native",
                rmsd_angstrom=1.25,
                affinity_kcal=-8.5
            )
            db.add(record)
            db.commit()
            
            assert record.id is not None
    
    def test_query_benchmark_by_type(self, test_db):
        """Test querying benchmarks by type"""
        with get_db() as db:
            # Create multiple benchmark records
            record1 = BenchmarkRecord(
                benchmark_type="casf2016",
                pdb_id="1A1E",
                ligand_name="Native",
                rmsd_angstrom=1.25,
                affinity_kcal=-8.5
            )
            record2 = BenchmarkRecord(
                benchmark_type="dude",
                pdb_id="ACHE",
                ligand_name="Active1",
                rmsd_angstrom=0.0,
                affinity_kcal=-9.2
            )
            db.add_all([record1, record2])
            db.commit()
            
            # Query CASF records
            casf_records = db.query(BenchmarkRecord).filter_by(benchmark_type="casf2016").all()
            assert len(casf_records) == 1
            assert casf_records[0].pdb_id == "1A1E"


class TestTransactionHandling:
    """Test database transaction rollback"""
    
    def test_rollback_on_error(self, test_db):
        """Test that transactions rollback on error"""
        try:
            with get_db() as db:
                session = DockingSession(
                    pdb_id="TEST",
                    ligand_name="Test",
                    ligand_smiles="C",
                    affinity_kcal=-5.0,
                    execution_duration_s=10.0
                )
                db.add(session)
                # Trigger an error to be handled by get_db() context manager
                raise ValueError("Test error")
        except ValueError:
            # Context manager should have rolled back
            pass
        
        # Verify no record was created
        with get_db() as db:
            sessions = db.query(DockingSession).filter_by(pdb_id="TEST").all()
            assert len(sessions) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
