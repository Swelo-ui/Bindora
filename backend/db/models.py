"""
ORM Models for Bindora Dock database
"""

from sqlalchemy import Column, String, Float, Integer, DateTime, Text
from sqlalchemy.dialects.postgresql import UUID as PGUUID
from backend.db.database import Base
import uuid
from datetime import datetime


class DockingSession(Base):
    """Record of a completed molecular docking simulation"""
    __tablename__ = "docking_sessions"
    
    id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    timestamp = Column(DateTime, default=datetime.utcnow, index=True, nullable=False)
    pdb_id = Column(String(10), index=True)
    ligand_name = Column(String(200))
    ligand_smiles = Column(Text)
    affinity_kcal = Column(Float)
    rmsd_lb = Column(Float, nullable=True)
    rmsd_ub = Column(Float, nullable=True)
    execution_duration_s = Column(Float)
    exhaustiveness = Column(Integer, default=8)
    num_modes = Column(Integer, default=9)
    receptor_pdbqt_path = Column(String(500), nullable=True)
    ligand_pdbqt_path = Column(String(500), nullable=True)
    top_pose_pdbqt_path = Column(String(500), nullable=True)
    
    def to_dict(self):
        """Convert to dictionary for JSON serialization"""
        return {
            'id': self.id,
            'timestamp': self.timestamp.isoformat() if self.timestamp else None,
            'pdb_id': self.pdb_id,
            'ligand_name': self.ligand_name,
            'ligand_smiles': self.ligand_smiles,
            'affinity_kcal': self.affinity_kcal,
            'rmsd_lb': self.rmsd_lb,
            'rmsd_ub': self.rmsd_ub,
            'execution_duration_s': self.execution_duration_s,
            'exhaustiveness': self.exhaustiveness,
            'num_modes': self.num_modes
        }


class BatchScreeningJob(Base):
    """Batch virtual screening job tracking"""
    __tablename__ = "batch_screening_jobs"
    
    id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    job_id = Column(String(100), unique=True, index=True, nullable=False)
    receptor_pdb_id = Column(String(10))
    total_ligands = Column(Integer, nullable=False)
    completed_ligands = Column(Integer, default=0)
    status = Column(String(20), default="queued", nullable=False)  # queued, running, completed, cancelled, interrupted
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    completed_at = Column(DateTime, nullable=True)
    top_hit_name = Column(String(200), nullable=True)
    top_affinity_kcal = Column(Float, nullable=True)
    
    def to_dict(self):
        """Convert to dictionary for JSON serialization"""
        return {
            'id': self.id,
            'job_id': self.job_id,
            'receptor_pdb_id': self.receptor_pdb_id,
            'total_ligands': self.total_ligands,
            'completed_ligands': self.completed_ligands,
            'status': self.status,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'completed_at': self.completed_at.isoformat() if self.completed_at else None,
            'top_hit_name': self.top_hit_name,
            'top_affinity_kcal': self.top_affinity_kcal
        }


class BenchmarkRecord(Base):
    """Benchmark validation results"""
    __tablename__ = "benchmark_records"
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    benchmark_type = Column(String(50), index=True, nullable=False)  # casf2016, dude, accuracy
    pdb_id = Column(String(10), nullable=False)
    ligand_name = Column(String(200))
    rmsd_angstrom = Column(Float)
    affinity_kcal = Column(Float)
    vinardo_affinity_kcal = Column(Float, nullable=True)
    timestamp = Column(DateTime, default=datetime.utcnow, index=True, nullable=False)
    
    def to_dict(self):
        """Convert to dictionary for JSON serialization"""
        return {
            'id': self.id,
            'benchmark_type': self.benchmark_type,
            'pdb_id': self.pdb_id,
            'ligand_name': self.ligand_name,
            'rmsd_angstrom': self.rmsd_angstrom,
            'affinity_kcal': self.affinity_kcal,
            'vinardo_affinity_kcal': self.vinardo_affinity_kcal,
            'timestamp': self.timestamp.isoformat() if self.timestamp else None
        }
