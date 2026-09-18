"""
Database initialization and session management
SQLAlchemy 2.0 with SQLite backend
"""

from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker, declarative_base
from contextlib import contextmanager
import os

# Global state
engine = None
SessionLocal = None
Base = declarative_base()


def init_db(database_url: str):
    """
    Initialize database engine and create tables
    
    Args:
        database_url: SQLAlchemy database URL (e.g., sqlite:///data/bindora.db)
    """
    global engine, SessionLocal
    
    # Create engine with SQLite optimizations
    engine = create_engine(
        database_url,
        echo=False,
        connect_args={"check_same_thread": False} if database_url.startswith("sqlite") else {},
        pool_pre_ping=True  # Verify connections before using
    )
    
    # Create scoped session factory for thread-safe access
    SessionLocal = scoped_session(
        sessionmaker(autocommit=False, autoflush=False, bind=engine)
    )
    
    # Create all tables
    Base.metadata.create_all(bind=engine)
    
    print(f"[DATABASE] Initialized at {database_url}")


@contextmanager
def get_db():
    """
    Get thread-safe database session context manager
    
    Usage:
        with get_db() as db:
            session = db.query(DockingSession).all()
    
    Yields:
        SQLAlchemy Session object
    """
    if SessionLocal is None:
        raise RuntimeError("Database not initialized. Call init_db() first.")
    
    db = SessionLocal()
    try:
        yield db
        db.commit()
    except Exception as e:
        db.rollback()
        raise e
    finally:
        db.close()


def get_db_session():
    """
    Get database session (non-context manager version)
    
    Returns:
        SQLAlchemy Session object
        
    Note: Caller must close the session manually
    """
    if SessionLocal is None:
        raise RuntimeError("Database not initialized. Call init_db() first.")
    
    return SessionLocal()


def close_db():
    """Close and dispose database connections and session pool."""
    global engine, SessionLocal
    if SessionLocal is not None:
        SessionLocal.remove()
    if engine is not None:
        engine.dispose()

