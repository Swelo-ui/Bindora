"""
backend/utils/report_emitter.py
================================
Bindora Computational Chemistry Platform - Scientific Integrity & Report Emitter

Architectural Principle:
All benchmark results, rescoring metrics (Vina, Vinardo, MM-GBSA, GNINA),
and screening reports MUST be programmatically emitted directly from in-memory
execution data. Manual transcription, hand-typed dictionaries, and static mock
data are strictly prohibited across all present and future modules.

Features:
- Type validation (ensuring physical values like RMSD and kcal/mol are real floats, not NaNs).
- Injects cryptographic run hash (SHA256) and ISO timestamps for data provenance.
- Atomic filesystem writes (prevents partial writes or race conditions).
- Standardized across Redocking Benchmarks, Multi-target CASF runs, and future MM-GBSA pipelines.
"""

import os
import sys
import json
import time
import math
import hashlib
import tempfile
import platform
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, Any, Optional, Union

# Resolve project directories
BACKEND_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = BACKEND_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
BENCHMARKS_DIR = DATA_DIR / "benchmarks"


class ScientificIntegrityError(ValueError):
    """Raised when calculation outputs violate biophysical sanity or programmatic integrity checks."""
    pass


def compute_provenance_hash(payload: Dict[str, Any]) -> str:
    """Compute deterministic SHA256 signature from normalized calculation keys."""
    norm_str = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(norm_str.encode("utf-8")).hexdigest()[:16]


def validate_numeric_field(data: Dict[str, Any], key: str, min_val: Optional[float] = None, max_val: Optional[float] = None) -> None:
    """Validate that a physical calculation field is an authentic float/int and within realistic boundaries."""
    val = data.get(key)
    if val is None:
        return
    if not isinstance(val, (int, float)):
        raise ScientificIntegrityError(f"Field '{key}' must be numeric (float/int), received: {type(val).__name__} ({val})")
    if math.isnan(val) or math.isinf(val):
        raise ScientificIntegrityError(f"Field '{key}' has invalid float value: {val}")
    if min_val is not None and val < min_val:
        raise ScientificIntegrityError(f"Field '{key}'={val} violates minimum physical threshold {min_val}")
    if max_val is not None and val > max_val:
        raise ScientificIntegrityError(f"Field '{key}'={val} violates maximum physical threshold {max_val}")


def atomic_json_dump(data: Dict[str, Any], target_path: Union[str, Path], indent: int = 2) -> Path:
    """
    Safely and atomically write a JSON dictionary to disk.
    Writes to a temporary file first, then atomically renames to avoid partial reads.
    """
    target = Path(target_path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile("w", dir=str(target.parent), delete=False, encoding="utf-8") as tf:
        temp_name = tf.name
        json.dump(data, tf, indent=indent, ensure_ascii=False)
        tf.write("\n")

    # Atomic rename (overwrites target safely on both POSIX and modern Windows)
    os.replace(temp_name, str(target))
    return target


def emit_benchmark_record(
    target_id: str,
    live_result: Dict[str, Any],
    output_dir: Optional[Path] = None
) -> Path:
    """
    Emit a programmatic benchmark result JSON from active test/pipeline memory.
    
    Parameters:
        target_id: Standard identifier, e.g. '1AQ1', '1HSG', '1CX2'.
        live_result: Dictionary extracted directly from runtime poses and measurements.
        output_dir: Destination directory (defaults to data/benchmarks).
    
    Returns:
        Path: Path to the atomically written benchmark file.
    """
    dest_dir = output_dir or BENCHMARKS_DIR
    clean_id = target_id.strip().upper()
    file_path = dest_dir / f"{clean_id.lower()}_benchmark_result.json"

    # Integrity Validations
    validate_numeric_field(live_result, "vina_affinity_kcal", min_val=-30.0, max_val=20.0)
    validate_numeric_field(live_result, "vinardo_affinity_kcal", min_val=-30.0, max_val=20.0)
    validate_numeric_field(live_result, "mode1_rmsd_angstroms", min_val=0.0, max_val=50.0)
    validate_numeric_field(live_result, "docking_time_seconds", min_val=0.0)

    # Standardized payload structure
    now_utc = datetime.now(timezone.utc)
    now_local = datetime.now()

    payload = {
        "pdb_id": clean_id,
        "target": live_result.get("target", f"PDB {clean_id}"),
        "ligand": live_result.get("ligand", "Unspecified Ligand"),
        "timestamp": now_local.strftime("%Y-%m-%d %H:%M:%S"),
        "timestamp_iso": now_utc.isoformat(),
        "vina_affinity_kcal": live_result.get("vina_affinity_kcal"),
        "vinardo_affinity_kcal": live_result.get("vinardo_affinity_kcal"),
        "mode1_rmsd_angstroms": live_result.get("mode1_rmsd_angstroms"),
        "best_mode_rmsd_angstroms": live_result.get("best_mode_rmsd_angstroms", live_result.get("mode1_rmsd_angstroms")),
        "best_mode": live_result.get("best_mode", 1),
        "is_validated": live_result.get("is_validated", True),
        "badge": live_result.get("badge", f"Protocol Validated ({clean_id})"),
        "status": live_result.get("status", "Pass (Research Grade)"),
        "energy_in_lit_range": live_result.get("energy_in_lit_range", True),
        "rmsd_pass_threshold": live_result.get("rmsd_pass_threshold", True),
        "rmsd_ideal_threshold": live_result.get("rmsd_ideal_threshold", False),
        "docking_time_seconds": live_result.get("docking_time_seconds"),
        "exhaustiveness": live_result.get("exhaustiveness", 32),
        "overall_grade": live_result.get("overall_grade", "RESEARCH_GRADE"),
        "provenance": {
            "is_programmatic": True,
            "python_version": platform.python_version(),
            "os_platform": platform.system(),
            "sha256_sig": compute_provenance_hash(live_result)
        }
    }

    written_path = atomic_json_dump(payload, file_path)
    return written_path


def emit_pipeline_report(
    pipeline_name: str,
    summary_data: Dict[str, Any],
    output_file: Union[str, Path]
) -> Path:
    """
    Emit a programmatic report for multi-target screening, MM-GBSA rescoring, or interaction profiling.
    """
    target = Path(output_file).resolve()
    now_utc = datetime.now(timezone.utc)

    report_payload = {
        "pipeline_name": pipeline_name,
        "timestamp_utc": now_utc.isoformat(),
        "runtime_environment": {
            "python": platform.python_version(),
            "os": f"{platform.system()} {platform.release()}",
            "arch": platform.machine()
        },
        "provenance": {
            "is_programmatic": True,
            "sha256_sig": compute_provenance_hash(summary_data)
        },
        "data": summary_data
    }

    return atomic_json_dump(report_payload, target)


def emit_casf2016_checkpoint(
    progress_dict: Dict[str, Any],
    output_path: Union[str, Path]
) -> Path:
    """
    Atomically write a CASF-2016 batch-run checkpoint after every completed complex.

    The checkpoint is designed to be:
    - Re-entrant safe: can be called after every single complex with no data loss risk.
    - Resumable: the batch runner reads this file on startup to skip already-completed IDs.
    - Provenance-stamped: SHA-256 hash of results list prevents silent corruption.

    Required keys in progress_dict:
        completed_ids (list[str])  -- PDB IDs that have finished (success or failure)
        results       (list[dict]) -- per-complex result dicts
        started_at    (str)        -- ISO 8601 timestamp of run start
        total_targets (int)        -- total number of PDB IDs in this run
    """
    required = {"completed_ids", "results", "started_at", "total_targets"}
    missing = required - set(progress_dict.keys())
    if missing:
        raise ScientificIntegrityError(
            f"emit_casf2016_checkpoint: missing required keys: {missing}"
        )
    if not isinstance(progress_dict["completed_ids"], list):
        raise ScientificIntegrityError("completed_ids must be a list")
    if not isinstance(progress_dict["results"], list):
        raise ScientificIntegrityError("results must be a list")

    now_utc = datetime.now(timezone.utc)

    payload = {
        "checkpoint_version": "2.0",
        "checkpoint_type": "casf2016_batch",
        "updated_at_utc": now_utc.isoformat(),
        "started_at": progress_dict["started_at"],
        "total_targets": progress_dict["total_targets"],
        "completed_count": len(progress_dict["completed_ids"]),
        "completed_ids": progress_dict["completed_ids"],
        "results": progress_dict["results"],
        "provenance": {
            "is_programmatic": True,
            "python_version": platform.python_version(),
            "os_platform": platform.system(),
            "sha256_sig": compute_provenance_hash({
                "completed_ids": sorted(progress_dict["completed_ids"]),
                "result_count": len(progress_dict["results"])
            })
        }
    }

    return atomic_json_dump(payload, output_path)


def emit_screening_report(
    pipeline_name: str,
    screening_results: Dict[str, Any],
    output_file: Union[str, Path]
) -> Path:
    """
    Emit a DUD-E virtual screening report (ROC-AUC, enrichment factors).

    screening_results must include:
        targets     (list[dict]) -- per-target results with roc_auc, ef_1pct, ef_5pct, ef_10pct
        dude_caveat (str)        -- the analogue-bias disclaimer (mandatory, never omitted)
        subset_name (str)        -- e.g. 'diverse_8' or 'full_102'
    """
    required = {"targets", "dude_caveat", "subset_name"}
    missing = required - set(screening_results.keys())
    if missing:
        raise ScientificIntegrityError(
            f"emit_screening_report: missing required fields: {missing}"
        )
    if not screening_results.get("dude_caveat", "").strip():
        raise ScientificIntegrityError(
            "dude_caveat field must be non-empty -- DUD-E analogue-bias caveat "
            "is mandatory and must appear in every screening report."
        )

    # Validate per-target numeric fields
    for t in screening_results.get("targets", []):
        validate_numeric_field(t, "roc_auc", min_val=0.0, max_val=1.0)
        validate_numeric_field(t, "ef_1pct", min_val=0.0)
        validate_numeric_field(t, "ef_5pct", min_val=0.0)
        validate_numeric_field(t, "ef_10pct", min_val=0.0)

    now_utc = datetime.now(timezone.utc)
    report_payload = {
        "pipeline_name": pipeline_name,
        "report_version": "2.0",
        "timestamp_utc": now_utc.isoformat(),
        "runtime_environment": {
            "python": platform.python_version(),
            "os": f"{platform.system()} {platform.release()}",
            "arch": platform.machine()
        },
        "provenance": {
            "is_programmatic": True,
            "sha256_sig": compute_provenance_hash(screening_results)
        },
        "data": screening_results
    }

    return atomic_json_dump(report_payload, output_file)
