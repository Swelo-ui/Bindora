# Bindora Scientific Integrity & Programmatic Data Standard

This document establishes the mandatory protocol for data persistence, benchmarking, scoring evaluations, and computational reporting within the **Bindora** platform.

---

## The Core Rule: Zero Manual Transcription

> **MANDATORY PRINCIPLE:**  
> **Under no circumstances should benchmark results, rescoring metrics, screening tables, or validation JSON files ever be hand-typed or written via manual dictionary scripts.**  
> Every report file must be emitted directly from active execution memory by the underlying calculation engine or test runner.

---

## 1. Architectural Architecture: `backend.utils.report_emitter`

All pipelines—current and future—must route result persistence through [`backend/utils/report_emitter.py`](backend/utils/report_emitter.py).

### Usage in Redocking & Benchmarks
```python
from backend.utils.report_emitter import emit_benchmark_record

# Extract live values directly from in-memory calculation objects
live_result = {
    "target": benchmark["target"],
    "ligand": benchmark["drug"],
    "vina_affinity_kcal": round(float(best_pose["affinity_kcal"]), 2),
    "vinardo_affinity_kcal": best_pose.get("vinardo_affinity_kcal"),
    "mode1_rmsd_angstroms": round(float(actual_rmsd), 2),
    "docking_time_seconds": round(time.time() - t_start, 2),
    "exhaustiveness": benchmark["exhaustiveness"],
    "overall_grade": "RESEARCH_GRADE"
}

# Programmatically validate types and emit atomically to data/benchmarks/
emit_benchmark_record(benchmark["pdb_id"], live_result)
```

### Usage in Rescoring (MM-GBSA, GNINA, Interaction Profiles)
```python
from backend.utils.report_emitter import emit_pipeline_report

summary_data = {
    "complex_id": "1AQ1",
    "delta_g_mmgbsa_kcal": -42.85,
    "coulombic_energy": -18.2,
    "vdW_energy": -34.6,
    "solvation_polar": 21.3,
    "solvation_apolar": -11.35,
    "sampling_time_seconds": 12.4
}

emit_pipeline_report(
    pipeline_name="MM-GBSA Rescoring Engine",
    summary_data=summary_data,
    output_file="data/reports/mmgbsa_1aq1.json"
)
```

---

## 2. Integrity Protections Enforced by `report_emitter`

1. **Biophysical Type & Bounds Checking:**
   - Numerical values (`vina_affinity_kcal`, `rmsd_angstroms`, etc.) are checked with `validate_numeric_field`.
   - Rejects `NaN`, infinite values, or strings mimicking numbers.
   - Enforces physical energy boundaries (e.g. affinities between $-30.0$ and $+20.0\text{ kcal/mol}$).

2. **Cryptographic Run Provenance:**
   - Every emitted file automatically receives a SHA-256 calculation signature (`provenance.sha256_sig`) calculated from the raw payload parameters.
   - Includes Python version, OS platform, and ISO 8601 UTC timestamps.

3. **Atomic Filesystem Writes:**
   - Files are written to temporary staging files first and then replaced via `os.replace`.
   - Prevents partial or corrupted JSON files during concurrent API reads or abrupt process termination.

---

## 3. Scope of Application

This standard applies strictly to:
- [x] **Astex Diverse Set & CASF Redocking Benchmarks** (`tests/test_research_grade.py`, `tests/benchmark_accuracy.py`)
- [x] **Virtual Screening Pipelines** (Multi-compound docking batches and ranking outputs)
- [x] **Advanced Rescoring Modules** (Vinardo empirical scoring, GNINA CNN affinity, MM-GBSA / MM-PBSA free energy calculations)
- [x] **ADMET & Bioactivity Summaries**
- [x] **API Cache & Export Endpoints**

By enforcing this discipline at the foundation, Bindora ensures that every benchmark figure, publication claim, and user-facing result is 100% reproducible and tamper-free.

