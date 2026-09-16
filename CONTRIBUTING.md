# Contributing to Bindora Dock

Thank you for your interest in contributing to Bindora Dock.
This document explains how to reproduce benchmarks, add new targets, and run CI checks.

---

## 1. Scientific Integrity Policy

All benchmark output is governed by [`SCIENTIFIC_INTEGRITY.md`](SCIENTIFIC_INTEGRITY.md):

- **Never write benchmark JSON by hand.** All output must go through
  `backend/utils/report_emitter.py` (SHA-256 provenance, atomic writes).
- **Never cherry-pick results.** Every benchmark script reports ALL complexes
  (successes and failures together). Do not add filter logic that hides failures.
- **Always include the DUD-E analogue-bias caveat.** `emit_screening_report()`
  raises `ScientificIntegrityError` if the caveat is missing.
- **Provenance hash.** The `sha256_sig` field is computed over the result payload.
  Editing results by hand breaks the hash.

---

## 2. Environment Setup

```bash
# Python 3.10+ required (tested on 3.11, 3.13)
pip install -r requirements.txt

# AutoDock Vina binary
# Windows: already at bin/vina.exe
# Linux/macOS:
# wget -O bin/vina https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64
# chmod +x bin/vina
```

---

## 3. Reproducing the Existing Benchmarks

Fast check (5 complexes, ~10 min):

```bash
python tests/benchmark_accuracy.py --fast --exhaustiveness 8 \
  --output-json data/benchmarks/my_run.json \
  --output-md   data/benchmarks/my_run.md
```

Full run (25 complexes, ~2-4 hours):

```bash
python tests/benchmark_accuracy.py --exhaustiveness 8 \
  --output-json data/benchmarks/full_run.json
```

---

## 4. Running the CASF-2016 Batch Runner

CASF-2016 (285 complexes) is the field-standard benchmark. The runner is resumable:

```bash
# Smoke test (10 complexes, ~15 min)
python tests/benchmark_casf2016.py --num 10 --exhaustiveness 4

# Full overnight run
python tests/benchmark_casf2016.py --exhaustiveness 8

# Resume after interruption
python tests/benchmark_casf2016.py --resume --exhaustiveness 8

# Use pre-downloaded PDBbind structures (faster)
python tests/benchmark_casf2016.py --pdb-dir path/to/structures/ --resume
```

Outputs: `data/benchmarks/casf2016_final_report.{json,md}`, `casf2016_errors.log`

Reference: Su M. et al. J. Chem. Inf. Model. 2019, 59(2), 895-913.
DOI: 10.1021/acs.jcim.8b00545

---

## 5. Running the DUD-E Screening Runner

DUD-E tests virtual-screening enrichment (ROC-AUC, EF1%, EF5%, EF10%).
SMILES are fetched on-demand from dude.docking.org and cached locally.

> **Mandatory caveat:** DUD-E results always report the analogue-bias disclaimer.
> `emit_screening_report()` raises `ScientificIntegrityError` if it is missing.

```bash
# Diverse 8-target subset
python tests/benchmark_screening.py --subset diverse --exhaustiveness 4

# Single target smoke test
python tests/benchmark_screening.py --target ache --exhaustiveness 4

# Quick test with limited compounds
python tests/benchmark_screening.py --target ache --max-actives 20 --max-decoys 100

# Resume interrupted run
python tests/benchmark_screening.py --subset diverse --resume
```

Reference: Mysinger MM et al. J. Med. Chem. 2012, 55(14), 6582-6594.
DOI: 10.1021/jm300687e

---

## 6. Running CI Locally

The GitHub Actions CI runs on every push to `main`. To replicate locally:

```bash
# Fast benchmark (same as CI)
python tests/benchmark_accuracy.py --fast --exhaustiveness 4 \
  --output-json data/benchmarks/ci_result.json

# Core unit tests
python -m pytest tests/test_docking.py tests/test_adme.py -v
```

The CI passes when at least 3 of 5 fast-benchmark complexes succeed.

---

## 7. Adding a New Benchmark Complex

Add an entry to `BENCHMARK_COMPLEXES` in `tests/benchmark_accuracy.py`:

```python
{
    "pdb_id": "XXXX",
    "target_name": "...",
    "target_class": "Kinase",
    "drug_name": "...",
    "ligand_resname": "LIG",
    "smiles": "...",
    "exp_affinity_type": "IC50",
    "exp_value_nm": 10.5,
    "exp_delta_g_kcal": -9.8,
    # Optional literature comparison fields:
    "literature_vina_rmsd": 1.2,
    "literature_r": 0.60,
    "literature_rmse": 1.8,
    "literature_ref": "DOI:10.xxxx/...",
}
```

Run the benchmark and verify the new complex appears in the full output table (not hidden).
Then open a PR -- CI will run automatically.

---

## 8. Reporting Issues

- **Docking failures on specific PDB IDs:** Open an issue with the PDB ID,
  the error from `casf2016_errors.log`, and your system details.
- **Wrong RMSD / energetics:** Include the JSON checkpoint and commit hash.
- **Security:** Do not open public issues for security vulnerabilities.

---

*Bindora Dock - research-grade molecular docking, openly benchmarked.*
