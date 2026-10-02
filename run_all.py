"""
run_all.py — Bindora Dock Research-Grade Benchmark Suite
Executes all held-out scientific benchmarks in sequence and produces a manifest.

Usage:
    python run_all.py [--skip-docking] [--skip-slow]

Stages:
    E1: PAINS benchmark (65 compounds, Baell & Holloway 2010)
    E2: P-gp substrate benchmark (66 held-out, Wang 2011)
    E3: BBB benchmark (7782 compounds, B3DB)
    E4: Decoy gating benchmark (9 compounds, EGFR/1M17)

Outputs written to: benchmarks/heldout/results/
Manifest written to: benchmarks/heldout/manifest.json
"""
import subprocess
import sys
import json
import hashlib
import argparse
from pathlib import Path
from datetime import datetime, timezone

REPO_ROOT = Path(__file__).resolve().parent
RESULTS_DIR = REPO_ROOT / "benchmarks" / "heldout" / "results"
SCRATCH_DIR = REPO_ROOT / "benchmarks" / "heldout" / "scratch"
MANIFEST_OUT = REPO_ROOT / "benchmarks" / "heldout" / "manifest.json"

SCRIPTS = [
    {
        "id": "E1",
        "name": "PAINS Benchmark (Baell & Holloway 2010, 65 compounds)",
        "script": "e1_pains_benchmark.py",
        "outputs": ["E1_pains_results.csv", "E1_pains_summary.json"],
        "slow": False,
        "requires_server": False,
    },
    {
        "id": "E2",
        "name": "P-gp Substrate Benchmark (Wang 2011, 66 held-out compounds)",
        "script": "e2_pgp_benchmark.py",
        "outputs": ["E2_pgp_results.csv", "E2_pgp_summary.json"],
        "slow": False,
        "requires_server": False,
    },
    {
        "id": "E3",
        "name": "BBB Permeability Benchmark (B3DB, 7782 compounds)",
        "script": "e3_bbb_benchmark.py",
        "outputs": ["E3_bbb_results.csv", "E3_bbb_summary.json"],
        "slow": False,
        "requires_server": False,
    },
    {
        "id": "E4",
        "name": "Decoy Gating Benchmark (live API, 9 compounds, EGFR/1M17)",
        "script": "e4_decoy_gating_benchmark.py",
        "outputs": ["E4_decoy_results.csv", "E4_decoy_summary.json"],
        "slow": True,
        "requires_server": True,
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def run_script(script_name: str) -> bool:
    script_path = SCRATCH_DIR / script_name
    if not script_path.exists():
        print(f"  ERROR: Script not found: {script_path}")
        return False
    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(REPO_ROOT),
        capture_output=False,
    )
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description="Bindora Benchmark Suite")
    parser.add_argument("--skip-docking", action="store_true", help="Skip live-server docking tests (E4)")
    parser.add_argument("--skip-slow", action="store_true", help="Skip slow benchmarks")
    args = parser.parse_args()

    print("=" * 70)
    print("  Bindora Dock — Research-Grade Benchmark Suite")
    print(f"  Started: {datetime.now(timezone.utc).isoformat()}")
    print("=" * 70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    manifest_entries = []

    for stage in SCRIPTS:
        sid = stage["id"]
        name = stage["name"]
        script = stage["script"]

        if args.skip_docking and stage.get("requires_server"):
            print(f"\n[{sid}] SKIPPED (--skip-docking): {name}")
            continue
        if args.skip_slow and stage.get("slow"):
            print(f"\n[{sid}] SKIPPED (--skip-slow): {name}")
            continue

        print(f"\n[{sid}] Running: {name}")
        print(f"  Script: {script}")

        ok = run_script(script)
        status = "SUCCESS" if ok else "FAILED"
        print(f"  Status: {status}")

        output_files = []
        for fname in stage["outputs"]:
            fpath = RESULTS_DIR / fname
            if fpath.exists():
                output_files.append({
                    "filename": fname,
                    "path": str(fpath.relative_to(REPO_ROOT)),
                    "size_bytes": fpath.stat().st_size,
                    "sha256": sha256(fpath),
                })
            else:
                output_files.append({"filename": fname, "path": None, "error": "Not generated"})

        manifest_entries.append({
            "stage_id": sid,
            "stage_name": name,
            "script": script,
            "status": status,
            "output_files": output_files,
        })

    # Write manifest
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "repo_root": str(REPO_ROOT),
        "python_version": sys.version,
        "stages": manifest_entries,
    }
    with open(MANIFEST_OUT, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"\n{'=' * 70}")
    print(f"  Manifest: {MANIFEST_OUT}")
    n_pass = sum(1 for e in manifest_entries if e["status"] == "SUCCESS")
    print(f"  Stages passed: {n_pass}/{len(manifest_entries)}")
    print(f"  Finished: {datetime.now(timezone.utc).isoformat()}")
    print("=" * 70)


if __name__ == "__main__":
    main()
