import os
import sys
import json
import datetime
import hashlib
import rdkit
import meeko
import numpy
import scipy
import sklearn

def get_file_sha256(path):
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def main():
    vina_path = "bin/vina.exe"
    vina_hash = get_file_sha256(vina_path)
    
    frozen_hashes_path = "benchmarks/heldout/FROZEN_HASHES.txt"
    frozen_hash = get_file_sha256(frozen_hashes_path)
    
    manifest_data = {
        "benchmark_audit_metadata": {
            "title": "Bindora Dock Independent Scientific Validation Audit",
            "audit_type": "Observation and Read-Only Empirical Held-Out Benchmarking",
            "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "random_seed": 42
        },
        "environment": {
            "python_version": sys.version,
            "os_platform": sys.platform,
            "packages": {
                "rdkit": rdkit.__version__,
                "meeko": meeko.__version__,
                "numpy": numpy.__version__,
                "scipy": scipy.__version__,
                "scikit-learn": sklearn.__version__
            },
            "binaries": {
                "vina_path": vina_path,
                "vina_sha256": vina_hash
            }
        },
        "frozen_datasets": {
            "frozen_hashes_file": frozen_hashes_path,
            "frozen_hashes_sha256": frozen_hash,
            "datasets": {
                "dev_set": {
                    "path": "benchmarks/heldout/dev_set.csv",
                    "sha256": get_file_sha256("benchmarks/heldout/dev_set.csv"),
                    "count": 42
                },
                "pgp_train": {
                    "path": "benchmarks/heldout/train.csv",
                    "sha256": get_file_sha256("benchmarks/heldout/train.csv"),
                    "count": 975
                },
                "pgp_test": {
                    "path": "benchmarks/heldout/test.csv",
                    "sha256": get_file_sha256("benchmarks/heldout/test.csv"),
                    "count": 244
                },
                "bbb_test": {
                    "path": "benchmarks/heldout/bbb_test.csv",
                    "sha256": get_file_sha256("benchmarks/heldout/bbb_test.csv"),
                    "count": 7782
                },
                "pains_dataset": {
                    "path": "benchmarks/heldout/pains_dataset.csv",
                    "sha256": get_file_sha256("benchmarks/heldout/pains_dataset.csv"),
                    "count": 67
                }
            }
        },
        "results_artifacts": {
            "macrocycle_eval_csv": "benchmarks/heldout/results/macrocycle_eval.csv",
            "pains_benchmark_csv": "benchmarks/heldout/results/pains_benchmark.csv",
            "pains_metrics_summary": "benchmarks/heldout/results/pains_metrics_summary.json",
            "pgp_benchmark_csv": "benchmarks/heldout/results/pgp_benchmark.csv",
            "pgp_metrics_summary": "benchmarks/heldout/results/pgp_metrics_summary.json",
            "bbb_benchmark_csv": "benchmarks/heldout/results/bbb_benchmark.csv",
            "bbb_metrics_summary": "benchmarks/heldout/results/bbb_metrics_summary.json"
        }
    }
    
    out_file = "benchmarks/heldout/manifest.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(manifest_data, f, indent=2)
        
    print(f"Manifest successfully created at {out_file}")

if __name__ == "__main__":
    main()
