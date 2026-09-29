import hashlib
import os
import datetime

FROZEN_FILES = [
    "benchmarks/heldout/dev_set.csv",
    "benchmarks/heldout/train.csv",
    "benchmarks/heldout/test.csv",
    "benchmarks/heldout/bbb_test.csv",
    "benchmarks/heldout/pains_dataset.csv"
]

def get_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def main():
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    lines = [
        f"# FROZEN BENCHMARK DATASET HASHES",
        f"# Generated: {timestamp}",
        f"# Rule: Once written, no modification to these files is permitted.",
        ""
    ]
    
    for fpath in FROZEN_FILES:
        if not os.path.exists(fpath):
            raise FileNotFoundError(f"Missing required file to freeze: {fpath}")
        size = os.path.getsize(fpath)
        sha = get_sha256(fpath)
        lines.append(f"{sha}  {fpath}  ({size} bytes)")
        print(f"{sha}  {fpath}  ({size} bytes)")
        
    out_file = "benchmarks/heldout/FROZEN_HASHES.txt"
    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        
    print(f"\nSaved frozen hashes to {out_file}")

if __name__ == "__main__":
    main()
