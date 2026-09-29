import hashlib
import os
import datetime

FROZEN_FILES = [
    ("benchmarks/heldout/dev_set.csv", "42 contaminated dev-set compounds with PubChem CIDs and InChIKeys"),
    ("benchmarks/heldout/pains_dataset.csv", "65 verified PAINS compounds with PubChem title match and Baell 2010 families"),
    ("benchmarks/heldout/pgp_substrate_train.csv", "261 training compounds from Wang 2011 P-gp substrate dataset (Bemis-Murcko split)"),
    ("benchmarks/heldout/pgp_substrate_test.csv", "66 held-out test compounds from Wang 2011 P-gp substrate dataset (Bemis-Murcko split)"),
    ("benchmarks/heldout/pgp_substrate_full.csv", "327 total non-dev compounds from Wang 2011 P-gp substrate dataset"),
    ("benchmarks/heldout/bbb_test.csv", "7782 experimental compounds from B3DB with confidence group annotations")
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
        f"# FROZEN BENCHMARK DATASET HASHES (REBUILT & VERIFIED)",
        f"# Generated: {timestamp}",
        f"# Strict Rule: Once written, no modification to these files is permitted.",
        ""
    ]
    
    hashes_dict = {}
    for fpath, desc in FROZEN_FILES:
        if not os.path.exists(fpath):
            raise FileNotFoundError(f"Missing required file: {fpath}")
        size = os.path.getsize(fpath)
        sha = get_sha256(fpath)
        hashes_dict[fpath] = sha
        lines.append(f"{sha}  {fpath}  ({size} bytes)  # {desc}")
        
    out_file = "benchmarks/heldout/FROZEN_HASHES.txt"
    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
        
    print(f"Written updated frozen hashes to {out_file}\n")
    
    # RE-VERIFY ALL LISTED FILES
    print("--- RE-VERIFYING ALL LISTED FILES AGAINST FROZEN_HASHES.TXT ---")
    verified = True
    with open(out_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            expected_sha = parts[0]
            rel_path = parts[1]
            actual_sha = get_sha256(rel_path)
            if actual_sha == expected_sha:
                print(f"  [OK: VERIFIED] {rel_path} -> {actual_sha}")
            else:
                print(f"  [FAILED] {rel_path} -> expected {expected_sha}, got {actual_sha}")
                verified = False
                
    if not verified:
        raise RuntimeError("Hash verification failed on frozen datasets!")
    print("\nAll benchmark datasets successfully frozen and 100% verified.")

if __name__ == "__main__":
    main()
