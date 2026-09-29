import subprocess
import hashlib
import os
import re

def get_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def run_cmd(cmd):
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=".")
    return res.stdout.strip()

def main():
    os.makedirs("benchmarks/heldout", exist_ok=True)
    
    rev = run_cmd("git rev-parse HEAD")
    status = run_cmd("git status --short")
    diff_stat = run_cmd("git diff HEAD --stat")
    diff_full = subprocess.run("git diff HEAD", shell=True, capture_output=True, text=True, cwd=".").stdout

    # Write 00_working_tree.patch
    with open("benchmarks/heldout/00_working_tree.patch", "w", encoding="utf-8") as f:
        f.write(diff_full)

    # Get untracked files
    untracked_out = run_cmd("git ls-files --others --exclude-standard")
    untracked_files = [line.strip() for line in untracked_out.splitlines() if line.strip()]
    
    untracked_hashes = []
    for fpath in untracked_files:
        if os.path.isfile(fpath):
            try:
                untracked_hashes.append(f"{get_sha256(fpath)}  {fpath}")
            except Exception as e:
                untracked_hashes.append(f"ERROR({e})  {fpath}")
        elif os.path.isdir(fpath):
            for root, dirs, files in os.walk(fpath):
                for fname in files:
                    full = os.path.join(root, fname).replace("\\", "/")
                    try:
                        untracked_hashes.append(f"{get_sha256(full)}  {full}")
                    except Exception as e:
                        untracked_hashes.append(f"ERROR({e})  {full}")

    # Write 00_repo_state.txt
    with open("benchmarks/heldout/00_repo_state.txt", "w", encoding="utf-8") as f:
        f.write(f"=== GIT REV-PARSE HEAD ===\n{rev}\n\n")
        f.write(f"=== GIT STATUS --SHORT ===\n{status}\n\n")
        f.write(f"=== GIT DIFF HEAD --STAT ===\n{diff_stat}\n\n")
        f.write(f"=== UNTRACKED FILES SHA-256 ===\n" + "\n".join(untracked_hashes) + "\n")

    print("Step 0 repo state written successfully.")

if __name__ == "__main__":
    main()
