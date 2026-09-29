import os
import sys
import subprocess
import time

def main():
    print("=" * 80)
    print("  BINODRA DOCK: INDEPENDENT HELD-OUT VALIDATION AUDIT SUITE")
    print("=" * 80)
    
    scripts = [
        ("Step 2: Macrocycle Perception", "benchmarks/heldout/step2_eval_macrocycles.py"),
        ("Step 3: PAINS Benchmark", "benchmarks/heldout/step3_pains_benchmark.py"),
        ("Step 4: P-gp Substrate Benchmark", "benchmarks/heldout/step4_pgp_benchmark.py"),
        ("Step 5: BBB Decoupling Benchmark", "benchmarks/heldout/step5_bbb_benchmark.py")
    ]
    
    t0 = time.time()
    for name, script_path in scripts:
        print(f"\n>>> Running {name} ({script_path})...")
        res = subprocess.run([sys.executable, script_path], capture_output=True, text=True, cwd=".")
        if res.returncode != 0:
            print(f"FAILED: {name}")
            print(res.stderr)
            sys.exit(res.returncode)
        else:
            print(res.stdout)
            
    total_time = time.time() - t0
    print("=" * 80)
    print(f"  ALL BENCHMARKS EXECUTED SUCCESSFULLY IN {total_time:.2f}s")
    print("=" * 80)

if __name__ == "__main__":
    main()
