import os
import re
import subprocess
import hashlib
from pathlib import Path

def get_file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

DRUG_NAMES = [
    "curcumin", "loperamide", "pentane", "gefitinib", "erlotinib",
    "haloperidol", "terfenadine", "nicotine", "metformin", "tacrolimus",
    "rapamycin", "cyclosporin", "lorlatinib", "vancomycin", "chalcone",
    "maleimide", "aspirin", "caffeine", "diazepam", "morphine", "isobutane",
    "aniline", "biotin", "benzamidine", "imatinib", "dasatinib", "ibrutinib",
    "nirmatrelvir", "afatinib", "indinavir", "atorvastatin", "fendiline",
    "fexofenadine", "methadone", "fentanyl", "risperidone", "verapamil"
]

# Known citations for standard thresholds
CITATION_MAP = {
    "mw > 500": "Lipinski et al. 1997, Adv. Drug Deliv. Rev. 23:3-25",
    "logp > 5": "Lipinski et al. 1997, Adv. Drug Deliv. Rev. 23:3-25",
    "hbd > 5": "Lipinski et al. 1997, Adv. Drug Deliv. Rev. 23:3-25",
    "hba > 10": "Lipinski et al. 1997, Adv. Drug Deliv. Rev. 23:3-25",
    "rotb > 10": "Veber et al. 2002, J. Med. Chem. 45:2615-2623",
    "tpsa > 140": "Veber et al. 2002, J. Med. Chem. 45:2615-2623",
    "160 <= mw <= 480": "Ghose et al. 1999, J. Comb. Chem. 1:55-68",
    "-0.4 <= logp <= 5.6": "Ghose et al. 1999, J. Comb. Chem. 1:55-68",
    "40 <= mr <= 130": "Ghose et al. 1999, J. Comb. Chem. 1:55-68",
    "20 <= total_atoms <= 70": "Ghose et al. 1999, J. Comb. Chem. 1:55-68",
    "sa_score < 3.0": "Ertl & Schuffenhauer 2009, J. Cheminform. 1:8",
    "sa_score <= 6.0": "Ertl & Schuffenhauer 2009, J. Cheminform. 1:8",
    "min_ring_size = 12": "Macrocycle definition: Villar et al. 2014 / ChemMedChem",
    "random_seed = 42": "Bindora heuristic (uncalibrated default seed)",
    "energy_window = 15.0": "Bindora heuristic (uncalibrated conformational window)",
    "rmsd_threshold = 0.5": "Bindora heuristic (uncalibrated RMSD clustering)",
    "threshold = -6.0": "Bindora heuristic (uncalibrated flat binding affinity gate)"
}

def scan_backend():
    findings_drugs = []
    findings_smiles = []
    findings_numeric = []
    
    backend_dir = Path("backend")
    py_files = sorted(list(backend_dir.rglob("*.py")))
    
    for fpath in py_files:
        rel_path = str(fpath).replace("\\", "/")
        with open(fpath, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
            
        for line_num, line in enumerate(lines, 1):
            stripped = line.strip()
            if not stripped:
                continue
                
            lower_line = line.lower()
            is_comment = stripped.startswith(("#", "//", "/*", "*", '"""', "'''"))
            
            # (i) Drug names
            for drug in DRUG_NAMES:
                if re.search(r'\b' + re.escape(drug) + r'\b', lower_line):
                    comment_flag = is_comment or ("#" in line and line.index("#") < lower_line.index(drug))
                    cat = "comment/description" if comment_flag else "in logic"
                    findings_drugs.append({
                        "file": rel_path,
                        "line": line_num,
                        "match": drug,
                        "category": cat,
                        "content": stripped
                    })
                    
            # (ii) SMILES / SMARTS literals
            # Find quoted strings that look like chemical strings or SMARTS
            smiles_candidates = re.finditer(r'["\']([A-Za-z0-9@\+\-\#\$\:\/\\\[\]\(\)\=\%~]{5,})["\']', line)
            for sm in smiles_candidates:
                candidate = sm.group(1)
                # Ignore common python strings, import paths, mime types
                if candidate in ["utf-8", "application/json", "error", "success", "message", "status", "warning"]:
                    continue
                # Heuristic for chemical SMARTS / SMILES
                has_chem_symbols = any(c in candidate for c in ["=", "#", "[", "]", "@", "/", "\\", "~", ":", "$"]) or (
                    re.search(r'[CNOFSIPClBr]', candidate) and re.search(r'\d', candidate)
                )
                if has_chem_symbols:
                    comment_flag = is_comment or ("#" in line and line.index("#") < sm.start())
                    cat = "comment/description" if comment_flag else "in logic"
                    findings_smiles.append({
                        "file": rel_path,
                        "line": line_num,
                        "match": candidate,
                        "category": cat,
                        "content": stripped
                    })
                    
            # (iii) Numeric hardcoded defaults and thresholds
            # Look for thresholds, seeds, cutoffs, window values
            num_patterns = [
                (r'\bseed\s*=\s*(\d+)', "seed"),
                (r'\brandom_seed\s*=\s*(\d+)', "random_seed"),
                (r'\bmin_ring_size\s*=\s*(\d+)', "min_ring_size"),
                (r'\benergy_window\s*=\s*([0-9\.]+)', "energy_window"),
                (r'\brmsd_threshold\s*=\s*([0-9\.]+)', "rmsd_threshold"),
                (r'\bmw\s*(>|<|>=|<=)\s*([0-9\.]+)', "mw_threshold"),
                (r'\blogp\s*(>|<|>=|<=)\s*([0-9\.]+)', "logp_threshold"),
                (r'\btpsa\s*(>|<|>=|<=)\s*([0-9\.]+)', "tpsa_threshold"),
                (r'\brotb\s*(>|<|>=|<=)\s*([0-9\.]+)', "rotb_threshold"),
                (r'(-[456789]\.[0-9]+)\b', "affinity_cutoff"),
                (r'\bsa_score\s*(<|<=|>|>=)\s*([0-9\.]+)', "sa_score_threshold")
            ]
            for pat, p_type in num_patterns:
                m = re.search(pat, line, re.IGNORECASE)
                if m:
                    comment_flag = is_comment or ("#" in line and line.index("#") < m.start())
                    cat = "comment/description" if comment_flag else "in logic"
                    
                    # Determine citation or uncalibrated label
                    citation = "Bindora heuristic (uncalibrated)"
                    for k, cit in CITATION_MAP.items():
                        if k in line.lower() or p_type in k:
                            citation = cit
                            break
                            
                    findings_numeric.append({
                        "file": rel_path,
                        "line": line_num,
                        "type": p_type,
                        "match": m.group(0),
                        "category": cat,
                        "citation": citation,
                        "content": stripped
                    })

    return findings_drugs, findings_smiles, findings_numeric

def main():
    os.makedirs("benchmarks/heldout", exist_ok=True)
    
    # 1. Write A1_repo_state.txt
    git_rev = subprocess.run("git rev-parse HEAD", shell=True, capture_output=True, text=True).stdout.strip()
    git_stat = subprocess.run("git status --short", shell=True, capture_output=True, text=True).stdout.strip()
    git_diff_stat = subprocess.run("git diff baseline-2026-09-29^..baseline-2026-09-29 --stat", shell=True, capture_output=True, text=True).stdout.strip()
    
    with open("benchmarks/heldout/A1_repo_state.txt", "w", encoding="utf-8") as f:
        f.write(f"=== GIT COMMIT HEAD ===\n{git_rev}\n\n")
        f.write(f"=== GIT TAG ===\nbaseline-2026-09-29\n\n")
        f.write(f"=== GIT STATUS ===\n{git_stat if git_stat else 'Clean'}\n\n")
        f.write(f"=== BASELINE COMMIT STATS ===\n{git_diff_stat}\n")
        
    print("A1 repo state written to benchmarks/heldout/A1_repo_state.txt")
    
    # 2. Run scan across backend/
    drugs, smiles, numerics = scan_backend()
    
    with open("benchmarks/heldout/A2_backend_audit_scan.txt", "w", encoding="utf-8") as f:
        f.write("================================================================================\n")
        f.write("  STAGE A2: BACKEND STATIC CODE AUDIT SCAN REPORT\n")
        f.write("  Scanned directory: backend/ (*.py)\n")
        f.write("================================================================================\n\n")
        
        f.write(f"--- 1. DRUG NAME OCCURRENCES ({len(drugs)} findings) ---\n")
        for d in drugs:
            f.write(f"[{d['category'].upper()}] {d['file']}:{d['line']} -> Match: '{d['match']}'\n")
            f.write(f"  Line: {d['content']}\n\n")
            
        f.write(f"\n--- 2. SMILES / SMARTS LITERALS ({len(smiles)} findings) ---\n")
        for s in smiles:
            f.write(f"[{s['category'].upper()}] {s['file']}:{s['line']} -> Match: '{s['match']}'\n")
            f.write(f"  Line: {s['content']}\n\n")
            
        f.write(f"\n--- 3. NUMERIC HARDCODED DEFAULTS & THRESHOLDS ({len(numerics)} findings) ---\n")
        for n in numerics:
            f.write(f"[{n['category'].upper()}] {n['file']}:{n['line']} -> {n['type']} ('{n['match']}')\n")
            f.write(f"  Source/Label: {n['citation']}\n")
            f.write(f"  Line: {n['content']}\n\n")
            
    print(f"A2 backend audit scan written to benchmarks/heldout/A2_backend_audit_scan.txt")
    print(f"Total findings: Drugs={len(drugs)}, SMILES/SMARTS={len(smiles)}, Numeric defaults={len(numerics)}")

if __name__ == "__main__":
    main()
