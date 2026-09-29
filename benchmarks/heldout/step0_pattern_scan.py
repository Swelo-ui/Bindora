import os
import re

DRUG_NAMES = [
    "curcumin", "loperamide", "pentane", "gefitinib", "erlotinib",
    "haloperidol", "terfenadine", "nicotine", "metformin", "tacrolimus",
    "rapamycin", "cyclosporin", "lorlatinib", "vancomycin", "chalcone",
    "maleimide", "aspirin", "caffeine", "diazepam", "morphine", "isobutane"
]

def scan_patch(patch_path):
    findings = []
    current_file = None
    line_num = 0
    
    with open(patch_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            if raw_line.startswith("diff --git"):
                parts = raw_line.split()
                current_file = parts[2].replace("a/", "") if len(parts) > 2 else "unknown"
                line_num = 0
                continue
            elif raw_line.startswith("@@"):
                m = re.search(r'\+(\d+)', raw_line)
                if m:
                    line_num = int(m.group(1)) - 1
                continue
            
            if raw_line.startswith("+") and not raw_line.startswith("+++"):
                line_num += 1
                line = raw_line[1:] # strip leading +
                
                # Check drug names
                lower_line = line.lower()
                for drug in DRUG_NAMES:
                    if re.search(r'\b' + re.escape(drug) + r'\b', lower_line):
                        is_comment = line.strip().startswith(("#", "//", "/*", "*", '"""', "'''")) or ("#" in line and line.index("#") < lower_line.index(drug))
                        cat = "comment/description" if is_comment else "in logic"
                        findings.append({
                            "type": "drug_name",
                            "match": drug,
                            "file": current_file,
                            "line": line_num,
                            "category": cat,
                            "content": line.strip()
                        })
                
                # Check SMILES / SMARTS literals
                # Strings enclosed in quotes with chemical symbols
                smiles_matches = re.finditer(r'["\']([A-Za-z0-9@\+\-\#\$\:\/\\\[\]\(\)\=\%]{7,})["\']', line)
                for sm in smiles_matches:
                    candidate = sm.group(1)
                    # Heuristic for smiles/smarts: contains typical atoms and ring digits/bonds/brackets
                    if any(c in candidate for c in ["=", "#", "[", "]", "@", "/", "\\"]) or (re.search(r'[CNOFSIPClBr]', candidate) and re.search(r'\d', candidate)):
                        is_comment = line.strip().startswith(("#", "//", "/*", "*", '"""', "'''")) or ("#" in line and line.index("#") < sm.start())
                        cat = "comment/description" if is_comment else "in logic"
                        findings.append({
                            "type": "smiles_literal",
                            "match": candidate,
                            "file": current_file,
                            "line": line_num,
                            "category": cat,
                            "content": line.strip()
                        })
            elif not raw_line.startswith("-"):
                line_num += 1

    return findings

def scan_files(file_list):
    findings = []
    for fpath in file_list:
        if not os.path.exists(fpath):
            continue
        with open(fpath, "r", encoding="utf-8", errors="ignore") as f:
            for idx, line in enumerate(f, 1):
                lower_line = line.lower()
                for drug in DRUG_NAMES:
                    if re.search(r'\b' + re.escape(drug) + r'\b', lower_line):
                        is_comment = line.strip().startswith(("#", "//", "/*", "*", '"""', "'''")) or ("#" in line and line.index("#") < lower_line.index(drug))
                        cat = "comment/description" if is_comment else "in logic"
                        findings.append({
                            "type": "drug_name",
                            "match": drug,
                            "file": fpath,
                            "line": idx,
                            "category": cat,
                            "content": line.strip()
                        })
                smiles_matches = re.finditer(r'["\']([A-Za-z0-9@\+\-\#\$\:\/\\\[\]\(\)\=\%]{7,})["\']', line)
                for sm in smiles_matches:
                    candidate = sm.group(1)
                    if any(c in candidate for c in ["=", "#", "[", "]", "@", "/", "\\"]) or (re.search(r'[CNOFSIPClBr]', candidate) and re.search(r'\d', candidate)):
                        is_comment = line.strip().startswith(("#", "//", "/*", "*", '"""', "'''")) or ("#" in line and line.index("#") < sm.start())
                        cat = "comment/description" if is_comment else "in logic"
                        findings.append({
                            "type": "smiles_literal",
                            "match": candidate,
                            "file": fpath,
                            "line": idx,
                            "category": cat,
                            "content": line.strip()
                        })
    return findings

def main():
    patch_file = "benchmarks/heldout/00_working_tree.patch"
    findings = []
    if os.path.exists(patch_file):
        findings.extend(scan_patch(patch_file))
    
    # Check untracked files
    untracked_py = []
    for root, dirs, files in os.walk("."):
        if ".git" in root or "heldout" in root:
            continue
        for file in files:
            if file.endswith(".py"):
                rel = os.path.relpath(os.path.join(root, file), ".").replace("\\", "/")
                # check if untracked
                untracked_py.append(rel)
                
    output_path = "benchmarks/heldout/00_pattern_scan.txt"
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"=== PATTERN SCAN REPORT ===\n")
        f.write(f"Scanned files/patch for drug names: {DRUG_NAMES}\n\n")
        f.write(f"Total findings: {len(findings)}\n\n")
        for item in findings:
            f.write(f"[{item['type'].upper()}] {item['file']}:{item['line']} ({item['category']})\n")
            f.write(f"  Match: {item['match']}\n")
            f.write(f"  Line:  {item['content']}\n\n")
            
    print(f"Pattern scan written to {output_path}. Total findings: {len(findings)}")

if __name__ == "__main__":
    main()
