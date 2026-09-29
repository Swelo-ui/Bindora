import urllib.request
import csv
import io
import os
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

def download_b3db():
    url = "https://raw.githubusercontent.com/theochem/B3DB/main/B3DB/B3DB_classification.tsv"
    print("Downloading B3DB classification dataset from GitHub...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        content = resp.read().decode("utf-8")
    return content

def main():
    os.makedirs("benchmarks/heldout", exist_ok=True)
    raw_content = download_b3db()
    reader = csv.DictReader(io.StringIO(raw_content), delimiter="\t")
    
    # Load dev set InChIKeys
    dev_inchikeys = set()
    dev_file = "benchmarks/heldout/dev_set.csv"
    if os.path.exists(dev_file):
        with open(dev_file, "r", encoding="utf-8") as f:
            dev_reader = csv.DictReader(f)
            for row in dev_reader:
                dev_inchikeys.add(row["inchikey"].strip())
    print(f"Loaded {len(dev_inchikeys)} dev-set InChIKeys to exclude")
    
    clean_rows = []
    excluded_count = 0
    parse_errors = 0
    
    for row in reader:
        smi = row.get("SMILES", "").strip()
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            parse_errors += 1
            continue
            
        can_smi = Chem.MolToSmiles(mol)
        inchikey = Chem.MolToInchiKey(mol)
        
        if inchikey in dev_inchikeys:
            excluded_count += 1
            continue
            
        label_raw = row.get("BBB+/BBB-", "").strip()
        if label_raw not in ["BBB+", "BBB-"]:
            continue
            
        binary_label = 1 if label_raw == "BBB+" else 0
        
        clean_rows.append({
            "NO": row.get("NO.", ""),
            "compound_name": row.get("compound_name", ""),
            "smiles": can_smi,
            "pubchem_cid": row.get("CID", ""),
            "logBB": row.get("logBB", ""),
            "bbb_label": label_raw,
            "y_binary": binary_label,
            "active_efflux": "NA_not_recorded_in_b3db",
            "group": row.get("group", ""),
            "inchikey": inchikey
        })
        
    print(f"B3DB processing complete: Total valid={len(clean_rows)}, Excluded contaminated={excluded_count}, Parse errors={parse_errors}")
    
    out_file = "benchmarks/heldout/bbb_test.csv"
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["NO", "compound_name", "smiles", "pubchem_cid", "logBB", "bbb_label", "y_binary", "active_efflux", "group", "inchikey"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(clean_rows)
        
    print(f"Saved {len(clean_rows)} compounds to {out_file}")

if __name__ == "__main__":
    main()
