import urllib.request
import csv
import io
import os
import hashlib
from collections import defaultdict
from random import Random
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

DATASET_METADATA = {
    "name": "Wang_2011_Pgp_Substrate",
    "description": "P-glycoprotein (ABCB1 / MDR1) active efflux substrate classification dataset",
    "citation": "Wang, Z. et al. (2011), 'P-glycoprotein substrate models using support vector machines based on a comprehensive data set', J. Chem. Inf. Model. 51(6):1447-1456",
    "doi": "10.1021/ci200057d",
    "license": "Apache License 2.0 (datagrok-ai/admetica repository)",
    "source_url": "https://raw.githubusercontent.com/datagrok-ai/admetica/main/ADMET/absorption/pgp-substrate/pgp-substrate.csv",
    "label_definition": "Binary active efflux substrate status: 1 = Substrate (actively transported by ABCB1), 0 = Non-substrate"
}

def get_sha256(data_bytes):
    return hashlib.sha256(data_bytes).hexdigest()

def scaffold_split(data_list, seed=42, frac=(0.8, 0.2)):
    scaffolds = defaultdict(set)
    for i, row in enumerate(data_list):
        smi = row["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol:
            try:
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                scaffolds[scaffold].add(i)
            except Exception:
                scaffolds[""].add(i)
        else:
            scaffolds[""].add(i)

    total_valid = len(data_list)
    train_size = int(total_valid * frac[0])
    
    index_sets = list(scaffolds.values())
    big_index_sets = [s for s in index_sets if len(s) > (total_valid - train_size) / 2]
    small_index_sets = [s for s in index_sets if len(s) <= (total_valid - train_size) / 2]
    
    rng = Random(seed)
    rng.shuffle(big_index_sets)
    rng.shuffle(small_index_sets)
    index_sets = big_index_sets + small_index_sets

    train_indices, test_indices = [], []
    for index_set in index_sets:
        if len(train_indices) + len(index_set) <= train_size:
            train_indices.extend(index_set)
        else:
            test_indices.extend(index_set)
            
    return train_indices, test_indices

def main():
    os.makedirs("benchmarks/heldout", exist_ok=True)
    
    print("Downloading Wang et al. 2011 P-gp substrate dataset from GitHub...")
    req = urllib.request.Request(DATASET_METADATA["source_url"], headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        raw_bytes = resp.read()
        
    raw_sha = get_sha256(raw_bytes)
    DATASET_METADATA["raw_file_sha256"] = raw_sha
    print(f"Downloaded {len(raw_bytes)} bytes. SHA-256: {raw_sha}")
    
    content = raw_bytes.decode("utf-8")
    reader = csv.DictReader(io.StringIO(content))
    
    # Load dev set InChIKeys
    dev_inchikeys = set()
    dev_file = "benchmarks/heldout/dev_set.csv"
    if os.path.exists(dev_file):
        with open(dev_file, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                dev_inchikeys.add(row["inchikey"].strip())
                
    valid_rows = []
    excluded_dev = []
    for idx, r in enumerate(reader, 1):
        smi = r["Smiles"].strip()
        y_val = int(r["Substrate"].strip())
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        can_smi = Chem.MolToSmiles(mol)
        inchikey = Chem.MolToInchiKey(mol)
        
        row_obj = {
            "compound_id": f"WANG2011_{idx:03d}",
            "smiles": can_smi,
            "y_substrate": y_val,
            "inchikey": inchikey
        }
        
        if inchikey in dev_inchikeys:
            excluded_dev.append(row_obj)
        else:
            valid_rows.append(row_obj)
            
    print(f"Processed molecules: Total={len(valid_rows) + len(excluded_dev)}, Valid non-dev={len(valid_rows)}, Excluded dev matches={len(excluded_dev)}")
    
    # Bemis-Murcko scaffold split 80/20
    train_idx, test_idx = scaffold_split(valid_rows, seed=42, frac=(0.8, 0.2))
    train_data = [valid_rows[i] for i in train_idx]
    test_data = [valid_rows[i] for i in test_idx]
    
    print(f"Bemis-Murcko Split: Train={len(train_data)} (pos={sum(r['y_substrate'] for r in train_data)}), Test={len(test_data)} (pos={sum(r['y_substrate'] for r in test_data)})")
    
    # Save train & test CSVs
    train_path = "benchmarks/heldout/pgp_substrate_train.csv"
    with open(train_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "smiles", "y_substrate", "inchikey"])
        w.writeheader()
        w.writerows(train_data)
        
    test_path = "benchmarks/heldout/pgp_substrate_test.csv"
    with open(test_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "smiles", "y_substrate", "inchikey"])
        w.writeheader()
        w.writerows(test_data)
        
    # Save full verified dataset
    full_path = "benchmarks/heldout/pgp_substrate_full.csv"
    with open(full_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "smiles", "y_substrate", "inchikey"])
        w.writeheader()
        w.writerows(valid_rows)
        
    # Save metadata
    import json
    with open("benchmarks/heldout/pgp_substrate_metadata.json", "w", encoding="utf-8") as f:
        json.dump(DATASET_METADATA, f, indent=2)
        
    print("Saved pgp_substrate_train.csv, pgp_substrate_test.csv, pgp_substrate_full.csv, and metadata.")

if __name__ == "__main__":
    main()
