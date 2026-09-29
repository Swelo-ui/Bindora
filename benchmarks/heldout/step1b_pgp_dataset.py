import urllib.request
import csv
import io
import os
from collections import defaultdict
from random import Random
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

def download_pgp_data():
    url = "https://dataverse.harvard.edu/api/access/datafile/4259597"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    print("Downloading TDC Pgp_Broccatelli from Harvard Dataverse...")
    with urllib.request.urlopen(req, timeout=30) as resp:
        content = resp.read().decode("utf-8")
    return content

def tdc_scaffold_split(data_list, seed=42, frac=(0.8, 0.2)):
    # data_list is a list of dicts with 'Drug_ID', 'Drug' (SMILES), 'Y'
    scaffolds = defaultdict(set)
    error_smiles = 0
    
    for i, row in enumerate(data_list):
        smi = row["Drug"]
        mol = Chem.MolFromSmiles(smi)
        if mol:
            try:
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                scaffolds[scaffold].add(i)
            except Exception:
                error_smiles += 1
        else:
            error_smiles += 1

    total_valid = len(data_list) - error_smiles
    train_size = int(total_valid * frac[0])
    test_size = total_valid - train_size
    
    index_sets = list(scaffolds.values())
    big_index_sets = []
    small_index_sets = []
    
    for index_set in index_sets:
        if len(index_set) > test_size / 2:
            big_index_sets.append(index_set)
        else:
            small_index_sets.append(index_set)
            
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
    raw_content = download_pgp_data()
    
    # Parse TSV
    reader = csv.DictReader(io.StringIO(raw_content), delimiter="\t")
    raw_rows = list(reader)
    print(f"Total raw rows in Pgp_Broccatelli: {len(raw_rows)}")
    
    # Process molecules, get canonical SMILES and InChIKeys
    valid_rows = []
    for r in raw_rows:
        smi = r["Drug"].strip().strip('"')
        mol = Chem.MolFromSmiles(smi)
        if mol:
            can_smi = Chem.MolToSmiles(mol)
            inchikey = Chem.MolToInchiKey(mol)
            valid_rows.append({
                "Drug_ID": r["Drug_ID"].strip().strip('"'),
                "Drug": can_smi,
                "Y": int(float(r["Y"].strip())),
                "InChIKey": inchikey
            })
            
    print(f"Valid molecules: {len(valid_rows)}")
    
    # Run TDC Bemis-Murcko scaffold split
    train_idx, test_idx = tdc_scaffold_split(valid_rows, seed=42, frac=(0.8, 0.2))
    print(f"Initial split sizes -> Train: {len(train_idx)}, Test: {len(test_idx)}")
    
    # Load dev set InChIKeys
    dev_inchikeys = set()
    dev_file = "benchmarks/heldout/dev_set.csv"
    if os.path.exists(dev_file):
        with open(dev_file, "r", encoding="utf-8") as f:
            dev_reader = csv.DictReader(f)
            for row in dev_reader:
                dev_inchikeys.add(row["inchikey"].strip())
    print(f"Loaded {len(dev_inchikeys)} dev-set InChIKeys to exclude")
    
    # Filter test set
    test_rows = [valid_rows[i] for i in test_idx]
    train_rows = [valid_rows[i] for i in train_idx]
    
    filtered_test = []
    excluded_test = []
    for r in test_rows:
        if r["InChIKey"] in dev_inchikeys:
            excluded_test.append(r)
        else:
            filtered_test.append(r)
            
    excluded_train = [r for r in train_rows if r["InChIKey"] in dev_inchikeys]
    print(f"Excluded from test set due to dev contamination: {len(excluded_test)} ({[r['Drug_ID'] for r in excluded_test]})")
    print(f"Found in train set matching dev contamination: {len(excluded_train)} ({[r['Drug_ID'] for r in excluded_train]})")
    print(f"Final sizes -> Train: {len(train_rows)}, Test: {len(filtered_test)}")
    
    # Save train.csv and test.csv
    with open("benchmarks/heldout/train.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Drug_ID", "Drug", "Y", "InChIKey"])
        w.writeheader()
        w.writerows(train_rows)
        
    with open("benchmarks/heldout/test.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Drug_ID", "Drug", "Y", "InChIKey"])
        w.writeheader()
        w.writerows(filtered_test)
        
    print("Saved benchmarks/heldout/train.csv and benchmarks/heldout/test.csv")

if __name__ == "__main__":
    main()
