import urllib.request
import urllib.parse
import json
import csv
import os
import time
from rdkit import Chem

DEV_COMPOUNDS = [
    # 50-run compounds
    {"name": "Indinavir", "cid": 5362440},
    {"name": "Biotin", "cid": 17154},
    {"name": "Benzamidine", "cid": 2332},
    {"name": "Erlotinib", "cid": 286311},
    {"name": "Imatinib", "cid": 5291},
    {"name": "Aspirin", "cid": 2244},
    {"name": "Atorvastatin", "cid": 60823},
    {"name": "Aniline", "cid": 6115},
    {"name": "Gefitinib", "cid": 123631},
    {"name": "Cyclosporine A", "cid": 5284373},
    {"name": "Pentane", "cid": 8003},
    {"name": "Ibrutinib", "cid": 24821094},
    {"name": "Nirmatrelvir", "cid": 155903259},
    {"name": "Afatinib", "cid": 24776445},
    {"name": "Dasatinib", "cid": 3062316},
    {"name": "Loperamide", "cid": 3955},
    {"name": "Curcumin", "cid": 969516},
    {"name": "5-Benzylidenerhodanine", "cid": 5354415},
    {"name": "Lorlatinib", "cid": 71731823},
    {"name": "Vancomycin", "cid": 14969},
    # test_generalization & prompt stress compounds
    {"name": "Chalcone", "cid": 637760},
    {"name": "Dibenzoylmethane", "cid": 14502},
    {"name": "Cinnamamide", "cid": 5372954},
    {"name": "Ferulic acid", "cid": 445858},
    {"name": "Capsaicin", "cid": 1548943},
    {"name": "Maleimide", "cid": 10935},
    {"name": "1,4-Benzoquinone", "cid": 4650},
    {"name": "Haloperidol", "cid": 3559},
    {"name": "Methadone", "cid": 4095},
    {"name": "Fentanyl", "cid": 3345},
    {"name": "Risperidone", "cid": 5073},
    {"name": "Verapamil", "cid": 2520},
    {"name": "Fexofenadine", "cid": 3410},
    {"name": "Tacrolimus", "cid": 445643},
    {"name": "Rapamycin", "cid": 5284616},
    {"name": "Nicotine", "cid": 89594},
    {"name": "Metformin", "cid": 4091},
    {"name": "Isobutane", "cid": 6360},
    {"name": "Hexane", "cid": 8058},
    {"name": "Caffeine", "cid": 2519},
    {"name": "Terfenadine", "cid": 5405},
    {"name": "Fendiline", "cid": 3348}
]

def fetch_pubchem(item):
    cid = item.get("cid")
    name = item["name"]
    if cid:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    else:
        q = urllib.parse.quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read().decode("utf-8"))
        props = data["PropertyTable"]["Properties"][0]
        # Prefer IsomericSMILES if present, else CanonicalSMILES / SMILES
        smi = props.get("IsomericSMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
        inchikey = props["InChIKey"]
        pub_cid = props["CID"]
        
        # Canonicalize via RDKit
        mol = Chem.MolFromSmiles(smi)
        rdkit_smi = Chem.MolToSmiles(mol) if mol else smi
        
        return {
            "name": name,
            "pubchem_cid": pub_cid,
            "smiles": rdkit_smi,
            "inchikey": inchikey
        }

def main():
    os.makedirs("benchmarks/heldout", exist_ok=True)
    out_file = "benchmarks/heldout/dev_set.csv"
    results = []
    
    print(f"Fetching {len(DEV_COMPOUNDS)} dev-set compounds from PubChem REST...")
    for item in DEV_COMPOUNDS:
        success = False
        for attempt in range(3):
            try:
                res = fetch_pubchem(item)
                results.append(res)
                print(f"  [OK] {item['name']} (CID {res['pubchem_cid']}) -> InChIKey: {res['inchikey']}")
                success = True
                break
            except Exception as e:
                print(f"  [RETRY {attempt+1}] {item['name']}: {e}")
                time.sleep(1)
        if not success:
            raise RuntimeError(f"PubChem fetch failed for {item['name']}")
        time.sleep(0.2) # polite rate limiting
        
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "pubchem_cid", "smiles", "inchikey"])
        writer.writeheader()
        writer.writerows(results)
        
    print(f"Successfully wrote {len(results)} compounds to {out_file}")

if __name__ == "__main__":
    main()
