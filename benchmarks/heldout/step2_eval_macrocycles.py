import os
import sys
sys.path.insert(0, os.path.abspath("."))
import csv
import json
import urllib.request
from rdkit import Chem
from backend.services.macrocycle import MacrocycleConformerEngine

MACROCYCLES = [
    {
        "name": "Cyclosporine A",
        "cid": 5284373,
        "lit_perimeter": 33,
        "lit_description": "33-membered monocyclic undecapeptide",
        "citation": "PMC3682974; DrugBank DB00091"
    },
    {
        "name": "Tacrolimus (FK506)",
        "cid": 445643,
        "lit_perimeter": 23,
        "lit_description": "23-membered macrolide lactone (IUPAC: 1,14-oxa-azabicyclo[19.3.1]pentacosane core)",
        "citation": "IUPAC / PMC3658888; DrugBank DB00864"
    },
    {
        "name": "Rapamycin (Sirolimus)",
        "cid": 5284616,
        "lit_perimeter": 31,
        "lit_description": "31-membered macrolide lactone (pipecolate-containing macrolide)",
        "citation": "PMC1805569; DrugBank DB00877"
    },
    {
        "name": "Lorlatinib",
        "cid": 71731823,
        "lit_perimeter": 12,
        "lit_description": "12-membered bridged kinase macrocycle (excluding fused pyrazole atoms) or 15-membered envelope",
        "citation": "PDB: 4CLI; J. Med. Chem. 2014, 57(11), 4720; DrugBank DB12130"
    },
    {
        "name": "Vancomycin",
        "cid": 14969,
        "lit_perimeter": 16,
        "lit_description": "Tricyclic glycopeptide core with 16-membered, 12-membered, and 12-membered crosslinked peptide/ether rings",
        "citation": "PMC6882670; DrugBank DB00512"
    }
]

def main():
    os.makedirs("benchmarks/heldout/results", exist_ok=True)
    results = []
    
    for item in MACROCYCLES:
        cid = item["cid"]
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/ConnectivitySMILES,InChIKey/JSON"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            props = data["PropertyTable"]["Properties"][0]
            smi = props["ConnectivitySMILES"]
            inchikey = props["InChIKey"]
            
        mol = Chem.MolFromSmiles(smi)
        Chem.FastFindRings(mol)
        ring_info = mol.GetRingInfo()
        raw_sssr = sorted([len(r) for r in ring_info.AtomRings()], reverse=True)
        raw_bond_rings = sorted([len(r) for r in ring_info.BondRings()], reverse=True)
        
        is_macro, bindora_sizes, bindora_rings = MacrocycleConformerEngine.is_macrocycle(mol)
        max_size = bindora_sizes[0] if bindora_sizes else 0
        lit_val = item["lit_perimeter"]
        lit_detected = lit_val in bindora_sizes
        
        results.append({
            "name": item["name"],
            "cid": cid,
            "inchikey": inchikey,
            "smiles": smi,
            "lit_perimeter": lit_val,
            "lit_description": item["lit_description"],
            "citation": item["citation"],
            "raw_sssr_sizes": str(raw_sssr),
            "bindora_is_macro": is_macro,
            "bindora_all_macro_sizes": str(bindora_sizes),
            "bindora_max_size": max_size,
            "lit_perimeter_in_bindora_sizes": lit_detected,
            "contracts_ester_amide_bridges": False,
            "cycle_detection_method": "SSSR + Pairwise Bond XOR Simple Cycles"
        })
        
    out_csv = "benchmarks/heldout/results/macrocycle_eval.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "name", "cid", "inchikey", "lit_perimeter", "lit_description", "citation",
            "raw_sssr_sizes", "bindora_is_macro", "bindora_max_size", "bindora_all_macro_sizes",
            "lit_perimeter_in_bindora_sizes", "contracts_ester_amide_bridges", "cycle_detection_method"
        ]
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
        
    print(f"Macrocycle evaluation complete. Saved to {out_csv}")
    print("\nSummary:")
    for r in results:
        print(f"  {r['name']:<25} | Lit: {r['lit_perimeter']:<2} | SSSR: {r['raw_sssr_sizes']:<20} | Bindora Sizes: {r['bindora_all_macro_sizes']:<35} | Lit In Sizes: {r['lit_perimeter_in_bindora_sizes']}")

if __name__ == "__main__":
    main()
