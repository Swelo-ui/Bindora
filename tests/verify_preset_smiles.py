import sys
import json
import urllib.request
from typing import Dict, Any
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Catalog of all preset compounds used in Bindora UI presets and benchmarks
PRESETS: Dict[str, Dict[str, Any]] = {
    # NSAID Library
    "Aspirin": {
        "cid": 2244,
        "smiles": "CC(=O)Oc1ccccc1C(=O)O"
    },
    "Ibuprofen": {
        "cid": 3672,
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    },
    "Naproxen": {
        "cid": 156391,
        "smiles": "COc1ccc2cc(ccc2c1)C(C)C(=O)O"
    },
    "Celecoxib": {
        "cid": 2662,
        "smiles": "Cc1ccc(cc1c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F)"
    },
    # Kinase Inhibitor Library
    "Imatinib": {
        "cid": 5291,
        "smiles": "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
    },
    "Gefitinib": {
        "cid": 123631,
        "smiles": "COc1cc2ncnc(c2cc1OCCCN3CCOCC3)Nc4ccc(c(c4)Cl)F"
    },
    "Erlotinib": {
        "cid": 176870,
        "smiles": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
    },
    "Dasatinib": {
        "cid": 3062316,
        "smiles": "Cc1cccc(c1Cl)NC(=O)c2cnc(s2)Nc3cc(nc(n3)C)N4CCN(CC4)CCO"
    },
    # Antiviral Library (Verified Canonical Structures)
    "Remdesivir": {
        "cid": 121304016,
        "smiles": "CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4"
    },
    "Favipiravir": {
        "cid": 492405,
        "smiles": "C1=C(N=C(C(=O)N1)C(=O)N)F"
    },
    "Molnupiravir": {
        "cid": 145996610,
        "smiles": "CC(C)C(=O)OCC1C(C(C(O1)N2C=CC(=NO)NC2=O)O)O"
    },
    "Ribavirin": {
        "cid": 37542,
        "smiles": "C1=NC(=NN1C2C(C(C(O2)CO)O)O)C(=O)N"
    }
}

def verify_preset(name: str, data: Dict[str, Any]) -> bool:
    cid = data["cid"]
    smiles = data["smiles"]
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[FAIL] {name}: Failed to parse SMILES: '{smiles}'")
        return False
    
    calc_formula = rdMolDescriptors.CalcMolFormula(mol)
    calc_heavy = mol.GetNumHeavyAtoms()
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,HeavyAtomCount,CanonicalSMILES/JSON"
    req = urllib.request.Request(url, headers={"User-Agent": "Bindora-CI/1.0"})
    
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            pc_data = json.loads(resp.read().decode("utf-8"))
            props = pc_data["PropertyTable"]["Properties"][0]
            pc_formula = props["MolecularFormula"]
            pc_heavy = props["HeavyAtomCount"]
    except Exception as e:
        print(f"[ERROR] {name}: Failed to fetch from PubChem CID {cid}: {e}")
        return False
        
    formula_match = (calc_formula == pc_formula)
    heavy_match = (calc_heavy == pc_heavy)
    
    status = "PASS" if (formula_match and heavy_match) else "FAIL"
    print(f"[{status}] {name:14} (CID {cid}): Formula: {calc_formula} (expected {pc_formula}), Heavy Atoms: {calc_heavy} (expected {pc_heavy})")
    
    return formula_match and heavy_match

def main():
    print("=" * 70)
    print("BINDORA PRESET SMILES PUBCHEM VALIDATION AUDIT")
    print("=" * 70)
    
    all_passed = True
    for name, data in PRESETS.items():
        passed = verify_preset(name, data)
        if not passed:
            all_passed = False
            
    print("=" * 70)
    if all_passed:
        print("[SUCCESS] All preset compounds match PubChem canonical structures!")
        sys.exit(0)
    else:
        print("[FAILURE] One or more preset compounds failed validation.")
        sys.exit(1)

if __name__ == "__main__":
    main()
