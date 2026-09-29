import os
import sys
sys.path.insert(0, os.path.abspath("."))
import csv
import json
import urllib.request
from rdkit import Chem
from backend.services.macrocycle import MacrocycleConformerEngine
from backend.services.docking import DockingEngine

MACROCYCLES = [
    {
        "name": "Cyclosporine A",
        "cid": 5284373,
        "lit_perimeter": 33,
        "citation": "PMC3682974; DrugBank DB00091; DOI: 10.1016/j.bmcl.2013.04.015",
        "lit_note": "33-membered monocyclic undecapeptide ring (no fused or bridged rings)"
    },
    {
        "name": "Tacrolimus",
        "cid": 445643,
        "lit_perimeter": 23,
        "citation": "PMC3658888; DrugBank DB00864; DOI: 10.1039/c2np20085a",
        "lit_note": "23-membered macrolide lactone (IUPAC: 1,14-oxa-azabicyclo[19.3.1]pentacosane core with fused pipecolate ring)"
    },
    {
        "name": "Rapamycin",
        "cid": 5284616,
        "lit_perimeter": 31,
        "citation": "PMC1805569; DrugBank DB00877; DOI: 10.1021/ja00188a037",
        "lit_note": "31-membered macrolide lactone (pipecolate-containing triene macrolide)"
    },
    {
        "name": "Lorlatinib",
        "cid": 71731823,
        "lit_perimeter": 12,
        "citation": "PDB: 4CLI; DrugBank DB12130; DOI: 10.1021/jm5001712",
        "lit_note": "12-membered bridged kinase macrocycle (excluding fused pyrazole ring atoms) or 15-membered envelope"
    },
    {
        "name": "Vancomycin",
        "cid": 14969,
        "lit_perimeter": 16,
        "citation": "PMC6882670; DrugBank DB00512; DOI: 10.1038/s41467-019-13318-z",
        "lit_note": "Tricyclic heptapeptide core with crosslinked aromatic ether/biphenyl rings forming 16-membered and 12-membered rings"
    }
]

def check_simple_cycle_connectivity(mol, atom_indices):
    """Verify if the sub-graph induced by atom_indices forms a valid, connected simple cycle."""
    sub_atoms = set(atom_indices)
    deg = {a: 0 for a in sub_atoms}
    cycle_bonds = []
    
    for b in mol.GetBonds():
        u = b.GetBeginAtomIdx()
        v = b.GetEndAtomIdx()
        if u in sub_atoms and v in sub_atoms:
            deg[u] += 1
            deg[v] += 1
            cycle_bonds.append(b.GetIdx())
            
    is_simple = len(sub_atoms) >= 3 and all(d == 2 for d in deg.values())
    is_connected = len(cycle_bonds) == len(sub_atoms) if is_simple else False
    
    return {
        "is_simple_cycle": is_simple,
        "is_connected": is_connected,
        "vertex_degrees": list(deg.values()),
        "bond_count": len(cycle_bonds)
    }

def main():
    os.makedirs("benchmarks/heldout/results", exist_ok=True)
    table_records = []
    cycle_details = {}
    
    for item in MACROCYCLES:
        cid = item["cid"]
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,InChIKey,ConnectivitySMILES/JSON"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            props = data["PropertyTable"]["Properties"][0]
            
        smi = props["ConnectivitySMILES"]
        inchikey = props["InChIKey"]
        title = props["Title"]
        
        mol = Chem.MolFromSmiles(smi)
        Chem.FastFindRings(mol)
        ring_info = mol.GetRingInfo()
        
        # 1. RDKit raw SSSR
        raw_atom_rings = [list(r) for r in ring_info.AtomRings()]
        raw_sssr_sizes = sorted([len(r) for r in raw_atom_rings], reverse=True)
        
        # 2. RDKit symmetrized SSSR (Chem.GetSymmSSSR)
        symm_rings = [list(r) for r in Chem.GetSymmSSSR(mol)]
        symm_sssr_sizes = sorted([len(r) for r in symm_rings], reverse=True)
        
        # 3. Bindora backend MacrocycleConformerEngine
        is_macro, backend_sizes, backend_rings = MacrocycleConformerEngine.is_macrocycle(mol)
        
        # 4. Docking.py detection test
        mol_h = Chem.AddHs(mol)
        docking_raw_detect = any(len(r) >= 12 for r in mol_h.GetRingInfo().AtomRings()) or (
            mol_h.GetNumHeavyAtoms() > 50 and Chem.Lipinski.NumRotatableBonds(mol_h) > 15
        )
        docking_passed_is_macro = is_macro
        
        # 5. Cycle connectivity and atom index breakdown
        compound_cycles = []
        for c_idx, r_atoms in enumerate(backend_rings):
            conn = check_simple_cycle_connectivity(mol, r_atoms)
            compound_cycles.append({
                "cycle_rank": c_idx + 1,
                "size": len(r_atoms),
                "is_simple_cycle": conn["is_simple_cycle"],
                "is_connected": conn["is_connected"],
                "atom_indices": r_atoms
            })
            
        cycle_details[item["name"]] = {
            "title": title,
            "cid": cid,
            "inchikey": inchikey,
            "smiles": smi,
            "literature_perimeter": item["lit_perimeter"],
            "literature_citation": item["citation"],
            "raw_sssr_rings": raw_sssr_sizes,
            "symm_sssr_rings": symm_sssr_sizes,
            "backend_macrocycle_sizes": backend_sizes,
            "cycles_breakdown": compound_cycles
        }
        
        table_records.append({
            "name": item["name"],
            "cid": cid,
            "inchikey": inchikey,
            "literature_perimeter": item["lit_perimeter"],
            "literature_note": item["lit_note"],
            "literature_citation": item["citation"],
            "rdkit_raw_sssr_sizes": str(raw_sssr_sizes),
            "rdkit_symm_sssr_sizes": str(symm_sssr_sizes),
            "backend_macrocycle_sizes": str(backend_sizes),
            "docking_default_detect": docking_raw_detect,
            "docking_with_backend_engine": docking_passed_is_macro,
            "all_backend_cycles_are_simple_cycles": all(c["is_simple_cycle"] for c in compound_cycles)
        })
        
    # Write CSV
    out_csv = "benchmarks/heldout/results/macrocycle_discrepancy_table.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "name", "cid", "inchikey", "literature_perimeter", "literature_note", "literature_citation",
            "rdkit_raw_sssr_sizes", "rdkit_symm_sssr_sizes", "backend_macrocycle_sizes",
            "docking_default_detect", "docking_with_backend_engine", "all_backend_cycles_are_simple_cycles"
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(table_records)
        
    # Write JSON details
    out_json = "benchmarks/heldout/results/macrocycle_cycle_breakdown.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(cycle_details, f, indent=2)
        
    print(f"Macrocycle discrepancy analysis complete. Saved to {out_csv} and {out_json}")
    for r in table_records:
        print(f"{r['name']:<18} | Lit: {r['literature_perimeter']:<2} | SSSR: {r['rdkit_raw_sssr_sizes']:<22} | Symm: {r['rdkit_symm_sssr_sizes']:<22} | Backend: {r['backend_macrocycle_sizes']}")

if __name__ == "__main__":
    main()
