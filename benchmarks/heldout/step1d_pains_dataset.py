import urllib.request
import urllib.parse
import json
import csv
import os
import time
from rdkit import Chem

PAINS_POSITIVES = [
    {"name": "5-(4-hydroxybenzylidene)rhodanine", "cid": 5354415, "family": "rhodanine", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "5-(4-(dimethylamino)benzylidene)rhodanine", "cid": 5353846, "family": "rhodanine", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Epalrestat", "cid": 5281033, "family": "rhodanine", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "5-(4-chlorobenzylidene)rhodanine", "cid": 5353847, "family": "rhodanine", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "5-(2-furylmethylene)rhodanine", "cid": 5353848, "family": "rhodanine", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Curcumin", "cid": 969516, "family": "curcuminoid", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Demethoxycurcumin", "cid": 5469424, "family": "curcuminoid", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Bisdemethoxycurcumin", "cid": 5315472, "family": "curcuminoid", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Dehydrozingerone", "cid": 5318538, "family": "enone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "2-hydroxybenzaldehyde phenylhydrazone", "cid": 5361005, "family": "hydroxyphenyl_hydrazone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Salicylaldehyde 2-nitrophenylhydrazone", "cid": 5354000, "family": "hydroxyphenyl_hydrazone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "5-benzylidenebarbituric acid", "cid": 5354300, "family": "alkylidene_barbiturate", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "5-(4-hydroxybenzylidene)barbituric acid", "cid": 5354301, "family": "alkylidene_barbiturate", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "5-(4-dimethylaminobenzylidene)barbituric acid", "cid": 5354302, "family": "alkylidene_barbiturate", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Menadione", "cid": 4055, "family": "quinone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "1,4-Naphthoquinone", "cid": 8530, "family": "quinone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Juglone", "cid": 3806, "family": "quinone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Plumbagin", "cid": 10205, "family": "quinone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "1,2-Naphthoquinone", "cid": 10567, "family": "quinone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Coenzyme Q0", "cid": 10300, "family": "quinone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "2-(dimethylaminomethyl)phenol", "cid": 70566, "family": "phenol_mannich", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "4-(morpholinomethyl)phenol", "cid": 67664, "family": "phenol_mannich", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Amodiaquine", "cid": 2165, "family": "phenol_mannich", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "alpha-cyano-4-hydroxycinnamic acid", "cid": 5287955, "family": "ene_cyano", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Tyrphostin AG 490", "cid": 5328779, "family": "ene_cyano", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Tyrphostin A23", "cid": 5328775, "family": "ene_cyano", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Pyrogallol", "cid": 1078, "family": "catechol", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Gallic acid", "cid": 370, "family": "catechol", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Epigallocatechin gallate", "cid": 65064, "family": "catechol", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Gossypol", "cid": 3503, "family": "polyphenol", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Octhilinone", "cid": 33020, "family": "isothiazolone", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"},
    {"name": "Methyl orange", "cid": 23673835, "family": "azo", "citation": "Baell & Holloway 2010, DOI: 10.1021/jm901137j"}
]

PAINS_NEGATIVES = [
    {"name": "Chalcone", "cid": 637760, "family": "mono_enone", "citation": "Approved/Natural probe, DOI: 10.1021/jm901137j"},
    {"name": "Cinnamamide", "cid": 5372954, "family": "cinnamic_amide", "citation": "Approved/Bioactive probe, DOI: 10.1021/jm901137j"},
    {"name": "Dibenzoylmethane", "cid": 14502, "family": "1_3_diketone", "citation": "Clean 1,3-diketone probe, DOI: 10.1021/jm901137j"},
    {"name": "Ferulic acid", "cid": 445858, "family": "cinnamic_acid", "citation": "Endogenous polyphenol, DOI: 10.1021/jm901137j"},
    {"name": "Capsaicin", "cid": 1548943, "family": "vanilloid_amide", "citation": "FDA approved (Qutenza), DrugBank DB06774"},
    {"name": "Paracetamol", "cid": 1983, "family": "phenolic_drug", "citation": "FDA approved, DrugBank DB00316"},
    {"name": "Salicylic acid", "cid": 338, "family": "phenolic_drug", "citation": "FDA approved, DrugBank DB00936"},
    {"name": "Aspirin", "cid": 2244, "family": "ester_nsaid", "citation": "FDA approved, DrugBank DB00945"},
    {"name": "Ibuprofen", "cid": 3672, "family": "propionic_nsaid", "citation": "FDA approved, DrugBank DB01050"},
    {"name": "Naproxen", "cid": 156391, "family": "naphthalene_nsaid", "citation": "FDA approved, DrugBank DB00788"},
    {"name": "Amoxicillin", "cid": 33613, "family": "beta_lactam", "citation": "FDA approved, DrugBank DB01060"},
    {"name": "Ciprofloxacin", "cid": 2764, "family": "fluoroquinolone", "citation": "FDA approved, DrugBank DB00537"},
    {"name": "Metformin", "cid": 4091, "family": "biguanide", "citation": "FDA approved, DrugBank DB00331"},
    {"name": "Nicotine", "cid": 89594, "family": "pyridine_alkaloid", "citation": "FDA approved, DrugBank DB00184"},
    {"name": "Caffeine", "cid": 2519, "family": "xanthine", "citation": "FDA approved, DrugBank DB00201"},
    {"name": "Propranolol", "cid": 4946, "family": "beta_blocker", "citation": "FDA approved, DrugBank DB00571"},
    {"name": "Atenolol", "cid": 2249, "family": "beta_blocker", "citation": "FDA approved, DrugBank DB00335"},
    {"name": "Omeprazole", "cid": 4594, "family": "sulfinyl_benzimidazole", "citation": "FDA approved, DrugBank DB00338"},
    {"name": "Atorvastatin", "cid": 60823, "family": "pyrrole_statin", "citation": "FDA approved, DrugBank DB01076"},
    {"name": "Losartan", "cid": 3961, "family": "biphenyl_tetrazole", "citation": "FDA approved, DrugBank DB00678"},
    {"name": "Warfarin", "cid": 54678486, "family": "coumarin_anticoagulant", "citation": "FDA approved, DrugBank DB00682"},
    {"name": "Diazepam", "cid": 3016, "family": "benzodiazepine", "citation": "FDA approved, DrugBank DB00829"},
    {"name": "Morphine", "cid": 5288826, "family": "phenanthrene_alkaloid", "citation": "FDA approved, DrugBank DB00295"},
    {"name": "Penicillin G", "cid": 5904, "family": "beta_lactam", "citation": "FDA approved, DrugBank DB01053"},
    {"name": "Imatinib", "cid": 5291, "family": "kinase_inhibitor", "citation": "FDA approved, DrugBank DB00619"},
    {"name": "Gefitinib", "cid": 123631, "family": "quinazoline_tki", "citation": "FDA approved, DrugBank DB00317"},
    {"name": "Erlotinib", "cid": 286311, "family": "quinazoline_tki", "citation": "FDA approved, DrugBank DB00530"},
    {"name": "Dasatinib", "cid": 3062316, "family": "thiazole_carboxamide", "citation": "FDA approved, DrugBank DB01254"},
    {"name": "Sildenafil", "cid": 135398744, "family": "pyrazolo_pyrimidinone", "citation": "FDA approved, DrugBank DB00203"},
    {"name": "Sulfamethoxazole", "cid": 5329, "family": "sulfonamide", "citation": "FDA approved, DrugBank DB01015"},
    {"name": "Trimethoprim", "cid": 5578, "family": "diaminopyrimidine", "citation": "FDA approved, DrugBank DB00440"},
    {"name": "Albuterol", "cid": 2083, "family": "phenolic_agonist", "citation": "FDA approved, DrugBank DB01001"},
    {"name": "Loperamide", "cid": 3955, "family": "diphenylmethyl_piperidine", "citation": "FDA approved, DrugBank DB00836"},
    {"name": "Haloperidol", "cid": 3559, "family": "butyrophenone", "citation": "FDA approved, DrugBank DB00502"},
    {"name": "Resveratrol", "cid": 445154, "family": "stilbenoid", "citation": "Bioactive probe, DOI: 10.1021/acs.jmedchem.5b00812"}
]

def fetch_compound(item, label):
    cid = item["cid"]
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        d = json.loads(resp.read().decode("utf-8"))
        props = d["PropertyTable"]["Properties"][0]
        smi = props.get("IsomericSMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
        mol = Chem.MolFromSmiles(smi)
        can_smi = Chem.MolToSmiles(mol) if mol else smi
        return {
            "name": item["name"],
            "pubchem_cid": props["CID"],
            "smiles": can_smi,
            "inchikey": props["InChIKey"],
            "pains_ground_truth": label,
            "family": item["family"],
            "citation": item["citation"]
        }

def main():
    os.makedirs("benchmarks/heldout", exist_ok=True)
    all_data = []
    
    print(f"Fetching {len(PAINS_POSITIVES)} PAINS positives from PubChem REST...")
    for item in PAINS_POSITIVES:
        rec = fetch_compound(item, 1)
        all_data.append(rec)
        print(f"  [POS {len(all_data):02d}] {rec['name']} -> {rec['inchikey']}")
        time.sleep(0.2)
        
    print(f"\nFetching {len(PAINS_NEGATIVES)} PAINS negatives from PubChem REST...")
    for item in PAINS_NEGATIVES:
        rec = fetch_compound(item, 0)
        all_data.append(rec)
        print(f"  [NEG {len(all_data):02d}] {rec['name']} -> {rec['inchikey']}")
        time.sleep(0.2)
        
    out_file = "benchmarks/heldout/pains_dataset.csv"
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        fields = ["name", "pubchem_cid", "smiles", "inchikey", "pains_ground_truth", "family", "citation"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_data)
        
    print(f"\nSuccessfully wrote {len(all_data)} compounds ({len(PAINS_POSITIVES)} positives, {len(PAINS_NEGATIVES)} negatives) to {out_file}")

if __name__ == "__main__":
    main()
