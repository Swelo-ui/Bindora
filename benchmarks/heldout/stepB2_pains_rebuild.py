import urllib.request
import urllib.parse
import json
import csv
import os
import time
from rdkit import Chem
from rdkit.Chem import FilterCatalog

# Curated, literature-anchored PAINS positives strictly mapped to Baell & Holloway (2010) SI families
PAINS_POSITIVES = [
    # Rhodanines / Thiazolidinediones (family: ene_rhod_A / rhodanine)
    {"query": "5-benzylidenerhodanine", "family": "ene_rhod_A", "doi": "10.1021/jm901137j"},
    {"query": "4-dimethylaminobenzylidenerhodanine", "family": "ene_rhod_A", "doi": "10.1021/jm901137j"},
    {"query": "5-(4-chlorobenzylidene)rhodanine", "family": "ene_rhod_A", "doi": "10.1021/jm901137j"},
    {"query": "epalrestat", "family": "ene_rhod_A", "doi": "10.1021/jm901137j"},
    {"query": "rhodanine", "family": "rhodanine", "doi": "10.1021/jm901137j"},

    # Catechols (family: catechol_A)
    {"query": "pyrogallol", "family": "catechol_A", "doi": "10.1021/jm901137j"},
    {"query": "catechol", "family": "catechol_A", "doi": "10.1021/jm901137j"},
    {"query": "gallic acid", "family": "catechol_A", "doi": "10.1021/jm901137j"},
    {"query": "dopamine", "family": "catechol_A", "doi": "10.1021/jm901137j"},
    {"query": "epinephrine", "family": "catechol_A", "doi": "10.1021/jm901137j"},

    # Quinones (family: quinone_A)
    {"query": "1,4-naphthoquinone", "family": "quinone_A", "doi": "10.1021/jm901137j"},
    {"query": "menadione", "family": "quinone_A", "doi": "10.1021/jm901137j"},
    {"query": "juglone", "family": "quinone_A", "doi": "10.1021/jm901137j"},
    {"query": "plumbagin", "family": "quinone_A", "doi": "10.1021/jm901137j"},
    {"query": "1,4-benzoquinone", "family": "quinone_A", "doi": "10.1021/jm901137j"},

    # Hydroxyphenyl Hydrazones (family: hzone_phenol_A)
    {"query": "salicylaldehyde phenylhydrazone", "family": "hzone_phenol_A", "doi": "10.1021/jm901137j"},
    {"query": "salicylaldehyde thiosemicarbazone", "family": "hzone_phenol_A", "doi": "10.1021/jm901137j"},

    # Alkylidene Barbiturates (family: keto_barbiturate_A)
    {"query": "5-benzylidenebarbituric acid", "family": "keto_barbiturate_A", "doi": "10.1021/jm901137j"},
    {"query": "5-cinnamylidenebarbituric acid", "family": "keto_barbiturate_A", "doi": "10.1021/jm901137j"},

    # Cyanoenones / Cyanovinyls (family: ene_cyano_A)
    {"query": "alpha-cyano-4-hydroxycinnamic acid", "family": "ene_cyano_A", "doi": "10.1021/jm901137j"},
    {"query": "tyrphostin A23", "family": "ene_cyano_A", "doi": "10.1021/jm901137j"},
    {"query": "tyrphostin AG 490", "family": "ene_cyano_A", "doi": "10.1021/jm901137j"},

    # Anilines / Alkylamino-aromatics (family: anil_di_alk_B / anil_no_alk)
    {"query": "4-(dimethylamino)benzaldehyde", "family": "anil_di_alk_B", "doi": "10.1021/jm901137j"},
    {"query": "crystal violet", "family": "anil_di_alk_A", "doi": "10.1021/jm901137j"},

    # Azo compounds (family: azo_A)
    {"query": "methyl orange", "family": "azo_A", "doi": "10.1021/jm901137j"},
    {"query": "sudan I", "family": "azo_A", "doi": "10.1021/jm901137j"},
    {"query": "phenazopyridine", "family": "azo_A", "doi": "10.1021/jm901137j"},

    # Curcuminoids (canonical PAINS / chemical aggregator; Baell & Walters 2014)
    {"query": "curcumin", "family": "curcuminoid", "doi": "10.1038/513481a"},
    {"query": "demethoxycurcumin", "family": "curcuminoid", "doi": "10.1038/513481a"},
    {"query": "bisdemethoxycurcumin", "family": "curcuminoid", "doi": "10.1038/513481a"}
]

# Hard Negatives: FDA-approved drugs or known clean non-PAINS probes with enones, phenols, diketones, etc.
PAINS_NEGATIVES = [
    {"query": "chalcone", "family": "mono_enone", "doi": "10.1021/jm901137j"},
    {"query": "cinnamamide", "family": "cinnamic_amide", "doi": "10.1021/jm901137j"},
    {"query": "dibenzoylmethane", "family": "1_3_diketone", "doi": "10.1021/jm901137j"},
    {"query": "ferulic acid", "family": "cinnamic_acid", "doi": "10.1021/jm901137j"},
    {"query": "capsaicin", "family": "vanilloid", "doi": "DrugBank:DB06774"},
    {"query": "paracetamol", "family": "phenol_approved_drug", "doi": "DrugBank:DB00316"},
    {"query": "salicylic acid", "family": "phenol_approved_drug", "doi": "DrugBank:DB00936"},
    {"query": "aspirin", "family": "ester_nsaid", "doi": "DrugBank:DB00945"},
    {"query": "ibuprofen", "family": "propionic_nsaid", "doi": "DrugBank:DB01050"},
    {"query": "naproxen", "family": "naphthalene_nsaid", "doi": "DrugBank:DB00788"},
    {"query": "amoxicillin", "family": "phenolic_beta_lactam", "doi": "DrugBank:DB01060"},
    {"query": "ciprofloxacin", "family": "fluoroquinolone", "doi": "DrugBank:DB00537"},
    {"query": "metformin", "family": "biguanide", "doi": "DrugBank:DB00331"},
    {"query": "nicotine", "family": "alkaloid", "doi": "DrugBank:DB00184"},
    {"query": "caffeine", "family": "xanthine", "doi": "DrugBank:DB00201"},
    {"query": "propranolol", "family": "beta_blocker", "doi": "DrugBank:DB00571"},
    {"query": "atenolol", "family": "beta_blocker", "doi": "DrugBank:DB00335"},
    {"query": "omeprazole", "family": "sulfinyl_benzimidazole", "doi": "DrugBank:DB00338"},
    {"query": "atorvastatin", "family": "pyrrole_statin", "doi": "DrugBank:DB01076"},
    {"query": "losartan", "family": "biphenyl_tetrazole", "doi": "DrugBank:DB00678"},
    {"query": "warfarin", "family": "coumarin_anticoagulant", "doi": "DrugBank:DB00682"},
    {"query": "diazepam", "family": "benzodiazepine", "doi": "DrugBank:DB00829"},
    {"query": "morphine", "family": "opioid_phenol", "doi": "DrugBank:DB00295"},
    {"query": "penicillin G", "family": "beta_lactam", "doi": "DrugBank:DB01053"},
    {"query": "imatinib", "family": "kinase_inhibitor", "doi": "DrugBank:DB00619"},
    {"query": "gefitinib", "family": "kinase_inhibitor", "doi": "DrugBank:DB00317"},
    {"query": "erlotinib", "family": "kinase_inhibitor", "doi": "DrugBank:DB00530"},
    {"query": "dasatinib", "family": "kinase_inhibitor", "doi": "DrugBank:DB01254"},
    {"query": "sildenafil", "family": "pde5_inhibitor", "doi": "DrugBank:DB00203"},
    {"query": "sulfamethoxazole", "family": "sulfonamide", "doi": "DrugBank:DB01015"},
    {"query": "trimethoprim", "family": "diaminopyrimidine", "doi": "DrugBank:DB00440"},
    {"query": "albuterol", "family": "phenolic_agonist", "doi": "DrugBank:DB01001"},
    {"query": "loperamide", "family": "diphenylmethyl_piperidine", "doi": "DrugBank:DB00836"},
    {"query": "haloperidol", "family": "butyrophenone", "doi": "DrugBank:DB00502"},
    {"query": "maleimide", "family": "covalent_electrophile", "doi": "Thiol-reactive probe"}
]

def fetch_pubchem_verified(query_name):
    q = urllib.parse.quote(query_name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/Title,InChIKey,ConnectivitySMILES/JSON"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        d = json.loads(resp.read().decode("utf-8"))
        props = d["PropertyTable"]["Properties"][0]
        return {
            "query": query_name,
            "cid": props["CID"],
            "title": props.get("Title", query_name),
            "smiles": props["ConnectivitySMILES"],
            "inchikey": props["InChIKey"]
        }

def main():
    os.makedirs("benchmarks/heldout/results", exist_ok=True)
    
    # Initialize RDKit PAINS_A/B/C catalog
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    pains_catalog = FilterCatalog.FilterCatalog(params)
    
    verified_records = []
    
    print("--- 1. Fetching & Verifying PAINS Positives from PubChem REST ---")
    for item in PAINS_POSITIVES:
        rec = fetch_pubchem_verified(item["query"])
        rec["ground_truth"] = 1
        rec["pains_family"] = item["family"]
        rec["doi_citation"] = item["doi"]
        
        # Test RDKit matches
        mol = Chem.MolFromSmiles(rec["smiles"])
        matches = [m.GetDescription() for m in pains_catalog.GetMatches(mol)] if mol else []
        rec["rdkit_matches"] = matches
        rec["rdkit_hit"] = 1 if len(matches) > 0 else 0
        
        verified_records.append(rec)
        print(f"  [POS {len(verified_records):02d}] {rec['title']} (CID {rec['cid']}) -> InChIKey: {rec['inchikey']} | RDKit: {matches}")
        time.sleep(0.15)
        
    print("\n--- 2. Fetching & Verifying Hard Negatives from PubChem REST ---")
    for item in PAINS_NEGATIVES:
        rec = fetch_pubchem_verified(item["query"])
        rec["ground_truth"] = 0
        rec["pains_family"] = item["family"]
        rec["doi_citation"] = item["doi"]
        
        mol = Chem.MolFromSmiles(rec["smiles"])
        matches = [m.GetDescription() for m in pains_catalog.GetMatches(mol)] if mol else []
        rec["rdkit_matches"] = matches
        rec["rdkit_hit"] = 1 if len(matches) > 0 else 0
        
        verified_records.append(rec)
        print(f"  [NEG {len(verified_records):02d}] {rec['title']} (CID {rec['cid']}) -> InChIKey: {rec['inchikey']} | RDKit: {matches}")
        time.sleep(0.15)
        
    # Write verified pains_dataset.csv
    out_csv = "benchmarks/heldout/pains_dataset.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fields = ["name", "pubchem_cid", "pubchem_title", "smiles", "inchikey", "pains_ground_truth", "family", "citation"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in verified_records:
            w.writerow({
                "name": r["query"],
                "pubchem_cid": r["cid"],
                "pubchem_title": r["title"],
                "smiles": r["smiles"],
                "inchikey": r["inchikey"],
                "pains_ground_truth": r["ground_truth"],
                "family": r["pains_family"],
                "citation": r["doi_citation"]
            })
            
    # Dump raw RDKit matches
    raw_matches_path = "benchmarks/heldout/results/pains_rdkit_raw_matches.csv"
    with open(raw_matches_path, "w", newline="", encoding="utf-8") as f:
        fields = ["query_name", "pubchem_cid", "pubchem_title", "ground_truth", "family", "rdkit_hit", "rdkit_raw_matches", "smiles", "inchikey"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in verified_records:
            w.writerow({
                "query_name": r["query"],
                "pubchem_cid": r["cid"],
                "pubchem_title": r["title"],
                "ground_truth": r["ground_truth"],
                "family": r["pains_family"],
                "rdkit_hit": r["rdkit_hit"],
                "rdkit_raw_matches": "; ".join(r["rdkit_matches"]),
                "smiles": r["smiles"],
                "inchikey": r["inchikey"]
            })
            
    # Print Curcumin Keto vs Enol diagnostic
    print("\n--- 3. Special Diagnostics: Curcumin Keto vs Enol & Reactive Alerts ---")
    curc_keto_smi = "COc1cc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(OC)c2)ccc1O"
    curc_enol_smi = "COc1cc(C=CC(=O)C=C(O)C=Cc2ccc(O)c(OC)c2)ccc1O"
    m_k = Chem.MolFromSmiles(curc_keto_smi)
    m_e = Chem.MolFromSmiles(curc_enol_smi)
    print("Curcumin Keto RDKit matches:", [m.GetDescription() for m in pains_catalog.GetMatches(m_k)])
    print("Curcumin Enol RDKit matches:", [m.GetDescription() for m in pains_catalog.GetMatches(m_e)])
    
    print(f"\nSuccessfully wrote verified PAINS dataset ({len(verified_records)} compounds) to {out_csv}")
    print(f"Dumped per-compound RDKit raw matches to {raw_matches_path}")

if __name__ == "__main__":
    main()
