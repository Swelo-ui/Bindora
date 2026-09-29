import os
import sys
sys.path.insert(0, os.path.abspath("."))
import csv
import math
from rdkit import Chem
from rdkit.Chem import FilterCatalog
from backend.services.adme import ADMEProfiler, _EXTENDED_PAINS_PATTERNS

def wilson_ci(k, n, z=1.95996):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denominator = 1 + (z ** 2) / n
    centre_adjusted_probability = p + (z ** 2) / (2 * n)
    adjusted_sd = math.sqrt((p * (1 - p) + (z ** 2) / (4 * n)) / n)
    lower_bound = (centre_adjusted_probability - z * adjusted_sd) / denominator
    upper_bound = (centre_adjusted_probability + z * adjusted_sd) / denominator
    return (max(0.0, round(lower_bound, 4)), min(1.0, round(upper_bound, 4)))

def compute_metrics(tp, fp, tn, fn):
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    sens_ci = wilson_ci(tp, tp + fn)
    spec_ci = wilson_ci(tn, tn + fp)
    
    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc = (tp * tn - fp * fn) / math.sqrt(denom) if denom > 0 else 0.0
    bacc = (sens + spec) / 2.0
    
    return {
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "Sensitivity": round(sens, 4),
        "Sensitivity_95CI": sens_ci,
        "Specificity": round(spec, 4),
        "Specificity_95CI": spec_ci,
        "Balanced_Accuracy": round(bacc, 4),
        "MCC": round(mcc, 4)
    }

def main():
    os.makedirs("benchmarks/heldout/results", exist_ok=True)
    
    # Initialize RDKit PAINS_A/B/C catalog
    pains_params = FilterCatalog.FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    pains_catalog = FilterCatalog.FilterCatalog(pains_params)
    
    # Read pains_dataset.csv
    dataset_file = "benchmarks/heldout/pains_dataset.csv"
    with open(dataset_file, "r", encoding="utf-8") as f:
        reader = list(csv.DictReader(f))
        
    records = []
    
    # Counters
    # rdkit alone
    rd_tp = rd_fp = rd_tn = rd_fn = 0
    # bindora full
    bf_tp = bf_fp = bf_tn = bf_fn = 0
    # bindora extended only
    be_tp = be_fp = be_tn = be_fn = 0
    
    from rdkit.Chem.MolStandardize import rdMolStandardize
    taut_enumerator = rdMolStandardize.TautomerEnumerator()

    for row in reader:
        name = row["name"]
        smi = row["smiles"]
        gt = int(row["pains_ground_truth"])
        fam = row["family"]
        
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
            
        # 1. RDKit Alone
        rd_matches = [m.GetDescription() for m in pains_catalog.GetMatches(mol)]
        rd_pred = 1 if len(rd_matches) > 0 else 0
        if gt == 1:
            if rd_pred == 1: rd_tp += 1
            else: rd_fn += 1
        else:
            if rd_pred == 1: rd_fp += 1
            else: rd_tn += 1
            
        # 2. Bindora Full
        adme_res = ADMEProfiler.calculate_adme(mol)
        bf_matches = adme_res["medicinal_chemistry_safety"]["pains_alerts"]["alerts"]
        bf_count = adme_res["medicinal_chemistry_safety"]["pains_alerts"]["count"]
        bf_pred = 1 if bf_count > 0 else 0
        if gt == 1:
            if bf_pred == 1: bf_tp += 1
            else: bf_fn += 1
        else:
            if bf_pred == 1: bf_fp += 1
            else: bf_tn += 1
            
        # 3. Bindora Extended Only
        try:
            can_mol = taut_enumerator.Canonicalize(mol)
        except Exception:
            can_mol = mol
        be_matches = []
        for ep in _EXTENDED_PAINS_PATTERNS:
            patt = Chem.MolFromSmarts(ep["smarts"])
            if patt and (mol.HasSubstructMatch(patt) or can_mol.HasSubstructMatch(patt)):
                be_matches.append(ep["description"])
        be_pred = 1 if len(be_matches) > 0 else 0
        if gt == 1:
            if be_pred == 1: be_tp += 1
            else: be_fn += 1
        else:
            if be_pred == 1: be_fp += 1
            else: be_tn += 1
            
        records.append({
            "name": name,
            "smiles": smi,
            "inchikey": row["inchikey"],
            "ground_truth": gt,
            "family": fam,
            "rdkit_alone_pred": rd_pred,
            "rdkit_alone_alerts": "; ".join(rd_matches),
            "bindora_full_pred": bf_pred,
            "bindora_full_alerts": "; ".join(bf_matches),
            "bindora_extended_pred": be_pred,
            "bindora_extended_alerts": "; ".join(be_matches)
        })
        
    # Write detailed CSV
    out_csv = "benchmarks/heldout/results/pains_benchmark.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "name", "inchikey", "ground_truth", "family",
            "rdkit_alone_pred", "rdkit_alone_alerts",
            "bindora_full_pred", "bindora_full_alerts",
            "bindora_extended_pred", "bindora_extended_alerts", "smiles"
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(records)
        
    m_rd = compute_metrics(rd_tp, rd_fp, rd_tn, rd_fn)
    m_bf = compute_metrics(bf_tp, bf_fp, bf_tn, bf_fn)
    m_be = compute_metrics(be_tp, be_fp, be_tn, be_fn)
    
    # Save metrics summary
    summary_path = "benchmarks/heldout/results/pains_metrics_summary.json"
    import json
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "RDKit_PAINS_ABC_Alone": m_rd,
            "Bindora_Full": m_bf,
            "Bindora_Extended_Only": m_be
        }, f, indent=2)
        
    print(f"PAINS benchmark complete. Saved to {out_csv}")
    print("\n=== METRICS SUMMARY ===")
    print("1. RDKit Alone (PAINS A/B/C):", m_rd)
    print("2. Bindora Full:", m_bf)
    print("3. Bindora Extended Only:", m_be)
    
    print("\n=== DISCREPANCIES (False Positives / False Negatives in Bindora Full) ===")
    for r in records:
        if r["ground_truth"] == 1 and r["bindora_full_pred"] == 0:
            print(f"  [FALSE NEGATIVE] {r['name']} (family: {r['family']}) -> No alert triggered")
        elif r["ground_truth"] == 0 and r["bindora_full_pred"] == 1:
            print(f"  [FALSE POSITIVE] {r['name']} (family: {r['family']}) -> Triggered: {r['bindora_full_alerts']}")

if __name__ == "__main__":
    main()
