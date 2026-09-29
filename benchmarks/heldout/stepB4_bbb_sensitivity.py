import os
import sys
sys.path.insert(0, os.path.abspath("."))
import csv
import json
import math
from rdkit import Chem
from rdkit.Chem import Descriptors
from backend.services.adme import _point_in_polygon, _BBB_COORDS, predict_pgp_substrate

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

def calc_metrics(tp, fp, tn, fn):
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    bacc = (sens + spec) / 2.0
    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc = (tp * tn - fp * fn) / math.sqrt(denom) if denom > 0 else 0.0
    return {
        "TP": int(tp), "FP": int(fp), "TN": int(tn), "FN": int(fn),
        "Sensitivity": round(sens, 4),
        "Sensitivity_95CI": wilson_ci(tp, tp + fn),
        "Specificity": round(spec, 4),
        "Specificity_95CI": wilson_ci(tn, tn + fp),
        "Balanced_Accuracy": round(bacc, 4),
        "MCC": round(mcc, 4)
    }

def evaluate_subset(rows):
    dec_tp = dec_fp = dec_tn = dec_fn = 0
    cpl_tp = cpl_fp = cpl_tn = cpl_fn = 0
    extra_fn = 0
    
    for r in rows:
        y_true = int(r["y_binary"])
        smi = r["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
            
        tpsa = round(Descriptors.TPSA(mol, includeSandP=True), 2)
        wlogp = round(Descriptors.MolLogP(mol), 2)
        
        inside_egg = _point_in_polygon(tpsa, wlogp, _BBB_COORDS)
        pred_decoupled = 1 if inside_egg else 0
        
        pgp_res = predict_pgp_substrate(mol)
        is_pgp = 1 if pgp_res["is_substrate"] else 0
        pred_coupled = 1 if (inside_egg and not is_pgp) else 0
        
        if y_true == 1 and pred_decoupled == 1 and pred_coupled == 0:
            extra_fn += 1
            
        if y_true == 1:
            if pred_decoupled == 1: dec_tp += 1
            else: dec_fn += 1
            if pred_coupled == 1: cpl_tp += 1
            else: cpl_fn += 1
        else:
            if pred_decoupled == 1: dec_fp += 1
            else: dec_tn += 1
            if pred_coupled == 1: cpl_fp += 1
            else: cpl_tn += 1
            
    return {
        "count": len(rows),
        "decoupled": calc_metrics(dec_tp, dec_fp, dec_tn, dec_fn),
        "coupled": calc_metrics(cpl_tp, cpl_fp, cpl_tn, cpl_fn),
        "extra_false_negatives": extra_fn,
        "fp_reduction": dec_fp - cpl_fp
    }

def main():
    os.makedirs("benchmarks/heldout/results", exist_ok=True)
    in_file = "benchmarks/heldout/bbb_test.csv"
    
    with open(in_file, "r", encoding="utf-8") as f:
        all_rows = list(csv.DictReader(f))
        
    subsets = {
        "Full_B3DB": all_rows,
        "Group_A_Only_Numerical_logBB": [r for r in all_rows if r.get("group") == "A"],
        "Non_Threshold_Groups_A_and_C": [r for r in all_rows if r.get("group") in ["A", "C"]],
        "Excluding_Uncertain_Group_D": [r for r in all_rows if r.get("group") != "D"]
    }
    
    results = {}
    csv_rows = []
    
    print("Running BBB sensitivity analysis across 4 confidence subsets...")
    for sub_name, rows in subsets.items():
        res = evaluate_subset(rows)
        results[sub_name] = res
        d = res["decoupled"]
        c = res["coupled"]
        print(f"\n--- {sub_name} (N = {res['count']}) ---")
        print(f"  Decoupled: Sens={d['Sensitivity']:.4f} | Spec={d['Specificity']:.4f} | BAcc={d['Balanced_Accuracy']:.4f} | MCC={d['MCC']:.4f}")
        print(f"  Coupled:   Sens={c['Sensitivity']:.4f} | Spec={c['Specificity']:.4f} | BAcc={c['Balanced_Accuracy']:.4f} | MCC={c['MCC']:.4f}")
        print(f"  Impact: Extra False Negatives = {res['extra_false_negatives']} | FP Reduction = {res['fp_reduction']}")
        
        csv_rows.append({
            "subset": sub_name,
            "sample_size": res["count"],
            "decoupled_sens": d["Sensitivity"],
            "decoupled_spec": d["Specificity"],
            "decoupled_bacc": d["Balanced_Accuracy"],
            "decoupled_mcc": d["MCC"],
            "coupled_sens": c["Sensitivity"],
            "coupled_spec": c["Specificity"],
            "coupled_bacc": c["Balanced_Accuracy"],
            "coupled_mcc": c["MCC"],
            "extra_false_negatives_from_coupling": res["extra_false_negatives"],
            "fp_reduction_from_coupling": res["fp_reduction"]
        })
        
    out_csv = "benchmarks/heldout/results/bbb_sensitivity_analysis.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "subset", "sample_size",
            "decoupled_sens", "decoupled_spec", "decoupled_bacc", "decoupled_mcc",
            "coupled_sens", "coupled_spec", "coupled_bacc", "coupled_mcc",
            "extra_false_negatives_from_coupling", "fp_reduction_from_coupling"
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(csv_rows)
        
    out_json = "benchmarks/heldout/results/bbb_sensitivity_summary.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print(f"\nSaved sensitivity analysis to {out_csv} and {out_json}")

if __name__ == "__main__":
    main()
