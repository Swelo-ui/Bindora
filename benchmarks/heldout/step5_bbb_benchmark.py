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

def calc_confusion_metrics(tp, fp, tn, fn):
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

def main():
    os.makedirs("benchmarks/heldout/results", exist_ok=True)
    in_file = "benchmarks/heldout/bbb_test.csv"
    
    with open(in_file, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
        
    print(f"Loaded {len(rows)} molecules from {in_file}")
    
    records = []
    # Counters for Decoupled (BOILED-Egg passive)
    dec_tp = dec_fp = dec_tn = dec_fn = 0
    # Counters for Coupled (BOILED-Egg passive AND NOT P-gp)
    cpl_tp = cpl_fp = cpl_tn = cpl_fn = 0
    
    extra_false_negatives = []

    for r in rows:
        y_true = int(r["y_binary"])
        smi = r["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
            
        tpsa = round(Descriptors.TPSA(mol, includeSandP=True), 2)
        wlogp = round(Descriptors.MolLogP(mol), 2)
        
        # 1. Passive BOILED-Egg BBB permeation
        inside_egg = _point_in_polygon(tpsa, wlogp, _BBB_COORDS)
        pred_decoupled = 1 if inside_egg else 0
        
        # 2. P-gp active efflux prediction
        pgp_res = predict_pgp_substrate(mol)
        is_pgp = 1 if pgp_res["is_substrate"] else 0
        
        # 3. Coupled rule: Permeant ONLY IF inside egg AND NOT P-gp substrate
        pred_coupled = 1 if (inside_egg and not is_pgp) else 0
        
        # Check if coupling created an extra false negative
        extra_fn = (y_true == 1 and pred_decoupled == 1 and pred_coupled == 0)
        if extra_fn:
            extra_false_negatives.append({
                "compound_name": r["compound_name"],
                "inchikey": r["inchikey"],
                "tpsa": tpsa,
                "wlogp": wlogp,
                "pgp_reason": pgp_res["reason"]
            })
            
        # Update decoupled stats
        if y_true == 1:
            if pred_decoupled == 1: dec_tp += 1
            else: dec_fn += 1
        else:
            if pred_decoupled == 1: dec_fp += 1
            else: dec_tn += 1
            
        # Update coupled stats
        if y_true == 1:
            if pred_coupled == 1: cpl_tp += 1
            else: cpl_fn += 1
        else:
            if pred_coupled == 1: cpl_fp += 1
            else: cpl_tn += 1
            
        records.append({
            "NO": r["NO"],
            "compound_name": r["compound_name"],
            "inchikey": r["inchikey"],
            "y_binary": y_true,
            "bbb_label": r["bbb_label"],
            "tpsa": tpsa,
            "wlogp": wlogp,
            "inside_boiled_egg_bbb": pred_decoupled,
            "is_pgp_substrate": is_pgp,
            "pred_decoupled_current": pred_decoupled,
            "pred_coupled_old": pred_coupled,
            "falsely_excluded_by_coupling": 1 if extra_fn else 0
        })
        
    out_csv = "benchmarks/heldout/results/bbb_benchmark.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "NO", "compound_name", "inchikey", "y_binary", "bbb_label",
            "tpsa", "wlogp", "inside_boiled_egg_bbb", "is_pgp_substrate",
            "pred_decoupled_current", "pred_coupled_old", "falsely_excluded_by_coupling"
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(records)
        
    m_dec = calc_confusion_metrics(dec_tp, dec_fp, dec_tn, dec_fn)
    m_cpl = calc_confusion_metrics(cpl_tp, cpl_fp, cpl_tn, cpl_fn)
    
    summary_path = "benchmarks/heldout/results/bbb_metrics_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "Total_Molecules_Evaluated": len(records),
            "Ground_Truth_BBB_Positive": dec_tp + dec_fn,
            "Ground_Truth_BBB_Negative": dec_tn + dec_fp,
            "Bindora_Current_Decoupled": m_dec,
            "Old_Coupled_Rule": m_cpl,
            "Extra_False_Negatives_From_Coupling": len(extra_false_negatives),
            "Sample_Falsely_Excluded_CNS_Drugs": extra_false_negatives[:20]
        }, f, indent=2)
        
    print(f"BBB benchmark complete. Saved to {out_csv}")
    print("\n=== BBB BENCHMARK COMPARISON ===")
    print("1. Bindora Current (Decoupled Passive BOILED-Egg):")
    print(f"   TP: {m_dec['TP']}, FP: {m_dec['FP']}, TN: {m_dec['TN']}, FN: {m_dec['FN']}")
    print(f"   Sensitivity: {m_dec['Sensitivity']:.4f} (95% CI: {m_dec['Sensitivity_95CI']})")
    print(f"   Specificity: {m_dec['Specificity']:.4f} (95% CI: {m_dec['Specificity_95CI']})")
    print(f"   Balanced Accuracy: {m_dec['Balanced_Accuracy']:.4f}, MCC: {m_dec['MCC']:.4f}")
    
    print("\n2. Old Coupled Rule (Passive BBB AND NOT P-gp):")
    print(f"   TP: {m_cpl['TP']}, FP: {m_cpl['FP']}, TN: {m_cpl['TN']}, FN: {m_cpl['FN']}")
    print(f"   Sensitivity: {m_cpl['Sensitivity']:.4f} (95% CI: {m_cpl['Sensitivity_95CI']})")
    print(f"   Specificity: {m_cpl['Specificity']:.4f} (95% CI: {m_cpl['Specificity_95CI']})")
    print(f"   Balanced Accuracy: {m_cpl['Balanced_Accuracy']:.4f}, MCC: {m_cpl['MCC']:.4f}")
    
    print(f"\nExtra True CNS Drugs Falsely Excluded as Non-Permeant by Coupling: {len(extra_false_negatives)}")
    print(f"Sensitivity dropped by {(m_dec['Sensitivity'] - m_cpl['Sensitivity'])*100:.2f} percentage points due to coupled P-gp overwrite!")

if __name__ == "__main__":
    main()
