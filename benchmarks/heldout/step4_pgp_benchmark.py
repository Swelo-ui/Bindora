import os
import sys
sys.path.insert(0, os.path.abspath("."))
import csv
import json
import math
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem
from backend.services.adme import predict_pgp_substrate

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

def calc_metrics(y_true, y_pred, y_score=None):
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    bacc = (sens + spec) / 2.0
    
    denom = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc = (tp * tn - fp * fn) / math.sqrt(denom) if denom > 0 else 0.0
    
    if y_score is not None:
        try:
            auc = roc_auc_score(y_true, y_score)
        except Exception:
            auc = bacc
    else:
        auc = bacc
        
    return {
        "TP": int(tp), "FP": int(fp), "TN": int(tn), "FN": int(fn),
        "Sensitivity": round(sens, 4),
        "Sensitivity_95CI": wilson_ci(tp, tp + fn),
        "Specificity": round(spec, 4),
        "Specificity_95CI": wilson_ci(tn, tn + fp),
        "Balanced_Accuracy": round(bacc, 4),
        "MCC": round(mcc, 4),
        "ROC_AUC": round(auc, 4)
    }

def get_morgan_fp(mol, radius=2, n_bits=1024):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.int8)
    for bit in fp.GetOnBits():
        arr[bit] = 1
    return arr

def main():
    os.makedirs("benchmarks/heldout/results", exist_ok=True)
    
    train_file = "benchmarks/heldout/train.csv"
    test_file = "benchmarks/heldout/test.csv"
    
    # Load Train
    with open(train_file, "r", encoding="utf-8") as f:
        train_rows = list(csv.DictReader(f))
    train_y = np.array([int(r["Y"]) for r in train_rows])
    train_mols = [Chem.MolFromSmiles(r["Drug"]) for r in train_rows]
    train_fps = np.array([get_morgan_fp(m) for m in train_mols])
    
    # Load Test
    with open(test_file, "r", encoding="utf-8") as f:
        test_rows = list(csv.DictReader(f))
    test_y = np.array([int(r["Y"]) for r in test_rows])
    test_mols = [Chem.MolFromSmiles(r["Drug"]) for r in test_rows]
    test_fps = np.array([get_morgan_fp(m) for m in test_mols])
    
    print(f"Dataset loaded: Train={len(train_y)} (pos={sum(train_y)}, neg={len(train_y)-sum(train_y)}), Test={len(test_y)} (pos={sum(test_y)}, neg={len(test_y)-sum(test_y)})")
    
    # 1. Majority Class Baseline
    maj_class = 1 if sum(train_y) >= len(train_y)/2 else 0
    pred_majority = np.full_like(test_y, maj_class)
    
    # 2. Simple Heuristic Rule (MW > 400 and LogP > 3)
    pred_heuristic = []
    for m in test_mols:
        mw = Descriptors.MolWt(m)
        logp = Descriptors.MolLogP(m)
        pred_heuristic.append(1 if (mw > 400.0 and logp > 3.0) else 0)
    pred_heuristic = np.array(pred_heuristic)
    
    # 3. Logistic Regression on Morgan FP
    clf = LogisticRegression(random_state=42, max_iter=1000)
    clf.fit(train_fps, train_y)
    pred_lr = clf.predict(test_fps)
    prob_lr = clf.predict_proba(test_fps)[:, 1]
    
    # 4. Bindora Rule-Based P-gp Substrate Classifier
    pred_bindora = []
    bindora_reasons = []
    bindora_motifs = []
    for m in test_mols:
        res = predict_pgp_substrate(m)
        pred_bindora.append(1 if res["is_substrate"] else 0)
        bindora_reasons.append(res["reason"])
        bindora_motifs.append("; ".join(res["motifs_identified"]))
    pred_bindora = np.array(pred_bindora)
    
    # Subgroup identification: Diphenylmethyl / Piperidine
    p_diphenyl = Chem.MolFromSmarts("[#6](c1ccccc1)(c2ccccc2)")
    p_pip = Chem.MolFromSmarts("N1CCC(CC1)")
    is_diphenyl_or_pip = np.array([
        bool((p_diphenyl and m.HasSubstructMatch(p_diphenyl)) or (p_pip and m.HasSubstructMatch(p_pip)))
        for m in test_mols
    ])
    print(f"Subgroup Diphenylmethyl / Piperidine count in test set: {sum(is_diphenyl_or_pip)} / {len(test_y)}")
    
    # Compute overall metrics
    m_maj = calc_metrics(test_y, pred_majority)
    m_heu = calc_metrics(test_y, pred_heuristic)
    m_lr = calc_metrics(test_y, pred_lr, prob_lr)
    m_bin = calc_metrics(test_y, pred_bindora)
    
    # Compute subgroup metrics
    sub_y = test_y[is_diphenyl_or_pip]
    sub_bin = pred_bindora[is_diphenyl_or_pip]
    sub_lr = pred_lr[is_diphenyl_or_pip]
    sub_lr_prob = prob_lr[is_diphenyl_or_pip]
    sub_heu = pred_heuristic[is_diphenyl_or_pip]
    
    other_y = test_y[~is_diphenyl_or_pip]
    other_bin = pred_bindora[~is_diphenyl_or_pip]
    other_lr = pred_lr[~is_diphenyl_or_pip]
    other_lr_prob = prob_lr[~is_diphenyl_or_pip]
    other_heu = pred_heuristic[~is_diphenyl_or_pip]
    
    m_sub_bin = calc_metrics(sub_y, sub_bin)
    m_sub_lr = calc_metrics(sub_y, sub_lr, sub_lr_prob)
    m_sub_heu = calc_metrics(sub_y, sub_heu)
    
    m_other_bin = calc_metrics(other_y, other_bin)
    m_other_lr = calc_metrics(other_y, other_lr, other_lr_prob)
    m_other_heu = calc_metrics(other_y, other_heu)
    
    # Save test predictions CSV
    out_csv = "benchmarks/heldout/results/pgp_benchmark.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fields = [
            "Drug_ID", "InChIKey", "Y_true", "has_diphenyl_or_piperidine",
            "pred_majority", "pred_heuristic_mw_logp", "pred_logistic_regression", "prob_logistic_regression",
            "pred_bindora", "bindora_motifs", "bindora_reason", "Drug"
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for idx in range(len(test_rows)):
            w.writerow({
                "Drug_ID": test_rows[idx]["Drug_ID"],
                "InChIKey": test_rows[idx]["InChIKey"],
                "Y_true": test_y[idx],
                "has_diphenyl_or_piperidine": int(is_diphenyl_or_pip[idx]),
                "pred_majority": int(pred_majority[idx]),
                "pred_heuristic_mw_logp": int(pred_heuristic[idx]),
                "pred_logistic_regression": int(pred_lr[idx]),
                "prob_logistic_regression": round(float(prob_lr[idx]), 4),
                "pred_bindora": int(pred_bindora[idx]),
                "bindora_motifs": bindora_motifs[idx],
                "bindora_reason": bindora_reasons[idx],
                "Drug": test_rows[idx]["Drug"]
            })
            
    # Save JSON summary
    summary_path = "benchmarks/heldout/results/pgp_metrics_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "Overall_Test_Set": {
                "Majority_Class": m_maj,
                "Heuristic_MW_LogP": m_heu,
                "Logistic_Regression_MorganFP": m_lr,
                "Bindora_Rule_Based": m_bin
            },
            "Subgroup_Diphenylmethyl_Piperidine": {
                "Count": int(sum(is_diphenyl_or_pip)),
                "Heuristic_MW_LogP": m_sub_heu,
                "Logistic_Regression_MorganFP": m_sub_lr,
                "Bindora_Rule_Based": m_sub_bin
            },
            "Subgroup_Other_Scaffolds": {
                "Count": int(sum(~is_diphenyl_or_pip)),
                "Heuristic_MW_LogP": m_other_heu,
                "Logistic_Regression_MorganFP": m_other_lr,
                "Bindora_Rule_Based": m_other_bin
            }
        }, f, indent=2)
        
    print(f"P-gp benchmark complete. Results saved to {out_csv}")
    print("\n=== OVERALL TEST SET PERFORMANCE ===")
    for model_name, metrics in [
        ("Majority Class", m_maj),
        ("Heuristic (MW>400, LogP>3)", m_heu),
        ("Logistic Regression (Morgan FP)", m_lr),
        ("Bindora Rule-Based", m_bin)
    ]:
        print(f"  {model_name:<32} | Sens: {metrics['Sensitivity']:.4f} | Spec: {metrics['Specificity']:.4f} | BAcc: {metrics['Balanced_Accuracy']:.4f} | MCC: {metrics['MCC']:.4f} | AUC: {metrics['ROC_AUC']:.4f}")
        
    print("\n=== SUBGROUP ANALYSIS ===")
    print("Diphenylmethyl / Piperidine Subgroup (N = %d):" % sum(is_diphenyl_or_pip))
    print(f"  Bindora Rule-Based: Sens: {m_sub_bin['Sensitivity']:.4f} | Spec: {m_sub_bin['Specificity']:.4f} | BAcc: {m_sub_bin['Balanced_Accuracy']:.4f} | MCC: {m_sub_bin['MCC']:.4f}")
    print(f"  Logistic Regression: Sens: {m_sub_lr['Sensitivity']:.4f} | Spec: {m_sub_lr['Specificity']:.4f} | BAcc: {m_sub_lr['Balanced_Accuracy']:.4f} | MCC: {m_sub_lr['MCC']:.4f}")
    print("Other Scaffolds (N = %d):" % sum(~is_diphenyl_or_pip))
    print(f"  Bindora Rule-Based: Sens: {m_other_bin['Sensitivity']:.4f} | Spec: {m_other_bin['Specificity']:.4f} | BAcc: {m_other_bin['Balanced_Accuracy']:.4f} | MCC: {m_other_bin['MCC']:.4f}")
    print(f"  Logistic Regression: Sens: {m_other_lr['Sensitivity']:.4f} | Spec: {m_other_lr['Specificity']:.4f} | BAcc: {m_other_lr['Balanced_Accuracy']:.4f} | MCC: {m_other_lr['MCC']:.4f}")

if __name__ == "__main__":
    main()
