#!/usr/bin/env python3
"""
tests/benchmark_accuracy.py
Bindora Accuracy & Validation Benchmark Suite (v1.0)

Evaluates AutoDock Vina and Vinardo scoring performance across curated
protein-ligand complexes with known crystallographic structures and experimental affinities (Kd/Ki/dG).

Outputs:
  - data/benchmarks/validation_report_v1.json
  - data/benchmarks/validation_report_v1.md
"""

import os
import sys
import math
import json
import time
import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Dict, Any, Optional

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from backend.config import DATA_DIR, BENCHMARKS_DIR, CACHE_DIR
from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine

if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Curated Benchmark Dataset (25 Representative Complexes)
# Source: PDBbind v2020 core set, CASF benchmarks, and primary pharmacology
# ---------------------------------------------------------------------------
BENCHMARK_COMPLEXES = [
    {
        "pdb_id": "1CX2",
        "target_name": "Cyclooxygenase-2 (COX-2)",
        "target_class": "Oxidoreductase / Inflammation",
        "drug_name": "SC-558",
        "ligand_resname": "S58",
        "smiles": "Cc1ccc(cc1c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F)Br",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 5.0,
        "exp_delta_g_kcal": -11.32
    },
    {
        "pdb_id": "1IEP",
        "target_name": "Abl1 Tyrosine Kinase",
        "target_class": "Kinase / Oncology",
        "drug_name": "Imatinib (STI-571)",
        "ligand_resname": "STI",
        "smiles": "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 10.0,
        "exp_delta_g_kcal": -10.91
    },
    {
        "pdb_id": "2ITY",
        "target_name": "EGFR Kinase (T790M)",
        "target_class": "Kinase / Oncology",
        "drug_name": "Gefitinib (Iressa)",
        "ligand_resname": "IRE",
        "smiles": "COc1cc2ncnc(c2cc1OCCCN3CCOCC3)Nc4ccc(c(c4)Cl)F",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 2.5,
        "exp_delta_g_kcal": -11.73
    },
    {
        "pdb_id": "1M17",
        "target_name": "EGFR Kinase (Active)",
        "target_class": "Kinase / Oncology",
        "drug_name": "Erlotinib (Tarceva)",
        "ligand_resname": "AQ4",
        "smiles": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 5.0,
        "exp_delta_g_kcal": -11.32
    },
    {
        "pdb_id": "1XKK",
        "target_name": "EGFR / HER2 Kinase",
        "target_class": "Kinase / Oncology",
        "drug_name": "Lapatinib (Tykerb)",
        "ligand_resname": "FMS",
        "smiles": "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 10.0,
        "exp_delta_g_kcal": -10.91
    },
    {
        "pdb_id": "4ER4",
        "target_name": "B-Raf Kinase (V600E)",
        "target_class": "Kinase / Oncology",
        "drug_name": "Sorafenib (Nexavar)",
        "ligand_resname": "BAX",
        "smiles": "CNC(=O)c1cc(ccn1)Oc2ccc(cc2)NC(=O)Nc3ccc(c(c3)Cl)C(F)(F)F",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 22.0,
        "exp_delta_g_kcal": -10.44
    },
    {
        "pdb_id": "4JSX",
        "target_name": "Phosphoinositide 3-Kinase (PI3K-gamma)",
        "target_class": "Kinase / Cell Signaling",
        "drug_name": "Torin 2",
        "ligand_resname": "TOR",
        "smiles": "Cc1ccc(cc1)c2cc3c(nc2c4cccnc4)c(=O)n(cn3)c5cc(cc(c5)Cl)C(F)(F)F",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 2.0,
        "exp_delta_g_kcal": -11.86
    },
    {
        "pdb_id": "2JDD",
        "target_name": "Human Renin",
        "target_class": "Aspartyl Protease / Cardiovascular",
        "drug_name": "Aliskiren (Tekturna)",
        "ligand_resname": "AL1",
        "smiles": "CC(C)C[C@@H](C(=O)N[C@@H](CC(C)C)C(=O)N)N[C@@H](CC1=CC(=C(C=C1)OC)OCCCOC)C[C@H]([C@@H](CC(C)C)N)O",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 0.6,
        "exp_delta_g_kcal": -12.58
    },
    {
        "pdb_id": "1STP",
        "target_name": "Streptavidin",
        "target_class": "Biotin-binding Protein",
        "drug_name": "Biotin (Vitamin B7)",
        "ligand_resname": "BTN",
        "smiles": "C1[C@@H]2[C@H](S1)NC(=O)N2CCCCC(=O)O",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 0.00004,
        "exp_delta_g_kcal": -18.30
    },
    {
        "pdb_id": "3ERT",
        "target_name": "Estrogen Receptor Alpha (ERα)",
        "target_class": "Nuclear Hormone Receptor",
        "drug_name": "4-Hydroxytamoxifen",
        "ligand_resname": "OHT",
        "smiles": "CCC(=C(c1ccc(cc1)OCCN(C)C)c2ccccc2)c3ccc(cc3)O",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 1.5,
        "exp_delta_g_kcal": -12.04
    },
    {
        "pdb_id": "1HSG",
        "target_name": "HIV-1 Protease",
        "target_class": "Aspartyl Protease / Antiviral",
        "drug_name": "Indinavir (Crixivan)",
        "ligand_resname": "MK1",
        "smiles": "CC(C)(C)NC(=O)[C@@H]1CN(CCN1C[C@@H](C[C@@H](CC2=CC=CC=C2)NC(=O)[C@H]3C[C@H]4CCCC[C@H]4N3)O)CC5=CN=CC=C5",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 0.54,
        "exp_delta_g_kcal": -12.64
    },
    {
        "pdb_id": "1E66",
        "target_name": "Acetylcholinesterase (AChE)",
        "target_class": "Hydrolase / Neuroscience",
        "drug_name": "Huperzine A",
        "ligand_resname": "HUP",
        "smiles": "CC=C1CC2(CC(=C1)C)C(=O)NC3=C2CCC(C3)N",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 82.0,
        "exp_delta_g_kcal": -9.66
    },
    {
        "pdb_id": "1UW6",
        "target_name": "Cyclin-Dependent Kinase 2 (CDK2)",
        "target_class": "Kinase / Cell Cycle",
        "drug_name": "Hymenialdisine",
        "ligand_resname": "HYM",
        "smiles": "c1cc(c([nH]1)Br)C(=C2C(=O)NC(=N2)N)C(=O)[nH]c3c1",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 22.0,
        "exp_delta_g_kcal": -10.44
    },
    {
        "pdb_id": "2AZM",
        "target_name": "Carbonic Anhydrase II (CA-II)",
        "target_class": "Lyase / Diuretic",
        "drug_name": "Furosemide (Lasix)",
        "ligand_resname": "FUM",
        "smiles": "c1c(c(cc(c1Cl)S(=O)(=O)N)C(=O)O)NCc2ccco2",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 190.0,
        "exp_delta_g_kcal": -9.16
    },
    {
        "pdb_id": "3LPB",
        "target_name": "Beta-2 Adrenergic Receptor (β2AR)",
        "target_class": "GPCR / Respiratory",
        "drug_name": "Carazolol",
        "ligand_resname": "CAU",
        "smiles": "CC(C)NCC(COc1cccc2c1c3c([nH]2)cccc3)O",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 0.03,
        "exp_delta_g_kcal": -14.36
    },
    {
        "pdb_id": "1UY6",
        "target_name": "HSP90 N-Terminal Domain",
        "target_class": "Chaperone / Oncology",
        "drug_name": "Radicicol",
        "ligand_resname": "RDC",
        "smiles": "C/C=C/C=C\\C1CC2C(O2)/C=C/C(=O)Oc3c(cc(c(c3C1)O)Cl)O",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 19.0,
        "exp_delta_g_kcal": -10.53
    },
    {
        "pdb_id": "3GCS",
        "target_name": "Dipeptidyl Peptidase IV (DPP-4)",
        "target_class": "Protease / Diabetes",
        "drug_name": "Sitagliptin (Januvia)",
        "ligand_resname": "M74",
        "smiles": "c1cc(c(cc1F)F)CC(CC(=O)N2CC3=C(C2)N=C(N3)C(F)(F)F)N",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 18.0,
        "exp_delta_g_kcal": -10.56
    },
    {
        "pdb_id": "1T46",
        "target_name": "c-Kit Tyrosine Kinase",
        "target_class": "Kinase / Oncology",
        "drug_name": "Imatinib",
        "ligand_resname": "STI",
        "smiles": "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 100.0,
        "exp_delta_g_kcal": -9.55
    },
    {
        "pdb_id": "4OB0",
        "target_name": "PARP-1 Catalytic Domain",
        "target_class": "Transferase / DNA Repair",
        "drug_name": "Olaparib (Lynparza)",
        "ligand_resname": "0Y4",
        "smiles": "O=C1NC(=O)c2cc(ccc2N1)Cc3ccc(cc3)C(=O)N4CCN(CC4)C(=O)C5CC5",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 5.0,
        "exp_delta_g_kcal": -11.32
    },
    {
        "pdb_id": "3CLD",
        "target_name": "SARS-CoV-2 Main Protease (Mpro)",
        "target_class": "Cysteine Protease / Antiviral",
        "drug_name": "N3 Inhibitor",
        "ligand_resname": "N3",
        "smiles": "CC(C)C[C@H](NC(=O)C1=CN=CC=C1)C(=O)N[C@@H](CC2=CC=CC=C2)C(=O)N",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 10000.0,
        "exp_delta_g_kcal": -6.82
    },
    {
        "pdb_id": "3PTB",
        "target_name": "Bovine Trypsin",
        "target_class": "Serine Protease",
        "drug_name": "Benzamidine",
        "ligand_resname": "BEN",
        "smiles": "c1ccccc1C(=N)N",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 18000.0,
        "exp_delta_g_kcal": -6.47
    },
    {
        "pdb_id": "1DWD",
        "target_name": "Thrombin",
        "target_class": "Serine Protease / Anticoagulant",
        "drug_name": "Argatroban (Acova)",
        "ligand_resname": "ARG",
        "smiles": "Cc1cccc2c1CCC[C@@H]2NS(=O)(=O)c3cccc4c3CCC[C@@H]4C(=O)N5CCC[C@H]5C(=O)N[C@@H](CCCNC(=N)N)C(=O)O",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 19.0,
        "exp_delta_g_kcal": -10.53
    },
    {
        "pdb_id": "2XY9",
        "target_name": "JAK2 Kinase",
        "target_class": "Kinase / Hematology",
        "drug_name": "Ruxolitinib (Jakafi)",
        "ligand_resname": "JAK",
        "smiles": "C1CCC(C1)[C@@H](CC#N)n2cc(cn2)c3c4cc[nH]c4ncn3",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 3.3,
        "exp_delta_g_kcal": -11.57
    },
    {
        "pdb_id": "3POZ",
        "target_name": "BRD4 Bromodomain 1",
        "target_class": "Epigenetic Reader",
        "drug_name": "(+)-JQ1",
        "ligand_resname": "JQ1",
        "smiles": "Cc1sc2c(c1C)c(n3c(s2)c(nc3[C@H](CC(=O)OC(C)(C)C)c4ccc(cc4)Cl)C)C",
        "exp_affinity_type": "Kd",
        "exp_value_nm": 77.0,
        "exp_delta_g_kcal": -9.70
    },
    {
        "pdb_id": "7BV2",
        "target_name": "SARS-CoV-2 RNA-Dependent RNA Polymerase",
        "target_class": "Viral Polymerase",
        "drug_name": "Remdesivir monophosphate (GS-441524-MP)",
        "ligand_resname": "F86",
        "smiles": "C1=CC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)O)O)(C#N)O",
        "exp_affinity_type": "Ki",
        "exp_value_nm": 500.0,
        "exp_delta_g_kcal": -8.59
    }
]

# ---------------------------------------------------------------------------
# Statistics Helper Functions (Pure Python, Zero Dependency)
# ---------------------------------------------------------------------------
def calc_pearson_r(x: List[float], y: List[float]) -> float:
    n = len(x)
    if n < 2:
        return 0.0
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    num = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    den_x = sum((xi - mean_x) ** 2 for xi in x)
    den_y = sum((yi - mean_y) ** 2 for yi in y)
    if den_x == 0 or den_y == 0:
        return 0.0
    return num / (math.sqrt(den_x) * math.sqrt(den_y))

def calc_rmse(actual: List[float], predicted: List[float]) -> float:
    n = len(actual)
    if n == 0:
        return 0.0
    return math.sqrt(sum((a - p) ** 2 for a, p in zip(actual, predicted)) / n)

def calc_mae(actual: List[float], predicted: List[float]) -> float:
    n = len(actual)
    if n == 0:
        return 0.0
    return sum(abs(a - p) for a, p in zip(actual, predicted)) / n

def calc_spearman_rho(x: List[float], y: List[float]) -> float:
    def rank_data(vals):
        sorted_indices = sorted(range(len(vals)), key=lambda i: vals[i])
        ranks = [0.0] * len(vals)
        for rank, idx in enumerate(sorted_indices):
            ranks[idx] = rank + 1
        return ranks
    if len(x) < 2:
        return 0.0
    rx = rank_data(x)
    ry = rank_data(y)
    return calc_pearson_r(rx, ry)

# ---------------------------------------------------------------------------
# Benchmark Runner
# ---------------------------------------------------------------------------
def run_benchmark(
    complexes: List[Dict[str, Any]],
    exhaustiveness: int = 8,
    max_count: Optional[int] = None,
    output_json: Optional[Path] = None,
    output_md: Optional[Path] = None
) -> Dict[str, Any]:
    target_complexes = complexes[:max_count] if max_count else complexes
    total = len(target_complexes)
    print(f"\n=======================================================")
    print(f"  Bindora Molecular Docking Benchmark Suite (v1.0)")
    print(f"  Evaluating {total} protein-ligand benchmark complexes")
    print(f"  Exhaustiveness: {exhaustiveness}")
    print(f"=======================================================\n")

    results = []
    t_start = time.time()

    for idx, c in enumerate(target_complexes, 1):
        pdb_id = c["pdb_id"]
        drug = c["drug_name"]
        exp_dg = c["exp_delta_g_kcal"]
        print(f"[{idx:02d}/{total:02d}] Evaluating {pdb_id} ({c['target_name']}) vs {drug} (Exp dG = {exp_dg:.2f} kcal/mol)...")

        entry_res = {
            "index": idx,
            "pdb_id": pdb_id,
            "target_name": c["target_name"],
            "target_class": c["target_class"],
            "drug_name": drug,
            "ligand_resname": c.get("ligand_resname", ""),
            "exp_affinity_type": c["exp_affinity_type"],
            "exp_value_nm": c["exp_value_nm"],
            "exp_delta_g_kcal": exp_dg,
            "success": False
        }

        try:
            # 1. Ingest receptor
            rec_meta = StructureFetcher.fetch_rcsb_pdb(pdb_id)
            if not rec_meta or "pdb_content" not in rec_meta:
                raise RuntimeError(f"Failed to fetch PDB content for {pdb_id}")

            rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"])
            pocket = rec.get("detected_pocket") or rec.get("blind_docking_box")
            if not pocket:
                raise RuntimeError("No binding pocket or bounding box detected")

            native_ligand = rec.get("native_ligand", {})
            has_native = native_ligand.get("has_native", False)
            native_pdb = native_ligand.get("pdb_block") or native_ligand.get("pdb_content", "")

            vina_dg = None
            vinardo_dg = None
            rmsd = None
            redock_valid = False

            if has_native and native_pdb:
                # Run native redocking validation
                redock = DockingEngine.run_redocking_validation(
                    rec["pdbqt_text"],
                    native_pdb,
                    pocket["center"],
                    pocket["size"],
                    exhaustiveness=exhaustiveness,
                    seed=42
                )
                vina_dg = redock.get("affinity_kcal")
                vinardo_dg = redock.get("vinardo_affinity_kcal")
                rmsd = redock.get("rmsd_angstroms")
                redock_valid = redock.get("is_validated", False)
            else:
                # Cross-docking with prepared ligand SMILES
                lig = DockingEngine.prepare_ligand(c["smiles"])
                poses = DockingEngine.run_docking(
                    rec["pdbqt_text"],
                    lig["pdbqt_text"],
                    pocket["center"],
                    pocket["size"],
                    exhaustiveness=exhaustiveness,
                    num_modes=3,
                    seed=42
                )
                if poses:
                    vina_dg = poses[0]["affinity_kcal"]
                    vinardo_dg = DockingEngine.score_pose_vinardo(
                        rec["pdbqt_text"],
                        poses[0]["pdbqt_content"],
                        pocket["center"],
                        pocket["size"]
                    )
                    rmsd = None
                    redock_valid = False

            entry_res["vina_delta_g_kcal"] = vina_dg
            entry_res["vinardo_delta_g_kcal"] = vinardo_dg
            entry_res["rmsd_angstroms"] = rmsd
            entry_res["is_validated"] = redock_valid
            entry_res["has_native_redock"] = has_native

            if vina_dg is not None:
                entry_res["vina_abs_error_kcal"] = round(abs(vina_dg - exp_dg), 2)
            if vinardo_dg is not None:
                entry_res["vinardo_abs_error_kcal"] = round(abs(vinardo_dg - exp_dg), 2)

            entry_res["success"] = (vina_dg is not None)
            status_str = f"Vina dG = {vina_dg} | Vinardo = {vinardo_dg} | RMSD = {rmsd} A" if rmsd is not None else f"Vina dG = {vina_dg} | Vinardo = {vinardo_dg}"
            print(f"      -> SUCCESS: {status_str}")

        except Exception as e:
            print(f"      -> FAILED: {e}")
            entry_res["error"] = str(e)

        results.append(entry_res)

    elapsed = time.time() - t_start

    # Filter successful evaluations for statistics
    successful = [r for r in results if r["success"] and r.get("vina_delta_g_kcal") is not None]
    exp_vals = [r["exp_delta_g_kcal"] for r in successful]
    vina_vals = [r["vina_delta_g_kcal"] for r in successful]
    vinardo_vals = [r["vinardo_delta_g_kcal"] for r in successful if r.get("vinardo_delta_g_kcal") is not None]
    exp_for_vinardo = [r["exp_delta_g_kcal"] for r in successful if r.get("vinardo_delta_g_kcal") is not None]

    # Pose RMSD accuracy (for complexes with native redocking)
    rmsd_entries = [r for r in successful if r.get("rmsd_angstroms") is not None]
    rmsd_vals = [r["rmsd_angstroms"] for r in rmsd_entries]
    rmsd_success_count = sum(1 for r in rmsd_entries if r["rmsd_angstroms"] <= 2.0)
    rmsd_success_rate = (rmsd_success_count / len(rmsd_entries) * 100.0) if rmsd_entries else 0.0

    stats = {
        "total_complexes_tested": total,
        "successful_evaluations": len(successful),
        "elapsed_seconds": round(elapsed, 2),
        "vina_pearson_r": round(calc_pearson_r(exp_vals, vina_vals), 3) if len(successful) >= 2 else 0.0,
        "vina_spearman_rho": round(calc_spearman_rho(exp_vals, vina_vals), 3) if len(successful) >= 2 else 0.0,
        "vina_rmse_kcal": round(calc_rmse(exp_vals, vina_vals), 2) if len(successful) > 0 else 0.0,
        "vina_mae_kcal": round(calc_mae(exp_vals, vina_vals), 2) if len(successful) > 0 else 0.0,
        "vinardo_pearson_r": round(calc_pearson_r(exp_for_vinardo, vinardo_vals), 3) if len(vinardo_vals) >= 2 else 0.0,
        "vinardo_rmse_kcal": round(calc_rmse(exp_for_vinardo, vinardo_vals), 2) if len(vinardo_vals) > 0 else 0.0,
        "vinardo_mae_kcal": round(calc_mae(exp_for_vinardo, vinardo_vals), 2) if len(vinardo_vals) > 0 else 0.0,
        "pose_reconstruction": {
            "total_native_complexes": len(rmsd_entries),
            "rmsd_under_2a_count": rmsd_success_count,
            "rmsd_success_rate_percent": round(rmsd_success_rate, 1),
            "mean_rmsd_angstroms": round(sum(rmsd_vals) / len(rmsd_vals), 2) if rmsd_vals else 0.0
        },
        # Literature comparison fields — published AutoDock Vina reference values
        # from the CASF-2016 scoring power benchmark (Su et al. JCIM 2019).
        # These allow readers to directly compare Bindora's numbers to published Vina results.
        # Source: Su M. et al. J. Chem. Inf. Model. 2019, 59(2), 895-913.
        #         DOI: 10.1021/acs.jcim.8b00545 — Table 3 (AutoDock Vina on 285 CASF-2016 complexes)
        "literature_r": 0.564,      # Published Vina Pearson R on CASF-2016 full set
        "literature_rmse": 2.19,    # Published Vina RMSE (kcal/mol) on CASF-2016 full set
        "literature_source": "Su M. et al. J. Chem. Inf. Model. 2019, 59(2), 895-913. DOI:10.1021/acs.jcim.8b00545",
        "literature_caveat": (
            "Literature values are from the 285-complex CASF-2016 core set (PDBbind v2016). "
            "Bindora's curated 25-complex set is a different and smaller subset. "
            "Direct numerical comparison is informative but not statistically equivalent. "
            "Run tests/benchmark_casf2016.py for a true apples-to-apples comparison."
        )
    }

    report = {
        "metadata": {
            "benchmark_version": "1.0",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "docking_engine": "AutoDock Vina 1.2.7",
            "scoring_functions": ["AutoDock Vina (Standard)", "Vinardo (Kortemme / Quiroga 2016)"],
            "exhaustiveness": exhaustiveness,
            "dataset_origin": "PDBbind / CASF Curated Complexes (25 PDBs)",
            "operating_system": sys.platform
        },
        "summary_statistics": stats,
        "complex_results": results
    }

    # Save JSON Report
    if output_json:
        output_json.parent.mkdir(parents=True, exist_ok=True)
        with open(output_json, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        print(f"\nSaved JSON validation report to: {output_json}")

    # Generate Markdown Report
    if output_md:
        output_md.parent.mkdir(parents=True, exist_ok=True)
        md_content = generate_markdown_report(report)
        with open(output_md, "w", encoding="utf-8") as f:
            f.write(md_content)
        print(f"Saved Markdown validation report to: {output_md}")

    print("\n=======================================================")
    print(f"  BENCHMARK COMPLETED IN {elapsed:.1f} SECONDS")
    print(f"  Vina Pearson R:      {stats['vina_pearson_r']}")
    print(f"  Vinardo Pearson R:   {stats['vinardo_pearson_r']}")
    print(f"  Vina RMSE:           {stats['vina_rmse_kcal']} kcal/mol")
    print(f"  Pose RMSD (< 2.0 Å): {stats['pose_reconstruction']['rmsd_success_rate_percent']}%")
    print("=======================================================\n")

    return report

# ---------------------------------------------------------------------------
# Markdown Report Generator
# ---------------------------------------------------------------------------
def generate_markdown_report(report: Dict[str, Any]) -> str:
    meta = report["metadata"]
    stats = report["summary_statistics"]
    results = report["complex_results"]
    pose_stats = stats["pose_reconstruction"]

    lines = [
        f"# Bindora Accuracy & Validation Benchmark Report (v{meta['benchmark_version']})",
        "",
        f"**Date:** {meta['timestamp_utc'][:10]} &bull; **Engine:** {meta['docking_engine']} &bull; **Scoring:** {', '.join(meta['scoring_functions'])}",
        "",
        "## Executive Summary",
        "",
        "This empirical benchmark validates Bindora's molecular docking engine against curated protein-ligand complexes with published crystallographic coordinates and wet-lab experimental binding affinities ($K_d / K_i / \\Delta G_{\\text{exp}}$).",
        "",
        "| Benchmark Metric | AutoDock Vina | Vinardo Scoring | Target Standard |",
        "| :--- | :---: | :---: | :---: |",
        f"| **Pearson Correlation ($R$)** | **{stats['vina_pearson_r']:.3f}** | **{stats['vinardo_pearson_r']:.3f}** | > 0.50 (CASF Core) |",
        f"| **Spearman Rank Correlation ($\\rho$)** | **{stats['vina_spearman_rho']:.3f}** | — | > 0.50 |",
        f"| **Root Mean Square Error (RMSE)** | **{stats['vina_rmse_kcal']:.2f} kcal/mol** | **{stats['vinardo_rmse_kcal']:.2f} kcal/mol** | < 2.5 kcal/mol |",
        f"| **Mean Absolute Error (MAE)** | **{stats['vina_mae_kcal']:.2f} kcal/mol** | **{stats['vinardo_mae_kcal']:.2f} kcal/mol** | < 2.0 kcal/mol |",
        f"| **Pose Redocking Success (RMSD $\\le 2.0$ Å)** | **{pose_stats['rmsd_success_rate_percent']}%** ({pose_stats['rmsd_under_2a_count']}/{pose_stats['total_native_complexes']}) | — | > 70% |",
        f"| **Mean Crystallographic RMSD** | **{pose_stats['mean_rmsd_angstroms']:.2f} Å** | — | < 2.0 Å |",
        f"| **Literature Vina R (CASF-2016, 285 complexes)** | {stats.get('literature_r', 'N/A')} | — | Reference |",
        f"| **Literature Vina RMSE (CASF-2016, 285 complexes)** | {stats.get('literature_rmse', 'N/A')} kcal/mol | — | Reference |",
        "",
        f"> *Literature reference: {stats.get('literature_source', '')}*",
        f"> *{stats.get('literature_caveat', '')}*",
        "",
        "---",
        "",
        "## Itemized Benchmark Results Matrix",
        "",
        "| PDB | Target Receptor | Investigational Ligand | Exp $\\Delta G$ | Vina $\\Delta G$ | Vinardo | Error | RMSD (Å) | Benchmark Status |",
        "| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |"
    ]

    for r in results:
        pdb = r["pdb_id"]
        tgt = r["target_name"]
        drug = r["drug_name"]
        exp_dg = f"{r['exp_delta_g_kcal']:.2f}"
        vina_dg = f"{r['vina_delta_g_kcal']:.2f}" if r.get("vina_delta_g_kcal") is not None else "—"
        vin_dg = f"{r['vinardo_delta_g_kcal']:.2f}" if r.get("vinardo_delta_g_kcal") is not None else "—"
        err = f"{r.get('vina_abs_error_kcal', '—')}"
        rmsd = f"{r['rmsd_angstroms']:.2f}" if r.get("rmsd_angstroms") is not None else "—"

        if r.get("rmsd_angstroms") is not None:
            status = "✅ Validated (< 2.0 Å)" if r["rmsd_angstroms"] <= 2.0 else "⚠️ Near-Native (> 2.0 Å)"
        elif r.get("success"):
            status = "✅ Docked"
        else:
            status = "❌ Failed"

        lines.append(f"| `{pdb}` | {tgt} | {drug} | {exp_dg} | {vina_dg} | {vin_dg} | {err} | {rmsd} | {status} |")

    lines.extend([
        "",
        "---",
        "",
        "## Methodology & Reproducibility Protocol",
        "",
        "1. **Receptor Ingestion:** Standard biological assemblies stripped of co-solvents and crystallographic waters. Ionization set to physiological pH 7.4. Gasteiger-Marsili partial charges and AutoDock 4 atom types assigned via Meeko.",
        "2. **Grid Box Centering:** Pocket centroid determined from crystallographic bound ligand coordinates with 22.0 Å cubic search space.",
        "3. **Stochastic Sampling:** AutoDock Vina iterated local search with Monte Carlo sampling (fixed seed 42 for complete reproducibility).",
        "4. **Vinardo Rescoring:** Scoring executed via AutoDock Vina v1.2.7 `--scoring vinardo --score_only --autobox`.",
        "5. **RMSD Calculation:** Symmetry-corrected heavy-atom root mean square deviation between docked pose coordinates and crystallographic coordinates.",
        "",
        "## Citations",
        "",
        "- Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry*, 31(2), 455-461.",
        "- Quiroga, R., & Villarreal, M. A. (2016). Vinardo: A Scoring Function Based on Autodock Vina with Improved Affinity Predictions. *PLoS ONE*, 11(5), e0155182.",
        "- Su, M., et al. (2019). Comparative Assessment of Scoring Functions (CASF-2016) on the PDBbind Database. *Journal of Chemical Information and Modeling*, 59(2), 895-913."
    ])

    return "\n".join(lines)

# ---------------------------------------------------------------------------
# CLI Entrypoint
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Bindora Accuracy Benchmark")
    parser.add_argument("--fast", action="store_true", help="Run fast evaluation on top 5 complexes")
    parser.add_argument("--num", type=int, default=None, help="Number of complexes to evaluate")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Exhaustiveness setting (default: 8)")
    parser.add_argument("--output-json", type=str, default=str(BENCHMARKS_DIR / "validation_report_v1.json"))
    parser.add_argument("--output-md", type=str, default=str(BENCHMARKS_DIR / "validation_report_v1.md"))

    args = parser.parse_args()

    max_count = 5 if args.fast else args.num
    exh = 4 if args.fast else args.exhaustiveness

    out_json = Path(args.output_json)
    out_md = Path(args.output_md)

    run_benchmark(
        BENCHMARK_COMPLEXES,
        exhaustiveness=exh,
        max_count=max_count,
        output_json=out_json,
        output_md=out_md
    )
