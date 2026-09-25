import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import time
import json
import math
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine
from backend.services.bioactivity import BioactivityService
from backend.utils.rmsd_calculator import calculate_rmsd

def run_1hsg_benchmark():
    print("=" * 60)
    print("STARTING 1HSG (HIV-1 Protease / Indinavir MK1) BENCHMARK")
    print("=" * 60)
    t0 = time.time()
    
    # 1. Fetch
    rec_meta = StructureFetcher.fetch_rcsb_pdb("1HSG")
    pdb_content = rec_meta["pdb_content"]
    
    # 2. Receptor Preparation (preserve C2 homodimer: Chains A & B)
    rec = DockingEngine.prepare_receptor(pdb_content)
    native = rec["native_ligand"]
    pocket_center = rec["detected_pocket"]["center"]
    pocket_size = rec["detected_pocket"]["size"]
    
    print(f"Receptor: HIV-1 Protease C2 Homodimer (Chains: {rec.get('chains', ['A', 'B'])})")
    print(f"Native Ligand: {native['name']} ({native['atom_count']} atoms)")
    print(f"Grid Center: x={pocket_center['x']:.2f}, y={pocket_center['y']:.2f}, z={pocket_center['z']:.2f}")
    print(f"Grid Box Size: x={pocket_size['x']:.2f}, y={pocket_size['y']:.2f}, z={pocket_size['z']:.2f}")
    
    # 3. Native Ligand Preparation (MK1)
    nat_prep = DockingEngine.prepare_native_ligand(native["pdb_block"])
    print(f"Prepared Heavy Atoms: {nat_prep['heavy_atom_count']}")
    print(f"Prepared Rotatable Bonds: {nat_prep['rotatable_bonds']}")
    
    # 4. Docking Execution
    seed = 42
    exhaustiveness = 8
    print(f"Running Vina (exhaustiveness={exhaustiveness}, seed={seed}, num_modes=9)...")
    td0 = time.time()
    poses = DockingEngine.run_docking(
        rec["pdbqt_text"],
        nat_prep["pdbqt_text"],
        pocket_center,
        pocket_size,
        exhaustiveness=exhaustiveness,
        num_modes=9,
        seed=seed
    )
    tdock = time.time() - td0
    print(f"Docking completed in {tdock:.2f} seconds. Poses generated: {len(poses)}")
    
    # 5. RMSD calculation against crystallographic reference
    cryst_ref_block = nat_prep.get("pdb_block", native["pdb_block"])
    
    pose_results = []
    min_rmsd = float("inf")
    best_mode = None
    
    for p in poses:
        rmsd_info = calculate_rmsd(p["pdbqt_content"], cryst_ref_block, return_details=True)
        rmsd_val = rmsd_info["rmsd"] if isinstance(rmsd_info, dict) else rmsd_info
        method = rmsd_info.get("method", "unknown") if isinstance(rmsd_info, dict) else "direct"
        autos = rmsd_info.get("automorphisms_tested", 1) if isinstance(rmsd_info, dict) else 1
        
        thermo = BioactivityService.calculate_thermodynamics(p["affinity_kcal"], nat_prep["heavy_atom_count"], 613.79)
        
        mode_data = {
            "mode": p["mode"],
            "affinity_kcal": p["affinity_kcal"],
            "derived_kd_nm": thermo["theoretical_kd_nm"],
            "ligand_efficiency": thermo["ligand_efficiency"]["value"],
            "rmsd_to_cryst": rmsd_val,
            "mapping_method": method,
            "automorphisms_tested": autos
        }
        pose_results.append(mode_data)
        if rmsd_val < min_rmsd:
            min_rmsd = rmsd_val
            best_mode = p["mode"]
            
        print(f"  Mode {p['mode']}: Score = {p['affinity_kcal']:6.2f} kcal/mol | "
              f"Kd-like = {thermo['theoretical_kd_nm']:8.1f} nM | "
              f"LE = {thermo['ligand_efficiency']['value']:.3f} | "
              f"Crystal RMSD = {rmsd_val:5.2f} Å ({method}, {autos} autos)")

    rank1 = pose_results[0]
    is_pass = rank1["rmsd_to_cryst"] <= 2.0 or min_rmsd <= 2.0
    status = "PASS" if is_pass else "FAIL"
    
    result_1hsg = {
        "system": "1HSG (HIV-1 Protease C2 Homodimer / Indinavir MK1)",
        "chains": rec.get("chains", ["A", "B"]),
        "ligand": "Indinavir (MK1)",
        "heavy_atoms": nat_prep["heavy_atom_count"],
        "rotatable_bonds": nat_prep["rotatable_bonds"],
        "grid_center": pocket_center,
        "grid_size": pocket_size,
        "seed": seed,
        "exhaustiveness": exhaustiveness,
        "num_modes": len(poses),
        "rank_1_score_kcal": rank1["affinity_kcal"],
        "rank_1_rmsd_angstroms": rank1["rmsd_to_cryst"],
        "rank_1_kd_like_nm": rank1["derived_kd_nm"],
        "rank_1_le": rank1["ligand_efficiency"],
        "best_mode": best_mode,
        "best_mode_rmsd_angstroms": min_rmsd,
        "atom_mapping_method": rank1["mapping_method"],
        "automorphisms_tested": rank1["automorphisms_tested"],
        "status": status,
        "docking_time_s": round(tdock, 2),
        "total_time_s": round(time.time() - t0, 2),
        "all_modes": pose_results
    }
    
    return result_1hsg

if __name__ == "__main__":
    res = run_1hsg_benchmark()
    with open("benchmark_1hsg_output.json", "w") as f:
        json.dump(res, f, indent=2)
    print("\n1HSG Benchmark Summary JSON saved to benchmark_1hsg_output.json")
