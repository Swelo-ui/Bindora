# scripts/run_fast_benchmark.py
import sys
import time
from pathlib import Path

# Ensure root package is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bindora.fetcher import PDBFetcherAndSplitter
from bindora.prep_engine import BiophysicalPrepEngine
from bindora.docking import BindoraDockingEngine
from bindora.report_emitter import (
    compute_symmetry_corrected_rmsd,
    check_posebusters_and_biophysics,
    emit_research_report,
)

BENCHMARK_TARGETS = [
    {
        "pdb_id": "1HSG",
        "ligand_resname": "MK1",
        "name": "HIV-1 Protease (Indinavir)",
        "smiles": "CC(C)(C)NC(=O)[C@@H]1CN(Cc2cccnc2)CCN1C[C@@H](O)C[C@@H](Cc1ccccc1)C(=O)N[C@@H]1[C@H](O)Cc2ccccc21",
        "target_chain": "ALL",
        "rmsd_threshold": 2.0,
    },
    {
        "pdb_id": "1AQ1",
        "ligand_resname": "STU",
        "name": "Human CDK2 Kinase (Staurosporine)",
        "smiles": "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)n1c3ccccc3c3c4CNC(=O)c4c4c5ccccc5n2c4c31",
        "target_chain": "A",
        "rmsd_threshold": 2.0,
    },
]


def run_benchmark_suite(work_dir: str = "./benchmark_workspace") -> int:
    work_path = Path(work_dir)
    fetcher = PDBFetcherAndSplitter(work_dir=str(work_path))
    prep = BiophysicalPrepEngine(work_dir=str(work_path), target_ph=7.4, seed=42)
    docker = BindoraDockingEngine(
        work_dir=str(work_path),
        scoring_function="vinardo",
        exhaustiveness=8,
        num_modes=9,
        seed=42,
    )

    print("=" * 78)
    print(" BINDORA DOCK - RESEARCH-GRADE FAST BENCHMARK SUITE (REDOCKING VALIDATION)")
    print("=" * 78)

    failures = 0
    for entry in BENCHMARK_TARGETS:
        t0 = time.time()
        pdb_id = entry["pdb_id"]
        lig_id = entry["ligand_resname"]
        print(f"\n[1/4] Fetching & splitting {pdb_id} ({entry['name']})...")

        raw_pdb = fetcher.fetch_pdb(pdb_id)
        split_data = fetcher.split_receptor_and_ligand(
            raw_pdb,
            ligand_resname=lig_id,
            target_chain=entry["target_chain"],
            retain_metals=True,
            retain_cofactors=True,
        )

        print(f"[2/4] Preparing Receptor & Ligand ({lig_id}) at pH 7.4...")
        rec_pdbqt, rec_log = prep.prepare_receptor_pdbqt(split_data["receptor_pdb"])
        lig_data = prep.prepare_ligand_from_smiles(entry["smiles"], ligand_name=f"{pdb_id}_{lig_id}")

        print(
            f"[3/4] Running Deterministic Vinardo Docking "
            f"(Center={split_data['pocket_center']}, Box={lig_data['optimal_box_size']})..."
        )
        dock_res = docker.dock(
            receptor_pdbqt=rec_pdbqt,
            ligand_pdbqt=lig_data["pdbqt_path"],
            center=split_data["pocket_center"],
            box_size=lig_data["optimal_box_size"],
            heavy_atoms=lig_data["heavy_atoms"],
        )

        print(f"[4/4] Evaluating Symmetry-Corrected RMSD & PoseBusters Validity...")
        top_mol = dock_res["pose_mols"][0]
        sym_rmsd = compute_symmetry_corrected_rmsd(
            top_mol, split_data["crystal_ligand_pdb"], reference_smiles=entry["smiles"]
        )
        val_res = check_posebusters_and_biophysics(
            top_mol, split_data["receptor_pdb"], split_data["crystal_ligand_pdb"]
        )

        reports = emit_research_report(
            pdb_id=pdb_id,
            ligand_id=lig_id,
            prep_info=lig_data,
            dock_info=dock_res,
            validation_info=val_res,
            rmsd_top1=sym_rmsd,
            output_dir=work_path / "reports",
        )

        elapsed = round(time.time() - t0, 1)
        top_aff = dock_res["top_pose"]["affinity_kcal_mol"]
        passed = (sym_rmsd <= entry["rmsd_threshold"]) and val_res["posebusters_valid"]
        status_str = "PASS" if passed else "FAIL"
        if not passed:
            failures += 1

        print(
            f"  -> Result [{status_str}] | Affinity: {top_aff} kcal/mol | "
            f"Sym-RMSD: {sym_rmsd} Å | PoseBusters: {val_res['posebusters_valid']} | "
            f"Time: {elapsed}s"
        )
        print(f"  -> Report saved: {reports['markdown_report']}")

    print("\n" + "=" * 78)
    if failures > 0:
        print(f"BENCHMARK SUITE COMPLETED WITH {failures} WARNING/FAILURE(S).")
        return 1
    print("ALL BENCHMARK TARGETS PASSED RESEARCH-GRADE THRESHOLDS (RMSD <= 2.0 Å).")
    return 0


if __name__ == "__main__":
    sys.exit(run_benchmark_suite())
