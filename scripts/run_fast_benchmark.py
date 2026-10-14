# scripts/run_fast_benchmark.py
import sys
import time
import argparse
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

CORE_BENCHMARK_TARGETS = [
    {
        "pdb_id": "1HSG",
        "ligand_resname": "MK1",
        "family": "Aspartyl Protease",
        "name": "HIV-1 Protease Homodimer",
        "smiles": "CC(C)(C)NC(=O)[C@@H]1CN(Cc2cccnc2)CCN1C[C@@H](O)C[C@@H](Cc1ccccc1)C(=O)N[C@@H]1[C@H](O)Cc2ccccc21",
        "target_chain": "ALL",
        "rmsd_threshold": 2.0,
    },
    {
        "pdb_id": "1AQ1",
        "ligand_resname": "STU",
        "family": "Ser/Thr Kinase",
        "name": "Human CDK2 Kinase",
        "smiles": "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)n1c3ccccc3c3c4CNC(=O)c4c4c5ccccc5n2c4c31",
        "target_chain": "A",
        "rmsd_threshold": 2.0,
    },
]

EXTENDED_TARGETS = [
    {
        "pdb_id": "3ERT",
        "ligand_resname": "OHT",
        "family": "Nuclear Receptor",
        "name": "Human Estrogen Receptor Alpha",
        "smiles": "CC/C(=C(\\c1ccccc1)/c1ccc(OCCN(C)C)cc1)/c1ccc(O)cc1",
        "target_chain": "A",
        "rmsd_threshold": 2.0,
    },
    {
        "pdb_id": "1BCD",
        "ligand_resname": "FMS",
        "family": "Zn2+ Metalloenzyme",
        "name": "Human Carbonic Anhydrase II (Zn2+)",
        "smiles": "NS(=O)(=O)c1ccc(F)cc1",
        "target_chain": "A",
        "rmsd_threshold": 2.0,
    },
]


def run_benchmark_suite(extended: bool = False, work_dir: str = "./benchmark_workspace") -> int:
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

    targets = list(CORE_BENCHMARK_TARGETS)
    if extended:
        targets.extend(EXTENDED_TARGETS)

    print("=" * 80)
    print(f" BINDORA DOCK - RESEARCH-GRADE BENCHMARK SUITE ({'EXTENDED 4-FAMILY' if extended else 'FAST 2-TARGET'} VALIDATION)")
    print("=" * 80)

    summary_rows = []
    failures = 0

    for entry in targets:
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
        top_pose = dock_res["top_pose"]
        top_aff = top_pose["affinity_kcal_mol"]
        ki_nm = top_pose["estimated_ki_nM"]
        le = top_pose["ligand_efficiency"]
        passed = (sym_rmsd <= entry["rmsd_threshold"]) and val_res["posebusters_valid"]
        status_str = "PASS" if passed else "FAIL"
        if not passed:
            failures += 1

        summary_rows.append({
            "pdb_id": pdb_id,
            "family": entry.get("family", "Other"),
            "name": entry["name"],
            "lig_id": lig_id,
            "heavy_atoms": lig_data["heavy_atoms"],
            "rotatable_bonds": lig_data["rotatable_bonds"],
            "affinity": top_aff,
            "ki_nm": ki_nm,
            "le": le,
            "rmsd": sym_rmsd,
            "posebusters": val_res["posebusters_valid"],
            "status": status_str,
            "time_s": elapsed,
        })

        print(
            f"  -> Result [{status_str}] | Affinity: {top_aff} kcal/mol | Ki: {ki_nm} nM | LE: {le} | "
            f"Sym-RMSD: {sym_rmsd} Å | PoseBusters: {val_res['posebusters_valid']} | Time: {elapsed}s"
        )
        print(f"  -> Report: {reports['markdown_report']}")

    # Write Markdown Summary Table
    reports_dir = work_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    summary_md_path = reports_dir / "BENCHMARK_SUMMARY.md"

    table_lines = [
        "# Bindora Dock: Research-Grade Redocking Benchmark Summary\n",
        f"**Search Engine:** AutoDock Vinardo (`seed=42`, `exhaustiveness=8`, dynamic $R_g$ box)\n",
        r"| Target Family | PDB ID | Target Protein | Ligand | Torsions | Heavy Atoms | Top-1 Affinity ($\Delta G$) | Estimated $K_i$ | Symmetry RMSD | PoseBusters | Status |",
        "| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for r in summary_rows:
        table_lines.append(
            f"| **{r['family']}** | `{r['pdb_id']}` | {r['name']} | `{r['lig_id']}` | "
            f"{r['rotatable_bonds']} | {r['heavy_atoms']} | **{r['affinity']} kcal/mol** | "
            f"{r['ki_nm']} nM | **{r['rmsd']:.3f} Å** | {'PASS' if r['posebusters'] else 'WARN'} | **{r['status']}** |"
        )
    table_lines.append("\n*All benchmarks evaluated ab initio strictly from 1D SMILES with zero coordinate leakage.*")
    summary_md_path.write_text("\n".join(table_lines) + "\n", encoding="utf-8")
    print(f"\n[SUMMARY] Formatted benchmark table saved to: {summary_md_path}")

    print("\n" + "=" * 80)
    if failures > 0:
        print(f"BENCHMARK SUITE COMPLETED WITH {failures} WARNING/FAILURE(S).")
        return 1
    print(f"ALL {len(targets)} BENCHMARK TARGET(S) PASSED RESEARCH-GRADE THRESHOLDS (RMSD <= 2.0 Å).")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bindora Dock Fast & Extended Benchmark Suite")
    parser.add_argument(
        "--extended",
        action="store_true",
        help="Run extended multi-family suite (Protease, Kinase, Nuclear Receptor 3ERT, Zn2+ Metalloenzyme 1BCD)",
    )
    parser.add_argument(
        "--work-dir",
        default="./benchmark_workspace",
        help="Workspace directory path for temporary pdb/pdbqt and outputs",
    )
    args = parser.parse_args()
    sys.exit(run_benchmark_suite(extended=args.extended, work_dir=args.work_dir))
