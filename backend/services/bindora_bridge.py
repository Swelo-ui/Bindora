# backend/services/bindora_bridge.py
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any

from bindora.fetcher import PDBFetcherAndSplitter
from bindora.prep_engine import BiophysicalPrepEngine
from bindora.docking import BindoraDockingEngine
from bindora.report_emitter import (
    compute_symmetry_corrected_rmsd,
    check_posebusters_and_biophysics,
    emit_research_report,
)


def run_bindora_scientific_pipeline(
    pdb_id: str,
    smiles: str,
    ligand_name: str = "LIG",
    co_crystal_resname: Optional[str] = None,
    active_site_residues: Optional[List[str]] = None,
    custom_center: Optional[Tuple[float, float, float]] = None,
    custom_box_size: Optional[Tuple[float, float, float]] = None,
    scoring_function: str = "vinardo",
    exhaustiveness: int = 8,
    target_ph: float = 7.4,
    seed: int = 42,
    work_dir: str = "./workspace",
) -> Dict[str, Any]:
    """
    Unified production entrypoint connecting the FastAPI backend to the
    research-grade Bindora biophysical docking and validation engine.
    """
    work_path = Path(work_dir)
    fetcher = PDBFetcherAndSplitter(work_dir=str(work_path))
    prep = BiophysicalPrepEngine(work_dir=str(work_path), target_ph=target_ph, seed=seed)
    docker = BindoraDockingEngine(
        work_dir=str(work_path),
        scoring_function=scoring_function,
        exhaustiveness=exhaustiveness,
        num_modes=9,
        seed=seed,
    )

    # 1. Fetch & Split Receptor / Pocket
    raw_pdb = fetcher.fetch_pdb(pdb_id)
    if co_crystal_resname:
        split_data = fetcher.split_receptor_and_ligand(
            raw_pdb,
            ligand_resname=co_crystal_resname,
            target_chain="ALL",
            retain_metals=True,
            retain_cofactors=True,
        )
    else:
        split_data = fetcher.extract_apo_or_residue_pocket(
            raw_pdb,
            active_site_residues=active_site_residues,
            target_chain="ALL",
        )

    # 2. Biophysical Receptor & Ligand Preparation
    rec_pdbqt, rec_prep_log = prep.prepare_receptor_pdbqt(split_data["receptor_pdb"])
    lig_data = prep.prepare_ligand_from_smiles(smiles, ligand_name=f"{pdb_id}_{ligand_name}")

    center = custom_center if custom_center is not None else split_data["pocket_center"]
    box_size = custom_box_size if custom_box_size is not None else lig_data["optimal_box_size"]

    # 3. Deterministic Docking Execution
    dock_res = docker.dock(
        receptor_pdbqt=rec_pdbqt,
        ligand_pdbqt=lig_data["pdbqt_path"],
        center=center,
        box_size=box_size,
        heavy_atoms=lig_data["heavy_atoms"],
    )

    # 4. Biophysical Validation (PoseBusters + Sym-RMSD if crystal reference exists)
    top_mol = dock_res["pose_mols"][0]
    sym_rmsd = None
    if split_data.get("crystal_ligand_pdb") is not None:
        sym_rmsd = compute_symmetry_corrected_rmsd(
            top_mol, split_data["crystal_ligand_pdb"], reference_smiles=smiles
        )

    val_res = check_posebusters_and_biophysics(
        top_mol, split_data["receptor_pdb"], split_data.get("crystal_ligand_pdb")
    )

    reports = emit_research_report(
        pdb_id=pdb_id,
        ligand_id=ligand_name,
        prep_info=lig_data,
        dock_info=dock_res,
        validation_info=val_res,
        rmsd_top1=sym_rmsd,
        output_dir=work_path / "reports",
    )

    return {
        "status": "success",
        "pdb_id": pdb_id,
        "ligand_id": ligand_name,
        "top_pose": dock_res["top_pose"],
        "all_poses": dock_res["poses"],
        "symmetry_rmsd_A": sym_rmsd,
        "posebusters_validation": val_res,
        "grid_center": center,
        "grid_box_size_A": box_size,
        "radius_of_gyration_A": lig_data["radius_of_gyration"],
        "receptor_prep_log": rec_prep_log,
        "artifacts": {
            "receptor_pdb": str(split_data["receptor_pdb"]),
            "docked_sdf": str(dock_res["docked_sdf"]),
            "docked_pdbqt": str(dock_res["docked_pdbqt"]),
            "json_report": str(reports["json_report"]),
            "methods_markdown": str(reports["markdown_report"]),
        },
    }
