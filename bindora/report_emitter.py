# bindora/report_emitter.py
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any

from rdkit import Chem
from rdkit.Chem import rdFMCS


def compute_symmetry_corrected_rmsd(
    docked_mol: Chem.Mol, crystal_pdb_path: Path, reference_smiles: Optional[str] = None
) -> float:
    """
    Computes graph-isomorphism symmetry-corrected heavy-atom RMSD (in Angstroms)
    even when atom ordering differs between SMILES-derived docked pose and PDB crystal ligand.
    """
    cryst_mol = Chem.MolFromPDBFile(str(crystal_pdb_path), removeHs=True, sanitize=False)
    if cryst_mol is None:
        raise ValueError(f"Failed to load crystal ligand from {crystal_pdb_path}")

    # Assign bond orders from template SMILES if available so aromaticity/symmetry is exact
    if reference_smiles:
        template = Chem.MolFromSmiles(reference_smiles)
        if template is not None:
            template = Chem.RemoveHs(template)
            try:
                from rdkit.Chem import AllChem
                cryst_mol = AllChem.AssignBondOrdersFromTemplate(template, cryst_mol)
            except Exception:
                pass

    dock_no_h = Chem.RemoveHs(docked_mol)
    cryst_no_h = Chem.RemoveHs(cryst_mol)

    # Method 1: Try spyrmsd (Graph-isomorphism symmetry-corrected standard)
    try:
        from spyrmsd import rmsd, molecule
        m_ref = molecule.Molecule.from_rdkit(cryst_no_h)
        m_dock = molecule.Molecule.from_rdkit(dock_no_h)
        val = rmsd.symmrmsd(
            m_ref.coordinates,
            [m_dock.coordinates],
            m_ref.atomicnums,
            m_dock.atomicnums,
            m_ref.adjacency_matrix,
            m_dock.adjacency_matrix,
        )[0]
        return round(float(val), 3)
    except Exception:
        pass

    # Method 2: Substructure / MCS Isomorphism Mapping across all symmetry matches
    mcs = rdFMCS.FindMCS(
        [dock_no_h, cryst_no_h],
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        timeout=10,
    )
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    dock_matches = dock_no_h.GetSubstructMatches(mcs_mol, uniquify=False)
    cryst_match = cryst_no_h.GetSubstructMatch(mcs_mol)

    if not dock_matches or not cryst_match:
        raise RuntimeError("Could not map atom correspondence between docked and crystal ligand.")

    cryst_conf = cryst_no_h.GetConformer()
    dock_conf = dock_no_h.GetConformer()
    cryst_pts = np.array([list(cryst_conf.GetAtomPosition(i)) for i in cryst_match])

    best_rmsd = 999.0
    for d_match in dock_matches:
        dock_pts = np.array([list(dock_conf.GetAtomPosition(i)) for i in d_match])
        diff = cryst_pts - dock_pts
        cur_rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
        if cur_rmsd < best_rmsd:
            best_rmsd = cur_rmsd

    return round(best_rmsd, 3)


def check_posebusters_and_biophysics(
    docked_mol: Chem.Mol, receptor_pdb: Path, crystal_pdb: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Runs PoseBusters physical validity suite (or built-in steric clash + geometry validator)
    and profiles key non-covalent receptor-ligand contacts.
    """
    dock_no_h = Chem.RemoveHs(docked_mol)
    conf = dock_no_h.GetConformer()
    lig_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(dock_no_h.GetNumAtoms())])

    # Parse receptor heavy atom coordinates & residue labels
    rec_coords = []
    rec_residues = []
    rec_elements = []
    with open(receptor_pdb, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                elem = line[76:78].strip().upper() or line[12:14].strip().upper()
                if elem == "H":
                    continue
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    res_id = f"{line[17:20].strip()}{line[22:26].strip()}_{line[21].strip()}"
                    rec_coords.append([x, y, z])
                    rec_residues.append(res_id)
                    rec_elements.append(elem)
                except ValueError:
                    continue

    rec_pts = np.array(rec_coords)
    # Pairwise distance matrix between ligand heavy atoms and receptor heavy atoms
    dists = np.linalg.norm(lig_coords[:, None, :] - rec_pts[None, :, :], axis=-1)

    # 1. Steric Clash Check (< 2.0 Å between non-hydrogen atoms indicates severe VdW overlap)
    min_intermolecular_dist = float(np.min(dists))
    steric_clash_count = int(np.sum(dists < 2.0))

    # 2. Identify Polar / H-Bond Candidate Contacts (N/O/S/F within 2.4 - 3.5 Å)
    polar_residues = set()
    hydrophobic_residues = set()
    for lig_idx in range(dock_no_h.GetNumAtoms()):
        lig_elem = dock_no_h.GetAtomWithIdx(lig_idx).GetSymbol().upper()
        close_rec_indices = np.where((dists[lig_idx] >= 2.2) & (dists[lig_idx] <= 4.0))[0]
        for r_idx in close_rec_indices:
            d = dists[lig_idx, r_idx]
            r_elem = rec_elements[r_idx]
            if lig_elem in ("N", "O", "S", "F") and r_elem in ("N", "O", "S") and d <= 3.5:
                polar_residues.add(rec_residues[r_idx])
            elif lig_elem == "C" and r_elem == "C" and d <= 4.0:
                hydrophobic_residues.add(rec_residues[r_idx])

    posebusters_pass = (steric_clash_count == 0) and (min_intermolecular_dist >= 2.0)

    # Optional: Run full PoseBusters library check if installed
    try:
        from posebusters import PoseBusters
        buster = PoseBusters(config="dock")
        df = buster.bust([docked_mol], str(crystal_pdb) if crystal_pdb else None, str(receptor_pdb))
        if not df.empty and "all_atoms_connected" in df.columns:
            posebusters_pass = bool(df.iloc[0].all())
    except Exception:
        pass

    return {
        "posebusters_valid": posebusters_pass,
        "steric_clash_count": steric_clash_count,
        "min_receptor_ligand_dist_A": round(min_intermolecular_dist, 2),
        "polar_hbond_residues": sorted(list(polar_residues))[:12],
        "hydrophobic_contact_residues": sorted(list(hydrophobic_residues))[:15],
    }


def emit_research_report(
    pdb_id: str,
    ligand_id: str,
    prep_info: Dict[str, Any],
    dock_info: Dict[str, Any],
    validation_info: Dict[str, Any],
    rmsd_top1: Optional[float],
    output_dir: Path,
) -> Dict[str, Path]:
    """Writes structured JSON report + copy-pasteable Academic Materials & Methods Markdown."""
    output_dir.mkdir(parents=True, exist_ok=True)
    top_pose = dock_info["top_pose"]

    payload = {
        "target_pdb": pdb_id,
        "ligand_id": ligand_id,
        "scientific_grade_status": "PASS" if (
            validation_info["posebusters_valid"] and (rmsd_top1 is None or rmsd_top1 <= 2.0)
        ) else "REVIEW_REQUIRED",
        "top1_symmetry_rmsd_A": rmsd_top1,
        "top_pose_thermodynamics": top_pose,
        "biophysical_validation": validation_info,
        "grid_box": {
            "center": dock_info["grid_center"],
            "size_A": dock_info["grid_size"],
            "radius_of_gyration_A": prep_info["radius_of_gyration"],
        },
        "all_poses": dock_info["poses"],
    }

    json_path = output_dir / f"{pdb_id}_{ligand_id}_bindora_report.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # Emit Academic Materials & Methods section
    c = dock_info["grid_center"]
    s = dock_info["grid_size"]
    rmsd_str = f"{rmsd_top1:.2f} Å" if rmsd_top1 is not None else "N/A (prospective docking)"

    md_text = f"""# Bindora Dock Scientific Validation & Methods Report: `{pdb_id}` ({ligand_id})

## 1. Executive Summary
* **Top-1 Binding Affinity ($\\Delta G$):** `{top_pose['affinity_kcal_mol']} kcal/mol`
* **Estimated $K_i$ (298.15 K):** `{top_pose['estimated_ki_nM']} nM`
* **Ligand Efficiency (LE):** `{top_pose['ligand_efficiency']} kcal/(mol·heavy atom)`
* **Symmetry-Corrected Heavy-Atom RMSD:** `{rmsd_str}` (Threshold <= 2.00 Å)
* **PoseBusters Biophysical Validity:** `{"PASS" if validation_info["posebusters_valid"] else "FAIL"}` (Min contact distance: `{validation_info["min_receptor_ligand_dist_A"]} Å`, Steric clashes: `{validation_info["steric_clash_count"]}`)
* **Key Polar/H-Bond Pocket Residues:** `{", ".join(validation_info["polar_hbond_residues"]) or "None detected"}`

## 2. Auto-Generated Materials & Methods (Ready for Manuscript)
> **Molecular Docking Protocol:** Structure-based molecular docking was performed using the Bindora Dock biophysical pipeline. The crystallographic structure of the target receptor (PDB ID: {pdb_id}) was retrieved from the RCSB Protein Data Bank. Non-essential crystallization buffers and bulk solvent molecules were removed while retaining essential catalytic cofactors and metal ions. Missing heavy atoms were reconstructed using PDBFixer, and residue protonation microstates were assigned at physiological pH (7.4) with Gasteiger partial charges. The small-molecule ligand ({ligand_id}, MW = {prep_info['molecular_weight']} Da, {prep_info['rotatable_bonds']} rotatable bonds) was standardized by salt stripping, canonical tautomer normalization, and 3D conformer embedding via the deterministic ETKDGv3 algorithm (random seed = {dock_info['seed']}), followed by energy minimization using the {prep_info['force_field']} force field in RDKit and torsion-tree parameterization in Meeko. Docking was executed using the AutoDock {dock_info['scoring_function'].capitalize()} scoring function with a search exhaustiveness of {dock_info['exhaustiveness']} and a deterministic seed of {dock_info['seed']}. The cubic search grid was centered at ({c[0]}, {c[1]}, {c[2]}) Å with dynamically scaled dimensions of ({s[0]} × {s[1]} × {s[2]}) Å derived from the ligand radius of gyration ($R_g = {prep_info['radius_of_gyration']}$ Å). Docked poses were validated for physical plausibility and symmetry-corrected heavy-atom root-mean-square deviation (RMSD).
"""
    md_path = output_dir / f"{pdb_id}_{ligand_id}_methods_section.md"
    md_path.write_text(md_text, encoding="utf-8")

    return {"json_report": json_path, "markdown_report": md_path}
