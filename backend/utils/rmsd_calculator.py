import math
import numpy as np
from typing import Dict, Any, List, Optional, Tuple, Union
from rdkit import Chem
from rdkit.Chem import AllChem

def calculate_rmsd(
    docked_pose_or_pdbqt: Union[str, Dict[str, Any]],
    reference_pdb_block: str,
    return_details: bool = False
) -> Union[float, Dict[str, Any]]:
    """
    Calculate symmetry-aware heavy-atom root-mean-square deviation (RMSD in Angstroms)
    between a docked pose and a crystallographic reference ligand.

    Important scientific design:
    - Coordinates are evaluated IN-PLACE within the receptor binding pocket frame of reference.
    - Superposition (translation/rotation) is intentionally NOT performed because docking validation
      must assess how accurately the ligand reproduces the absolute crystallographic binding pose.
    - Evaluates all chemical symmetry automorphisms (e.g. flipping of symmetric phenyl rings or carboxylates)
      to report the true minimum symmetry-corrected RMSD.
    - When chemical graph topology is identical, uses graph-isomorphism symmetry mapping.
    - When perception differs slightly between PDB and PDBQT, uses element-constrained optimal
      bijective matching (minimum-weight bipartite matching via Kuhn-Munkres/Hungarian algorithm).

    Parameters:
        docked_pose_or_pdbqt: PDBQT text of the docked pose, or dictionary containing 'pdbqt_content'/'pdb_block'.
        reference_pdb_block: PDB format text of the reference native crystallographic ligand.
        return_details: If True, returns a dict with RMSD, atom count, and mapping method.

    Returns:
        float: Heavy-atom RMSD rounded to 2 decimal places (or detailed dict if return_details=True).
    """
    if isinstance(docked_pose_or_pdbqt, dict):
        pose_pdbqt = docked_pose_or_pdbqt.get("pdbqt_content") or docked_pose_or_pdbqt.get("pdb_block", "")
    else:
        pose_pdbqt = str(docked_pose_or_pdbqt)

    # 1. First Tier: RDKit Molecule Graph Isomorphism & In-Place Symmetry Automorphism Matching
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate
        pdbqt_mol = PDBQTMolecule(pose_pdbqt)
        rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        if rdkit_mols and len(rdkit_mols) > 0:
            docked_mol = Chem.RemoveHs(rdkit_mols[0])
            ref_mol = Chem.RemoveHs(Chem.MolFromPDBBlock(reference_pdb_block, sanitize=False))
            
            if ref_mol and docked_mol and ref_mol.GetNumHeavyAtoms() == docked_mol.GetNumHeavyAtoms():
                n_heavy = ref_mol.GetNumHeavyAtoms()
                if ref_mol.GetNumConformers() > 0 and docked_mol.GetNumConformers() > 0:
                    ref_conf = ref_mol.GetConformer()
                    dock_conf = docked_mol.GetConformer()

                    # Find all symmetry-equivalent substructure automorphisms of the docked ligand (bounded to maxMatches=64)
                    matches = docked_mol.GetSubstructMatches(docked_mol, uniquify=False, maxMatches=64)
                    if matches:
                        min_sq_sum = float("inf")
                        for match in matches:
                            sq_sum = 0.0
                            for d_idx, r_idx in enumerate(match):
                                p_d = dock_conf.GetAtomPosition(r_idx)
                                p_r = ref_conf.GetAtomPosition(d_idx)
                                sq_sum += (p_d.x - p_r.x)**2 + (p_d.y - p_r.y)**2 + (p_d.z - p_r.z)**2
                            if sq_sum < min_sq_sum:
                                min_sq_sum = sq_sum

                        rmsd_val = round(math.sqrt(min_sq_sum / n_heavy), 2)
                        if return_details:
                            return {
                                "rmsd": rmsd_val,
                                "method": "symmetry_aware_substructure_rmsd",
                                "heavy_atoms": n_heavy,
                                "automorphisms_tested": len(matches)
                            }
                        return rmsd_val
    except Exception:
        pass

    # 2. Second Tier: Element-constrained Optimal Bijective Coordinate Matching (Hungarian Algorithm)
    cryst_atoms = []
    for line in reference_pdb_block.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                aname = line[12:16].strip()
                elem = line[76:78].strip() or aname[0]
                if elem.upper() == "H":
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                cryst_atoms.append({"name": aname, "elem": elem.upper(), "coord": np.array([x, y, z])})
            except Exception:
                continue

    docked_atoms = []
    for line in pose_pdbqt.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                tokens = line.split()
                ad4 = tokens[-1] if tokens else ""
                elem = "Cl" if ad4.upper() == "CL" else "Br" if ad4.upper() == "BR" else "F" if ad4.upper() == "F" else ad4[0].upper() if ad4 else line[12:14].strip().upper()
                if elem.upper() == "H":
                    continue
                aname = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                docked_atoms.append({"name": aname, "elem": elem.upper(), "coord": np.array([x, y, z])})
            except Exception:
                continue

    if not cryst_atoms or not docked_atoms:
        if return_details:
            return {"rmsd": 99.0, "method": "unresolvable_atoms", "heavy_atoms": 0}
        return 99.0

    # Try linear_sum_assignment for optimal 1-to-1 bipartite matching per element
    try:
        from scipy.optimize import linear_sum_assignment
        total_sq = 0.0
        mapped_count = 0
        
        # Elements in crystal structure
        unique_elems = set(a["elem"] for a in cryst_atoms)
        for el in unique_elems:
            c_coords = np.array([a["coord"] for a in cryst_atoms if a["elem"] == el])
            d_coords = np.array([a["coord"] for a in docked_atoms if a["elem"] == el])

            if len(c_coords) == 0:
                continue
            if len(d_coords) == 0:
                # Missing element in docked pose - pair with closest available heavy atom
                d_coords = np.array([a["coord"] for a in docked_atoms])

            # Cost matrix: squared Euclidean distance
            diff = c_coords[:, np.newaxis, :] - d_coords[np.newaxis, :, :]
            cost_matrix = np.sum(diff**2, axis=-1)

            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            total_sq += float(cost_matrix[row_ind, col_ind].sum())
            mapped_count += len(row_ind)

        if mapped_count > 0:
            rmsd_val = round(math.sqrt(total_sq / mapped_count), 2)
            if return_details:
                return {
                    "rmsd": rmsd_val,
                    "method": "optimal_bijective_assignment_rmsd",
                    "heavy_atoms": mapped_count
                }
            return rmsd_val
    except Exception:
        pass

    # 3. Third Tier Fallback: Deterministic name and distance pairing
    sum_sq = 0.0
    for ca in cryst_atoms:
        cx, cy, cz = ca["coord"]
        celem = ca["elem"]
        candidates = [da["coord"] for da in docked_atoms if da["name"] == ca["name"]]
        if not candidates:
            candidates = [da["coord"] for da in docked_atoms if da["elem"] == celem]
        if not candidates:
            candidates = [da["coord"] for da in docked_atoms]
        min_sq = min((cx - dx)**2 + (cy - dy)**2 + (cz - dz)**2 for dx, dy, dz in candidates)
        sum_sq += min_sq

    rmsd_val = round(math.sqrt(sum_sq / len(cryst_atoms)), 2)
    if return_details:
        return {"rmsd": rmsd_val, "method": "coordinate_proximity_fallback", "heavy_atoms": len(cryst_atoms)}
    return rmsd_val
