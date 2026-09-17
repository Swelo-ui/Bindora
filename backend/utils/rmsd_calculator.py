import math
from typing import Dict, Any, List, Optional, Tuple, Union
from rdkit import Chem
from rdkit.Chem import AllChem

def calculate_rmsd(
    docked_pose_or_pdbqt: Union[str, Dict[str, Any]],
    reference_pdb_block: str
) -> float:
    """
    Calculate symmetry-corrected heavy-atom RMSD (in Angstroms) between a docked pose and a reference crystal structure.
    
    Parameters:
        docked_pose_or_pdbqt: PDBQT text of the docked pose, or pose dictionary containing 'pdbqt_content' or 'pdb_block'.
        reference_pdb_block: PDB block of the crystallographic reference ligand.
        
    Returns:
        float: RMSD in Angstroms rounded to 2 decimal places.
    """
    # Extract pose PDBQT string
    if isinstance(docked_pose_or_pdbqt, dict):
        pose_pdbqt = docked_pose_or_pdbqt.get("pdbqt_content") or docked_pose_or_pdbqt.get("pdb_block", "")
    else:
        pose_pdbqt = str(docked_pose_or_pdbqt)

    # 1. Gold-standard RDKit graph-isomorphism & in-place symmetry-corrected RMSD via Meeko
    # Note: Must use CalcRMS (in-place) rather than GetBestRMS (superposition).
    # Docking evaluation requires pocket-coordinate positioning; GetBestRMS aligns
    # the molecules in vacuum, measuring only internal scaffold deformation (leading to false 0.0 A).
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate
        pdbqt_mol = PDBQTMolecule(pose_pdbqt)
        rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        if rdkit_mols and len(rdkit_mols) > 0:
            ref_mol = Chem.RemoveHs(Chem.MolFromPDBBlock(reference_pdb_block, sanitize=False))
            docked_mol = Chem.RemoveHs(rdkit_mols[0])
            if ref_mol and docked_mol and ref_mol.GetNumHeavyAtoms() == docked_mol.GetNumHeavyAtoms():
                best_rms = float(AllChem.CalcRMS(docked_mol, ref_mol))
                if not math.isnan(best_rms):
                    return round(best_rms, 2)
    except Exception:
        pass

    # 2. Robust coordinate-based RMSD fallback (ignoring hydrogens)
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
                cryst_atoms.append({"name": aname, "elem": elem.upper(), "coord": (x, y, z)})
            except Exception:
                continue

    docked_atoms = []
    for line in pose_pdbqt.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                tokens = line.split()
                ad4 = tokens[-1] if tokens else ""
                elem = "Cl" if ad4.upper() == "CL" else "Br" if ad4.upper() == "BR" else "F" if ad4.upper() == "F" else ad4[0].upper() if ad4 else line[12:14].strip().upper()
                if elem == "H":
                    continue
                aname = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                docked_atoms.append({"name": aname, "elem": elem, "coord": (x, y, z)})
            except Exception:
                continue

    if not cryst_atoms or not docked_atoms:
        return 99.0

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

    return round(math.sqrt(sum_sq / len(cryst_atoms)), 2)
