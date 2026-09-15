import re
import math
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski
from meeko import MoleculePreparation, PDBQTWriterLegacy
from backend.config import VINA_EXE
from backend.utils.vina_setup import ensure_vina

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}

SOLVENTS_AND_IONS = {
    "HOH", "WAT", "DOD", "TIP", "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE",
    "SO4", "PO4", "GOL", "EDO", "DMS", "ACT", "FMT"
}

class DockingEngine:
    """Service to handle receptor & ligand preparation, Vina docking execution, and intermolecular contact analysis."""

    @staticmethod
    def prepare_receptor(pdb_content: str, target_chain: Optional[str] = None) -> Dict[str, Any]:
        """Clean receptor PDB, identify binding site / co-crystallized ligand, and produce PDBQT."""
        lines = pdb_content.splitlines()
        protein_lines = []
        hetatm_ligand_atoms = []
        all_ca_coords = []
        chains_found = set()

        for line in lines:
            if line.startswith("ATOM  "):
                res_name = line[17:20].strip()
                chain = line[21:22].strip()
                chains_found.add(chain)
                if target_chain and chain != target_chain:
                    continue
                if res_name in STANDARD_AMINO_ACIDS:
                    protein_lines.append(line)
                    atom_name = line[12:16].strip()
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        if atom_name == "CA":
                            all_ca_coords.append((x, y, z))
                    except Exception:
                        pass
            elif line.startswith("HETATM"):
                res_name = line[17:20].strip()
                chain = line[21:22].strip()
                if res_name not in SOLVENTS_AND_IONS:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        hetatm_ligand_atoms.append({
                            "name": res_name,
                            "chain": chain,
                            "coord": (x, y, z)
                        })
                    except Exception:
                        pass

        if not protein_lines:
            raise ValueError("No standard amino acid protein atoms found in receptor PDB file.")

        # Determine grid box center
        has_co_ligand = len(hetatm_ligand_atoms) > 5
        co_ligand_name = hetatm_ligand_atoms[0]["name"] if has_co_ligand else None

        if has_co_ligand:
            coords = [a["coord"] for a in hetatm_ligand_atoms]
            center_x = sum(c[0] for c in coords) / len(coords)
            center_y = sum(c[1] for c in coords) / len(coords)
            center_z = sum(c[2] for c in coords) / len(coords)
            size_x, size_y, size_z = 22.0, 22.0, 22.0
            pocket_desc = f"Auto-centered on co-crystallized ligand pocket ({co_ligand_name})"
        elif all_ca_coords:
            center_x = sum(c[0] for c in all_ca_coords) / len(all_ca_coords)
            center_y = sum(c[1] for c in all_ca_coords) / len(all_ca_coords)
            center_z = sum(c[2] for c in all_ca_coords) / len(all_ca_coords)
            size_x, size_y, size_z = 26.0, 26.0, 26.0
            pocket_desc = "Auto-centered on protein geometric center"
        else:
            center_x, center_y, center_z = 0.0, 0.0, 0.0
            size_x, size_y, size_z = 24.0, 24.0, 24.0
            pocket_desc = "Default center"

        # Generate PDBQT format lines for protein
        pdbqt_lines = []
        for line in protein_lines:
            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            elem = line[76:78].strip() or atom_name[0]
            ad4_type = elem
            if elem == "C":
                ad4_type = "A" if res_name in ("PHE", "TYR", "TRP", "HIS") else "C"
            elif elem == "O":
                ad4_type = "OA"
            elif elem == "N":
                ad4_type = "NA" if res_name in ("HIS", "TRP") else "N"
            elif elem == "S":
                ad4_type = "SA"
            elif elem == "H":
                ad4_type = "HD"

            charge = 0.00
            pdbqt_line = f"{line[:54]:<54}{0.00:>6.2f}{0.00:>6.2f}    {charge:>6.3f} {ad4_type:<2}"
            pdbqt_lines.append(pdbqt_line)

        cleaned_pdb = "\n".join(protein_lines) + "\nEND\n"
        pdbqt_text = "\n".join(pdbqt_lines) + "\nTER\nEND\n"

        return {
            "cleaned_pdb": cleaned_pdb,
            "pdbqt_text": pdbqt_text,
            "atom_count": len(protein_lines),
            "chains": sorted(list(chains_found)),
            "detected_pocket": {
                "has_co_crystallized_ligand": has_co_ligand,
                "co_ligand_name": co_ligand_name,
                "description": pocket_desc,
                "center": {"x": round(center_x, 2), "y": round(center_y, 2), "z": round(center_z, 2)},
                "size": {"x": round(size_x, 1), "y": round(size_y, 1), "z": round(size_z, 1)}
            }
        }

    @staticmethod
    def prepare_ligand(input_data: str, is_sdf: bool = False) -> Dict[str, Any]:
        """Convert ligand SMILES or SDF into 3D conformer, energy-minimize, and produce PDBQT via Meeko."""
        if is_sdf or "\n" in input_data:
            mol = Chem.MolFromMolBlock(input_data)
            if not mol:
                mol = Chem.MolFromMolBlock(input_data, sanitize=False)
                if mol:
                    Chem.SanitizeMol(mol)
        else:
            smiles = input_data.strip()
            mol = Chem.MolFromSmiles(smiles)

        if not mol:
            raise ValueError("Failed to parse ligand structure from provided input.")

        # Ensure hydrogens are present
        mol_h = Chem.AddHs(mol)

        # Generate 3D coordinates using ETKDGv3
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        embed_result = AllChem.EmbedMolecule(mol_h, params)
        if embed_result != 0:
            # Fallback with random coordinates
            AllChem.EmbedMolecule(mol_h, useRandomCoords=True)

        # Energy minimization with MMFF94
        try:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
        except Exception:
            pass

        # Prepare PDBQT using Meeko
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol_h)
        if not mol_setups:
            raise RuntimeError("Meeko could not construct ligand flexible torsions setup.")

        pdbqt_str, is_ok, err_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            raise RuntimeError(f"Meeko PDBQT conversion failed: {err_msg}")

        # Also prepare PDB block for 3Dmol.js viewer
        pdb_block = Chem.MolToPDBBlock(mol_h)
        canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol_h))

        return {
            "pdbqt_text": pdbqt_str,
            "pdb_block": pdb_block,
            "canonical_smiles": canonical_smiles,
            "heavy_atom_count": mol.GetNumHeavyAtoms(),
            "rotatable_bonds": Lipinski.NumRotatableBonds(mol)
        }

    @staticmethod
    def run_docking(
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        center: Dict[str, float],
        size: Dict[str, float],
        exhaustiveness: int = 8,
        num_modes: int = 9
    ) -> List[Dict[str, Any]]:
        """Run AutoDock Vina on the prepared receptor and ligand PDBQT files."""
        vina_path = ensure_vina()

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            rec_file = tmp_path / "receptor.pdbqt"
            lig_file = tmp_path / "ligand.pdbqt"
            out_file = tmp_path / "docked_out.pdbqt"

            rec_file.write_text(receptor_pdbqt, encoding="utf-8")
            lig_file.write_text(ligand_pdbqt, encoding="utf-8")

            cmd = [
                str(vina_path),
                "--receptor", str(rec_file),
                "--ligand", str(lig_file),
                "--center_x", str(center["x"]),
                "--center_y", str(center["y"]),
                "--center_z", str(center["z"]),
                "--size_x", str(size["x"]),
                "--size_y", str(size["y"]),
                "--size_z", str(size["z"]),
                "--exhaustiveness", str(exhaustiveness),
                "--num_modes", str(num_modes),
                "--out", str(out_file)
            ]

            process = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            if process.returncode != 0:
                err_msg = process.stderr or process.stdout
                raise RuntimeError(f"AutoDock Vina execution error: {err_msg}")

            if not out_file.exists():
                raise RuntimeError("AutoDock Vina finished without generating an output PDBQT file.")

            output_text = out_file.read_text(encoding="utf-8")

        # Parse poses from output PDBQT
        poses = []
        current_mode = None
        current_affinity = None
        current_rmsd_lb = 0.0
        current_rmsd_ub = 0.0
        current_lines = []

        for line in output_text.splitlines():
            if line.startswith("MODEL"):
                parts = line.split()
                current_mode = int(parts[1]) if len(parts) > 1 else len(poses) + 1
                current_lines = [line]
            elif "REMARK VINA RESULT:" in line:
                # e.g.: REMARK VINA RESULT:    -8.4      0.000      0.000
                m = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                if len(m) >= 1:
                    current_affinity = float(m[0])
                    current_rmsd_lb = float(m[1]) if len(m) > 1 else 0.0
                    current_rmsd_ub = float(m[2]) if len(m) > 2 else 0.0
                current_lines.append(line)
            elif line.startswith("ENDMDL"):
                current_lines.append(line)
                pose_pdbqt = "\n".join(current_lines)
                # Convert PDBQT lines to standard PDB format for 3Dmol.js
                pdb_lines = []
                for pline in current_lines:
                    if pline.startswith(("ATOM", "HETATM")):
                        # Standard PDB format line
                        pdb_lines.append(pline[:66])
                pdb_block = "\n".join(pdb_lines) + "\nEND\n"

                poses.append({
                    "mode": current_mode or len(poses) + 1,
                    "affinity_kcal": current_affinity or 0.0,
                    "rmsd_lb": current_rmsd_lb,
                    "rmsd_ub": current_rmsd_ub,
                    "pdbqt_content": pose_pdbqt,
                    "pdb_block": pdb_block
                })
                current_lines = []
            else:
                if current_mode is not None:
                    current_lines.append(line)

        return poses

    @staticmethod
    def analyze_interactions(
        receptor_pdb: str,
        ligand_pdb_or_pdbqt: str,
        hbond_cutoff: float = 3.5,
        hydrophobic_cutoff: float = 4.0
    ) -> Dict[str, Any]:
        """Compute atomic contacts: hydrogen bonds and hydrophobic interactions between receptor and docked pose."""
        # 1. Parse receptor atoms
        receptor_atoms = []
        for line in receptor_pdb.splitlines():
            if line.startswith("ATOM  "):
                try:
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    chain = line[21:22].strip()
                    atom_name = line[12:16].strip()
                    elem = line[76:78].strip() or atom_name[0]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    receptor_atoms.append({
                        "res_name": res_name,
                        "res_num": res_num,
                        "chain": chain,
                        "atom_name": atom_name,
                        "elem": elem.upper(),
                        "coord": (x, y, z)
                    })
                except Exception:
                    continue

        # 2. Parse ligand atoms
        ligand_atoms = []
        for line in ligand_pdb_or_pdbqt.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    elem = line[76:78].strip() or atom_name[0]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ligand_atoms.append({
                        "atom_name": atom_name,
                        "elem": elem.upper(),
                        "coord": (x, y, z)
                    })
                except Exception:
                    continue

        hbonds = []
        hydrophobics = []
        contact_residues = set()

        hbond_donors_acceptors = {"N", "O"}

        for latom in ligand_atoms:
            lx, ly, lz = latom["coord"]
            lelem = latom["elem"]

            for ratom in receptor_atoms:
                rx, ry, rz = ratom["coord"]
                relem = ratom["elem"]

                # Quick bounding box filter
                if abs(lx - rx) > hydrophobic_cutoff or abs(ly - ry) > hydrophobic_cutoff or abs(lz - rz) > hydrophobic_cutoff:
                    continue

                dist = math.sqrt((lx - rx)**2 + (ly - ry)**2 + (lz - rz)**2)
                res_id = f"{ratom['res_name']} {ratom['res_num']}:{ratom['chain']}"

                # Hydrogen bond check
                if dist <= hbond_cutoff and lelem in hbond_donors_acceptors and relem in hbond_donors_acceptors:
                    hbonds.append({
                        "type": "Hydrogen Bond",
                        "distance": round(dist, 2),
                        "residue": res_id,
                        "res_name": ratom["res_name"],
                        "res_num": ratom["res_num"],
                        "chain": ratom["chain"],
                        "receptor_atom": ratom["atom_name"],
                        "ligand_atom": latom["atom_name"],
                        "start_coord": [lx, ly, lz],
                        "end_coord": [rx, ry, rz]
                    })
                    contact_residues.add(res_id)

                # Hydrophobic contact check (between carbon atoms)
                elif dist <= hydrophobic_cutoff and lelem == "C" and relem == "C":
                    hydrophobics.append({
                        "type": "Hydrophobic Contact",
                        "distance": round(dist, 2),
                        "residue": res_id,
                        "res_name": ratom["res_name"],
                        "res_num": ratom["res_num"],
                        "chain": ratom["chain"],
                        "receptor_atom": ratom["atom_name"],
                        "ligand_atom": latom["atom_name"],
                        "start_coord": [lx, ly, lz],
                        "end_coord": [rx, ry, rz]
                    })
                    contact_residues.add(res_id)

        # Deduplicate to top closest interactions per residue
        unique_hbonds = {}
        for hb in hbonds:
            key = f"{hb['residue']}_{hb['ligand_atom']}"
            if key not in unique_hbonds or unique_hbonds[key]["distance"] > hb["distance"]:
                unique_hbonds[key] = hb

        unique_hydrophobics = {}
        for hp in hydrophobics:
            key = f"{hp['residue']}_{hp['ligand_atom']}"
            if key not in unique_hydrophobics or unique_hydrophobics[key]["distance"] > hp["distance"]:
                unique_hydrophobics[key] = hp

        return {
            "hydrogen_bonds": list(unique_hbonds.values()),
            "hydrophobic_contacts": list(unique_hydrophobics.values())[:12],
            "interacting_residues": sorted(list(contact_residues)),
            "total_hbond_count": len(unique_hbonds),
            "total_hydrophobic_count": len(unique_hydrophobics)
        }
