# bindora/fetcher.py
import os
import requests
import numpy as np
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any

# Biophysical Whitelists & Blacklists
ESSENTIAL_METALS = {"ZN", "MG", "MN", "FE", "CA", "CU", "CO", "NI"}
ESSENTIAL_COFACTORS = {"HEM", "HEC", "FAD", "FMN", "NAD", "NAP", "NDP", "ATP", "ADP", "SAM", "SAH", "PLP"}
CRYSTAL_ARTIFACTS = {
    "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "DMS", "ACT", "ACE",
    "CL", "NA", "K", "BR", "IOD", "FMT", "TRS", "MES", "EPE", "CIT", "MPD"
}

class PDBFetcherAndSplitter:
    """Fetches PDB structures and cleanly separates receptor, cofactors, and reference ligands."""

    def __init__(self, work_dir: str = "./workspace"):
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def fetch_pdb(self, pdb_id: str) -> Path:
        pdb_id = pdb_id.strip().upper()
        out_path = self.work_dir / f"{pdb_id}_raw.pdb"
        if out_path.exists() and out_path.stat().st_size > 100:
            return out_path

        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            raise RuntimeError(f"Failed to fetch PDB {pdb_id} from RCSB (HTTP {resp.status_code})")
        
        out_path.write_text(resp.text, encoding="utf-8")
        return out_path

    def split_receptor_and_ligand(
        self,
        pdb_path: Path,
        ligand_resname: str,
        target_chain: Optional[str] = None,
        retain_metals: bool = True,
        retain_cofactors: bool = True,
        retain_bridging_waters: bool = False,
    ) -> Dict[str, Any]:
        """
        Splits raw PDB into sanitized receptor PDB and reference crystal ligand PDB,
        while calculating exact binding site center and dimensions.
        """
        ligand_resname = ligand_resname.strip().upper()
        receptor_lines: List[str] = []
        ligand_lines: List[str] = []
        water_lines: List[Tuple[str, np.ndarray]] = []

        receptor_coords: List[np.ndarray] = []
        ligand_coords: List[np.ndarray] = []
        selected_lig_chain: Optional[str] = None if target_chain == "ALL" else target_chain
        selected_lig_resseq: Optional[str] = None

        with open(pdb_path, "r", encoding="utf-8") as f:
            for line in f:
                record = line[:6].strip()
                if record not in ("ATOM", "HETATM"):
                    continue

                # Handle Alternate Locations (keep only ' ' or 'A')
                alt_loc = line[16]
                if alt_loc not in (" ", "A", "1"):
                    continue

                resname = line[17:20].strip().upper()
                chain_id = line[21].strip()
                resseq = line[22:26].strip()

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord = np.array([x, y, z], dtype=float)
                except ValueError:
                    continue

                # Standard Protein Residues
                if record == "ATOM":
                    if target_chain is None or chain_id == target_chain or target_chain == "ALL":
                        clean_line = line[:16] + " " + line[17:]
                        receptor_lines.append(clean_line)
                        receptor_coords.append(coord)
                    continue

                # HETATM handling
                if resname == ligand_resname:
                    # Lock onto the first matching ligand instance if multiple copies exist in asymmetric unit
                    if selected_lig_chain is None:
                        selected_lig_chain = chain_id
                    if selected_lig_resseq is None and chain_id == selected_lig_chain:
                        selected_lig_resseq = resseq

                    if chain_id == selected_lig_chain and resseq == selected_lig_resseq:
                        clean_line = line[:16] + " " + line[17:]
                        ligand_lines.append(clean_line)
                        ligand_coords.append(coord)
                    continue

                if resname in ("HOH", "WAT"):
                    water_lines.append((line[:16] + " " + line[17:], coord))
                    continue

                if retain_metals and resname in ESSENTIAL_METALS:
                    receptor_lines.append(line[:16] + " " + line[17:])
                    receptor_coords.append(coord)
                    continue

                if retain_cofactors and resname in ESSENTIAL_COFACTORS:
                    receptor_lines.append(line[:16] + " " + line[17:])
                    receptor_coords.append(coord)
                    continue

                # Ignore crystallization artifacts (SO4, GOL, EDO, etc.)
                if resname in CRYSTAL_ARTIFACTS:
                    continue

        if not ligand_coords:
            raise ValueError(
                f"Co-crystallized ligand '{ligand_resname}' not found in {pdb_path.name}."
            )

        lig_pts = np.vstack(ligand_coords)
        rec_pts = np.vstack(receptor_coords)

        # Optional: Retain conserved structural waters bridging receptor & ligand (< 3.2 Å to both)
        bridging_waters_kept = 0
        if retain_bridging_waters and water_lines:
            for w_line, w_coord in water_lines:
                dist_to_lig = np.min(np.linalg.norm(lig_pts - w_coord, axis=1))
                dist_to_rec = np.min(np.linalg.norm(rec_pts - w_coord, axis=1))
                if dist_to_lig <= 3.2 and dist_to_rec <= 3.2:
                    receptor_lines.append(w_line)
                    bridging_waters_kept += 1

        stem = pdb_path.stem.replace("_raw", "")
        rec_out = self.work_dir / f"{stem}_receptor_clean.pdb"
        lig_out = self.work_dir / f"{stem}_{ligand_resname}_crystal.pdb"

        rec_out.write_text("".join(receptor_lines) + "END\n", encoding="utf-8")
        lig_out.write_text("".join(ligand_lines) + "END\n", encoding="utf-8")

        pocket_center = np.mean(lig_pts, axis=0)
        pocket_span = np.ptp(lig_pts, axis=0)  # max - min along x, y, z

        return {
            "receptor_pdb": rec_out,
            "crystal_ligand_pdb": lig_out,
            "pocket_center": tuple(float(x) for x in np.round(pocket_center, 3)),
            "ligand_span": tuple(float(x) for x in np.round(pocket_span, 3)),
            "bridging_waters_kept": bridging_waters_kept,
        }

    def extract_apo_or_residue_pocket(
        self,
        pdb_path: Path,
        active_site_residues: Optional[List[str]] = None,
        target_chain: Optional[str] = None,
        retain_metals: bool = True,
        retain_cofactors: bool = True,
    ) -> Dict[str, Any]:
        """
        Prepares receptor and computes binding pocket center when no co-crystallized
        ligand is present (using specified active-site residues or geometric interior centroid).
        active_site_residues format: ['A:25', 'B:25'] or ['25', '27', '48']
        """
        receptor_lines: List[str] = []
        receptor_coords: List[np.ndarray] = []
        active_coords: List[np.ndarray] = []

        target_res_set = set()
        if active_site_residues:
            target_res_set = {r.strip().upper() for r in active_site_residues}

        with open(pdb_path, "r", encoding="utf-8") as f:
            for line in f:
                record = line[:6].strip()
                if record not in ("ATOM", "HETATM"):
                    continue
                alt_loc = line[16]
                if alt_loc not in (" ", "A", "1"):
                    continue

                resname = line[17:20].strip().upper()
                chain_id = line[21].strip().upper()
                resseq = line[22:26].strip()

                try:
                    coord = np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                        dtype=float,
                    )
                except ValueError:
                    continue

                if record == "ATOM":
                    if target_chain in (None, "ALL", chain_id):
                        clean_line = line[:16] + " " + line[17:]
                        receptor_lines.append(clean_line)
                        receptor_coords.append(coord)
                        if (
                            f"{chain_id}:{resseq}" in target_res_set
                            or resseq in target_res_set
                            or f"{resname}{resseq}" in target_res_set
                        ):
                            active_coords.append(coord)
                elif record == "HETATM":
                    if (retain_metals and resname in ESSENTIAL_METALS) or (
                        retain_cofactors and resname in ESSENTIAL_COFACTORS
                    ):
                        receptor_lines.append(line[:16] + " " + line[17:])
                        receptor_coords.append(coord)

        if not receptor_coords:
            raise ValueError(f"No valid receptor ATOM records found in {pdb_path.name}")

        rec_pts = np.vstack(receptor_coords)
        pocket_pts = np.vstack(active_coords) if active_coords else rec_pts

        stem = pdb_path.stem.replace("_raw", "")
        rec_out = self.work_dir / f"{stem}_receptor_clean.pdb"
        rec_out.write_text("".join(receptor_lines) + "END\n", encoding="utf-8")

        mean_center = np.mean(pocket_pts, axis=0)
        return {
            "receptor_pdb": rec_out,
            "crystal_ligand_pdb": None,
            "pocket_center": tuple(float(x) for x in np.round(mean_center, 3)),
            "pocket_mode": "residue_guided" if active_coords else "geometric_centroid",
        }

