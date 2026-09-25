import os
import re
import shutil
import tempfile
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional
import numpy as np
from scipy.spatial import Voronoi, cKDTree
from backend.config import get_subprocess_kwargs

HYDROPHOBIC_RESIDUES = {
    "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"
}

class PocketDetectionService:
    """Service to perform blind pocket detection using fpocket or Voronoi alpha-sphere cavity clustering."""

    @staticmethod
    def detect_pockets(pdb_content: str, max_pockets: int = 5) -> List[Dict[str, Any]]:
        """Detect candidate binding pockets on a receptor structure.
        
        Attempts to run fpocket binary if available on system; otherwise seamlessly
        executes pure-Python Voronoi alpha-sphere tessellation (Le Guilloux et al. 2009).
        """
        if not pdb_content or not pdb_content.strip():
            return []

        # 1. Attempt fpocket executable if installed
        fpocket_path = shutil.which("fpocket")
        if fpocket_path:
            try:
                pockets = PocketDetectionService._run_fpocket(fpocket_path, pdb_content, max_pockets)
                if pockets:
                    return pockets
            except Exception as e:
                print(f"[POCKET DETECTION] fpocket execution failed, falling back to Voronoi: {e}")

        # 2. Check for fpocket in WSL if on Windows
        if os.name == "nt" and shutil.which("wsl"):
            try:
                wsl_check = subprocess.run(["wsl", "which", "fpocket"], capture_output=True, text=True, timeout=5, **get_subprocess_kwargs())
                if wsl_check.returncode == 0 and wsl_check.stdout.strip():
                    pockets = PocketDetectionService._run_fpocket_wsl(pdb_content, max_pockets)
                    if pockets:
                        return pockets
            except Exception:
                pass

        # 3. Pure-Python Voronoi alpha-sphere cavity clustering (native fallback)
        return PocketDetectionService._run_voronoi_alpha_spheres(pdb_content, max_pockets)

    @staticmethod
    def _run_fpocket(fpocket_bin: str, pdb_content: str, max_pockets: int) -> List[Dict[str, Any]]:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            pdb_file = tmp_path / "receptor.pdb"
            pdb_file.write_text(pdb_content, encoding="utf-8")

            cmd = [fpocket_bin, "-f", str(pdb_file)]
            res = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True, timeout=60, **get_subprocess_kwargs())
            if res.returncode != 0:
                return []

            out_dir = tmp_path / "receptor_out"
            info_file = out_dir / "receptor_info.txt"
            if not info_file.exists():
                return []

            return PocketDetectionService._parse_fpocket_info(info_file.read_text(encoding="utf-8"), out_dir, max_pockets)

    @staticmethod
    def _run_fpocket_wsl(pdb_content: str, max_pockets: int) -> List[Dict[str, Any]]:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            pdb_file = tmp_path / "receptor.pdb"
            pdb_file.write_text(pdb_content, encoding="utf-8")

            # Convert Windows path to WSL path
            wsl_path = subprocess.run(["wsl", "wslpath", "-a", str(pdb_file).replace("\\", "/")], capture_output=True, text=True, **get_subprocess_kwargs()).stdout.strip()

            cmd = ["wsl", "fpocket", "-f", wsl_path]
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=60, **get_subprocess_kwargs())
            if res.returncode != 0:
                return []

            out_dir = tmp_path / "receptor_out"
            info_file = out_dir / "receptor_info.txt"
            if not info_file.exists():
                return []

            return PocketDetectionService._parse_fpocket_info(info_file.read_text(encoding="utf-8"), out_dir, max_pockets)

    @staticmethod
    def _parse_fpocket_info(info_text: str, out_dir: Path, max_pockets: int) -> List[Dict[str, Any]]:
        pockets = []
        blocks = info_text.split("Pocket ")
        for b in blocks[1:]:
            lines = b.splitlines()
            if not lines:
                continue
            try:
                pocket_id = int(lines[0].split(":")[0].strip())
            except Exception:
                pocket_id = len(pockets) + 1

            score = 0.5
            as_count = 0
            volume = 500.0
            for l in lines:
                if "Druggability Score" in l:
                    try:
                        score = float(l.split(":")[1].strip())
                    except Exception:
                        pass
                elif "Number of Alpha Spheres" in l:
                    try:
                        as_count = int(l.split(":")[1].strip())
                    except Exception:
                        pass
                elif "Volume" in l:
                    try:
                        volume = float(l.split(":")[1].strip())
                    except Exception:
                        pass

            # Read centroid from pocket PDB if present
            pocket_pdb = out_dir / f"pockets/pocket{pocket_id}_vert.pqr"
            cx, cy, cz = 0.0, 0.0, 0.0
            coords = []
            if pocket_pdb.exists():
                for pl in pocket_pdb.read_text(encoding="utf-8", errors="ignore").splitlines():
                    if pl.startswith(("ATOM", "HETATM")):
                        try:
                            coords.append([float(pl[30:38]), float(pl[38:46]), float(pl[46:54])])
                        except Exception:
                            pass
            if coords:
                ca = np.mean(coords, axis=0)
                cx, cy, cz = float(ca[0]), float(ca[1]), float(ca[2])
                cmin = np.min(coords, axis=0)
                cmax = np.max(coords, axis=0)
                sx = float(min(32.0, max(18.0, (cmax[0] - cmin[0]) + 8.0)))
                sy = float(min(32.0, max(18.0, (cmax[1] - cmin[1]) + 8.0)))
                sz = float(min(32.0, max(18.0, (cmax[2] - cmin[2]) + 8.0)))
            else:
                sx, sy, sz = 22.0, 22.0, 22.0

            pockets.append({
                "pocket_id": pocket_id,
                "rank": len(pockets) + 1,
                "center": {"x": round(cx, 2), "y": round(cy, 2), "z": round(cz, 2)},
                "size": {"x": round(sx, 1), "y": round(sy, 1), "z": round(sz, 1)},
                "druggability_score": round(score, 3),
                "alpha_spheres": as_count,
                "volume_a3": round(volume, 1),
                "method": "fpocket",
                "description": f"Blind pocket #{len(pockets) + 1} (Druggability: {score:.2f}, {as_count} alpha-spheres)"
            })
            if len(pockets) >= max_pockets:
                break
        return pockets

    @staticmethod
    def _run_voronoi_alpha_spheres(pdb_content: str, max_pockets: int) -> List[Dict[str, Any]]:
        """Scientific Voronoi-tessellation / alpha-sphere cavity clustering engine."""
        atom_res = []
        coords = []
        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")) and not line[17:20].strip() in ("HOH", "WAT", "DOD", "TIP"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                    atom_res.append(line[17:20].strip())
                except Exception:
                    pass

        if len(coords) < 30:
            return []

        coords = np.array(coords)
        atom_tree = cKDTree(coords)

        # 1. Voronoi Tessellation
        try:
            vor = Voronoi(coords)
        except Exception as ve:
            print(f"[POCKET VORONOI ERROR] {ve}")
            return []

        # 2. Distance query to find alpha-spheres
        dists, _ = atom_tree.query(vor.vertices)

        # Alpha-sphere radii between 2.8Å and 5.0Å correspond to pocket clefts
        mask = (dists >= 2.8) & (dists <= 5.0)
        alpha_pts = vor.vertices[mask]
        if len(alpha_pts) < 15:
            return []

        # 3. Filter for buriedness (surrounded by at least 12 protein atoms within 6Å)
        num_neighbors = atom_tree.query_ball_point(alpha_pts, r=6.0)
        buried_mask = np.array([len(n) >= 12 for n in num_neighbors])
        alpha_pts = alpha_pts[buried_mask]
        if len(alpha_pts) < 15:
            return []

        # 4. Spatial clustering of alpha spheres within 2.8Å
        sphere_tree = cKDTree(alpha_pts)
        pairs = sphere_tree.query_pairs(r=2.8)
        adj = {i: [] for i in range(len(alpha_pts))}
        for i, j in pairs:
            adj[i].append(j)
            adj[j].append(i)

        visited = set()
        clusters = []
        for i in range(len(alpha_pts)):
            if i not in visited:
                comp = []
                queue = [i]
                visited.add(i)
                while queue:
                    node = queue.pop(0)
                    comp.append(node)
                    for neighbor in adj[node]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
                if len(comp) >= 15:
                    clusters.append(comp)

        if not clusters:
            return []

        # 5. Score and rank pockets
        scored_pockets = []
        for idx, comp in enumerate(clusters):
            pts = alpha_pts[comp]
            cx, cy, cz = pts.mean(axis=0)
            cmin = pts.min(axis=0)
            cmax = pts.max(axis=0)
            sx = float(min(32.0, max(20.0, (cmax[0] - cmin[0]) + 8.0)))
            sy = float(min(32.0, max(20.0, (cmax[1] - cmin[1]) + 8.0)))
            sz = float(min(32.0, max(20.0, (cmax[2] - cmin[2]) + 8.0)))

            # Approximate volume
            vol = float(len(comp) * 11.5)

            # Hydrophobicity of neighboring residues
            contact_indices = atom_tree.query_ball_point([cx, cy, cz], r=7.0)
            hydrophobic_count = sum(1 for ai in contact_indices if atom_res[ai] in HYDROPHOBIC_RESIDUES)
            hydro_ratio = hydrophobic_count / max(1, len(contact_indices))

            # Druggability score (size + hydrophobicity)
            size_norm = min(1.0, len(comp) / 60.0)
            drug_score = min(0.98, max(0.15, 0.45 * size_norm + 0.55 * hydro_ratio))

            scored_pockets.append({
                "pocket_id": idx + 1,
                "center": {"x": round(float(cx), 2), "y": round(float(cy), 2), "z": round(float(cz), 2)},
                "size": {"x": round(sx, 1), "y": round(sy, 1), "z": round(sz, 1)},
                "druggability_score": round(float(drug_score), 3),
                "alpha_spheres": len(comp),
                "volume_a3": round(vol, 1),
                "method": "voronoi_alpha_sphere"
            })

        # Sort descending by druggability score
        scored_pockets.sort(key=lambda p: (p["druggability_score"], p["alpha_spheres"]), reverse=True)

        # Assign final ranks
        for rank, p in enumerate(scored_pockets[:max_pockets], 1):
            p["rank"] = rank
            p["description"] = f"Blind pocket #{rank} (Druggability: {p['druggability_score']:.2f}, {p['alpha_spheres']} alpha-spheres)"

        return scored_pockets[:max_pockets]
