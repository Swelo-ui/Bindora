# bindora/docking.py
import os
import re
import math
import shutil
import subprocess
from pathlib import Path
from typing import Tuple, List, Dict, Optional, Any

from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

R_GAS_KCAL = 1.98720425864083e-3  # kcal / (mol * K)
TEMP_KELVIN = 298.15              # 25 °C standard temperature


def calculate_thermodynamics(delta_g: float, heavy_atoms: int) -> Dict[str, float]:
    """Calculates dissociation constant Ki (in nanomolar) and Ligand Efficiency (LE)."""
    # Ki = exp(dG / (R * T)) in Molar, convert to nM (1e9)
    ki_molar = math.exp(delta_g / (R_GAS_KCAL * TEMP_KELVIN))
    ki_nm = round(ki_molar * 1e9, 2)
    le = round(-delta_g / max(1, heavy_atoms), 3)
    return {"affinity_kcal_mol": round(delta_g, 3), "estimated_ki_nM": ki_nm, "ligand_efficiency": le}


class BindoraDockingEngine:
    """Executes deterministic AutoDock Vina / Vinardo docking with full thermodynamic readout."""

    def __init__(
        self,
        work_dir: str = "./workspace",
        scoring_function: str = "vinardo",
        exhaustiveness: int = 8,
        num_modes: int = 9,
        seed: int = 42,
        cpu_threads: int = 0,
    ):
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.scoring = scoring_function.lower()
        if self.scoring not in ("vina", "vinardo", "ad4"):
            raise ValueError("scoring_function must be 'vina', 'vinardo', or 'ad4'")
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.seed = seed
        self.cpu_threads = cpu_threads

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        center: Tuple[float, float, float],
        box_size: Tuple[float, float, float],
        heavy_atoms: int,
    ) -> Dict[str, Any]:
        out_pdbqt = self.work_dir / f"{ligand_pdbqt.stem}_docked_{self.scoring}.pdbqt"
        out_sdf = self.work_dir / f"{ligand_pdbqt.stem}_docked_{self.scoring}.sdf"

        # Prefer Python bindings if available, fallback to CLI binary (Windows/Linux compatible)
        poses_energies = self._run_vina_python_or_cli(
            receptor_pdbqt, ligand_pdbqt, out_pdbqt, center, box_size
        )

        # Convert output PDBQT poses back to bond-order-preserved RDKit Molecules via Meeko
        pdbqt_mol = PDBQTMolecule.from_file(str(out_pdbqt), skip_typing=True)
        rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

        pose_records: List[Dict[str, Any]] = []
        sdf_writer = Chem.SDWriter(str(out_sdf))
        pose_mols: List[Chem.Mol] = []

        # Extract conformers from Meeko output RDKit mol
        base_mol = rdkit_mols[0] if rdkit_mols else None
        if base_mol is not None:
            for idx, (affinity, rmsd_lb, rmsd_ub) in enumerate(poses_energies):
                if idx >= base_mol.GetNumConformers():
                    break
                thermo = calculate_thermodynamics(affinity, heavy_atoms)
                
                # Create standalone single-conformer RDKit Mol for this pose
                single_mol = Chem.Mol(base_mol)
                single_mol.RemoveAllConformers()
                single_mol.AddConformer(base_mol.GetConformer(idx), assignId=True)

                single_mol.SetProp("_Name", f"Pose_{idx + 1}")
                single_mol.SetProp("Docking_Engine", f"AutoDock_Vina_{self.scoring}")
                single_mol.SetProp("Affinity_kcal_mol", str(thermo["affinity_kcal_mol"]))
                single_mol.SetProp("Estimated_Ki_nM", str(thermo["estimated_ki_nM"]))
                single_mol.SetProp("Ligand_Efficiency", str(thermo["ligand_efficiency"]))
                sdf_writer.write(single_mol)
                pose_mols.append(single_mol)

                pose_records.append({
                    "pose_rank": idx + 1,
                    **thermo,
                    "rmsd_lb": rmsd_lb,
                    "rmsd_ub": rmsd_ub,
                })

        sdf_writer.close()

        return {
            "docked_pdbqt": out_pdbqt,
            "docked_sdf": out_sdf,
            "pose_mols": pose_mols,
            "poses": pose_records,
            "top_pose": pose_records[0] if pose_records else None,
            "grid_center": center,
            "grid_size": box_size,
            "scoring_function": self.scoring,
            "exhaustiveness": self.exhaustiveness,
            "seed": self.seed,
        }

    def _run_vina_python_or_cli(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        out_pdbqt: Path,
        center: Tuple[float, float, float],
        box_size: Tuple[float, float, float],
    ) -> List[Tuple[float, float, float]]:
        try:
            from vina import Vina
            v = Vina(sf_name=self.scoring, cpu=self.cpu_threads, seed=self.seed, verbosity=0)
            v.set_receptor(str(receptor_pdbqt))
            v.set_ligand_from_file(str(ligand_pdbqt))
            v.compute_vina_maps(center=list(center), box_size=list(box_size))
            v.dock(exhaustiveness=self.exhaustiveness, n_poses=self.num_modes)
            v.write_poses(str(out_pdbqt), n_poses=self.num_modes, overwrite=True)
            energies = v.energies(n_poses=self.num_modes)
            return [(float(row[0]), float(row[1]), float(row[2])) for row in energies]
        except ImportError:
            pass

        # CLI Binary Fallback (Windows vina.exe or Linux vina)
        vina_candidates = [
            os.environ.get("VINA_BIN"),
            shutil.which("vina"),
            shutil.which("vina.exe"),
            str(Path("./bin/vina.exe").resolve()),
            str(Path("./bin/vina").resolve()),
            str(Path("H:/AnuDock/bin/vina.exe")),
            "H:/Autodock/vina.exe.EXE",
        ]
        vina_bin = next((c for c in vina_candidates if c and Path(c).is_file()), None)
        if not vina_bin:
            raise RuntimeError(
                "Neither Python 'vina' package nor 'vina' CLI executable found in PATH."
            )

        cmd = [
            vina_bin,
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(box_size[0]),
            "--size_y", str(box_size[1]),
            "--size_z", str(box_size[2]),
            "--scoring", self.scoring,
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--seed", str(self.seed),
            "--out", str(out_pdbqt),
        ]
        if self.cpu_threads > 0:
            cmd.extend(["--cpu", str(self.cpu_threads)])

        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not out_pdbqt.exists():
            raise RuntimeError(f"Vina CLI execution failed:\n{res.stderr}\n{res.stdout}")

        # Parse energies from PDBQT REMARK VINA RESULT lines
        energies: List[Tuple[float, float, float]] = []
        with open(out_pdbqt, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.replace("REMARK VINA RESULT:", "").split()
                    if len(parts) >= 3:
                        energies.append((float(parts[0]), float(parts[1]), float(parts[2])))
        return energies
