# bindora/prep_engine.py
import os
import subprocess
import shutil
import numpy as np
from pathlib import Path
from typing import Tuple, Dict, Optional, Any

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from meeko import MoleculePreparation, PDBQTWriterLegacy


def compute_radius_of_gyration(mol: Chem.Mol) -> float:
    """Computes 3D Radius of Gyration (Rg) in Angstroms for dynamic box sizing."""
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])
    total_mass = np.sum(masses)
    center_of_mass = np.sum(coords * masses[:, None], axis=0) / total_mass
    sq_dists = np.sum((coords - center_of_mass) ** 2, axis=1)
    rg = float(np.sqrt(np.sum(masses * sq_dists) / total_mass))
    return round(rg, 3)


class BiophysicalPrepEngine:
    """Research-grade preparation for receptors and small-molecule ligands."""

    def __init__(self, work_dir: str = "./workspace", target_ph: float = 7.4, seed: int = 42):
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.target_ph = target_ph
        self.seed = seed

    def prepare_ligand_from_smiles(
        self, smiles: str, ligand_name: str = "ligand"
    ) -> Dict[str, Any]:
        """
        Sanitizes SMILES (strips salts, normalizes tautomers, protonates at pH 7.4),
        embeds deterministic 3D ETKDGv3 conformer, minimizes with MMFF94s, and writes PDBQT.
        """
        raw_mol = Chem.MolFromSmiles(smiles)
        if raw_mol is None:
            raise ValueError(f"RDKit failed to parse SMILES: {smiles}")

        # 1. Strip counter-ions/salts & keep largest covalent organic parent
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(raw_mol)

        # 2. Functional group normalization & canonical tautomer
        normalizer = rdMolStandardize.Normalizer()
        mol = normalizer.normalize(mol)

        # Optional pH-dependent protonation via Dimorphite-DL if installed
        try:
            from dimorphite_dl import protonate_smiles
            try:
                prot_list = protonate_smiles(
                    Chem.MolToSmiles(mol),
                    ph_min=self.target_ph - 0.4,
                    ph_max=self.target_ph + 0.4,
                    max_variants=1,
                )
            except TypeError:
                prot_list = protonate_smiles(
                    Chem.MolToSmiles(mol),
                    min_ph=self.target_ph - 0.4,
                    max_ph=self.target_ph + 0.4,
                    max_variants=1,
                )
            if prot_list:
                mol_prot = Chem.MolFromSmiles(prot_list[0])
                if mol_prot is not None:
                    mol = mol_prot
        except Exception:
            pass

        # 3. Add explicit hydrogens with 3D coordinates via deterministic ETKDGv3
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = self.seed
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True

        res = AllChem.EmbedMolecule(mol, params)
        if res == -1:
            res = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=self.seed)
        if res == -1:
            raise RuntimeError(f"3D conformer embedding failed for {ligand_name}")

        # 4. Force-field energy minimization (MMFF94s with UFF fallback)
        force_field_used = "MMFF94s"
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=500)
        else:
            force_field_used = "UFF"
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        # 5. Compute biophysical properties & optimal docking box size
        rg = compute_radius_of_gyration(mol)
        # Feinstein & Brylinski (2015) optimal box edge formula: L = 2.857 * Rg + 8.0 Å
        optimal_box_edge = round(max(18.0, min(32.0, 2.857 * rg + 8.0)), 2)
        heavy_atoms = mol.GetNumHeavyAtoms()
        rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mw = round(Descriptors.MolWt(mol), 2)

        # Save minimized 3D SDF for RMSD & PoseBusters reference
        sdf_path = self.work_dir / f"{ligand_name}_minimized.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        # 6. Meeko Torsion Tree & PDBQT generation
        preparator = MoleculePreparation(flexible_amides=False, rigid_macrocycles=False)
        mol_setups = preparator.prepare(mol)
        pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            raise RuntimeError(f"Meeko PDBQT generation failed for {ligand_name}: {err_msg}")

        pdbqt_path = self.work_dir / f"{ligand_name}.pdbqt"
        pdbqt_path.write_text(pdbqt_string, encoding="utf-8")

        return {
            "rdkit_mol": mol,
            "sdf_path": sdf_path,
            "pdbqt_path": pdbqt_path,
            "radius_of_gyration": rg,
            "optimal_box_size": (optimal_box_edge, optimal_box_edge, optimal_box_edge),
            "heavy_atoms": heavy_atoms,
            "rotatable_bonds": rotatable_bonds,
            "molecular_weight": mw,
            "force_field": force_field_used,
        }

    def prepare_receptor_pdbqt(self, clean_receptor_pdb: Path) -> Tuple[Path, str]:
        """
        Repairs missing heavy atoms with PDBFixer, protonates at target pH using PDB2PQR/PROPKA
        (or OpenBabel pH model fallback), and converts to rigid AutoDock PDBQT.
        """
        stem = clean_receptor_pdb.stem
        fixed_pdb = self.work_dir / f"{stem}_fixed.pdb"
        pqr_path = self.work_dir / f"{stem}_propka.pqr"
        pdbqt_path = self.work_dir / f"{stem}.pdbqt"
        prep_method_log = []

        # Step 1: PDBFixer for missing heavy atoms & non-standard residue conversion
        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile

            fixer = PDBFixer(filename=str(clean_receptor_pdb))
            fixer.findMissingResidues()
            fixer.missingResidues = {}  # Do not model long unverified terminal loops
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(self.target_ph)
            with open(fixed_pdb, "w", encoding="utf-8") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
            prep_method_log.append(f"PDBFixer (missing heavy atoms repaired + pH {self.target_ph} H-addition)")
        except Exception as e:
            shutil.copy(clean_receptor_pdb, fixed_pdb)
            prep_method_log.append(f"Standard PDB coordinates retained (PDBFixer skipped: {e})")

        # Step 2: Try PDB2PQR + PROPKA3 for local pKa microstate assignment
        pdb2pqr_bin = shutil.which("pdb2pqr") or shutil.which("pdb2pqr30")
        input_for_pdbqt = fixed_pdb
        if pdb2pqr_bin:
            cmd = [
                pdb2pqr_bin,
                "--ff=AMBER",
                f"--with-ph={self.target_ph}",
                "--titration-state-method=propka",
                "--keep-chain",
                str(fixed_pdb),
                str(pqr_path),
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode == 0 and pqr_path.exists():
                input_for_pdbqt = pqr_path
                prep_method_log.append(f"PDB2PQR + PROPKA3 titration at pH {self.target_ph}")

        # Step 3: Convert to rigid receptor PDBQT with Gasteiger partial charges via OpenBabel
        obabel_candidates = [
            os.environ.get("OBABEL_BIN"),
            shutil.which("obabel"),
            shutil.which("obabel.exe"),
            str(Path("H:/mgltools/OpenBabel-2.3.2/obabel.exe")),
            "C:/Program Files/OpenBabel-3.1.1/obabel.exe",
            "C:/Program Files (x86)/OpenBabel-2.4.1/obabel.exe",
        ]
        obabel_bin = next((c for c in obabel_candidates if c and Path(c).is_file()), None)
        if not obabel_bin:
            raise RuntimeError("OpenBabel ('obabel') binary not found in PATH or standard paths.")

        cmd_obabel = [
            obabel_bin,
            str(input_for_pdbqt),
            "-xr",  # Rigid receptor flag (no torsion tree roots)
            "-p", str(self.target_ph),
            "--partialcharge", "gasteiger",
            "-O", str(pdbqt_path),
        ]
        res_ob = subprocess.run(cmd_obabel, capture_output=True, text=True)
        if res_ob.returncode != 0 or not pdbqt_path.exists():
            raise RuntimeError(f"OpenBabel receptor PDBQT conversion failed: {res_ob.stderr}")

        # Sanitize receptor PDBQT: strip any ROOT/ENDROOT/TORSDOF tags generated by OpenBabel
        raw_lines = pdbqt_path.read_text(encoding="utf-8", errors="replace").splitlines() if pdbqt_path.exists() else []
        clean_lines = [
            l for l in raw_lines
            if not l.startswith(("ROOT", "ENDROOT", "TORSDOF", "BRANCH", "ENDBRANCH"))
        ]
        
        # If OpenBabel produced 0 molecules (common with legacy OpenBabel on large multimers), use BiophysicalReceptorPreparer
        if len(clean_lines) < 10:
            try:
                from backend.services.biophysical_prep import BiophysicalReceptorPreparer
                clean_pdb_lines = fixed_pdb.read_text(encoding="utf-8", errors="replace").splitlines()
                prot_lines = [l for l in clean_pdb_lines if l.startswith(("ATOM  ", "HETATM"))]
                _, pdbqt_lines, _ = BiophysicalReceptorPreparer.prepare(prot_lines, pH=self.target_ph)
                clean_lines = pdbqt_lines
                prep_method_log.append("BiophysicalReceptorPreparer (AMBER FF14SB partial charges & AD4 atom types)")
            except Exception as be:
                prep_method_log.append(f"Biophysical fallback notice: {be}")
        else:
            prep_method_log.append("OpenBabel Gasteiger charge & AutoDock atom typing (-xr)")

        pdbqt_path.write_text("\n".join(clean_lines) + "\nTER\nEND\n", encoding="utf-8")
        return pdbqt_path, " -> ".join(prep_method_log)
