# Bindora Dock: Biophysical Structure-Based Molecular Docking & Validation Platform

[![Bindora Fast Benchmark Suite](https://github.com/Swelo-ui/Bindora/actions/workflows/fast_benchmark.yml/badge.svg)](https://github.com/Swelo-ui/Bindora/actions/workflows/fast_benchmark.yml)
![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-blue.svg)
![Scoring: Vinardo / Vina](https://img.shields.io/badge/Engine-AutoDock%20Vina%20%7C%20Vinardo-success.svg)
![Validation: PoseBusters](https://img.shields.io/badge/Validation-PoseBusters%20%2B%20spyrmsd-purple.svg)

**Bindora Dock** is an automated, end-to-end biophysical molecular docking and validation platform designed for reproducible *in silico* structure-based drug design (SBDD). Unlike naive docking wrappers, Bindora enforces strict chemical and structural preparation standards prior to conformational search and validates every docked pose for physical plausibility and symmetry-corrected heavy-atom RMSD.

---

## Key Scientific Capabilities

1. **Biophysical Receptor Sanitization (`bindora/fetcher.py` & `bindora/prep_engine.py`)**
   - Automatically strips crystallization artifacts (`SO4`, `GOL`, `EDO`, `PEG`, bulk solvent) while **preserving essential catalytic metal ions** ($\text{Zn}^{2+}$, $\text{Mg}^{2+}$, $\text{Fe}^{2+}$) and **cofactors** (`HEM`, `FAD`, `NAD`, `ATP`).
   - Reconstructs missing heavy atoms via **PDBFixer (OpenMM)** and assigns pH-dependent protonation microstates ($\text{pH } 7.4$) using **PDB2PQR / PROPKA3** and **OpenBabel** Gasteiger partial charges.
   - Supports **Co-Crystallized Pocket Extraction**, **Residue-Guided Active Site Targeting** (e.g., `['A:25', 'B:25']`), and **Apo Protein Centroid Detection**.

2. **Deterministic Ligand Standardization & Torsion Tree Generation**
   - Strips counter-ions/salts (`LargestFragmentChooser`), normalizes functional group tautomers, and enumerates physiological ionization states via **Dimorphite-DL**.
   - Generates deterministic 3D conformers using **RDKit ETKDGv3 (`seed=42`)**, relaxes internal strain energy via **MMFF94s / UFF** force-field minimization, and builds flexible torsion trees with **Meeko**.

3. **Radius-of-Gyration ($R_g$) Dynamic Box Scaling & Vinardo/Vina Search (`bindora/docking.py`)**
   - Eliminates arbitrary grid box truncation by dynamically scaling cubic search dimensions from the ligand's 3D Radius of Gyration ($L = 2.857 \times R_g + 8.0\text{ \AA}$, Feinstein & Brylinski heuristic).
   - Computes full thermodynamic profiles for every binding mode:
     - **Binding Free Energy ($\Delta G$)** in $\text{kcal/mol}$
     - **Estimated Inhibition Constant ($K_i$)** at $298.15\text{ K}$ ($\text{nM}$)
     - **Ligand Efficiency ($\text{LE} = -\Delta G / N_{\text{heavy}}$)**

4. **Symmetry-Corrected RMSD & PoseBusters Physical Validation (`bindora/report_emitter.py`)**
   - Evaluates redocking accuracy from 1D SMILES without coordinate leakage using graph-isomorphism symmetry-corrected heavy-atom RMSD (**`spyrmsd`** + MCS isomorphism fallback).
   - Screens docked poses with **PoseBusters** and intermolecular Van der Waals clash checks (minimum heavy-atom contact distance $\ge 2.0\text{ \AA}$).
   - Automatically emits a **manuscript-ready Materials & Methods section** (`.md`) and structured JSON telemetry (`.json`).

---

## Verified Redocking Benchmark Performance

All benchmarks are executed strictly **ab initio from 1D SMILES strings** (`seed=42`, `exhaustiveness=8`, `scoring=vinardo`), ensuring zero crystal coordinate bias:

| PDB ID | Target Protein | Ligand (`ResName`) | Rotatable Bonds | Grid Box ($R_g$-Scaled) | Top-1 Affinity ($\Delta G$) | Estimated $K_i$ ($298\text{ K}$) | Symmetry-Corrected RMSD | PoseBusters Validity | Status |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`1HSG`** | HIV-1 Protease Homodimer | Indinavir (`MK1`) | 14 | $23.88 \times 23.88 \times 23.88\text{ \AA}$ | **$-8.65\text{ kcal/mol}$** | $455.3\text{ nM}$ | **0.901 Å** | **PASS (`True`)** | **PASS (< 2.0 Å)** |
| **`1AQ1`** | Human CDK2 Kinase | Staurosporine (`STU`) | 2 | $18.80 \times 18.80 \times 18.80\text{ \AA}$ | **$-9.00\text{ kcal/mol}$** | $250.1\text{ nM}$ | **1.369 Å** | **PASS (`True`)** | **PASS (< 2.0 Å)** |

---

## Quickstart

### 1. Create the Scientific Environment
```bash
mamba env create -f environment.yml
mamba activate bindora-env
```

### 2. Run the Automated Redocking Benchmark Suite

```bash
python scripts/run_fast_benchmark.py
```

### 3. Programmatic Python API Usage

```python
from backend.services.bindora_bridge import run_bindora_scientific_pipeline

result = run_bindora_scientific_pipeline(
    pdb_id="1HSG",
    smiles="CC(C)(C)NC(=O)[C@@H]1CN(Cc2cccnc2)CCN1C[C@@H](O)C[C@@H](Cc1ccccc1)C(=O)N[C@@H]1[C@H](O)Cc2ccccc21",
    ligand_name="MK1",
    co_crystal_resname="MK1",
    scoring_function="vinardo",
    exhaustiveness=8,
    seed=42,
)

print("Top-1 Affinity:", result["top_pose"]["affinity_kcal_mol"], "kcal/mol")
print("Symmetry RMSD:", result["symmetry_rmsd_A"], "Å")
print("PoseBusters Valid:", result["posebusters_validation"]["posebusters_valid"])
```
