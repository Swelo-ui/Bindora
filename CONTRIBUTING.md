# Contributing to Bindora Dock

Thank you for your interest in contributing to **Bindora Dock**, an open-source, research-grade biophysical molecular docking, induced-fit modeling, and computational pharmacology platform developed under **NexPharmaTech**.

Whether you are fixing a bug, contributing new crystallographic benchmark complexes, enhancing scoring physics, expanding chemoinformatics rules, or improving the 3D WebGL studio, these guidelines ensure scientific rigor, reproducibility, and high software quality.

---

## 1. Scientific Integrity Policy

Bindora Dock is designed for publication-level computational pharmacology. All benchmark and analytical data must adhere to strict scientific integrity standards:

1. **Zero Synthetic / Cherry-Picked Data:**
   - Every reported affinity metric ($\Delta G$), estimated $K_i$, and RMSD must be computed through live execution of AutoDock Vina / Vinardo or biophysical engines.
   - Hand-authored, fabricated, or cherry-picked benchmark JSON/Markdown outputs are strictly prohibited.
   - All benchmark scripts report **all** complexes evaluated (both successes and failures) without filtering.

2. **Ab Initio 1D SMILES Conformation Generation:**
   - Docking simulations and redocking benchmarks must generate initial 3D conformers purely from 1D SMILES strings (`RDKit ETKDGv3` / `Meeko`).
   - Coordinate leakage from holo crystal ligands into the search starting point is forbidden.

3. **Topology-Aware Chemical Symmetry RMSD:**
   - Conformational fidelity must be evaluated using graph-isomorphic automorphism search (`spyrmsd`) and template bond-order assignment (`AssignBondOrdersFromTemplate`).
   - Symmetry RMSD must be computed strictly in place within the binding pocket coordinate frame without unconstrained rigid-body alignment.

4. **Cryptographic Provenance:**
   - Benchmark outputs must be serialized through `backend/utils/report_emitter.py` or `bindora/report_emitter.py`, producing atomic records with hardware metadata and SHA-256 digital signatures.
   - Modifying benchmark output files invalidates the checksum signature.

5. **Mandatory DUD-E Analogue-Bias Caveat:**
   - Virtual screening evaluations against the DUD-E dataset must display the standard analogue-bias caveat. `emit_screening_report()` enforces this requirement programmatically.

---

## 2. Environment Setup

### Option A: Conda / Mamba Environment (Recommended for Scientific Workflows)

We provide a streamlined `environment.yml` anchored to `conda-forge`:

```bash
# Clone the repository
git clone https://github.com/Swelo-ui/Bindora.git
cd Bindora

# Create environment via micromamba or conda
micromamba env create -f environment.yml
micromamba activate bindora-env
```

### Option B: Python Virtual Environment (Standard Pip)

Ensure Python 3.10+ (tested on Python 3.11, 3.12, and 3.13) is installed:

```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# Linux / macOS:
source .venv/bin/activate

pip install --upgrade pip
pip install -r requirements.txt
```

### AutoDock Vina Binary Setup

- **Windows:** The official 64-bit binary is pre-bundled in `bin/vina.exe`.
- **Linux / macOS:** Download the official Scripps binary:
  ```bash
  mkdir -p bin
  curl -L -o bin/vina https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
  chmod +x bin/vina
  export PATH="$PWD/bin:$PATH"
  ```

---

## 3. Reproducing the Benchmark Suites

### 3.1. Fast Benchmark Suite (CI Validation)
The deterministic fast benchmark suite evaluates HIV-1 Protease (`1HSG`) and Human CDK2 Kinase (`1AQ1`) ab initio from 1D SMILES, reporting symmetry-corrected RMSD, PoseBusters biophysical checks, and free energy estimates:

```bash
python scripts/run_fast_benchmark.py
```

### 3.2. Extended Multi-Family Suite
Runs a 4-family validation suite covering Aspartyl Protease (`1HSG`), Ser/Thr Kinase (`1AQ1`), Nuclear Receptor (`3ERT`), and $\text{Zn}^{2+}$ Metalloenzyme (`1BCD`):

```bash
python scripts/run_fast_benchmark.py --extended
```

### 3.3. CASF-2016 Core Set Benchmark (285 Complexes)
The internationally recognized CASF-2016 benchmark evaluates redocking success rate across 285 diverse crystallographic targets:

```bash
# Smoke test (5 complexes)
python tests/benchmark_casf2016.py --max-complexes 5 --exhaustiveness 4

# Full benchmark run (resumable)
python tests/benchmark_casf2016.py --exhaustiveness 8 --resume
```

Outputs are automatically saved to `data/benchmarks/casf2016_final_report.json` and `casf2016_final_report.md`.

### 3.4. DUD-E Virtual Screening Enrichment
Evaluates enrichment factors (EF1%, EF5%, EF10%) and receiver operating characteristic curves (ROC-AUC):

```bash
# Diverse 8-target subset
python tests/benchmark_screening.py --subset diverse --exhaustiveness 4

# Single target smoke test
python tests/benchmark_screening.py --target ache --exhaustiveness 4
```

---

## 4. Running Tests & Quality Assurance

Before submitting a Pull Request, verify that all automated unit and integration tests pass:

```bash
# Run core scientific integrity tests
python -m pytest tests/test_scientific_integrity.py -v

# Run docking and ADME test suites
python -m pytest tests/test_docking.py tests/test_adme.py -v

# Run full test suite
python -m pytest tests/ -v
```

---

## 5. Contributing Code & Guidelines

### 5.1. Pull Request Checklist
When submitting a Pull Request:
- [ ] **Tests Pass:** All tests pass with zero errors (`pytest tests/`).
- [ ] **Fast Benchmark Pass:** Fast benchmark finishes with status PASS for all targets (`python scripts/run_fast_benchmark.py`).
- [ ] **Type Annotations:** Public service functions and classes include Python type annotations (`typing`).
- [ ] **Deterministic Seeds:** Do not hardcode silent fixed seeds in production services. Use `secrets.randbelow()` for production runs or explicit parameters for reproducible benchmarks.
- [ ] **Clean Markdown:** Ensure all Ångström symbols in documentation use standard UTF-8 characters (`Å`, `Å²`, `Å³`) rather than raw LaTeX `\AA` to ensure flawless rendering across GitHub and Web browsers.
- [ ] **Preserve Comments:** Keep scientific rationale, physical units, and methodology citations documented in source code.

### 5.2. Adding a New Benchmark Target
To add a new benchmark target to `scripts/run_fast_benchmark.py` or `tests/benchmark_accuracy.py`:
```python
{
    "pdb_id": "XXXX",
    "ligand_resname": "LIG",
    "family": "Kinase",
    "name": "Target Full Name",
    "smiles": "Canonical_1D_SMILES_String",
    "target_chain": "A",
    "rmsd_threshold": 2.0,
}
```
Verify the benchmark script executes and passes before committing.

### 5.3. Web Studio & Frontend Contributions
- The frontend uses vanilla HTML5, CSS3 with modern glassmorphism design tokens, and WebGL via `3Dmol.js`.
- Ensure all UI changes are responsive, support both Dark and Light color schemes seamlessly, and do not introduce unstyled image assets or broken layout paddings.

---

## 6. Reporting Issues & Vulnerabilities

### Bug Reports
- Open a GitHub Issue specifying:
  - Operating system and Python version.
  - PDB ID or ligand SMILES string causing the error.
  - Full traceback or log output from `logs/` or console.
  - Expected versus observed behavior.

### Security Disclosures
- For responsible security disclosures (e.g. arbitrary code execution, server vulnerabilities), please contact the engineering team directly at **`sharmaji.pharmatech.info@gmail.com`** instead of filing a public issue.

---

## 7. License

By contributing to Bindora Dock, you agree that your contributions are licensed under the **Apache License, Version 2.0**. See the [LICENSE](LICENSE) file for terms and conditions.

*Bindora Dock — Research-grade molecular docking and computational pharmacology, openly benchmarked.*
