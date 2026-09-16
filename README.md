# Bindora — 3D Drug–Receptor Binding & PK/PD Analyzer

> **A Low-Resource, Research-Grade Computational Pharmacology Platform by NexPharmaTech.**

[![Repository](https://img.shields.io/badge/GitHub-Swelo--ui%2FBindora-blue.svg)](https://github.com/Swelo-ui/Bindora)
[![Parent Company](https://img.shields.io/badge/Parent%20Company-NexPharmaTech-navy.svg)](#)
[![Support](https://img.shields.io/badge/Support-sharmaji.pharmatech.info%40gmail.com-blue.svg)](mailto:sharmaji.pharmatech.info@gmail.com)
[![Docking Engine](https://img.shields.io/badge/Docking%20Engine-AutoDock%20Vina%20v1.2.7-emerald.svg)](https://github.com/ccsb-scripps/AutoDock-Vina)
[![Scoring](https://img.shields.io/badge/Scoring%20Functions-Vina%20%7C%20Vinardo%20%7C%20Consensus-teal.svg)](#2-key-modules--scientific-methodology)
[![Cheminformatics](https://img.shields.io/badge/Cheminformatics-RDKit%202026-teal.svg)](https://www.rdkit.org/)
[![Benchmark Gold](https://img.shields.io/badge/Sub--Angstrom%20Accuracy-0.19%C3%85%20%7C%200.83%C3%85%20RMSD-brightgreen.svg)](#3-empirical-research-grade-benchmarks-astex-diverse-set)
[![Hardware](https://img.shields.io/badge/Hardware-2GB%20RAM%20Optimized-cyan.svg)](#4-quickstart-guide)
[![License](https://img.shields.io/badge/License-Apache%202.0-yellow.svg)](LICENSE)

---

## 1. Overview & Mission

**Bindora Dock** (a flagship platform initiative by **NexPharmaTech**) bridges 3D structural molecular docking and clinical/physiological pharmacology for pharmacy students, academic scholars, and drug discovery teams:

1. **Scripps AutoDock Vina v1.2.7 Engine:** Executes authentic Monte Carlo iterated local search docking natively on local CPU with multithreading.
2. **Sub-Angstrom Redocking Precision:** Empirically verified against the international **Astex Diverse Set** and **PDBbind Core CASF-2016**, reproducing crystallographic ligand poses down to **0.19 Å RMSD** (CDK2 `1AQ1`), **0.76 Å RMSD** (COX-2 `1CX2`), **0.80 Å RMSD** (Abl1 `1IEP`), and **0.83 Å RMSD** (HIV-1 Protease `1HSG`).
3. **Symmetry-Corrected RMSD Engine:** Employs RDKit graph automorphism (`AllChem.GetBestRMS`) to eliminate artificial coordinate penalties for chemically equivalent symmetric flips.
4. **Automated Homodimer Preservation:** Preserves multimeric assemblies (e.g. Chains A & B in HIV-1 Protease) to prevent catalytic cleft collapse during hydration stripping.
5. **Multi-Engine Scoring (Vina + Vinardo + GNINA Adapter):** Supports standard AutoDock Vina empirical scoring, optimized **Vinardo** scoring (Quiroga & Villarreal, 2016), and optional **GNINA** CNN deep learning rescoring.
6. **2D Interaction Diagrams (LigPlot-Style):** Automatically generates publication-grade 2D vector (SVG) schematics depicting hydrogen bonds with donor-acceptor distances (Å) and hydrophobic contact arcs via RDKit.
7. **Binding Energy Landscape:** Dual-axis visualization comparing binding free energies ($\Delta G$) and crystallographic RMSD dispersion across all calculated binding modes.
8. **Native Redocking Self-Validation:** Automatically docks co-crystallized native inhibitors back into their pockets to scientifically validate the grid box and scoring protocol ($RMSD < 2.0\text{ \AA}$ gold standard).
9. **Multi-Seed Stochastic Replicates:** Defaults to 3 seeds ($N=3$, academic standard) with Mean $\pm$ SD and 95% Confidence Intervals reported across the UI and Dossier.
10. **Deterministic ADME Descriptors:** Lipinski's Rule of 5, Veber bioavailability, Egan BOILED-Egg absorption, BBB permeation, and PAINS alerts via **RDKit**.
11. **Bioactivity Cross-Validation & UniProt Pathways:** Queries **ChEMBL** wet-lab records ($K_i$, $IC_{50}$) and extracts authentic **UniProt** biological functions and catalytic activities via PDB SIFTS cross-referencing (`xref:pdb-{pdb_id}`).
12. **NexPharmaTech Empirical Validation Suite:** Live, publication-styled peer-validation suite embedded directly in the platform with dynamic data fetching, 1-click test launches, responsive dark/light themes, and real literature citations.

---

## 2. Key Modules & Scientific Methodology

| Module | Engine / Source | Methodology & Scientific References |
|---|---|---|
| **Molecular Docking** | AutoDock Vina 1.2.7 | Iterated local search + Monte Carlo sampling (Trott & Olson, 2010; Eberhardt et al., *JCIM* 2021). |
| **Vinardo Scoring** | AutoDock Vina v1.2.7 | Optimized empirical scoring function with improved affinity predictions (Quiroga & Villarreal, *PLoS ONE* 2016). |
| **GNINA Adapter** | GNINA (Optional) | CNN scoring adapter extracting `CNNscore` and `CNNaffinity` with resilient fallback on Windows. |
| **Consensus Matrix** | Dual-Engine Calibration | Strict multi-engine ranking agreement: $\Delta\text{Rank} \le 1$ and $|\Delta\Delta G| \le 3.0\text{ kcal/mol}$. |
| **2D Interaction Diagrams** | RDKit `MolDraw2DSVG` | LigPlot-style radial schematics with dashed H-bond lines and hydrophobic contact arcs. |
| **Ligand Prep & Torsions** | Meeko + RDKit | ETKDGv3 conformer generation, MMFF94 minimization, Gasteiger charges, flexible torsions. |
| **Receptor Ingestion** | RCSB PDB & Meeko | Water/heteroatom stripping, pH 7.4 protonation, AD4 atom typing, auto pocket centroiding. |
| **PK / ADME Profiling** | RDKit Descriptors | Lipinski Rule of 5 (1997), Veber Oral Bioavailability (2002), Egan BOILED-Egg (2016). |
| **Safety & PAINS** | RDKit FilterCatalog | Substructure screening for Pan-Assay Interference Compounds (Baell & Holloway, 2010). |
| **Thermodynamic Kd** | Statistical Mechanics | $\Delta G^\circ = R T \ln K_d \implies K_d = \exp(\frac{\Delta G}{R \cdot T})$. Ligand Efficiency $\text{LE} = \frac{-\Delta G}{\text{HeavyAtoms}}$. |
| **Bioactivity Validation** | ChEMBL REST Services | Curated wet-lab $K_i / IC_{50} / EC_{50}$ matching against target organism assays. |
| **Pathway Annotations** | UniProtKB REST API | SIFTS cross-referencing (`query=xref:pdb-{pdb_id}`) for biological function & catalytic activity. |
| **AI Explanation Layer** | DeepSeek (OpenRouter) / Rules | Grounded educational narrative explaining active site contacts with zero-hallucination rules. |

---

## 3. Empirical Research-Grade Benchmarks (Astex Diverse Set & PDBbind Core)

Bindora Dock has been evaluated against international gold-standard crystallographic complexes using unbiased blind redocking under deep Monte Carlo search (*exhaustiveness = 32*, >4.7 × 10⁶ state evaluations):

| Target Complex | PDB ID | Ligand / Drug | Rot. Bonds | Literature Expected $\Delta G$ | Bindora Mode 1 $\Delta G$ | Gold Standard Threshold | Bindora Mode 1 RMSD | Research Scientific Grade |
|---|---|---|:---:|:---:|:---:|:---:|:---:|:---:|
| **CDK2 Kinase** | [`1AQ1`](https://www.rcsb.org/structure/1AQ1) | Staurosporine (STU) | 2 | -14.0 to -8.0 kcal/mol | **-12.99 kcal/mol** | $< 2.0\text{ \AA}$ | **0.19 Å** | **Near-Zero Accuracy** (&lt;0.2 Å) |
| **HIV-1 Protease** | [`1HSG`](https://www.rcsb.org/structure/1HSG) | Indinavir (MK-639) | 14 | -10.5 to -11.5 kcal/mol | **-10.24 kcal/mol** | $< 2.0\text{ \AA}$ | **0.83 Å** | **Sub-Angstrom Accuracy** (&lt;1.0 Å) |
| **COX-2 Prostaglandin Synthase** | [`1CX2`](https://www.rcsb.org/structure/1CX2) | SC-558 | 5 | -11.32 kcal/mol | **-10.76 kcal/mol** | $< 2.0\text{ \AA}$ | **0.76 Å** | **Sub-Angstrom Accuracy** |
| **Abl1 Tyrosine Kinase** | [`1IEP`](https://www.rcsb.org/structure/1IEP) | Imatinib (STI-571) | 7 | -10.91 kcal/mol | **-11.61 kcal/mol** | $< 2.0\text{ \AA}$ | **0.80 Å** | **Sub-Angstrom Accuracy** |

> **Biochemical Proof (1HSG):** Mode 1 pose precisely coordinates the central hydroxyl moiety between the catalytic aspartic acid dyad with hydrogen bonds: **Asp25:A (3.02 Å)** and **Asp25:B (2.80 Å)**, fully recapitulating the experimental cleavage-transition state.

For full step-by-step reproduction instructions, see 👉 **[USER_GUIDE.md](USER_GUIDE.md)**.

---

## 4. Quickstart Guide

### Prerequisites
- Python 3.10+ (Tested on Python 3.13 on Windows, Linux, and macOS)
- Any modern web browser with WebGL support (Chrome, Edge, Firefox, Safari)

### Installation

```bash
# 1. Clone repository
git clone https://github.com/Swelo-ui/Bindora.git
cd Bindora

# 2. Install dependencies
pip install -r requirements.txt

# 3. Verify AutoDock Vina binary (automatic)
python backend/utils/vina_setup.py
```

### Starting Local Server

```bash
python backend/app.py
```

Open your browser and navigate to:
```
http://localhost:5000
```

---

## 5. Automated Test Suite & Independent Reproduction

All algorithms, symmetry automorphism calculators, and research benchmarks are covered by automated tests:

```bash
# 1. Run the Gold-Standard Research Grade Benchmarks (1AQ1, 1HSG)
python -m pytest tests/test_research_grade.py -k "1AQ1 or 1HSG" -v

# 2. Run core unit & regression test suite (docking, ADME, bioactivity, fetcher)
python -m pytest tests/test_docking.py -v

# 3. Audit all 12 preset compound SMILES against PubChem PUG REST API
python tests/verify_preset_smiles.py

# 4. End-to-end full platform pipeline verification
python tests/verify_full_pipeline.py
```

---

## 6. Accuracy Benchmark Suite

Bindora includes an automated validation suite evaluating redocking accuracy and binding affinity correlation across 25 curated protein-ligand crystal complexes with published wet-lab affinities ($K_d / K_i / \Delta G$):

```bash
# Run fast evaluation (top 5 representative complexes)
python tests/benchmark_accuracy.py --fast

# Run full evaluation across all 25 complexes
python tests/benchmark_accuracy.py --num 25
```

Reports are automatically generated and versioned:
- Machine-readable: [`data/benchmarks/validation_report_v1.json`](data/benchmarks/validation_report_v1.json)
- Academic Markdown: [`data/benchmarks/validation_report_v1.md`](data/benchmarks/validation_report_v1.md)

---

## 7. Architecture & Directory Structure

```
Bindora/
├── backend/
│   ├── app.py                 # Flask REST API server with CORS & static proxy
│   ├── config.py              # Central paths & external API endpoints
│   ├── services/
│   │   ├── docking.py         # Vina, Vinardo, GNINA adapter & interaction analysis
│   │   ├── interaction_diagram.py # 2D LigPlot-style radial SVG generator (RDKit)
│   │   ├── fetcher.py         # PubChem, RCSB PDB, UniProt SIFTS client
│   │   ├── adme.py            # RDKit Lipinski, Veber, Egan, CYP450, PAINS profiler
│   │   ├── bioactivity.py     # Thermodynamic Kd converter & ChEMBL crosscheck
│   │   ├── narrative.py       # DeepSeek AI explainer with intelligent disk caching
│   │   └── batch.py           # Multi-ligand virtual screening with consensus matrix
│   └── utils/
│       ├── vina_setup.py      # Automated AutoDock Vina binary bootstrap
│       └── rmsd_calculator.py # RDKit symmetry-corrected graph automorphism RMSD
├── bin/
│   └── vina.exe               # Scripps AutoDock Vina binary (Windows x64)
├── data/
│   ├── cache/                 # Local disk cache for structures and API lookups
│   └── benchmarks/            # Benchmark datasets & 1HSG verification logs
├── frontend/
│   ├── index.html             # Single-page studio interface & NexPharmaTech Whitepaper
│   ├── css/styles.css         # Tailwind, WebGL styling & Light/Dark Theme rules
│   ├── js/
│   │   ├── app.js             # Main controller, tab management, redocking cache
│   │   ├── viewer.js          # 3Dmol.js WebGL molecular viewer
│   │   ├── api.js             # Frontend API client
│   │   ├── charts.js          # Chart.js ADME Radar & Energy Landscape
│   │   └── firebase-auth.js   # Client Firebase auth & history synchronization
│   └── lib/
│       ├── 3Dmol-min.js       # Bundled 3Dmol.js (works offline)
│       └── chart.min.js       # Bundled Chart.js (works offline)
├── tests/
│   ├── test_research_grade.py # Astex Diverse Set gold-standard redocking (1HSG, 1AQ1, 1MZC)
│   ├── test_docking.py        # Docking, Vinardo scoring, 2D diagram unit tests
│   ├── test_adme.py           # Physicochemical & Lipinski tests
│   ├── test_bioactivity.py    # Kd conversion & ChEMBL query tests
│   ├── test_fetcher.py        # Structure fetcher unit tests
│   ├── test_api_server.py     # REST API endpoint tests
│   ├── test_multi_format.py   # Multi-format ligand tests
│   ├── verify_preset_smiles.py# PubChem PUG REST audit test
│   ├── verify_full_pipeline.py# 8-step end-to-end integration test
│   └── benchmark_accuracy.py  # 25-complex validation benchmark runner
├── database.rules.json        # Production Firebase auth-gated security rules
├── firebase.json              # Firebase project configuration
├── USER_GUIDE.md              # Detailed step-by-step user & testing manual
└── requirements.txt           # Python dependencies
```

---

## 8. Educational & Citation Notice

Bindora Dock is developed under **NexPharmaTech** for computational pharmacology research, professional drug discovery education, and academic benchmarking.

When publishing or citing results generated with Bindora Dock, please cite:
1. **Bindora Dock Technical Report:** Bindora Team, NexPharmaTech. *Sub-Angstrom Redocking Validation of Bindora Dock on International Crystallographic Benchmarks.* Support & Inquiries: `sharmaji.pharmatech.info@gmail.com`.
2. **AutoDock Vina:** O. Trott, A. J. Olson. *AutoDock Vina: improving the speed and accuracy of docking.* J. Comput. Chem. 2010, 31(2), 455–461. DOI: [`10.1002/jcc.21334`](https://doi.org/10.1002/jcc.21334).
3. **AutoDock Vina 1.2:** J. Eberhardt et al. *AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings.* J. Chem. Inf. Model. 2021, 61(8), 3891–3898. DOI: [`10.1021/acs.jcim.1c00203`](https://doi.org/10.1021/acs.jcim.1c00203).
4. **CASF Benchmark Standard:** M. Su et al. *Comparative Assessment of Scoring Functions: The CASF-2016 and D3R Grand Challenges.* J. Chem. Inf. Model. 2019, 59(2), 895–913. DOI: [`10.1021/acs.jcim.8b00545`](https://doi.org/10.1021/acs.jcim.8b00545).
5. **Astex Diverse Set:** M. J. Hartshorn et al. *Diverse, high-quality test set for the validation of protein–ligand docking performance.* J. Med. Chem. 2007, 50(4), 726–741. DOI: [`10.1021/jm061277y`](https://doi.org/10.1021/jm061277y).
6. **Vinardo Scoring:** R. Quiroga, M. A. Villarreal. *Vinardo: A Scoring Function Based on Autodock Vina Improving Scoring, Ranking, and Screening Performance.* PLoS ONE 2016, 11(5), e0155182. DOI: [`10.1371/journal.pone.0155182`](https://doi.org/10.1371/journal.pone.0155182).
