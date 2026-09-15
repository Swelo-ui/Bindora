# Bindora — 3D Drug–Receptor Binding & PK/PD Analyzer

> **A Low-Resource, Zero-Hallucination Computational Pharmacology Platform for Pharmacy Students & Researchers.**

[![Repository](https://img.shields.io/badge/GitHub-Swelo--ui%2FBindora-blue.svg)](https://github.com/Swelo-ui/Bindora)
[![Docking Engine](https://img.shields.io/badge/Docking%20Engine-AutoDock%20Vina%20v1.2.7-emerald.svg)](https://github.com/ccsb-scripps/AutoDock-Vina)
[![Cheminformatics](https://img.shields.io/badge/Cheminformatics-RDKit%202026-teal.svg)](https://www.rdkit.org/)
[![AI Engine](https://img.shields.io/badge/OpenRouter-deepseek--v4--flash--0731-purple.svg)](https://openrouter.ai/)
[![Hardware](https://img.shields.io/badge/Hardware-2GB%20RAM%20Optimized-cyan.svg)](#low-resource-optimization)
[![License](https://img.shields.io/badge/License-Apache%202.0-yellow.svg)](LICENSE)

---

## 1. Overview & Mission

**Bindora** bridges 3D structural molecular docking and clinical/physiological pharmacology for pharmacy students, academic scholars, and independent drug researchers without wet-lab access:

1. **Genuine Molecular Docking:** Executes Scripps CCSB **AutoDock Vina v1.2.7** natively on server CPU.
2. **Deterministic ADME Descriptors:** Calculates Lipinski's Rule of 5, Veber oral bioavailability, Egan BOILED-Egg intestinal absorption, Blood-Brain Barrier (BBB) permeation, and PAINS safety alerts via **RDKit**.
3. **Bioactivity Cross-Validation:** Automatically queries EMBL-EBI **ChEMBL** curated wet-lab records to validate computational predictions against real $K_i$, $IC_{50}$, and $EC_{50}$ measurements.
4. **Interactive 3D WebGL Visualization:** Real-time protein ribbons, ligand sticks, active pocket surfaces, and dashed hydrogen-bond contact lines with live distance measurement tools via **3Dmol.js** — specifically optimized for **2 GB RAM** hardware.
5. **AI Research Pharmacologist:** Powered by `deepseek/deepseek-v4-flash-0731` via OpenRouter with strict anti-hallucination guardrails and intelligent disk caching to prevent excess API consumption, with a 100% offline deterministic fallback.

---

## 2. Key Modules & Scientific Methodology

| Module | Engine / Source | Methodology & Scientific References |
|---|---|---|
| **Molecular Docking** | AutoDock Vina 1.2.7 | Empirical scoring function + iterated local search (Eberhardt et al., *J. Chem. Inf. Model.* 2021). |
| **Ligand Prep & Torsions** | Meeko + RDKit | ETKDGv3 3D conformer generation, MMFF94 minimization, Gasteiger charges, flexible rotatable bond assignment. |
| **Structure Ingestion** | PubChem & RCSB PDB REST | Automated 2D/3D structure acquisition with keyless, CC0 / Public Domain public APIs. |
| **PK / ADME Profiling** | RDKit Descriptors | Lipinski Rule of Five (1997), Veber Oral Bioavailability (2002), Ghose Filter (1999), Egan BOILED-Egg (2016). |
| **Safety & PAINS** | RDKit FilterCatalog | Substructure screening for Pan-Assay Interference Compounds (Baell & Holloway, *J. Med. Chem.* 2010) and Brenk toxicophores (2008). |
| **Thermodynamic Conversion** | Statistical Mechanics | $\Delta G^\circ = R T \ln K_d \implies K_d = \exp(\frac{\Delta G \times 1000}{R \cdot T})$. Ligand Efficiency $\text{LE} = \frac{-\Delta G}{\text{HeavyAtoms}}$ (Hopkins et al., 2004). |
| **Bioactivity Validation** | ChEMBL REST Services | Curated wet-lab $K_i / IC_{50} / EC_{50}$ matching against target organism assays. |
| **AI Explanation Layer** | DeepSeek (OpenRouter) / Offline Rules | Grounded educational narrative explaining active site contacts, oral drug-likeness, and certainty levels. |

---

## 3. Preloaded Benchmark Case Studies

Bindora includes instant 1-click exploration for classic drug-target systems:
1. **Imatinib vs BCR-ABL1 Kinase (`1IEP`)**: Type II kinase inhibitor in Chronic Myeloid Leukemia (ChEMBL $IC_{50} \approx 38\text{ nM}$).
2. **Aspirin vs Cyclooxygenase-2 (`1CX2`)**: Classic NSAID targeting active site Ser530 in inflammation and pain.
3. **Gefitinib vs EGFR Kinase (`2ITY`)**: 4-Anilinoquinazoline ATP-competitive inhibitor in Non-Small Cell Lung Cancer (ChEMBL $IC_{50} \approx 33\text{ nM}$).
4. **Remdesivir (GS-441524) vs SARS-CoV-2 RdRp (`7BV2`)**: Antiviral nucleoside analog inhibiting viral RNA synthesis.

---

## 4. Quickstart Guide

### Prerequisites
- Python 3.10+ (Tested on Python 3.13 on Windows / Linux / macOS)
- Web browser with WebGL support

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/Swelo-ui/Bindora.git
cd Bindora

# 2. Install Python dependencies
pip install -r requirements.txt

# 3. Setup environment secrets (optional for AI narration)
cp .env.example .env
# Edit .env to add your OPENROUTER_API_KEY

# 4. Setup and verify AutoDock Vina binary (automatic)
python backend/utils/vina_setup.py
```

### Running the Application

```bash
python backend/app.py
```

Open your browser and navigate to:
```
http://localhost:5000
```

---

## 5. Running the Automated Test Suite

All algorithms and API endpoints are covered by automated unit and integration tests:

```bash
python -m pytest tests/ -v
```

Output:
```
tests/test_adme.py::test_adme_aspirin PASSED
tests/test_adme.py::test_adme_invalid_smiles PASSED
tests/test_bioactivity.py::test_thermodynamics_conversion PASSED
tests/test_bioactivity.py::test_chembl_crosscheck_imatinib PASSED
tests/test_docking.py::test_ligand_preparation PASSED
tests/test_docking.py::test_receptor_preparation_and_docking PASSED
tests/test_fetcher.py::test_pubchem_search_aspirin PASSED
tests/test_fetcher.py::test_rcsb_pdb_fetch_1cx2 PASSED
tests/test_fetcher.py::test_uniprot_search PASSED
tests/test_api_server.py::test_api_health PASSED
tests/test_api_server.py::test_api_benchmarks PASSED
tests/test_api_server.py::test_api_pubchem_search PASSED
tests/test_api_server.py::test_api_adme_calculation PASSED
============================== 13 passed in 5.4s ==============================
```

---

## 6. Architecture & Directory Tree

```
Bindora/
├── backend/
│   ├── app.py                 # Flask entrypoint with CORS & REST API routes
│   ├── config.py              # Configuration & public API endpoints
│   ├── services/
│   │   ├── fetcher.py         # PubChem, RCSB PDB, UniProt API clients
│   │   ├── docking.py         # AutoDock Vina execution & contact analyzer
│   │   ├── adme.py            # RDKit Lipinski, Veber, Egan, BBB, PPB, PAINS profiler
│   │   ├── bioactivity.py     # Thermodynamic Kd converter & ChEMBL crosscheck
│   │   ├── narrative.py       # DeepSeek AI explainer with smart disk caching
│   │   └── batch.py           # Multi-ligand comparative virtual screening
│   └── utils/
│       └── vina_setup.py      # Automated Scripps Vina binary installer
├── bin/
│   └── vina.exe               # Scripps AutoDock Vina binary (Windows x64)
├── data/
│   ├── cache/                 # Local disk cache for API responses & structures
│   └── benchmarks/            # Pre-indexed benchmark case studies
├── frontend/
│   ├── index.html             # Single-page application interface
│   ├── css/
│   │   └── styles.css         # Styling & WebGL canvas viewport setup
│   ├── js/
│   │   ├── app.js             # Main frontend controller & tab state
│   │   ├── viewer.js          # 3Dmol.js WebGL molecular viewer controller
│   │   ├── api.js             # Client API communication wrapper
│   │   └── charts.js          # Radar & ADME visualization with Chart.js
│   └── lib/
│       ├── 3Dmol-min.js       # Bundled 3Dmol.js (works offline)
│       └── chart.min.js       # Bundled Chart.js (works offline)
├── tests/                     # Automated pytest suite
├── .env.example               # Environment variables template
├── .gitignore                 # Strict rules to prevent secret and binary leaks
└── requirements.txt           # Python dependencies
```

---

## 7. Educational & Citation Notice

Bindora is developed for **research, training, and educational use** (B.Pharm, M.Pharm, Pharm.D, and computational pharmacology curricula).

When publishing or citing results generated with Bindora, please cite:
1. **AutoDock Vina:** J. Eberhardt et al., *J. Chem. Inf. Model.* 2021, 61, 8, 3891–3898.
2. **Meeko:** Scripps Research Center for Computational Structural Biology (CCSB).
3. **RDKit:** Open-source cheminformatics toolkit (`https://www.rdkit.org`).
4. **ChEMBL:** European Bioinformatics Institute (EMBL-EBI), *Nucleic Acids Res.* 2019.
