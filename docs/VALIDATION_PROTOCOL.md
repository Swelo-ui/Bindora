# Bindora Dock — Scientific Validation, Thermodynamic Standardization & Benchmarking Protocol

**Version:** 1.0.0 (March 2026)  
**Authors:** Senior Computational Chemistry Software Engineering & Molecular Docking Team  
**Platform Stack:** AutoDock Vina 1.2.5 (Scripps CCSB) | RDKit v2026.03.6 | Meeko | SciPy | PLIP Non-Covalent Engine  

---

## 1. Executive Summary & Foundational Principles

Bindora Dock is an advanced educational and translational molecular docking and drug-discovery platform. It couples AutoDock Vina conformational search with RDKit cheminformatics profiling, automated active-site pocket detection, biophysical interaction fingerprinting, and LLM-synthesized pharmacological narratives.

### Core Scientific Grounding Mandate
To ensure academic publication grade and scientific defensibility:
1. **Never conflate computational predictions with physical experimental measurements.**
2. **Never hard-code benchmark numbers** or artificially force AutoDock Vina 1.2.5 to match literature numbers produced under different empirical scoring functions (e.g., AutoDock 4.2 grid-based electrostatic and desolvation scoring).
3. **Explicitly categorize every reported metric** into one of five mutually distinct data tiers:
   - **Tier A (Predicted):** AutoDock Vina docking score in kcal/mol (empirical scoring function output).
   - **Tier B (Derived / Model-based):** Affinity-Derived Kd-like estimate derived via standard isothermal equilibrium model (T = 298.15 K, RT ≈ 0.592485 kcal/mol). Explicitly marked as a model-derived estimate, not an experimental thermodynamic constant.
   - **Tier C (Calculated):** Deterministic cheminformatics descriptors (RDKit molecular weight, MolLogP, HBD, HBA, TPSA, Rotatable bonds, and normalized Ligand Efficiency).
   - **Tier D (Validation):** Crystallographic redocking heavy-atom coordinate RMSD (≤ 2.0 Å) evaluating protocol pose reproduction.
   - **Tier E (Experimental):** Curated physical wet-lab assay bioactivities (ChEMBL / BindingDB Ki, Kd, IC₅₀).

---

## 2. Standardized Biophysical & Thermodynamic Formulations

### 2.1 AutoDock Vina Docking Score
- **Terminology:** Strictly designated as **"AutoDock Vina Docking Score (kcal/mol)"** or **"Predicted Binding Score"**.
- **Nature of the Metric:** Vina uses an empirical, knowledge-based scoring function parameterized against the PDBbind refined set, combining steric interactions (Gauss 1, Gauss 2), repulsion, hydrophobic contacts, and directional hydrogen bonding.
- **Limitation:** It is an empirical free energy estimator (ΔG_score), not a rigorous path-integral thermodynamic free energy of binding (ΔG°_bind). It omits explicit solvent polarization, finite-temperature receptor conformational entropy, and ion-solvation equilibria.

### 2.2 Affinity-Derived Kd-like Estimate (Model-Derived)
The affinity-derived Kd-like estimate is computed from the standard equilibrium state relation:

```text
ΔG° = R · T · ln(Kd)   ⟹   Kd = exp(ΔG° / (R · T))
```

Where:
- **Temperature (T):** 298.15 K (25.0°C, standard thermodynamic state).
- **Molar Gas Constant (R):** 0.0019872041 kcal/(mol·K).
- **Thermal Energy Product (RT):**
  ```text
  RT = 298.15 × 0.0019872041 ≈ 0.5924849 kcal/mol
  ```
- **Concentration Standard State:** 1.0 M.

#### Conversion Formulas:
```text
Kd [M] = exp(Docking Score / 0.5924849)
Kd [nM] = Kd [M] × 10⁹
Kd [µM] = Kd [M] × 10⁶
pKd = -log₁₀(Kd [M])
```

> **Mandatory Scientific Disclaimer:**
> *"This value is mathematically derived from the docking score and is not an experimentally measured or rigorously calculated thermodynamic Kd."*
> It is reported strictly as a model-derived indicator of energetic magnitude, **never as a substitute for wet-lab in vitro assay constants (Ki / IC₅₀)**.

### 2.3 Ligand Efficiency (LE)
Ligand Efficiency measures binding energy contribution per non-hydrogen (heavy) atom:

```text
LE = |AutoDock Vina Docking Score| / N_heavy
```

- **Units:** kcal/(mol·heavy atom).
- **Lead Discovery Benchmark:** LE ≥ 0.30 kcal/(mol·HA) indicates an efficient ligand scaffold.
- **Size-Independent Ligand Efficiency (SILE):**
  ```text
  SILE = |Score| / (N_heavy)^0.3
  ```
- **Binding Efficiency Index (BEI):**
  ```text
  BEI = pKd / MW [kDa]
  ```

---

## 3. Crystallographic Redocking Self-Validation Protocol

### 3.1 Purpose & Scientific Scope
Crystallographic redocking is the gold standard method to evaluate whether a molecular docking engine and parameter set can accurately predict the crystallographic binding pose of a known ligand inside its native cognate receptor pocket.

### 3.2 Evaluation Threshold
- **Success Criterion:** Heavy-atom coordinate RMSD ≤ 2.0 Å relative to the experimentally determined crystallographic coordinates.
- **Status Classification:**
  - RMSD ≤ 2.0 Å ⟹ **PASS** (Crystallographic Reproduction Validated)
  - RMSD > 2.0 Å ⟹ **FAIL** (Pose Deviation Exceeds Acceptance Threshold)

### 3.3 True In-Place Coordinate RMSD vs Spatial Superposition
1. **No Superposition Permitted:** Docked poses are assessed directly within the receptor binding pocket coordinate frame without translation or rotation. Structural alignment (superposition) is strictly prohibited during redocking validation because it conceals translational/rotational shifts within the pocket.
2. **Symmetry-Aware Graph Isomorphism (Tier 1):**
   Molecules often feature chemically indistinguishable atoms with differing coordinate labels (e.g., 180° flipping of symmetric phenyl rings, carboxylate oxygens, or symmetric branching). Bindora utilizes RDKit graph isomorphism (`GetSubstructMatches(uniquify=False)`) to evaluate all valid symmetry automorphisms without superposition, calculating the true minimum coordinate RMSD:
   ```text
   RMSD_symmetry = min_{m ∈ Automorphisms} √[ (1/N) · Σ ||x_i^docked - x_m(i)^crystal||² ]
   ```
3. **Element-Constrained Optimal Bijective Matching (Tier 2):**
   When format conversion (PDB to PDBQT) alters atom naming or hydrogen assignment, Bindora applies the **Kuhn-Munkres algorithm (Hungarian method)** via `scipy.optimize.linear_sum_assignment` under strict element identity constraints:
   ```text
   C_ij = ||x_i^docked - x_j^crystal||²   if Elem(i) == Elem(j)
        = ∞                              if Elem(i) != Elem(j)
   ```
   This completely eliminates nearest-neighbor double-mapping artifacts and ensures an exact 1-to-1 bijective correspondence.

### 3.4 Disambiguation: Pose vs Rank 1 RMSD vs Crystal Redocking RMSD
- **Pose vs Rank 1 RMSD (`rmsd_lb` / `rmsd_ub`):** Generated by AutoDock Vina's internal clustering algorithm. It quantifies the conformational divergence between alternative modes (modes 2–9) and the top-ranked mode (mode 1). Mode 1 is always 0.000 Å. **This is NOT a redocking validation metric.**
- **Crystallographic Redocking RMSD:** The coordinate distance between the docked mode and the independently resolved crystallographic reference ligand coordinates extracted from the PDB structure.

---

## 4. Docking Experiment Mode Taxonomy

Bindora Dock automatically classifies docking simulations into four distinct experimental regimes:

| Experiment Mode | Receptor Context | Ligand Identity | Validation Applicability | Scientific Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **Native Redocking (Self-Validation)** | Co-crystallized holo structure (e.g. 1M17) | Cognate crystal ligand (e.g. Erlotinib) | **Directly Applicable** (RMSD ≤ 2.0 Å) | Confirms protocol capacity to locate crystallographic binding mode in the cognate pocket. |
| **Cross-Docking / Benchmark Docking** | Holo structure crystallized with a different ligand (e.g. 1M17) | Non-native investigational ligand (e.g. Gefitinib) | **Non-Native Benchmark** (Relative scoring) | Compares relative scoring against literature. Scores may differ from literature due to receptor induced-fit adaptation and scoring function differences (AutoDock 4.2 vs Vina 1.2.5). |
| **Targeted Pocket Docking** | Apo structure or structure without co-crystallized ligand | Investigational compound | **Not Applicable** (No crystal reference) | Docking targeted to computationally predicted cavities or catalytic residues. |
| **Blind Docking** | Full macromolecular surface | Investigational compound | **Exploratory** (Cavity discovery) | Global cavity search across large search volume. |

---

## 5. Non-Covalent Interaction Analysis Standards

Bindora's interaction engine executes rigorous 3D spatial and geometric analysis based on Protein-Ligand Interaction Profiler (PLIP) criteria:

### 5.1 Interaction Criteria Matrix
| Interaction Type | Cutoff Distance | Geometric & Angle Criteria | Participating Elements / Chemical Groups |
| :--- | :--- | :--- | :--- |
| **Hydrogen Bond** | ≤ 3.5 Å | Donor-H···Acceptor ≥ 115.0° | Donor: O, N, S (with polar H); Acceptor: O, N, S |
| **Salt Bridge** | ≤ 4.2 Å | Charged center centroid distance | Cation: Lys NZ, Arg guanidinium; Anion: Asp/Glu carboxylate |
| **π-π Stacking (Parallel)** | ≤ 5.5 Å | Ring normal angle θ ≤ 30.0° | Aromatic rings (Phe, Tyr, Trp, His and ligand aromatics) |
| **π-π Stacking (T-Shaped)** | ≤ 5.5 Å | Ring normal angle 60.0° ≤ θ ≤ 90.0° | Aromatic rings (edge-to-face geometry) |
| **π-Cation** | ≤ 4.5 Å | Cation to aromatic centroid distance | Cation (Lys, Arg) to aromatic ring |
| **Classical Halogen Bond (σ-hole)** | ≤ 3.8 Å | C-X···Acceptor ≥ 130.0° | X ∈ {Cl, Br, I} interacting with Lewis base (O, N, S) |
| **Fluorine Polar Contact** | ≤ 3.5 Å | Non-directional multipolar contact | Organic fluorine (F lacks a classical σ-hole; electrostatically polar) |
| **Hydrophobic Contact** | ≤ 4.0 Å | Inter-atomic distance | Non-polar carbon–carbon contacts |

### 5.2 Physical Distinction: σ-Hole Halogen Bonds vs. Fluorine Contacts
- **Chlorine, Bromine, Iodine:** Possess significant anisotropic electron distribution leading to an electropositive region along the extension of the C-X covalent bond (the σ-hole). These form highly directional halogen bonds (θ ≥ 130°, optimal at 180°).
- **Fluorine:** Due to extremely high electronegativity and low polarizability, fluorine does **not** exhibit a classical electrophilic σ-hole in organic molecules. Its interactions are categorized transparently as **"Fluorine Polar / Multipolar Contacts"** to preserve physical accuracy.

---

## 6. Cheminformatics (ADME & Drug-Likeness) Integrity

#### 6.1 RDKit Descriptor Reliability
All physicochemical descriptors are generated deterministically using RDKit v2026.03.6. The software guarantees numeric integer and float population without missing keys (`"—"`):
- **Molecular Weight (MW):** Exact monoisotopic and average molecular weight (≤ 500 Da for Lipinski compliance).
- **MolLogP:** Wildman-Crippen calculated partition coefficient (≤ 5.0).
- **Hydrogen Bond Donors (HBD):** `rdMolDescriptors.CalcNumLipinskiHBD` (≤ 5).
- **Hydrogen Bond Acceptors (HBA):** `rdMolDescriptors.CalcNumLipinskiHBA` (≤ 10).
- **Topological Polar Surface Area (TPSA):** Calculated polar surface area (≤ 140 Å² for Veber oral bioavailability).
- **Rotatable Bonds (RotB):** Number of non-terminal, non-ring single bonds (≤ 10).

---

## 7. LLM Mechanistic Narrative Anti-Hallucination Guardrails

When synthesizing pharmacological dossier briefings via Gemini or OpenRouter LLM engines:
1. **Zero Number Mutation:** The model is prohibited from modifying, inventing, or hallucinating numerical values, scores, or physical constants.
2. **Explicit Data Separation:**
   - Computational predictions must be identified as such.
   - Derived values (e.g. theoretical Kd) must state their thermodynamic origin.
   - Literature benchmarks from ChEMBL must be clearly labeled as independent wet-lab assays.
3. **Prohibition of In Vivo Claims:** The model is strictly barred from asserting that in silico docking poses establish in vivo clinical efficacy, pharmacodynamic safety, or therapeutic potency.
4. **Deterministic Fallback Engine:** If API access is offline, the platform seamlessly deploys an offline, rule-based pharmacology reasoning engine that produces 100% grounded dossiers.

---

## 8. Verification Benchmark Case Studies (Genuine Empirical Outputs)

### Benchmark 1: Erlotinib Native Redocking into EGFR Kinase (PDB: 1M17)
- **Target Receptor:** Epidermal Growth Factor Receptor (EGFR) Kinase Domain (PDB ID: `1M17`, Chain A).
- **Native Crystallographic Ligand:** Erlotinib (AQ4, N_heavy = 29, Rotatable Bonds = 11).
- **Search Space (Bounding Box):** Center = (22.01, 0.25, 52.79), Dimensions = 22.0 × 22.0 × 22.0 Å.
- **Docking Engine Configuration:** AutoDock Vina 1.2.5, Exhaustiveness = 8, Seed = 42, Modes = 9.
- **Docking Execution Time:** 51.61 s.
- **AutoDock Vina Docking Score (Rank 1):** **-7.10 kcal/mol**.
- **Crystallographic Redocking RMSD (Rank 1):** **1.52 Å** (≤ 2.0 Å ⟹ **PASS**).
- **Atom Mapping Engine:** `topological_symmetry_graph_isomorphism` (Exact chemical graph isomorphism via template bond-order assignment).
- **Affinity-Derived Kd-like Estimate:** 6289.2 nM (6.29 µM) *(Model-derived estimate; not an experimental thermodynamic Kd)*.
- **Ligand Efficiency (LE):** 0.245 kcal/(mol·HA) (|-7.10| / 29).
- **Pose Distribution Across Modes:**
  - Mode 1: -7.10 kcal/mol, RMSD = 1.52 Å (Native active pose)
  - Mode 2: -6.98 kcal/mol, RMSD = 2.40 Å
  - Mode 3: -6.97 kcal/mol, RMSD = 8.70 Å
- **Validation Outcome:** **PASS** (Crystallographic Reproduction Validated).
- **ChEMBL Comparison:** Experimental wet-lab biochemical Ki for Erlotinib against human EGFR is ≈ 2.1 nM (ChEMBL assay records). The difference from the in silico derived estimate (6.29 µM) reflects the empirical nature of scoring functions and reinforces the necessity of explicit terminology separation.

### Benchmark 2: Indinavir Native Redocking into HIV-1 Protease C2 Homodimer (PDB: 1HSG)
- **Target Receptor:** HIV-1 Protease Homodimer (PDB ID: `1HSG`, preserving both catalytic chains A & B).
- **Native Crystallographic Ligand:** Indinavir (MK1, N_heavy = 45, Rotatable Bonds = 13).
- **Search Space (Bounding Box):** Center = (13.07, 22.47, 5.56), Dimensions = 22.0 × 22.0 × 22.0 Å.
- **Docking Engine Configuration:** AutoDock Vina 1.2.5, Exhaustiveness = 8, Seed = 42, Modes = 9.
- **Docking Execution Time:** 257.57 s across 4 parallel CPU worker threads.
- **AutoDock Vina Docking Score (Rank 1):** **-10.54 kcal/mol**.
- **Crystallographic Redocking RMSD (Rank 1):** **0.46 Å** (≤ 2.0 Å ⟹ **PASS — Sub-Angstrom Accuracy**).
- **Atom Mapping Engine:** `topological_symmetry_graph_isomorphism` evaluating 12 chemical symmetry automorphisms.
- **Affinity-Derived Kd-like Estimate:** 19.0 nM (0.019 µM) *(Model-derived estimate; not an experimental thermodynamic Kd)*.
- **Ligand Efficiency (LE):** 0.234 kcal/(mol·HA) (|-10.54| / 45).
- **Pose Distribution Across Modes:**
  - Mode 1: -10.54 kcal/mol, RMSD = 0.46 Å (Near-perfect crystallographic reproduction)
  - Mode 2: -10.24 kcal/mol, RMSD = 10.54 Å (Inverted homodimer binding orientation)
  - Mode 3: -10.12 kcal/mol, RMSD = 10.67 Å
- **Validation Outcome:** **PASS** (Gold-standard crystallographic pose reproduced with sub-Angstrom precision).

### Benchmark 3: Gefitinib Cross-Docking into EGFR Kinase (PDB: 1M17)
- **Target Receptor:** EGFR Kinase Domain (PDB ID: `1M17`, crystallized with native Erlotinib).
- **Docked Ligand:** Gefitinib (Iressa, N_heavy = 31).
- **Experiment Mode:** **Cross-Docking / Benchmark Docking** (Pocket crystallized around Erlotinib, not Gefitinib).
- **AutoDock Vina 1.2.5 Docking Score:** -7.99 kcal/mol.
- **Affinity-Derived Kd-like Estimate:** 1390 nM (1.39 µM).
- **Ligand Efficiency (LE):** 0.258 kcal/(mol·HA) (|-7.99| / 31).
- **Literature Comparison Context:** AutoDock 4.2 literature benchmarks report -9.0 to -10.5 kcal/mol using AutoDock 4.2 semi-empirical force fields with Amber-based electrostatic charges. AutoDock Vina 1.2.5 uses an independent empirical scoring function; reporting -7.99 kcal/mol is scientifically honest and must never be altered to fake a match with literature.
