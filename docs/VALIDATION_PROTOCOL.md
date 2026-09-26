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
   - **Tier B (Derived / Model-based):** Affinity-Derived $K_d$-like estimate derived via standard isothermal equilibrium model ($T = 298.15\text{ K}$, $RT \approx 0.592485\text{ kcal/mol}$). Explicitly marked as a model-derived estimate, not an experimental thermodynamic constant.
   - **Tier C (Calculated):** Deterministic cheminformatics descriptors (RDKit molecular weight, MolLogP, HBD, HBA, TPSA, Rotatable bonds, and normalized Ligand Efficiency).
   - **Tier D (Validation):** Crystallographic redocking heavy-atom coordinate RMSD ($\le 2.0\text{ \AA}$) evaluating protocol pose reproduction.
   - **Tier E (Experimental):** Curated physical wet-lab assay bioactivities (ChEMBL / BindingDB $K_i$, $K_d$, $\text{IC}_{50}$).

---

## 2. Standardized Biophysical & Thermodynamic Formulations

### 2.1 AutoDock Vina Docking Score
- **Terminology:** Strictly designated as **"AutoDock Vina Docking Score (kcal/mol)"** or **"Predicted Binding Score"**.
- **Nature of the Metric:** Vina uses an empirical, knowledge-based scoring function parameterized against the PDBbind refined set, combining steric interactions (Gauss 1, Gauss 2), repulsion, hydrophobic contacts, and directional hydrogen bonding.
- **Limitation:** It is an empirical free energy estimator ($\Delta G_\text{score}$), not a rigorous path-integral thermodynamic free energy of binding ($\Delta G^\circ_\text{bind}$). It omits explicit solvent polarization, finite-temperature receptor conformational entropy, and ion-solvation equilibria.

### 2.2 Affinity-Derived $K_d$-like Estimate (Model-Derived)
The affinity-derived $K_d$-like estimate is computed from the standard equilibrium state relation:
$$\Delta G^\circ = R \cdot T \cdot \ln(K_d) \quad \implies \quad K_d = \exp\left(\frac{\Delta G^\circ}{R \cdot T}\right)$$

Where:
- **Temperature ($T$):** $298.15\text{ K}$ ($25.0^\circ\text{C}$, standard thermodynamic state).
- **Molar Gas Constant ($R$):** $0.0019872041\text{ kcal}/(\text{mol}\cdot\text{K})$.
- **Thermal Energy Product ($RT$):**
  $$RT = 298.15 \times 0.0019872041 \approx 0.5924849\text{ kcal/mol}$$
- **Concentration Standard State:** $1.0\text{ M}$.

#### Conversion Formulas:
$$K_d\text{ [M]} = \exp\left(\frac{\text{Docking Score}}{0.5924849}\right)$$
$$K_d\text{ [nM]} = K_d\text{ [M]} \times 10^9$$
$$K_d\text{ [}\mu\text{M]} = K_d\text{ [M]} \times 10^6$$
$$pK_d = -\log_{10}(K_d\text{ [M]})$$

> **Mandatory Scientific Disclaimer:**
> *"This value is mathematically derived from the docking score and is not an experimentally measured or rigorously calculated thermodynamic $K_d$."*
> It is reported strictly as a model-derived indicator of energetic magnitude, **never as a substitute for wet-lab in vitro assay constants ($K_i / \text{IC}_{50}$)**.

### 2.3 Ligand Efficiency (LE)
Ligand Efficiency measures binding energy contribution per non-hydrogen (heavy) atom:
$$\text{LE} = \frac{|\text{AutoDock Vina Docking Score}|}{N_\text{heavy}}$$

- **Units:** $\text{kcal}/(\text{mol}\cdot\text{heavy atom})$.
- **Lead Discovery Benchmark:** $\text{LE} \ge 0.30\text{ kcal}/(\text{mol}\cdot\text{HA})$ indicates an efficient ligand scaffold.
- **Size-Independent Ligand Efficiency (SILE):**
  $$\text{SILE} = \frac{|\text{Score}|}{(N_\text{heavy})^{0.3}}$$
- **Binding Efficiency Index (BEI):**
  $$\text{BEI} = \frac{pK_d}{\text{MW [kDa]}}$$

---

## 3. Crystallographic Redocking Self-Validation Protocol

### 3.1 Purpose & Scientific Scope
Crystallographic redocking is the gold standard method to evaluate whether a molecular docking engine and parameter set can accurately predict the crystallographic binding pose of a known ligand inside its native cognate receptor pocket.

### 3.2 Evaluation Threshold
- **Success Criterion:** Heavy-atom coordinate $\text{RMSD} \le 2.0\text{ \AA}$ relative to the experimentally determined crystallographic coordinates.
- **Status Classification:**
  - $\text{RMSD} \le 2.0\text{ \AA} \implies$ **PASS** (Crystallographic Reproduction Validated)
  - $\text{RMSD} > 2.0\text{ \AA} \implies$ **FAIL** (Pose Deviation Exceeds Acceptance Threshold)

### 3.3 True In-Place Coordinate RMSD vs Spatial Superposition
1. **No Superposition Permitted:** Docked poses are assessed directly within the receptor binding pocket coordinate frame without translation or rotation. Structural alignment (superposition) is strictly prohibited during redocking validation because it conceals translational/rotational shifts within the pocket.
2. **Symmetry-Aware Graph Isomorphism (Tier 1):**
   Molecules often feature chemically indistinguishable atoms with differing coordinate labels (e.g., $180^\circ$ flipping of symmetric phenyl rings, carboxylate oxygens, or symmetric branching). Bindora utilizes RDKit graph isomorphism (`GetSubstructMatches(uniquify=False)`) to evaluate all valid symmetry automorphisms without superposition, calculating the true minimum coordinate RMSD:
   $$\text{RMSD}_\text{symmetry} = \min_{m \in \text{Automorphisms}} \sqrt{\frac{1}{N}\sum_{i=1}^{N} \|\mathbf{x}_i^\text{docked} - \mathbf{x}_{m(i)}^\text{crystal}\|^2}$$
3. **Element-Constrained Optimal Bijective Matching (Tier 2):**
   When format conversion (PDB to PDBQT) alters atom naming or hydrogen assignment, Bindora applies the **Kuhn-Munkres algorithm (Hungarian method)** via `scipy.optimize.linear_sum_assignment` under strict element identity constraints:
   $$C_{ij} = \begin{cases} \|\mathbf{x}_i^\text{docked} - \mathbf{x}_j^\text{crystal}\|^2 & \text{if } \text{Elem}(i) = \text{Elem}(j) \\ \infty & \text{if } \text{Elem}(i) \ne \text{Elem}(j) \end{cases}$$
   This completely eliminates nearest-neighbor double-mapping artifacts and ensures an exact 1-to-1 bijective correspondence.

### 3.4 Disambiguation: Pose vs Rank 1 RMSD vs Crystal Redocking RMSD
- **Pose vs Rank 1 RMSD (`rmsd_lb` / `rmsd_ub`):** Generated by AutoDock Vina's internal clustering algorithm. It quantifies the conformational divergence between alternative modes (modes 2–9) and the top-ranked mode (mode 1). Mode 1 is always $0.000\text{ \AA}$. **This is NOT a redocking validation metric.**
- **Crystallographic Redocking RMSD:** The coordinate distance between the docked mode and the independently resolved crystallographic reference ligand coordinates extracted from the PDB structure.

---

## 4. Docking Experiment Mode Taxonomy

Bindora Dock automatically classifies docking simulations into four distinct experimental regimes:

| Experiment Mode | Receptor Context | Ligand Identity | Validation Applicability | Scientific Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **Native Redocking (Self-Validation)** | Co-crystallized holo structure (e.g. 1M17) | Cognate crystal ligand (e.g. Erlotinib) | **Directly Applicable** ($\text{RMSD} \le 2.0\text{ \AA}$) | Confirms protocol capacity to locate crystallographic binding mode in the cognate pocket. |
| **Cross-Docking / Benchmark Docking** | Holo structure crystallized with a different ligand (e.g. 1M17) | Non-native investigational ligand (e.g. Gefitinib) | **Non-Native Benchmark** (Relative scoring) | Compares relative scoring against literature. Scores may differ from literature due to receptor induced-fit adaptation and scoring function differences (AutoDock 4.2 vs Vina 1.2.5). |
| **Targeted Pocket Docking** | Apo structure or structure without co-crystallized ligand | Investigational compound | **Not Applicable** (No crystal reference) | Docking targeted to computationally predicted cavities or catalytic residues. |
| **Blind Docking** | Full macromolecular surface | Investigational compound | **Exploratory** (Cavity discovery) | Global cavity search across large search volume. |

---

## 5. Non-Covalent Interaction Analysis Standards

Bindora's interaction engine executes rigorous 3D spatial and geometric analysis based on Protein-Ligand Interaction Profiler (PLIP) criteria:

### 5.1 Interaction Criteria Matrix
| Interaction Type | Cutoff Distance | Geometric & Angle Criteria | Participating Elements / Chemical Groups |
| :--- | :--- | :--- | :--- |
| **Hydrogen Bond** | $\le 3.5\text{ \AA}$ | $\text{Donor-H}\cdots\text{Acceptor} \ge 115.0^\circ$ | Donor: O, N, S (with polar H); Acceptor: O, N, S |
| **Salt Bridge** | $\le 4.2\text{ \AA}$ | Charged center centroid distance | Cation: Lys $\text{NZ}$, Arg guanidinium; Anion: Asp/Glu carboxylate |
| **$\pi$-$\pi$ Stacking (Parallel)** | $\le 5.5\text{ \AA}$ | Ring normal angle $\theta \le 30.0^\circ$ | Aromatic rings (Phe, Tyr, Trp, His and ligand aromatics) |
| **$\pi$-$\pi$ Stacking (T-Shaped)** | $\le 5.5\text{ \AA}$ | Ring normal angle $60.0^\circ \le \theta \le 90.0^\circ$ | Aromatic rings (edge-to-face geometry) |
| **$\pi$-Cation** | $\le 4.5\text{ \AA}$ | Cation to aromatic centroid distance | Cation (Lys, Arg) to aromatic ring |
| **Classical Halogen Bond ($\sigma$-hole)** | $\le 3.8\text{ \AA}$ | $\text{C-X}\cdots\text{Acceptor} \ge 130.0^\circ$ | $\text{X} \in \{\text{Cl}, \text{Br}, \text{I}\}$ interacting with Lewis base (O, N, S) |
| **Fluorine Polar Contact** | $\le 3.5\text{ \AA}$ | Non-directional multipolar contact | Organic fluorine (F lacks a classical $\sigma$-hole; electrostatically polar) |
| **Hydrophobic Contact** | $\le 4.0\text{ \AA}$ | Inter-atomic distance | Non-polar carbon–carbon contacts |

### 5.2 Physical Distinction: $\sigma$-Hole Halogen Bonds vs. Fluorine Contacts
- **Chlorine, Bromine, Iodine:** Possess significant anisotropic electron distribution leading to an electropositive region along the extension of the $\text{C-X}$ covalent bond (the $\sigma$-hole). These form highly directional halogen bonds ($\theta \ge 130^\circ$, optimal at $180^\circ$).
- **Fluorine:** Due to extremely high electronegativity and low polarizability, fluorine does **not** exhibit a classical electrophilic $\sigma$-hole in organic molecules. Its interactions are categorized transparently as **"Fluorine Polar / Multipolar Contacts"** to preserve physical accuracy.

---

## 6. Cheminformatics (ADME & Drug-Likeness) Integrity

### 6.1 RDKit Descriptor Reliability
All physicochemical descriptors are generated deterministically using RDKit v2026.03.6. The software guarantees numeric integer and float population without missing keys (`"—"`):
- **Molecular Weight (MW):** Exact monoisotopic and average molecular weight ($\le 500\text{ Da}$ for Lipinski compliance).
- **MolLogP:** Wildman-Crippen calculated partition coefficient ($\le 5.0$).
- **Hydrogen Bond Donors (HBD):** `rdMolDescriptors.CalcNumLipinskiHBD` ($\le 5$).
- **Hydrogen Bond Acceptors (HBA):** `rdMolDescriptors.CalcNumLipinskiHBA` ($\le 10$).
- **Topological Polar Surface Area (TPSA):** Calculated polar surface area ($\le 140\text{ \AA}^2$ for Veber oral bioavailability).
- **Rotatable Bonds (RotB):** Number of non-terminal, non-ring single bonds ($\le 10$).

---

## 7. LLM Mechanistic Narrative Anti-Hallucination Guardrails

When synthesizing pharmacological dossier briefings via Gemini or OpenRouter LLM engines:
1. **Zero Number Mutation:** The model is prohibited from modifying, inventing, or hallucinating numerical values, scores, or physical constants.
2. **Explicit Data Separation:**
   - Computational predictions must be identified as such.
   - Derived values (e.g. theoretical $K_d$) must state their thermodynamic origin.
   - Literature benchmarks from ChEMBL must be clearly labeled as independent wet-lab assays.
3. **Prohibition of In Vivo Claims:** The model is strictly barred from asserting that in silico docking poses establish in vivo clinical efficacy, pharmacodynamic safety, or therapeutic potency.
4. **Deterministic Fallback Engine:** If API access is offline, the platform seamlessly deploys an offline, rule-based pharmacology reasoning engine that produces 100% grounded dossiers.

---

## 8. Verification Benchmark Case Studies (Genuine Empirical Outputs)

### Benchmark 1: Erlotinib Native Redocking into EGFR Kinase (PDB: 1M17)
- **Target Receptor:** Epidermal Growth Factor Receptor (EGFR) Kinase Domain (PDB ID: `1M17`, Chain A).
- **Native Crystallographic Ligand:** Erlotinib (AQ4, $N_\text{heavy} = 29$, Rotatable Bonds = 11).
- **Search Space (Bounding Box):** Center = $(22.01, 0.25, 52.79)$, Dimensions = $22.0 \times 22.0 \times 22.0\text{ \AA}$.
- **Docking Engine Configuration:** AutoDock Vina 1.2.5, Exhaustiveness = 8, Seed = 42, Modes = 9.
- **Docking Execution Time:** $51.61\text{ s}$.
- **AutoDock Vina Docking Score (Rank 1):** **$-7.10\text{ kcal/mol}$**.
- **Crystallographic Redocking RMSD (Rank 1):** **$1.52\text{ \AA}$** ($\le 2.0\text{ \AA} \implies$ **PASS**).
- **Atom Mapping Engine:** `topological_symmetry_graph_isomorphism` (Exact chemical graph isomorphism via template bond-order assignment).
- **Affinity-Derived $K_d$-like Estimate:** $6289.2\text{ nM}$ ($6.29\text{ }\mu\text{M}$) *(Model-derived estimate; not an experimental thermodynamic $K_d$)*.
- **Ligand Efficiency (LE):** $0.245\text{ kcal}/(\text{mol}\cdot\text{HA})$ ($|-7.10| / 29$).
- **Pose Distribution Across Modes:**
  - Mode 1: $-7.10\text{ kcal/mol}$, $\text{RMSD} = 1.52\text{ \AA}$ (Native active pose)
  - Mode 2: $-6.98\text{ kcal/mol}$, $\text{RMSD} = 2.40\text{ \AA}$
  - Mode 3: $-6.97\text{ kcal/mol}$, $\text{RMSD} = 8.70\text{ \AA}$
- **Validation Outcome:** **PASS** (Crystallographic Reproduction Validated).
- **ChEMBL Comparison:** Experimental wet-lab biochemical $K_i$ for Erlotinib against human EGFR is $\approx 2.1\text{ nM}$ (ChEMBL assay records). The difference from the in silico derived estimate ($6.29\text{ }\mu\text{M}$) reflects the empirical nature of scoring functions and reinforces the necessity of explicit terminology separation.

### Benchmark 2: Indinavir Native Redocking into HIV-1 Protease C2 Homodimer (PDB: 1HSG)
- **Target Receptor:** HIV-1 Protease Homodimer (PDB ID: `1HSG`, preserving both catalytic chains A & B).
- **Native Crystallographic Ligand:** Indinavir (MK1, $N_\text{heavy} = 45$, Rotatable Bonds = 13).
- **Search Space (Bounding Box):** Center = $(13.07, 22.47, 5.56)$, Dimensions = $22.0 \times 22.0 \times 22.0\text{ \AA}$.
- **Docking Engine Configuration:** AutoDock Vina 1.2.5, Exhaustiveness = 8, Seed = 42, Modes = 9.
- **Docking Execution Time:** $257.57\text{ s}$ across 4 parallel CPU worker threads.
- **AutoDock Vina Docking Score (Rank 1):** **$-10.54\text{ kcal/mol}$**.
- **Crystallographic Redocking RMSD (Rank 1):** **$0.46\text{ \AA}$** ($\le 2.0\text{ \AA} \implies$ **PASS — Sub-Angstrom Accuracy**).
- **Atom Mapping Engine:** `topological_symmetry_graph_isomorphism` evaluating 12 chemical symmetry automorphisms.
- **Affinity-Derived $K_d$-like Estimate:** $19.0\text{ nM}$ ($0.019\text{ }\mu\text{M}$) *(Model-derived estimate; not an experimental thermodynamic $K_d$)*.
- **Ligand Efficiency (LE):** $0.234\text{ kcal}/(\text{mol}\cdot\text{HA})$ ($|-10.54| / 45$).
- **Pose Distribution Across Modes:**
  - Mode 1: $-10.54\text{ kcal/mol}$, $\text{RMSD} = 0.46\text{ \AA}$ (Near-perfect crystallographic reproduction)
  - Mode 2: $-10.24\text{ kcal/mol}$, $\text{RMSD} = 10.54\text{ \AA}$ (Inverted homodimer binding orientation)
  - Mode 3: $-10.12\text{ kcal/mol}$, $\text{RMSD} = 10.67\text{ \AA}$
- **Validation Outcome:** **PASS** (Gold-standard crystallographic pose reproduced with sub-Angstrom precision).

### Benchmark 3: Gefitinib Cross-Docking into EGFR Kinase (PDB: 1M17)
- **Target Receptor:** EGFR Kinase Domain (PDB ID: `1M17`, crystallized with native Erlotinib).
- **Docked Ligand:** Gefitinib (Iressa, $N_\text{heavy} = 31$).
- **Experiment Mode:** **Cross-Docking / Benchmark Docking** (Pocket crystallized around Erlotinib, not Gefitinib).
- **AutoDock Vina 1.2.5 Docking Score:** $-7.99\text{ kcal/mol}$.
- **Affinity-Derived $K_d$-like Estimate:** $1390\text{ nM}$ ($1.39\text{ }\mu\text{M}$).
- **Ligand Efficiency (LE):** $0.258\text{ kcal}/(\text{mol}\cdot\text{HA})$ ($|-7.99| / 31$).
- **Literature Comparison Context:** AutoDock 4.2 literature benchmarks report $-9.0\text{ to }-10.5\text{ kcal/mol}$ using AutoDock 4.2 semi-empirical force fields with Amber-based electrostatic charges. AutoDock Vina 1.2.5 uses an independent empirical scoring function; reporting $-7.99\text{ kcal/mol}$ is scientifically honest and must never be altered to fake a match with literature.
