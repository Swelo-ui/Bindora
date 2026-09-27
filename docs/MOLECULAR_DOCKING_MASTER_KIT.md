# Bindora Dock — Complete Molecular Docking Master Handbook & Scientific Kit
**Author:** Bindora Computational Pharmacology Studio  
**Engine:** Scripps AutoDock Vina v1.2.5 & RDKit Open-Source Cheminformatics  
**Target Audience:** Students, Researchers, Medicinal Chemists, & Molecular Modelers  

---

## Table of Contents
1. [Introduction: Molecular Docking Kya Hai aur Kyu Zaroori Hai?](#1-introduction-molecular-docking-kya-hai)
2. [Complete Dictionary of Abbreviations & Scientific Terms (Har Word Ka Matlab)](#2-complete-dictionary-of-abbreviations--scientific-terms)
3. [Biological & 3D Structural Anatomy (Structure Ko Kaise Samajhein?)](#3-biological--3d-structural-anatomy)
4. [Non-Covalent Chemical Interactions (3D Bonds Ki Pehchan)](#4-non-covalent-chemical-interactions)
5. [The Science Behind the "Cloud" (Molecular Surface & Pocket Cavity)](#5-the-science-behind-the-cloud-molecular-surface)
6. [Bindora Dock: Step-by-Step Practical Docking Workflow](#6-bindora-dock-step-by-step-practical-docking-workflow)
7. [Grid Box, Exhaustiveness, & Native CPU Multi-Threading](#7-grid-box-exhaustiveness--native-cpu-multi-threading)
8. [Flexible Side Chains (Induced-Fit Simulation)](#8-flexible-side-chains-induced-fit-simulation)
9. [Redocking Self-Validation & RMSD Benchmark](#9-redocking-self-validation--rmsd-benchmark)
10. [ADME, BOILED-Egg, & Pharmacokinetics Interpretation](#10-adme-boiled-egg--pharmacokinetics-interpretation)
11. [Troubleshooting & Pro Tips for Publishing Research Papers](#11-troubleshooting--pro-tips-for-publishing-research-papers)

---

## 1. Introduction: Molecular Docking Kya Hai?

### 1.1 Simple Real-World Intuition
Sochiye aapke paas ek **Taala (Lock)** hai — ye hamara **Target Protein / Receptor** (jaise koi bimari failane wala enzyme ya virus ka protein) hai.  
Aapke paas hajaron alag-alag aakaar ki **Chabiyan (Keys)** hain — ye hamare **Drug Molecules / Ligands** hain.

**Molecular Docking** computer par chalne wali wo simulation hai jo do sawalon ka jawab deti hai:
1. **Pose Prediction:** Kya chabi taale ke keyhole (binding pocket) ke andar ghus sakti hai, aur ghusne ke baad kis angle par baithti hai?
2. **Affinity Estimation:** Chabi taale ke kitne tight fit baithti hai? (Binding Free Energy, $\Delta G$).

```
[Target Protein / Lock] + [Drug Candidate / Key] 
              ↓ (AutoDock Vina Global Search)
    [Protein-Ligand Complex in Optimal 3D Conformation]
```

### 1.2 Two Classical Models of Binding
1. **Fischer's Lock-and-Key Model (1894):** Protein ka pocket rigid (sakht) hota hai aur ligand aakar exact fit hota hai. (Standard rigid docking is concept par chalti hai).
2. **Koshland's Induced-Fit Model (1958):** Jaise haath dastane (glove) me jata hai to dastana haath ke aakaar ke anusaar thoda sa adjust hota hai, waise hi jab drug protein pocket me ghusti hai to pocket ke amino acids thode move karte hain. (Bindora ka **Flexible Side Chains** feature is Induced-Fit ko simulate karta hai).

---

## 2. Complete Dictionary of Abbreviations & Scientific Terms

Yeh dictionary docking me aane wale har technical word aur shortform ka complete encyclopedia hai:

| Term / Shortform | Full Form | Scientific Meaning (Asal Matlab) | Good vs Bad Value |
| :--- | :--- | :--- | :--- |
| **$\Delta G$ (Delta G)** | Gibbs Free Energy of Binding | Ligand aur protein judne par kitni energy release hoti hai ($\text{kcal/mol}$). Jitni zyada negative, utna strong bond. | **Good:** $\le -7.0\text{ kcal/mol}$<br>**Weak:** $> -5.0\text{ kcal/mol}$ |
| **$K_d$ / $K_i$** | Dissociation / Inhibition Constant | Dawa aur protein ke complex ko todne ke liye kitni concentration chahiye ($\text{nM}$ ya $\mu\text{M}$). | **Good:** $< 100\text{ nM}$ (Tight)<br>**Weak:** $> 10\text{ }\mu\text{M}$ |
| **RMSD** | Root Mean Square Deviation | Do 3D structures ke beech ka physical distance difference ($\text{\AA}$). Redocking validation me use hota hai. | **Success:** $< 2.0\text{ \AA}$<br>**Fail:** $> 2.0\text{ \AA}$ |
| **PDB** | Protein Data Bank | Protein ke har atom ke $X, Y, Z$ 3D coordinates store karne wala standard format. | File extension: `.pdb` |
| **PDBQT** | Protein Data Bank + Charges (Q) + Torsions (T) | AutoDock Vina ka required format. Isme har atom par **Partial Charge ($Q$)** aur ghumnay wale bonds (**Torsion Tree, $T$**) hote hain. | File extension: `.pdbqt` |
| **SMILES** | Simplified Molecular Input Line Entry System | Kisi chemical molecule ka 1-line text representation (e.g., Ethanol = `CCO`, Aspirin = `CC(=O)Oc1ccccc1C(=O)O`). | Standard chemical string |
| **LE** | Ligand Efficiency | Binding energy divided by heavy atom count ($-\Delta G / N_{\text{heavy}}$). Batata hai ki molecule ka har atom kitna productive hai. | **Good:** $\ge 0.30\text{ kcal/mol/atom}$ |
| **LipE / LLE** | Lipophilic Ligand Efficiency | $\text{pIC}_{50} - \text{cLogP}$. Batata hai ki binding sirf grease (fat) ki wajah se to nahi aa rahi, balki specific chemical bonds se aa rahi hai. | **Target:** $\ge 5.0$ |
| **cLogP / WLOGP** | Octanol-Water Partition Coefficient | Molecule kitna oily/fat-soluble (lipophilic) hai vs water-soluble (hydrophilic). | **Ideal:** $1.0 \text{ to } 3.0$<br>**Poor:** $> 5.0$ (Too greasy) |
| **TPSA** | Topological Polar Surface Area | Molecule ke polar atoms (Oxygen, Nitrogen, Polar H) ka total surface area ($\text{\AA}^2$). | **Oral Drug:** $\le 140\text{ \AA}^2$<br>**Brain (BBB):** $\le 90\text{ \AA}^2$ |
| **HBD** | Hydrogen Bond Donors | Wo Hydrogen atoms jo electronegative N ya O se jude hote hain (e.g., $-\text{OH}, -\text{NH}_2$). | **Rule of 5:** $\le 5$ |
| **HBA** | Hydrogen Bond Acceptors | Wo Nitrogen ya Oxygen jinke paas lone pair hota hai jo H ko attract karta hai. | **Rule of 5:** $\le 10$ |
| **RotB** | Rotatable Bonds | Molecule ke single non-ring bonds jo freely 360 degree ghum sakte hain. Jyada hone par molecule loose/floppy ho jata hai. | **Ideal:** $\le 10$ |
| **MW** | Molecular Weight | Molecule ka molecular mass ($\text{g/mol}$). | **Lipinski:** $\le 500\text{ Da}$ |
| **PAINS** | Pan-Assay Interference Compounds | Aise chemical groups jo lab testing me fake/false positive results dete hain (e.g., rhodanines, catechols). | **Result:** Zero alerts required |
| **QED** | Quantitative Estimate of Drug-likeness | 0 se 1 ke scale par dawa jese gun hone ka overall composite index. | **Good:** $> 0.67$ |
| **SA Score** | Synthetic Accessibility Score | 1 (chemistry lab me banana bohot aasan) se 10 (banana lagbhag namumkin) tak ka score. | **Lead candidate:** $< 4.0$ |
| **CYP450** | Cytochrome P450 Enzymes | Liver ke 5 mukhya enzymes (CYP1A2, 2C9, 2C19, 2D6, 3A4) jo dawa ko metabolize karke body se bahar nikalte hain. | Inhibition check zaroori hai |
| **hERG** | Human Ether-à-go-go-Related Gene | Dil (heart) ka potassium channel. Agar dawa ise block kare to cardiac arrhythmia (heart attack) ho sakta hai. | Non-inhibitor hona chahiye |
| **Ames Mutagenicity** | Ames Bacterial Reverse Mutation Test | Batata hai ki kya molecule DNA me mutation ya Cancer create kar sakta hai. | Negative hona zaroori hai |

---

## 3. Biological & 3D Structural Anatomy

Protein aur Ligand ko 3D viewer me pehchanne ka tarika:

```
               PROTEIN ARCHITECTURE
┌────────────────────────────────────────────────────────┐
│  Primary:    Amino Acid Sequence (Met-Ala-His-Leu...)   │
│  Secondary:  Alpha-Helices (Springs) & Beta-Sheets     │
│  Tertiary:   Folded 3D Globular Structure              │
│  Quaternary: Multi-subunit Complex (e.g. Hemoglobin)   │
└────────────────────────────────────────────────────────┘
```

### 3.1 Amino Acids Ke 4 Parivaar (The 20 Building Blocks)
Protein 20 amino acids se banta hai. Inke gun yaad rakhna docking samajhne ki chabi hai:

1. **Hydrophobic / Greasy (Non-Polar):**
   - **Members:** `LEU` (Leucine), `ILE` (Isoleucine), `VAL` (Valine), `PHE` (Phenylalanine), `MET` (Methionine), `ALA` (Alanine), `PRO` (Proline), `TRP` (Tryptophan).
   - **Role:** Ye pocket ke andar paani se door chupkar baithte hain. Drug ke aromatic rings aur carbon chains ko pakadte hain.
2. **Positively Charged / Basic:**
   - **Members:** `LYS` (Lysine, $+\text{NH}_3^+$), `ARG` (Arginine, Guanidinium group), `HIS` (Histidine, partial positive).
   - **Role:** Drug ke negative groups ($-\text{COO}^-$, $-\text{SO}_3^-$) ke sath majboot **Salt Bridges** banate hain.
3. **Negatively Charged / Acidic:**
   - **Members:** `ASP` (Aspartate, $-\text{COO}^-$), `GLU` (Glutamate, $-\text{COO}^-$).
   - **Role:** Drug ke basic Nitrogen ($-\text{NH}_2, -\text{NH}^+$) ko khinchkar Ionic interaction karte hain.
4. **Polar Neutral (Hydrogen Bonders):**
   - **Members:** `SER` (Serine, $-\text{OH}$), `THR` (Threonine, $-\text{OH}$), `TYR` (Tyrosine, Phenolic $-\text{OH}$), `ASN` (Asparagine), `GLN` (Glutamine), `CYS` (Cysteine, $-\text{SH}$).
   - **Role:** Directional Hydrogen Bonds bana kar drug ki orientation fix karte hain.

### 3.2 3D Element Color Code (CPK Standard)
Viewer me har atom ka rang international standard par based hota hai:
- **Carbon (C):** Cyan (`cyanCarbon`) ya Green (`greenCarbon`) ya Gray (`whiteCarbon`).
- **Oxygen (O):** Red (Lal).
- **Nitrogen (N):** Dark Blue (Neela).
- **Hydrogen (H):** White (Safed).
- **Sulfur (S):** Yellow (Peela).
- **Phosphorus (P):** Orange (Narangi).
- **Fluorine (F):** Light Green.
- **Chlorine (Cl):** Bright Green.
- **Bromine (Br):** Dark Red / Brown.
- **Iodine (I):** Purple.

---

## 4. Non-Covalent Chemical Interactions (3D Bonds Ki Pehchan)

Dawa protein se Fevicol ki tarah permanent nahi judti; wo **kamjor lekin hajaron non-covalent forces** ke sum se judti hai. Bindora inko 3D viewer aur 2D LigPlot schematic me live draw karta hai:

```
┌────────────────────────────────────────────────────────────────────────┐
│                    NON-COVALENT FORCES HIERARCHY                       │
│                                                                        │
│  1. Salt Bridge (Ionic)   :  ~3.0 - 5.0 kcal/mol  [Strongest]          │
│  2. Hydrogen Bond         :  ~1.0 - 3.0 kcal/mol  [Directional]        │
│  3. Pi-Cation / Pi-Pi     :  ~1.0 - 2.5 kcal/mol  [Aromatic rings]     │
│  4. Halogen Bond          :  ~0.5 - 2.0 kcal/mol  [Sigma-hole]         │
│  5. Hydrophobic / VDW     :  ~0.5 - 1.0 kcal/mol  [Bulk surface pack]  │
└────────────────────────────────────────────────────────────────────────┘
```

1. **Hydrogen Bonds (H-Bonds — Peeli/Yellow Dashed Lines):**
   - **Distance:** $1.8\text{ \AA} \text{ se } 3.2\text{ \AA}$.
   - **Rule:** Ek Donor ($-\text{N}-\text{H}$ ya $-\text{O}-\text{H}$) apna proton Acceptor ($=\text{O}$ ya $:\text{N}$) ko share karta hai. Angle $120^\circ \text{ se } 180^\circ$ hona chahiye.
2. **Salt Bridges (Ionic Bonds — Gulabi/Magenta Dashed Lines):**
   - **Distance:** $\le 4.0\text{ \AA}$.
   - **Rule:** Full formal negative charge (`ASP`/`GLU`) aur full formal positive charge (`ARG`/`LYS`) ke beech ka electrostatic attraction.
3. **$\pi$-$\pi$ Stacking (Aromatic Interaction — Cyan/Aasmani Dashed Lines):**
   - **Distance:** $\le 4.5\text{ \AA}$.
   - **Rule:** Do aromatic benzene rings (jaise drug ki ring aur protein ke `PHE`, `TYR`, ya `TRP`) ek dusre ke upar parallel (Face-to-Face) ya T-shape (Edge-to-Face) baithti hain.
4. **$\pi$-Cation Interaction (Narangi/Orange Lines):**
   - **Distance:** $\le 4.5\text{ \AA}$.
   - **Rule:** Aromatic electron cloud aur positive charge (`LYS`, `ARG`, ya tertiary amine) ka attraction.
5. **Halogen Bonds (Green Lines):**
   - **Distance:** $\le 3.5\text{ \AA}$.
   - **Rule:** Halogen atom (Cl, Br, I) ke aage ek electropositive "$\sigma$-hole" hota hai jo Oxygen lone pair se judta hai.

---

## 5. The Science Behind the "Cloud" (Molecular Surface)

Bindora ke 3D viewer me **"Pocket Surface (5Å)"** ya **"Surface"** toggle on karne par jo translucent grey/colored cloud dikhta hai, uska exact scientific formula yeh hai:

$$\text{Surface Selection} = \left\{ \text{Atom}_i \in \text{Receptor} \mid \min_{j \in \text{Ligand}} \|\mathbf{r}_i - \mathbf{r}_j\| \le 5.5\text{ \AA} \right\}$$

### Kyu Ye Kisi Par Banta Hai aur Kisi Par Adha/Nahi Banta?
1. **Gufa (Buried Pocket) vs. Darwaza (Solvent Front):**
   - **Pichla Hissa (Cloud Covered):** Drug ka pichla hissa active site ke andar gehraai me hota hai. Waha protein ke amino acids use charo taraf se gherte hain, isliye waha protein ki deewaar (cloud) banti hai.
   - **Aage Ka Hissa (Khula / No Cloud):** Dawa ka aage ka hissa (jaise Gefitinib ka morpholine group) cellular cytoplasm/paani ki taraf bahar jhankta hai taaki dawa blood me dissolve ho sake. Darwaze par protein atoms maujood nahi hote, isliye waha cloud physically ban hi nahi sakta!
2. **Cloud Colors Ka Scientific Code:**
   - **White/Grey:** Hydrophobic non-polar pocket walls (`LEU`, `VAL`, `PHE`).
   - **Lal (Red):** Electronegative Oxygen (`ASP`, `GLU`).
   - **Neela (Blue):** Electropositive Nitrogen (`LYS`, `ARG`).
   - **Peela/Green:** Gatekeeper Methionine/Cysteine Sulfur (`MET`, `CYS`).

---

## 6. Bindora Dock: Step-by-Step Practical Docking Workflow

Docking run karne ka standard scientific protocol:

```
[Tab 1: Target Prep]  ──→  [Tab 2: Molecular Docking]  ──→  [Tab 3: ADME Profile]
   • RCSB PDB Download         • SMILES to 3D Conformer        • Lipinski Rule of 5
   • Water Strip               • Grid Box Definition           • BOILED-Egg Brain/Gut
   • Gasteiger Charges         • Vina Multi-Core Search        • Bioavailability Radar
                               • Redocking Self-Validation
                                        ↓
[Tab 7: Dossier / Audit] ←── [Tab 5: AI Report] ←── [Tab 4: ChEMBL Cross-Check]
   • Scripps Vina SHA-256      • Publication Paper Synth       • Experimental IC50 / Ki
   • Reproducibility Seeds     • Lead Optimization Advice      • Wet-Lab Assays Match
```

### Step 1: Target Protein Prepare Karna (Tab 1)
1. RCSB PDB ID enter karein (e.g., `1CX2` for COX-2, `2ITY` for EGFR, `1IEP` for BCR-ABL Kinase, `1HSG` for HIV Protease).
2. Bindora automatic clean-up run karta hai:
   - Bulk crystallographic waters (`HOH`) ko remove karta hai taaki pocket khali ho sake.
   - Crystallization buffer ions (`SO4`, `GOL`, `CL`, `NA`) ko filter karta hai.
   - PDBQT format me convert karke polar Hydrogens aur **Gasteiger Partial Charges** assign karta hai.

### Step 2: Ligand Prepare Karna (Tab 2)
1. Ligand input karein: SMILES paste karein ya 2D Molecule Ketcher me draw karein.
2. Bindora internally **RDKit ETKDGv3** algorithm se 2D representation ko lowest-energy 3D conformation me convert karta hai aur **MMFF94 force-field** se geometry minimize karta hai.
3. Meeko engine se ligand ki **Torsion Tree** define hoti hai (konsa bond ghumega, konsa rigid ring hai).

### Step 3: Grid Box Set Karna
- **Pocket Centroid Button:** Automatically co-crystallized drug ke coordinates ka geometric center $(X, Y, Z)$ calculate karta hai.
- **Blind Docking Button:** Grid box ko $60\text{ \AA} \times 60\text{ \AA} \times 60\text{ \AA}$ tak expand karta hai agar aapko pocket ki location pata na ho.
- **Ideal Box Size:** Usually $20\text{ \AA} \times 20\text{ \AA} \times 20\text{ \AA}$ se $24\text{ \AA} \times 24\text{ \AA} \times 24\text{ \AA}$ sabse accurate results deta hai.

### Step 4: Execute Docking
- Click **"Execute 3D Molecular Docking"**. AutoDock Vina local search chalayega aur best 9 binding poses calculate karega.

---

## 7. Grid Box, Exhaustiveness, & Native CPU Multi-Threading

### 7.1 Exhaustiveness Kya Hai?
AutoDock Vina ka global search algorithm **Iterated Local Search (ILS) with Monte Carlo** par chalta hai.  
**Exhaustiveness** ka matlab hai: *Algorithm ne energy landscape par kitni independent search trajectories chalayi hain.*

| Exhaustiveness Level | Independent Searches | CPU Time | Accuracy & Use Case |
| :--- | :--- | :--- | :--- |
| **4 (Screening Fast)** | 4 Trajectories | ~2 to 5 seconds | Huge compound libraries screening (rough filtering). |
| **8 (Standard Academic)** | 8 Trajectories | ~8 to 15 seconds | College research, preliminary binding pose check. |
| **16 (Publication Grade)** | 16 Trajectories | ~25 to 45 seconds | Peer-reviewed journals, high-confidence binding modes. |
| **32 (Deep Exploration)** | 32 Trajectories | ~1 to 2 minutes | Highly flexible molecules ($RotB > 8$) & benchmark redocking. |

### 7.2 Native CPU Hardware Acceleration
Bindora web-browser me chalne ke bawajood client PC ke **physical & logical CPU cores** ko real-time utilize karta hai:
- Agar aapke laptop me **4 Cores / 4 Threads** hain (e.g. Intel Core i5), to Vina `--cpu 4` command ke sath 4 parallel C++ threads dispatch karta hai.
- 8-Core PC par Vina `--cpu 8` par chalta hai, jisse calculation 8 guna fast ho jaati hai.
- Header me **"Native CPU"** badge par click karke aap har ek core ka live utilization load monitor kar sakte hain.

---

## 8. Flexible Side Chains (Induced-Fit Simulation)

### 8.1 Rigid Receptor Ki Kami
Standard docking me protein ke amino acids pathar ki deewaar ki tarah freeze hote hain. Real biology me jab drug pocket me ghusti hai, to active site ke amino acids apna rasta badalte hain. Agar koi amino acid rasta block kar raha ho, to standard docking me drug waha ghus hi nahi paati (steric clash penalty).

### 8.2 Bindora Flexible Residue Engine
Bindora Meeko polymer partitioning use karta hai:
1. Target active site ke flexible candidates (e.g., `A:TYR:456`, `A:MET:766`) select karein.
2. Engine protein structure ko do hisson me tod deta hai:
   - **Rigid Receptor PDBQT:** Baaki poora protein static rehta hai.
   - **Flex PDBQT:** Selected amino acids ke single bonds freely rotate ho sakte hain.
3. Vina docking ke dauran ligand aur selected side-chains **dono ko ek sath move karta hai**. Isse true **Induced-Fit Docking** prapt hoti hai!

---

## 9. Redocking Self-Validation & RMSD Benchmark

Science me kisi computational prediction ko tab tak valid nahi maana jata jab tak wo **known crystallographic benchmark** ko reproduce na kar sake.

```
       REDOCKING VALIDATION WORKFLOW
       
       Crystallographic Co-crystal (RCSB PDB)
                     │
       ┌─────────────┴─────────────┐
       ▼                           ▼
[Extract Native Ligand]     [Clean Receptor PDBQT]
       │                           │
       └─────────────┬─────────────┘
                     ▼
          [Dock Native Ligand Back]
                     │
                     ▼
         [Calculate RMSD Alignment]
    ┌────────────────┬────────────────┐
    ▼                                 ▼
RMSD ≤ 2.0 Å                      RMSD > 2.0 Å
[VALIDATED: PASS]                 [INVALID: FAIL]
```

### 9.1 RMSD Ka Formula
Root Mean Square Deviation crystal pose aur redocked pose ke har heavy atom ke beech ke distance ka quadratic mean hota hai:

$$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \|\mathbf{r}_i^{\text{docked}} - \mathbf{r}_i^{\text{crystal}}\|^2}$$

### 9.2 Validation Criteria
- **$\text{RMSD} \le 1.0\text{ \AA}$ (Sub-Angstrom):** Exceptional crystallographic match (Gold Standard).
- **$1.0\text{ \AA} < \text{RMSD} \le 2.0\text{ \AA}$:** Validated research benchmark. Paper me publish karne yogya.
- **$\text{RMSD} > 2.0\text{ \AA}$:** Failed validation. Iska matlab grid box galat hai ya exhaustiveness kam hai.

---

## 10. ADME, BOILED-Egg, & Pharmacokinetics Interpretation

Dawa chahe pocket me kitni bhi majbooti se dock ho jaye (chahe $\Delta G = -15\text{ kcal/mol}$ ho), agar wo pet me pachegi nahi ya liver me toxic ban jayegi, to wo mariz ko theek nahi kar sakti. Isliye **ADME Profiling** zaroori hai:

```
            BOILED-Egg PHARMACOKINETIC MODEL
┌─────────────────────────────────────────────────────────┐
│  WLOGP (Lipophilicity)                                  │
│   ▲                                                     │
│   │       ┌────────────────────────┐                    │
│   │       │   WHITE (HIA)          │                    │
│   │       │   High Gastrointestinal│                    │
│   │       │   Absorption           │                    │
│   │       │       ┌─────────────┐  │                    │
│   │       │       │ YELLOW YOLK │  │                    │
│   │       │       │ Blood-Brain │  │                    │
│   │       │       │ Barrier BBB │  │                    │
│   │       │       └─────────────┘  │                    │
│   │       └────────────────────────┘                    │
│   │                                                     │
│   │   GREY: Poor absorption / Outside Druggable Space   │
│   └───────────────────────────────────────────────►     │
│                                           TPSA (Polar)  │
└─────────────────────────────────────────────────────────┘
```

1. **BOILED-Egg White Zone (HIA):** Human Intestinal Absorption. Dawa goli (tablet) ke roop me khane par aanto me aasani se sokh li jayegi.
2. **BOILED-Egg Yellow Yolk (BBB):** Blood-Brain Barrier penetration.
   - **Depression / Alzheimer / Epilepsy Drugs:** Yolk ke andar honi chahiye.
   - **Heart / Diabetes / Cancer Drugs:** Yolk se bahar honi chahiye taaki dimaagi side-effects (dizziness, sedation) na hon.
3. **P-glycoprotein Substrate ($PGP^+$):** Blue dot ka matlab P-gp pump dawa ko brain se bahar phenk dega.
4. **Lipinski's Rule of 5:**
   - Molecular Weight $\le 500\text{ Da}$
   - $\text{cLogP} \le 5$
   - $\text{HBD} \le 5$
   - $\text{HBA} \le 10$
   - *Max 1 violation allowed for oral drugs.*

---

## 11. Troubleshooting & Pro Tips for Publishing Research Papers

1. **Binding Energy Bohot Kam ($>-5.0\text{ kcal/mol}$) Kyu Aa Rahi Hai?**
   - Molecule bohot chota (fragment) hai. Heavy atoms badhayein ya aromatic rings jodkar hydrophobic contacts banayein.
2. **Redocking Fail Kyu Ho Rahi Hai?**
   - Grid box size badhayein ($22\text{ \AA}$ karein).
   - Exhaustiveness ko 8 se badhakar 32 karein.
   - Pocket Centroid coordinate check karein.
3. **Research Paper Me Kaise Report Karein?**
   - Software citation: *AutoDock Vina 1.2.5 (Scripps Research Institute)*, *RDKit (Open-source cheminformatics)*, *Meeko 0.5.x*.
   - Binding Affinity ko hamesha **$\text{kcal/mol}$** me likhein aur saath me **$K_i$ calculate** karein.
   - 2D Interaction Diagram (Hydrogen bonds, Salt bridges, RMSD) ko figure bana kar insert karein.

---
*Bindora Dock — Precision Preclinical Pharmacology & Molecular Simulation System.*
