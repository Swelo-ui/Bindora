# Bindora Dock — 50-Run Research-Grade QA Protocol Audit Report
**Execution Target:** Live Flask Backend (`http://127.0.0.1:5000`)  
**Timestamp:** 2026-09-29  
**Platform:** Windows 10 (AMD64), Python 3.13.1, AutoDock Vina 1.2.5, RDKit 2026.03.6, OpenMM 8.x (OpenCL / Reference / CPU), Gemmi 0.7.5, Meeko 0.5.x  
**External Binaries:** AutoDock Vina (`H:\AnuDock\bin\vina.exe` — available); GNINA (not installed / unavailable); fpocket (not installed / unavailable)  
**Standard Execution Parameters:** Exhaustiveness = 8 (strictly maintained across all docking evaluations), Default Seed = 42  
**Consolidated Raw Logs Directory:** `data/benchmarks/qa_50run_raw_logs/` (50 JSON run records preserved)

---

## Executive Summary

A comprehensive 50-run observational quality assurance audit was conducted live against Bindora Dock running on localhost. Across 16 distinct feature categories, the platform demonstrated strong performance in native crystallographic redocking (5/5 PASS, all RMSD ≤ 1.50 Å), pocket auto-detection (100% active-site centering), interaction fingerprint detection (Asp189 salt bridge in trypsin, Met769/793 hinge H-bond in EGFR, 9-H-bond network in streptavidin), and ChEMBL bioactivity cross-referencing.

However, the audit revealed several critical anomalies, including:
1. An unhandled `HTTP 500` crash across all Induced-Fit Docking runs due to an internal keyword argument mismatch (`ligand_sdf_or_pdbqt` vs `ligand_smiles_or_pdbqt`).
2. A failure of the PAINS filter to detect curcumin (a canonical PAINS prototype), while correctly flagging rhodanines.
3. Overprediction of blood-brain barrier permeability for loperamide due to the absence of a P-glycoprotein active efflux transporter model in the BOILED-Egg implementation.
4. False-positive `HIGH_CONFIDENCE` consensus scoring of a volatile alkane decoy (pentane) over a nanomolar nanomolar drug (erlotinib).
5. Lack of determinism in `/api/docking/run` due to omission of the `seed` parameter when calling `DockingEngine.run_docking`, causing inter-run pose drift up to 1.30 Å.
6. Severe timeouts (>450–720s) on large macrocycles (>80 heavy atoms: Cyclosporine A, Vancomycin) under standard exhaustiveness.
7. AI Narrative factual misattributions where metrics from low-ranked decoy poses were erroneously attributed to the top-ranked mode.

In accordance with the audit protocol rules, **zero files in the Bindora codebase were modified, patched, or improved**. All findings below represent pure observational empirical data.

---

## 1. Feature Category Results Tables (Runs 1 to 50)

### Category 1: Native Redocking Validation (5 runs)
*Protocol: Redock native co-crystallized ligand into its crystallographic pocket at exhaustiveness = 8, seed = 42.*

| Run # | Target PDB | Native Ligand | Observed Vina Score | Rank 1 RMSD (Å) | Best Mode RMSD (Å) | Benchmark Status | Literature Expectation | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **1** | 1HSG | Indinavir (MK1) | −10.90 kcal/mol | 0.45 Å | 0.45 Å (Mode 1) | PASS | Sub-2.0 Å (often < 1.5 Å) | **MATCH** | Exceptional near-zero RMSD; 173.4s runtime reflects 14 rotatable bonds |
| **2** | 1STP | Biotin (BTN) | −6.97 kcal/mol | 0.62 Å | 0.62 Å (Mode 1) | PASS | Well under 1.0 Å | **MATCH** | Reconstructed sub-Ångström pose; tight pocket lock |
| **3** | 3PTB | Benzamidine (BEN) | −6.08 kcal/mol | 0.39 Å | 0.37 Å (Mode 2) | PASS | Sub-1.0 Å; rigid ligand | **MATCH** | Mode 1 RMSD = 0.39 Å; salt bridge with Asp189 recovered |
| **4** | 1M17 | Erlotinib (AQ4) | −7.13 kcal/mol | 1.50 Å | 1.50 Å (Mode 1) | PASS | RMSD ≤ 2.0 Å; hinge binder | **MATCH** | Hinge H-bond to Met769/793 formed; Mode 1 satisfies standard criterion |
| **5** | 1IEP | Imatinib (STI) | −12.65 kcal/mol | 0.34 Å | 0.34 Å (Mode 1) | PASS | Rigid redocking known hard case; large RMSD expected in cross-conformation | **MATCH (SURPASS)** | Native redocking into crystallized DFG-out conformation achieved 0.34 Å |

---

### Category 2: Cross-Docking Validation (3 runs)
*Protocol: Dock non-native ligands into unrelated/different receptors at exhaustiveness = 8, seed = 42.*

| Run # | Investigational Compound | Target Receptor | Observed Vina Score | Ligand Efficiency (LE) | Consensus Rank | Literature Expectation | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **6** | Aspirin | 1HSG (HIV-1 Protease) | −5.55 kcal/mol | 0.277 kcal/mol/atom | 1 | Weak, non-specific binding (> −6.0 kcal/mol) | **MATCH** | Score above −6.0 kcal/mol threshold; flagged as weak binder; LE suboptimal (< 0.30) |
| **7** | Atorvastatin | 3PTB (Trypsin) | −7.07 kcal/mol | 0.354 kcal/mol/atom | 1 | Poor steric fit into small S1 protease pocket; low LE | **PARTIAL MATCH** | Score is −7.07 kcal/mol but statin is forced onto outer cleft; runtime 153.5s due to 10 rotatable bonds |
| **8** | Erlotinib | 1IEP (Abl Kinase) | −7.02 kcal/mol | 0.351 kcal/mol/atom | 1 | Moderate affinity plausible (both ATP kinases), but clearly distinguishable from native | **MATCH** | Score (−7.02 kcal/mol) is markedly weaker than native Imatinib (−12.65 kcal/mol in Run 5) |

---

### Category 3: MM-GBSA Physics Consistency (3 runs)
*Protocol: EGFR Kinase (1M17) targeted with three ligands of increasing documented potency at exhaustiveness = 8, seed = 42.*

| Run # | Compound | Potency Profile | Observed Vina Score | MM-GBSA ΔG (docked) | MM-GBSA ΔG (refined) | Expected Directional Trend | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **9** | Aniline | Non-binder / fragment | −4.59 kcal/mol | −11.54 kcal/mol | −11.54 kcal/mol | A clearly worse than B & C | **MATCH** | Both Vina and MM-GBSA score fragment clearly lower than active kinase inhibitors |
| **10** | Gefitinib | Nanomolar EGFR binder | −7.82 kcal/mol | −18.13 kcal/mol | −18.13 kcal/mol | High affinity; similar to C | **MATCH** | Physics-based MM-GBSA ΔG (−18.13 kcal/mol) reflects strong van der Waals & desolvation burial |
| **11** | Erlotinib | Nanomolar EGFR binder | −7.17 kcal/mol | −16.14 kcal/mol | −16.14 kcal/mol | High affinity; similar to B | **MATCH** | Scores in identical range to Gefitinib (within 0.65 kcal Vina, 1.99 kcal MM-GBSA) |

---

### Category 4: Ligand Strain Differentiation (3 runs)
*Protocol: Evaluate intramolecular strain on docked 3D coordinates via `/api/docking/refine`.*

| Run # | Compound | Target | Conformational Flexibility | Initial Strain | Minimized Strain | Strain Relaxation | Literature Expectation | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **12** | Benzamidine | 3PTB | Rigid amidine fragment (0 active rotatable bonds) | 19.40 kcal/mol | 19.40 kcal/mol | **0.00 kcal/mol** | Low/negligible strain relaxation | **MATCH** | Zero relaxation energy; confirmed rigid conformation |
| **13** | Atorvastatin | 3PTB | Highly flexible statin (10 rotatable bonds) | 66.61 kcal/mol | 62.17 kcal/mol | **4.44 kcal/mol** | Moderate/variable strain depending on docked pose | **MATCH** | Notable 4.44 kcal/mol relaxation upon local vacuum minimization |
| **14** | Cyclosporine A | 1M17 | Large cyclic undecapeptide (85 heavy atoms) | N/A | N/A | N/A | Stress test: measure stability or failure | **TIMEOUT / CRASH** | Request timed out after **722.92 seconds** during exhaustiveness=8 docking |

---

### Category 5: Consensus Scoring Concordance (3 runs)
*Protocol: Target 1M17; examine multi-engine consensus matrix, rank aggregation, and CNN status.*

| Run # | Compound | Role / Nature | Observed Vina Score | Consensus Score | Consensus Confidence | GNINA / CNN Status | Literature Expectation | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **15** | Erlotinib | True clinical binder | −7.17 kcal/mol | −0.786 | `MODERATE_CONFIDENCE` | Unavailable (`bin/` missing) | Expect `HIGH_CONFIDENCE` | **MISMATCH** | Scored only MODERATE due to mode rank dispersion across internal terms |
| **16** | Pentane | Simple alkane decoy | −2.89 kcal/mol | −1.320 | `HIGH_CONFIDENCE` | Unavailable (`bin/` missing) | Expect `DECOY` or low confidence | **CRITICAL MISMATCH** | Decoy assigned `HIGH_CONFIDENCE` because single mode had zero rank spread across available engines |
| **17** | Aspirin | Borderline non-binder | −5.71 kcal/mol | −1.619 | `MODERATE_CONFIDENCE` | Unavailable (`bin/` missing) | Marginal / Low confidence | **BORDERLINE** | Assigned MODERATE; GNINA binary absent, surrogate formula in effect |

---

### Category 6: Covalent Warhead Detection (3 runs)
*Protocol: Evaluate electrophilic warheads and active-site nucleophile geometry via `/api/covalent/evaluate`.*

| Run # | Ligand | Receptor PDB | Warhead Chemotype | Detected Warhead Type | Target Nucleophile | Reactive Distance | Literature Mechanism | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **18** | Ibrutinib | 5P9J (BTK) | Acrylamide | `Michael Acceptor (Acrylamide / Enone)` | CYS 481:A (SG) | 4.14 Å | Conjugate 1,4-Addition to Cys481 | **MATCH** | Warhead recognized; Cys481 SG located at 4.14 Å (`PERMISSIVE_COVALENT_PROXIMITY`, bonus −1.47 kcal) |
| **19** | Nirmatrelvir | 7RFS (Mpro) | Nitrile | `Electrophilic Nitrile / Cyanamide` | CYS 145:A (SG) | 7.73 Å | Reversible nucleophilic addition to Cys145 | **MATCH (Detection)** | Warhead recognized; catalytic dyad Cys145 and His41 identified; rigid docked pose had Cys145 at 7.73 Å |
| **20** | Afatinib | 4I22 (EGFR) | Acrylamide | `Michael Acceptor (Acrylamide / Enone)` | CYS 797:A (SG) | 4.01 Å | Conjugate 1,4-Addition to Cys797 | **MATCH** | Cys797 SG located at 4.01 Å, attack angle 127.5°, feasibility score 0.629 (`PERMISSIVE_PROXIMITY`) |

---

### Category 7: Induced-Fit Docking (3 runs)
*Protocol: Execute Monte Carlo loop and backbone phi/psi induced-fit docking via `/api/docking/induced-fit`.*

| Run # | Compound | Receptor PDB | Target Feature | HTTP Status | Response Error Message | Literature Expectation | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **21** | Imatinib | 1IEP (Abl) | DFG-out loop plastic rearrangement | **500** | `unexpected keyword argument 'ligand_sdf_or_pdbqt'` | Large backbone adjustment | **SYSTEM ERROR** | Backend crash: `backend/app.py` passes `ligand_sdf_or_pdbqt` but service expects `ligand_smiles_or_pdbqt` |
| **22** | Dasatinib | 1IEP (Abl) | Alternative Abl inhibitor chemotype | **500** | `unexpected keyword argument 'ligand_sdf_or_pdbqt'` | Moderate loop adjustment | **SYSTEM ERROR** | Same backend parameter mismatch crash |
| **23** | Erlotinib | 1M17 (EGFR) | Subtle hinge-region induced fit | **500** | `unexpected keyword argument 'ligand_sdf_or_pdbqt'` | Smaller backbone adjustment than Abl | **SYSTEM ERROR** | Same backend parameter mismatch crash |

---

### Category 8: ADME / Drug-Likeness Plausibility (3 runs)
*Protocol: Evaluate pharmacokinetic predictions via `/api/pkpd/adme` against established pharmacology.*

| Run # | Compound | MW (Da) | WLogP | TPSA (Å²) | Predicted GI Absorption | Predicted BBB Permeation | Known Clinical Pharmacokinetics | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **24** | Aspirin | 180.16 | 1.31 | 63.60 | High | Permeant | High GI absorption; crosses BBB | **MATCH** | Correctly predicted High GI and BBB permeation; 0 Rule-of-5 violations |
| **25** | Atorvastatin | 558.65 | 6.31 | 111.79 | Moderate / Low | Non-permeant | Poor BBB penetration; high hepatic uptake | **MATCH** | Accurately predicted non-permeant; flags MW > 500 and LogP > 5.0 |
| **26** | Loperamide | 477.05 | 5.09 | 43.78 | High | **Permeant (Likely crosses BBB)** | Peripheral opioid; **NO CNS penetration** due to P-gp efflux | **CRITICAL MISMATCH** | BOILED-Egg fails because active P-gp efflux is omitted; predicts CNS permeation based on low TPSA |

---

### Category 9: PAINS / Brenk Alert Detection (3 runs)
*Protocol: Substructure alert scanning via `/api/pkpd/adme` against Baell & Holloway (2010) and Brenk (2008).*

| Run # | Compound | Test Profile | Observed PAINS Alerts | Observed Brenk Alerts | Literature Expectation | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **27** | Curcumin | Canonical PAINS benchmark | **0 alerts (Clear)** | 2 alerts (`beta-keto/anhydride`, `Michael_acceptor_1`) | **MUST trigger PAINS** (Baell & Holloway 2010) | **CRITICAL MISMATCH** | PAINS filter completely missed curcumin; Brenk filter caught 2 structural alerts |
| **28** | Aspirin | Clean NSAID drug | **0 alerts (Clear)** | 1 alert (`phenol_ester`) | Clean drug; 0 PAINS alerts | **MATCH** | Clean on PAINS; Brenk ester alert accurately reflects hydrolyzable acetyl ester |
| **29** | 5-benzylidenerhodanine | Classic Rhodanine PAINS | **1 alert (`rhod_sat_A(33)`)** | 1 alert (`Thiocarbonyl_group`) | **MUST trigger PAINS** | **MATCH** | Correctly flagged with specific rhodanine SMARTS catalog match |

---

### Category 10: Interaction Fingerprint Accuracy (3 runs)
*Protocol: Contact analysis via `/api/docking/interactions` on crystallographic / native redocked complexes.*

| Run # | Complex | Target Active Site | H-Bonds Detected | Salt Bridges Detected | Key Literature Contact Checked | Observed Distance | Match / Mismatch | Notes & Anomalies |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **30** | 3PTB / Benzamidine | Trypsin S1 pocket | 3 | 1 | **Asp189 Salt Bridge** | 3.26 Å | **MATCH** | Asp189:A carboxylate to amidine cation salt bridge captured at 3.26 Å |
| **31** | 1M17 / Erlotinib | EGFR ATP cleft | 2 | 0 | **Met769/793 Hinge H-Bond** | 3.23 Å | **MATCH** | Backbone amide N of Met769 (Met793 canonical) to quinazoline N captured at 3.23 Å, 155.7° |
| **32** | 1STP / Biotin | Streptavidin tetramer | 9 | 1 | **Extensive H-Bond Network** | Multiple (2.36–3.44 Å) | **MATCH** | Dense network: Ser88, Asn49 (x2), Ser45 (x2), Trp92, Asn23 (x2), Ser27; Asp128 salt bridge |

---

### Category 11: Batch / Virtual Screening Pipeline (3 runs)
*Protocol: Screen 5 compounds (Gefitinib, Erlotinib, Aspirin, Benzamidine, Aniline) against 1M17 across 3 independent iterations at exhaustiveness = 8.*

| Run # | Iteration | Runtime | Rank 1 | Rank 2 | Rank 3 | Rank 4 | Rank 5 | Ranking Stability | Score Drift Across Runs |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **33** | Batch 1 | 92.50s | Gefitinib (−7.36) | Erlotinib (−7.23) | Aspirin (−5.71) | Benzamidine (−5.31) | Aniline (−4.60) | **100% Concordant** | Sub-0.01 kcal for rigid (Aspirin, Benzamidine) |
| **34** | Batch 2 | 91.29s | Gefitinib (−7.81) | Erlotinib (−7.36) | Aspirin (−5.71) | Benzamidine (−5.30) | Aniline (−4.55) | **100% Concordant** | ~0.45 kcal for flexible Gefitinib (unseeded) |
| **35** | Batch 3 | 90.08s | Gefitinib (−7.85) | Erlotinib (−7.15) | Aspirin (−5.71) | Benzamidine (−5.31) | Aniline (−4.54) | **100% Concordant** | Ranking order strictly preserved (1→5) |

---

### Category 12: ChEMBL Bioactivity Cross-Validation (3 runs)
*Protocol: Automated wet-lab assay retrieval via `/api/pkpd/crosscheck`.*

| Run # | Drug Searched | Target Searched | Molecule ChEMBL ID | Target ChEMBL ID | Status Badge | Experimental Records | Top Wet-Lab Measurement | Match / Mismatch |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **36** | Imatinib | ABL1 | CHEMBL941 | CHEMBL3099 | Experimentally Corroborated | 3 records | IC50 = 10.8 nM (Abl radiometric kinase assay) | **MATCH** |
| **37** | Aspirin | PTGS1 | CHEMBL25 | CHEMBL221 | ChEMBL Records Available | 5 records | IC50 = 11.4 µM (COX1 radioimmunoassay) | **MATCH** |
| **38** | ZINC999999NonexistentMol | ABL1 | None | CHEMBL3099 | No ChEMBL Records | 0 records | Graceful "no data" reported; zero hallucinations | **MATCH** |

---

### Category 13: Pocket Auto-Detection Accuracy (3 runs)
*Protocol: Active-site cavity auto-detection via `/api/structure/receptor` compared against co-crystallized native ligand centroid.*

| Run # | Target PDB | Target Biology | Auto-Detected Center (x, y, z) | Native Ligand Centroid (x, y, z) | Centroid Offset (Å) | Detected Pocket Description | Match / Mismatch |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **39** | 1HSG | HIV-1 Protease | (13.07, 22.47, 5.56) | (13.07, 22.47, 5.56) | **0.00 Å** | Auto-centered on MK1 native pocket (Chain B) | **MATCH** |
| **40** | 3PTB | Trypsin S1 Pocket | (−1.76, 14.46, 16.92) | (−1.76, 14.46, 16.92) | **0.00 Å** | Auto-centered on BEN native pocket (Chain A) | **MATCH** |
| **41** | 1M17 | EGFR ATP Cleft | (22.01, 0.25, 52.79) | (22.01, 0.25, 52.79) | **0.00 Å** | Auto-centered on AQ4 native pocket (Chain A) | **MATCH** |

---

### Category 14: Macrocycle Handling (3 runs)
*Protocol: Conformational ensemble generation via `/api/macrocycle/sample` (srETKDGv3 / MMFF94) and docking at exhaustiveness = 8.*

| Run # | Compound | Macrocycle Profile | Ring Size | Heavy Atoms | Conformers Sampled | Docking Affinity | Runtime | Observation & Verdict |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **42** | Cyclosporine A | Cyclic undecapeptide | 33 | 85 | N/A | N/A | 458.37s | **TIMEOUT**: AutoDock Vina Monte Carlo search cannot converge at exhaustiveness=8 |
| **43** | Lorlatinib | Synthetic kinase macrocycle | 14 | 30 | 4 retained | −8.86 kcal/mol | 50.62s | **SUCCESS**: srETKDGv3 samples low-strain ring; docks cleanly with high affinity |
| **44** | Vancomycin | Glycopeptide antibiotic | 16 | 99 | N/A | N/A | 528.61s | **TIMEOUT**: Massively exceeds practical small-molecule docking limits (MW 1449 Da) |

---

### Category 15: Reproducibility / Determinism (3 runs)
*Protocol: Erlotinib docked against EGFR (1M17) with fixed seed = 42, exhaustiveness = 8, back-to-back under identical runtime conditions.*

| Run # | Iteration | Vina Docking Score | MM-GBSA ΔG | Consensus Rank | Consensus Confidence | RMSD to Run 45 (Å) | Bit-Identical? | Drift Severity |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **45** | Run 1 | −7.129 kcal/mol | −15.44 kcal/mol | 1 | `MODERATE_CONFIDENCE` | Baseline (0.00 Å) | No | Baseline |
| **46** | Run 2 | −7.230 kcal/mol | −17.75 kcal/mol | 1 | `HIGH_CONFIDENCE` | **1.30 Å** | No | Moderate drift (score Δ = 0.10 kcal; RMSD = 1.30 Å; tier changed) |
| **47** | Run 3 | −7.286 kcal/mol | −17.64 kcal/mol | 1 | `HIGH_CONFIDENCE` | **0.79 Å** | No | Moderate drift (score Δ = 0.16 kcal; RMSD = 0.79 Å; tier changed) |

---

### Category 16: AI Narrative Self-Consistency Check (3 runs)
*Protocol: Cross-check AI Pharmacologist narrative text against raw JSON data generated by backend.*

| Run # | Test Case | Source Run ID | Stated Narrative Claim | Underlying Raw JSON Value | Consistency Status | Discrepancy Description |
|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **48** | Erlotinib / 1M17 | Run 11 | "pose carries a high intramolecular strain of 15.38 kcal/mol... consensus scoring is DISCORDANT" | Mode 1 strain = **8.08 kcal/mol**; Mode 1 confidence = `MODERATE_CONFIDENCE` | **FACTUAL MISATTRIBUTION** | LLM pulled 15.38 kcal/mol from **Mode 9** (`poses[8]`, `DECOY_HIGH_STRAIN`) and attributed it to top pose Mode 1 |
| **48** | Erlotinib / 1M17 | Run 11 | "MM-GBSA ΔG: −9.42 kcal/mol" | Mode 1 MM-GBSA = **−16.14 kcal/mol** | **FACTUAL MISATTRIBUTION** | LLM pulled −9.42 kcal/mol from **Mode 8** (`poses[7]`) instead of Mode 1 |
| **49** | Aspirin / 1HSG | Run 6 | "RMSD relative to the native reference (9.52/10.99 Å)" | Native reference comparison was **Not Computed** | **METRIC MISINTERPRETATION** | LLM misread Vina's internal pose-dispersion values (`rmsd_lb`=9.52 Å, `rmsd_ub`=10.99 Å from Mode 1) as native crystal RMSD |
| **50** | Pentane / 1M17 | Run 16 | "Heavy atoms: 20... Ligand efficiency (LE): 0.145" | Pentane has **5 heavy atoms**; server defaulted `heavy_atoms: 20` | **PROPAGATED ERROR** | Narrative propagated backend default heavy atom count of 20, calculating false LE of 0.145 instead of 0.579 |

---

## 2. Anomalies Found

During the 50-run live protocol, seven distinct anomalies were discovered:

### Anomaly 1: Induced-Fit Docking Endpoint Crash (HTTP 500)
- **Affected Runs:** Runs 21, 22, 23 (Category 7)
- **Exact Error:** `HTTP Error 500: INTERNAL SERVER ERROR`
- **Backend Traceback:**
  ```text
  TypeError: InducedFitService.run_induced_fit_docking() got an unexpected keyword argument 'ligand_sdf_or_pdbqt'.
  Did you mean 'ligand_smiles_or_pdbqt'?
  ```
- **Analysis:** In `backend/app.py` line 898, the endpoint calls `run_induced_fit_docking(..., ligand_sdf_or_pdbqt=...)`, but the method definition in `backend/services/induced_fit.py` line 185 defines the parameter as `ligand_smiles_or_pdbqt`. This keyword mismatch prevents induced-fit docking from executing through the REST API.

### Anomaly 2: False-Negative PAINS Filter on Curcumin
- **Affected Run:** Run 27 (Category 9)
- **Observed Result:** `pains_alerts_count: 0`, `status: "Clear"`
- **Expected Result:** Curcumin is one of the classic prototypical PAINS compounds documented by Baell & Holloway (J. Med. Chem. 2010). It MUST trigger a PAINS alert.
- **Analysis:** While Bindora's Brenk catalog correctly flagged curcumin with 2 alerts (`beta-keto/anhydride` and `Michael_acceptor_1`), the PAINS SMARTS library failed to detect the diferuloylmethane scaffold. In contrast, 5-benzylidenerhodanine was correctly flagged (`rhod_sat_A(33)`).

### Anomaly 3: False-Positive Blood-Brain Barrier Penetration for Loperamide
- **Affected Run:** Run 26 (Category 8)
- **Observed Result:** `bbb_permeation: "Permeant (Likely crosses BBB)"`
- **Known Pharmacology:** Loperamide (Imodium) is a peripheral μ-opioid receptor agonist that does not produce central opioid effects at therapeutic doses because it is an avid substrate for active P-glycoprotein (MDR1) efflux at the blood-brain barrier.
- **Analysis:** Bindora's BOILED-Egg implementation classifies loperamide as BBB permeant purely because its physicochemical coordinates (WLogP 5.09, TPSA 43.78 Å²) fall within the yolk ellipse. The platform lacks an active P-gp substrate classifier or efflux correction, creating a false-positive CNS penetration prediction for a known clinical non-CNS drug.

### Anomaly 4: Inverted Consensus Confidence on Decoy vs True Binder
- **Affected Runs:** Run 15 vs Run 16 (Category 5)
- **Observed Result:** Pentane (`CCCCC`, Vina = −2.897 kcal/mol) was classified as `HIGH_CONFIDENCE`, whereas Erlotinib (nanomolar EGFR inhibitor, Vina = −7.173 kcal/mol) was classified as `MODERATE_CONFIDENCE`.
- **Analysis:** In `backend/services/consensus.py`, the classification rule awards `HIGH_CONFIDENCE` when `rank_spread <= 2 and mean_rank <= 2.5`. For pentane, only a single mode was returned, resulting in a trivial rank spread of 0 across all available engines (Vina, Vinardo, MM-GBSA). Because GNINA was unavailable and pentane did not trigger the high-strain cutoff, the algorithm treated concordant poor binding across engines as "high confidence".

### Anomaly 5: Missing Seed Propagation in Standard Docking API
- **Affected Runs:** Runs 45, 46, 47 (Category 15)
- **Observed Result:** Three consecutive identical runs with `seed=42` produced drifting scores (−7.129, −7.230, −7.286 kcal/mol) and pose RMSDs up to 1.30 Å.
- **Analysis:** In `backend/app.py` line 331, `/api/docking/run` fails to extract `data.get("seed")` and pass it to `DockingEngine.run_docking()`. AutoDock Vina consequently initializes its random number generator from the system clock, destroying reproducible determinism for standard docking runs.

### Anomaly 6: Timeouts on Large Macrocycles (>80 Heavy Atoms)
- **Affected Runs:** Run 14 (Category 4), Runs 42 & 44 (Category 14)
- **Observed Result:** Cyclosporine A (85 heavy atoms, 33-membered ring) and Vancomycin (99 heavy atoms, MW 1449 Da) consistently timed out after 450–722 seconds.
- **Analysis:** AutoDock Vina's Iterated Local Search Monte Carlo algorithm scales poorly with high degrees of conformational freedom. While synthetic kinase macrocycles with constrained ring sizes (Lorlatinib, 30 heavy atoms) dock successfully within 50.6 seconds, cyclic undecapeptides and glycopeptides exceed practical limits for standard Vina without rigidifying peptide backbones.

### Anomaly 7: AI Narrative Pose Misattribution and Metric Conflation
- **Affected Runs:** Runs 48, 49, 50 (Category 16)
- **Observed Discrepancies:**
  - In Run 48, the AI narrative claimed the top pose had a high strain of 15.38 kcal/mol and MM-GBSA of −9.42 kcal/mol, whereas Mode 1 had a strain of 8.08 kcal/mol and MM-GBSA of −16.14 kcal/mol. The LLM conflated Mode 1 with Mode 8 and Mode 9.
  - In Run 49, the narrative reported "RMSD relative to the native reference (9.52/10.99 Å)", mistaking Vina's internal pose-dispersion metrics (`rmsd_lb` and `rmsd_ub` from Vina Mode 1) for a crystallographic comparison to native Indinavir.
  - In Run 50, the narrative propagated the default parameter `heavy_atoms: 20` for pentane (which has only 5 heavy atoms), reporting a false ligand efficiency of 0.145 kcal/mol/atom.

---

## 3. Reproducibility Verdict

Based on Category 15 (Runs 45, 46, 47), Bindora Dock currently **FAILS** the strict determinism test for standard docking runs executed via `/api/docking/run`:

1. **Numeric Drift:** Under identical user parameters (`target: 1M17`, `ligand: Erlotinib`, `exhaustiveness: 8`, `seed: 42`), the binding scores drifted from −7.129 kcal/mol to −7.286 kcal/mol (a variance of 0.157 kcal/mol).
2. **Pose Dispersion:** The top-ranked pose in Run 46 deviated by **1.30 Å heavy-atom RMSD** from Run 45. Run 47 deviated by **0.79 Å RMSD** from Run 45.
3. **Confidence Tier Oscillation:** Run 45 was classified as `MODERATE_CONFIDENCE`, while Runs 46 and 47 were classified as `HIGH_CONFIDENCE`.
4. **Root Cause:** The `seed` parameter is correctly accepted and passed in `/api/docking/redock-validate` (which is why Category 1 redocking is deterministic), but is dropped in `/api/docking/run`.

**Verdict:** For publication-grade or regulatory academic workflows, Bindora Dock cannot currently guarantee bit-identical pose reproducibility when using the primary `/api/docking/run` endpoint until the user-specified seed is forwarded to the underlying Vina process.

---

*Report compiled strictly from live localhost test execution data. No codebase files were modified.*
