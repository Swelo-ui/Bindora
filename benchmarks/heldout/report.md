# Bindora Dock: Independent Scientific Validation Audit Report
**Protocol:** Observation and Read-Only Empirical Benchmarking on Frozen Held-Out Datasets
**Execution Timestamp:** 2026-09-29T07:31:14.559468+00:00
**Audit Mode:** Read-Only Observer. Zero modifications to existing codebase files.

---

## 1. Executive Summary & Verification of Protocol Adherence

This independent audit was conducted under strict scientific observation constraints:
1. **Codebase Immutability:** No tracked files in the repository were altered, patched, or refactored during this audit. All evaluation scripts and data artifacts reside exclusively in `benchmarks/heldout/`.
2. **No Data from Memory:** Every molecular structure, canonical SMILES, and InChIKey was resolved via PubChem REST API lookup or downloaded official dataset files (TDC, B3DB).
3. **Freeze First, Run Second:** Benchmark datasets (`dev_set.csv`, `train.csv`, `test.csv`, `bbb_test.csv`, `pains_dataset.csv`) were frozen and hashed prior to running any validation scripts.
4. **Script-Generated Metrics:** All tables and confusion matrices reported below are generated programmatically from raw CSV outputs.

### Environment & Toolchain Metadata
- **Python:** `3.13.1 (tags/v3.13.1:0671451, Dec  3 2024, 19:06:28) [MSC v.1942 64 bit (AMD64)]`
- **RDKit:** `2026.03.6`
- **Meeko:** `0.8.0`
- **NumPy:** `2.5.3`
- **SciPy:** `1.18.1`
- **Scikit-learn:** `1.9.1`
- **AutoDock Vina Binary (`bin/vina.exe`) SHA-256:** `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5`

---

## 2. Step 0: Repository State & Static Pattern Audit

### 2.1 Git Working Tree Status
```
=== GIT REV-PARSE HEAD ===
325868c7a190a600bcb307d3a5d7ecd4e8759d39

=== GIT STATUS --SHORT ===
M backend/app.py
 M backend/services/adme.py
 M backend/services/consensus.py
 M backend/services/docking.py
 M backend/services/induced_fit.py
 M backend/services/macrocycle.py
 M backend/services/narrative.py
?? benchmarks/
?? data/benchmarks/qa_50run_raw_logs/
?? data/benchmarks/qa_50run_report.md

=== GIT DIFF HEAD --STAT ===
backend/app.py                  |  36 +++++++++--
 backend/services/adme.py        | 128 ++++++++++++++++++++++++++++++++++++++--
 backend/services/consensus.py   |  15 ++++-
 backend/services/docking.py     |  34 ++++++++---
 backend/services/induced_fit.py |  42 +++++++++----
 backend/services/macrocycle.py  |  30 ++++++++--
 backend/services/narrative.py   |  31 +++++++++-
 7 files changed, 279 insertions(+), 37 deletions(-)

=== UNTRACKED FILES SHA-256 ===
50128800ba5670c1c78bd247b356241e9d5ac3ea0f3b5577894e880b03f5bb41  benchmarks/heldout/00_working_tree.patch
fe203cf22b376dd6bdcdca24a2e49e4ffe8621539e26bb75a91715afc3e2e6e5  benchmarks/heldout/step0_repo_state.py
235603b74751bb784986ac6d70fb97957426a11a900c5b06de5790b6a9dd41c9  data/benchmarks/qa_50run_raw_logs/run_01_Native_Redocking_Indinavir_(MK1).json
944a295835b8291b3a39eaf7e0e64297dbe58373e6aad1d109a2ef6dc9234af2  data/benchmarks/qa_50run_raw_logs/run_02_Native_Redocking_Biotin_(BTN).json
58311956331cefbcebccbd89743deec97e7b8063587963c4de80dc5421aa564c  data/benchmarks/qa_50run_raw_logs/run_03_Native_Redocking_Benzamidine_(BEN).json
764d14b2b7efb16ce53d7ada137650099839bed3f56872c3486dfad6c0e8371f  data/benchmarks/qa_50run_raw_logs/run_04_Native_Redocking_Erlotinib_(AQ4).json
57746f960623800cf50e489f8ce4f2d3b86b26a04a30fc230d4d10025734691c  data/benchmarks/qa_50run_raw_logs/run_05_Native_Redocking_Imatinib_(STI).json
e15890cb28520dcb42806de4d21f171e48fe29a73706927e8eb23f8f44a1dde0  data/benchmarks/qa_50run_raw_logs/run_06_Cross_Docking_Aspirin.json
ee279c8234c4a0c5ee6413c2fc54015f04c48153f53cc253daeb5e7ee80c0c6a  data/benchmarks/qa_50run_raw_logs/run_07_Cross_Docking_Atorvastatin.json
95d8f12c4afd47734d8848d4e74011501f64572bb5639bd57361b058293f730d  data/benchmarks/qa_50run_raw_logs/run_08_Cross_Docking_Erlotinib.json
11ac28003b2ab33838d5f0d50e1c0e2cfc3a518047da179cd46612f5f71171d7  data/benchmarks/qa_50run_raw_logs/run_09_MMGBSA_Physics_Aniline.json
8896bd4c6f354830275fbe2e47348cd2e3c8c6b320c9f181fe3b419185a96e06  data/benchmarks/qa_50run_raw_logs/run_10_MMGBSA_Physics_Gefitinib.json
c3812aad4eb67bbdfa527690b1fc941ba056c752fcbebe2a0aa6e615c0c842aa  data/benchmarks/qa_50run_raw_logs/run_11_MMGBSA_Physics_Erlotinib.json
7259f7e55c0e8abb7f5340a46509aec575f09f16c127470fb5c7e12e6dc54293  data/benchmarks/qa_50run_raw_logs/run_12_Ligand_Strain_Benzamidine.json
2aa1891abe7b75ba200eaedf7512acd111c60d6c2aa7dc269f4800c58dd76ce9  data/benchmarks/qa_50run_raw_logs/run_13_Ligand_Strain_Atorvastatin.json
9a918be8c1e1e4d2f074eb82e73ec19019f41bc8224777a50e062aa0df468a1d  data/benchmarks/qa_50run_raw_logs/run_14_Ligand_Strain_Cyclosporine_A.json
638f38d3098ef55e10c250d343a00aeb1da0dacf960d219bdfbac76744e6401a  data/benchmarks/qa_50run_raw_logs/run_15_Consensus_Concordance_Erlotinib.json
cfdd54167b94157167f350a7f9d54cbc5a08668e17b5a04d8de5b3ae6bf42bd5  data/benchmarks/qa_50run_raw_logs/run_16_Consensus_Concordance_Pentane.json
6dfc8f60e64b6bce385ad906057294972d4fe4c5477f37458d25bb9b40674581  data/benchmarks/qa_50run_raw_logs/run_17_Consensus_Concordance_Aspirin.json
b3c0e91b45c5423d844c1355dd76823c6b9ece96ddb2b56587b5766872395151  data/benchmarks/qa_50run_raw_logs/run_18_Covalent_Warhead_Ibrutinib.json
77536b2ace627b249164213aa923fc28c436d95d21863258c8c80c85482afd06  data/benchmarks/qa_50run_raw_logs/run_19_Covalent_Warhead_Nirmatrelvir.json
e41694c1d73375fb59ba7d3573260aa7a4f99036a199b257a886e3c347c5b629  data/benchmarks/qa_50run_raw_logs/run_20_Covalent_Warhead_Afatinib.json
7f4c34be480624cf1034db7b6fe594b16cb7817d5fd96e746ef448c2cc481fe9  data/benchmarks/qa_50run_raw_logs/run_21_Induced_Fit_Imatinib.json
ed7f242eb2934ccfc1a14b6caab1cfc1cf0431b5d702a4c046119e7b451f00f6  data/benchmarks/qa_50run_raw_logs/run_22_Induced_Fit_Dasatinib.json
73845ed4b02278916a0bf9e2843b5717b37aea3425db4348d7097c1976250bf6  data/benchmarks/qa_50run_raw_logs/run_23_Induced_Fit_Erlotinib.json
9866ac7afd2971a7f962a12ef2a43a1707a040c38f1ed82171eb62b539e76968  data/benchmarks/qa_50run_raw_logs/run_24_ADME_Plausibility_Aspirin.json
7a90e283eed2f173cd6840bbc7c6d509ffeb82975ec1256a4df4e0f9de648dc8  data/benchmarks/qa_50run_raw_logs/run_25_ADME_Plausibility_Atorvastatin.json
c7c5f7ae32d99abc51e516881c75a87ac3847e642d17aeac91c20cba74cfbf2b  data/benchmarks/qa_50run_raw_logs/run_26_ADME_Plausibility_Loperamide.json
2724fbd1fe095ac7e1238788fc17447132eef7105f1aadfbe96a87f08cb31919  data/benchmarks/qa_50run_raw_logs/run_27_PAINS_Detection_Curcumin.json
226e49f883da651bc4e74f3ecd7b01e8f35be37506c41a6993aba74d2c201a3b  data/benchmarks/qa_50run_raw_logs/run_28_PAINS_Detection_Aspirin.json
0488f8565fd27fcaa1de621ad78e6e255c6a025c13fa0c0659cc484ac3809e77  data/benchmarks/qa_50run_raw_logs/run_29_PAINS_Detection_5-benzylidenerhodanine.json
0ee06d4b8305dc2087362b6b3ae40b072fb13f16747dd6d9d33665ec11c0fa96  data/benchmarks/qa_50run_raw_logs/run_30_Interaction_Fingerprint_Benzamidine.json
d129d1e0a4bcae5c753b1d99a24683423a22f76e0b24041e59fe1545cb5a4611  data/benchmarks/qa_50run_raw_logs/run_31_Interaction_Fingerprint_Erlotinib.json
ce5d0411eaeab462f074aa21f99534227cb82b5a0cc67f45f30edc7a3c1c0a67  data/benchmarks/qa_50run_raw_logs/run_32_Interaction_Fingerprint_Biotin.json
613da9cd8babb961923747f79d6c1c5eadadb5c5daad8c5be2f2744a7bcd0f30  data/benchmarks/qa_50run_raw_logs/run_33_Batch_Screening_Iteration_1.json
8f22b02ae9d8162c8716dcc5a8036e7d071775a0a614c164e077f3909978614a  data/benchmarks/qa_50run_raw_logs/run_34_Batch_Screening_Iteration_2.json
ed1007f8789d0332714e0f3bfe67bc10e4e90bcafdaec1ee777c2e3851c15d08  data/benchmarks/qa_50run_raw_logs/run_35_Batch_Screening_Iteration_3.json
d869dd2687fb25e47ea7d0c3ec4440cf0822f3a3779100762db2f1e8e60edfb0  data/benchmarks/qa_50run_raw_logs/run_36_ChEMBL_Crosscheck_Imatinib.json
f7b358c1ec0e8489d4c0a70969e75025a50fc3a0744e0523b3a2da1eda4357ff  data/benchmarks/qa_50run_raw_logs/run_37_ChEMBL_Crosscheck_Aspirin.json
c2869b18247a6b93d8072699d97dbc32db4ffdc0258a33c3be87c40f0afa5538  data/benchmarks/qa_50run_raw_logs/run_38_ChEMBL_Crosscheck_ZINC999999NonexistentMol.json
eb1140f2a708068d10a734535e13e44acb93b71e8b3d344600ae9654d59a1e5d  data/benchmarks/qa_50run_raw_logs/run_39_Pocket_Detection_1HSG.json
68190f2f879c153da3c62e25facc9e5db86ac53c8218e47599917a0dfe10c398  data/benchmarks/qa_50run_raw_logs/run_40_Pocket_Detection_3PTB.json
24b2c3a0309bccef0d3a98b1d77f92369cfc528e3de3048fe4f776eb913dac90  data/benchmarks/qa_50run_raw_logs/run_41_Pocket_Detection_1M17.json
c3be5f5616ce2c52aff80655540520eec33b41227b301632da26b1e4a0a0cfd1  data/benchmarks/qa_50run_raw_logs/run_42_Macrocycle_Handling_Cyclosporine_A.json
39696cb56cd4819ee1af67599966fc68f88b3ae8ba851415e53d30bea8e026e8  data/benchmarks/qa_50run_raw_logs/run_43_Macrocycle_Handling_Lorlatinib.json
9d6307c05c7f5482fe945c0f37e2560f55c2ad5640d6801566f8a5a6983aff8e  data/benchmarks/qa_50run_raw_logs/run_44_Macrocycle_Handling_Vancomycin.json
c42d263515194715da7d5730f0308cee04fcdd006f260531f3e34a1492bef397  data/benchmarks/qa_50run_raw_logs/run_45_Reproducibility_Determinism_Erlotinib_Rep_1.json
27978e3cb7f47bbac2d44c8aa5bb5c5078844b536655696c08d1923c11910663  data/benchmarks/qa_50run_raw_logs/run_46_Reproducibility_Determinism_Erlotinib_Rep_2.json
c4100da925f5a81ecbab132a8620978cc39523f2b38df7dc3b19f0f0cfff4eba  data/benchmarks/qa_50run_raw_logs/run_47_Reproducibility_Determinism_Erlotinib_Rep_3.json
8407e5df0ca7d6a3a4649fe11bab9fa72194a5ea1631a69bdec1823b5adceb0a  data/benchmarks/qa_50run_raw_logs/run_48_AI_Narrative_Consistency_Run_11_Erlotinib_1M17.json
7528ea389e1933ed726b78f97a2bbc4f21d63201b7f0123b7ab3446d7d421323  data/benchmarks/qa_50run_raw_logs/run_48_AI_Narrative_Consistency_Run_4_Erlotinib_1M17.json
dacfcdcb52191674397247dd680d59600b25f938f6fb723b3df79fabca2f7038  data/benchmarks/qa_50run_raw_logs/run_49_AI_Narrative_Consistency_Run_06_Aspirin_1HSG.json
55fe7cf67fbdd4b1e69a0c4c91bd796dddfeea0a1d10d39e0367b071925754be  data/benchmarks/qa_50run_raw_logs/run_49_AI_Narrative_Consistency_Run_1_Indinavir_1HSG.json
d3b410d4a095c7f20c080e498a0abf89bf4b2193928c676cb69534f9f7391bbf  data/benchmarks/qa_50run_raw_logs/run_50_AI_Narrative_Consistency_Run_16_Pentane_1M17.json
42f01aa0e904ba6c8d40d374844aa83467f523ec45424e7b516a86c703143d8d  data/benchmarks/qa_50run_report.md
```

### 2.2 Static Pattern Scan (Drug Names & SMILES Literals in Patch)
A regex scan across the uncommitted working tree patch was executed searching for drug names (`curcumin`, `loperamide`, `pentane`, `gefitinib`, `erlotinib`, `haloperidol`, `terfenadine`, etc.) and chemical literals:
```
=== PATTERN SCAN REPORT ===
Scanned files/patch for drug names: ['curcumin', 'loperamide', 'pentane', 'gefitinib', 'erlotinib', 'haloperidol', 'terfenadine', 'nicotine', 'metformin', 'tacrolimus', 'rapamycin', 'cyclosporin', 'lorlatinib', 'vancomycin', 'chalcone', 'maleimide', 'aspirin', 'caffeine', 'diazepam', 'morphine', 'isobutane']

Total findings: 12

[SMILES_LITERAL] backend/services/adme.py:96 (in logic)
  Match: C(=O)-C=C-C(=O)
  Line:  "smarts": "C(=O)-C=C-C(=O)",

[SMILES_LITERAL] backend/services/adme.py:101 (in logic)
  Match: c1cc(O)c(CN)cc1
  Line:  "smarts": "c1cc(O)c(CN)cc1",

[DRUG_NAME] backend/services/adme.py:119 (comment/description)
  Match: loperamide
  Line:  # 1. Gem-diphenyl / diphenylmethyl group (e.g. Loperamide, Terfenadine, Fendiline)

[DRUG_NAME] backend/services/adme.py:119 (comment/description)
  Match: terfenadine
  Line:  # 1. Gem-diphenyl / diphenylmethyl group (e.g. Loperamide, Terfenadine, Fendiline)

[SMILES_LITERAL] backend/services/adme.py:120 (in logic)
  Match: [#6](c1ccccc1)(c2ccccc2)
  Line:  p_diphenyl = Chem.MolFromSmarts("[#6](c1ccccc1)(c2ccccc2)")

[DRUG_NAME] backend/services/adme.py:124 (comment/description)
  Match: loperamide
  Line:  # 2. 4-Arylpiperidine / piperazine basic system (e.g. Loperamide, Haloperidol)

[DRUG_NAME] backend/services/adme.py:124 (comment/description)
  Match: haloperidol
  Line:  # 2. 4-Arylpiperidine / piperazine basic system (e.g. Loperamide, Haloperidol)

[SMILES_LITERAL] backend/services/adme.py:125 (in logic)
  Match: c1ccccc1-C1CCNCC1
  Line:  p_aryl_pip = Chem.MolFromSmarts("c1ccccc1-C1CCNCC1")

[SMILES_LITERAL] backend/services/adme.py:126 (in logic)
  Match: N1CCC(CC1)
  Line:  p_pip_ring = Chem.MolFromSmarts("N1CCC(CC1)")

[DRUG_NAME] backend/services/adme.py:138 (comment/description)
  Match: tacrolimus
  Line:  # 4. Large macrocyclic polyketide / cyclic peptide (e.g. Cyclosporine, Ivermectin, Tacrolimus)

[DRUG_NAME] backend/services/macrocycle.py:60 (comment/description)
  Match: tacrolimus
  Line:  # 2. Algebraic cycle combination for fused/bridged macrocycles (e.g. Tacrolimus, Rapamycin)

[DRUG_NAME] backend/services/macrocycle.py:60 (comment/description)
  Match: rapamycin
  Line:  # 2. Algebraic cycle combination for fused/bridged macrocycles (e.g. Tacrolimus, Rapamycin)
```
**Audit Finding on Static Patterns:**
- **Drug Names:** Drug names in the patch occur exclusively inside non-executable code comments or docstrings as illustrative pharmacophoric examples (e.g. `# 1. Gem-diphenyl / diphenylmethyl group (e.g. Loperamide, Terfenadine, Fendiline)` in `adme.py:119`). No drug names are present in branch logic or conditional checks.
- **SMARTS Literals:** Concrete SMARTS patterns exist in executable logic for generalized structural motifs: `C(=O)-C=C-C(=O)` (bis-enone / ene-dione), `c1cc(O)c(CN)cc1` (phenol-Mannich), `[#6](c1ccccc1)(c2ccccc2)` (diphenylmethyl core), `c1ccccc1-C1CCNCC1` (4-arylpiperidine), and `N1CCC(CC1)` (piperidine heterocycle).

---

## 3. Step 1: Frozen Benchmark Datasets & Integrity Hashes

Prior to evaluation, 5 held-out datasets were constructed, resolved via PubChem REST, and frozen in `benchmarks/heldout/FROZEN_HASHES.txt`:

| Dataset File | Record Count | Description | SHA-256 Hash |
| :--- | :--- | :--- | :--- |
| `benchmarks/heldout/dev_set.csv` | 42 | dev_set | `3fce1b013eb3770a54cd51d9208dbf2f85ae7d99a87563bb2b6e0beec1e5bc58` |
| `benchmarks/heldout/train.csv` | 975 | pgp_train | `3d7ad170299233164a551d33a5e6ecda43420dee97fd2f35c43c0e291904cc8d` |
| `benchmarks/heldout/test.csv` | 244 | pgp_test | `ed85c4c64124f9c02435de70d5506d38d689c91699cb000b11a98cbfc68a396f` |
| `benchmarks/heldout/bbb_test.csv` | 7782 | bbb_test | `579a470255319d5203ec7dd5e97adf607507f3302c69f40f74d25dc5b0037c78` |
| `benchmarks/heldout/pains_dataset.csv` | 67 | pains_dataset | `fdb7da68f7f187b5edb65761ec0ef8754d8e0a078568493fe8ec3c33ba6e6e47` |

**Data Curation Details:**
- **`dev_set.csv` (42 compounds):** Compounds present in earlier 50-run QA suites and prompt stress tests (Indinavir, Aspirin, Erlotinib, Gefitinib, Curcumin, Loperamide, Haloperidol, Tacrolimus, Rapamycin, etc.), each resolved with official PubChem CID, canonical SMILES, and InChIKey.
- **`train.csv` (975 compounds) & `test.csv` (244 compounds):** TDC `Pgp_Broccatelli` (1,219 valid compounds) partitioned via Bemis-Murcko scaffold split (80/20, `random_state=42`). 7 dev-set compounds were naturally isolated to the training split; zero dev-set compounds remained in the test split.
- **`bbb_test.csv` (7,782 compounds):** Full experimental B3DB classification dataset downloaded from GitHub (`theochem/B3DB`). 23 dev-set compounds were excluded.
- **`pains_dataset.csv` (67 compounds):** 32 literature PAINS positives across 8 chemical families (Baell & Holloway 2010) and 35 hard negatives (approved drugs and non-interfering analogues). Every entry was retrieved directly from PubChem REST.

---

## 4. Step 2: Macrocycle Perception Evaluation

Evaluation of 5 complex macrocycles: Cyclosporine A, Tacrolimus, Rapamycin, Lorlatinib, and Vancomycin.

| Compound | PubChem InChIKey | Literature Ring Size | RDKit Raw SSSR | Bindora All Macro Sizes (>=12) | Lit Size in Bindora? | Bridge Contraction? |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Cyclosporine A | `PMATZTZNYRCHOR-CGLBZJNRSA-N` | 33 (33-membered monocyclic undecapeptide) | `[33]` | `[33]` | **True** | False |
| Tacrolimus (FK506) | `QJJXYPPXXYFBGM-LFZNUXCKSA-N` | 23 (23-membered macrolide lactone (IUPAC: 1,14-oxa-azabicyclo[19.3.1]pentacosane core)) | `[27, 6, 6, 6]` | `[33, 27, 25, 23, 12]` | **True** | False |
| Rapamycin (Sirolimus) | `QFJCIRLUMZQUOT-HPLJOQBZSA-N` | 31 (31-membered macrolide lactone (pipecolate-containing macrolide)) | `[35, 33, 6, 6]` | `[41, 39, 35, 33, 31, 29, 12]` | **True** | False |
| Lorlatinib | `IIXWYSCJSQVBQM-LLVKDONJSA-N` | 12 (12-membered bridged kinase macrocycle (excluding fused pyrazole atoms) or 15-membered envelope) | `[15, 6, 6, 5]` | `[19, 17, 15, 12]` | **True** | False |
| Vancomycin | `MYPYJXKWCTUITO-LYRMYLQWSA-N` | 16 (Tricyclic glycopeptide core with 16-membered, 12-membered, and 12-membered crosslinked peptide/ether rings) | `[36, 16, 12, 6, 6, 6, 6, 6, 6, 6]` | `[42, 40, 38, 36, 34, 30, 28, 22, 18, 16, 14, 12]` | **True** | False |

### Mechanistic Analysis of Macrocycle Detection:
1. **SSSR Limitations in Complex Fused Systems:** Standard RDKit SSSR (`GetRingInfo().AtomRings()`) constructs an arbitrary minimal cycle basis. In Tacrolimus, the 6-membered pipecolic acid ring is fused to the macrolactone, causing SSSR to report `[27, 6, 6, 6]`. In Rapamycin, SSSR reports `[35, 33, 6, 6]`. Thus, raw SSSR alone fails to isolate the true macrolide lactone perimeter.
2. **Algebraic Cycle Combination (`edge_basis[i] ^ edge_basis[j]`):** Bindora combines basis cycle edges pairwise. When an edge combination forms a simple cycle (all vertices have degree == 2), it is retained. Consequently:
   - For Tacrolimus: Bindora detects `[33, 27, 25, 23, 12]`. The official 23-membered macrolide lactone is captured.
   - For Rapamycin: Bindora detects `[41, 39, 35, 33, 31, 29, 12]`. The official 31-membered macrolide lactone is captured.
   - For Lorlatinib: SSSR reports `[15, 6, 6, 5]`, while Bindora detects `[19, 17, 15, 12]`, capturing both the 12-membered bridged perimeter and the 15-membered pyrazole envelope.
   - For Cyclosporine A: Unfused monocyclic peptide reports exactly `[33]`.
3. **Bridge Contraction:** Bindora does **NOT** contract ester or amide bonds into virtual single edges; detection is purely topological cycle algebra on the atomic graph.
4. **Envelope vs Macrocycle Identification:** Because pairwise XOR produces both the smaller perimeter and outer composite cycles, `max_ring_size` reports the composite outer envelope (e.g. 33 for Tacrolimus, 41 for Rapamycin), while `macrocycle_ring_sizes` contains the true lactone perimeters.

---

## 5. Step 3: PAINS Benchmark Evaluation

Evaluated on frozen `pains_dataset.csv` (N = 67: 32 literature positives, 35 hard negatives).

### 5.1 Performance Comparison Table
| Model | TP | FP | TN | FN | Sensitivity [95% Wilson CI] | Specificity [95% Wilson CI] | Balanced Accuracy | MCC |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **RDKit_PAINS_ABC_Alone** | 11 | 0 | 35 | 21 | 0.3438 [0.2041, 0.5169] | 1.0000 [0.9011, 1.0000] | 0.6719 | 0.4635 |
| **Bindora_Full** | 14 | 0 | 35 | 18 | 0.4375 [0.2817, 0.6067] | 1.0000 [0.9011, 1.0000] | 0.7188 | 0.5375 |
| **Bindora_Extended_Only** | 8 | 0 | 35 | 24 | 0.2500 [0.1325, 0.4211] | 1.0000 [0.9011, 1.0000] | 0.6250 | 0.3851 |

### 5.2 False Positives & False Negatives Analysis
- **False Positives (FP = 0):** Across all 35 hard negatives (including Chalcone, Cinnamamide, Dibenzoylmethane, Ferulic acid, Capsaicin, Aspirin, Ibuprofen, Paracetamol, Metformin, Caffeine, etc.), **zero false positives** were generated by RDKit alone, Bindora Extended, or Bindora Full. Specificity is 1.0000 (95% CI: 0.9011 - 1.0000).
- **False Negatives (FN = 18 in Bindora Full):**
  Bindora Full missed 18 PAINS compounds:
  - `5-(4-hydroxybenzylidene)rhodanine` (Family: *rhodanine*) — RDKit alert: none, Extended alert: none.
  - `5-(4-(dimethylamino)benzylidene)rhodanine` (Family: *rhodanine*) — RDKit alert: none, Extended alert: none.
  - `Epalrestat` (Family: *rhodanine*) — RDKit alert: none, Extended alert: none.
  - `5-(4-chlorobenzylidene)rhodanine` (Family: *rhodanine*) — RDKit alert: none, Extended alert: none.
  - `5-(2-furylmethylene)rhodanine` (Family: *rhodanine*) — RDKit alert: none, Extended alert: none.
  - `Dehydrozingerone` (Family: *enone*) — RDKit alert: none, Extended alert: none.
  - `2-hydroxybenzaldehyde phenylhydrazone` (Family: *hydroxyphenyl_hydrazone*) — RDKit alert: none, Extended alert: none.
  - `Salicylaldehyde 2-nitrophenylhydrazone` (Family: *hydroxyphenyl_hydrazone*) — RDKit alert: none, Extended alert: none.
  - `5-benzylidenebarbituric acid` (Family: *alkylidene_barbiturate*) — RDKit alert: none, Extended alert: none.
  - `5-(4-hydroxybenzylidene)barbituric acid` (Family: *alkylidene_barbiturate*) — RDKit alert: none, Extended alert: none.
  - `5-(4-dimethylaminobenzylidene)barbituric acid` (Family: *alkylidene_barbiturate*) — RDKit alert: none, Extended alert: none.
  - `1,2-Naphthoquinone` (Family: *quinone*) — RDKit alert: none, Extended alert: none.
  - `Coenzyme Q0` (Family: *quinone*) — RDKit alert: none, Extended alert: none.
  - `2-(dimethylaminomethyl)phenol` (Family: *phenol_mannich*) — RDKit alert: none, Extended alert: none.
  - `4-(morpholinomethyl)phenol` (Family: *phenol_mannich*) — RDKit alert: none, Extended alert: none.
  - `alpha-cyano-4-hydroxycinnamic acid` (Family: *ene_cyano*) — RDKit alert: none, Extended alert: none.
  - `Pyrogallol` (Family: *catechol*) — RDKit alert: none, Extended alert: none.
  - `Octhilinone` (Family: *isothiazolone*) — RDKit alert: none, Extended alert: none.

**Root Cause of False Negatives:**
1. *Rhodanines:* Standard RDKit FilterCatalog and Bindora's `rhodanine_expanded` SMARTS (`O=C1CSC(=[S,O])N1`) check for saturated C5 (`CSC`), whereas active 5-benzylidenerhodanines have an exocyclic alkene at C5 (`C(=C)SC`), causing them to evade the SMARTS pattern.
2. *Phenol-Mannich Bases:* Bindora's pattern `c1cc(O)c(CN)cc1` requires a primary amine (`CN`), but classical Mannich bases (e.g. 2-(dimethylaminomethyl)phenol) contain tertiary amines (`CN(C)C`).
3. *Hydrazones & Barbiturates:* RDKit's built-in PAINS catalog in FilterCatalog implements a limited subset of the original 480 Wehi/Baell SMARTS.

---

## 6. Step 4: P-gp Substrate Classifier Benchmark

Evaluated on frozen TDC `Pgp_Broccatelli` held-out test split (`test.csv`, N = 244: 152 substrates, 92 non-substrates).

### 6.1 Overall Test Set Performance
| Model | TP | FP | TN | FN | Sensitivity [95% CI] | Specificity [95% CI] | Balanced Accuracy | MCC | ROC-AUC |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Majority_Class** | 152 | 92 | 0 | 0 | 1.0000 [0.9754, 1.0000] | 0.0000 [0.0000, 0.0401] | 0.5000 | 0.0000 | 0.5000 |
| **Heuristic_MW_LogP** | 74 | 2 | 90 | 78 | 0.4868 [0.4087, 0.5656] | 0.9783 [0.9242, 0.9940] | 0.7326 | 0.4867 | 0.7326 |
| **Logistic_Regression_MorganFP** | 131 | 26 | 66 | 21 | 0.8618 [0.7980, 0.9078] | 0.7174 [0.6181, 0.7992] | 0.7896 | 0.5861 | 0.8942 |
| **Bindora_Rule_Based** | 58 | 3 | 89 | 94 | 0.3816 [0.3082, 0.4608] | 0.9674 [0.9085, 0.9888] | 0.6745 | 0.3906 | 0.6745 |

### 6.2 Subgroup Analysis: Diphenylmethyl / Piperidine Chemotypes vs Other Scaffolds

Test set partitioned into:
1. **Diphenylmethyl / Piperidine Chemotypes (N = 41):** Scaffolds matching `[#6](c1ccccc1)(c2ccccc2)` or `N1CCC(CC1)` (e.g. Loperamide-like / Haloperidol-like systems).
2. **Other Scaffolds (N = 203):** All other chemical scaffolds in the held-out test set.

| Subgroup | Model | Sensitivity | Specificity | Balanced Accuracy | MCC |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Diphenylmethyl / Piperidine (N=41)** | Bindora Rule-Based | **0.7727** | **0.8947** | **0.8337** | **0.6675** |
| | Logistic Regression (Morgan FP) | 0.7273 | 0.5263 | 0.6268 | 0.2592 |
| | Heuristic (MW>400, LogP>3) | 0.5000 | 0.9474 | 0.7237 | 0.4903 |
| **Other Scaffolds (N=203)** | Bindora Rule-Based | **0.3154** | **0.9863** | **0.6508** | **0.3574** |
| | Logistic Regression (Morgan FP) | **0.8846** | **0.7671** | **0.8259** | **0.6558** |
| | Heuristic (MW>400, LogP>3) | 0.4846 | 0.9863 | 0.7355 | 0.4864 |

### Inductive Bias & Generalization Findings:
1. **High Specificity, Low Global Sensitivity:** The Bindora rule-based model operates with high specificity (0.9674 on test set, 0.9863 on other scaffolds), meaning it rarely generates false positive substrate calls.
2. **Targeted Inductive Bias:** On the diphenylmethyl / piperidine sub-family, Bindora achieves **0.7727 Sensitivity** and **0.8337 Balanced Accuracy** (MCC 0.6675), outperforming Logistic Regression.
3. **Generalization Gap:** On general scaffolds outside this motif, Bindora's sensitivity drops sharply to **0.3154** (missing 68.5% of true substrates), whereas Logistic Regression maintains **0.8846 Sensitivity** (Balanced Accuracy 0.8259, MCC 0.6558).

---

## 7. Step 5: BBB Decoupling Benchmark

Evaluated on frozen `bbb_test.csv` (N = 7,782 experimental compounds from B3DB: 4,942 BBB+, 2,840 BBB-).

### 7.1 Confusion Matrix & Comparative Performance

| Model Architecture | TP | FP | TN | FN | Sensitivity [95% CI] | Specificity [95% CI] | Balanced Accuracy | MCC |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Bindora Current (Decoupled Passive BOILED-Egg)** | 3009 | 486 | 2354 | 1933 | 0.6089 [0.5952, 0.6224] | 0.8289 [0.8146, 0.8423] | **0.7189** | **0.4237** |
| **Old Coupled Rule (`if PGP: BBB=False`)** | 2543 | 437 | 2403 | 2399 | 0.5146 [0.5006, 0.5285] | 0.8461 [0.8324, 0.8589] | 0.6803 | 0.3572 |

### 7.2 Quantitative Impact of Coupling
- **True CNS Drugs Falsely Excluded by Coupling:** Exactly **466 true permeant compounds** in B3DB were erroneously flipped from true positive to false negative when P-gp substrate predictions were allowed to overwrite BBB permeability.
- **Sensitivity Penalty:** Coupling drops sensitivity by **9.43 percentage points** (from 60.89% down to 51.46%).
- **MCC Degradation:** Matthews Correlation drops from 0.4237 to 0.3572.

### 7.3 Representative True CNS Permeants Falsely Excluded by Old Coupled Rule
| Compound Name | InChIKey | TPSA (Å²) | LogP | Triggered P-gp Exclusion Rule |
| :--- | :--- | :--- | :--- | :--- |
| proc-19m | `RLJKFXRSANDEJS-UHFFFAOYSA-N` | 53.09 | 4.73 | High lipophilicity (LogP 4.7 > 2.8), MW (477.6 > 400 Da), and presence of Piperidine basic heterocycle, Basic aliphatic nitrogen center. |
| - | `LIVAIIMPQYHZOR-UHFFFAOYSA-N` | 70.26 | 4.01 | High lipophilicity (LogP 4.0 > 2.8), MW (517.7 > 400 Da), and presence of Diphenylmethyl / bis-aryl lipophilic core, Piperidine basic heterocycle, Basic aliphatic nitrogen center. |
| chembl1185227 | `XMKPWAJOBWNRGY-UHFFFAOYSA-N` | 73.4 | 3.42 | High lipophilicity (LogP 3.4 > 2.8), MW (516.6 > 400 Da), and presence of Piperidine basic heterocycle, Basic aliphatic nitrogen center. |
| 161105-56-2 | `FBOIFJGFNUAMIT-UHFFFAOYSA-N` | 44.73 | 5.31 | High lipophilicity (LogP 5.3 > 2.8), MW (470.7 > 400 Da), and presence of Diphenylmethyl / bis-aryl lipophilic core, Piperidine basic heterocycle, Basic aliphatic nitrogen center. |
| domperidone | `FGXWKSZFVQUSTL-UHFFFAOYSA-N` | 78.82 | 3.35 | High lipophilicity (LogP 3.4 > 2.8), MW (425.9 > 400 Da), and presence of Piperidine basic heterocycle, Basic aliphatic nitrogen center. |
| risperidone-9-oh | `PNKSVGFMOKGXFK-UHFFFAOYSA-N` | 71.5 | 3.69 | High lipophilicity (LogP 3.7 > 2.8), MW (425.5 > 400 Da), and presence of Piperidine basic heterocycle, Basic aliphatic nitrogen center. |
| levocabastine hydrochloride | `ZCGOMHNNNFPNMX-UHFFFAOYSA-N` | 64.33 | 4.89 | High lipophilicity (LogP 4.9 > 2.8), MW (420.5 > 400 Da), and presence of 4-Arylpiperidine basic pharmacophore, Basic aliphatic nitrogen center. |
| dexverapamil | `SGTNSNPWRIOYBX-HHHXNRCGSA-N` | 63.95 | 5.09 | High lipophilicity (LogP 5.1 > 2.8), MW (454.6 > 400 Da), and presence of Basic aliphatic nitrogen center. |
| chembl239099 | `IFFOAHBQMKKYRJ-UHFFFAOYSA-N` | 76.88 | 3.44 | High lipophilicity (LogP 3.4 > 2.8), MW (458.6 > 400 Da), and presence of Piperidine basic heterocycle, Basic aliphatic nitrogen center. |
| ketoconazole | `XMAYWYJOQHXEEK-ANWICMFUSA-N` | 69.06 | 4.21 | High lipophilicity (LogP 4.2 > 2.8), MW (531.4 > 400 Da), and presence of Basic aliphatic nitrogen center. |

### Scientific Conclusion on BBB Decoupling:
Passive membrane permeability (governed by physicochemical factors: TPSA, WLogP, polar/apolar balance) and active efflux (governed by ABCB1 transporter binding) represent distinct physical mechanisms. Conflating them into a single binary boolean creates severe systematic under-prediction of CNS drug candidates.

---

## 8. Summary of Objective Findings & Limitations

1. **No Evidence of Name-Based Hardcoding:** String and SMARTS audits show no drug-name conditional branching or target-specific overrides. Drug names appear exclusively as illustrative comments.
2. **Generalization Gap in Rule-Based P-gp:** The P-gp predictor is strongly effective on diphenylmethyl/4-arylpiperidine chemotypes (Balanced Acc: 83.4%, MCC: 0.668), but exhibits significant under-recall on broader chemical scaffolds (Sensitivity: 31.5%). A statistical machine learning model (e.g. Logistic Regression on Morgan fingerprints) achieves far superior generalization (Balanced Acc: 82.6%, MCC: 0.656).
3. **Macrocycle Logic Operates on Cycle Algebra:** The algebraic bond XOR approach successfully recovers the true perimeter ring sizes of fused macrocycles (Tacrolimus 23, Rapamycin 31, Lorlatinib 12), resolving the prior SSSR under-reporting bug, while also capturing outer envelope cycles.
4. **Decoupling Validated by B3DB:** Decoupling passive BOILED-Egg BBB permeation from active P-gp efflux is quantitatively substantiated, rescuing 466 true CNS drugs from false negative classification.

---
*End of Independent Scientific Validation Audit Report.*