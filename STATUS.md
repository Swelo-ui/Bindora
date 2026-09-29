# Bindora Dock: Research-Grade Engineering & Validation Status

**Repository:** `h:\AnuDock`  
**Current Branch:** `research-grade-fixes`  
**Baseline Tag:** `baseline-2026-09-29`  
**Updated:** 2026-09-29T07:40:00Z  

---

## Stage Progress & Gates

| Stage | Gate Description | Status | Evidence / Notes |
| :--- | :--- | :--- | :--- |
| **STAGE A** | Branch `research-grade-fixes` + tag `baseline-2026-09-29` exist; full `backend/` scan completed | **PASSED** | Branch & tag verified; `A1_repo_state.txt` and `A2_backend_audit_scan.txt` generated (Drugs=24, SMILES=135, Numeric defaults=48). |
| **STAGE B** | Benchmark harness rebuilt: P-gp substrate dataset verified, PAINS verified with PubChem titles/synonyms & raw RDKit matches, Macrocycle discrepancy table & cycle graph connectivity, BBB sensitivity, FROZEN_HASHES re-verified | **PASSED** | Wang 2011 P-gp substrate dataset (327 valid, 80/20 scaffold split); 65 verified PAINS compounds with PubChem title matches & raw RDKit dump; macrocycle discrepancy & simple-cycle connectivity check; BBB sensitivity across 4 tiers; all 6 frozen files 100% re-verified. |
| **STAGE C** | Code fixes C1–C15 committed individually with failing-then-passing pytest tests and CHANGELOG | *IN PROGRESS* | C1 (Seed randomness), C2 (RDKit ligand descriptors), and C3 (Endpoint validation, alias normalization, structured JSON errors) completed and passing tests. |
| **STAGE D** | Docking regression suite (Erlotinib, CsA, Lorlatinib, IFD, Narrative) with raw logs and manifests | *PENDING* | Awaiting Gate C. |
| **STAGE E** | Scientific benchmarks executed on frozen held-out sets (PAINS, P-gp, BBB, DUD-E decoy gate) | *PENDING* | Awaiting Gate D. |
| **STAGE F** | Final deliverables: run_all.py, CSVs, manifest.json, STATUS.md, report.md, CHANGELOG.md | *PENDING* | Awaiting Gate E. |

---

## Unverified Claims & Critical Observations Log
- **P-gp Broccatelli Dataset Label:** User correctly observed that TDC `Pgp_Broccatelli` is an **inhibition** dataset (measuring whether a compound inhibits ABCB1 ATPase or calcein-AM efflux), NOT a substrate dataset. The previous audit benchmark evaluated Bindora's substrate logic against inhibitor labels. Stage B1 must obtain a true documented substrate dataset.
- **RDKit PAINS Catalog:** The hypothesis that RDKit FilterCatalog contains only a "limited subset" was investigated. FilterCatalog loads `PAINS_A`, `PAINS_B`, and `PAINS_C` (480 total patterns). The missed rhodanines and catechol/pyrogallol require exact SMARTS matching inspection in Stage B2.
