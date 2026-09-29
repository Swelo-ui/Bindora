# Bindora Dock: Research-Grade Engineering Changelog

All notable technical fixes, scientific upgrades, and defect remediations are documented here.
Format: `[Fix ID] What + Why + Test added`.

---

## [Fix C1] Elimination of Silent Default Seed 42 and Dynamic Random Seed Generation
- **What:** Removed hardcoded silent default `seed = 42` across `DockingEngine.run_docking`, `DockingEngine.run_redocking_validation`, and API endpoints (`/api/docking/run`, `/api/docking/redock-validate`, `/api/docking/induced-fit`). When `seed` is absent or `None`, the engine generates an explicit random integer seed using `secrets.randbelow(2147483647) + 1`, passes it to Vina, and reports it in the response payload.
- **Why:** Silently defaulting to seed 42 hid stochastic sampling variance and created artificial reproducibility across docking runs, violating independent scientific auditing standards.
- **Test Added:** `tests/test_audit_fixes.py::test_c1_seed_not_default_42_when_absent`, `test_c1_seed_generation_in_docking`.

## [Fix C2] Strict RDKit-Driven Ligand Descriptors & 400 On Parse Failure
- **What:** Heavy atoms count and molecular weight (MW) are now calculated strictly and dynamically from sanitized RDKit molecules (`mol.GetNumHeavyAtoms()`, `Descriptors.MolWt(mol)`). Eliminated hardcoded fallback defaults (`20` heavy atoms and `300.0/350.0` MW). Any invalid or unparseable ligand SMILES causes an immediate HTTP 400 rejection with explicit error details and `null` descriptor fields.
- **Why:** Fallback numbers and trusting unvalidated client inputs masked malformed chemical structures and produced misleading thermodynamic estimates.
- **Test Added:** `tests/test_audit_fixes.py::test_c2_ligand_descriptors_computed_strictly_from_rdkit`, `test_c2_app_endpoint_parse_failure_returns_400`.

## [Fix C3] Endpoint Validation, Alias Normalization, and Structured Error Handling
- **What:** Audited all API endpoints (`/api/docking/run`, `/api/docking/classify-experiment`, `/api/docking/analyze-interactions`, `/api/docking/interaction-diagram`, `/api/docking/refine`, `/api/docking/redock-validate`, `/api/adme/profile`, `/api/induced-fit`, `/api/batch-screen`, `/api/narrative`). Implemented centralized JSON payload parsing via `get_request_json()` which validates JSON syntax and enforces dictionary payload structure, returning structured HTTP 400 JSON on malformed JSON syntax or non-dict payloads (e.g. lists). Added route aliases `/api/adme/profile`, `/api/induced-fit`, `/api/batch-screen`, and `/api/narrative`. Registered global structured JSON handlers for `HTTPException` and `Exception` on `/api/` routes, ensuring no endpoint ever leaks unhandled 500 HTML error pages.
- **Why:** Malformed JSON syntax, non-dict JSON bodies, and missing parameters caused uncaught `AttributeError` 500s or returned default HTML error pages, crashing client integrations.
- **Test Added:** `tests/test_audit_fixes.py::test_c3_endpoint_aliases_exist`, `test_c3_invalid_json_syntax_returns_400_not_500`, `test_c3_non_dict_json_returns_400_not_500`, `test_c3_empty_or_missing_parameters_return_400`.

## [Fix C4] Induced-Fit Docking Improvements, Dynamic Random Seed, and Honest Backbone Reporting
- **What:** In `InducedFitService`, eliminated default `seed = 42` across sampling and docking, generating dynamic non-deterministic integer seeds (`secrets.randbelow`) and reporting `seed_used` in responses. Upgraded receptor perturbation analysis to explicitly separate `backbone_rmsd_angstroms` from `sidechain_rmsd_angstroms`. Added transparent biophysical movement assessment (`backbone_movement_assessment`) explicitly flagging when backbone displacement is minimal (< 0.25 Å) and adaptation is driven primarily by sidechain breathing. Tested against real QA targets (Imatinib/1IEP, Dasatinib/1IEP, Erlotinib/1M17) using verified PubChem SMILES.
- **Why:** Masking minimal backbone displacement as large-scale induced fit produced misleading scientific claims; hardcoded seeds concealed Monte Carlo sampling variability.
- **Test Added:** `tests/test_audit_fixes.py::test_c4_seed_signature_and_generation`, `test_c4_honest_backbone_reporting`, `test_c4_real_qa_targets_sampling`.


