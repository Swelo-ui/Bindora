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
