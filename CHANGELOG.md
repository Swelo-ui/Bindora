# Bindora Dock: Research-Grade Engineering Changelog

All notable technical fixes, scientific upgrades, and defect remediations are documented here.
Format: `[Fix ID] What + Why + Test added`.

---

## [Fix C1] Elimination of Silent Default Seed 42 and Dynamic Random Seed Generation
- **What:** Removed hardcoded silent default `seed = 42` across `DockingEngine.run_docking`, `DockingEngine.run_redocking_validation`, and API endpoints (`/api/dock`, `/api/docking/redock-validate`, `/api/induced-fit`). When `seed` is absent or `None`, the engine generates an explicit random integer seed using `secrets.randbelow(2147483647) + 1`, passes it to Vina, and reports it in the response payload.
- **Why:** Silently defaulting to seed 42 hid stochastic sampling variance and created artificial reproducibility across docking runs, violating independent scientific auditing standards.
- **Test Added:** `tests/test_audit_fixes.py::test_c1_seed_not_default_42_when_absent`, `test_c1_seed_generation_in_docking`.
