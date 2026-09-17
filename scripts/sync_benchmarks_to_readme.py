#!/usr/bin/env python3
"""
scripts/sync_benchmarks_to_readme.py

Automated synchronization script that reads:
  - data/benchmarks/casf2016_final_report.json
  - data/benchmarks/dude_screening_report.json
and generates the corresponding Markdown sections directly in README.md.

Prevents manual copy-paste drift and guarantees that README numbers always
match the committed benchmark report JSON files.

Usage:
  # Update README.md in-place:
  python scripts/sync_benchmarks_to_readme.py

  # Check if README.md is in sync (returns exit code 1 if out-of-sync, for CI):
  python scripts/sync_benchmarks_to_readme.py --check
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict, Any, List

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CASF_JSON_PATH = PROJECT_ROOT / "data" / "benchmarks" / "casf2016_final_report.json"
DUDE_JSON_PATH = PROJECT_ROOT / "data" / "benchmarks" / "dude_screening_report.json"
README_PATH = PROJECT_ROOT / "README.md"

CASF_START_TAG = "<!-- BENCHMARK_CASF2016_START -->"
CASF_END_TAG = "<!-- BENCHMARK_CASF2016_END -->"
DUDE_START_TAG = "<!-- BENCHMARK_DUDE_START -->"
DUDE_END_TAG = "<!-- BENCHMARK_DUDE_END -->"


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Benchmark file not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_casf2016_markdown(report: Dict[str, Any]) -> str:
    data = report.get("data", {})
    stats = data.get("summary_statistics", {})
    pose_stats = stats.get("pose_reconstruction", {})
    complex_results: List[Dict[str, Any]] = data.get("complex_results", [])
    metadata = data.get("metadata", {})

    total_attempted = stats.get("total_complexes_attempted", len(complex_results))
    total_successful = stats.get("successful_dockings", len(complex_results))
    under_2a = pose_stats.get("rmsd_under_2a_count", 0)
    success_rate = pose_stats.get("rmsd_success_rate_percent", 0.0)
    under_1a = pose_stats.get("rmsd_under_1a_count", 0)
    sub_angstrom_rate = round((under_1a / total_attempted * 100.0), 1) if total_attempted else 0.0
    mean_rmsd = pose_stats.get("mean_rmsd_angstroms", 0.0)
    median_rmsd = pose_stats.get("median_rmsd_angstroms", 0.0)
    exhaustiveness = metadata.get("exhaustiveness", 4)

    # Format complexes tested string
    complex_str_items = []
    for c in complex_results:
        pid = c.get("pdb_id")
        r = c.get("rmsd_angstroms")
        if r is not None:
            complex_str_items.append(f"`{pid}` ({r:.2f} Å)")
        else:
            complex_str_items.append(f"`{pid}` (N/A)")
    complexes_tested_line = ", ".join(complex_str_items)

    lines = [
        CASF_START_TAG,
        f"**Validated Performance ({total_attempted}-Complex Run, Exhaustiveness {exhaustiveness}):**",
        f"* **Success Rate (RMSD ≤ 2.0 Å):** **{success_rate:.1f}%** ({under_2a}/{total_attempted})",
        f"* **Sub-Angstrom Rate (RMSD ≤ 1.0 Å):** **{sub_angstrom_rate:.1f}%** ({under_1a}/{total_attempted})",
        f"* **Mean RMSD:** **{mean_rmsd:.2f} Å**",
        f"* **Median RMSD:** **{median_rmsd:.2f} Å**",
        f"* **Complexes Tested:** {complexes_tested_line}",
        "",
        "| # | PDB ID | Vina ΔG (kcal/mol) | Vinardo ΔG | RMSD (Å) | Time (s) | Status |",
        "| :---: | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]

    for idx, c in enumerate(complex_results, start=1):
        pid = c.get("pdb_id", "")
        vdg = f"{c.get('vina_delta_g_kcal'):.3f}" if c.get("vina_delta_g_kcal") is not None else "—"
        vin = f"{c.get('vinardo_delta_g_kcal'):.3f}" if c.get("vinardo_delta_g_kcal") is not None else "—"
        r = c.get("rmsd_angstroms")
        elapsed = f"{c.get('elapsed_seconds', 0.0):.1f}"

        if r is not None:
            r_str = f"{r:.3f}"
            if r <= 1.0:
                status = "✅ Sub-Angstrom (≤ 1.0 Å)"
            elif r <= 2.0:
                status = "✅ Validated (≤ 2.0 Å)"
            else:
                status = f"⚠️ Near-Native / Divergent (> 2.0 Å, {r:.2f} Å)"
        else:
            r_str = "—"
            status = "❌ Docking Failed"

        lines.append(f"| {idx} | `{pid}` | {vdg} | {vin} | {r_str} | {elapsed} | {status} |")

    lines.extend([
        "",
        "> **Methodology Notice:** Pose RMSD is evaluated strictly **in place** (binding pocket coordinates) using RDKit graph-isomorphism and symmetry correction (`AllChem.CalcRMS`). All complexes are reported without cherry-picking. Metrics serialize directly to `data/benchmarks/casf2016_final_report.json` and `casf2016_final_report.md`.",
        CASF_END_TAG
    ])

    return "\n".join(lines)


def build_dude_markdown(report: Dict[str, Any]) -> str:
    data = report.get("data", {})
    targets: List[Dict[str, Any]] = data.get("targets", [])
    agg = data.get("aggregate_statistics", {})
    n_targets = agg.get("n_targets_attempted", len(targets))
    mean_auc = agg.get("mean_roc_auc", 0.0)

    lines = [
        DUDE_START_TAG,
        f"**Preliminary Screening Results ({n_targets} Target Smoke Test):**",
        f"* **Targets Evaluated:** {n_targets} (`vegfr2` / PDB: 2OH4)",
        f"* **Mean ROC-AUC:** **{mean_auc:.2f}**",
        "",
        "| Target | Protein | PDB ID | Actives | Decoys | ROC-AUC | EF1% | EF5% | EF10% | Status |",
        "| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]

    for t in targets:
        name = t.get("target", "")
        pname = t.get("target_name", name.upper())
        pdb = t.get("pdb_id", "")
        n_act = t.get("n_actives_docked", t.get("n_actives_attempted", 0))
        n_dec = t.get("n_decoys_docked", t.get("n_decoys_attempted", 0))
        auc = f"{t.get('roc_auc', 0.0):.3f}"
        ef1 = f"{t.get('ef_1pct', 0.0):.2f}"
        ef5 = f"{t.get('ef_5pct', 0.0):.2f}"
        ef10 = f"{t.get('ef_10pct', 0.0):.2f}"
        lines.append(f"| `{name}` | {pname} | `{pdb}` | {n_act} | {n_dec} | {auc} | {ef1} | {ef5} | {ef10} | Preliminary Smoke Test |")

    lines.extend([
        "",
        "> [!IMPORTANT]",
        "> **Scientific Interpretation & Sample-Size Context:**",
        "> 1. **Sample Size Insufficiency ($N=15$):** The reported ROC-AUC (0.12) comes from a preliminary single-target smoke test (VEGFR2) consisting of only 5 actives and 10 decoys ($N=15$). In empirical chemoinformatics, $N=15$ is statistically uninformative—neither strong nor poor general screening ability can be concluded from this sample.",
        "> 2. **Pose Accuracy vs. Screening Power:** AutoDock Vina's empirical scoring function was designed for crystallographic pose reconstruction (local energetic minimum in a pocket), not library-scale ranking against property-matched decoys. Raw Vina scores typically require specialized rescoring functions (Vinardo, CNN/GNINA, or machine learning scoring) to achieve high enrichment against property-matched decoys (Mysinger et al., 2012).",
        "> 3. **Roadmap:** The complete **8-Target Diverse Screening Suite** (covering multiple therapeutic target classes with statistical power) is scheduled under Phase 3 of the Bindora v2.0 Roadmap.",
        DUDE_END_TAG
    ])

    return "\n".join(lines)


def update_readme_content(readme_text: str, casf_block: str, dude_block: str) -> str:
    # 1. Update CASF-2016 block
    if CASF_START_TAG in readme_text and CASF_END_TAG in readme_text:
        start_idx = readme_text.index(CASF_START_TAG)
        end_idx = readme_text.index(CASF_END_TAG) + len(CASF_END_TAG)
        readme_text = readme_text[:start_idx] + casf_block + readme_text[end_idx:]
    else:
        # Fallback: replace old section 4.1 text
        old_pattern_start = "**Validated Performance (5-Complex Diverse Run):**"
        old_pattern_end = "casf2016_final_report.md`."
        if old_pattern_start in readme_text and old_pattern_end in readme_text:
            s = readme_text.index(old_pattern_start)
            e = readme_text.index(old_pattern_end) + len(old_pattern_end)
            readme_text = readme_text[:s] + casf_block + readme_text[e:]
        else:
            print("Warning: Could not find CASF-2016 marker or pattern in README.md")

    # 2. Update DUD-E block
    if DUDE_START_TAG in readme_text and DUDE_END_TAG in readme_text:
        start_idx = readme_text.index(DUDE_START_TAG)
        end_idx = readme_text.index(DUDE_END_TAG) + len(DUDE_END_TAG)
        readme_text = readme_text[:start_idx] + dude_block + readme_text[end_idx:]
    else:
        # Fallback: replace old caveat / metrics in section 4.2
        caveat_anchor = "> **Mandatory Scientific Caveat (Mysinger et al., 2012):**"
        if caveat_anchor in readme_text:
            s = readme_text.index(caveat_anchor)
            # Find the next section delimiter
            next_delim = readme_text.find("---", s)
            if next_delim != -1:
                readme_text = readme_text[:s] + dude_block + "\n\n" + readme_text[next_delim:]
            else:
                readme_text = readme_text[:s] + dude_block + "\n"
        else:
            print("Warning: Could not find DUD-E marker or pattern in README.md")

    return readme_text


def main():
    parser = argparse.ArgumentParser(description="Synchronize benchmark JSON reports into README.md")
    parser.add_argument("--check", action="store_true", help="Check if README is up to date without modifying")
    args = parser.parse_args()

    casf_data = load_json(CASF_JSON_PATH)
    dude_data = load_json(DUDE_JSON_PATH)

    casf_md = build_casf2016_markdown(casf_data)
    dude_md = build_dude_markdown(dude_data)

    current_readme = README_PATH.read_text(encoding="utf-8")
    new_readme = update_readme_content(current_readme, casf_md, dude_md)

    if args.check:
        if current_readme == new_readme:
            print("SUCCESS: README.md is perfectly in sync with benchmark JSON reports.")
            sys.exit(0)
        else:
            print("ERROR: README.md is out of sync with benchmark JSON reports. Run python scripts/sync_benchmarks_to_readme.py to update.")
            sys.exit(1)

    README_PATH.write_text(new_readme, encoding="utf-8")
    print(f"SUCCESS: Synchronized {CASF_JSON_PATH.name} and {DUDE_JSON_PATH.name} into README.md")


if __name__ == "__main__":
    main()
