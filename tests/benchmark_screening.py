#!/usr/bin/env python3
"""
tests/benchmark_screening.py
Bindora DUD-E Virtual Screening Benchmark (v2.0)

Evaluates virtual-screening enrichment power using DUD-E.
Metrics: ROC-AUC, EF1%, EF5%, EF10%.

MANDATORY CAVEAT (always reported with every DUD-E result):
  DUD-E decoys are property-matched but NOT topologically diversified
  from actives. Analogue bias can inflate enrichment scores.
  Mysinger et al. J. Med. Chem. 2012, 55(14), 6582.

Data: http://dude.docking.org/  (UCSF, free, no registration)

Usage:
  python tests/benchmark_screening.py --subset diverse
  python tests/benchmark_screening.py --target ache --resume
  python tests/benchmark_screening.py --subset full --exhaustiveness 4 --resume
"""

import os
import sys
import ssl
import math
import json
import time
import logging
import argparse
import traceback
import urllib.request
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Dict, Any, Optional, Set

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from backend.config import DATA_DIR, BENCHMARKS_DIR
from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine
from backend.utils.report_emitter import (
    emit_screening_report, emit_casf2016_checkpoint,
    atomic_json_dump, ScientificIntegrityError,
)

if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# Paths
DUDE_CACHE_DIR        = DATA_DIR / "dude"
SCREENING_CHECKPOINT  = BENCHMARKS_DIR / "dude_progress.json"
SCREENING_REPORT_JSON = BENCHMARKS_DIR / "dude_screening_report.json"
SCREENING_REPORT_MD   = BENCHMARKS_DIR / "dude_screening_report.md"
SCREENING_ERROR_LOG   = BENCHMARKS_DIR / "dude_errors.log"
DUDE_BASE_URL = "http://dude.docking.org/targets"

# Mandatory DUD-E analogue-bias caveat -- NEVER omitted from any report
DUDE_ANALOGUE_BIAS_CAVEAT = (
    "DUD-E decoys are property-matched (MW, cLogP, HBA, HBD, rotatable bonds) "
    "but NOT topologically diversified from actives. Some decoys carry structural "
    "similarity to actives (analogue bias), inflating enrichment scores for some "
    "methods. This caveat applies to all published DUD-E benchmarks "
    "(Mysinger et al. 2012). Results must be interpreted alongside "
    "pose-accuracy data (CASF-2016 redocking), not in isolation."
)

# DUD-E target definitions
DUDE_DIVERSE_8 = [
    dict(target="ace",    name="Angiotensin Converting Enzyme",  pdb_id="1O86", cls="Metalloprotease"),
    dict(target="ache",   name="Acetylcholinesterase",           pdb_id="1EVE", cls="Hydrolase"),
    dict(target="gcr",    name="Glucocorticoid Receptor",        pdb_id="1M2Z", cls="Nuclear Receptor"),
    dict(target="hivpr",  name="HIV-1 Protease",                 pdb_id="1HSG", cls="Aspartyl Protease"),
    dict(target="hivrt",  name="HIV-1 Reverse Transcriptase",    pdb_id="1RT1", cls="Polymerase"),
    dict(target="kif11",  name="Kinesin Eg5",                    pdb_id="1II6", cls="Motor Protein"),
    dict(target="src",    name="Src Tyrosine Kinase",            pdb_id="2SRC", cls="Kinase"),
    dict(target="vegfr2", name="VEGFR2",                         pdb_id="2OH4", cls="RTK"),
]

DUDE_EXTENDED = DUDE_DIVERSE_8 + [
    dict(target="abl1",  name="Abl Tyrosine Kinase",   pdb_id="2HYY", cls="Kinase"),
    dict(target="bace1", name="Beta-Secretase 1",      pdb_id="2QP8", cls="Aspartyl Protease"),
    dict(target="cah2",  name="Carbonic Anhydrase 2",  pdb_id="1OKL", cls="Lyase"),
    dict(target="cdk2",  name="CDK2",                  pdb_id="1H00", cls="Kinase"),
    dict(target="dpp4",  name="DPP-4",                 pdb_id="1X70", cls="Protease"),
    dict(target="esr1",  name="Estrogen Receptor A",   pdb_id="3ERT", cls="Nuclear Receptor"),
    dict(target="fa10",  name="Coagulation Factor X",  pdb_id="1FJS", cls="Serine Protease"),
    dict(target="hsp90", name="HSP90",                 pdb_id="1UYG", cls="Chaperone"),
    dict(target="jak2",  name="JAK2",                  pdb_id="2XA4", cls="Kinase"),
    dict(target="parp1", name="PARP-1",                pdb_id="3L3M", cls="Transferase"),
    dict(target="pgh2",  name="COX-2",                 pdb_id="3LN1", cls="Oxidoreductase"),
    dict(target="pparg", name="PPAR Gamma",            pdb_id="2GTK", cls="Nuclear Receptor"),
    dict(target="reni",  name="Human Renin",           pdb_id="2IKO", cls="Aspartyl Protease"),
    dict(target="thb",   name="Thrombin",              pdb_id="1VZQ", cls="Serine Protease"),
    dict(target="try1",  name="Trypsin",               pdb_id="2FTL", cls="Serine Protease"),
    dict(target="urok",  name="Urokinase (uPA)",       pdb_id="1GJA", cls="Serine Protease"),
    dict(target="wee1",  name="Wee1 Kinase",           pdb_id="3CQU", cls="Kinase"),
]

SUBSETS = {
    "diverse":  DUDE_DIVERSE_8,
    "extended": DUDE_EXTENDED,
    "full":     DUDE_EXTENDED,
}


def setup_logger():
    logger = logging.getLogger("dude_screening")
    logger.setLevel(logging.DEBUG)
    if not logger.handlers:
        fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        ch.setFormatter(fmt)
        logger.addHandler(ch)
        BENCHMARKS_DIR.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(SCREENING_ERROR_LOG), encoding="utf-8", mode="a")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def roc_auc(labels: List[int], scores: List[float]) -> float:
    """Compute ROC-AUC via trapezoid rule. Higher score = more likely active."""
    if not labels or sum(labels) == 0 or sum(labels) == len(labels):
        return 0.5
    paired = sorted(zip(scores, labels), key=lambda x: x[0], reverse=True)
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    tp, area = 0, 0.0
    for _, lbl in paired:
        if lbl == 1:
            tp += 1
        else:
            area += tp
    return round(area / (n_pos * n_neg), 4) if n_pos * n_neg > 0 else 0.5


def enrichment_factor(labels: List[int], scores: List[float], frac: float) -> float:
    """EF at top (frac*100)% of ranked library. Random baseline = 1.0."""
    if not labels:
        return 0.0
    n_active = sum(labels)
    if n_active == 0:
        return 0.0
    n_top = max(1, int(len(labels) * frac))
    paired = sorted(zip(scores, labels), key=lambda x: x[0], reverse=True)
    hits = sum(lbl for _, lbl in paired[:n_top])
    expected = n_active * frac
    return round(hits / expected, 3) if expected > 0 else 0.0


# ---------------------------------------------------------------------------
# ChEMBL target IDs for DUD-E diverse-8 + extended targets
# Used as fallback when dude.docking.org is unreachable
# ---------------------------------------------------------------------------
_CHEMBL_TARGET_MAP: Dict[str, str] = {
    "ace":    "CHEMBL1808",   # Angiotensin Converting Enzyme
    "ache":   "CHEMBL220",    # Acetylcholinesterase
    "gcr":    "CHEMBL2034",   # Glucocorticoid Receptor
    "hivpr":  "CHEMBL247",    # HIV-1 Protease
    "hivrt":  "CHEMBL4153",   # HIV-1 Reverse Transcriptase
    "kif11":  "CHEMBL1836",   # Kinesin Eg5
    "src":    "CHEMBL267",    # Src Tyrosine Kinase
    "vegfr2": "CHEMBL279",    # VEGFR2
    "abl1":   "CHEMBL1862",   # Abl Tyrosine Kinase
    "bace1":  "CHEMBL4822",   # Beta-Secretase 1
    "cah2":   "CHEMBL205",    # Carbonic Anhydrase 2
    "cdk2":   "CHEMBL301",    # CDK2
    "dpp4":   "CHEMBL284",    # DPP-4
    "esr1":   "CHEMBL206",    # Estrogen Receptor A
    "fa10":   "CHEMBL244",    # Coagulation Factor X
    "hsp90":  "CHEMBL3880",   # HSP90
    "jak2":   "CHEMBL2971",   # JAK2
    "parp1":  "CHEMBL3105",   # PARP-1
    "pgh2":   "CHEMBL230",    # COX-2
    "pparg":  "CHEMBL235",    # PPAR Gamma
    "reni":   "CHEMBL286",    # Human Renin
    "thb":    "CHEMBL204",    # Thrombin
    "try1":   "CHEMBL209",    # Trypsin
    "urok":   "CHEMBL3286",   # Urokinase
    "wee1":   "CHEMBL3955",   # Wee1 Kinase
}


def _ssl_no_verify_opener() -> urllib.request.OpenerDirector:
    """Return a URL opener that skips SSL certificate verification."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    https_handler = urllib.request.HTTPSHandler(context=ctx)
    return urllib.request.build_opener(https_handler)


def _fetch_from_chembl(
    target: str,
    kind: str,
    n: int = 500,
    logger: Optional[logging.Logger] = None,
) -> Optional[List[str]]:
    """
    Fetch real published bioactivity SMILES from ChEMBL REST API.

    Actives  (kind contains 'actives'):
        IC50 <= 1000 nM, sorted most-potent-first.
    Inactives (kind contains 'decoys'):
        IC50 >= 50000 nM (>50 µM), representing weak/non-binders.

    Data source: ChEMBL database (https://www.ebi.ac.uk/chembl/)
    Citation: Mendez D et al. ChEMBL: towards direct deposition of bioassay
              data. Nucleic Acids Res 2019, 47(D1):D930-D940.

    Returns deduplicated list of canonical SMILES sorted by potency,
    or None if ChEMBL is unreachable or target not mapped.
    """
    chembl_id = _CHEMBL_TARGET_MAP.get(target.lower())
    if not chembl_id:
        if logger:
            logger.warning("[%s] No ChEMBL ID mapped — cannot use ChEMBL fallback", target)
        return None
    try:
        is_active = "actives" in kind
        if is_active:
            # Potent: IC50 <= 1000 nM, units explicitly nM to prevent µM/µg/mL contamination
            filt = (
                "standard_type=IC50"
                "&standard_relation=%3D"
                "&standard_value__lte=1000"
                "&standard_units=nM"
                "&order_by=standard_value"
            )
        else:
            # Inactives: IC50 >= 50000 nM (weak binders / non-binders), units explicitly nM
            filt = (
                "standard_type=IC50"
                "&standard_relation=%3D"
                "&standard_value__gte=50000"
                "&standard_units=nM"
            )
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/activity?"
            + filt
            + "&target_chembl_id=" + chembl_id
            + "&assay_type=B"            # Binding assays only
            + "&limit=" + str(n)
            + "&format=json"
        )
        opener = _ssl_no_verify_opener()
        if logger:
            logger.info(
                "[%s] Fetching from ChEMBL (%s, %s) — %s",
                target, chembl_id, "actives" if is_active else "inactives", url,
            )
        with opener.open(url, timeout=45) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))

        # Deduplicate by canonical SMILES; preserve potency order
        seen: Set[str] = set()
        smiles_list: List[str] = []
        for act in data.get("activities", []):
            smi = act.get("canonical_smiles")
            if isinstance(smi, str) and smi.strip() and smi not in seen:
                seen.add(smi)
                smiles_list.append(smi.strip())

        total = data.get("page_meta", {}).get("total_count", "?")
        if logger:
            logger.info(
                "[%s] ChEMBL returned %d unique SMILES (total available: %s) "
                "[source: ChEMBL %s, IC50 %s]",
                target, len(smiles_list), total, chembl_id,
                "<= 1 µM" if is_active else ">= 50 µM",
            )
        return smiles_list if smiles_list else None
    except Exception as exc:
        if logger:
            logger.warning("[%s] ChEMBL fetch failed: %s", target, exc)
        return None


def fetch_dude_smiles(
    target: str,
    kind: str,
    logger: Optional[logging.Logger] = None,
) -> Optional[List[str]]:
    """
    Fetch actives_final.ism or decoys_final.ism for a DUD-E target.

    Data-source priority (real published data only — no hardcoding):

      Layer 1 — Local disk cache (data/dude/<target>_<kind>)
                 Written by any successful layer; persists across runs.

      Layer 2 — DUD-E server (dude.docking.org) with SSL-bypass
                 Original source; SSL cert may be expired but data is authentic.
                 Multiple URL patterns attempted.

      Layer 3 — ChEMBL REST API (www.ebi.ac.uk/chembl)
                 Real published IC50 bioactivity records.
                 Actives  : IC50 <= 1 µM, sorted most-potent-first.
                 Inactives: IC50 >= 50 µM (weak/non-binders as decoy proxy).
                 Citation : Mendez D et al. Nucleic Acids Res 2019, 47(D1):D930.

    Returns list of SMILES strings, or None if ALL real-data layers fail.
    Raises RuntimeError if target exists in ChEMBL map but all APIs are down.
    """
    DUDE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = DUDE_CACHE_DIR / (target + "_" + kind)

    # --- Layer 1: local cache ---
    if cache_path.exists():
        try:
            lines = cache_path.read_text(encoding="utf-8").splitlines()
            smiles = [
                ln.split()[0]
                for ln in lines
                if ln.strip() and not ln.startswith("#")
            ]
            if smiles:
                if logger:
                    logger.info(
                        "[%s] Using cached %s (%d compounds, %s)",
                        target, kind, len(smiles), cache_path,
                    )
                return smiles
        except Exception:
            pass

    # --- Layer 2: DUD-E server (multiple URL patterns, SSL bypass) ---
    dude_urls = [
        DUDE_BASE_URL + "/" + target + "/" + kind,
        "https://dude.docking.org/targets/" + target + "/" + kind,
    ]
    opener = _ssl_no_verify_opener()
    for url in dude_urls:
        try:
            if logger:
                logger.info("[%s] Trying DUD-E: %s", target, url)
            req = urllib.request.Request(url, headers={"User-Agent": "Bindora-Research/2.0"})
            with opener.open(req, timeout=30) as resp:
                content = resp.read().decode("utf-8", errors="replace")
            if content.strip():
                smiles = [
                    ln.split()[0]
                    for ln in content.splitlines()
                    if ln.strip() and not ln.startswith("#")
                ]
                if smiles:
                    # Prepend source metadata as comment line
                    tagged = (
                        "# source=dude.docking.org target=" + target + " kind=" + kind + "\n"
                        + content
                    )
                    cache_path.write_text(tagged, encoding="utf-8")
                    if logger:
                        logger.info(
                            "[%s] DUD-E: fetched %d SMILES from %s",
                            target, len(smiles), url,
                        )
                    return smiles
        except Exception as exc:
            if logger:
                logger.debug("[%s] DUD-E URL failed (%s): %s", target, url, exc)
            continue

    # --- Layer 3: ChEMBL REST API ---
    chembl_smiles = _fetch_from_chembl(target, kind, n=500, logger=logger)
    if chembl_smiles:
        # Cache with provenance header
        chembl_id = _CHEMBL_TARGET_MAP.get(target.lower(), "?")
        header = (
            "# source=ChEMBL target=" + target
            + " chembl_id=" + chembl_id
            + " kind=" + kind
            + " fetched=" + datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            + "\n"
        )
        cache_path.write_text(header + "\n".join(chembl_smiles), encoding="utf-8")
        return chembl_smiles

    # All real-data layers exhausted
    if logger:
        logger.error(
            "[%s] All data sources failed for %s.\n"
            "  - DUD-E server (dude.docking.org): unreachable or 404.\n"
            "  - ChEMBL API: unreachable or no data for %s.\n"
            "  Fix: check internet connection, or pre-download DUD-E files to %s",
            target, kind, _CHEMBL_TARGET_MAP.get(target.lower(), "unknown ChEMBL ID"),
            cache_path,
        )
    return None


def screen_target(
    ti: Dict[str, Any],
    exh: int,
    max_a: Optional[int],
    max_d: Optional[int],
    logger: logging.Logger,
) -> Dict[str, Any]:
    """Dock actives + decoys for one DUD-E target. Never raises."""
    tgt = ti["target"]
    pdb_id = ti["pdb_id"]
    result: Dict[str, Any] = {
        "target": tgt, "target_name": ti["name"],
        "pdb_id": pdb_id, "target_class": ti.get("cls", "Unknown"),
        "success": False,
        "n_actives_attempted": 0, "n_actives_docked": 0,
        "n_decoys_attempted": 0, "n_decoys_docked": 0,
        "roc_auc": None, "ef_1pct": None, "ef_5pct": None, "ef_10pct": None,
        "error": None, "elapsed_seconds": 0.0,
    }
    t0 = time.time()
    try:
        act = fetch_dude_smiles(tgt, "actives_final.ism", logger=logger)
        dec = fetch_dude_smiles(tgt, "decoys_final.ism", logger=logger)
        if not act:
            raise RuntimeError("Failed to fetch actives for DUD-E target " + repr(tgt))
        if not dec:
            raise RuntimeError("Failed to fetch decoys for DUD-E target " + repr(tgt))

        # Record data source provenance in result
        cache_act = DUDE_CACHE_DIR / (tgt + "_actives_final.ism")
        if cache_act.exists():
            first_line = cache_act.read_text(encoding="utf-8").splitlines()[0] if cache_act.stat().st_size > 0 else ""
            result["data_source"] = (
                "ChEMBL" if "source=ChEMBL" in first_line else
                "DUD-E"  if "source=dude"  in first_line else
                "cache"
            )
        else:
            result["data_source"] = "unknown"

        if max_a: act = act[:max_a]
        if max_d: dec = dec[:max_d]
        result["n_actives_attempted"] = len(act)
        result["n_decoys_attempted"] = len(dec)
        logger.info(
            "  [%s] Actives=%d, Decoys=%d  [data_source=%s]",
            tgt, len(act), len(dec), result.get("data_source", "?"),
        )
        # Prepare receptor
        rec_meta = StructureFetcher.fetch_rcsb_pdb(pdb_id)
        if not rec_meta or "pdb_content" not in rec_meta:
            raise RuntimeError("Cannot fetch receptor PDB " + pdb_id)
        rec = DockingEngine.prepare_receptor(rec_meta["pdb_content"])
        pocket = rec.get("detected_pocket") or rec.get("blind_docking_box")
        if not pocket:
            raise RuntimeError("No binding pocket detected for " + pdb_id)
        # Dock all actives + decoys
        scores: List[float] = []
        labels: List[int] = []
        fails = 0
        all_ligs = [(s, 1) for s in act] + [(s, 0) for s in dec]
        for i, (smi, lbl) in enumerate(all_ligs):
            if (i + 1) % 100 == 0:
                logger.info("    [" + tgt + "] " + str(i+1) + "/" + str(len(all_ligs)) + " screened, " + str(fails) + " failures")
            try:
                lig = DockingEngine.prepare_ligand(smi)
                poses = DockingEngine.run_docking(
                    rec["pdbqt_text"], lig["pdbqt_text"],
                    pocket["center"], pocket["size"],
                    exhaustiveness=exh, num_modes=1, seed=42,
                )
                if poses and poses[0].get("affinity_kcal") is not None:
                    scores.append(-poses[0]["affinity_kcal"])
                    labels.append(lbl)
                    if lbl == 1: result["n_actives_docked"] += 1
                    else: result["n_decoys_docked"] += 1
                else:
                    fails += 1
            except Exception:
                fails += 1
        if not scores:
            raise RuntimeError("No ligands successfully docked for " + tgt)
        result["roc_auc"] = roc_auc(labels, scores)
        result["ef_1pct"]  = enrichment_factor(labels, scores, 0.01)
        result["ef_5pct"]  = enrichment_factor(labels, scores, 0.05)
        result["ef_10pct"] = enrichment_factor(labels, scores, 0.10)
        result["success"] = True
    except Exception as exc:
        result["error"] = str(exc)
        result["success"] = False
        logger.warning("[" + tgt + "] FAILED: " + str(exc))
        logger.debug(traceback.format_exc())
    result["elapsed_seconds"] = round(time.time() - t0, 2)
    return result


def gen_md(results: List[Dict], stats: Dict, meta: Dict) -> str:
    """Generate Markdown -- always shows ALL targets, no cherry-picking."""
    lines = [
        "# Bindora DUD-E Virtual Screening Report (v2.0)",
        "",
        "**Date:** " + meta["timestamp_utc"][:10] + " | Subset: " + meta["subset_name"] + " | Exh: " + str(meta["exhaustiveness"]),
        "",
        "> **Mandatory Caveat (always reported alongside DUD-E results):**",
        "> " + DUDE_ANALOGUE_BIAS_CAVEAT,
        "",
        "## Aggregate Statistics",
        "",
        "| Metric | Value |",
        "| :--- | :---: |",
        "| Targets attempted | " + str(stats["n_targets_attempted"]) + " |",
        "| Targets successful | " + str(stats["n_targets_successful"]) + " |",
        "| Mean ROC-AUC | " + str(stats["mean_roc_auc"]) + " |",
        "| Mean EF1% | " + str(stats["mean_ef_1pct"]) + " |",
        "| Mean EF5% | " + str(stats["mean_ef_5pct"]) + " |",
        "| Mean EF10% | " + str(stats["mean_ef_10pct"]) + " |",
        "",
        "## Per-Target Results (Full Distribution)",
        "",
        "| Target | Name | Data Source | Actives | Decoys | ROC-AUC | EF1% | EF5% | EF10% | Status |",
        "| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]
    for r in results:
        au = ("{:.3f}".format(r["roc_auc"])  if r.get("roc_auc")  is not None else "--")
        e1 = ("{:.2f}".format(r["ef_1pct"])  if r.get("ef_1pct")  is not None else "--")
        e5 = ("{:.2f}".format(r["ef_5pct"])  if r.get("ef_5pct")  is not None else "--")
        e10= ("{:.2f}".format(r["ef_10pct"]) if r.get("ef_10pct") is not None else "--")
        ds = r.get("data_source", "DUD-E")
        if not r.get("success"):
            st = "Failed: " + (r.get("error") or "Unknown")[:40]
        elif (r.get("roc_auc") or 0) >= 0.7:
            st = "Good (AUC >= 0.7)"
        else:
            st = "Marginal"
        lines.append(
            "| " + r["target"] + " | " + r["target_name"] + " | " + ds + " | " +
            str(r["n_actives_docked"]) + " | " + str(r["n_decoys_docked"]) +
            " | " + au + " | " + e1 + " | " + e5 + " | " + e10 + " | " + st + " |"
        )
    lines += [
        "",
        "## Data Sources & References",
        "",
        "- **DUD-E**: Mysinger MM et al. J. Med. Chem. 2012, 55(14), 6582. DOI: 10.1021/jm300687e",
        "- **ChEMBL** (fallback when DUD-E server unavailable): Mendez D et al. "
        "Nucleic Acids Res 2019, 47(D1):D930-D940. DOI: 10.1093/nar/gky1075",
        "  - Actives: binding IC50 <= 1 µM (most potent first)",
        "  - Inactives (decoy proxy): binding IC50 >= 50 µM",
    ]

    return "\n".join(lines)


def run_screening(targets, exh, resume, subset_name, max_a=None, max_d=None):
    logger = setup_logger()
    completed: Set[str] = set()
    existing: List[Dict] = []
    if resume and SCREENING_CHECKPOINT.exists():
        try:
            cp = json.loads(SCREENING_CHECKPOINT.read_text(encoding="utf-8"))
            completed = set(cp.get("completed_targets", []))
            existing = cp.get("results", [])
            logger.info("Resume: " + str(len(completed)) + " targets done")
        except Exception:
            pass
    remaining = [t for t in targets if t["target"] not in completed]
    started_at = datetime.now(timezone.utc).isoformat()
    all_results = list(existing)
    print("\nDUD-E: " + subset_name + ", " + str(len(remaining)) + " targets, exh=" + str(exh))
    for idx, ti in enumerate(remaining, 1):
        print("[{:02d}/{:02d}] {} ({}) ...".format(idx, len(remaining), ti["target"], ti["name"]), flush=True)
        r = screen_target(ti, exh, max_a, max_d, logger)
        if r["success"]:
            print("  AUC={:.3f} EF1%={:.2f} EF5%={:.2f} ({:.1f}s)".format(r["roc_auc"], r["ef_1pct"], r["ef_5pct"], r["elapsed_seconds"]))
        else:
            print("  FAILED: " + (r.get("error") or "?")[:70])
        all_results.append(r)
        completed.add(ti["target"])
        try:
            emit_casf2016_checkpoint(
                {"completed_ids": list(completed), "completed_targets": list(completed),
                 "results": all_results, "started_at": started_at, "total_targets": len(targets)},
                SCREENING_CHECKPOINT,
            )
        except Exception as ce:
            logger.error("Checkpoint failed: " + str(ce))
    succ = [r for r in all_results if r.get("success")]
    def _m(vals): return round(sum(vals)/len(vals), 3) if vals else None
    agg = {
        "n_targets_attempted":  len(all_results),
        "n_targets_successful": len(succ),
        "mean_roc_auc":  _m([r["roc_auc"]  for r in succ if r.get("roc_auc")  is not None]),
        "mean_ef_1pct":  _m([r["ef_1pct"]  for r in succ if r.get("ef_1pct")  is not None]),
        "mean_ef_5pct":  _m([r["ef_5pct"]  for r in succ if r.get("ef_5pct")  is not None]),
        "mean_ef_10pct": _m([r["ef_10pct"] for r in succ if r.get("ef_10pct") is not None]),
    }
    meta = {
        "benchmark_version": "2.0", "benchmark_type": "DUD-E_virtual_screening",
        "subset_name": subset_name, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exhaustiveness": exh, "operating_system": sys.platform,
    }
    sd = {
        "targets": all_results, "aggregate_statistics": agg, "subset_name": subset_name,
        "dude_caveat": DUDE_ANALOGUE_BIAS_CAVEAT, "metadata": meta,
    }
    try:
        emit_screening_report("DUD-E Virtual Screening Benchmark", sd, SCREENING_REPORT_JSON)
        print("Report: " + str(SCREENING_REPORT_JSON))
    except ScientificIntegrityError as sie:
        print("INTEGRITY ERROR: " + str(sie), file=sys.stderr)
        raise
    except Exception as e:
        logger.error("emit failed: " + str(e))
        atomic_json_dump(sd, SCREENING_REPORT_JSON)
    try:
        SCREENING_REPORT_MD.write_text(gen_md(all_results, agg, meta), encoding="utf-8")
    except Exception as e:
        logger.error("MD failed: " + str(e))
    # Auto-synchronize README.md with newly generated report
    try:
        sync_script = PROJECT_ROOT / "scripts" / "sync_benchmarks_to_readme.py"
        if sync_script.exists():
            import subprocess
            subprocess.run([sys.executable, str(sync_script)], capture_output=True, text=True)
    except Exception:
        pass
    print("Done. AUC=" + str(agg["mean_roc_auc"]) + " EF1%=" + str(agg["mean_ef_1pct"]))
    return sd


def main():
    p = argparse.ArgumentParser(description="Bindora DUD-E Screening Benchmark v2.0")
    p.add_argument("--subset", choices=list(SUBSETS.keys()), default="diverse")
    p.add_argument("--target", type=str, default=None)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--exhaustiveness", type=int, default=4)
    p.add_argument("--max-actives", type=int, default=None)
    p.add_argument("--max-decoys",  type=int, default=None)
    args = p.parse_args()
    if args.target:
        matched = [t for t in DUDE_EXTENDED if t["target"].lower() == args.target.lower()]
        if not matched:
            print("Unknown target: " + repr(args.target), file=sys.stderr)
            sys.exit(1)
        targets, sn = matched, "single_" + args.target
    else:
        targets = SUBSETS[args.subset]
        sn = args.subset + "_" + str(len(targets))
    run_screening(targets, args.exhaustiveness, args.resume, sn, args.max_actives, args.max_decoys)


if __name__ == "__main__":
    main()
