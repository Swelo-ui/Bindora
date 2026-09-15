import math
import json
import urllib.request
import urllib.parse
from typing import Dict, Any, List, Optional
from backend.config import CACHE_DIR, CHEMBL_BASE_URL

HEADERS = {"User-Agent": "AnuDock-Research-Tool/1.0 (academic; +https://github.com/AnuDock)"}

def _get_json(url: str, timeout: int = 15) -> Optional[Dict[str, Any]]:
    try:
        req = urllib.request.Request(url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception:
        return None

class BioactivityService:
    """Service for thermodynamic conversion of docking affinities and cross-checking against ChEMBL/BindingDB."""

    @staticmethod
    def calculate_thermodynamics(affinity_kcal: float, heavy_atoms: int, mw_da: float) -> Dict[str, Any]:
        """Convert binding free energy (delta G) to theoretical Kd, pKd, and Ligand Efficiency."""
        # delta G = R * T * ln(Kd)
        # R = 1.9872036 cal/(mol*K), T = 298.15 K => RT = 592.4847 cal/mol = 0.5924847 kcal/mol
        rt = 0.5924847
        kd_molar = math.exp(affinity_kcal / rt)
        kd_nm = kd_molar * 1e9
        kd_um = kd_molar * 1e6
        pkd = -math.log10(kd_molar) if kd_molar > 0 else 0.0

        # Ligand Efficiency (LE) = -delta G / heavy_atoms
        le = round(-affinity_kcal / heavy_atoms, 3) if heavy_atoms > 0 else 0.0

        # Binding Efficiency Index (BEI) = pKd / (MW in kDa)
        mw_kda = mw_da / 1000.0 if mw_da > 0 else 1.0
        bei = round(pkd / mw_kda, 2)

        # Size-Independent Ligand Efficiency (SILE) = -delta G / (heavy_atoms ** 0.3)
        sile = round(-affinity_kcal / (heavy_atoms ** 0.3), 3) if heavy_atoms > 0 else 0.0

        # Fit Quality (FQ) = LE / LE_max, where LE_max ~ 0.073 + 5.14 / heavy_atoms (Reynolds et al.)
        le_max = 0.073 + (5.14 / heavy_atoms) if heavy_atoms > 0 else 1.0
        fit_quality = round(le / le_max, 3) if le_max > 0 else 0.0

        # Classification of predicted affinity & Weak-Binder Alerting
        is_weak = affinity_kcal > -6.0
        if affinity_kcal <= -9.0:
            potency_class = "Sub-micromolar to Nanomolar (High Affinity Candidate)"
            weak_warning = None
        elif affinity_kcal <= -7.0:
            potency_class = "Low Micromolar (Moderate Affinity Hit)"
            weak_warning = None
        elif affinity_kcal <= -6.0:
            potency_class = "High Micromolar (Weak Screening Hit)"
            weak_warning = None
        else:
            potency_class = "Sub-threshold / Marginal Binding (High Micromolar/Millimolar Kd)"
            weak_warning = (
                f"Sub-threshold / Weak Binding Alert (ΔG = {round(affinity_kcal, 2)} kcal/mol, Kd ≈ {round(kd_um, 1)} µM): "
                "Calculated affinity falls above the -6.0 kcal/mol threshold. Such weak interactions typically reflect "
                "superficial surface adhesion or numerical artifacts rather than biologically meaningful active-site inhibition. "
                "Treat docking pose strictly as hypothesis-generating and interpret with extreme caution."
            )

        return {
            "binding_affinity_kcal": round(affinity_kcal, 2),
            "theoretical_kd_nm": round(kd_nm, 2) if kd_nm < 1e6 else round(kd_nm, 0),
            "theoretical_kd_um": round(kd_um, 3),
            "pkd": round(pkd, 2),
            "is_weak_binder": is_weak,
            "weak_binder_warning": weak_warning,
            "ligand_efficiency": {
                "value": le,
                "unit": "kcal/mol/heavy atom",
                "quality": "Optimal (> 0.30)" if le >= 0.3 else "Suboptimal (< 0.30)",
                "citation": "Hopkins et al., Drug Discov. Today 2004"
            },
            "size_independent_le": {
                "value": sile,
                "citation": "Nissink, ChemMedChem 2009"
            },
            "fit_quality": {
                "value": fit_quality,
                "quality": "Optimal (>= 0.8)" if fit_quality >= 0.8 else "Suboptimal (< 0.8)",
                "citation": "Reynolds et al., J. Med. Chem. 2008"
            },
            "binding_efficiency_index": {
                "value": bei,
                "unit": "pKd / kDa",
                "citation": "Abad-Zapatero & Metz, Drug Discov. Today 2005"
            },
            "potency_class": potency_class
        }

    @staticmethod
    def crosscheck_chembl(drug_name: str, target_name: str) -> Dict[str, Any]:
        """Search ChEMBL for curated experimental Ki, IC50, Kd or EC50 measurements for drug-target pair."""
        cache_key = f"chembl_xcheck_{urllib.parse.quote_plus(drug_name.lower())}_{urllib.parse.quote_plus(target_name.lower())}.json"
        cache_file = CACHE_DIR / cache_key
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                pass

        # 1. Search compound in ChEMBL
        mol_url = f"{CHEMBL_BASE_URL}/molecule/search?q={urllib.parse.quote(drug_name)}&format=json"
        mol_data = _get_json(mol_url)
        mol_chembl_id = None
        mol_pref_name = drug_name
        if mol_data and mol_data.get("molecules"):
            mol_obj = mol_data["molecules"][0]
            mol_chembl_id = mol_obj.get("molecule_chembl_id")
            mol_pref_name = mol_obj.get("pref_name") or drug_name

        # 2. Search target in ChEMBL
        target_url = f"{CHEMBL_BASE_URL}/target/search?q={urllib.parse.quote(target_name)}&format=json"
        target_data = _get_json(target_url)
        target_chembl_id = None
        target_pref_name = target_name
        target_organism = "Unknown"
        if target_data and target_data.get("targets"):
            # Prefer single protein targets
            single_targets = [t for t in target_data["targets"] if t.get("target_type") == "SINGLE PROTEIN"]
            target_obj = single_targets[0] if single_targets else target_data["targets"][0]
            target_chembl_id = target_obj.get("target_chembl_id")
            target_pref_name = target_obj.get("pref_name") or target_name
            target_organism = target_obj.get("organism") or "Homo sapiens"

        experimental_activities = []
        if mol_chembl_id and target_chembl_id:
            act_url = (
                f"{CHEMBL_BASE_URL}/activity?molecule_chembl_id={mol_chembl_id}&"
                f"target_chembl_id={target_chembl_id}&limit=10&format=json"
            )
            act_data = _get_json(act_url)
            if act_data and act_data.get("activities"):
                for act in act_data["activities"]:
                    std_type = act.get("standard_type")
                    std_val = act.get("standard_value")
                    std_units = act.get("standard_units")
                    if std_type in ("IC50", "Ki", "Kd", "EC50") and std_val is not None:
                        experimental_activities.append({
                            "type": std_type,
                            "value": float(std_val),
                            "units": std_units or "nM",
                            "relation": act.get("standard_relation", "="),
                            "assay_description": act.get("assay_description", "Bioactivity assay"),
                            "assay_type": act.get("assay_type", "B"),
                            "pubmed_id": act.get("document_chembl_id", ""),
                            "source": "ChEMBL Curated Database"
                        })

        has_experimental = len(experimental_activities) > 0
        status_label = "Experimentally Corroborated" if has_experimental else "Computational Prediction Only"
        status_color = "emerald" if has_experimental else "amber"
        summary_note = (
            f"Cross-checked with ChEMBL. Found {len(experimental_activities)} experimental record(s) "
            f"for {mol_pref_name} against {target_pref_name} ({target_organism})."
            if has_experimental else
            f"No direct wet-lab assay records found in ChEMBL for {mol_pref_name} + {target_name}. "
            f"All values are computational estimations and must be treated as hypothesis-generating."
        )

        result = {
            "is_cross_checked": has_experimental,
            "status_badge": status_label,
            "status_color": status_color,
            "summary_note": summary_note,
            "drug_searched": drug_name,
            "target_searched": target_name,
            "molecule_chembl_id": mol_chembl_id,
            "target_chembl_id": target_chembl_id,
            "target_organism": target_organism,
            "experimental_records": experimental_activities,
            "chembl_target_url": f"https://www.ebi.ac.uk/chembl/target_report_card/{target_chembl_id}/" if target_chembl_id else None,
            "chembl_compound_url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{mol_chembl_id}/" if mol_chembl_id else None
        }

        # Cache result
        try:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2)
        except Exception:
            pass

        return result
