import math
from typing import Dict, Any, List, Optional
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, FilterCatalog

from backend.utils import sascorer
from pathlib import Path
import json

# Initialize PAINS, Brenk, NIH, and ZINC structural alert catalogs once for efficiency
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

_brenk_params = FilterCatalog.FilterCatalogParams()
_brenk_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
_brenk_catalog = FilterCatalog.FilterCatalog(_brenk_params)

_nih_params = FilterCatalog.FilterCatalogParams()
_nih_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
_nih_catalog = FilterCatalog.FilterCatalog(_nih_params)

_zinc_params = FilterCatalog.FilterCatalogParams()
_zinc_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC)
_zinc_catalog = FilterCatalog.FilterCatalog(_zinc_params)

# Load BOILED-Egg 101-point polygon coordinates from Daina & Zoete (2016) ChemMedChem SI
_BOILED_EGG_FILE = Path(__file__).resolve().parent / "boiled_egg_coords.json"
if _BOILED_EGG_FILE.exists():
    with open(_BOILED_EGG_FILE, "r", encoding="utf-8") as _f:
        _egg_data = json.load(_f)
        _GIA_COORDS = _egg_data.get("gia_coords", [])
        _BBB_COORDS = _egg_data.get("bbb_coords", [])
else:
    _GIA_COORDS = []
    _BBB_COORDS = []

def _point_in_polygon(x: float, y: float, poly: List[List[float]]) -> bool:
    """Ray-casting algorithm to determine if a point (x, y) is inside a polygon."""
    if not poly:
        return False
    n = len(poly)
    inside = False
    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside

# Common CYP structural alert SMARTS
_CYP_ALERTS = {
    "CYP1A2_inhibition": {
        "smarts": "[#6]1:[#6]:[#6]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:12", # Polyaromatic planar ring
        "name": "Planar polycyclic aromatic system (associated with CYP1A2 binding)"
    },
    "CYP2D6_inhibition": {
        "smarts": "[NX3;H2,H1,H0;!$(NC=O)]~[#6]~[#6]~[#6]~c1ccccc1", # Basic nitrogen approx 5-7Å from aromatic ring
        "name": "Basic nitrogen with hydrophobic aryl feature (CYP2D6 pharmacophore)"
    },
    "CYP3A4_inhibition": {
        "smarts": "c1cncn1", # Imidazole / azole ring
        "name": "Azole ring nitrogen (coordinates directly with heme iron in CYP3A4)"
    }
}

class ADMEProfiler:
    """Service for deterministic, calculated physicochemical and pharmacokinetic (ADME) profiling."""

    @staticmethod
    def calculate_adme(smiles_or_mol: Any) -> Dict[str, Any]:
        """Compute comprehensive ADME descriptors, drug-likeness rules, and toxicity alerts."""
        if isinstance(smiles_or_mol, str):
            mol = Chem.MolFromSmiles(smiles_or_mol)
        else:
            mol = smiles_or_mol

        if not mol:
            return {"error": "Invalid chemical structure or SMILES string"}

        # Basic physicochemical properties
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rotb = int(Lipinski.NumRotatableBonds(mol))
        tpsa = round(Descriptors.TPSA(mol), 2)
        mr = round(Crippen.MolMR(mol), 2)
        heavy_atoms = int(mol.GetNumHeavyAtoms())
        aromatic_rings = int(Lipinski.NumAromaticRings(mol))
        fsp3 = round(Descriptors.FractionCSP3(mol), 2)

        # 1. Lipinski's Rule of 5 (Pfizer criteria for oral bioavailability)
        lipinski_violations = []
        if mw > 500:
            lipinski_violations.append(f"MW > 500 Da ({mw})")
        if logp > 5.0:
            lipinski_violations.append(f"LogP > 5.0 ({logp})")
        if hbd > 5:
            lipinski_violations.append(f"H-Bond Donors > 5 ({hbd})")
        if hba > 10:
            lipinski_violations.append(f"H-Bond Acceptors > 10 ({hba})")

        lipinski_status = "Pass" if len(lipinski_violations) <= 1 else ("Borderline" if len(lipinski_violations) == 2 else "Fail")

        # 2. Veber's Rule (GSK criteria for oral bioavailability)
        veber_violations = []
        if rotb > 10:
            veber_violations.append(f"Rotatable bonds > 10 ({rotb})")
        if tpsa > 140.0:
            veber_violations.append(f"TPSA > 140 Å² ({tpsa})")
        veber_pass = len(veber_violations) == 0

        # 3. Ghose Filter
        ghose_violations = []
        if not (160 <= mw <= 480):
            ghose_violations.append(f"MW outside 160-480 Da ({mw})")
        if not (-0.4 <= logp <= 5.6):
            ghose_violations.append(f"LogP outside -0.4 to 5.6 ({logp})")
        if not (40 <= mr <= 130):
            ghose_violations.append(f"MR outside 40-130 ({mr})")
        total_atoms = mol.GetNumAtoms()
        if not (20 <= total_atoms <= 70):
            ghose_violations.append(f"Total atoms outside 20-70 ({total_atoms})")
        ghose_pass = len(ghose_violations) == 0

        # 4. Pharmacokinetic estimations (GI Absorption & BBB Permeation via published BOILED-Egg model)
        # Exact 101-point ellipse coordinates from Daina & Zoete, ChemMedChem 2016
        tpsa_sandp = Descriptors.TPSA(mol, includeSandP=True)
        wlogp = round(Descriptors.MolLogP(mol), 2)
        gi_high = _point_in_polygon(tpsa_sandp, wlogp, _GIA_COORDS)
        gi_absorption = "High" if gi_high else "Moderate / Low"

        bbb_permeant = _point_in_polygon(tpsa_sandp, wlogp, _BBB_COORDS)
        bbb_status = "Permeant (Likely crosses BBB)" if bbb_permeant else "Non-permeant (Low CNS penetration likelihood)"

        # 5. Synthetic Accessibility Score (SAScore: 1.0 easy to 10.0 difficult)
        sa_score_raw = sascorer.calculateScore(mol)
        sa_score = round(sa_score_raw, 2) if sa_score_raw is not None else None
        if sa_score is not None:
            if sa_score < 3.0:
                sa_desc = "Easy synthetic accessibility"
            elif sa_score <= 6.0:
                sa_desc = "Moderate synthetic accessibility"
            else:
                sa_desc = "Difficult / complex synthesis"
        else:
            sa_desc = "N/A"

        # Plasma Protein Binding (PPB) heuristic
        if logp > 3.5:
            ppb_tier = "High (> 90%)"
            ppb_note = "Extensive binding to albumin / alpha-1-acid glycoprotein driven by lipophilicity."
        elif logp >= 1.5:
            ppb_tier = "Moderate (50% - 90%)"
            ppb_note = "Balanced free fraction available for target interaction and tissue distribution."
        else:
            ppb_tier = "Low (< 50%)"
            ppb_note = "High free fraction in plasma; rapid distribution and predominantly renal clearance."

        # Cytochrome P450 (CYP) inhibition alerts
        cyp_alerts = []
        for cyp_name, alert in _CYP_ALERTS.items():
            patt = Chem.MolFromSmarts(alert["smarts"])
            if patt and mol.HasSubstructMatch(patt):
                cyp_alerts.append({
                    "cyp": cyp_name.replace("_inhibition", ""),
                    "description": alert["name"],
                    "status": "Potential Inhibitor"
                })

        # Toxicity & Screening Alerts: PAINS (Pan-Assay Interference Compounds)
        pains_matches = _pains_catalog.GetMatches(mol)
        pains_list = [entry.GetDescription() for entry in pains_matches]

        # Brenk structural alerts (reactive / unstable / toxicophores)
        brenk_matches = _brenk_catalog.GetMatches(mol)
        brenk_list = [entry.GetDescription() for entry in brenk_matches]

        # NIH Screening alert catalog (reactive / unwanted groups)
        nih_matches = _nih_catalog.GetMatches(mol)
        nih_list = [entry.GetDescription() for entry in nih_matches]

        # ZINC Screening alert catalog (problematic / aggregator chemotypes)
        zinc_matches = _zinc_catalog.GetMatches(mol)
        zinc_list = [entry.GetDescription() for entry in zinc_matches]

        total_alerts = len(pains_list) + len(brenk_list) + len(nih_list) + len(zinc_list)

        return {
            "physicochemical": {
                "molecular_weight": {"value": mw, "unit": "g/mol", "description": "Molecular Weight", "source": f"RDKit v{rdkit.__version__}"},
                "logp": {"value": logp, "unit": "unitless", "description": "Wildman-Crippen MolLogP", "source": f"RDKit v{rdkit.__version__}"},
                "hbd": {"value": hbd, "unit": "count", "description": "Hydrogen Bond Donors (Lipinski)", "source": f"RDKit v{rdkit.__version__}"},
                "hba": {"value": hba, "unit": "count", "description": "Hydrogen Bond Acceptors (Lipinski)", "source": f"RDKit v{rdkit.__version__}"},
                "h_bond_donors": {"value": hbd, "unit": "count", "description": "Hydrogen Bond Donors (Lipinski)", "source": f"RDKit v{rdkit.__version__}"},
                "h_bond_acceptors": {"value": hba, "unit": "count", "description": "Hydrogen Bond Acceptors (Lipinski)", "source": f"RDKit v{rdkit.__version__}"},
                "tpsa": {"value": tpsa, "unit": "Å²", "description": "Topological Polar Surface Area", "source": f"RDKit v{rdkit.__version__}"},
                "rotatable_bonds": {"value": rotb, "unit": "count", "description": "Rotatable single bonds", "source": f"RDKit v{rdkit.__version__}"},
                "molar_refractivity": {"value": mr, "unit": "cm³/mol", "description": "Molar Refractivity", "source": f"RDKit v{rdkit.__version__}"},
                "heavy_atoms": {"value": heavy_atoms, "unit": "count", "description": "Non-hydrogen heavy atoms", "source": f"RDKit v{rdkit.__version__}"},
                "aromatic_rings": {"value": aromatic_rings, "unit": "count", "description": "Aromatic ring systems", "source": f"RDKit v{rdkit.__version__}"},
                "fsp3": {"value": fsp3, "unit": "ratio", "description": "Carbon saturation index (sp3 carbons / total carbons)", "source": f"RDKit v{rdkit.__version__}"},
                "formal_charge": {"value": Chem.GetFormalCharge(mol), "unit": "elementary charge", "description": "Net formal charge", "source": f"RDKit v{rdkit.__version__}"},
                "sascore": {
                    "value": sa_score,
                    "unit": "scale 1-10",
                    "interpretation": sa_desc,
                    "description": "Calculated Synthetic Accessibility Score (1=easy, 10=very difficult)",
                    "method": "Ertl & Schuffenhauer fragment-based penalty",
                    "citation": "Ertl & Schuffenhauer, J. Cheminform. 2009",
                    "source": f"RDKit Contrib sascorer"
                },
                "metadata": {
                    "cheminformatics_engine": f"RDKit v{rdkit.__version__}",
                    "logp_algorithm": "Wildman-Crippen (J. Chem. Inf. Comput. Sci. 1999, 39, 868-873)",
                    "tpsa_algorithm": "Prasanna & Doerksen / Ertl et al. (J. Med. Chem. 2000, 43, 3714-3717)",
                    "hbd_hba_definition": "Lipinski standard (O/N atoms with attached H for HBD; O/N atoms for HBA)"
                }
            },
            "drug_likeness": {
                "lipinski": {
                    "rule": "Rule of Five (Pfizer)",
                    "status": lipinski_status,
                    "violations_count": len(lipinski_violations),
                    "violations": lipinski_violations,
                    "mw_value": mw,
                    "logp_value": logp,
                    "hbd_value": hbd,
                    "hba_value": hba,
                    "citation": "Lipinski et al., Adv. Drug Deliv. Rev. 1997"
                },
                "veber": {
                    "rule": "Veber Oral Bioavailability (GSK)",
                    "status": "Pass" if veber_pass else "Fail",
                    "violations_count": len(veber_violations),
                    "violations": veber_violations,
                    "citation": "Veber et al., J. Med. Chem. 2002"
                },
                "ghose": {
                    "rule": "Ghose Drug-likeness Filter",
                    "status": "Pass" if ghose_pass else "Fail",
                    "violations_count": len(ghose_violations),
                    "violations": ghose_violations,
                    "citation": "Ghose et al., J. Comb. Chem. 1999"
                }
            },
            "pharmacokinetics": {
                "gi_absorption": {
                    "level": gi_absorption,
                    "model": "BOILED-Egg published ellipse model (WLogP & TPSA)",
                    "citation": "Daina & Zoete, ChemMedChem 2016, 11, 1117-1121"
                },
                "bbb_permeation": {
                    "status": bbb_status,
                    "model": "BOILED-Egg yolk ellipse model (WLogP & TPSA)",
                    "citation": "Daina & Zoete, ChemMedChem 2016, 11, 1117-1121"
                },
                "plasma_protein_binding": {
                    "tier": ppb_tier,
                    "rationale": ppb_note
                },
                "cyp450_inhibition": {
                    "alerts": cyp_alerts,
                    "count": len(cyp_alerts),
                    "methodology": "heuristic_smarts",
                    "validated": False,
                    "disclaimer": "Exploratory SMARTS substructure heuristics — Not a validated quantitative predictor",
                    "status": f"{len(cyp_alerts)} heuristic alert(s)" if cyp_alerts else "No structural alerts"
                }
            },
            "medicinal_chemistry_safety": {
                "total_alerts_count": total_alerts,
                "overall_status": "Clean (No alerts)" if total_alerts == 0 else f"{total_alerts} alert(s) across 4 catalogs",
                "pains_alerts": {
                    "count": len(pains_list),
                    "alerts": pains_list,
                    "status": "Clear" if len(pains_list) == 0 else f"{len(pains_list)} PAINS alert(s) detected",
                    "citation": "Baell & Holloway, J. Med. Chem. 2010"
                },
                "brenk_alerts": {
                    "count": len(brenk_list),
                    "alerts": brenk_list[:5],
                    "status": "Clear" if len(brenk_list) == 0 else f"{len(brenk_list)} structural alert(s)",
                    "citation": "Brenk et al., ChemMedChem 2008"
                },
                "nih_alerts": {
                    "count": len(nih_list),
                    "alerts": nih_list[:5],
                    "status": "Clear" if len(nih_list) == 0 else f"{len(nih_list)} NIH alert(s)",
                    "citation": "NIH Molecular Libraries Program Clinical Alerts"
                },
                "zinc_alerts": {
                    "count": len(zinc_list),
                    "alerts": zinc_list[:5],
                    "status": "Clear" if len(zinc_list) == 0 else f"{len(zinc_list)} ZINC alert(s)",
                    "citation": "Irwin & Shoichet, J. Chem. Inf. Model. 2005"
                }
            }
        }
