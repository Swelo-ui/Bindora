import math
from typing import Dict, Any, List, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, FilterCatalog

# Initialize PAINS and Brenk structural alert catalogs once for efficiency
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

_brenk_params = FilterCatalog.FilterCatalogParams()
_brenk_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
_brenk_catalog = FilterCatalog.FilterCatalog(_brenk_params)

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

        # 4. Pharmacokinetic estimations (GI Absorption & BBB Permeation via Egan BOILED-Egg model)
        # BOILED-Egg ellipse criteria: TPSA <= 131.6 and -0.4 <= LogP <= 5.6 indicates high GI absorption
        gi_high = (tpsa <= 131.6) and (-0.4 <= logp <= 5.6)
        gi_absorption = "High" if gi_high else "Moderate / Low"

        # BBB permeability: TPSA < 90 Å², MW < 400 Da, 1.2 <= LogP <= 3.2
        bbb_permeant = (tpsa < 90.0) and (mw < 400.0) and (1.0 <= logp <= 3.5)
        bbb_status = "Permeant (Likely crosses BBB)" if bbb_permeant else "Non-permeant (Low CNS penetration likelihood)"

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

        return {
            "physicochemical": {
                "molecular_weight": {"value": mw, "unit": "g/mol", "description": "Molecular Weight"},
                "logp": {"value": logp, "unit": "unitless", "description": "Wildman-Crippen calculated lipophilicity (MolLogP)"},
                "hbd": {"value": hbd, "unit": "count", "description": "Hydrogen Bond Donors"},
                "hba": {"value": hba, "unit": "count", "description": "Hydrogen Bond Acceptors"},
                "tpsa": {"value": tpsa, "unit": "Å²", "description": "Topological Polar Surface Area"},
                "rotatable_bonds": {"value": rotb, "unit": "count", "description": "Rotatable single bonds"},
                "molar_refractivity": {"value": mr, "unit": "cm³/mol", "description": "Molar Refractivity"},
                "heavy_atoms": {"value": heavy_atoms, "unit": "count", "description": "Non-hydrogen heavy atoms"},
                "aromatic_rings": {"value": aromatic_rings, "unit": "count", "description": "Aromatic ring systems"},
                "fsp3": {"value": fsp3, "unit": "ratio", "description": "Carbon saturation index (sp3 carbons / total carbons)"}
            },
            "drug_likeness": {
                "lipinski": {
                    "rule": "Rule of Five (Pfizer)",
                    "status": lipinski_status,
                    "violations_count": len(lipinski_violations),
                    "violations": lipinski_violations,
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
                    "model": "BOILED-Egg / Egan Ellipse criteria (WLogP & TPSA)",
                    "citation": "Daina & Zoete, ChemMedChem 2016; Egan et al., J. Med. Chem. 2000"
                },
                "bbb_permeation": {
                    "status": bbb_status,
                    "model": "Clark's polar surface area (<90 Å²) & molecular size heuristic",
                    "citation": "Clark, Drug Discov Today 2003"
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
                "pains_alerts": {
                    "count": len(pains_list),
                    "alerts": pains_list,
                    "status": "Clear" if len(pains_list) == 0 else f"{len(pains_list)} PAINS alert(s) detected",
                    "citation": "Baell & Holloway, J. Med. Chem. 2010"
                },
                "brenk_alerts": {
                    "count": len(brenk_list),
                    "alerts": brenk_list[:5], # top 5
                    "status": "Clear" if len(brenk_list) == 0 else f"{len(brenk_list)} structural alert(s)",
                    "citation": "Brenk et al., ChemMedChem 2008"
                }
            }
        }
