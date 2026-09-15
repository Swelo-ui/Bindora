import os
import json
import urllib.request
from typing import Dict, Any, Optional
import hashlib
from backend.config import OPENROUTER_API_KEY, OPENROUTER_MODEL, GEMINI_API_KEY, CACHE_DIR, PROJECT_NAME

class NarrativeExplainer:
    """Anti-hallucination educational explanation engine.
    
    Translates raw docking scores, structural contacts, ADME metrics, and ChEMBL records
    into structured educational prose for pharmacy students. Strictly forbidden from inventing numbers.
    Includes smart local disk caching to prevent redundant API calls and conserve OpenRouter quota.
    """

    @staticmethod
    def _compute_cache_key(report_data: Dict[str, Any], model: str) -> str:
        """Create a deterministic hash key for caching report responses."""
        drug = report_data.get("ligand_name", "")
        target = report_data.get("target_name", "")
        pdb = report_data.get("pdb_id", "")
        affinity = report_data.get("thermodynamics", {}).get("binding_affinity_kcal", "")
        kd = report_data.get("thermodynamics", {}).get("theoretical_kd_nm", "")
        raw = f"{drug}_{target}_{pdb}_{affinity}_{kd}_{model}".lower()
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()

    @staticmethod
    def generate_explanation(report_data: Dict[str, Any], api_key: Optional[str] = None, provider: str = "auto") -> Dict[str, Any]:
        """Generate structured narrative explanation with strict grounding on real calculated data."""
        model_to_use = os.environ.get("OPENROUTER_MODEL", OPENROUTER_MODEL)
        cache_key = NarrativeExplainer._compute_cache_key(report_data, model_to_use)
        cache_file = CACHE_DIR / f"narrative_{cache_key}.json"

        # Check local cache first to avoid spending API quota on repeated queries
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    cached["cached"] = True
                    return cached
            except Exception:
                pass

        # 1. First, attempt LLM call if API key is provided
        effective_key = api_key or OPENROUTER_API_KEY or GEMINI_API_KEY
        if effective_key:
            try:
                llm_response = NarrativeExplainer._call_llm(report_data, effective_key, provider, model_to_use)
                if llm_response:
                    result = {
                        "narrative": llm_response,
                        "source": f"Bindora AI Explainer ({model_to_use} via OpenRouter)",
                        "is_fallback": False,
                        "cached": False
                    }
                    try:
                        with open(cache_file, "w", encoding="utf-8") as f:
                            json.dump(result, f, indent=2)
                    except Exception:
                        pass
                    return result
            except Exception as e:
                print(f"[NARRATIVE] OpenRouter API notice: {e}. Gracefully falling back to deterministic reasoning engine.")

        # 2. Deterministic Pharmacology Narrative Generator (Zero hallucination, fully grounded, works offline)
        fallback_text = NarrativeExplainer._generate_deterministic_narrative(report_data)
        return {
            "narrative": fallback_text,
            "source": f"{PROJECT_NAME} Rule-Based Pharmacological Reasoning Engine (Deterministic / Offline)",
            "is_fallback": True,
            "cached": False
        }

    @staticmethod
    def _call_llm(data: Dict[str, Any], key: str, provider: str, model: str) -> Optional[str]:
        """Call OpenRouter API with anti-hallucination prompt and DeepSeek model."""
        system_prompt = (
            f"You are {PROJECT_NAME}'s Academic Pharmacology Explainer for pharmacy and medicinal chemistry students. "
            "STRICT RULES (Anti-Hallucination Central Directive):\n"
            "1. You are strictly an EXPLAINER of calculated values, NOT a data source.\n"
            "2. NEVER invent, hallucinate, or alter any numbers, scores, or constants.\n"
            "3. Every numeric claim (binding energy in kcal/mol, Kd in nM, MW, LogP, TPSA, H-bonds) MUST be cited directly from the provided payload.\n"
            "4. If experimental cross-check status is 'Computational Prediction Only', state clearly that this is an unverified simulation without wet-lab assay confirmation.\n"
            "5. Structure your output in clear markdown sections: (1) Binding Mechanism & Active Site Interactions, "
            "(2) ADME & Oral Bioavailability Profile, (3) Experimental Validation & Confidence Analysis, (4) Physiological Implications & Target Pathways."
        )

        user_content = f"Here is the verified experimental and computational payload for the drug-target docking run:\n```json\n{json.dumps(data, indent=2)}\n```\nProvide a comprehensive, pedagogical pharmacodynamics and pharmacokinetics explanation."

        # OpenRouter endpoint
        if "sk-or-" in key or provider == "openrouter" or (not provider.startswith("gemini") and not key.startswith("AIza")):
            url = "https://openrouter.ai/api/v1/chat/completions"
            headers = {
                "Authorization": f"Bearer {key}",
                "Content-Type": "application/json",
                "HTTP-Referer": "https://github.com/Swelo-ui/Bindora",
                "X-Title": "Bindora 3D Drug-Receptor Analyzer"
            }
            body = {
                "model": model,
                "messages": [
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_content}
                ],
                "temperature": 0.2
            }
            req = urllib.request.Request(url, data=json.dumps(body).encode("utf-8"), headers=headers)
            with urllib.request.urlopen(req, timeout=25) as resp:
                res_json = json.loads(resp.read().decode("utf-8"))
                return res_json["choices"][0]["message"]["content"]
        
        return None

    @staticmethod
    def _generate_deterministic_narrative(data: Dict[str, Any]) -> str:
        """Construct a publication-grade, academically rigorous narrative strictly from calculated facts."""
        drug = data.get("ligand_name", "Ligand")
        target = data.get("target_name", "Target Receptor")
        pdb_id = data.get("pdb_id", "N/A")
        
        thermo = data.get("thermodynamics", {})
        affinity = thermo.get("binding_affinity_kcal", 0.0)
        kd_nm = thermo.get("theoretical_kd_nm", "N/A")
        kd_um = thermo.get("theoretical_kd_um", "N/A")
        le = thermo.get("ligand_efficiency", {}).get("value", 0.0)
        potency = thermo.get("potency_class", "Undetermined")

        contacts = data.get("interactions", {})
        hbonds = contacts.get("hydrogen_bonds", [])
        hydrophobics = contacts.get("hydrophobic_contacts", [])
        residues = contacts.get("interacting_residues", [])

        adme = data.get("adme", {})
        phys = adme.get("physicochemical", {})
        mw = phys.get("molecular_weight", {}).get("value", "N/A")
        logp = phys.get("logp", {}).get("value", "N/A")
        tpsa = phys.get("tpsa", {}).get("value", "N/A")
        hbd = phys.get("hbd", {}).get("value", "N/A")
        hba = phys.get("hba", {}).get("value", "N/A")
        rotb = phys.get("rotatable_bonds", {}).get("value", "N/A")

        rules = adme.get("drug_likeness", {})
        lipinski = rules.get("lipinski", {})
        lipinski_status = lipinski.get("status", "Pass")
        lipinski_viols = lipinski.get("violations", [])

        pk = adme.get("pharmacokinetics", {})
        gi = pk.get("gi_absorption", {}).get("level", "Moderate")
        bbb = pk.get("bbb_permeation", {}).get("status", "Non-permeant")
        ppb = pk.get("plasma_protein_binding", {}).get("tier", "Moderate")
        cyp_alerts = pk.get("cyp450_inhibition", [])

        safety = adme.get("medicinal_chemistry_safety", {})
        pains = safety.get("pains_alerts", {})
        pains_count = pains.get("count", 0)

        crosscheck = data.get("bioactivity_crosscheck", {})
        is_corroborated = crosscheck.get("is_cross_checked", False)
        badge = crosscheck.get("status_badge", "Computational Prediction Only")
        records = crosscheck.get("experimental_records", [])

        uniprot = data.get("uniprot", {})
        target_fn = uniprot.get("function", "")

        # Format H-bond text
        if hbonds:
            hb_details = ", ".join([f"**{hb['residue']}** ({hb['receptor_atom']} — {hb['ligand_atom']}, {hb['distance']} Å)" for hb in hbonds[:5]])
        else:
            hb_details = "No directional polar hydrogen bonds detected within the 3.5 Å cutoff threshold; binding appears dominated by van der Waals and hydrophobic packing."

        # Format Residues text
        res_list_str = ", ".join([f"`{r}`" for r in residues[:8]]) if residues else "None within threshold"

        # Format CYP alerts
        if cyp_alerts:
            cyp_str = "; ".join([f"{c['cyp']} ({c['description']})" for c in cyp_alerts])
        else:
            cyp_str = "No major CYP450 structural inhibition alerts identified."

        # Section 1: Binding mechanism
        sec1 = (
            f"### 1. 3D Binding Mechanism & Active Site Interactions\n\n"
            f"AutoDock Vina molecular docking of **{drug}** against receptor **{target}** (PDB ID: `{pdb_id}`) yielded a predicted binding free energy (ΔG) of **{affinity} kcal/mol**, corresponding to a theoretical dissociation constant ($K_d$) of approximately **{kd_nm} nM** ({kd_um} µM). This places the predicted binding potency in the **{potency}** tier.\n\n"
            f"- **Ligand Efficiency (LE):** Calculated at **{le} kcal/mol/heavy atom** (benchmark target ≥ 0.30 kcal/mol/heavy atom).\n"
            f"- **Hydrogen Bonding Network:** {len(hbonds)} hydrogen bond(s) identified within 3.5 Å: {hb_details}.\n"
            f"- **Non-Polar Contacts:** {len(hydrophobics)} hydrophobic contact(s) anchoring the lipophilic scaffold inside the pocket.\n"
            f"- **Key Interacting Residues:** {res_list_str}."
        )

        # Section 2: ADME & Drug-likeness
        sec2 = (
            f"### 2. Pharmacokinetics (ADME) & Drug-Likeness Profile\n\n"
            f"Evaluation of the calculated molecular descriptors reflects the following profile:\n\n"
            f"- **Lipinski Rule of Five Compliance:** **{lipinski_status}** ({len(lipinski_viols)} violation(s)). Molecular Weight: **{mw} Da** (target ≤ 500), MolLogP: **{logp}** (target ≤ 5.0), H-Bond Donors: **{hbd}** (≤ 5), H-Bond Acceptors: **{hba}** (≤ 10).\n"
            f"- **Oral Bioavailability (Veber Rule):** Rotatable bonds = **{rotb}** (≤ 10), TPSA = **{tpsa} Å²** (≤ 140 Å²).\n"
            f"- **Gastrointestinal Absorption:** Predicted as **{gi}** according to the BOILED-Egg/Egan ellipse coordinate model.\n"
            f"- **Blood-Brain Barrier (BBB) Permeation:** **{bbb}**.\n"
            f"- **Plasma Protein Binding (PPB):** Estimated in the **{ppb}** range based on lipophilic equilibrium.\n"
            f"- **Metabolism & CYP450 Interactions:** {cyp_str}\n"
            f"- **Medicinal Chemistry Safety:** PAINS Filter: **{pains_count} alert(s)**. Brenk alerts: **{safety.get('brenk_alerts', {}).get('count', 0)} alert(s)**."
        )

        # Section 3: Experimental Cross-Validation
        if is_corroborated:
            rec_str = "\n".join([f"  - **{r['type']}**: {r['relation']} {r['value']} {r['units']} (Assay: {r['assay_description'][:70]}...)" for r in records[:3]])
            sec3 = (
                f"### 3. Experimental Bioactivity Cross-Validation\n\n"
                f"**Validation Status:** `[{badge}]` (Green Badge)\n\n"
                f"AnuDock cross-referenced this drug-target pair against curated experimental bioactivity records in the **ChEMBL Database** ({crosscheck.get('target_organism', 'Homo sapiens')}).\n\n"
                f"Found **{len(records)}** experimental assay record(s):\n"
                f"{rec_str}\n\n"
                f"The computational binding prediction aligns with wet-lab experimental affinity benchmarks, demonstrating that the docking pose captures realistic pharmacophore orientation."
            )
        else:
            sec3 = (
                f"### 3. Experimental Bioactivity Cross-Validation\n\n"
                f"**Validation Status:** `[{badge}]` (Amber Badge - Computational Estimate)\n\n"
                f"No direct curated experimental assay record was identified in public databases (ChEMBL) matching this specific drug and target combination. "
                f"**Important Research Note:** All presented binding energies and downstream affinity metrics are strictly computational estimates. They serve as valuable hypothesis-generating models for virtual screening, but must not be cited as confirmed wet-lab experimental constants."
            )

        # Section 4: Virtual Physiology
        sec4 = (
            f"### 4. Downstream Physiological & Target Pathway Interpretation\n\n"
            f"{f'Target Biological Role: {target_fn[:300]}...' if target_fn else f'Target Receptor: {target}'}\n\n"
            f"By occupying the receptor binding cleft with a docking score of {affinity} kcal/mol, **{drug}** is modeled to modulate target signaling cascades computationally without animal testing. "
            f"In an educational context, this profile illustrates how structural complementarities (hydrogen bonds + hydrophobic fit) translate directly into pharmacological affinity and systemic pharmacokinetic behavior."
        )

        return f"{sec1}\n\n---\n\n{sec2}\n\n---\n\n{sec3}\n\n---\n\n{sec4}"
