import os
import json
import requests
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
    def _is_valid_narrative(content: str) -> bool:
        """Validate that LLM response is complete, well-formed, and contains all required sections."""
        if not content or not isinstance(content, str):
            return False
        trimmed = content.strip()
        if len(trimmed) < 600:
            return False
        
        # Check that it contains required sections
        has_sections = ("1." in trimmed and "2." in trimmed and "3." in trimmed and "4." in trimmed)
        if not has_sections:
            return False

        # Reject if model returned scratchpad thoughts
        if "The user wants" in trimmed or "Let me structure" in trimmed or trimmed.startswith("Key data points:"):
            return False

        # Check for unclosed or truncated sentences
        words = trimmed.split()
        if words:
            last_word = words[-1].lower().rstrip(".*_#`")
            dangling_words = {"and", "or", "is", "the", "a", "an", "with", "at", "by", "of", "to", "in", "for", "as", "values"}
            if last_word in dangling_words:
                return False
        if trimmed.endswith((",", "-", "(", "/", ":", "—", ";")):
            return False

        return True

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
                    narr = cached.get("narrative", "")
                    # Only return from cache if it is a complete, valid narrative!
                    if NarrativeExplainer._is_valid_narrative(narr) or cached.get("is_fallback", False):
                        cached["cached"] = True
                        return cached
                    else:
                        # Stale or corrupted cache file - delete it
                        try:
                            cache_file.unlink()
                        except Exception:
                            pass
            except Exception:
                pass

        # 1. First, attempt LLM call if API key is provided
        effective_key = api_key or OPENROUTER_API_KEY or GEMINI_API_KEY
        if effective_key:
            try:
                llm_result = NarrativeExplainer._call_llm(report_data, effective_key, provider, model_to_use)
                if llm_result:
                    llm_response, actual_model = llm_result
                    result = {
                        "narrative": llm_response,
                        "source": f"Bindora AI Explainer ({actual_model} via OpenRouter)",
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
    def _call_llm(data: Dict[str, Any], key: str, provider: str, model: str) -> Optional[tuple]:
        """Call OpenRouter API with anti-hallucination prompt, disabled reasoning overhead, and multi-model fallback."""
        system_prompt = (
            "You are a Senior Computational Chemist and Molecular Docking Scientist providing an academic research dossier briefing.\n"
            "CRITICAL SCIENTIFIC & TERMINOLOGY CONSTRAINTS:\n"
            "1. Output ONLY the final Markdown formatted briefing. Do NOT output any internal chain-of-thought, planning notes, or meta-comments.\n"
            "2. NEVER invent, hallucinate, or alter any numbers, scores, or constants. Use the exact values provided in the JSON data.\n"
            "3. TERMINOLOGY HONESTY:\n"
            "   - Refer to AutoDock Vina output strictly as 'AutoDock Vina Docking Score (kcal/mol)' or 'predicted binding score'. Do NOT describe it as an experimentally measured binding free energy (ΔG°).\n"
            "   - Refer to Kd strictly as 'affinity-derived Kd-like estimate (model-derived)', derived via standard isothermal thermodynamic approximation (T = 298.15 K, RT ≈ 0.592 kcal/mol). Explicitly clarify that this value is mathematically derived from the docking score and is not an experimentally measured or rigorously calculated thermodynamic Kd.\n"
            "   - Clearly separate: (A) Docking-derived predictions, (B) Calculated RDKit cheminformatics descriptors, (C) Derived thermodynamic estimates, and (D) Experimental wet-lab assay data from ChEMBL (if available).\n"
            "4. CAUTIOUS SCIENTIFIC TONE:\n"
            "   - Use cautious, publication-grade academic prose ('in silico docking predicts', 'computationally modeled interaction', 'theoretical estimate').\n"
            "   - NEVER claim that a docking score proves nanomolar in vivo efficacy, clinical potency, or therapeutic safety.\n"
            "   - Never claim an interaction is 'experimentally validated' solely based on a computational docking score.\n"
            "5. DOCKING CLASSIFICATION:\n"
            "   - If 'experiment_classification' is present in the data, explicitly state whether the run is 'Native Redocking (Self-Validation)' or 'Cross-Docking / Benchmark Docking'.\n"
            "   - If Cross-Docking, note that scoring differences relative to literature may arise from receptor conformational adaptation (induced fit) or scoring function differences (e.g., Vina 1.2.5 vs AutoDock 4.2).\n"
            "6. Markdown Tables: When summarizing ADME properties or molecular descriptors, use clean GitHub-flavored markdown tables with standard pipes and headers.\n"
            "7. If experimental cross-check status is 'Computational Prediction Only', state clearly that this is an in silico estimation without deposited wet-lab binding assays in ChEMBL.\n"
            "8. Structure your output into four clean sections with markdown headers:\n"
            "### 1. 3D Binding Mechanism & Active Site Interactions\n"
            "### 2. Pharmacokinetics (ADME) & Drug-Likeness Profile\n"
            "### 3. Bioactivity & Experimental Cross-Validation\n"
            "### 4. Physiological Implications & Target Context"
        )

        user_content = (
            f"Here is the verified experimental and computational payload for the drug-target docking run:\n"
            f"```json\n{json.dumps(data, indent=2)}\n```\n"
            f"Provide a concise, publication-grade pharmacological evaluation in 350-500 words. Begin directly with the report."
        )

        # OpenRouter endpoint
        if "sk-or-" in key or provider == "openrouter" or (not provider.startswith("gemini") and not key.startswith("AIza")):
            url = "https://openrouter.ai/api/v1/chat/completions"
            headers = {
                "Authorization": f"Bearer {key}",
                "Content-Type": "application/json",
                "HTTP-Referer": "https://github.com/Swelo-ui/Bindora",
                "X-Title": "Bindora 3D Drug-Receptor Analyzer"
            }

            # Try requested model first, then fallback to high-reliability models if truncated
            candidate_models = [model, "google/gemini-2.0-flash-001", "meta-llama/llama-3.3-70b-instruct"]
            seen_models = set()

            for cand_model in candidate_models:
                if not cand_model or cand_model in seen_models:
                    continue
                seen_models.add(cand_model)

                body = {
                    "model": cand_model,
                    "messages": [
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": user_content}
                    ],
                    "temperature": 0.2,
                    "max_tokens": 3500,
                    # Suppress reasoning token consumption so the completion is never truncated
                    "reasoning": {"max_tokens": 0}
                }

                try:
                    resp = requests.post(url, json=body, headers=headers, timeout=15)
                    if resp.status_code == 200:
                        res_json = resp.json()
                        choices = res_json.get("choices", [])
                        if choices:
                            finish_reason = choices[0].get("finish_reason")
                            # If finished due to length, it was truncated; skip
                            if finish_reason == "length":
                                print(f"[NARRATIVE] Model {cand_model} truncated by token limit. Trying fallback...")
                                continue

                            msg = choices[0].get("message", {})
                            content = msg.get("content") or ""

                            if NarrativeExplainer._is_valid_narrative(content):
                                return content.strip(), cand_model
                            else:
                                print(f"[NARRATIVE] Model {cand_model} output incomplete or missing sections. Trying fallback...")
                    else:
                        print(f"[NARRATIVE] OpenRouter HTTP {resp.status_code} for {cand_model}: {resp.text[:150]}")
                        # If unauthorized or insufficient credits, do not waste time retrying other models with same key
                        if resp.status_code in (401, 402, 403):
                            break
                except Exception as ex:
                    print(f"[NARRATIVE] Request error with model {cand_model}: {ex}")

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
        cyp_data = pk.get("cyp450_inhibition", [])
        cyp_alerts = cyp_data.get("alerts", []) if isinstance(cyp_data, dict) else cyp_data

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
            cyp_str = "; ".join([f"{c['cyp']} ({c['description']})" for c in cyp_alerts]) + " *(Exploratory SMARTS Heuristic — Not a Validated Predictor)*"
        else:
            cyp_str = "No major CYP450 structural inhibition alerts identified *(Exploratory SMARTS Heuristic)*."

        exp_class = data.get("experiment_classification", {})
        docking_mode = exp_class.get("docking_mode") or data.get("docking_mode", "")
        exp_desc = exp_class.get("description", "")
        mode_line = f"- **Docking Classification:** **{docking_mode}**{f' — {exp_desc}' if exp_desc else ''}\n" if docking_mode else ""

        # Section 1: Binding mechanism
        sec1 = (
            f"### 1. 3D Binding Mechanism & Active Site Interactions\n\n"
            f"AutoDock Vina molecular docking of **{drug}** against receptor **{target}** (PDB ID: `{pdb_id}`) yielded a docking score of **{affinity} kcal/mol**. "
            f"Based on standard isothermal thermodynamic approximation ($T = 298.15\\text{{ K}}$, $RT \\approx 0.592\\text{{ kcal/mol}}$), this corresponds to an affinity-derived Kd-like estimate of approximately **{kd_nm} nM** ({kd_um} µM) *(Model-derived estimate: $K_d = \\exp(\\text{{score}}/RT)$; not an experimental thermodynamic $K_d$)*, placing the predicted score in the **{potency}** tier.\n\n"
            f"{mode_line}"
            f"- **Ligand Efficiency (LE):** Calculated at **{le} kcal/mol/heavy atom** (benchmark target ≥ 0.30 kcal/mol/heavy atom; $|\\text{{score}}| / \\text{{heavy atoms}}$).\n"
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
                f"### 3. Bioactivity & Experimental Cross-Validation\n\n"
                f"**Validation Status:** `[{badge}]` (Green Badge — Literature Assay Records Available)\n\n"
                f"Curated experimental bioactivity records from the **ChEMBL Database** ({crosscheck.get('target_organism', 'Homo sapiens')}) for this drug-target pair include:\n\n"
                f"Found **{len(records)}** experimental assay record(s):\n"
                f"{rec_str}\n\n"
                f"**Methodological Note:** AutoDock Vina scores reflect empirical scoring function approximations, whereas ChEMBL records document physical in vitro biological assays. Literature benchmarks provide comparative biological context rather than direct free energy equivalences."
            )
        else:
            sec3 = (
                f"### 3. Bioactivity & Experimental Cross-Validation\n\n"
                f"**Validation Status:** `[{badge}]` (Amber Badge — In Silico Prediction Only)\n\n"
                f"No direct curated experimental assay record was identified in public databases (ChEMBL) matching this specific drug and target combination. "
                f"**Important Research Note:** All presented binding energies and downstream affinity metrics are strictly computational estimates. They serve as valuable hypothesis-generating models for virtual screening, but must not be cited as confirmed wet-lab experimental constants."
            )

        # Section 4: Virtual Physiology
        sec4 = (
            f"### 4. Downstream Physiological & Target Pathway Interpretation\n\n"
            f"{f'Target Biological Role: {target_fn[:300]}...' if target_fn else f'Target Receptor: {target}'}\n\n"
            f"The predicted docking pose and energetic score of {affinity} kcal/mol reflect structural shape complementarity and key residue interactions within the pocket. "
            f"In an educational and drug discovery context, these findings provide structural hypotheses for lead optimization and understanding receptor-ligand pharmacophoric determinants."
        )

        return f"{sec1}\n\n---\n\n{sec2}\n\n---\n\n{sec3}\n\n---\n\n{sec4}"
