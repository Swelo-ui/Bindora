import os
import json
import urllib.request
import urllib.parse
from pathlib import Path
from typing import List, Dict, Any, Optional

from rdkit import Chem, RDConfig
from rdkit.Chem import rdMolChemicalFeatures

from backend.config import CACHE_DIR, CHEMBL_BASE_URL

HEADERS = {"User-Agent": "Bindora-Research-Tool/1.0 (academic; +https://github.com/Swelo-ui/Bindora)"}

# Initialize RDKit Chemical Feature Factory
_FDEF_PATH = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
_FEATURE_FACTORY = None
if os.path.exists(_FDEF_PATH):
    try:
        _FEATURE_FACTORY = rdMolChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)
    except Exception as e:
        print(f"[PHARMACOPHORE INIT WARNING] Feature factory build failed: {e}")

class PharmacophoreService:
    """Service to derive consensus ligand pharmacophore profiles from known actives and screen candidates."""

    @staticmethod
    def fetch_target_actives(
        target_name: str = "",
        chembl_target_id: Optional[str] = None,
        pdb_id: Optional[str] = None,
        uniprot_accession: Optional[str] = None,
        max_actives: int = 10
    ) -> List[Dict[str, Any]]:
        """Fetch curated high-affinity active molecules (IC50 <= 1000 nM) from ChEMBL for this target."""
        target_name_clean = (target_name or "").strip()
        pdb_clean = (pdb_id or "").strip().upper()
        acc_clean = (uniprot_accession or "").strip().upper()
        chembl_tid = (chembl_target_id or "").strip().upper()

        if not target_name_clean and not chembl_tid and not pdb_clean and not acc_clean:
            return []

        ident_key = chembl_tid or acc_clean or pdb_clean or urllib.parse.quote_plus(target_name_clean[:40].lower())
        cache_key = f"chembl_actives_{ident_key}.json"
        cache_file = CACHE_DIR / cache_key
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    if cached and len(cached) >= 3:
                        return cached
            except Exception:
                pass

        tid = chembl_tid

        # 1. Try UniProt accession direct match (most accurate for biological targets)
        if not tid and acc_clean:
            try:
                u_url = f"{CHEMBL_BASE_URL}/target?target_components__accession={acc_clean}&format=json"
                req = urllib.request.Request(u_url, headers=HEADERS)
                with urllib.request.urlopen(req, timeout=10) as resp:
                    udata = json.loads(resp.read().decode("utf-8"))
                    targets = udata.get("targets", [])
                    if targets:
                        tid = targets[0].get("target_chembl_id")
            except Exception as e:
                print(f"[PHARMACOPHORE UNIPROT LOOKUP WARNING] {e}")

        # 2. Try PDB ID search in ChEMBL (ChEMBL indexes RCSB PDB complexes)
        if not tid and pdb_clean:
            try:
                p_url = f"{CHEMBL_BASE_URL}/target/search?q={pdb_clean}&format=json"
                req = urllib.request.Request(p_url, headers=HEADERS)
                with urllib.request.urlopen(req, timeout=10) as resp:
                    pdata = json.loads(resp.read().decode("utf-8"))
                    targets = pdata.get("targets", [])
                    if targets:
                        tid = targets[0].get("target_chembl_id")
            except Exception as e:
                print(f"[PHARMACOPHORE PDB LOOKUP WARNING] {e}")

        # 3. Try Target Name search in ChEMBL
        if not tid and target_name_clean:
            search_queries = [target_name_clean]
            if len(target_name_clean) > 30:
                for kw in ["HIV-1 PROTEASE", "HIV PROTEASE", "PROTEASE", "KINASE", "REVERSE TRANSCRIPTASE", "POLYPROTEIN"]:
                    if kw in target_name_clean.upper():
                        search_queries.append(kw.title())
                        break
                words = target_name_clean.split()
                if len(words) > 3:
                    search_queries.append(" ".join(words[:3]))

            for q in search_queries:
                try:
                    turl = f"{CHEMBL_BASE_URL}/target/search?q={urllib.parse.quote(q)}&format=json"
                    req = urllib.request.Request(turl, headers=HEADERS)
                    with urllib.request.urlopen(req, timeout=10) as resp:
                        tdata = json.loads(resp.read().decode("utf-8"))
                        targets = tdata.get("targets", [])
                        if targets:
                            tid = targets[0].get("target_chembl_id")
                            break
                except Exception as e:
                    print(f"[PHARMACOPHORE TARGET LOOKUP WARNING] {e}")

        if not tid:
            return []

        # Query activities for target: first try potent IC50 <= 1000 nM, if < 3 try up to 10000 nM
        actives = []
        seen_smiles = set()

        for lte_val in [1000, 10000]:
            try:
                act_url = (
                    f"{CHEMBL_BASE_URL}/activity?target_chembl_id={tid}&standard_type=IC50&"
                    f"standard_value__lte={lte_val}&standard_units=nM&limit=30&format=json"
                )
                req = urllib.request.Request(act_url, headers=HEADERS)
                with urllib.request.urlopen(req, timeout=12) as resp:
                    act_data = json.loads(resp.read().decode("utf-8"))
                    for act in act_data.get("activities", []):
                        smi = act.get("canonical_smiles")
                        if smi and smi not in seen_smiles:
                            m = Chem.MolFromSmiles(smi)
                            if m and m.GetNumHeavyAtoms() >= 6:
                                seen_smiles.add(smi)
                                actives.append({
                                    "chembl_id": act.get("molecule_chembl_id"),
                                    "smiles": smi,
                                    "ic50_nm": float(act.get("standard_value", 0.0)),
                                    "assay_type": act.get("standard_type", "IC50")
                                })
                        if len(actives) >= max_actives:
                            break
                if len(actives) >= 3:
                    break
            except Exception as e:
                print(f"[PHARMACOPHORE ACTIVES FETCH WARNING] {e}")

        if actives:
            try:
                with open(cache_file, "w", encoding="utf-8") as f:
                    json.dump(actives, f, indent=2)
            except Exception:
                pass

        return actives

    @staticmethod
    def extract_features(smiles_or_mol: Any) -> Dict[str, Any]:
        """Extract chemical features (Donors, Acceptors, Aromatic rings, Hydrophobes) using RDKit BaseFeatures."""
        global _FEATURE_FACTORY
        if _FEATURE_FACTORY is None:
            if os.path.exists(_FDEF_PATH):
                _FEATURE_FACTORY = rdMolChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)
            else:
                return {}

        if isinstance(smiles_or_mol, str):
            mol = Chem.MolFromSmiles(smiles_or_mol)
        else:
            mol = smiles_or_mol

        if not mol:
            return {}

        counts = {
            "Donor": 0,
            "Acceptor": 0,
            "Aromatic": 0,
            "Hydrophobe": 0,
            "PosIonizable": 0,
            "NegIonizable": 0
        }

        num_feats = _FEATURE_FACTORY.GetNumMolFeatures(mol)
        for i in range(num_feats):
            feat = _FEATURE_FACTORY.GetMolFeature(mol, i)
            fam = feat.GetFamily()
            if fam in counts:
                counts[fam] += 1
            elif fam == "LumpedHydrophobe":
                counts["Hydrophobe"] += 1

        return counts

    @staticmethod
    def build_consensus_profile(actives: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Build consensus pharmacophore requirements from a collection of known target actives."""
        if not actives:
            return {"error": "No active compounds provided"}

        all_counts = []
        for act in actives:
            smi = act.get("smiles") or act.get("canonical_smiles", "")
            feats = PharmacophoreService.extract_features(smi)
            if feats:
                all_counts.append(feats)

        if not all_counts:
            return {"error": "Failed to extract features from actives"}

        consensus_features = {}
        for fam in ["Donor", "Acceptor", "Aromatic", "Hydrophobe", "PosIonizable", "NegIonizable"]:
            vals = [c[fam] for c in all_counts]
            non_zeros = sum(1 for v in vals if v > 0)
            presence_ratio = non_zeros / len(all_counts)

            # Feature is part of consensus if present in >= 50% of actives
            if presence_ratio >= 0.5:
                # Median count
                vals_sorted = sorted(vals)
                median_val = vals_sorted[len(vals_sorted) // 2]
                consensus_features[fam] = {
                    "required": max(1, median_val),
                    "presence_frequency": round(presence_ratio, 2),
                    "min": min(vals),
                    "max": max(vals)
                }

        return {
            "actives_count": len(all_counts),
            "consensus_features": consensus_features,
            "core_requirements": {fam: d["required"] for fam, d in consensus_features.items()}
        }

    @staticmethod
    def match_candidate(candidate_smiles: str, consensus_profile: Dict[str, Any]) -> Dict[str, Any]:
        """Match a candidate compound against the consensus pharmacophore profile."""
        if not candidate_smiles or "consensus_features" not in consensus_profile:
            return {"score": 0.0, "status": "Invalid profile"}

        cand_feats = PharmacophoreService.extract_features(candidate_smiles)
        if not cand_feats:
            return {"score": 0.0, "status": "Invalid molecule"}

        reqs = consensus_profile["consensus_features"]
        if not reqs:
            return {"score": 100.0, "status": "No consensus requirements"}

        total_weight = len(reqs)
        match_points = 0.0
        details = {}

        for fam, req_info in reqs.items():
            req_count = req_info["required"]
            actual_count = cand_feats.get(fam, 0)
            if actual_count >= req_count:
                match_points += 1.0
                details[fam] = f"Satisfied ({actual_count} >= {req_count})"
            elif actual_count > 0:
                match_points += (actual_count / req_count) * 0.7
                details[fam] = f"Partial ({actual_count}/{req_count})"
            else:
                details[fam] = f"Missing (0/{req_count})"

        score_pct = round((match_points / total_weight) * 100.0, 1)
        if score_pct >= 80.0:
            status = "Strong Pharmacophore Match"
        elif score_pct >= 50.0:
            status = "Moderate Pharmacophore Match"
        else:
            status = "Weak / Low Match"

        return {
            "match_score_pct": score_pct,
            "status": status,
            "candidate_features": cand_feats,
            "feature_alignment": details
        }
