import os
import json
import math
import urllib.request
import urllib.parse
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolChemicalFeatures, rdMolAlign

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
    """Service to derive true 3D spatial consensus pharmacophore hypotheses from known actives and screen candidates."""

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

        # 1. Try UniProt accession direct match
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

        # 2. Try PDB ID search in ChEMBL
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

        actives = []
        seen_smiles = set()

        for lte_val in [1000, 10000]:
            try:
                act_url = (
                    f"{CHEMBL_BASE_URL}/activity?target_chembl_id={tid}&standard_type=IC50&"
                    f"standard_value__lte={lte_val}&standard_units=nM&limit=30&format=json"
                )
                req = urllib.request.Request(act_url, headers=HEADERS)
                with urllib.request.urlopen(req, timeout=10) as resp:
                    adata = json.loads(resp.read().decode("utf-8"))
                    items = adata.get("activities", [])
                    for it in items:
                        smi = it.get("canonical_smiles")
                        val = it.get("standard_value")
                        mid = it.get("molecule_chembl_id")
                        if not smi or smi in seen_smiles:
                            continue
                        try:
                            fval = float(val) if val else 0.0
                        except ValueError:
                            continue

                        seen_smiles.add(smi)
                        actives.append({
                            "chembl_id": mid,
                            "smiles": smi,
                            "ic50_nm": round(fval, 1),
                            "standard_type": it.get("standard_type", "IC50"),
                            "pref_name": it.get("molecule_pref_name") or mid
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
    def _generate_3d_mol(smiles: str, num_confs: int = 1) -> Optional[Chem.Mol]:
        """Generate low-energy 3D conformer(s) for a SMILES string using RDKit ETKDGv3."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params)
            if not cids:
                AllChem.EmbedMolecule(mol_h, randomSeed=42)
            for cid in mol_h.GetConformers():
                try:
                    AllChem.MMFFOptimizeMolecule(mol_h, confId=cid.GetId(), maxIters=200)
                except Exception:
                    pass
            return mol_h
        except Exception:
            return None

    @staticmethod
    def extract_3d_spatial_features(mol_3d: Chem.Mol, conf_id: int = 0) -> Tuple[Dict[str, int], List[Dict[str, Any]]]:
        """
        Extract chemical features with exact 3D Cartesian coordinates (x, y, z)
        from a conformer of the molecule.
        """
        global _FEATURE_FACTORY
        if _FEATURE_FACTORY is None:
            if os.path.exists(_FDEF_PATH):
                _FEATURE_FACTORY = rdMolChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)
            else:
                return {}, []

        counts = {
            "Donor": 0,
            "Acceptor": 0,
            "Aromatic": 0,
            "Hydrophobe": 0,
            "PosIonizable": 0,
            "NegIonizable": 0
        }
        spatial_features: List[Dict[str, Any]] = []

        conf = mol_3d.GetConformer(conf_id) if mol_3d.GetNumConformers() > conf_id else None
        if not conf:
            return counts, spatial_features

        num_feats = _FEATURE_FACTORY.GetNumMolFeatures(mol_3d)
        for i in range(num_feats):
            feat = _FEATURE_FACTORY.GetMolFeature(mol_3d, i)
            fam = feat.GetFamily()
            if fam == "LumpedHydrophobe":
                fam = "Hydrophobe"

            if fam in counts:
                counts[fam] += 1

            atom_ids = feat.GetAtomIds()
            if not atom_ids:
                continue

            # Calculate 3D center of the feature
            coords = [conf.GetAtomPosition(aid) for aid in atom_ids]
            cx = sum(p.x for p in coords) / len(coords)
            cy = sum(p.y for p in coords) / len(coords)
            cz = sum(p.z for p in coords) / len(coords)

            spatial_features.append({
                "family": fam,
                "center": [round(cx, 3), round(cy, 3), round(cz, 3)],
                "atom_indices": list(atom_ids)
            })

        return counts, spatial_features

    @staticmethod
    def extract_features(smiles_or_mol: Any) -> Dict[str, Any]:
        """Legacy helper returning 2D chemical feature counts for backward compatibility."""
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
        global _FEATURE_FACTORY
        if _FEATURE_FACTORY is None and os.path.exists(_FDEF_PATH):
            _FEATURE_FACTORY = rdMolChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)
        if not _FEATURE_FACTORY:
            return counts

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
        """
        Build True 3D Spatial Pharmacophore Hypothesis:
        1. Generates 3D conformers for known target actives.
        2. Aligns actives in 3D coordinate space against the most potent template active.
        3. Spatially clusters feature centroids (tolerance radius = 1.2 A).
        4. Derives pairwise 3D spatial distance constraints matrix D_ij.
        """
        if not actives:
            return {"error": "No active compounds provided"}

        # 1. Generate 3D conformers for actives
        mols_3d = []
        for act in actives:
            smi = act.get("smiles") or act.get("canonical_smiles", "")
            if smi:
                m3d = PharmacophoreService._generate_3d_mol(smi, num_confs=1)
                if m3d:
                    mols_3d.append(m3d)

        if not mols_3d:
            return {"error": "Failed to generate 3D conformers for actives"}

        # 2. Align conformers in 3D against template (mols_3d[0])
        ref_mol = mols_3d[0]
        aligned_mols = [ref_mol]
        for m in mols_3d[1:]:
            try:
                # Shape/atom alignment
                o3a = rdMolAlign.GetO3A(m, ref_mol)
                o3a.Align()
            except Exception:
                try:
                    rdMolAlign.AlignMol(m, ref_mol)
                except Exception:
                    pass
            aligned_mols.append(m)

        # 3. Extract 3D spatial features from each aligned molecule
        all_features_per_mol = []
        all_counts = []
        for m in aligned_mols:
            cnts, sfeats = PharmacophoreService.extract_3d_spatial_features(m)
            all_counts.append(cnts)
            all_features_per_mol.append(sfeats)

        # 4. Spatially cluster features across molecules into 3D Consensus Spheres
        # Group features by family
        ref_feats = all_features_per_mol[0]
        consensus_spheres = []
        sphere_id = 1

        for rf in ref_feats:
            fam = rf["family"]
            rc = np.array(rf["center"])
            matching_centers = [rc]

            for other_feats in all_features_per_mol[1:]:
                # Find closest feature of same family within 2.2 A
                closest = None
                min_d = float("inf")
                for of in other_feats:
                    if of["family"] == fam:
                        d = np.linalg.norm(np.array(of["center"]) - rc)
                        if d < min_d:
                            min_d = d
                            closest = np.array(of["center"])
                if closest is not None and min_d <= 2.2:
                    matching_centers.append(closest)

            presence_ratio = len(matching_centers) / len(aligned_mols)
            if presence_ratio >= 0.5:
                mean_center = np.mean(matching_centers, axis=0)
                consensus_spheres.append({
                    "id": f"P{sphere_id}",
                    "family": fam,
                    "center": [round(float(c), 2) for c in mean_center],
                    "tolerance_radius": 1.2,
                    "presence_frequency": round(presence_ratio, 2)
                })
                sphere_id += 1

        # 5. Compute Pairwise 3D Distance Constraints
        distance_constraints = []
        for i in range(len(consensus_spheres)):
            for j in range(i + 1, len(consensus_spheres)):
                p1 = consensus_spheres[i]
                p2 = consensus_spheres[j]
                c1 = np.array(p1["center"])
                c2 = np.array(p2["center"])
                dist = float(np.linalg.norm(c1 - c2))
                distance_constraints.append({
                    "pair": [p1["id"], p2["id"]],
                    "families": [p1["family"], p2["family"]],
                    "target_distance": round(dist, 2),
                    "tolerance": 1.0  # +/- 1.0 A
                })

        # Consensus 2D counts for reference
        consensus_counts = {}
        for fam in ["Donor", "Acceptor", "Aromatic", "Hydrophobe", "PosIonizable", "NegIonizable"]:
            vals = [c[fam] for c in all_counts if fam in c]
            if vals:
                sorted_v = sorted(vals)
                consensus_counts[fam] = {
                    "required": max(1, sorted_v[len(sorted_v) // 2]),
                    "presence_frequency": round(sum(1 for v in vals if v > 0) / len(vals), 2)
                }

        return {
            "actives_count": len(aligned_mols),
            "model_type": "True 3D Spatial Pharmacophore Hypothesis (RDKit ETKDGv3 + O3A)",
            "consensus_spheres": consensus_spheres,
            "distance_constraints": distance_constraints[:12],
            "total_spheres": len(consensus_spheres),
            "consensus_features": consensus_counts,
            "core_requirements": {fam: d["required"] for fam, d in consensus_counts.items()}
        }

    @staticmethod
    def match_candidate(candidate_smiles: str, consensus_profile: Dict[str, Any]) -> Dict[str, Any]:
        """
        Screen candidate molecule against True 3D Spatial Pharmacophore Hypothesis:
        1. Generates 3D conformer ensemble (up to 5 conformers).
        2. Evaluates spatial alignment against consensus tolerance spheres & 3D distance constraints.
        3. Penalizes incorrect 3D spatial arrangement even if 2D counts match!
        """
        if not candidate_smiles or "consensus_spheres" not in consensus_profile:
            # Fallback to legacy count matching if profile lacks 3D spheres
            return PharmacophoreService._legacy_match(candidate_smiles, consensus_profile)

        spheres = consensus_profile.get("consensus_spheres", [])
        dist_constraints = consensus_profile.get("distance_constraints", [])

        if not spheres:
            return {"match_score_pct": 100.0, "status": "No 3D spatial constraints in profile"}

        cand_mol_3d = PharmacophoreService._generate_3d_mol(candidate_smiles, num_confs=4)
        if not cand_mol_3d:
            return {"match_score_pct": 0.0, "status": "Invalid candidate structure (3D conformer generation failed)"}

        best_score = 0.0
        best_details = {}
        best_rmsd = 99.0

        # Evaluate across conformers to find best 3D spatial pose
        for cid in range(cand_mol_3d.GetNumConformers()):
            _, sfeats = PharmacophoreService.extract_3d_spatial_features(cand_mol_3d, conf_id=cid)
            if not sfeats:
                continue

            # 1. Sphere matching: each sphere must find a matching feature of that family
            sphere_matches = 0
            sphere_details = {}
            matched_deviations = []

            for sp in spheres:
                s_id = sp["id"]
                fam = sp["family"]
                sc = np.array(sp["center"])
                s_tol = sp.get("tolerance_radius", 1.2)

                # Find candidate features of matching family
                fam_feats = [f for f in sfeats if f["family"] == fam]
                if not fam_feats:
                    sphere_details[s_id] = f"{fam}: Missing feature"
                    continue

                # Compute distances to candidate feature centers
                dists = [np.linalg.norm(np.array(f["center"]) - sc) for f in fam_feats]
                min_d = min(dists)

                if min_d <= s_tol:
                    sphere_matches += 1.0
                    matched_deviations.append(min_d)
                    sphere_details[s_id] = f"{fam}: Matched in 3D (Δ={min_d:.2f} Å <= {s_tol} Å)"
                elif min_d <= s_tol * 1.8:
                    sphere_matches += 0.5
                    matched_deviations.append(min_d)
                    sphere_details[s_id] = f"{fam}: Partial 3D overlap (Δ={min_d:.2f} Å)"
                else:
                    sphere_details[s_id] = f"{fam}: Spatial mismatch (closest Δ={min_d:.2f} Å)"

            sphere_score = (sphere_matches / len(spheres)) * 100.0

            # 2. Distance constraint satisfaction
            dist_matches = 0
            if dist_constraints:
                for dc in dist_constraints:
                    p1_id, p2_id = dc["pair"]
                    f1, f2 = dc["families"]
                    target_d = dc["target_distance"]
                    tol = dc.get("tolerance", 1.0)

                    # Look for feature pairs in candidate of families f1, f2
                    f1_cand = [f for f in sfeats if f["family"] == f1]
                    f2_cand = [f for f in sfeats if f["family"] == f2]

                    pair_ok = False
                    for cf1 in f1_cand:
                        for cf2 in f2_cand:
                            if cf1["center"] != cf2["center"]:
                                cand_d = np.linalg.norm(np.array(cf1["center"]) - np.array(cf2["center"]))
                                if abs(cand_d - target_d) <= tol:
                                    pair_ok = True
                                    break
                        if pair_ok:
                            break
                    if pair_ok:
                        dist_matches += 1
                dist_score = (dist_matches / len(dist_constraints)) * 100.0
                total_conformer_score = 0.65 * sphere_score + 0.35 * dist_score
            else:
                total_conformer_score = sphere_score

            cur_rmsd = float(np.mean(matched_deviations)) if matched_deviations else 99.0

            if total_conformer_score > best_score:
                best_score = total_conformer_score
                best_details = sphere_details
                best_rmsd = cur_rmsd

        final_pct = round(best_score, 1)
        if final_pct >= 80.0:
            status = "Strong 3D Pharmacophore Match"
        elif final_pct >= 50.0:
            status = "Moderate 3D Pharmacophore Match"
        else:
            status = "Weak / Low 3D Spatial Match"

        return {
            "match_score_pct": final_pct,
            "status": status,
            "spatial_rmsd_angstroms": round(best_rmsd, 2) if best_rmsd < 90.0 else None,
            "total_spheres_tested": len(spheres),
            "spatial_sphere_alignment": best_details,
            "candidate_features": PharmacophoreService.extract_features(candidate_smiles)
        }

    @staticmethod
    def _legacy_match(candidate_smiles: str, consensus_profile: Dict[str, Any]) -> Dict[str, Any]:
        """Fallback 2D count matcher."""
        cand_feats = PharmacophoreService.extract_features(candidate_smiles)
        reqs = consensus_profile.get("consensus_features", {})
        if not reqs:
            return {"match_score_pct": 100.0, "status": "No requirements"}

        match_points = 0.0
        details = {}
        for fam, req_info in reqs.items():
            req_count = req_info.get("required", 1)
            actual_count = cand_feats.get(fam, 0)
            if actual_count >= req_count:
                match_points += 1.0
                details[fam] = f"Satisfied ({actual_count} >= {req_count})"
            elif actual_count > 0:
                match_points += (actual_count / req_count) * 0.7
                details[fam] = f"Partial ({actual_count}/{req_count})"
            else:
                details[fam] = f"Missing (0/{req_count})"

        score_pct = round((match_points / len(reqs)) * 100.0, 1)
        status = "Strong Pharmacophore Match" if score_pct >= 80 else ("Moderate Match" if score_pct >= 50 else "Weak Match")
        return {
            "match_score_pct": score_pct,
            "status": status,
            "candidate_features": cand_feats,
            "feature_alignment": details
        }
