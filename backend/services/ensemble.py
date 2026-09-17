import json
import math
import urllib.request
import urllib.parse
from pathlib import Path
from typing import List, Dict, Any, Optional

from backend.config import CACHE_DIR, UNIPROT_BASE_URL
from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine
from backend.services.bioactivity import BioactivityService

HEADERS = {"User-Agent": "Bindora-Research-Tool/1.0 (academic; +https://github.com/Swelo-ui/Bindora)"}

class EnsembleDockingService:
    """Service to discover multiple target structures via UniProt and run ensemble cross-conformation docking."""

    @staticmethod
    def fetch_ensemble_structures(uniprot_accession: str, limit: int = 6) -> List[Dict[str, Any]]:
        """Query UniProt for deposited PDB structures mapped to this protein accession."""
        acc = uniprot_accession.strip().upper()
        if not acc:
            return []

        if len(acc) == 4:
            try:
                from backend.services.fetcher import StructureFetcher
                u_info = StructureFetcher.fetch_uniprot_by_pdb_id(acc)
                if u_info and u_info.get("accession"):
                    acc = u_info["accession"].strip().upper()
            except Exception:
                pass

        cache_file = CACHE_DIR / f"ensemble_uniprot_{acc}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    if cached:
                        return cached[:limit]
            except Exception:
                pass

        url = f"{UNIPROT_BASE_URL}/{acc}.json"
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=12) as resp:
                data = json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            print(f"[ENSEMBLE FETCH WARNING] UniProt lookup failed for {acc}: {e}")
            return []

        structures = []
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "PDB":
                props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
                res_str = props.get("Resolution", "N/A")
                res_val = 99.0
                try:
                    res_val = float(res_str.replace("A", "").strip())
                except Exception:
                    pass

                structures.append({
                    "pdb_id": xref.get("id"),
                    "method": props.get("Method", "Unknown"),
                    "resolution": res_str,
                    "resolution_val": res_val,
                    "chains": props.get("Chains", "")
                })

        # Sort primarily by resolution (highest quality crystal structures first)
        structures.sort(key=lambda x: x["resolution_val"])

        cleaned = []
        for s in structures:
            cleaned.append({
                "pdb_id": s["pdb_id"],
                "method": s["method"],
                "resolution": s["resolution"],
                "chains": s["chains"]
            })

        try:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(cleaned, f, indent=2)
        except Exception:
            pass

        return cleaned[:limit]

    @staticmethod
    def run_ensemble_docking(
        pdb_ids: List[str],
        ligand_smiles: str,
        exhaustiveness: int = 4
    ) -> Dict[str, Any]:
        """Dock a fixed candidate ligand across multiple macromolecular conformations of the target."""
        if not pdb_ids or not ligand_smiles:
            return {"error": "Missing pdb_ids or ligand_smiles for ensemble docking"}

        # Prepare ligand once
        lig_prep = DockingEngine.prepare_ligand(ligand_smiles)
        lig_pdbqt = lig_prep.get("pdbqt_text")
        heavy_atoms = lig_prep.get("heavy_atom_count", 20)
        mw = 300.0
        if "adme" in lig_prep and "physicochemical" in lig_prep["adme"]:
            mw = float(lig_prep["adme"]["physicochemical"].get("molecular_weight", {}).get("value", 300.0))

        results = []
        affinities = []

        for pdb_id in pdb_ids:
            pdb_id = pdb_id.strip().upper()
            try:
                meta = StructureFetcher.fetch_rcsb_pdb(pdb_id)
                if not meta or not meta.get("pdb_content"):
                    results.append({
                        "pdb_id": pdb_id,
                        "title": f"Structure {pdb_id}",
                        "resolution": "N/A",
                        "status": "Failed (PDB unavailable)",
                        "affinity_kcal": None
                    })
                    continue

                rec_prep = DockingEngine.prepare_receptor(meta["pdb_content"])
                pocket = rec_prep.get("detected_pocket", {})
                center = pocket.get("center", {"x": 0.0, "y": 0.0, "z": 0.0})
                size = pocket.get("size", {"x": 22.0, "y": 22.0, "z": 22.0})

                poses = DockingEngine.run_docking(
                    receptor_pdbqt=rec_prep["pdbqt_text"],
                    ligand_pdbqt=lig_pdbqt,
                    center=center,
                    size=size,
                    exhaustiveness=exhaustiveness,
                    num_modes=1
                )

                if not poses:
                    results.append({
                        "pdb_id": pdb_id,
                        "title": meta.get("title", f"Structure {pdb_id}"),
                        "resolution": meta.get("resolution", "N/A"),
                        "status": "No pose found",
                        "affinity_kcal": None
                    })
                    continue

                best = poses[0]
                aff = best["affinity_kcal"]
                affinities.append(aff)

                # Vinardo rescoring
                vinardo = DockingEngine.score_pose_vinardo(
                    rec_prep["pdbqt_text"], best["pdbqt_content"], center, size
                )

                # Interactions
                contacts = DockingEngine.analyze_interactions(rec_prep["cleaned_pdb"], best["pdbqt_content"])
                hbond_cnt = contacts.get("total_hbond_count", 0)

                # Thermodynamics
                thermo = BioactivityService.calculate_thermodynamics(aff, heavy_atoms, mw)

                results.append({
                    "pdb_id": pdb_id,
                    "title": meta.get("title", f"Structure {pdb_id}")[:60],
                    "resolution": meta.get("resolution", "N/A"),
                    "affinity_kcal": round(aff, 2),
                    "vinardo_score": vinardo,
                    "kd_nanomolar": thermo.get("kd_nanomolar"),
                    "ligand_efficiency": thermo.get("ligand_efficiency"),
                    "hbond_count": hbond_cnt,
                    "status": "Completed",
                    "top_pose_pdbqt": best["pdbqt_content"]
                })
            except Exception as ex:
                results.append({
                    "pdb_id": pdb_id,
                    "title": f"Structure {pdb_id}",
                    "resolution": "N/A",
                    "status": f"Error: {str(ex)[:40]}",
                    "affinity_kcal": None
                })

        # Calculate ensemble summary statistics
        summary = {
            "total_structures": len(pdb_ids),
            "successful_runs": len(affinities),
            "mean_affinity": round(sum(affinities) / len(affinities), 2) if affinities else None,
            "std_affinity": round(math.sqrt(sum((a - (sum(affinities)/len(affinities)))**2 for a in affinities) / len(affinities)), 2) if len(affinities) > 1 else 0.0,
            "best_affinity": min(affinities) if affinities else None,
            "spread_kcal": round(max(affinities) - min(affinities), 2) if affinities else None
        }

        # Sort results: successful runs first, ordered by binding affinity (most negative first)
        results.sort(key=lambda r: (0 if r.get("affinity_kcal") is not None else 1, r.get("affinity_kcal") or 999.0))

        return {
            "structures": results,
            "summary": summary
        }
