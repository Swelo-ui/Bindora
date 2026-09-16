from typing import List, Dict, Any
from rdkit import Chem
from backend.services.docking import DockingEngine
from backend.services.adme import ADMEProfiler
from backend.services.bioactivity import BioactivityService

class BatchScreeningService:
    """Service to execute multi-ligand batch docking, SMILES sanitization, and comparative ranking against a single target."""

    @staticmethod
    def run_batch(
        receptor_pdbqt: str,
        receptor_pdb: str,
        pocket_center: Dict[str, float],
        pocket_size: Dict[str, float],
        ligand_list: List[Dict[str, str]], # [{"name": "Aspirin", "smiles": "CC(=O)..."}, ...]
        exhaustiveness: int = 4
    ) -> List[Dict[str, Any]]:
        """Dock multiple candidate ligands against the same receptor pocket and rank them with consensus scoring."""
        valid_results = []
        failed_results = []

        for item in ligand_list:
            name = item.get("name", "Unknown").strip()
            smiles = item.get("smiles", "").strip()
            if not smiles:
                continue

            # 1. Strict SMILES Sanitization Pre-flight
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                failed_results.append({
                    "name": name,
                    "smiles": smiles,
                    "valid": False,
                    "status": "Failed: Syntax Error",
                    "error": "Malformed SMILES string (unparseable chemical syntax)",
                    "affinity_kcal": 0.0,
                    "consensus_score": 0.0,
                    "rank": "—"
                })
                continue

            try:
                Chem.SanitizeMol(mol)
            except Exception as e:
                failed_results.append({
                    "name": name,
                    "smiles": smiles,
                    "valid": False,
                    "status": "Failed: Sanitization Error",
                    "error": f"Invalid valence or aromaticity: {str(e)}",
                    "affinity_kcal": 0.0,
                    "consensus_score": 0.0,
                    "rank": "—"
                })
                continue

            try:
                # 2. ADME profile
                adme = ADMEProfiler.calculate_adme(smiles)
                if "error" in adme:
                    failed_results.append({
                        "name": name,
                        "smiles": smiles,
                        "valid": False,
                        "status": "Failed: ADME Calculation",
                        "error": adme.get("error", "Cheminformatics descriptor computation failed"),
                        "affinity_kcal": 0.0,
                        "consensus_score": 0.0,
                        "rank": "—"
                    })
                    continue

                phys = adme.get("physicochemical", {})
                mw = phys.get("molecular_weight", {}).get("value", 0.0)
                heavy_atoms = phys.get("heavy_atoms", {}).get("value", 1)

                # 3. Prepare ligand 3D conformer & PDBQT
                lig_prep = DockingEngine.prepare_ligand(smiles)

                # 4. Docking run
                poses = DockingEngine.run_docking(
                    receptor_pdbqt,
                    lig_prep["pdbqt_text"],
                    pocket_center,
                    pocket_size,
                    exhaustiveness=exhaustiveness,
                    num_modes=3
                )

                if not poses:
                    failed_results.append({
                        "name": name,
                        "smiles": smiles,
                        "valid": False,
                        "status": "Failed: No Poses",
                        "error": "Vina returned no docked conformations within bounding volume",
                        "affinity_kcal": 0.0,
                        "consensus_score": 0.0,
                        "rank": "—"
                    })
                    continue

                best_pose = poses[0]
                affinity = best_pose["affinity_kcal"]

                # 5. Thermodynamic conversion
                thermo = BioactivityService.calculate_thermodynamics(affinity, heavy_atoms, mw)
                le = thermo["ligand_efficiency"]["value"]

                # 6. Intermolecular contacts for top pose
                contacts = DockingEngine.analyze_interactions(receptor_pdb, best_pose["pdbqt_content"])
                hbond_cnt = contacts.get("total_hbond_count", 0)

                # 7. Multi-engine Rescoring & Consensus Confidence Evaluation
                vinardo_aff = best_pose.get("vinardo_affinity_kcal")
                if vinardo_aff is None:
                    vinardo_aff = DockingEngine.score_pose_vinardo(receptor_pdbqt, best_pose["pdbqt_content"])
                
                gnina_info = best_pose.get("gnina") or DockingEngine.score_pose_gnina(receptor_pdbqt, best_pose["pdbqt_content"])
                gnina_score = gnina_info.get("cnn_score")

                # Balances: |ΔG| (50%), Ligand Efficiency (scaled x10), and H-bond contact reward (+0.25 each, up to +1.0)
                consensus = round(abs(affinity) * 0.5 + (le * 10.0) + min(hbond_cnt, 4) * 0.25, 2)

                valid_results.append({
                    "name": name,
                    "smiles": smiles,
                    "valid": True,
                    "status": "Screened",
                    "affinity_kcal": affinity,
                    "vinardo_affinity_kcal": vinardo_aff,
                    "gnina_score": gnina_score,
                    "consensus_confidence": "Pending",
                    "consensus_score": consensus,
                    "theoretical_kd_nm": thermo["theoretical_kd_nm"],
                    "ligand_efficiency": le,
                    "potency_class": thermo["potency_class"],
                    "is_weak_binder": thermo.get("is_weak_binder", False),
                    "hbond_count": hbond_cnt,
                    "hydrophobic_count": contacts["total_hydrophobic_count"],
                    "lipinski_status": adme["drug_likeness"]["lipinski"]["status"],
                    "lipinski_violations": adme["drug_likeness"]["lipinski"]["violations_count"],
                    "gi_absorption": adme["pharmacokinetics"]["gi_absorption"]["level"],
                    "pains_count": adme["medicinal_chemistry_safety"]["pains_alerts"]["count"],
                    "top_pose_pdb": best_pose["pdb_block"]
                })
            except Exception as e:
                failed_results.append({
                    "name": name,
                    "smiles": smiles,
                    "valid": False,
                    "status": "Failed: Pipeline Exception",
                    "error": str(e),
                    "affinity_kcal": 0.0,
                    "consensus_score": 0.0,
                    "rank": "—"
                })

        # Multi-Engine Consensus Ranking Calibration:
        # 1. Sort by Vina affinity to assign Vina ranks
        valid_results.sort(key=lambda x: x["affinity_kcal"])
        for idx, item in enumerate(valid_results):
            item["vina_rank"] = idx + 1
            item["rank"] = idx + 1

        # 2. Sort copy by Vinardo affinity to assign Vinardo ranks
        vinardo_valid = [item for item in valid_results if item.get("vinardo_affinity_kcal") is not None]
        vinardo_valid.sort(key=lambda x: x["vinardo_affinity_kcal"])
        vinardo_ranks = {item["name"]: idx + 1 for idx, item in enumerate(vinardo_valid)}

        # 3. Evaluate Concrete Consensus Confidence:
        # - High-Confidence: Rank difference <= 1 AND energy difference <= 3.0 kcal/mol
        # - Moderate-Confidence: Rank difference <= 2 AND energy difference <= 4.5 kcal/mol
        # - Divergent — Low Confidence: Rank difference > 2 OR energy difference > 4.5 kcal/mol
        for item in valid_results:
            vinardo_rank = vinardo_ranks.get(item["name"])
            vinardo_aff = item.get("vinardo_affinity_kcal")
            if vinardo_rank is not None and vinardo_aff is not None:
                rank_delta = abs(item["vina_rank"] - vinardo_rank)
                energy_delta = round(abs(item["affinity_kcal"] - vinardo_aff), 2)
                item["rank_delta"] = rank_delta
                item["energy_delta_kcal"] = energy_delta
                item["consensus_criteria"] = f"|ΔRank|={rank_delta}, |ΔΔG|={energy_delta} kcal/mol"

                if rank_delta <= 1 and energy_delta <= 3.0:
                    item["consensus_confidence"] = "High-Confidence"
                elif rank_delta <= 2 and energy_delta <= 4.5:
                    item["consensus_confidence"] = "Moderate-Confidence"
                else:
                    item["consensus_confidence"] = "Divergent — Low Confidence"
            else:
                item["consensus_confidence"] = "Standard (Single Engine)"
                item["consensus_criteria"] = "Single Engine"

        # Combine valid and failed (failed listed at bottom with rank "—")
        return valid_results + failed_results
