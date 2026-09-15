from typing import List, Dict, Any
from backend.services.docking import DockingEngine
from backend.services.adme import ADMEProfiler
from backend.services.bioactivity import BioactivityService

class BatchScreeningService:
    """Service to execute multi-ligand batch docking and comparative ranking against a single target."""

    @staticmethod
    def run_batch(
        receptor_pdbqt: str,
        receptor_pdb: str,
        pocket_center: Dict[str, float],
        pocket_size: Dict[str, float],
        ligand_list: List[Dict[str, str]], # [{"name": "Aspirin", "smiles": "CC(=O)..."}, ...]
        exhaustiveness: int = 4
    ) -> List[Dict[str, Any]]:
        """Dock multiple candidate ligands against the same receptor pocket and rank them."""
        results = []

        for item in ligand_list:
            name = item.get("name", "Unknown")
            smiles = item.get("smiles", "")
            if not smiles:
                continue

            try:
                # 1. ADME profile
                adme = ADMEProfiler.calculate_adme(smiles)
                if "error" in adme:
                    continue

                phys = adme.get("physicochemical", {})
                mw = phys.get("molecular_weight", {}).get("value", 0.0)
                heavy_atoms = phys.get("heavy_atoms", {}).get("value", 1)

                # 2. Prepare ligand
                lig_prep = DockingEngine.prepare_ligand(smiles)

                # 3. Docking run
                poses = DockingEngine.run_docking(
                    receptor_pdbqt,
                    lig_prep["pdbqt_text"],
                    pocket_center,
                    pocket_size,
                    exhaustiveness=exhaustiveness,
                    num_modes=3
                )

                if not poses:
                    continue

                best_pose = poses[0]
                affinity = best_pose["affinity_kcal"]

                # 4. Thermodynamic conversion
                thermo = BioactivityService.calculate_thermodynamics(affinity, heavy_atoms, mw)

                # 5. Intermolecular contacts for top pose
                contacts = DockingEngine.analyze_interactions(receptor_pdb, best_pose["pdbqt_content"])

                results.append({
                    "name": name,
                    "smiles": smiles,
                    "affinity_kcal": affinity,
                    "theoretical_kd_nm": thermo["theoretical_kd_nm"],
                    "ligand_efficiency": thermo["ligand_efficiency"]["value"],
                    "potency_class": thermo["potency_class"],
                    "hbond_count": contacts["total_hbond_count"],
                    "hydrophobic_count": contacts["total_hydrophobic_count"],
                    "lipinski_status": adme["drug_likeness"]["lipinski"]["status"],
                    "lipinski_violations": adme["drug_likeness"]["lipinski"]["violations_count"],
                    "gi_absorption": adme["pharmacokinetics"]["gi_absorption"]["level"],
                    "pains_count": adme["medicinal_chemistry_safety"]["pains_alerts"]["count"],
                    "top_pose_pdb": best_pose["pdb_block"]
                })
            except Exception as e:
                # Log error and skip failed candidate
                continue

        # Sort by binding affinity (most negative is strongest binding)
        results.sort(key=lambda x: x["affinity_kcal"])

        # Add rank
        for idx, item in enumerate(results):
            item["rank"] = idx + 1

        return results
