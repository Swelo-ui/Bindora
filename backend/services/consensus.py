"""
Multi-Engine Deep Learning & Consensus Ranking Matrix for Bindora Dock.

Combines empirical, physics-based, and deep learning scoring metrics into a unified consensus matrix:
1. AutoDock Vina Empirical Score
2. Vinardo Empirical Free Energy
3. MM-GBSA Continuum Physics Solvation Free Energy
4. GNINA / CNN Shape Complementarity & Affinity
5. Ligand Intramolecular Strain Penalty Filtering
6. Targeted Covalent Enthalpy Bonus (when applicable)
"""

import math
from typing import Dict, Any, List, Optional, Tuple


class ConsensusScoringService:
    """
    Consensus Scoring and Multi-Engine Rank Aggregator.
    Eliminates false-positive decoys by requiring concordant agreement between
    empirical potentials (Vina/Vinardo), physics mechanics (MM-GBSA), and deep learning (CNN).
    """

    DEFAULT_WEIGHTS = {
        "vina": 0.30,
        "vinardo": 0.20,
        "mmgbsa": 0.30,
        "cnn": 0.20
    }

    STRAIN_PENALTY_THRESHOLD = 8.0   # kcal/mol: strain above this begins penalty
    HIGH_STRAIN_CUTOFF = 15.0        # kcal/mol: flags high-strain decoy (> 15.0 kcal/mol)

    @classmethod
    def compute_pose_consensus(
        cls,
        poses: List[Dict[str, Any]],
        weights: Optional[Dict[str, float]] = None,
        covalent_bonus_kcal: float = 0.0
    ) -> List[Dict[str, Any]]:
        """
        Evaluate and re-rank a collection of docked poses using multi-engine consensus matrix.

        Parameters:
            poses: List of pose dictionaries with 'affinity_kcal', optional 'vinardo_affinity_kcal',
                   'mmgbsa_delta_g_kcal', 'gnina', and 'ligand_strain_kcal'.
            weights: Custom weighting dict for scoring terms.
            covalent_bonus_kcal: Optional covalent enthalpy bonus from CovalentDockingService.

        Returns:
            Re-ranked list of poses augmented with consensus scores, ranks, and confidence ratings.
        """
        if not poses:
            return []

        w = dict(cls.DEFAULT_WEIGHTS)
        if weights:
            w.update(weights)

        n = len(poses)

        # 1. Extract values for each metric
        vina_vals = []
        vinardo_vals = []
        mmgbsa_vals = []
        cnn_vals = []
        strain_vals = []

        for p in poses:
            # Vina affinity (kcal/mol, negative is better)
            v_aff = p.get("affinity_kcal")
            vina_vals.append(v_aff if v_aff is not None else 0.0)

            # Vinardo affinity
            vin_aff = p.get("vinardo_affinity_kcal")
            vinardo_vals.append(vin_aff if vin_aff is not None else (v_aff if v_aff is not None else 0.0))

            # MM-GBSA delta G
            gbsa = p.get("mmgbsa_delta_g_kcal")
            if gbsa is None and "mmgbsa" in p and isinstance(p["mmgbsa"], dict):
                gbsa = p["mmgbsa"].get("mmgbsa_delta_g_kcal")
            mmgbsa_vals.append(gbsa if gbsa is not None else (v_aff if v_aff is not None else 0.0))

            # CNN affinity or score
            gnina_dict = p.get("gnina", {})
            cnn_score = None
            if isinstance(gnina_dict, dict) and gnina_dict.get("available"):
                cnn_score = gnina_dict.get("cnn_score")
            # If CNN not available, estimate contact-based surrogate from normalized Vina/strain
            if cnn_score is None:
                # Surrogate CNN score in [0.0, 1.0] based on binding efficiency
                s_est = max(0.0, min(1.0, 0.5 - (v_aff or -7.0) / 20.0))
                cnn_score = round(s_est, 3)
            cnn_vals.append(cnn_score)

            # Strain energy
            st = p.get("ligand_strain_kcal")
            strain_vals.append(st if st is not None else 0.0)

        # 2. Compute Z-score normalizations
        def z_score(arr: List[float]) -> List[float]:
            if len(arr) <= 1:
                return [0.0] * len(arr)
            mean = sum(arr) / len(arr)
            std = math.sqrt(sum((x - mean)**2 for x in arr) / len(arr))
            if std < 1e-5:
                return [0.0] * len(arr)
            return [(x - mean) / std for x in arr]

        z_vina = z_score(vina_vals)
        z_vinardo = z_score(vinardo_vals)
        z_mmgbsa = z_score(mmgbsa_vals)
        # For CNN, higher is better, so negate z-score so that negative is better
        z_cnn = [-z for z in z_score(cnn_vals)]

        # 3. Calculate Individual Ranks
        def get_ranks(arr: List[float], ascending: bool = True) -> List[int]:
            indexed = list(enumerate(arr))
            indexed.sort(key=lambda x: x[1], reverse=not ascending)
            ranks = [0] * len(arr)
            for r, (idx, _) in enumerate(indexed, 1):
                ranks[idx] = r
            return ranks

        ranks_vina = get_ranks(vina_vals, ascending=True)
        ranks_vinardo = get_ranks(vinardo_vals, ascending=True)
        ranks_mmgbsa = get_ranks(mmgbsa_vals, ascending=True)
        ranks_cnn = get_ranks(cnn_vals, ascending=False)  # Higher CNN is better

        # 4. Compute composite consensus score for each pose
        augmented_poses = []
        for i in range(n):
            p_copy = dict(poses[i])

            # Strain penalty: zero below threshold, linear above
            st = strain_vals[i]
            strain_pen = 0.0
            if st > cls.STRAIN_PENALTY_THRESHOLD:
                strain_pen = (st - cls.STRAIN_PENALTY_THRESHOLD) * 0.75

            # Per-pose covalent energy bonus
            p_cov_bonus = (
                poses[i].get("covalent_energy_bonus_kcal") or
                poses[i].get("covalent", {}).get("covalent_energy_bonus_kcal") or
                covalent_bonus_kcal
            )

            # Composite Z-score (lower is better)
            comp_z = (
                w["vina"] * z_vina[i] +
                w["vinardo"] * z_vinardo[i] +
                w["mmgbsa"] * z_mmgbsa[i] +
                w["cnn"] * z_cnn[i] +
                strain_pen +
                p_cov_bonus
            )

            # Mean rank (Borda count rank aggregation)
            mean_rank = (ranks_vina[i] + ranks_vinardo[i] + ranks_mmgbsa[i] + ranks_cnn[i]) / 4.0

            # Rank variance (concordance measure)
            ind_ranks = [ranks_vina[i], ranks_vinardo[i], ranks_mmgbsa[i], ranks_cnn[i]]
            rank_spread = max(ind_ranks) - min(ind_ranks)

            # Confidence classification & SBDD Decoy Gate:
            # 1. High Strain Decoy (> 15.0 kcal/mol)
            # 2. Grease-Ball Decoy (MM-GBSA confirmed zero polar contacts & negligible electrostatics)
            # 3. High/Moderate/Discordant Confidence
            is_unfavorable_desolv = mmgbsa_vals[i] > 2.0
            decoy_flag = p_copy.get("decoy_filter_flag") or p_copy.get("mmgbsa", {}).get("decoy_filter_verdict")
            is_grease_decoy = (
                decoy_flag in ("FLAGGED_GREASY_DECOY", "FLAGGED_LIPOPHILIC_AGGREGATOR") or
                bool(p_copy.get("mmgbsa", {}).get("is_grease_ball_decoy", False))
            )

            if is_grease_decoy:
                confidence = "DECOY_GREASE_BALL"
                decoy_why = p_copy.get("mmgbsa", {}).get("decoy_reason") or "Lacks polar active site complementarity (opportunistic grease)"
                conf_desc = f"Pose flagged as false-positive decoy: {decoy_why}"
            elif st > cls.HIGH_STRAIN_CUTOFF:
                confidence = "DECOY_HIGH_STRAIN"
                conf_desc = f"Pose flagged as high-strain decoy ({st:.1f} kcal/mol > {cls.HIGH_STRAIN_CUTOFF} kcal/mol cutoff)."
            elif rank_spread <= 2 and mean_rank <= 2.5:
                confidence = "HIGH_CONFIDENCE"
                conf_desc = "Strong multi-engine concordance across empirical, physics (MM-GBSA), and DL scoring."
            elif rank_spread <= 4 or mean_rank <= 4.0:
                confidence = "MODERATE_CONFIDENCE"
                conf_desc = "Acceptable agreement across scoring engines with moderate rank dispersion."
            else:
                confidence = "DISCORDANT_SCORING"
                conf_desc = "Significant disagreement between empirical Vina and physics/CNN scores; verify solvation contacts."

            p_copy["consensus_score"] = round(comp_z, 3)
            p_copy["borda_mean_rank"] = round(mean_rank, 2)
            p_copy["rank_spread"] = rank_spread
            p_copy["consensus_confidence"] = confidence
            p_copy["consensus_confidence_description"] = conf_desc
            p_copy["individual_ranks"] = {
                "vina": ranks_vina[i],
                "vinardo": ranks_vinardo[i],
                "mmgbsa": ranks_mmgbsa[i],
                "cnn": ranks_cnn[i]
            }
            p_copy["scoring_breakdown"] = {
                "vina_kcal": vina_vals[i],
                "vinardo_kcal": vinardo_vals[i],
                "mmgbsa_kcal": mmgbsa_vals[i],
                "cnn_score": cnn_vals[i],
                "strain_penalty_kcal": round(strain_pen, 2),
                "covalent_bonus_kcal": round(covalent_bonus_kcal, 2)
            }

            augmented_poses.append(p_copy)

        # 5. Sort by consensus score: decoys automatically relegated to bottom
        augmented_poses.sort(key=lambda x: (
            x["consensus_confidence"] in ("DECOY_HIGH_STRAIN", "DECOY_GREASE_BALL"),
            x["consensus_score"]
        ))

        # Assign final consensus ranks
        for final_rank, p in enumerate(augmented_poses, 1):
            p["consensus_rank"] = final_rank

        return augmented_poses
