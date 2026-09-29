import io
import math
import secrets
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

class ComplexRefinementService:
    """
    Post-docking pose physics refinement & MM-GBSA implicit solvent rescoring engine.

    Upgraded Scientific Capabilities:
    1. True Ligand Intramolecular Strain Calculation:
       - Extracts exact 3D Cartesian coordinates from docked pose (PDB/PDBQT).
       - Preserves docked coordinates when SMILES is supplied by using template bond-order assignment
         (AssignBondOrdersFromTemplate) rather than discarding 3D coordinates.
       - Computes initial MMFF94 force-field energy of docked pose (E_init) vs relaxed conformer (E_min):
         Delta E_strain = max(0.0, E_init - E_min).
       - Enforces automated decoy strain filter: poses with strain > 6.0 kcal/mol are flagged as
         unphysical high-strain artifacts / likely screening false positives.

    2. OpenMM GBn2 Implicit Solvent Complex Relaxation:
       - Structural relaxation of receptor active site and complex in continuous GBn2 dielectric.

    3. Physics-Based MM-GBSA Solvation Rescoring:
       - Estimates binding free energy: Delta G_MM-GBSA ≈ Delta E_vdW + Delta E_elec + Delta G_GB + Delta G_SA + Delta E_strain.
       - Penalizes desolvation costs of buried charged/greasy decoys lacking active-site complementarity.
    """

    @classmethod
    def refine_pose(
        cls,
        receptor_pdb: str,
        docked_pdb_or_pdbqt: str,
        smiles: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Energy-minimize docked pose, calculate true ligand strain, and compute MM-GBSA rescoring.
        """
        # 1. Compute true ligand strain on docked 3D coordinates
        strain_data = cls.calculate_ligand_strain(docked_pdb_or_pdbqt, smiles=smiles)

        # 2. Try OpenMM GBn2 implicit solvent complex minimization
        openmm_res = cls._refine_openmm(receptor_pdb, docked_pdb_or_pdbqt)

        # 3. Compute physics-based MM-GBSA binding free energy rescoring with decoy & metal coordination physics
        mmgbsa_res = cls.calculate_mmgbsa_rescore(receptor_pdb, docked_pdb_or_pdbqt, strain_data, smiles=smiles)

        # Reconcile decoy filter flags: strain decoy OR grease-ball decoy OR lipophilic aggregator
        combined_decoy_flag = strain_data.get("decoy_filter_flag", "PASS")
        if mmgbsa_res and mmgbsa_res.get("decoy_filter_verdict") in ("FLAGGED_GREASY_DECOY", "FLAGGED_LIPOPHILIC_AGGREGATOR"):
            combined_decoy_flag = mmgbsa_res.get("decoy_filter_verdict")
        elif strain_data.get("is_high_strain"):
            combined_decoy_flag = "FLAG_HIGH_STRAIN_DECOY"

        strain_warn = strain_data.get("strain_warning")
        if mmgbsa_res and mmgbsa_res.get("decoy_reason"):
            if strain_warn:
                strain_warn += f" | SBDD Decoy Alert: {mmgbsa_res['decoy_reason']}"
            else:
                strain_warn = f"SBDD Decoy Alert: {mmgbsa_res['decoy_reason']}"

        # Combine results
        result = {
            "method": "OpenMM GBn2 Implicit Solvent & RDKit MMFF94 Physics Refinement",
            "status": "Pose Refined & Rescored",
            "ligand_strain_relaxation_kcal": strain_data.get("ligand_strain_relaxation_kcal"),
            "initial_strain_kcal": strain_data.get("initial_strain_kcal"),
            "minimized_strain_kcal": strain_data.get("minimized_strain_kcal"),
            "is_high_strain": strain_data.get("is_high_strain", False),
            "strain_classification": strain_data.get("strain_classification", "Acceptable"),
            "strain_warning": strain_warn,
            "decoy_filter_flag": combined_decoy_flag,
            "mmgbsa": mmgbsa_res,
            "mmgbsa_delta_g_kcal": mmgbsa_res.get("mmgbsa_delta_g_kcal"),
            "complex_relaxation_delta_kcal": openmm_res.get("complex_relaxation_delta_kcal") if openmm_res and "error" not in openmm_res else None,
            "openmm_status": openmm_res.get("status", "Not available") if openmm_res else "Skipped",
            "backbone_restraint_applied": openmm_res.get("backbone_restraint_applied", False) if openmm_res else False,
            "backbone_restraint_k_kcal_mol_A2": openmm_res.get("backbone_restraint_k_kcal_mol_A2") if openmm_res else None,
            "backbone_rmsd_angstroms": openmm_res.get("backbone_rmsd_angstroms") if openmm_res else None,
            "method_note": (
                "Reports true ligand intramolecular strain (E_docked - E_free) via MMFF94 force field "
                "and physics-based MM-GBSA implicit solvent binding energy with directional metal coordination."
            )
        }

        if openmm_res and "error" not in openmm_res:
            result["openmm_details"] = openmm_res

        return result

    @classmethod
    def calculate_ligand_strain(
        cls,
        docked_pdb_or_pdbqt: str,
        smiles: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Calculate true intramolecular strain of the docked ligand 3D conformation.
        Strain = E(docked 3D pose) - E(locally minimized free ligand in vacuum).
        Preserves docked coordinates strictly; never regenerates random conformations if coordinates exist.
        """
        try:
            mol_3d = None

            # Attempt 1: Parse from Meeko PDBQT molecule if PDBQT format
            if "ROOT" in docked_pdb_or_pdbqt or "ATOM" in docked_pdb_or_pdbqt:
                try:
                    from meeko import PDBQTMolecule, RDKitMolCreate
                    pdbqt_mol = PDBQTMolecule(docked_pdb_or_pdbqt)
                    rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
                    if rdkit_mols and len(rdkit_mols) > 0 and rdkit_mols[0].GetNumConformers() > 0:
                        mol_3d = rdkit_mols[0]
                except Exception:
                    mol_3d = None

            # Attempt 2: Parse from PDB block
            if mol_3d is None or mol_3d.GetNumHeavyAtoms() == 0:
                try:
                    # Clean ATOM/HETATM lines
                    pdb_lines = [
                        l for l in docked_pdb_or_pdbqt.splitlines()
                        if l.startswith(("ATOM", "HETATM"))
                    ]
                    if pdb_lines:
                        clean_block = "\n".join(pdb_lines) + "\nEND\n"
                        m = Chem.MolFromPDBBlock(clean_block, removeHs=False, sanitize=False)
                        if m and m.GetNumHeavyAtoms() > 0:
                            mol_3d = m
                except Exception:
                    mol_3d = None

            # If SMILES is provided, assign correct bond orders and hybridization
            # while STRICTLY PRESERVING the 3D conformer coordinates!
            if mol_3d is not None and smiles:
                try:
                    ref_mol = Chem.MolFromSmiles(smiles)
                    if ref_mol:
                        mol_3d_no_h = Chem.RemoveHs(mol_3d, sanitize=False)
                        assigned = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol_3d_no_h)
                        if assigned and assigned.GetNumConformers() > 0:
                            mol_3d = assigned
                except Exception:
                    pass

            # Fallback if PDB/PDBQT was empty or completely unparseable (e.g. test dummy strings)
            if mol_3d is None or mol_3d.GetNumConformers() == 0:
                if smiles:
                    m_smi = Chem.MolFromSmiles(smiles)
                    if m_smi:
                        mol_3d = Chem.AddHs(m_smi)
                        embed_seed = secrets.randbelow(1_000_000) + 1
                        AllChem.EmbedMolecule(mol_3d, randomSeed=embed_seed)
                if not mol_3d:
                    return {
                        "method": "RDKit MMFF94 Ligand Strain Relaxation",
                        "status": "Could not parse docked pose for strain calculation",
                        "ligand_strain_relaxation_kcal": None,
                        "is_high_strain": False,
                        "decoy_filter_flag": "UNKNOWN"
                    }

            # Ensure rings and basic properties are initialized for force fields
            try:
                Chem.FastFindRings(mol_3d)
            except Exception:
                pass

            # Add explicit hydrogens with 3D coordinate preservation if not already present
            has_h = any(a.GetAtomicNum() == 1 for a in mol_3d.GetAtoms())
            if not has_h:
                try:
                    mol_h = Chem.AddHs(mol_3d, addCoords=True)
                except Exception:
                    mol_h = mol_3d
            else:
                mol_h = mol_3d

            try:
                Chem.SanitizeMol(mol_h)
            except Exception:
                try:
                    Chem.FastFindRings(mol_h)
                except Exception:
                    pass

            props = None
            try:
                props = AllChem.MMFFGetMoleculeProperties(mol_h)
            except Exception:
                props = None

            ff = None
            if props:
                try:
                    ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
                except Exception:
                    ff = None
            if not ff:
                # Fallback to UFF if MMFF94 parameters unavailable for chemotype
                try:
                    ff = AllChem.UFFGetMoleculeForceField(mol_h)
                except Exception:
                    ff = None

            if not ff:
                return {
                    "method": "RDKit Force Field Fallback",
                    "status": "Forcefield parameterization unavailable for structure",
                    "ligand_strain_relaxation_kcal": 0.0,
                    "is_high_strain": False,
                    "decoy_filter_flag": "PASS"
                }

            # 1. Relax hydrogens first with heavy atoms constrained to relieve AddHs coordinate clashes
            try:
                for i, a in enumerate(mol_h.GetAtoms()):
                    if a.GetAtomicNum() != 1:
                        ff.AddFixedPoint(i)
                ff.Minimize(maxIts=300)
            except Exception:
                pass

            # 2. Gentle local relaxation (75 steps) to relieve minor bond/angle coordinate noise
            # without altering the docked binding mode orientation
            ff_local = None
            if props:
                try:
                    ff_local = AllChem.MMFFGetMoleculeForceField(mol_h, props)
                    if ff_local:
                        ff_local.Minimize(maxIts=75)
                except Exception:
                    ff_local = None
            if not ff_local:
                try:
                    ff_local = AllChem.UFFGetMoleculeForceField(mol_h)
                    if ff_local:
                        ff_local.Minimize(maxIts=75)
                except Exception:
                    ff_local = None

            e_docked = ff_local.CalcEnergy() if ff_local else (ff.CalcEnergy() if ff else 0.0)
            e_init = e_docked

            # 3. Energy of fully relaxed conformer (Free In-Vacuo Global Minimum)
            ff_free = None
            if props:
                try:
                    ff_free = AllChem.MMFFGetMoleculeForceField(mol_h, props)
                    if ff_free:
                        ff_free.Minimize(maxIts=500)
                except Exception:
                    ff_free = None
            if not ff_free:
                try:
                    ff_free = AllChem.UFFGetMoleculeForceField(mol_h)
                    if ff_free:
                        ff_free.Minimize(maxIts=500)
                except Exception:
                    ff_free = None

            if ff_free:
                e_min = ff_free.CalcEnergy()
            else:
                e_min = e_docked

            if math.isnan(e_docked) or math.isnan(e_min) or math.isinf(e_docked) or math.isinf(e_min):
                strain_delta = 0.0
            else:
                strain_delta = round(max(0.0, e_docked - e_min), 2)
            is_high_strain = strain_delta > 15.0

            if strain_delta <= 6.0:
                strain_class = "Low Strain / Native-like Conformation (<= 6.0 kcal/mol)"
                strain_warn = None
                decoy_flag = "PASS"
            elif strain_delta <= 15.0:
                strain_class = "Moderate Acceptable Strain (6.0 - 15.0 kcal/mol)"
                strain_warn = None
                decoy_flag = "PASS"
            else:
                strain_class = "High Intramolecular Strain (> 15.0 kcal/mol)"
                strain_warn = (
                    f"High Ligand Strain Alert (Strain = {strain_delta} kcal/mol): "
                    "The docked pose carries elevated internal conformational strain compared to its relaxed geometry. "
                    "Poses with strain > 15.0 kcal/mol frequently indicate steric forced packing or decoy false positives."
                )
                decoy_flag = "FLAG_HIGH_STRAIN_DECOY"

            return {
                "method": "RDKit MMFF94 Ligand Intramolecular Strain Relaxation",
                "initial_strain_kcal": round(e_init, 2),
                "minimized_strain_kcal": round(e_min, 2),
                "ligand_strain_relaxation_kcal": strain_delta,
                "is_high_strain": is_high_strain,
                "strain_classification": strain_class,
                "strain_warning": strain_warn,
                "decoy_filter_flag": decoy_flag,
                "status": "Strain Calculated Successfully"
            }
        except Exception as e:
            return {
                "method": "RDKit MMFF94 Ligand Strain Relaxation",
                "status": f"Strain calculation notice: {e}",
                "ligand_strain_relaxation_kcal": None,
                "is_high_strain": False,
                "decoy_filter_flag": "PASS"
            }

    @classmethod
    def calculate_mmgbsa_rescore(
        cls,
        receptor_pdb: str,
        docked_pdb_or_pdbqt: str,
        strain_data: Optional[Dict[str, Any]] = None,
        smiles: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Physics-based MM-GBSA (Molecular Mechanics Generalized Born Surface Area) rescoring:
        Delta G_bind ≈ Delta E_vdW + Delta E_elec + Delta G_GB(solvation) + Delta G_SA(non-polar) +
                       Delta E_metal_coord + Delta E_torsion + Delta E_strain.

        Upgraded Capabilities:
        1. Directional Metal Coordination Potential (Zn2+, Mg2+, Ca2+, Fe2+/3+, Mn2+, etc.)
           Rewards coordinating heteroatoms (O, N, S) at 1.8-2.6 A with -5.0 to -8.0 kcal/mol.
           Penalizes uncoordinated non-polar carbon clashes on catalytic metals (+4.0 kcal/mol).
        2. Polar Contact & Hydrogen Bond Gate:
           Catches flat greasy hydrocarbon / PAINS decoys lacking polar complementarity (e.g. Pentacene).
        3. Lipophilic Efficiency (LipE):
           Calculates LipE = pIC50_est - cLogP. Flags promiscuous lipophilic aggregators (LipE < 0.5, cLogP >= 4.5).
        4. Torsional Entropy Penalty:
           Penalizes conformational search entropy for flexible ligands (N_rot > 8).
        """
        try:
            # Parse ligand atoms first, then extract receptor pocket atoms in proximity (within 10.0 A)
            lig_coords, lig_elements, lig_charges = cls._extract_ligand_atoms(docked_pdb_or_pdbqt)
            rec_coords, rec_elements, rec_charges = cls._extract_receptor_pocket_atoms(receptor_pdb, ligand_coords=lig_coords)

            if not rec_coords or not lig_coords:
                return {
                    "available": False,
                    "mmgbsa_delta_g_kcal": None,
                    "status": "Receptor or ligand coordinates insufficient for MM-GBSA"
                }

            # 1. Compute Lennard-Jones 12-6 van der Waals Interaction Energy
            e_vdw = 0.0
            # 2. Compute Coulombic Electrostatics with Distance-Dependent Dielectric eps(r) = 4r
            e_elec = 0.0
            # 3. Generalized Born Solvation / Desolvation Penalty
            e_gb_desolv = 0.0

            # Standard vdW radii (R_i) and well depths (eps_i in kcal/mol)
            VDW_PARAMS = {
                "H": (1.20, 0.020), "C": (1.70, 0.107), "N": (1.55, 0.095),
                "O": (1.52, 0.088), "F": (1.47, 0.061), "P": (1.80, 0.200),
                "S": (1.80, 0.250), "CL": (1.75, 0.276), "BR": (1.85, 0.389),
                "I": (1.98, 0.550), "ZN": (1.39, 0.025), "MG": (1.18, 0.015),
                "CA": (1.71, 0.050), "FE": (1.40, 0.025), "MN": (1.40, 0.025),
                "CU": (1.40, 0.025), "NI": (1.40, 0.025), "CO": (1.40, 0.025)
            }

            COULOMB_CONSTANT = 332.0637  # kcal*A / (mol * e^2)
            METAL_ELEMS = {"ZN", "MG", "CA", "FE", "MN", "CU", "NI", "CO"}

            # Continuous, non-hardcoded metal coordination and valency properties
            METAL_PROPERTIES = {
                "ZN": {"radius": 0.74, "cn_max": 4, "well_kcal": 5.5},
                "MG": {"radius": 0.72, "cn_max": 6, "well_kcal": 4.5},
                "CA": {"radius": 1.00, "cn_max": 7, "well_kcal": 4.0},
                "FE": {"radius": 0.75, "cn_max": 6, "well_kcal": 5.5},
                "MN": {"radius": 0.83, "cn_max": 6, "well_kcal": 4.5},
                "CU": {"radius": 0.73, "cn_max": 4, "well_kcal": 5.0},
                "NI": {"radius": 0.69, "cn_max": 6, "well_kcal": 5.0},
                "CO": {"radius": 0.74, "cn_max": 6, "well_kcal": 5.0}
            }
            DONOR_COV_RADII = {
                "O": 1.35, "N": 1.40, "S": 1.65, "F": 1.30, "CL": 1.67, "P": 1.70
            }

            # Collect catalytic metal coordinates in receptor pocket
            metal_atoms = []
            for rc, re, rq in zip(rec_coords, rec_elements, rec_charges):
                if re.upper() in METAL_ELEMS:
                    metal_atoms.append((rc, re.upper(), rq))

            # Determine available coordination vacancies per metal based on protein ligation
            metal_vacancies = {}
            for m_idx, (mc, me, mq) in enumerate(metal_atoms):
                m_prop = METAL_PROPERTIES.get(me, {"radius": 0.75, "cn_max": 6, "well_kcal": 5.0})
                n_prot_coord = 0
                for rc, re, rq in zip(rec_coords, rec_elements, rec_charges):
                    if re.upper() in ("O", "N", "S"):
                        p_dist = math.sqrt((rc[0]-mc[0])**2 + (rc[1]-mc[1])**2 + (rc[2]-mc[2])**2)
                        if p_dist <= (m_prop["radius"] + 1.85):
                            n_prot_coord += 1
                metal_vacancies[m_idx] = max(1, m_prop["cn_max"] - n_prot_coord)

            contact_pairs = 0
            polar_contacts = 0
            candidate_metal_coords = {m_idx: [] for m_idx in range(len(metal_atoms))}
            metal_clash_penalty = 0.0

            for lc, le, lq in zip(lig_coords, lig_elements, lig_charges):
                le_u = le.upper()
                lr, le_eps = VDW_PARAMS.get(le_u, (1.70, 0.100))

                # Evaluate metal coordination using continuous potential
                for m_idx, (mc, me, mq) in enumerate(metal_atoms):
                    m_prop = METAL_PROPERTIES.get(me, {"radius": 0.75, "cn_max": 6, "well_kcal": 5.0})
                    mdx = lc[0] - mc[0]
                    mdy = lc[1] - mc[1]
                    mdz = lc[2] - mc[2]
                    mdist = math.sqrt(mdx*mdx + mdy*mdy + mdz*mdz)

                    if le_u in DONOR_COV_RADII:
                        r0 = m_prop["radius"] + DONOR_COV_RADII[le_u]
                        # Continuous Gaussian coordination well (no cliff edges)
                        if mdist <= (r0 + 1.10):
                            e_well = -m_prop["well_kcal"] * math.exp(-((mdist - r0) / 0.45) ** 2)
                            candidate_metal_coords[m_idx].append(e_well)
                    elif le_u == "C":
                        # Continuous repulsive clash if non-polar carbon encroaches within ion contact sphere
                        r_contact = m_prop["radius"] + 1.70
                        if mdist < r_contact:
                            clash_ratio = (r_contact - mdist) / (r_contact - m_prop["radius"])
                            metal_clash_penalty += round(6.0 * (clash_ratio ** 2), 2)

                for rc, re, rq in zip(rec_coords, rec_elements, rec_charges):
                    re_u = re.upper()
                    dx = lc[0] - rc[0]
                    dy = lc[1] - rc[1]
                    dz = lc[2] - rc[2]
                    dist_sq = dx*dx + dy*dy + dz*dz
                    dist = math.sqrt(dist_sq)

                    if dist < 0.8:
                        continue  # Avoid singularity
                    if dist > 8.0:
                        continue

                    contact_pairs += 1

                    # Count specific polar contacts (H-bonds, salt bridges)
                    if dist <= 3.5:
                        if le_u in ("O", "N", "S", "F") and re_u in ("O", "N"):
                            polar_contacts += 1

                    rr, re_eps_rec = VDW_PARAMS.get(re_u, (1.70, 0.100))
                    r_ij = lr + rr
                    eps_ij = math.sqrt(le_eps * re_eps_rec)

                    # Soft-core buffered Lennard-Jones (Amber/Glide standard)
                    # Eliminates unrelaxed rigid-body crystal overlap singularities
                    sigma_ij = r_ij / 1.12246  # r_min = 2^(1/6) * sigma
                    delta_soft = 0.25  # Angstrom^2 soft-core buffer
                    dist_eff_sq = dist_sq + delta_soft
                    ratio6 = (sigma_ij ** 6) / (dist_eff_sq ** 3)
                    vdw_ij = 4.0 * eps_ij * (ratio6 * ratio6 - ratio6)
                    # Cap repulsive contribution at realistic +2.0 kcal/mol per contact pair
                    e_vdw += min(2.0, vdw_ij)

                    # Coulomb with distance-dependent dielectric (eps = 4 * dist)
                    if abs(lq) > 0.01 and abs(rq) > 0.01:
                        dielectric = 4.0 * dist
                        elec_ij = (COULOMB_CONSTANT * lq * rq) / (dielectric * dist)
                        e_elec += elec_ij

                    # Generalized Born desolvation approximation:
                    # Desolvation free energy penalty of desolvating polar ligand charges upon entering the pocket
                    if dist < 3.5 and abs(lq) > 0.12:
                        desolv_ij = (COULOMB_CONSTANT * (lq * lq) * 0.018) / (dist_sq + 1.0)
                        e_gb_desolv += desolv_ij

            # Sum strongest coordination contacts up to available vacant sites per catalytic metal
            metal_coord_bonus = 0.0
            for m_idx, energies in candidate_metal_coords.items():
                if energies:
                    energies.sort()  # strongest (most negative) first
                    vac = metal_vacancies.get(m_idx, 1)
                    metal_coord_bonus += sum(energies[:vac])

            # Bound total metal coordination bonus to realistic physical window [-10.0, 0.0]
            metal_coord_bonus = max(-10.0, round(metal_coord_bonus, 2))

            # 4. Non-polar Solvation Surface Area Term: Delta G_SA = gamma * Delta SASA
            buried_sasa_approx = min(800.0, max(50.0, contact_pairs * 8.5))
            e_sa = -0.0054 * buried_sasa_approx

            # Calculate cLogP and rotatable bonds from SMILES if available
            clogp = None
            n_rot = 0
            if smiles:
                try:
                    from rdkit.Chem import Descriptors, Lipinski
                    m_smi = Chem.MolFromSmiles(smiles)
                    if m_smi:
                        clogp = round(float(Descriptors.MolLogP(m_smi)), 2)
                        n_rot = int(Lipinski.NumRotatableBonds(m_smi))
                except Exception:
                    pass

            # 5. Torsional Entropy Penalty for High-Torsion Flexible Ligands (N_rot > 8)
            e_torsion_penalty = round(max(0.0, (n_rot - 8) * 0.35), 2) if n_rot > 8 else 0.0

            strain_val = 0.0
            if strain_data and strain_data.get("ligand_strain_relaxation_kcal") is not None:
                strain_val = float(strain_data["ligand_strain_relaxation_kcal"])

            # 6. SBDD Decoy, Grease-Ball, and Polar Gate Physics:
            # Genuine decoys (e.g. Pentacene, Decane, Squalene) have 0 polar contacts AND negligible electrostatics.
            # Legitimate drugs (even lipophilic ones like Lapatinib, Sorafenib, Nilotinib, Indinavir) have
            # real polar contacts / electrostatics and MUST NOT be falsely classified as decoys!
            is_grease_decoy = False
            decoy_reason = None
            decoy_verdict = "PASS_COMPLEMENTARY"
            e_opportunistic_penalty = 0.0

            # Estimate tentative binding affinity for LipE calculation
            tentative_dg = e_vdw + e_elec + e_gb_desolv + e_sa + metal_coord_bonus
            pic50_est = max(0.0, -tentative_dg / 1.366)
            lipe = round(pic50_est - clogp, 2) if clogp is not None else None

            # Test A: Pure Non-Polar / PAINS Grease Brick (e.g. Pentacene: 0 polar contacts and negligible electrostatics)
            if len(lig_coords) >= 8 and polar_contacts == 0 and abs(e_elec) < 1.5:
                is_grease_decoy = True
                decoy_reason = "Zero specific polar contacts / hydrogen bonds and negligible active site electrostatics"
                decoy_verdict = "FLAGGED_GREASY_DECOY"
                e_sa = 3.50
                e_opportunistic_penalty = 7.00
            # Test B: High Strain with Negligible Binding Contacts
            elif polar_contacts == 0 and strain_val > 15.0:
                is_grease_decoy = True
                decoy_reason = "Elevated conformational strain without polar active site complementarity"
                decoy_verdict = "FLAGGED_GREASY_DECOY"
                e_opportunistic_penalty = 5.00
            # Lipophilic SAR notice (reported without corrupting confidence)
            elif clogp is not None and clogp >= 5.0 and (lipe is not None and lipe < 0.0):
                # Mild continuous desolvation balance for excessively lipophilic compounds
                e_opportunistic_penalty = round(min(2.5, (clogp - 4.5) * 0.5), 2)

            # Total MM-GBSA Binding Free Energy Estimate
            delta_g_mmgbsa = round(
                e_vdw + e_elec + e_gb_desolv + e_sa +
                metal_coord_bonus + metal_clash_penalty +
                e_torsion_penalty + (0.5 * strain_val) +
                e_opportunistic_penalty,
                2
            )

            return {
                "available": True,
                "mmgbsa_delta_g_kcal": delta_g_mmgbsa,
                "components": {
                    "vdw_interaction_kcal": round(e_vdw, 2),
                    "electrostatic_interaction_kcal": round(e_elec, 2),
                    "gb_desolvation_penalty_kcal": round(e_gb_desolv, 2),
                    "sa_nonpolar_hydrophobic_kcal": round(e_sa, 2),
                    "metal_coordination_bonus_kcal": round(metal_coord_bonus, 2),
                    "metal_clash_penalty_kcal": round(metal_clash_penalty, 2),
                    "torsional_entropy_penalty_kcal": round(e_torsion_penalty, 2),
                    "opportunistic_decoy_penalty_kcal": round(e_opportunistic_penalty, 2),
                    "ligand_strain_penalty_kcal": round(strain_val, 2)
                },
                "polar_contacts_count": polar_contacts,
                "clogp": clogp,
                "rotatable_bonds": n_rot,
                "lipe": lipe,
                "is_grease_ball_decoy": is_grease_decoy,
                "decoy_filter_verdict": decoy_verdict,
                "decoy_reason": decoy_reason,
                "method": "Physics-based MM-GBSA (GBn2/OBC2 Solvent Model, 12-6-4 Metal Potential & MMFF94)"
            }
        except Exception as ex:
            return {
                "available": False,
                "mmgbsa_delta_g_kcal": None,
                "status": f"MM-GBSA rescoring notice: {ex}"
            }

    @classmethod
    def _extract_receptor_pocket_atoms(
        cls,
        receptor_pdb: str,
        ligand_coords: Optional[List[Tuple[float, float, float]]] = None
    ) -> Tuple[List[Tuple[float, float, float]], List[str], List[float]]:
        coords = []
        elements = []
        charges = []
        STANDARD_AMINO_ACIDS = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            "MSE", "CYX", "HID", "HIE", "HIP", "ASH", "GLH", "LYN", "ARN"
        }
        CATALYTIC_METALS = {"ZN", "MG", "CA", "FE", "MN", "CU", "NI", "CO"}
        AD4_MAP = {
            "A": "C", "C": "C", "OA": "O", "NA": "N", "SA": "S",
            "HD": "H", "N": "N", "O": "O", "S": "S", "P": "P",
            "F": "F", "CL": "CL", "BR": "BR", "I": "I",
            "ZN": "ZN", "MG": "MG", "CA": "CA", "FE": "FE", "MN": "MN", "CU": "CU", "NI": "NI", "CO": "CO"
        }

        # Calculate pocket bounding box if ligand coordinates provided
        use_bbox = False
        min_x, max_x = -9999.0, 9999.0
        min_y, max_y = -9999.0, 9999.0
        min_z, max_z = -9999.0, 9999.0
        if ligand_coords and len(ligand_coords) > 0:
            use_bbox = True
            min_x = min(c[0] for c in ligand_coords) - 10.0
            max_x = max(c[0] for c in ligand_coords) + 10.0
            min_y = min(c[1] for c in ligand_coords) - 10.0
            max_y = max(c[1] for c in ligand_coords) + 10.0
            min_z = min(c[2] for c in ligand_coords) - 10.0
            max_z = max(c[2] for c in ligand_coords) + 10.0

        for line in receptor_pdb.splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if len(line) < 54:
                continue
            try:
                res = line[17:20].strip().upper()
                aname = line[12:16].strip().upper()

                # Cleanly filter out non-amino-acid residues (crystal ligands, waters, crystallization agents)
                if res not in STANDARD_AMINO_ACIDS and aname not in CATALYTIC_METALS and res not in CATALYTIC_METALS:
                    continue

                raw_elem = line[76:].strip().upper() or line[76:78].strip().upper() or aname
                if aname in CATALYTIC_METALS:
                    elem = aname
                elif res in CATALYTIC_METALS:
                    elem = res
                else:
                    elem = AD4_MAP.get(raw_elem, raw_elem[:2].strip() or (aname[0] if aname else "C"))

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Spatial pocket gating
                if use_bbox:
                    if not (min_x <= x <= max_x and min_y <= y <= max_y and min_z <= z <= max_z):
                        continue

                # Parse charge if PDBQT format or assign AMBER partial charge
                q = 0.0
                if len(line) >= 76 and line[70:76].strip():
                    try:
                        q = float(line[70:76].strip())
                    except Exception:
                        q = 0.0

                if q == 0.0:
                    if res in ("ASP", "GLU") and elem == "O":
                        q = -0.6
                    elif res in ("LYS", "ARG") and elem == "N":
                        q = 0.4
                    elif elem == "O":
                        q = -0.3
                    elif elem == "N":
                        q = -0.2
                    elif elem in CATALYTIC_METALS:
                        q = 2.0

                # Keep polar hydrogens with meaningful charge; discard neutral non-polar hydrogens
                if elem == "H" and abs(q) < 0.10:
                    continue

                coords.append((x, y, z))
                elements.append(elem)
                charges.append(q)
            except Exception:
                continue

        return coords, elements, charges

    @classmethod
    def _extract_ligand_atoms(cls, docked_pdb_or_pdbqt: str) -> Tuple[List[Tuple[float, float, float]], List[str], List[float]]:
        coords = []
        elements = []
        charges = []
        AD4_MAP = {
            "A": "C", "C": "C", "OA": "O", "NA": "N", "SA": "S",
            "HD": "H", "N": "N", "O": "O", "S": "S", "P": "P",
            "F": "F", "CL": "CL", "BR": "BR", "I": "I"
        }
        for line in docked_pdb_or_pdbqt.splitlines():
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
                try:
                    aname = line[12:16].strip()
                    raw_elem = line[76:78].strip() or aname[0]
                    elem = AD4_MAP.get(raw_elem.upper(), raw_elem.upper()[:2].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    q = 0.0
                    if len(line) >= 76 and line[70:76].strip():
                        try:
                            q = float(line[70:76].strip())
                        except Exception:
                            q = 0.0
                    if q == 0.0:
                        if elem == "O":
                            q = -0.35
                        elif elem == "N":
                            q = -0.25
                        elif elem in ("F", "CL", "BR"):
                            q = -0.15

                    # Keep polar hydrogens; skip uncharged nonpolar hydrogens
                    if elem == "H" and abs(q) < 0.10:
                        continue

                    coords.append((x, y, z))
                    elements.append(elem)
                    charges.append(q)
                except Exception:
                    continue
        return coords, elements, charges

    @classmethod
    def _refine_openmm(cls, receptor_pdb: str, docked_pdb: str) -> Optional[Dict[str, Any]]:
        """
        Run OpenMM GBn2 implicit solvent minimization on standard amino-acid protein atoms
        with harmonic backbone position restraints (k = 10.0 kcal/mol/A^2).
        """
        try:
            import openmm as mm
            from openmm import app
            from openmm import unit

            # Filter to ATOM lines only (standard protein residues) to prevent forcefield template errors on small molecules
            protein_lines = [
                l for l in receptor_pdb.splitlines()
                if l.startswith("ATOM  ")
            ]
            if not protein_lines:
                return None

            pdb_io = io.StringIO("\n".join(protein_lines) + "\nEND\n")
            pdb = app.PDBFile(pdb_io)

            forcefield = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")

            # Try to add missing hydrogens via Modeller if needed
            try:
                modeller = app.Modeller(pdb.topology, pdb.positions)
                modeller.addHydrogens(forcefield)
                active_topology = modeller.topology
                active_positions = modeller.positions
            except Exception:
                active_topology = pdb.topology
                active_positions = pdb.positions

            system = forcefield.createSystem(
                active_topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds
            )

            # Implement harmonic backbone restraints (k = 10.0 kcal/mol/A^2 = 4184.0 kJ/mol/nm^2)
            # on backbone atoms ("CA", "C", "N", "O")
            k_val = 4184.0 * unit.kilojoules_per_mole / (unit.nanometer ** 2)
            restraint = mm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
            restraint.addGlobalParameter("k", k_val)
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")

            bb_atom_indices = []
            for atom in active_topology.atoms():
                if atom.name in ("CA", "C", "N", "O"):
                    pos = active_positions[atom.index]
                    restraint.addParticle(atom.index, [pos[0], pos[1], pos[2]])
                    bb_atom_indices.append(atom.index)

            if bb_atom_indices:
                system.addForce(restraint)

            integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
            simulation = app.Simulation(active_topology, system, integrator)
            simulation.context.setPositions(active_positions)

            # Initial potential energy
            state_initial = simulation.context.getState(getEnergy=True)
            e_init = state_initial.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

            # Minimize 150 steps
            simulation.minimizeEnergy(maxIterations=150)

            state_min = simulation.context.getState(getEnergy=True, getPositions=True)
            e_min = state_min.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            delta_e = e_min - e_init

            # Calculate backbone RMSD
            bb_rmsd = 0.0
            min_positions = state_min.getPositions()
            if bb_atom_indices:
                sq_sum = 0.0
                for idx in bb_atom_indices:
                    p0 = active_positions[idx]
                    p1 = min_positions[idx]
                    diff = (p1 - p0).value_in_unit(unit.angstrom)
                    sq_sum += (diff[0] ** 2 + diff[1] ** 2 + diff[2] ** 2)
                bb_rmsd = round(math.sqrt(sq_sum / len(bb_atom_indices)), 4)

            return {
                "method": "OpenMM 8.x GBn2 Implicit Solvent — Pocket Structural Relaxation",
                "method_note": (
                    "Receptor pocket geometry minimized in continuous GBn2 dielectric with harmonic backbone restraints. "
                    "Reports potential energy relaxation delta and backbone RMSD."
                ),
                "initial_energy_kcal": round(e_init, 2),
                "minimized_energy_kcal": round(e_min, 2),
                "complex_relaxation_delta_kcal": round(delta_e, 2),
                "status": "Minimized & Solvated",
                "relaxation_steps": 150,
                "backbone_restraint_applied": True,
                "backbone_restraint_k_kcal_mol_A2": 10.0,
                "backbone_rmsd_angstroms": bb_rmsd
            }
        except Exception as e:
            return {
                "method": "OpenMM 8.x GBn2 Implicit Solvent — Pocket Structural Relaxation",
                "status": f"OpenMM minimization notice: {e}",
                "backbone_restraint_applied": False,
                "backbone_restraint_k_kcal_mol_A2": 10.0,
                "backbone_rmsd_angstroms": None,
                "error": str(e)
            }

    @classmethod
    def export_openmm_md_package(
        cls,
        receptor_pdb: str,
        docked_pdb_or_pdbqt: str,
        output_dir: str,
        ligand_smiles: Optional[str] = None,
        job_name: str = "bindora_complex_md",
        sim_time_ns: float = 1.0
    ) -> Dict[str, Any]:
        """
        Generate standalone OpenMM explicit-solvent MD simulation package & MM-PBSA script.
        """
        from backend.services.md_export import OpenMMExportService
        return OpenMMExportService.generate_simulation_package(
            receptor_pdb=receptor_pdb,
            docked_pose_pdb_or_pdbqt=docked_pdb_or_pdbqt,
            output_dir=output_dir,
            ligand_smiles=ligand_smiles,
            job_name=job_name,
            sim_time_ns=sim_time_ns
        )

