"""
Macrocyclic Conformational Sampling Module for Bindora Dock.

Addresses the conformational sampling bottleneck for macrocycles (ring size >= 12, cyclic peptides,
macrolides, kinase macrocycles like Lorlatinib).

Provides:
- Macrocycle ring size and topology detection
- Distance Geometry with srETKDGv3 / ETKDGv3 macrocycle parameters
- MMFF94 / UFF force field relaxation of ring closures
- Relative energy window filtering & RMSD-based conformational pruning
- Multi-conformer ensemble preparation for flexible docking
"""

import math
from typing import Dict, Any, List, Optional, Tuple, Union
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdMolAlign


class MacrocycleConformerEngine:
    """
    Specialized conformational sampling engine for macrocyclic ligands (ring size >= 12).
    Standard distance geometry frequently fails to sample native ring topologies; this engine
    uses experimental torsion distance geometry (srETKDGv3) with random coordinates and force-field
    gradient minimization to identify low-strain ring conformations.
    """

    DEFAULT_MIN_RING_SIZE = 12

    @classmethod
    def is_macrocycle(
        cls,
        mol_or_smiles: Union[Chem.Mol, str],
        min_ring_size: int = DEFAULT_MIN_RING_SIZE
    ) -> Tuple[bool, List[int], List[List[int]]]:
        """
        Check if a molecule contains any simple macrocyclic ring of size >= min_ring_size.
        Uses symmetrized SSSR (Chem.GetSymmSSSR) with rigorous 2-connected cycle graph verification,
        plus topological perimeter combination for fused macrocyclic lactones (e.g. Tacrolimus, Rapamycin).
        Strictly excludes disconnected ring unions and multi-component pseudo-cycles.

        Returns:
            (is_macro, ring_sizes, ring_atom_indices)
        """
        mol = None
        if isinstance(mol_or_smiles, str):
            mol = Chem.MolFromSmiles(mol_or_smiles)
        else:
            mol = mol_or_smiles

        if mol is None:
            return False, [], []

        # 1. Symmetrized SSSR (Chem.GetSymmSSSR)
        symm_rings = [list(r) for r in Chem.GetSymmSSSR(mol)]
        
        # Edge sets for each basis ring
        ring_edges = []
        for r in symm_rings:
            atoms = set(r)
            bonds = {
                b.GetIdx() for b in mol.GetBonds()
                if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
            }
            ring_edges.append(bonds)

        # 2. Candidate edge sets: SymmSSSR rings + fused ring symmetric difference
        candidate_edge_sets = list(ring_edges)
        for i in range(len(ring_edges)):
            for j in range(i + 1, len(ring_edges)):
                shared = ring_edges[i] & ring_edges[j]
                # Rings must share at least one bond (fused/bridged)
                if shared:
                    candidate_edge_sets.append(ring_edges[i] ^ ring_edges[j])

        valid_cycles = []
        seen_atom_tuples = set()

        for edges in candidate_edge_sets:
            adj = {}
            for b_idx in edges:
                b = mol.GetBondWithIdx(b_idx)
                u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                adj.setdefault(u, []).append(v)
                adj.setdefault(v, []).append(u)

            atoms = set(adj.keys())
            if len(atoms) < min_ring_size:
                continue
            if len(edges) != len(atoms):
                continue

            # Rigorous simple cycle check: the induced subgraph in the full molecule
            # must contain exactly len(atoms) bonds (no cross-ring internal chords).
            induced_bonds = [
                b.GetIdx() for b in mol.GetBonds()
                if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
            ]
            if len(induced_bonds) != len(atoms):
                continue

            # Every vertex in a simple cycle must have degree exactly 2
            if not all(len(neighbors) == 2 for neighbors in adj.values()):
                continue

            # Verify single 2-connected simple cycle traversal (no disconnected components)
            start = next(iter(atoms))
            visited = set()
            curr = start
            prev = None
            while curr not in visited:
                visited.add(curr)
                nbrs = adj[curr]
                next_node = nbrs[0] if nbrs[0] != prev else nbrs[1]
                prev = curr
                curr = next_node

            if visited == atoms and curr == start:
                atom_tuple = tuple(sorted(atoms))
                if atom_tuple not in seen_atom_tuples:
                    seen_atom_tuples.add(atom_tuple)
                    valid_cycles.append(list(atoms))

        # Sort cycles descending by size
        valid_cycles.sort(key=lambda c: len(c), reverse=True)
        macro_sizes = [len(c) for c in valid_cycles]
        # Unique sizes preserving descending order
        unique_sizes = []
        for s in macro_sizes:
            if s not in unique_sizes:
                unique_sizes.append(s)

        return len(unique_sizes) > 0, unique_sizes, valid_cycles

    @classmethod
    def sample_macrocycle_conformers(
        cls,
        mol_or_smiles: Union[Chem.Mol, str],
        num_confs: int = 25,
        energy_window: float = 15.0,
        rmsd_threshold: float = 0.5,
        random_seed: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        Sample low-strain conformational ensemble for macrocyclic molecule.

        Parameters:
            mol_or_smiles: RDKit Mol or SMILES string
            num_confs: Number of initial conformers to sample via Distance Geometry
            energy_window: Energy cutoff (kcal/mol) relative to global minimum
            rmsd_threshold: Minimum RMSD (Angstroms) to consider conformers distinct
            random_seed: Reproducible random seed (if None, dynamic random seed is generated)

        Returns:
            Dictionary containing best conformer mol, energies, and sampling metadata.
        """
        if isinstance(mol_or_smiles, str):
            mol = Chem.MolFromSmiles(mol_or_smiles)
        else:
            mol = Chem.Mol(mol_or_smiles)

        if mol is None:
            raise ValueError("Invalid molecule provided for macrocycle sampling.")

        # Ensure dynamic random seed when None is provided
        if random_seed is None:
            import secrets
            random_seed = int(secrets.randbelow(2147483647) + 1)
        else:
            random_seed = int(random_seed)

        # Scale conformer count for large/complex macrocycles (> 50 heavy atoms, e.g. CsA)
        # to prevent process lockups while guaranteeing conformational sampling
        heavy_count = mol.GetNumHeavyAtoms()
        if heavy_count > 50 and num_confs > 3:
            num_confs = 3

        is_macro, ring_sizes, ring_rings = cls.is_macrocycle(mol)
        mol_h = Chem.AddHs(mol)

        # 1. Setup Distance Geometry Parameters (Tier 1: srETKDGv3 macrocycle parameters)
        params = None
        engine_name = "srETKDGv3"
        try:
            if hasattr(AllChem, "srETKDGv3"):
                params = AllChem.srETKDGv3()
            else:
                params = AllChem.ETKDGv3()
                engine_name = "ETKDGv3"
        except Exception:
            params = AllChem.ETKDGv3()
            engine_name = "ETKDGv3"

        params.randomSeed = random_seed
        params.useRandomCoords = True
        params.boxSizeMult = 2.0
        if hasattr(params, "useMacrocycleTorsions"):
            params.useMacrocycleTorsions = True
        if hasattr(params, "useMacrocycle14config"):
            params.useMacrocycle14config = True
        if hasattr(params, "boundsMatForceScaling"):
            params.boundsMatForceScaling = 1.0
        params.clearConfs = True
        params.numThreads = 0

        # 2. Embed multiple conformers with multi-tier fallback
        conf_ids = []
        try:
            conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params))
        except Exception:
            conf_ids = []

        # Tier 2 fallback: standard ETKDGv3 without specialized macrocycle torsion constraints
        if len(conf_ids) == 0:
            engine_name = "ETKDGv3_Fallback"
            try:
                fb_v3 = AllChem.ETKDGv3()
                fb_v3.randomSeed = random_seed
                fb_v3.useRandomCoords = True
                fb_v3.numThreads = 0
                conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=fb_v3))
            except Exception:
                conf_ids = []

        # Tier 3 fallback: basic ETKDG with random coordinates
        if len(conf_ids) == 0:
            engine_name = "ETKDG_Random_Fallback"
            try:
                fb_params = AllChem.ETKDG()
                fb_params.randomSeed = random_seed
                fb_params.useRandomCoords = True
                fb_params.useExpTorsionAnglePrefs = False
                fb_params.useBasicKnowledge = False
                fb_params.numThreads = 0
                conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=fb_params))
            except Exception:
                conf_ids = []

        # Tier 4 fallback: single standard embedding
        if len(conf_ids) == 0:
            try:
                res_embed = AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=random_seed)
                conf_ids = [0] if res_embed == 0 and mol_h.GetNumConformers() > 0 else []
            except Exception:
                conf_ids = []

        if len(conf_ids) == 0:
            raise RuntimeError("Distance geometry failed to generate 3D coordinates for macrocycle.")

        # 3. Energy minimization of all conformers using MMFF94 (fallback to UFF)
        ff_type = "MMFF94"
        raw_energies: List[Tuple[int, float]] = []

        # Check if MMFF parameters exist
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol_h)
        if mmff_props is not None:
            results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=300, numThreads=0)
            for cid, (converged, energy) in zip(conf_ids, results):
                raw_energies.append((cid, energy))
        else:
            ff_type = "UFF"
            results = AllChem.UFFOptimizeMoleculeConfs(mol_h, maxIters=300, numThreads=0)
            for cid, (converged, energy) in zip(conf_ids, results):
                raw_energies.append((cid, energy))

        if not raw_energies:
            raise RuntimeError("Force field optimization failed on macrocycle conformers.")

        # 4. Sort by absolute energy
        raw_energies.sort(key=lambda x: x[1])
        min_energy = raw_energies[0][1]

        # 5. Energy window filtering
        within_window: List[Tuple[int, float, float]] = []
        for cid, energy in raw_energies:
            rel_e = energy - min_energy
            if rel_e <= energy_window:
                within_window.append((cid, energy, rel_e))

        # 6. RMSD Pruning to eliminate duplicate conformers
        retained_confs: List[Tuple[int, float, float]] = []
        for item in within_window:
            cid, energy, rel_e = item
            if not retained_confs:
                retained_confs.append(item)
                continue

            # Check RMSD against all already retained conformers using fast coordinate alignment
            is_duplicate = False
            for ret_cid, _, _ in retained_confs:
                try:
                    rmsd = rdMolAlign.AlignMol(mol_h, mol_h, prbCid=cid, refCid=ret_cid)
                    if rmsd < rmsd_threshold:
                        is_duplicate = True
                        break
                except Exception:
                    pass

            if not is_duplicate:
                retained_confs.append(item)

        # 7. Construct best conformer mol (with lowest energy coordinates)
        best_cid = retained_confs[0][0]
        best_mol = Chem.Mol(mol_h)
        # Keep only best conformer
        for c in list(best_mol.GetConformers()):
            if c.GetId() != best_cid:
                best_mol.RemoveConformer(c.GetId())

        # Construct list of separate Mol objects for retained ensemble
        ensemble_mols = []
        for r_cid, r_e, r_rel in retained_confs:
            conf_mol = Chem.Mol(mol_h)
            for c in list(conf_mol.GetConformers()):
                if c.GetId() != r_cid:
                    conf_mol.RemoveConformer(c.GetId())
            conf_mol.SetProp("_Energy_kcal", f"{r_e:.2f}")
            conf_mol.SetProp("_RelEnergy_kcal", f"{r_rel:.2f}")
            ensemble_mols.append(conf_mol)

        return {
            "is_macrocycle": is_macro,
            "max_ring_size": max(ring_sizes) if ring_sizes else 0,
            "macrocycle_ring_sizes": ring_sizes,
            "sampling_engine": engine_name,
            "force_field": ff_type,
            "initial_conformers_sampled": len(conf_ids),
            "conformers_within_energy_window": len(within_window),
            "conformers_retained_after_rmsd_pruning": len(retained_confs),
            "global_min_energy_kcal": round(min_energy, 2),
            "energy_window_cutoff_kcal": energy_window,
            "rmsd_pruning_cutoff_angstroms": rmsd_threshold,
            "relative_energies_kcal": [round(x[2], 2) for x in retained_confs],
            "best_mol": best_mol,
            "ensemble_mols": ensemble_mols,
            "seed_used": random_seed
        }

    @classmethod
    def prepare_macrocycle_ensemble_for_docking(
        cls,
        mol: Chem.Mol,
        max_docking_confs: int = 5
    ) -> List[Chem.Mol]:
        """
        Produce top N diverse, low-energy macrocycle conformers formatted for docking.
        """
        sampled = cls.sample_macrocycle_conformers(
            mol,
            num_confs=30,
            energy_window=12.0,
            rmsd_threshold=0.75
        )
        ensemble = sampled.get("ensemble_mols", [])
        return ensemble[:max_docking_confs]
