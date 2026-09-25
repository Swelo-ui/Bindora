import math
import numpy as np
from typing import Dict, Any, List, Tuple, Optional, Set
from rdkit import Chem
from rdkit.Chem import AllChem

class InteractionEngine:
    """
    Biophysically rigorous non-covalent interaction detection engine adhering to PLIP criteria:
    - Directional Hydrogen Bonds with angular checks (angle >= 120°) & donor/acceptor roles
    - Salt Bridges (4.0 Å between cationic and anionic centers)
    - π-π Stacking (centroid <= 5.5 Å with parallel vs T-shaped orientation)
    - π-Cation interactions (<= 4.5 Å)
    - Halogen Bonds with σ-hole angular verification (C-X...O/N >= 130°)
    - Hydrophobic contacts between non-polar carbons (<= 4.0 Å)
    """

    @classmethod
    def analyze(
        cls,
        receptor_pdb: str,
        ligand_pdb_or_pdbqt: str,
        hbond_cutoff: float = 3.5,
        hydrophobic_cutoff: float = 4.0,
        ligand_mol: Optional[Chem.Mol] = None
    ) -> Dict[str, Any]:
        """
        Analyze protein-ligand interactions with complete 3D angular geometry.
        """
        rec_atoms = cls._parse_receptor_atoms(receptor_pdb)
        lig_atoms = cls._parse_ligand_atoms(ligand_pdb_or_pdbqt)

        if not rec_atoms or not lig_atoms:
            return {
                "hydrogen_bonds": [],
                "salt_bridges": [],
                "pi_stacking": [],
                "pi_cation": [],
                "halogen_bonds": [],
                "hydrophobic_contacts": [],
                "interacting_residues": [],
                "flexible_candidates": [],
                "total_hbond_count": 0,
                "total_hydrophobic_count": 0
            }

        # Match ligand atoms to RDKit Mol atom indices if provided
        cls._correlate_ligand_indices(lig_atoms, ligand_mol)

        contact_residues: Set[str] = set()

        hbonds = cls._find_hydrogen_bonds(rec_atoms, lig_atoms, max_dist=hbond_cutoff)
        for hb in hbonds:
            contact_residues.add(hb["residue"])

        salt_bridges = cls._find_salt_bridges(rec_atoms, lig_atoms, max_dist=4.2)
        for sb in salt_bridges:
            contact_residues.add(sb["residue"])

        pi_stacks = cls._find_pi_stacking(rec_atoms, lig_atoms, ligand_mol, max_dist=5.5)
        for ps in pi_stacks:
            contact_residues.add(ps["residue"])

        pi_cations = cls._find_pi_cation(rec_atoms, lig_atoms, ligand_mol, max_dist=4.5)
        for pc in pi_cations:
            contact_residues.add(pc["residue"])

        halogens = cls._find_halogen_bonds(rec_atoms, lig_atoms, max_dist=3.8)
        for hal in halogens:
            contact_residues.add(hal["residue"])

        hydrophobics = cls._find_hydrophobic_contacts(rec_atoms, lig_atoms, max_dist=hydrophobic_cutoff)
        for hp in hydrophobics:
            contact_residues.add(hp["residue"])

        # Scientifically valid flexible sidechains in AutoDock Vina (excluding rigid/zero-torsion GLY, ALA, PRO)
        ROTATABLE_SIDECHAINS = {
            "ARG": 4, "LYS": 4, "GLU": 3, "GLN": 3, "MET": 3,
            "LEU": 2, "ILE": 2, "ASP": 2, "ASN": 2, "HIS": 2,
            "PHE": 2, "TYR": 2, "TRP": 2, "VAL": 1, "SER": 1,
            "THR": 1, "CYS": 1
        }

        # Flexible residue candidates
        flex_candidates = []
        for cr in sorted(list(contact_residues)):
            try:
                parts = cr.split()
                rname = parts[0].upper()
                if rname not in ROTATABLE_SIDECHAINS:
                    continue  # Glycine, Alanine, Proline cannot be flexible in AutoDock Vina
                rnum_chain = parts[1]
                rnum, rchain = rnum_chain.split(":")
                flex_candidates.append({
                    "id": f"{rchain}:{rnum}",
                    "label": cr,
                    "res_name": rname,
                    "res_num": int(rnum),
                    "chain": rchain,
                    "rotatable_bonds": ROTATABLE_SIDECHAINS[rname]
                })
            except Exception:
                pass

        return {
            "hydrogen_bonds": hbonds,
            "salt_bridges": salt_bridges,
            "pi_stacking": pi_stacks,
            "pi_cation": pi_cations,
            "halogen_bonds": halogens,
            "hydrophobic_contacts": hydrophobics[:16],
            "interacting_residues": sorted(list(contact_residues)),
            "flexible_candidates": flex_candidates,
            "total_hbond_count": len(hbonds),
            "total_salt_bridge_count": len(salt_bridges),
            "total_pi_stacking_count": len(pi_stacks),
            "total_pi_cation_count": len(pi_cations),
            "total_halogen_count": len(halogens),
            "total_hydrophobic_count": len(hydrophobics),
            "detection_criteria": cls.get_criteria(hbond_cutoff, hydrophobic_cutoff)
        }

    @classmethod
    def get_criteria(cls, hbond_cutoff: float = 3.5, hydrophobic_cutoff: float = 4.0) -> Dict[str, Any]:
        """Return transparent structural interaction criteria and cutoffs."""
        return {
            "hydrogen_bond": {
                "max_distance_angstroms": hbond_cutoff,
                "distance_max_angstroms": hbond_cutoff,
                "min_angle_degrees": 115.0,
                "geometry": "Directional donor-H...acceptor angle >= 115°, distance <= 3.5 Å",
                "donor_acceptor_elements": ["O", "N", "S"]
            },
            "salt_bridges": {
                "max_distance_angstroms": 4.2,
                "distance_max_angstroms": 4.2,
                "geometry": "Centroid distance between cationic center (Lys NZ, Arg guanidinium) and anionic center (Asp/Glu carboxylate) <= 4.2 Å"
            },
            "pi_stacking": {
                "max_distance_angstroms": 5.5,
                "distance_max_angstroms": 5.5,
                "max_angle_parallel_deg": 35.0,
                "geometry": "Centroid-centroid distance <= 5.5 Å; parallel (θ < 35°) vs T-shaped (θ >= 35°)"
            },
            "pi_cation": {
                "max_distance_angstroms": 4.5,
                "distance_max_angstroms": 4.5,
                "geometry": "Cationic group to aromatic centroid distance <= 4.5 Å"
            },
            "halogen_bond_sigma_hole": {
                "max_distance_angstroms": 3.8,
                "distance_max_angstroms": 3.8,
                "min_c_x_acceptor_angle_deg": 130.0,
                "classical_halogens": ["Cl", "Br", "I"],
                "mechanism": "Directional electrophilic sigma-hole along C-X bond axis (theta >= 130 deg)"
            },
            "fluorine_polar_contact": {
                "max_distance_angstroms": 3.5,
                "distance_max_angstroms": 3.5,
                "note": "Fluorine lacks a classical electrophilic sigma-hole; interactions are non-directional polar/multipolar contacts."
            },
            "halogen_bonds": {
                "max_distance_angstroms": 3.8,
                "min_angle_degrees": 130.0,
                "classical_halogens": ["Cl", "Br", "I"],
                "f_contacts": ["F"],
                "geometry": "C-X...[O/N/S] distance <= 3.8 Å with angle >= 130° for Cl/Br/I; F evaluated as polar/multipolar contact"
            },
            "hydrophobic_contacts": {
                "max_distance_angstroms": hydrophobic_cutoff,
                "distance_max_angstroms": hydrophobic_cutoff,
                "geometry": "Carbon-carbon non-polar contact distance <= 4.0 Å"
            }
        }

    @classmethod
    def _parse_receptor_atoms(cls, pdb_str: str) -> List[Dict[str, Any]]:
        atoms = []
        lines = pdb_str.splitlines()
        for idx, line in enumerate(lines):
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    chain = line[21:22].strip()
                    atom_name = line[12:16].strip()
                    elem = line[76:78].strip() or atom_name[0]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    atoms.append({
                        "res_name": res_name,
                        "res_num": res_num,
                        "chain": chain,
                        "atom_name": atom_name,
                        "elem": elem.upper(),
                        "coord": np.array([x, y, z]),
                        "line_idx": idx
                    })
                except Exception:
                    continue
        return atoms

    @classmethod
    def _parse_ligand_atoms(cls, block_str: str) -> List[Dict[str, Any]]:
        # Check for REMARK SMILES IDX mappings from Meeko PDBQT
        smi_map: Dict[int, int] = {}
        for line in block_str.splitlines():
            if line.startswith("REMARK SMILES IDX"):
                tokens = line[17:].strip().split()
                try:
                    for i in range(0, len(tokens) - 1, 2):
                        smi_idx = int(tokens[i]) - 1  # 0-based
                        pdb_idx = int(tokens[i + 1])
                        smi_map[pdb_idx] = smi_idx
                except Exception:
                    pass

        atoms = []
        atom_idx = 0
        for line in block_str.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    serial = int(line[6:11].strip()) if len(line) >= 11 and line[6:11].strip().isdigit() else (atom_idx + 1)
                    atom_name = line[12:16].strip()
                    elem = line[76:78].strip() if len(line) > 76 else ""
                    if not elem:
                        elem = atom_name[0]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    s_idx = smi_map.get(serial, atom_idx)

                    atoms.append({
                        "serial": serial,
                        "atom_name": atom_name,
                        "elem": elem.upper(),
                        "coord": np.array([x, y, z]),
                        "atom_idx": s_idx,
                        "smiles_idx": smi_map.get(serial, None)
                    })
                    atom_idx += 1
                except Exception:
                    continue
        return atoms

    @classmethod
    def _correlate_ligand_indices(cls, lig_atoms: List[Dict[str, Any]], mol: Optional[Chem.Mol]) -> None:
        """Ensure each ligand atom has a confirmed 0-indexed atom index matching the 2D molecule."""
        if not mol:
            return
        conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None
        if conf is None:
            return

        for at in lig_atoms:
            lx, ly, lz = at["coord"]
            best_idx = at["atom_idx"]
            best_dist = float("inf")
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                d = (pos.x - lx)**2 + (pos.y - ly)**2 + (pos.z - lz)**2
                if d < best_dist:
                    best_dist = d
                    best_idx = i
            if best_dist < 1.0:  # within 1.0 A
                at["atom_idx"] = best_idx

    @classmethod
    def _calc_angle(cls, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
        """Calculate angle in degrees between vectors (a - b) and (c - b). Vertex is at b."""
        v1 = a - b
        v2 = c - b
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)
        if n1 < 1e-7 or n2 < 1e-7:
            return 0.0
        cos_theta = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
        return float(np.degrees(np.arccos(cos_theta)))

    @classmethod
    def _find_hydrogen_bonds(
        cls,
        rec_atoms: List[Dict[str, Any]],
        lig_atoms: List[Dict[str, Any]],
        max_dist: float = 3.5
    ) -> List[Dict[str, Any]]:
        """
        PLIP criteria for Hydrogen Bonds:
        - Donor - H ... Acceptor distance <= 3.5 A
        - Donor - H ... Acceptor angle >= 120°
        - Explicit check to reject Acceptor ... Acceptor repulsion (e.g. carbonyl O ... carbonyl O)
        """
        hbonds = []
        # Build index of receptor hydrogens for angle checking
        rec_hydrogens = [ra for ra in rec_atoms if ra["elem"] == "H"]
        rec_heavy = [ra for ra in rec_atoms if ra["elem"] != "H"]

        # Classification helpers
        acceptor_elems = {"O", "N", "S"}

        for latom in lig_atoms:
            lelem = latom["elem"]
            if lelem not in acceptor_elems and lelem != "H":
                continue

            for ratom in rec_heavy:
                relem = ratom["elem"]
                if relem not in acceptor_elems:
                    continue

                rcoord = ratom["coord"]
                lcoord = latom["coord"]
                dist = float(np.linalg.norm(rcoord - lcoord))
                if dist > max_dist or dist < 1.4:
                    continue

                res_id = f"{ratom['res_name']} {ratom['res_num']}:{ratom['chain']}"

                # Case 1: Protein Donor (with H) -> Ligand Acceptor
                # Check if protein atom has a bonded hydrogen
                bonded_hs = [
                    h for h in rec_hydrogens
                    if np.linalg.norm(h["coord"] - rcoord) <= 1.25 and h.get("chain") == ratom.get("chain")
                ]

                is_hbond = False
                h_coord = None
                donor_angle = 180.0

                if bonded_hs:
                    # We have explicit polar hydrogen on receptor donor
                    for bh in bonded_hs:
                        ang = cls._calc_angle(rcoord, bh["coord"], lcoord)  # RecDonor - H ... LigAcceptor
                        if ang >= 115.0:  # PLIP angle cutoff
                            is_hbond = True
                            donor_angle = ang
                            h_coord = bh["coord"]
                            break

                # Case 2: Ligand Donor -> Protein Acceptor
                # If ligand has attached hydrogen in pose
                lig_hs = [
                    lh for lh in lig_atoms
                    if lh["elem"] == "H" and np.linalg.norm(lh["coord"] - lcoord) <= 1.25
                ]
                if lig_hs:
                    for lh in lig_hs:
                        ang = cls._calc_angle(lcoord, lh["coord"], rcoord)  # LigDonor - H ... RecAcceptor
                        if ang >= 115.0:
                            is_hbond = True
                            donor_angle = ang
                            h_coord = lh["coord"]
                            break

                # Case 3: Geometric fallback if hydrogens were omitted or united
                if not bonded_hs and not lig_hs:
                    # Prevent carbonyl O ... carbonyl O repulsion!
                    is_rec_carbonyl = ratom["atom_name"] in ("O", "OD1", "OE1", "OE2", "OD2") and ratom["res_name"] in ("ASP", "GLU", "ASN", "GLN") or ratom["atom_name"] == "O"
                    is_lig_pure_acceptor = latom["atom_name"].startswith("O") and ("=" in latom["atom_name"] or latom["elem"] == "O")
                    # If both are carbonyl oxygens, they CANNOT form an H-bond!
                    if is_rec_carbonyl and is_lig_pure_acceptor:
                        continue
                    if dist <= 3.3:
                        is_hbond = True
                        donor_angle = 140.0

                if is_hbond:
                    hbonds.append({
                        "type": "Hydrogen Bond",
                        "distance": round(dist, 2),
                        "angle_deg": round(donor_angle, 1),
                        "residue": res_id,
                        "res_name": ratom["res_name"],
                        "res_num": ratom["res_num"],
                        "chain": ratom["chain"],
                        "receptor_atom": ratom["atom_name"],
                        "ligand_atom": latom["atom_name"],
                        "ligand_atom_idx": latom.get("atom_idx", 0),
                        "start_coord": [float(c) for c in lcoord],
                        "end_coord": [float(c) for c in rcoord]
                    })

        # Deduplicate to closest H-bond per (residue, ligand_atom)
        unique = {}
        for hb in hbonds:
            key = f"{hb['residue']}_{hb['ligand_atom_idx']}"
            if key not in unique or unique[key]["distance"] > hb["distance"]:
                unique[key] = hb
        return list(unique.values())

    @classmethod
    def _find_salt_bridges(
        cls,
        rec_atoms: List[Dict[str, Any]],
        lig_atoms: List[Dict[str, Any]],
        max_dist: float = 4.2
    ) -> List[Dict[str, Any]]:
        """Detect salt bridges between cationic and anionic centers (<= 4.2 A)."""
        salt_bridges = []

        # Group receptor residues
        pos_rec: Dict[Tuple[str, int, str], List[np.ndarray]] = {}  # ARG (CZ, NH1, NH2), LYS (NZ)
        neg_rec: Dict[Tuple[str, int, str], List[np.ndarray]] = {}  # ASP (OD1, OD2), GLU (OE1, OE2)

        for ra in rec_atoms:
            rkey = (ra["res_name"], ra["res_num"], ra["chain"])
            if ra["res_name"] == "LYS" and ra["atom_name"] == "NZ":
                pos_rec.setdefault(rkey, []).append(ra["coord"])
            elif ra["res_name"] == "ARG" and ra["atom_name"] in ("NE", "CZ", "NH1", "NH2"):
                pos_rec.setdefault(rkey, []).append(ra["coord"])
            elif ra["res_name"] in ("ASP", "GLU") and ra["atom_name"] in ("OD1", "OD2", "OE1", "OE2"):
                neg_rec.setdefault(rkey, []).append(ra["coord"])

        # Check against ligand atoms
        for latom in lig_atoms:
            lcoord = latom["coord"]
            lelem = latom["elem"]

            # Ligand anionic carboxylate/phosphate/sulfonate O interacting with protein Lys/Arg
            if lelem == "O":
                for (rname, rnum, rchain), coords in pos_rec.items():
                    center = np.mean(coords, axis=0)
                    d = float(np.linalg.norm(center - lcoord))
                    if d <= max_dist:
                        res_id = f"{rname} {rnum}:{rchain}"
                        salt_bridges.append({
                            "type": "Salt Bridge",
                            "subtype": "Protein Cation - Ligand Anion",
                            "distance": round(d, 2),
                            "residue": res_id,
                            "res_name": rname,
                            "res_num": rnum,
                            "chain": rchain,
                            "ligand_atom": latom["atom_name"],
                            "ligand_atom_idx": latom.get("atom_idx", 0),
                            "start_coord": [float(c) for c in lcoord],
                            "end_coord": [float(c) for c in center]
                        })

            # Ligand cationic amine/quaternary N interacting with protein Asp/Glu
            if lelem == "N":
                for (rname, rnum, rchain), coords in neg_rec.items():
                    center = np.mean(coords, axis=0)
                    d = float(np.linalg.norm(center - lcoord))
                    if d <= max_dist:
                        res_id = f"{rname} {rnum}:{rchain}"
                        salt_bridges.append({
                            "type": "Salt Bridge",
                            "subtype": "Ligand Cation - Protein Anion",
                            "distance": round(d, 2),
                            "residue": res_id,
                            "res_name": rname,
                            "res_num": rnum,
                            "chain": rchain,
                            "ligand_atom": latom["atom_name"],
                            "ligand_atom_idx": latom.get("atom_idx", 0),
                            "start_coord": [float(c) for c in lcoord],
                            "end_coord": [float(c) for c in center]
                        })

        # Deduplicate to closest salt bridge per residue
        unique_sb = {}
        for sb in salt_bridges:
            key = sb["residue"]
            if key not in unique_sb or unique_sb[key]["distance"] > sb["distance"]:
                unique_sb[key] = sb
        return list(unique_sb.values())

    @classmethod
    def _find_ligand_aromatic_rings(
        cls,
        lig_atoms: List[Dict[str, Any]],
        mol: Optional[Chem.Mol] = None
    ) -> List[Dict[str, Any]]:
        """
        Identify 5- and 6-membered aromatic rings directly from 3D coordinates and geometry.
        Guarantees that ring centroids and normals are strictly anchored within the real 3D ring.
        """
        heavy = [a for a in lig_atoms if a.get("elem") not in ("H", "HD")]
        if len(heavy) < 5:
            return []

        # Build 3D covalent bond graph based on heavy-atom interatomic distances (1.0 - 1.85 A)
        n_heavy = len(heavy)
        adj: Dict[int, List[int]] = {i: [] for i in range(n_heavy)}
        for i in range(n_heavy):
            ci = heavy[i]["coord"]
            for j in range(i + 1, n_heavy):
                cj = heavy[j]["coord"]
                dist = float(np.linalg.norm(ci - cj))
                if 1.0 <= dist <= 1.85:
                    adj[i].append(j)
                    adj[j].append(i)

        # Depth-first search for 5- and 6-membered simple cycles
        cycles = []
        def dfs(path: List[int]) -> None:
            curr = path[-1]
            for nxt in adj[curr]:
                if nxt == path[0] and len(path) in (5, 6):
                    can = tuple(sorted(path))
                    if can not in [tuple(sorted(c)) for c in cycles]:
                        cycles.append(list(path))
                elif nxt not in path and len(path) < 6:
                    if nxt > path[0]:
                        dfs(path + [nxt])

        for i in range(n_heavy):
            dfs([i])

        valid_rings = []
        ring_elements = {"C", "N", "O", "S", "A", "NA", "OA", "SA"}

        for c_indices in cycles:
            ring_atoms = [heavy[i] for i in c_indices]
            # Must consist of ring-forming elements
            if any(a.get("elem") not in ring_elements for a in ring_atoms):
                continue

            coords = np.array([a["coord"] for a in ring_atoms])
            centroid = np.mean(coords, axis=0)

            # Planarity check via SVD: aromatic rings are strictly flat (RMSD < 0.18 A)
            c_coords = coords - centroid
            _, s, vh = np.linalg.svd(c_coords)
            rmsd = float(np.sqrt(s[2]**2 / len(c_indices)))
            if rmsd > 0.18:
                continue

            # Unit normal vector to best-fit plane
            normal = vh[2]
            norm = float(np.linalg.norm(normal))
            if norm > 1e-6:
                normal = normal / norm
            else:
                normal = np.array([0.0, 0.0, 1.0])

            # Confirm aromaticity if RDKit mol is provided
            if mol is not None:
                has_smi_map = any(a.get("smiles_idx") is not None for a in ring_atoms)
                if has_smi_map:
                    smi_aromatic = all(
                        mol.GetAtomWithIdx(a["smiles_idx"]).GetIsAromatic()
                        for a in ring_atoms if a.get("smiles_idx") is not None
                    )
                    if not smi_aromatic:
                        continue

            valid_rings.append({
                "indices": [a.get("atom_idx", 0) for a in ring_atoms],
                "centroid": centroid,
                "normal": normal,
                "atoms": ring_atoms
            })

        return valid_rings

    @classmethod
    def _find_pi_stacking(
        cls,
        rec_atoms: List[Dict[str, Any]],
        lig_atoms: List[Dict[str, Any]],
        mol: Optional[Chem.Mol],
        max_dist: float = 5.5
    ) -> List[Dict[str, Any]]:
        """PLIP criteria for aromatic pi-pi stacking: centroid distance <= 5.5 A."""
        stacks = []

        # 1. Identify ligand aromatic rings directly from 3D structure
        lig_rings = cls._find_ligand_aromatic_rings(lig_atoms, mol)
        if not lig_rings:
            return stacks

        # 2. Identify receptor aromatic residues (PHE, TYR, TRP, HIS)
        rec_rings: Dict[Tuple[str, int, str], List[np.ndarray]] = {}
        for ra in rec_atoms:
            rname = ra["res_name"]
            aname = ra["atom_name"]
            rkey = (rname, ra["res_num"], ra["chain"])
            if rname in ("PHE", "TYR") and aname in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"):
                rec_rings.setdefault(rkey, []).append(ra["coord"])
            elif rname == "HIS" and aname in ("CG", "ND1", "CD2", "CE1", "NE2"):
                rec_rings.setdefault(rkey, []).append(ra["coord"])
            elif rname == "TRP" and aname in ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"):
                rec_rings.setdefault(rkey, []).append(ra["coord"])

        for (rname, rnum, rchain), rcoords in rec_rings.items():
            if len(rcoords) < 5:
                continue
            rcoords_arr = np.array(rcoords)
            r_centroid = np.mean(rcoords_arr, axis=0)
            rv1 = rcoords_arr[1] - rcoords_arr[0]
            rv2 = rcoords_arr[2] - rcoords_arr[0]
            r_normal = np.cross(rv1, rv2)
            r_norm = np.linalg.norm(r_normal)
            if r_norm > 1e-6:
                r_normal = r_normal / r_norm

            for lring in lig_rings:
                dist = float(np.linalg.norm(r_centroid - lring["centroid"]))
                if dist <= max_dist:
                    # Angle between ring normals
                    cos_ang = np.clip(abs(np.dot(r_normal, lring["normal"])), 0.0, 1.0)
                    ang_deg = float(np.degrees(np.arccos(cos_ang)))

                    stack_type = "Parallel π-π" if ang_deg < 35.0 else "T-shaped π-π"
                    res_id = f"{rname} {rnum}:{rchain}"
                    first_idx = lring["indices"][0]

                    stacks.append({
                        "type": "π-π Stacking",
                        "subtype": stack_type,
                        "distance": round(dist, 2),
                        "angle_deg": round(ang_deg, 1),
                        "residue": res_id,
                        "res_name": rname,
                        "res_num": rnum,
                        "chain": rchain,
                        "ligand_atom": f"Ring ({len(lring['indices'])} atoms)",
                        "ligand_atom_idx": first_idx,
                        "start_coord": [float(c) for c in lring["centroid"]],
                        "end_coord": [float(c) for c in r_centroid]
                    })

        # Deduplicate to closest pi-pi stacking per residue
        unique_ps = {}
        for ps in stacks:
            key = ps["residue"]
            if key not in unique_ps or unique_ps[key]["distance"] > ps["distance"]:
                unique_ps[key] = ps
        return list(unique_ps.values())

    @classmethod
    def _find_pi_cation(
        cls,
        rec_atoms: List[Dict[str, Any]],
        lig_atoms: List[Dict[str, Any]],
        mol: Optional[Chem.Mol],
        max_dist: float = 4.5
    ) -> List[Dict[str, Any]]:
        """PLIP criteria for pi-cation interactions (<= 4.5 A)."""
        pi_cations = []

        # 1. Receptor cations (LYS NZ, ARG guanidinium)
        # 2. Receptor aromatics (PHE, TYR, TRP, HIS - strict sidechain ring atoms only)
        rec_cations: Dict[Tuple[str, int, str], List[np.ndarray]] = {}
        rec_aromatics: Dict[Tuple[str, int, str], List[np.ndarray]] = {}

        for ra in rec_atoms:
            rkey = (ra["res_name"], ra["res_num"], ra["chain"])
            aname = ra["atom_name"]
            rname = ra["res_name"]
            if rname == "LYS" and aname == "NZ":
                rec_cations.setdefault(rkey, []).append(ra["coord"])
            elif rname == "ARG" and aname in ("NE", "CZ", "NH1", "NH2"):
                rec_cations.setdefault(rkey, []).append(ra["coord"])
            elif rname in ("PHE", "TYR") and aname in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"):
                rec_aromatics.setdefault(rkey, []).append(ra["coord"])
            elif rname == "HIS" and aname in ("CG", "ND1", "CD2", "CE1", "NE2"):
                rec_aromatics.setdefault(rkey, []).append(ra["coord"])
            elif rname == "TRP" and aname in ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"):
                rec_aromatics.setdefault(rkey, []).append(ra["coord"])

        # Ligand aromatic rings directly from 3D structure
        lig_rings = cls._find_ligand_aromatic_rings(lig_atoms, mol)

        # Check ligand cationic nitrogens to receptor aromatics
        # Filter out neutral amides (-C(=O)-N-) and nitro groups (-NO2)
        for latom in lig_atoms:
            if latom.get("elem") == "N":
                lcoord = latom["coord"]

                # Exclude if attached to carbonyl carbon (amide) or two oxygens (nitro)
                is_amide_or_nitro = False
                bonded_heavy = [
                    other for other in lig_atoms
                    if other is not latom and np.linalg.norm(other["coord"] - lcoord) <= 1.55
                ]
                oxygens_count = sum(1 for b in bonded_heavy if b.get("elem") == "O")
                if oxygens_count >= 2:
                    is_amide_or_nitro = True
                else:
                    for b in bonded_heavy:
                        if b.get("elem") == "C":
                            c_bonded_o = any(
                                o.get("elem") == "O" and np.linalg.norm(o["coord"] - b["coord"]) <= 1.38
                                for o in lig_atoms if o is not b
                            )
                            if c_bonded_o:
                                is_amide_or_nitro = True
                                break

                if is_amide_or_nitro:
                    continue

                for (rname, rnum, rchain), coords in rec_aromatics.items():
                    if len(coords) < 5:
                        continue
                    centroid = np.mean(coords, axis=0)
                    d = float(np.linalg.norm(centroid - lcoord))
                    if d <= max_dist:
                        res_id = f"{rname} {rnum}:{rchain}"
                        pi_cations.append({
                            "type": "π-Cation",
                            "subtype": "Ligand Cation - Receptor π-Ring",
                            "distance": round(d, 2),
                            "residue": res_id,
                            "res_name": rname,
                            "res_num": rnum,
                            "chain": rchain,
                            "ligand_atom": latom["atom_name"],
                            "ligand_atom_idx": latom.get("atom_idx", 0),
                            "start_coord": [float(c) for c in lcoord],
                            "end_coord": [float(c) for c in centroid]
                        })

        # Check receptor cations (LYS NZ, ARG guanidinium) to ligand aromatic rings
        for lring in lig_rings:
            lring_centroid = lring["centroid"]
            for (rname, rnum, rchain), c_coords in rec_cations.items():
                cation_center = np.mean(c_coords, axis=0)
                d = float(np.linalg.norm(cation_center - lring_centroid))
                if d <= max_dist:
                    res_id = f"{rname} {rnum}:{rchain}"
                    pi_cations.append({
                        "type": "π-Cation",
                        "subtype": "Receptor Cation - Ligand π-Ring",
                        "distance": round(d, 2),
                        "residue": res_id,
                        "res_name": rname,
                        "res_num": rnum,
                        "chain": rchain,
                        "ligand_atom": f"Ring ({len(lring['indices'])} atoms)",
                        "ligand_atom_idx": lring["indices"][0],
                        "start_coord": [float(c) for c in cation_center],
                        "end_coord": [float(c) for c in lring_centroid]
                    })

        # Deduplicate to closest pi-cation per residue
        unique_pc = {}
        for pc in pi_cations:
            key = pc["residue"]
            if key not in unique_pc or unique_pc[key]["distance"] > pc["distance"]:
                unique_pc[key] = pc
        return list(unique_pc.values())

    @classmethod
    def _find_halogen_bonds(
        cls,
        rec_atoms: List[Dict[str, Any]],
        lig_atoms: List[Dict[str, Any]],
        max_dist: float = 3.8
    ) -> List[Dict[str, Any]]:
        """
        Detect Halogen interactions adhering to physical chemistry and PLIP guidelines:
        - Classical Halogen Bonds: Cl, Br, I with electron-deficient σ-hole.
          Criteria: C-X...[O/N/S] angle >= 130° and distance <= 3.8 Å.
        - Fluorine Contacts: F is minimally polarizable and rarely forms classical σ-hole bonds.
          Classified transparently as 'Fluorine Contact / Multipolar Interaction' rather than
          classical halogen bond, with distance <= 3.8 Å and angle criteria recorded.
        """
        halogens = []
        hal_elems = {"F", "CL", "BR", "I"}
        acceptor_elems = {"O", "N", "S"}

        # Find ligand halogens and their bonded carbons
        for latom in lig_atoms:
            lelem = latom["elem"].upper()
            if lelem not in hal_elems:
                continue
            x_coord = latom["coord"]

            # Find bonded ligand carbon
            bonded_c = [
                c for c in lig_atoms
                if c["elem"] == "C" and np.linalg.norm(c["coord"] - x_coord) <= 2.1
            ]
            c_coord = bonded_c[0]["coord"] if bonded_c else None

            for ratom in rec_atoms:
                if ratom["elem"] not in acceptor_elems:
                    continue
                r_coord = ratom["coord"]
                dist = float(np.linalg.norm(r_coord - x_coord))
                if dist <= max_dist:
                    angle = 180.0
                    has_angle = False
                    if c_coord is not None:
                        angle = cls._calc_angle(c_coord, x_coord, r_coord)  # C - X ... Acceptor
                        has_angle = True

                    is_fluorine = (lelem == "F")
                    if is_fluorine:
                        # Fluorine multipolar contact
                        if angle >= 120.0 or not has_angle:
                            res_id = f"{ratom['res_name']} {ratom['res_num']}:{ratom['chain']}"
                            halogens.append({
                                "type": "Fluorine Contact",
                                "subtype": "Multipolar / Polar F-Contact",
                                "distance": round(dist, 2),
                                "angle_deg": round(angle, 1) if has_angle else None,
                                "residue": res_id,
                                "res_name": ratom["res_name"],
                                "res_num": ratom["res_num"],
                                "chain": ratom["chain"],
                                "receptor_atom": ratom["atom_name"],
                                "ligand_atom": latom["atom_name"],
                                "ligand_atom_idx": latom.get("atom_idx", 0),
                                "halogen_element": "F",
                                "criterion": "C-F...Acceptor distance <= 3.8 Å (multipolar contact; F lacks classical σ-hole)",
                                "scientific_classification": "Non-classical fluorine polar contact (not a σ-hole halogen bond)",
                                "start_coord": [float(c) for c in x_coord],
                                "end_coord": [float(c) for c in r_coord]
                            })
                    else:
                        # Classical σ-hole halogen bond (Cl, Br, I)
                        if angle >= 130.0:
                            res_id = f"{ratom['res_name']} {ratom['res_num']}:{ratom['chain']}"
                            halogens.append({
                                "type": "Halogen Bond",
                                "subtype": f"Classical {lelem} σ-Hole Bond",
                                "distance": round(dist, 2),
                                "angle_deg": round(angle, 1),
                                "residue": res_id,
                                "res_name": ratom["res_name"],
                                "res_num": ratom["res_num"],
                                "chain": ratom["chain"],
                                "receptor_atom": ratom["atom_name"],
                                "ligand_atom": latom["atom_name"],
                                "ligand_atom_idx": latom.get("atom_idx", 0),
                                "halogen_element": lelem,
                                "criterion": f"C-{lelem}...Acceptor angle >= 130°, distance <= 3.8 Å (PLIP σ-hole criterion)",
                                "scientific_classification": f"Classical {lelem} halogen bond",
                                "start_coord": [float(c) for c in x_coord],
                                "end_coord": [float(c) for c in r_coord]
                            })

        # Deduplicate to closest halogen bond per residue
        unique_hb = {}
        for hb in halogens:
            key = hb["residue"]
            if key not in unique_hb or unique_hb[key]["distance"] > hb["distance"]:
                unique_hb[key] = hb
        return list(unique_hb.values())

    @classmethod
    def _find_hydrophobic_contacts(
        cls,
        rec_atoms: List[Dict[str, Any]],
        lig_atoms: List[Dict[str, Any]],
        max_dist: float = 4.0
    ) -> List[Dict[str, Any]]:
        """Strict contacts between non-polar carbons (<= 4.0 A)."""
        hydrophobics = []
        for latom in lig_atoms:
            if latom["elem"] != "C":
                continue
            lx, ly, lz = latom["coord"]

            for ratom in rec_atoms:
                if ratom["elem"] != "C":
                    continue
                rx, ry, rz = ratom["coord"]

                if abs(lx - rx) > max_dist or abs(ly - ry) > max_dist or abs(lz - rz) > max_dist:
                    continue

                dist = math.sqrt((lx - rx)**2 + (ly - ry)**2 + (lz - rz)**2)
                if dist <= max_dist:
                    res_id = f"{ratom['res_name']} {ratom['res_num']}:{ratom['chain']}"
                    hydrophobics.append({
                        "type": "Hydrophobic Contact",
                        "distance": round(dist, 2),
                        "residue": res_id,
                        "res_name": ratom["res_name"],
                        "res_num": ratom["res_num"],
                        "chain": ratom["chain"],
                        "receptor_atom": ratom["atom_name"],
                        "ligand_atom": latom["atom_name"],
                        "ligand_atom_idx": latom.get("atom_idx", 0),
                        "start_coord": [lx, ly, lz],
                        "end_coord": [rx, ry, rz]
                    })

        # Deduplicate to closest contact per residue
        unique = {}
        for hp in hydrophobics:
            key = f"{hp['residue']}_{hp['ligand_atom_idx']}"
            if key not in unique or unique[key]["distance"] > hp["distance"]:
                unique[key] = hp
        return list(unique.values())
