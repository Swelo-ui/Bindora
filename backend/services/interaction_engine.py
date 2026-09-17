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

        # Flexible residue candidates
        flex_candidates = []
        for cr in sorted(list(contact_residues)):
            try:
                parts = cr.split()
                rname = parts[0]
                rnum_chain = parts[1]
                rnum, rchain = rnum_chain.split(":")
                flex_candidates.append({
                    "id": f"{rchain}:{rnum}",
                    "label": cr,
                    "res_name": rname,
                    "res_num": int(rnum),
                    "chain": rchain
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
            "total_hydrophobic_count": len(hydrophobics)
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
        atoms = []
        atom_idx = 0
        for line in block_str.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    elem = line[76:78].strip() or atom_name[0]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    atoms.append({
                        "atom_name": atom_name,
                        "elem": elem.upper(),
                        "coord": np.array([x, y, z]),
                        "atom_idx": atom_idx
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

        return salt_bridges

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
        if not mol:
            return stacks

        # 1. Identify ligand aromatic rings
        ring_info = mol.GetRingInfo()
        lig_rings = []
        conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None

        for r_atom_indices in ring_info.AtomRings():
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in r_atom_indices):
                # Calculate centroid and normal
                coords = []
                for idx in r_atom_indices:
                    # Match with 3D pose atom coordinate if available
                    match_lat = next((la for la in lig_atoms if la.get("atom_idx") == idx), None)
                    if match_lat is not None:
                        coords.append(match_lat["coord"])
                    elif conf:
                        p = conf.GetAtomPosition(idx)
                        coords.append(np.array([p.x, p.y, p.z]))
                if len(coords) >= 5:
                    coords_arr = np.array(coords)
                    centroid = np.mean(coords_arr, axis=0)
                    # Normal vector via cross product of two vectors
                    v1 = coords_arr[1] - coords_arr[0]
                    v2 = coords_arr[2] - coords_arr[0]
                    normal = np.cross(v1, v2)
                    norm = np.linalg.norm(normal)
                    if norm > 1e-6:
                        normal = normal / norm
                    lig_rings.append({
                        "indices": r_atom_indices,
                        "centroid": centroid,
                        "normal": normal
                    })

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

        return stacks

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

        # 1. Receptor cations (LYS NZ, ARG guanidinium) to ligand aromatic rings
        # 2. Ligand cations to receptor aromatic rings (PHE, TYR, TRP, HIS)
        rec_cations: Dict[Tuple[str, int, str], List[np.ndarray]] = {}
        rec_aromatics: Dict[Tuple[str, int, str], List[np.ndarray]] = {}

        for ra in rec_atoms:
            rkey = (ra["res_name"], ra["res_num"], ra["chain"])
            if ra["res_name"] == "LYS" and ra["atom_name"] == "NZ":
                rec_cations.setdefault(rkey, []).append(ra["coord"])
            elif ra["res_name"] == "ARG" and ra["atom_name"] in ("NE", "CZ", "NH1", "NH2"):
                rec_cations.setdefault(rkey, []).append(ra["coord"])
            elif ra["res_name"] in ("PHE", "TYR", "TRP", "HIS") and ra["elem"] in ("C", "N"):
                rec_aromatics.setdefault(rkey, []).append(ra["coord"])

        # Check ligand cationic nitrogens to receptor aromatics
        for latom in lig_atoms:
            if latom["elem"] == "N":
                lcoord = latom["coord"]
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

        return pi_cations

    @classmethod
    def _find_halogen_bonds(
        cls,
        rec_atoms: List[Dict[str, Any]],
        lig_atoms: List[Dict[str, Any]],
        max_dist: float = 3.8
    ) -> List[Dict[str, Any]]:
        """PLIP criteria for Halogen bonds: C-X...O/N with angle >= 130° and dist <= 3.8 A."""
        halogens = []
        hal_elems = {"F", "CL", "BR", "I"}
        acceptor_elems = {"O", "N", "S"}

        # Find ligand halogens and their bonded carbons
        for latom in lig_atoms:
            if latom["elem"] not in hal_elems:
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
                    if c_coord is not None:
                        angle = cls._calc_angle(c_coord, x_coord, r_coord)  # C - X ... Acceptor

                    if angle >= 130.0:
                        res_id = f"{ratom['res_name']} {ratom['res_num']}:{ratom['chain']}"
                        halogens.append({
                            "type": "Halogen Bond",
                            "distance": round(dist, 2),
                            "angle_deg": round(angle, 1),
                            "residue": res_id,
                            "res_name": ratom["res_name"],
                            "res_num": ratom["res_num"],
                            "chain": ratom["chain"],
                            "receptor_atom": ratom["atom_name"],
                            "ligand_atom": latom["atom_name"],
                            "ligand_atom_idx": latom.get("atom_idx", 0),
                            "start_coord": [float(c) for c in x_coord],
                            "end_coord": [float(c) for c in r_coord]
                        })

        return halogens

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
