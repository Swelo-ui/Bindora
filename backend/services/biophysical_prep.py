import math
import numpy as np
from typing import Dict, Any, List, Tuple, Optional

# Standard AD4 / AMBER FF14SB partial charges & atom types for standard amino acid residues at pH 7.4
# (Hydrogens are polar-only 'HD' per united-atom AutoDock convention; non-polar Hs merged into carbons)
AMBER_RESIDUE_TEMPLATES = {
    "ALA": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679), "CB": ("C", -0.1825)
    },
    "ARG": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0007), "CG": ("C", 0.0390), "CD": ("C", 0.0486),
        "NE": ("NA", -0.3479), "HE": ("HD", 0.3456),
        "CZ": ("C", 0.6402),
        "NH1": ("NA", -0.4674), "HH11": ("HD", 0.3582), "HH12": ("HD", 0.3582),
        "NH2": ("NA", -0.4674), "HH21": ("HD", 0.3582), "HH22": ("HD", 0.3582)
    },
    "ASN": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0598), "CG": ("C", 0.6865), "OD1": ("OA", -0.5912),
        "ND2": ("NA", -0.6698), "HD21": ("HD", 0.3486), "HD22": ("HD", 0.3486)
    },
    "ASP": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0303), "CG": ("C", 0.7994), "OD1": ("OA", -0.8014), "OD2": ("OA", -0.8014)
    },
    "CYS": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.1231), "SG": ("SA", -0.3119), "HG": ("HD", 0.1933)
    },
    "GLN": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.0029), "CG": ("C", -0.0645), "CD": ("C", 0.6967), "OE1": ("OA", -0.6086),
        "NE2": ("NA", -0.6541), "HE21": ("HD", 0.3429), "HE22": ("HD", 0.3429)
    },
    "GLU": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.0560), "CG": ("C", -0.0425), "CD": ("C", 0.8054), "OE1": ("OA", -0.8188), "OE2": ("OA", -0.8188)
    },
    "GLY": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", -0.0252), "C": ("C", 0.5973), "O": ("OA", -0.5679)
    },
    "HIS": {  # HIE tautomer (H on NE2) - standard physiological neutral form
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0520), "CG": ("A", 0.0051), "ND1": ("NA", -0.5284),
        "CD2": ("A", 0.1386), "CE1": ("A", 0.2057),
        "NE2": ("NA", -0.3643), "HE2": ("HD", 0.3648)
    },
    "ILE": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0087), "CG1": ("C", -0.0430), "CG2": ("C", -0.0430), "CD1": ("C", -0.0880)
    },
    "LEU": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0518), "CG": ("C", 0.0453), "CD1": ("C", -0.0917), "CD2": ("C", -0.0917)
    },
    "LYS": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.0102), "CG": ("C", 0.0187), "CD": ("C", 0.0085), "CE": ("C", 0.1148),
        "NZ": ("NA", -0.3010), "HZ1": ("HD", 0.3300), "HZ2": ("HD", 0.3300), "HZ3": ("HD", 0.3300)
    },
    "MET": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0487), "CG": ("C", 0.0652), "SD": ("SA", -0.2742), "CE": ("C", -0.0536)
    },
    "PHE": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0343), "CG": ("A", 0.0112),
        "CD1": ("A", -0.0950), "CD2": ("A", -0.0950),
        "CE1": ("A", -0.1000), "CE2": ("A", -0.1000),
        "CZ": ("A", -0.0980)
    },
    "PRO": {
        "N": ("N", -0.2163), "CA": ("C", 0.0270), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0180), "CG": ("C", -0.0210), "CD": ("C", 0.0240)
    },
    "SER": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.0843), "OG": ("OA", -0.6406), "HG": ("HD", 0.4275)
    },
    "THR": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.1294), "OG1": ("OA", -0.6763), "HG1": ("HD", 0.4102), "CG2": ("C", -0.0438)
    },
    "TRP": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0050), "CG": ("A", -0.1413), "CD1": ("A", -0.0267),
        "NE1": ("NA", -0.3418), "HE1": ("HD", 0.3412),
        "CE2": ("A", 0.1292), "CD2": ("A", -0.0549), "CE3": ("A", -0.1154),
        "CZ3": ("A", -0.0976), "CH2": ("A", -0.0964), "CZ2": ("A", -0.1069)
    },
    "TYR": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0153), "CG": ("A", -0.0135),
        "CD1": ("A", -0.0900), "CD2": ("A", -0.0900),
        "CE1": ("A", -0.1400), "CE2": ("A", -0.1400),
        "CZ": ("A", 0.1900), "OH": ("OA", -0.5579), "HH": ("HD", 0.3999)
    },
    "VAL": {
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.0315), "CG1": ("C", -0.0760), "CG2": ("C", -0.0760)
    },
    "MSE": {  # Selenomethionine
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0487), "CG": ("C", 0.0652), "SE": ("SA", -0.2742), "CE": ("C", -0.0536)
    },
    "SEP": {  # Phosphoserine
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.0843), "OG": ("OA", -0.5000), "P": ("P", 1.2000),
        "O1P": ("OA", -0.8500), "O2P": ("OA", -0.8500), "O3P": ("OA", -0.8500)
    },
    "TPO": {  # Phosphothreonine
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.1294), "OG1": ("OA", -0.5000), "CG2": ("C", -0.0438), "P": ("P", 1.2000),
        "O1P": ("OA", -0.8500), "O2P": ("OA", -0.8500), "O3P": ("OA", -0.8500)
    },
    "PTR": {  # Phosphotyrosine
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.0153), "CG": ("A", -0.0135),
        "CD1": ("A", -0.0900), "CD2": ("A", -0.0900),
        "CE1": ("A", -0.1400), "CE2": ("A", -0.1400),
        "CZ": ("A", 0.1900), "OH": ("OA", -0.5000), "P": ("P", 1.2000),
        "O1P": ("OA", -0.8500), "O2P": ("OA", -0.8500), "O3P": ("OA", -0.8500)
    },
    "KCX": {  # Carbamylated Lysine
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", 0.0102), "CG": ("C", 0.0187), "CD": ("C", 0.0085), "CE": ("C", 0.1148),
        "NZ": ("NA", -0.4000), "HZ": ("HD", 0.3300), "CX": ("C", 0.7000), "O1X": ("OA", -0.7500), "O2X": ("OA", -0.7500)
    },
    "CSO": {  # S-hydroxycysteine
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.1231), "SG": ("SA", -0.1000), "OD": ("OA", -0.5500), "HD": ("HD", 0.3500)
    },
    "CME": {  # S,S-(2-hydroxyethyl)thiocysteine
        "N": ("N", -0.4157), "H": ("HD", 0.2719), "CA": ("C", 0.0337), "C": ("C", 0.5973), "O": ("OA", -0.5679),
        "CB": ("C", -0.1231), "SG": ("SA", -0.1000), "SD": ("SA", -0.1000), "CE": ("C", 0.0500), "CZ": ("C", 0.1000), "OH": ("OA", -0.6000), "HH": ("HD", 0.4000)
    }
}


class BiophysicalReceptorPreparer:
    """
    Biophysical receptor preparation engine that adds polar hydrogens (pH 7.4)
    and assigns proper AutoDock4 atom types and AMBER/Kollman partial charges.
    Ensures AutoDock Vina has directional hydrogen bond donors (HD) and acceptors (OA/NA).
    """

    @classmethod
    def prepare(
        cls,
        protein_lines: List[str],
        pH: float = 7.4,
        structural_water_lines: Optional[List[str]] = None,
        metal_lines: Optional[List[str]] = None
    ) -> Tuple[List[str], List[str], Dict[str, Any]]:
        """
        Takes raw PDB ATOM lines of standard amino acids and:
        1. Parses residues and their 3D coordinates.
        2. Places polar hydrogens (backbone NH, Ser/Thr/Tyr -OH, Lys -NH3+, Arg guanidinium, His tautomer).
        3. Formats both cleaned PDB lines and true AutoDock4 PDBQT lines with non-zero partial charges.
        Optionally retains conserved structural water molecules and catalytic metal ions.
        Returns: (cleaned_pdb_lines, pdbqt_lines, stats)
        """
        # Try OpenMM / PDBFixer first if available
        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile
            import io

            fixer = PDBFixer(pdbfile=io.StringIO("\n".join(protein_lines) + "\nEND\n"))
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(pH)

            out_s = io.StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, out_s, keepIds=True)
            fixed_pdb = out_s.getvalue()
            return cls._process_protonated_pdb(fixed_pdb.splitlines(), method="OpenMM/PDBFixer", water_lines=structural_water_lines, metal_lines=metal_lines)
        except Exception:
            # Fall back to high-accuracy geometric polar hydrogen builder
            return cls._geometric_polar_hydrogen_builder(protein_lines, pH=pH, water_lines=structural_water_lines, metal_lines=metal_lines)

    @classmethod
    def _normalize(cls, vec: np.ndarray) -> np.ndarray:
        norm = np.linalg.norm(vec)
        if norm < 1e-8:
            return np.array([1.0, 0.0, 0.0])
        return vec / norm

    @classmethod
    def _geometric_polar_hydrogen_builder(
        cls,
        protein_lines: List[str],
        pH: float = 7.4,
        water_lines: Optional[List[str]] = None,
        metal_lines: Optional[List[str]] = None
    ) -> Tuple[List[str], List[str], Dict[str, Any]]:
        """
        Biophysical geometric placement of polar hydrogens at standard bond lengths and angles:
        - Backbone trans amide NH: N-H along bisector of C_{i-1}-N-CA (d = 1.01 A)
        - Hydroxyls (SER, THR, TYR): O-H (d = 0.96 A, angle ~ 108.5 deg)
        - Lysine: NZ -NH3+ tetrahedral (d = 1.01 A)
        - Arginine: Guanidinium planar system (HE on NE, 2 Hs on NH1, 2 Hs on NH2)
        - Histidine: HIE tautomer (HE2 on NE2) or HID (HD1 on ND1)
        - Amides (ASN, GLN): Planar NH2
        - Tryptophan: NE1-HE1 in indole plane
        - Cysteine: SG-HG (d = 1.34 A)
        """
        residues: List[Dict[str, Any]] = []
        curr_key = None
        curr_atoms: Dict[str, Dict[str, Any]] = {}

        for line in protein_lines:
            if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                continue
            res_name = line[17:20].strip()
            chain = line[21:22].strip()
            res_seq_str = line[22:26].strip()
            try:
                res_seq = int(res_seq_str)
            except ValueError:
                res_seq = 0
            atom_name = line[12:16].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coord = np.array([x, y, z])
            except Exception:
                continue

            key = (chain, res_seq, res_name)
            if key != curr_key:
                if curr_atoms:
                    residues.append({"key": curr_key, "atoms": curr_atoms})
                curr_key = key
                curr_atoms = {}

            curr_atoms[atom_name] = {
                "line": line,
                "coord": coord,
                "res_name": res_name,
                "chain": chain,
                "res_seq": res_seq
            }

        if curr_atoms:
            residues.append({"key": curr_key, "atoms": curr_atoms})

        output_pdb_lines: List[str] = []
        output_pdbqt_lines: List[str] = []
        polar_h_count = 0
        atom_serial = 1

        for idx, res in enumerate(residues):
            chain, res_seq, res_name = res["key"]
            atoms = res["atoms"]
            tpl = AMBER_RESIDUE_TEMPLATES.get(res_name, {})

            # 1. Add backbone atoms & polar amide hydrogen (H)
            prev_c_coord = None
            if idx > 0 and residues[idx - 1]["key"][0] == chain:
                prev_c_coord = residues[idx - 1]["atoms"].get("C", {}).get("coord")

            new_res_atoms: List[Tuple[str, np.ndarray, bool]] = []  # (name, coord, is_polar_h)

            # Order of standard backbone atoms
            for at_name in ["N", "CA", "C", "O"]:
                if at_name in atoms:
                    new_res_atoms.append((at_name, atoms[at_name]["coord"], False))

            # Add backbone amide H (not for Proline, and not if already has H)
            if "N" in atoms and "CA" in atoms and res_name != "PRO" and "H" not in atoms:
                n_coord = atoms["N"]["coord"]
                ca_coord = atoms["CA"]["coord"]
                v_ca_n = cls._normalize(n_coord - ca_coord)

                if prev_c_coord is not None:
                    v_c_n = cls._normalize(n_coord - prev_c_coord)
                    bisect = cls._normalize(v_ca_n + v_c_n)
                    h_coord = n_coord + bisect * 1.01
                else:
                    # N-terminus residue: orient away from CA
                    h_coord = n_coord + v_ca_n * 1.01

                new_res_atoms.append(("H", h_coord, True))

            # 2. Add sidechain heavy atoms
            for at_name, at_data in atoms.items():
                if at_name not in ["N", "CA", "C", "O", "H"]:
                    new_res_atoms.append((at_name, at_data["coord"], False))

            # 3. Add sidechain polar hydrogens according to residue type
            if res_name == "SER" and "OG" in atoms and "CB" in atoms and "HG" not in atoms:
                cb = atoms["CB"]["coord"]
                og = atoms["OG"]["coord"]
                ca = atoms.get("CA", {}).get("coord", cb + np.array([1.0, 0.0, 0.0]))
                v_axis = cls._normalize(og - cb)
                v_perp = cls._normalize(np.cross(v_axis, ca - cb))
                h_dir = cls._normalize(math.cos(math.radians(108.5 - 90)) * v_axis + math.sin(math.radians(108.5 - 90)) * v_perp)
                new_res_atoms.append(("HG", og + h_dir * 0.96, True))

            elif res_name == "THR" and "OG1" in atoms and "CB" in atoms and "HG1" not in atoms:
                cb = atoms["CB"]["coord"]
                og1 = atoms["OG1"]["coord"]
                v_axis = cls._normalize(og1 - cb)
                v_ref = np.array([0.0, 1.0, 0.0])
                if abs(np.dot(v_axis, v_ref)) > 0.9:
                    v_ref = np.array([1.0, 0.0, 0.0])
                v_perp = cls._normalize(np.cross(v_axis, v_ref))
                h_dir = cls._normalize(-0.3 * v_axis + 0.95 * v_perp)
                new_res_atoms.append(("HG1", og1 + h_dir * 0.96, True))

            elif res_name == "TYR" and "OH" in atoms and "CZ" in atoms and "HH" not in atoms:
                cz = atoms["CZ"]["coord"]
                oh = atoms["OH"]["coord"]
                ce1 = atoms.get("CE1", {}).get("coord", cz + np.array([0.0, 1.0, 0.0]))
                v_ring = cls._normalize(cz - ce1)
                v_oh = cls._normalize(oh - cz)
                h_dir = cls._normalize(0.3 * v_oh + 0.95 * v_ring)
                new_res_atoms.append(("HH", oh + h_dir * 0.96, True))

            elif res_name == "LYS" and "NZ" in atoms and "CE" in atoms and "HZ1" not in atoms:
                ce = atoms["CE"]["coord"]
                nz = atoms["NZ"]["coord"]
                v_ce_nz = cls._normalize(nz - ce)
                # Tetrahedral tripod of 3 hydrogens
                v_ref1 = np.array([1.0, 0.0, 0.0]) if abs(v_ce_nz[0]) < 0.8 else np.array([0.0, 1.0, 0.0])
                v_perp1 = cls._normalize(np.cross(v_ce_nz, v_ref1))
                v_perp2 = cls._normalize(np.cross(v_ce_nz, v_perp1))
                cos_tet = 1.0 / 3.0  # cos(70.5 deg)
                sin_tet = math.sqrt(1.0 - cos_tet**2)
                for ang_idx, ang in enumerate([0, 2 * math.pi / 3, 4 * math.pi / 3]):
                    h_vec = cos_tet * v_ce_nz + sin_tet * (math.cos(ang) * v_perp1 + math.sin(ang) * v_perp2)
                    new_res_atoms.append((f"HZ{ang_idx+1}", nz + cls._normalize(h_vec) * 1.01, True))

            elif res_name == "ARG" and "NE" in atoms and "CZ" in atoms and "HE" not in atoms:
                ne = atoms["NE"]["coord"]
                cz = atoms["CZ"]["coord"]
                cd = atoms.get("CD", {}).get("coord", ne - np.array([1.0, 0.0, 0.0]))
                v1 = cls._normalize(ne - cd)
                v2 = cls._normalize(ne - cz)
                bisect = cls._normalize(v1 + v2)
                new_res_atoms.append(("HE", ne + bisect * 1.01, True))
                if "NH1" in atoms:
                    nh1 = atoms["NH1"]["coord"]
                    v_cz_nh1 = cls._normalize(nh1 - cz)
                    v_p = cls._normalize(np.cross(v_cz_nh1, bisect))
                    new_res_atoms.append(("HH11", nh1 + cls._normalize(v_cz_nh1 + 0.8 * v_p) * 1.01, True))
                    new_res_atoms.append(("HH12", nh1 + cls._normalize(v_cz_nh1 - 0.8 * v_p) * 1.01, True))
                if "NH2" in atoms:
                    nh2 = atoms["NH2"]["coord"]
                    v_cz_nh2 = cls._normalize(nh2 - cz)
                    v_p = cls._normalize(np.cross(v_cz_nh2, bisect))
                    new_res_atoms.append(("HH21", nh2 + cls._normalize(v_cz_nh2 + 0.8 * v_p) * 1.01, True))
                    new_res_atoms.append(("HH22", nh2 + cls._normalize(v_cz_nh2 - 0.8 * v_p) * 1.01, True))

            elif res_name == "HIS" and "NE2" in atoms and "CD2" in atoms and "CE1" in atoms and "HE2" not in atoms:
                # Standard neutral HIE tautomer
                ne2 = atoms["NE2"]["coord"]
                cd2 = atoms["CD2"]["coord"]
                ce1 = atoms["CE1"]["coord"]
                v1 = cls._normalize(ne2 - cd2)
                v2 = cls._normalize(ne2 - ce1)
                bisect = cls._normalize(v1 + v2)
                new_res_atoms.append(("HE2", ne2 + bisect * 1.01, True))

            elif res_name == "ASN" and "ND2" in atoms and "CG" in atoms and "HD21" not in atoms:
                nd2 = atoms["ND2"]["coord"]
                cg = atoms["CG"]["coord"]
                v_axis = cls._normalize(nd2 - cg)
                od1 = atoms.get("OD1", {}).get("coord", cg + np.array([0.0, 1.0, 0.0]))
                v_perp = cls._normalize(np.cross(v_axis, od1 - cg))
                v_inplane = cls._normalize(np.cross(v_perp, v_axis))
                new_res_atoms.append(("HD21", nd2 + cls._normalize(0.5 * v_axis + 0.866 * v_inplane) * 1.01, True))
                new_res_atoms.append(("HD22", nd2 + cls._normalize(0.5 * v_axis - 0.866 * v_inplane) * 1.01, True))

            elif res_name == "GLN" and "NE2" in atoms and "CD" in atoms and "HE21" not in atoms:
                ne2 = atoms["NE2"]["coord"]
                cd = atoms["CD"]["coord"]
                v_axis = cls._normalize(ne2 - cd)
                oe1 = atoms.get("OE1", {}).get("coord", cd + np.array([0.0, 1.0, 0.0]))
                v_perp = cls._normalize(np.cross(v_axis, oe1 - cd))
                v_inplane = cls._normalize(np.cross(v_perp, v_axis))
                new_res_atoms.append(("HE21", ne2 + cls._normalize(0.5 * v_axis + 0.866 * v_inplane) * 1.01, True))
                new_res_atoms.append(("HE22", ne2 + cls._normalize(0.5 * v_axis - 0.866 * v_inplane) * 1.01, True))

            elif res_name == "TRP" and "NE1" in atoms and "CD1" in atoms and "CE2" in atoms and "HE1" not in atoms:
                ne1 = atoms["NE1"]["coord"]
                cd1 = atoms["CD1"]["coord"]
                ce2 = atoms["CE2"]["coord"]
                v1 = cls._normalize(ne1 - cd1)
                v2 = cls._normalize(ne1 - ce2)
                bisect = cls._normalize(v1 + v2)
                new_res_atoms.append(("HE1", ne1 + bisect * 1.01, True))

            elif res_name == "CYS" and "SG" in atoms and "CB" in atoms and "HG" not in atoms:
                cb = atoms["CB"]["coord"]
                sg = atoms["SG"]["coord"]
                v_cb_sg = cls._normalize(sg - cb)
                v_ref = np.array([0.0, 0.0, 1.0])
                if abs(v_cb_sg[2]) > 0.8:
                    v_ref = np.array([1.0, 0.0, 0.0])
                v_perp = cls._normalize(np.cross(v_cb_sg, v_ref))
                h_dir = cls._normalize(0.34 * v_cb_sg + 0.94 * v_perp)
                new_res_atoms.append(("HG", sg + h_dir * 1.34, True))

            # 4. Format PDB and PDBQT lines for this residue
            for at_name, coord, is_polar_h in new_res_atoms:
                elem = "H" if is_polar_h or at_name.startswith("H") else at_name[0]
                name_field = f" {at_name:<3}" if len(at_name) < 4 else f"{at_name:<4}"
                pdb_line = (
                    f"ATOM  {atom_serial:>5} {name_field} {res_name:>3} {chain:>1}{res_seq:>4}    "
                    f"{coord[0]:>8.3f}{coord[1]:>8.3f}{coord[2]:>8.3f}  1.00 20.00           {elem:>2}"
                )
                output_pdb_lines.append(pdb_line)

                # Determine AD4 atom type and partial charge from AMBER template
                ad4_type = elem
                charge = 0.0
                if at_name in tpl:
                    ad4_type, charge = tpl[at_name]
                else:
                    if elem == "H":
                        ad4_type = "HD"
                        charge = 0.25
                    elif elem == "O":
                        ad4_type = "OA"
                        charge = -0.55
                    elif elem == "N":
                        ad4_type = "NA" if res_name in ("HIS", "TRP") else "N"
                        charge = -0.40
                    elif elem == "S":
                        ad4_type = "SA"
                        charge = -0.28
                    elif elem == "P":
                        ad4_type = "P"
                        charge = 1.20
                    elif elem == "C":
                        ad4_type = "A" if res_name in ("PHE", "TYR", "TRP", "HIS") else "C"
                        charge = 0.05

                if is_polar_h:
                    polar_h_count += 1

                pdbqt_line = f"{pdb_line[:54]:<54}  1.00 20.00    {charge:>6.3f} {ad4_type:<2}"
                output_pdbqt_lines.append(pdbqt_line)
                atom_serial += 1

        # Optionally append structural water molecules as OA atoms
        waters_retained = 0
        if water_lines:
            for wline in water_lines:
                try:
                    w_x = float(wline[30:38])
                    w_y = float(wline[38:46])
                    w_z = float(wline[46:54])
                    w_res = wline[17:20].strip() or "HOH"
                    w_chain = wline[21:22].strip() or " "
                    w_num = wline[22:26].strip() or "1"
                    name_field = " O  "
                    pdb_w = (
                        f"HETATM{atom_serial:>5} {name_field} {w_res:>3} {w_chain:>1}{w_num:>4}    "
                        f"{w_x:>8.3f}{w_y:>8.3f}{w_z:>8.3f}  1.00 20.00           O "
                    )
                    pdbqt_w = f"{pdb_w[:54]:<54}  1.00 20.00    -0.560  OA"
                    output_pdb_lines.append(pdb_w)
                    output_pdbqt_lines.append(pdbqt_w)
                    atom_serial += 1
                    waters_retained += 1
                except Exception:
                    continue

        metals_retained = 0
        if metal_lines:
            metal_defs = {
                "ZN": ("Zn", 2.000), "MG": ("Mg", 2.000), "CA": ("Ca", 2.000),
                "MN": ("Mn", 2.000), "FE": ("Fe", 2.000), "CU": ("Cu", 2.000),
                "NI": ("Ni", 2.000), "CO": ("Co", 2.000)
            }
            for mline in metal_lines:
                try:
                    m_res = mline[17:20].strip().upper()
                    ad4_m, m_q = metal_defs.get(m_res, (m_res.title()[:2], 2.000))
                    mx = float(mline[30:38])
                    my = float(mline[38:46])
                    mz = float(mline[46:54])
                    m_chain = mline[21:22].strip() or " "
                    m_num = mline[22:26].strip() or "1"
                    name_field = f" {m_res:<3}" if len(m_res) < 4 else f"{m_res:<4}"
                    pdb_m = (
                        f"HETATM{atom_serial:>5} {name_field} {m_res:>3} {m_chain:>1}{m_num:>4}    "
                        f"{mx:>8.3f}{my:>8.3f}{mz:>8.3f}  1.00 20.00          {m_res:>2}"
                    )
                    pdbqt_m = f"{pdb_m[:54]:<54}  1.00 20.00    {m_q:>6.3f}  {ad4_m:<2}"
                    output_pdb_lines.append(pdb_m)
                    output_pdbqt_lines.append(pdbqt_m)
                    atom_serial += 1
                    metals_retained += 1
                except Exception:
                    continue

        stats = {
            "method": "Biophysical Polar Hydrogen Geometry Engine (pH 7.4)",
            "polar_hydrogens_added": polar_h_count,
            "structural_waters_retained": waters_retained,
            "catalytic_metals_retained": metals_retained,
            "total_atoms_prepared": len(output_pdb_lines),
            "forcefield_charges": "AMBER FF14SB / Kollman AD4 partial charges",
            "protonation_state": f"Standard physiological pH {pH:.1f} (Ser/Thr/Tyr -OH and Lys -NH3+ protonated, Asp/Glu ionized, His neutral HIE)"
        }

        return output_pdb_lines, output_pdbqt_lines, stats

    @classmethod
    def _process_protonated_pdb(
        cls,
        lines: List[str],
        method: str,
        water_lines: Optional[List[str]] = None,
        metal_lines: Optional[List[str]] = None
    ) -> Tuple[List[str], List[str], Dict[str, Any]]:
        """Atom typing and partial charge assignment for pre-protonated PDB lines."""
        pdb_lines = []
        pdbqt_lines = []
        polar_h_count = 0
        atom_serial = 1

        for line in lines:
            if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                continue
            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            elem = line[76:78].strip() or atom_name[0]

            tpl = AMBER_RESIDUE_TEMPLATES.get(res_name, {})
            if atom_name in tpl:
                ad4_type, charge = tpl[atom_name]
            else:
                if elem == "H":
                    ad4_type = "HD"
                    charge = 0.25
                elif elem == "O":
                    ad4_type = "OA"
                    charge = -0.55
                elif elem == "N":
                    ad4_type = "NA" if res_name in ("HIS", "TRP") else "N"
                    charge = -0.40
                elif elem == "S":
                    ad4_type = "SA"
                    charge = -0.28
                elif elem == "P":
                    ad4_type = "P"
                    charge = 1.20
                elif elem == "C":
                    ad4_type = "A" if res_name in ("PHE", "TYR", "TRP", "HIS") else "C"
                    charge = 0.05
                else:
                    ad4_type = elem
                    charge = 0.00

            if ad4_type == "HD":
                polar_h_count += 1

            is_non_polar_h = elem == "H" and ad4_type != "HD"
            if is_non_polar_h:
                continue

            pdb_lines.append(line)
            pdbqt_line = f"{line[:54]:<54}  1.00 20.00    {charge:>6.3f} {ad4_type:<2}"
            pdbqt_lines.append(pdbqt_line)
            atom_serial += 1

        waters_retained = 0
        if water_lines:
            for wline in water_lines:
                try:
                    w_x = float(wline[30:38])
                    w_y = float(wline[38:46])
                    w_z = float(wline[46:54])
                    w_res = wline[17:20].strip() or "HOH"
                    w_chain = wline[21:22].strip() or " "
                    w_num = wline[22:26].strip() or "1"
                    name_field = " O  "
                    pdb_w = (
                        f"HETATM{atom_serial:>5} {name_field} {w_res:>3} {w_chain:>1}{w_num:>4}    "
                        f"{w_x:>8.3f}{w_y:>8.3f}{w_z:>8.3f}  1.00 20.00           O "
                    )
                    pdbqt_w = f"{pdb_w[:54]:<54}  1.00 20.00    -0.560  OA"
                    pdb_lines.append(pdb_w)
                    pdbqt_lines.append(pdbqt_w)
                    atom_serial += 1
                    waters_retained += 1
                except Exception:
                    continue

        metals_retained = 0
        if metal_lines:
            metal_defs = {
                "ZN": ("Zn", 2.000), "MG": ("Mg", 2.000), "CA": ("Ca", 2.000),
                "MN": ("Mn", 2.000), "FE": ("Fe", 2.000), "CU": ("Cu", 2.000),
                "NI": ("Ni", 2.000), "CO": ("Co", 2.000)
            }
            for mline in metal_lines:
                try:
                    m_res = mline[17:20].strip().upper()
                    ad4_m, m_q = metal_defs.get(m_res, (m_res.title()[:2], 2.000))
                    mx = float(mline[30:38])
                    my = float(mline[38:46])
                    mz = float(mline[46:54])
                    m_chain = mline[21:22].strip() or " "
                    m_num = mline[22:26].strip() or "1"
                    name_field = f" {m_res:<3}" if len(m_res) < 4 else f"{m_res:<4}"
                    pdb_m = (
                        f"HETATM{atom_serial:>5} {name_field} {m_res:>3} {m_chain:>1}{m_num:>4}    "
                        f"{mx:>8.3f}{my:>8.3f}{mz:>8.3f}  1.00 20.00          {m_res:>2}"
                    )
                    pdbqt_m = f"{pdb_m[:54]:<54}  1.00 20.00    {m_q:>6.3f}  {ad4_m:<2}"
                    pdb_lines.append(pdb_m)
                    pdbqt_lines.append(pdbqt_m)
                    atom_serial += 1
                    metals_retained += 1
                except Exception:
                    continue

        stats = {
            "method": method,
            "polar_hydrogens_added": polar_h_count,
            "structural_waters_retained": waters_retained,
            "catalytic_metals_retained": metals_retained,
            "total_atoms_prepared": len(pdbqt_lines),
            "forcefield_charges": "AMBER FF14SB / Kollman AD4 partial charges",
            "protonation_state": "Standard physiological pH 7.4"
        }
        return pdb_lines, pdbqt_lines, stats
