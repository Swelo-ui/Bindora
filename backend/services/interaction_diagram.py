import math
from typing import Dict, Any, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

class InteractionDiagramGenerator:
    """Generates 2D LigPlot-style protein-ligand interaction diagrams as scalable, dual-theme SVG."""

    @staticmethod
    def generate_diagram_svg(
        smiles_or_mol: Any,
        interactions: Dict[str, Any],
        width: int = 720,
        height: int = 540
    ) -> str:
        """
        Render a 2D interaction diagram showing the central ligand surrounded by
        contact residues, hydrogen bonds (dashed lines with distance), salt bridges,
        pi-stacking, pi-cation, halogen bonds, and hydrophobic contacts.
        Supports both Dark Mode and Light Mode seamlessly via CSS variables,
        with multi-tier collision avoidance, zero-overlap residue badges, and interactive toggles.
        """
        if isinstance(smiles_or_mol, str):
            mol = Chem.MolFromSmiles(smiles_or_mol)
        else:
            mol = smiles_or_mol

        if not mol:
            return (
                f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">'
                f'<rect width="100%" height="100%" fill="#090a0f" rx="8"/>'
                f'<text x="50%" y="50%" fill="#94a3b8" font-family="system-ui, sans-serif" font-size="14" text-anchor="middle">'
                f'Structure unavailable for 2D diagram</text></svg>'
            )

        # Ensure 2D coordinates exist
        mol_copy = Chem.Mol(mol)
        try:
            AllChem.Compute2DCoords(mol_copy)
        except Exception:
            pass

        # Prepare RDKit SVG Drawer with ample perimeter margin (30%)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.padding = 0.30
        opts.bondLineWidth = 2.4
        opts.minFontSize = 11
        opts.maxFontSize = 16

        drawer.DrawMolecule(mol_copy)
        drawer.FinishDrawing()
        base_svg = drawer.GetDrawingText()

        # Canvas center and bottom legend offset
        legend_y = height - 20
        cx = width / 2.0
        cy = 245.0

        # Extract precise 2D atom coordinates
        num_atoms = mol_copy.GetNumAtoms()
        atom_coords = []
        for i in range(num_atoms):
            try:
                pt = drawer.GetDrawCoords(i)
                atom_coords.append((pt.x, pt.y))
            except Exception:
                atom_coords.append((cx, cy))

        if atom_coords:
            min_atom_x = min(p[0] for p in atom_coords)
            max_atom_x = max(p[0] for p in atom_coords)
            min_atom_y = min(p[1] for p in atom_coords)
            max_atom_y = max(p[1] for p in atom_coords)
            mol_w = max(10.0, max_atom_x - min_atom_x)
            mol_h = max(10.0, max_atom_y - min_atom_y)
            lig_cx = (min_atom_x + max_atom_x) / 2.0
            lig_cy = (min_atom_y + max_atom_y) / 2.0
        else:
            min_atom_x, max_atom_x = cx - 50, cx + 50
            min_atom_y, max_atom_y = cy - 50, cy + 50
            mol_w, mol_h = 100.0, 100.0
            lig_cx, lig_cy = cx, cy

        hbonds: List[Dict[str, Any]] = interactions.get("hydrogen_bonds", [])
        salt_bridges: List[Dict[str, Any]] = interactions.get("salt_bridges", [])
        pi_stacks: List[Dict[str, Any]] = interactions.get("pi_stacking", [])
        pi_cations: List[Dict[str, Any]] = interactions.get("pi_cation", [])
        halogens: List[Dict[str, Any]] = interactions.get("halogen_bonds", [])
        hydrophobics: List[Dict[str, Any]] = interactions.get("hydrophobic_contacts", [])

        # Color definitions for all 6 interaction types (Scientific standard palette)
        interaction_styles = {
            "hbond": {
                "name": "Hydrogen Bond",
                "color": "#facc15",          # Gold / Yellow
                "line_dash": "5,4",
                "line_width": "2.0",
                "pill_bg": "var(--bd-hb-bg, #2a2004)",
                "pill_border": "var(--bd-hb-border, #ca8a04)",
                "pill_text": "var(--bd-hb-text, #fde047)",
                "badge_border": "#ca8a04"
            },
            "salt_bridge": {
                "name": "Salt Bridge",
                "color": "#ec4899",          # Magenta / Deep Pink
                "line_dash": "4,4",
                "line_width": "2.0",
                "pill_bg": "#380a24",
                "pill_border": "#db2777",
                "pill_text": "#f472b6",
                "badge_border": "#db2777"
            },
            "pi_stack": {
                "name": "π-π Stacking",
                "color": "#10b981",          # Emerald / Forest Green
                "line_dash": "4,4",
                "line_width": "2.0",
                "pill_bg": "#062c1d",
                "pill_border": "#059669",
                "pill_text": "#6ee7b7",
                "badge_border": "#059669"
            },
            "pi_cation": {
                "name": "π-Cation",
                "color": "#f97316",          # Warm Amber / Orange
                "line_dash": "4,4",
                "line_width": "1.8",
                "pill_bg": "#381604",
                "pill_border": "#ea580c",
                "pill_text": "#fdba74",
                "badge_border": "#ea580c"
            },
            "halogen": {
                "name": "Halogen Bond",
                "color": "#a855f7",          # Violet / Purple
                "line_dash": "4,4",
                "line_width": "1.8",
                "pill_bg": "#2a0845",
                "pill_border": "#9333ea",
                "pill_text": "#d8b4fe",
                "badge_border": "#9333ea"
            },
            "hydrophobic": {
                "name": "Hydrophobic",
                "color": "#94a3b8",          # Neutral Cool Slate
                "line_dash": "2,3",
                "line_width": "1.5",
                "pill_bg": "#1e293b",
                "pill_border": "#64748b",
                "pill_text": "#cbd5e1",
                "badge_border": "#64748b"
            }
        }

        # Step 1: Collect and group interactions by unique residue (res_name, res_num, chain)
        grouped_residues: Dict[str, Dict[str, Any]] = {}

        def get_target_pt(latom_idx: Optional[int], fallback_idx: int) -> tuple:
            if latom_idx is not None and 0 <= latom_idx < len(atom_coords):
                return atom_coords[latom_idx]
            elif atom_coords:
                return atom_coords[fallback_idx % len(atom_coords)]
            return (lig_cx, lig_cy)

        def add_contact(itype: str, data: Dict[str, Any], fallback_idx: int, default_dist: Optional[float] = None):
            res_name = data.get('res_name', 'RES')
            res_num = data.get('res_num', '')
            chain = data.get('chain', 'A')
            dist = data.get('distance', default_dist)
            latom_idx = data.get('ligand_atom_idx')
            target_pt = get_target_pt(latom_idx, fallback_idx)

            key = f"{res_name}_{res_num}_{chain}"
            if key not in grouped_residues:
                grouped_residues[key] = {
                    "key": key,
                    "res_name": res_name,
                    "res_num": res_num,
                    "chain": chain,
                    "label": f"{res_name} {res_num}:{chain}",
                    "contacts": [],
                    "types_set": set(),
                    "target_pts": []
                }

            # If hydrophobic and the residue already has directional polar/aromatic contacts, skip
            if itype == "hydrophobic" and len(grouped_residues[key]["types_set"]) > 0:
                return

            # For hydrophobic residues, strictly enforce AT MOST ONE contact (closest one) to avoid clutter
            if itype == "hydrophobic" and "hydrophobic" in grouped_residues[key]["types_set"]:
                existing = [c for c in grouped_residues[key]["contacts"] if c["type"] == "hydrophobic"][0]
                if dist is not None and (existing["distance"] is None or dist < existing["distance"]):
                    existing["distance"] = dist
                    existing["target_pt"] = target_pt
                    existing["latom_idx"] = latom_idx
                return

            # For directional contacts, limit to at most 2 contacts per residue
            if itype != "hydrophobic" and len([c for c in grouped_residues[key]["contacts"] if c["type"] == itype]) >= 2:
                return

            grouped_residues[key]["contacts"].append({
                "type": itype,
                "distance": dist,
                "target_pt": target_pt,
                "latom_idx": latom_idx
            })
            grouped_residues[key]["types_set"].add(itype)
            grouped_residues[key]["target_pts"].append(target_pt)

        # Populate from top contacts
        for idx, hb in enumerate(hbonds[:8]):
            add_contact("hbond", hb, idx, default_dist=3.0)

        for idx, sb in enumerate(salt_bridges[:4]):
            add_contact("salt_bridge", sb, idx, default_dist=3.8)

        for idx, ps in enumerate(pi_stacks[:4]):
            add_contact("pi_stack", ps, idx, default_dist=4.5)

        for idx, pc in enumerate(pi_cations[:3]):
            add_contact("pi_cation", pc, idx, default_dist=4.2)

        for idx, hal in enumerate(halogens[:3]):
            add_contact("halogen", hal, idx, default_dist=3.4)

        for idx, hp in enumerate(hydrophobics[:8]):
            add_contact("hydrophobic", hp, idx, default_dist=None)

        if not grouped_residues:
            return (
                f'<svg class="bindora-svg-root" width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">'
                f'<rect width="100%" height="100%" fill="#090a0f" rx="8"/>'
                f'<text x="50%" y="50%" fill="#94a3b8" font-family="system-ui, sans-serif" font-size="14" text-anchor="middle">'
                f'No contacts within threshold for 2D diagram</text></svg>'
            )

        # Step 2: Calculate initial radial placement and multi-tier staggering
        residue_list = list(grouped_residues.values())

        for r in residue_list:
            pts = r["target_pts"]
            avg_x = sum(p[0] for p in pts) / len(pts)
            avg_y = sum(p[1] for p in pts) / len(pts)
            r["avg_target"] = (avg_x, avg_y)

            dx = avg_x - lig_cx
            dy = avg_y - lig_cy
            r["angle"] = math.atan2(dy, dx)

            priority_order = ["salt_bridge", "hbond", "pi_stack", "pi_cation", "halogen", "hydrophobic"]
            primary_type = "hydrophobic"
            for ptype in priority_order:
                if ptype in r["types_set"]:
                    primary_type = ptype
                    break
            r["primary_type"] = primary_type

            # Dynamic badge dimensions with generous text and dot separation
            n_dots = min(len(r["types_set"]), 4)
            text_width = len(r["label"]) * 7.2
            if n_dots > 0:
                badge_w = max(94.0, 24.0 + (n_dots - 1) * 9.0 + 10.0 + text_width + 14.0)
            else:
                badge_w = max(76.0, text_width + 24.0)
            r["badge_w"] = badge_w
            r["half_w"] = badge_w / 2.0

        # Sort circularly by angle around the ligand center
        residue_list.sort(key=lambda r: r["angle"])
        n_res = len(residue_list)

        # Enforce minimum angular separation
        if n_res > 1:
            min_angular_sep = (2.0 * math.pi) / max(n_res + 1, 10)
            for i in range(1, n_res):
                diff = residue_list[i]["angle"] - residue_list[i-1]["angle"]
                if diff < min_angular_sep:
                    residue_list[i]["angle"] = residue_list[i-1]["angle"] + min_angular_sep

        # Multi-tier radial distribution: ray from ligand center outward past molecular boundary
        tier_offsets = [0.0, 32.0] if n_res <= 8 else [0.0, 28.0, 52.0]

        nodes = []
        for idx, r in enumerate(residue_list):
            ang = r["angle"]
            tier = idx % len(tier_offsets)
            r_tier = tier_offsets[tier]

            # Ray distance from ligand center to molecular envelope along angle
            r_lig = math.sqrt((mol_w * 0.5 * math.cos(ang))**2 + (mol_h * 0.5 * math.sin(ang))**2) + 56.0
            r_total = r_lig + r_tier

            px = lig_cx + r_total * math.cos(ang)
            py = lig_cy + r_total * math.sin(ang)

            # Guaranteed canvas insetting: badge border can NEVER clip SVG boundary
            min_x = r["half_w"] + 14.0
            max_x = width - r["half_w"] - 14.0
            min_y = 24.0
            max_y = legend_y - 26.0

            px = max(min_x, min(max_x, px))
            py = max(min_y, min(max_y, py))

            nodes.append({
                "res": r,
                "x": px,
                "y": py,
                "min_x": min_x,
                "max_x": max_x,
                "min_y": min_y,
                "max_y": max_y
            })

        # Step 3: Exact 2D AABB Box Collision Relaxation for Badges (35 passes)
        for _ in range(35):
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    dx = nodes[i]["x"] - nodes[j]["x"]
                    dy = nodes[i]["y"] - nodes[j]["y"]
                    req_w = nodes[i]["res"]["half_w"] + nodes[j]["res"]["half_w"] + 10.0
                    req_h = 30.0
                    overlap_x = req_w - abs(dx)
                    overlap_y = req_h - abs(dy)

                    if overlap_x > 0 and overlap_y > 0:
                        if overlap_x < overlap_y:
                            shift = overlap_x / 2.0 + 1.0
                            sgn = 1.0 if dx >= 0 else -1.0
                            nodes[i]["x"] += shift * sgn
                            nodes[j]["x"] -= shift * sgn
                        else:
                            shift = overlap_y / 2.0 + 1.0
                            sgn = 1.0 if dy >= 0 else -1.0
                            nodes[i]["y"] += shift * sgn
                            nodes[j]["y"] -= shift * sgn

                        nodes[i]["x"] = max(nodes[i]["min_x"], min(nodes[i]["max_x"], nodes[i]["x"]))
                        nodes[i]["y"] = max(nodes[i]["min_y"], min(nodes[i]["max_y"], nodes[i]["y"]))
                        nodes[j]["x"] = max(nodes[j]["min_x"], min(nodes[j]["max_x"], nodes[j]["x"]))
                        nodes[j]["y"] = max(nodes[j]["min_y"], min(nodes[j]["max_y"], nodes[j]["y"]))

        # Step 4: Distance Pill Optimizer & Annotation Rendering
        contact_elements_svg = []
        badge_elements_svg = []
        placed_pills = []

        def find_best_distance_pill_position(res_x: float, res_y: float, half_w: float, target_pt: tuple, latom_idx: Optional[int]):
            vx = target_pt[0] - res_x
            vy = target_pt[1] - res_y
            line_len = math.hypot(vx, vy)
            if line_len < 32.0:
                return None

            dir_x = vx / line_len
            dir_y = vy / line_len
            norm_x = -dir_y
            norm_y = dir_x

            best_cand = None
            best_penalty = 1e9

            # Evaluate 63 candidate positions around the line (7 longitudinal x 9 lateral offsets)
            for t in [0.20, 0.28, 0.38, 0.48, 0.58, 0.68, 0.78]:
                base_x = res_x + dir_x * (line_len * t)
                base_y = res_y + dir_y * (line_len * t)

                for lat in [0.0, 18.0, -18.0, 28.0, -28.0, 38.0, -38.0, 50.0, -50.0]:
                    cpx = base_x + norm_x * lat
                    cpy = base_y + norm_y * lat

                    # Check canvas boundaries with padding
                    if cpx < 22.0 or cpx > width - 22.0 or cpy < 18.0 or cpy > legend_y - 18.0:
                        continue

                    penalty = 0.0

                    # 1. Box collision test with the residue badge
                    if abs(cpx - res_x) < (half_w + 12.0) and abs(cpy - res_y) < 22.0:
                        penalty += 60000.0

                    # 2. Distance to target interacting atom (must be >= 22 px)
                    d_tpt = math.hypot(cpx - target_pt[0], cpy - target_pt[1])
                    if d_tpt < 22.0:
                        penalty += 4000.0 * (22.0 - d_tpt)

                    # 3. Distance to ALL other ligand atoms (strictly prevent landing inside aromatic rings or on bonds)
                    for a_idx, ap in enumerate(atom_coords):
                        if a_idx == latom_idx:
                            continue
                        d_atom = math.hypot(cpx - ap[0], cpy - ap[1])
                        if d_atom < 22.0:
                            penalty += 5000.0 * (22.0 - d_atom)

                    # 4. Box clearance from all previously placed distance pills (must be >= 38px X and >= 18px Y)
                    for prev in placed_pills:
                        dx = abs(cpx - prev['x'])
                        dy = abs(cpy - prev['y'])
                        if dx < 38.0 and dy < 18.0:
                            penalty += 80000.0 * (1.0 + (38.0 - dx) + (18.0 - dy))

                    # Minor preference for centered t and smaller lateral displacement
                    penalty += abs(lat) * 1.5 + abs(t - 0.45) * 10.0

                    if penalty < best_penalty:
                        best_penalty = penalty
                        best_cand = (cpx, cpy)

            # If all candidates have collision penalties, suppress visual pill to prevent overlapping numbers
            if best_penalty >= 20000.0:
                return None

            return best_cand

        for node in nodes:
            r = node["res"]
            res_x = node["x"]
            res_y = node["y"]
            half_w = r["half_w"]
            res_label = r["label"]
            primary_style = interaction_styles[r["primary_type"]]
            types_class_str = " ".join([f"itype-{t}" for t in r["types_set"]])
            types_data_str = " ".join(r["types_set"])

            for c_idx, c in enumerate(r["contacts"]):
                ctype = c["type"]
                dist = c["distance"]
                target_pt = c["target_pt"]
                cstyle = interaction_styles[ctype]

                # Vector from ligand atom to residue badge center
                vx = res_x - target_pt[0]
                vy = res_y - target_pt[1]
                line_len = math.hypot(vx, vy)
                if line_len < 1.0:
                    line_len = 1.0

                dir_x = vx / line_len
                dir_y = vy / line_len

                # Line starts right outside the ligand atom circle
                lx1 = target_pt[0] + dir_x * 12.0
                ly1 = target_pt[1] + dir_y * 12.0

                # Line terminates precisely on the residue badge box border
                dx_to_target = target_pt[0] - res_x
                dy_to_target = target_pt[1] - res_y
                if abs(dx_to_target) > 1e-4 and abs(dy_to_target) > 1e-4:
                    scale_to_edge = min((half_w + 2.0) / abs(dx_to_target), 14.0 / abs(dy_to_target))
                    lx2 = res_x + dx_to_target * scale_to_edge
                    ly2 = res_y + dy_to_target * scale_to_edge
                else:
                    lx2 = res_x
                    ly2 = res_y

                dist_str = f"{dist:.2f}Å" if dist is not None else ""
                tooltip_title = f"{cstyle['name']}: {dist_str} to {res_label}" if dist_str else f"{cstyle['name']} to {res_label}"

                line_svg = (
                    f'<line class="interaction-line" data-res-key="{r["key"]}" data-target-x="{target_pt[0]:.1f}" data-target-y="{target_pt[1]:.1f}" '
                    f'x1="{lx1:.1f}" y1="{ly1:.1f}" x2="{lx2:.1f}" y2="{ly2:.1f}" '
                    f'stroke="{cstyle["color"]}" stroke-width="{cstyle["line_width"]}" '
                    f'stroke-dasharray="{cstyle["line_dash"]}" opacity="0.9">'
                    f'<title>{tooltip_title}</title></line>'
                )

                pill_svg = ""
                # Distance pills are displayed ONLY on directional/polar bonds (H-bonds, salt bridges, halogens, pi)
                # Hydrophobic contacts remain clean non-polar dashed lines with hover tooltip (matching LigPlot & PoseView)
                is_directional = ctype in ("hbond", "salt_bridge", "halogen", "pi_stack", "pi_cation")
                if is_directional and dist is not None and line_len >= 34.0:
                    best_pos = find_best_distance_pill_position(res_x, res_y, half_w, target_pt, c.get("latom_idx"))
                    if best_pos:
                        px, py = best_pos
                        placed_pills.append({"x": px, "y": py})
                        pill_svg = (
                            f'<g class="interaction-distance-pill" data-dist-type="{ctype}" data-res-key="{r["key"]}" data-target-x="{target_pt[0]:.1f}" data-target-y="{target_pt[1]:.1f}">\n'
                            f'  <rect x="{px - 17.0:.1f}" y="{py - 7.5:.1f}" width="34" height="15" rx="3.5" '
                            f'fill="{cstyle["pill_bg"]}" stroke="{cstyle["pill_border"]}" stroke-width="1"/>\n'
                            f'  <text x="{px:.1f}" y="{py + 3.5:.1f}" fill="{cstyle["pill_text"]}" '
                            f'font-family="system-ui, monospace" font-size="8.5" font-weight="bold" text-anchor="middle">{dist_str}</text>\n'
                            f'  <title>{cstyle["name"]}: {dist_str}</title>\n'
                            f'</g>'
                        )

                contact_elements_svg.append(
                    f'<g class="interaction-contact itype-{ctype}" data-type="{ctype}" data-res-key="{r["key"]}">\n'
                    f'  {line_svg}\n'
                    f'  {pill_svg}\n'
                    f'</g>'
                )

            # Build residue badge with multi-type interaction indicators
            dots_svg = []
            priority_order = ["salt_bridge", "hbond", "pi_stack", "pi_cation", "halogen", "hydrophobic"]
            sorted_types = sorted(list(r["types_set"]), key=lambda t: priority_order.index(t) if t in priority_order else 99)

            n_dots = min(len(sorted_types), 4)
            badge_w = r["badge_w"]
            half_w = r["half_w"]

            if n_dots > 0:
                dot_start_x = -half_w + 12.0
                for d_idx, dtype in enumerate(sorted_types[:4]):
                    d_color = interaction_styles[dtype]["color"]
                    dx_pos = dot_start_x + (d_idx * 9.0)
                    dots_svg.append(f'<circle class="badge-dot" data-dot-type="{dtype}" cx="{dx_pos:.1f}" cy="0" r="3.2" fill="{d_color}"/>')

                # Guaranteed clean 10px separation after the last dot so text NEVER sticks to or intrudes into dots
                last_dot_right = dot_start_x + (n_dots - 1) * 9.0 + 3.2
                text_x = last_dot_right + 10.0
                text_anchor = "start"
            else:
                text_x = 0.0
                text_anchor = "middle"

            badge_svg = (
                f'<g class="interaction-badge-node {types_class_str}" data-res-key="{r["key"]}" data-types="{types_data_str}" data-primary-type="{r["primary_type"]}" transform="translate({res_x:.1f},{res_y:.1f})" style="cursor: grab;">\n'
                f'  <rect class="badge-bg" x="{-half_w:.1f}" y="-12" width="{badge_w:.1f}" height="24" rx="6" '
                f'fill="{primary_style["pill_bg"]}" stroke="{primary_style["badge_border"]}" stroke-width="1.4"/>\n'
                f'  {"".join(dots_svg)}\n'
                f'  <text class="badge-text" x="{text_x:.1f}" y="3.8" fill="{primary_style["pill_text"]}" '
                f'font-family="system-ui, -apple-system, sans-serif" font-size="9.5" font-weight="bold" text-anchor="{text_anchor}">{res_label}</text>\n'
                f'  <title>{res_label} ({", ".join(r["types_set"])})</title>\n'
                f'</g>'
            )
            badge_elements_svg.append(badge_svg)

        legend_y = height - 20
        legend_svg = (
            f'<g id="bindora-diagram-legend" transform="translate(15, {legend_y})">'
            f'<rect x="-5" y="-13" width="{width - 20}" height="26" rx="5" '
            f'fill="var(--bd-legend-bg, #121319)" stroke="var(--bd-legend-border, #22242e)" stroke-width="1"/>'
            f'<line x1="8" y1="0" x2="22" y2="0" stroke="#facc15" stroke-width="2" stroke-dasharray="4,3"/>'
            f'<text x="26" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9"><title>Hydrogen Bond</title>H-Bond</text>'
            f'<line x1="82" y1="0" x2="96" y2="0" stroke="#ec4899" stroke-width="2" stroke-dasharray="4,3"/>'
            f'<text x="100" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9">Salt Bridge</text>'
            f'<line x1="168" y1="0" x2="182" y2="0" stroke="#10b981" stroke-width="2" stroke-dasharray="4,3"/>'
            f'<text x="186" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9">π-π Stack</text>'
            f'<line x1="250" y1="0" x2="264" y2="0" stroke="#f97316" stroke-width="1.8" stroke-dasharray="4,3"/>'
            f'<text x="268" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9">π-Cation</text>'
            f'<line x1="324" y1="0" x2="338" y2="0" stroke="#a855f7" stroke-width="1.8" stroke-dasharray="4,3"/>'
            f'<text x="342" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9">Halogen</text>'
            f'<line x1="398" y1="0" x2="412" y2="0" stroke="#94a3b8" stroke-width="1.5" stroke-dasharray="2,3"/>'
            f'<text x="416" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9">Hydrophobic (Slate)</text>'
            f'<text x="{width - 32}" y="3.5" fill="var(--bd-legend-muted, #64748b)" font-family="system-ui, sans-serif" font-size="9" text-anchor="end">LigPlot-style 2D Map</text>'
            f'</g>'
        )

        style_defs = (
            '<defs>\n'
            '<style>\n'
            '  .bindora-svg-root {\n'
            '    --bd-bg: #090a0f;\n'
            '    --bd-bond: #e2e8f0;\n'
            '    --bd-text: #f1f5f9;\n'
            '    --bd-subtext: #cbd5e1;\n'
            '    --bd-legend-muted: #64748b;\n'
            '    --bd-legend-bg: #121319;\n'
            '    --bd-legend-border: #22242e;\n'
            '    --bd-hb-bg: #2a2004;\n'
            '    --bd-hb-border: #ca8a04;\n'
            '    --bd-hb-text: #fef08a;\n'
            '    --bd-hb-line: #facc15;\n'
            '    --bd-hp-bg: #1e293b;\n'
            '    --bd-hp-border: #475569;\n'
            '    --bd-hp-text: #e2e8f0;\n'
            '    --bd-hp-line: #94a3b8;\n'
            '  }\n'
            '  .bindora-svg-root.hide-hbond .itype-hbond { display: none !important; }\n'
            '  .bindora-svg-root.hide-salt_bridge .itype-salt_bridge { display: none !important; }\n'
            '  .bindora-svg-root.hide-pi_stack .itype-pi_stack { display: none !important; }\n'
            '  .bindora-svg-root.hide-pi_cation .itype-pi_cation { display: none !important; }\n'
            '  .bindora-svg-root.hide-halogen .itype-halogen { display: none !important; }\n'
            '  .bindora-svg-root.hide-hydrophobic .itype-hydrophobic { display: none !important; }\n'
            '  html.light .bindora-svg-root,\n'
            '  .theme-light .bindora-svg-root,\n'
            '  .bindora-svg-root.theme-light {\n'
            '    --bd-bg: #ffffff;\n'
            '    --bd-bond: #1e293b;\n'
            '    --bd-text: #0f172a;\n'
            '    --bd-subtext: #334155;\n'
            '    --bd-legend-muted: #94a3b8;\n'
            '    --bd-legend-bg: #f8fafc;\n'
            '    --bd-legend-border: #cbd5e1;\n'
            '    --bd-hb-bg: #fefce8;\n'
            '    --bd-hb-border: #ca8a04;\n'
            '    --bd-hb-text: #713f12;\n'
            '    --bd-hb-line: #ca8a04;\n'
            '    --bd-hp-bg: #f1f5f9;\n'
            '    --bd-hp-border: #94a3b8;\n'
            '    --bd-hp-text: #334155;\n'
            '    --bd-hp-line: #64748b;\n'
            '  }\n'
            '  @media print {\n'
            '    .bindora-svg-root {\n'
            '      --bd-bg: #ffffff !important;\n'
            '      --bd-bond: #1e293b !important;\n'
            '      --bd-text: #0f172a !important;\n'
            '      --bd-subtext: #334155 !important;\n'
            '      --bd-legend-muted: #94a3b8 !important;\n'
            '      --bd-legend-bg: #f8fafc !important;\n'
            '      --bd-legend-border: #cbd5e1 !important;\n'
            '      --bd-hb-bg: #fefce8 !important;\n'
            '      --bd-hb-border: #ca8a04 !important;\n'
            '      --bd-hb-text: #713f12 !important;\n'
            '      --bd-hb-line: #ca8a04 !important;\n'
            '      --bd-hp-bg: #f1f5f9 !important;\n'
            '      --bd-hp-border: #94a3b8 !important;\n'
            '      --bd-hp-text: #334155 !important;\n'
            '      --bd-hp-line: #64748b !important;\n'
            '    }\n'
            '  }\n'
            '</style>\n'
            '</defs>\n'
        )

        header_end = base_svg.find('<!-- END OF HEADER -->')
        svg_end = base_svg.rfind('</svg>')
        if header_end != -1 and svg_end != -1:
            molecule_paths = base_svg[header_end + len('<!-- END OF HEADER -->'):svg_end].strip()
        else:
            molecule_paths = base_svg

        molecule_paths = molecule_paths.replace('stroke:#000000', 'stroke:var(--bd-bond, #e2e8f0)')

        all_contacts = "\n".join(contact_elements_svg)
        all_badges = "\n".join(badge_elements_svg)

        final_svg = (
            f'<svg class="bindora-svg-root" id="bindora-2d-interaction-svg" '
            f'xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'viewBox="0 0 {width} {height}" width="100%" height="100%" '
            f'style="background-color: var(--bd-bg, #090a0f); border-radius: 8px; display: block; overflow: hidden;">\n'
            f'{style_defs}\n'
            f'  <!-- Zoomable and Pannable Diagram Content -->\n'
            f'  <g id="bindora-diagram-content" transform="matrix(1 0 0 1 0 0)">\n'
            f'    <g id="bindora-ligand-layer">\n{molecule_paths}\n    </g>\n'
            f'    <g id="bindora-contacts-layer">\n{all_contacts}\n    </g>\n'
            f'    <g id="bindora-badges-layer">\n{all_badges}\n    </g>\n'
            f'  </g>\n'
            f'  <!-- Pinned Legend Bar (Always in fixed position) -->\n'
            f'  {legend_svg}\n'
            f'</svg>'
        )

        return final_svg
