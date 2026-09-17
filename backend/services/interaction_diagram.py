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
        width: int = 700,
        height: int = 520
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
                f'<rect width="100%" height="100%" fill="#091428" rx="8"/>'
                f'<text x="50%" y="50%" fill="#94a3b8" font-family="system-ui, sans-serif" font-size="14" text-anchor="middle">'
                f'Structure unavailable for 2D diagram</text></svg>'
            )

        # Ensure 2D coordinates exist
        mol_copy = Chem.Mol(mol)
        try:
            AllChem.Compute2DCoords(mol_copy)
        except Exception:
            pass

        # Prepare RDKit SVG Drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.padding = 0.28  # Ensures ample margin around molecule for interaction badges
        opts.bondLineWidth = 2.4
        opts.minFontSize = 11
        opts.maxFontSize = 16

        drawer.DrawMolecule(mol_copy)
        drawer.FinishDrawing()
        base_svg = drawer.GetDrawingText()

        # Extract canvas center
        cx = width / 2.0
        cy = (height - 35) / 2.0  # Slightly offset upward for bottom legend bar

        # Get atom drawing positions
        num_atoms = mol_copy.GetNumAtoms()
        atom_coords = []
        for i in range(num_atoms):
            try:
                pt = drawer.GetDrawCoords(i)
                atom_coords.append((pt.x, pt.y))
            except Exception:
                atom_coords.append((cx, cy))

        hbonds: List[Dict[str, Any]] = interactions.get("hydrogen_bonds", [])
        salt_bridges: List[Dict[str, Any]] = interactions.get("salt_bridges", [])
        pi_stacks: List[Dict[str, Any]] = interactions.get("pi_stacking", [])
        pi_cations: List[Dict[str, Any]] = interactions.get("pi_cation", [])
        halogens: List[Dict[str, Any]] = interactions.get("halogen_bonds", [])
        hydrophobics: List[Dict[str, Any]] = interactions.get("hydrophobic_contacts", [])

        # Color definitions for all 6 interaction types (Scientific standard palette)
        interaction_styles = {
            "hbond": {
                "name": "H-Bond",
                "color": "#facc15",          # Gold / Yellow (Universally recognized for H-bonds)
                "line_dash": "5,4",
                "line_width": "2.0",
                "pill_bg": "var(--bd-hb-bg, #2a2004)",
                "pill_border": "var(--bd-hb-border, #ca8a04)",
                "pill_text": "var(--bd-hb-text, #fde047)",
                "badge_border": "#ca8a04"
            },
            "salt_bridge": {
                "name": "Salt Bridge",
                "color": "#ec4899",          # Magenta / Deep Pink (Ionic ion-pair standard)
                "line_dash": "4,4",
                "line_width": "2.0",
                "pill_bg": "#380a24",
                "pill_border": "#db2777",
                "pill_text": "#f472b6",
                "badge_border": "#db2777"
            },
            "pi_stack": {
                "name": "π-π Stacking",
                "color": "#10b981",          # Emerald / Forest Green (Aromatic stacking standard)
                "line_dash": "4,4",
                "line_width": "2.0",
                "pill_bg": "#062c1d",
                "pill_border": "#059669",
                "pill_text": "#6ee7b7",
                "badge_border": "#059669"
            },
            "pi_cation": {
                "name": "π-Cation",
                "color": "#f97316",          # Warm Amber / Orange (Cation-aromatic standard)
                "line_dash": "4,4",
                "line_width": "1.8",
                "pill_bg": "#381604",
                "pill_border": "#ea580c",
                "pill_text": "#fdba74",
                "badge_border": "#ea580c"
            },
            "halogen": {
                "name": "Halogen Bond",
                "color": "#a855f7",          # Violet / Purple (Halogen sigma-hole standard)
                "line_dash": "4,4",
                "line_width": "1.8",
                "pill_bg": "#2a0845",
                "pill_border": "#9333ea",
                "pill_text": "#d8b4fe",
                "badge_border": "#9333ea"
            },
            "hydrophobic": {
                "name": "Hydrophobic",
                "color": "#94a3b8",          # Neutral Cool Slate / Charcoal Gray (Non-polar carbon standard)
                "line_dash": "2,3",
                "line_width": "1.5",
                "pill_bg": "#1e293b",
                "pill_border": "#64748b",
                "pill_text": "#cbd5e1",
                "badge_border": "#64748b"
            }
        }

        # Step 1: Collect and group interactions by unique residue (res_name, res_num, chain)
        # In structural biology (PoseView, LigPlot+), contacts to the same amino acid converge on that residue
        grouped_residues: Dict[str, Dict[str, Any]] = {}

        def get_target_pt(latom_idx: Optional[int], fallback_idx: int) -> tuple:
            if latom_idx is not None and 0 <= latom_idx < len(atom_coords):
                return atom_coords[latom_idx]
            elif atom_coords:
                return atom_coords[fallback_idx % len(atom_coords)]
            return (cx, cy)

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
                    "res_name": res_name,
                    "res_num": res_num,
                    "chain": chain,
                    "label": f"{res_name} {res_num}:{chain}",
                    "contacts": [],
                    "types_set": set(),
                    "target_pts": []
                }

            grouped_residues[key]["contacts"].append({
                "type": itype,
                "distance": dist,
                "target_pt": target_pt,
                "latom_idx": latom_idx
            })
            grouped_residues[key]["types_set"].add(itype)
            grouped_residues[key]["target_pts"].append(target_pt)

        # Populate from top contacts (prioritizing strong directional bonds first)
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
                f'<rect width="100%" height="100%" fill="#091428" rx="8"/>'
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

            dx = avg_x - cx
            dy = avg_y - cy
            angle = math.atan2(dy, dx)
            r["angle"] = angle

            # Hierarchy priority for primary badge theme:
            priority_order = ["salt_bridge", "hbond", "pi_stack", "pi_cation", "halogen", "hydrophobic"]
            primary_type = "hydrophobic"
            for ptype in priority_order:
                if ptype in r["types_set"]:
                    primary_type = ptype
                    break
            r["primary_type"] = primary_type

        # Sort circularly by angle around the ligand
        residue_list.sort(key=lambda r: r["angle"])
        n_res = len(residue_list)

        # Enforce angular separation to prevent initial overlap
        if n_res > 1:
            min_angular_sep = (2.0 * math.pi) / max(n_res + 1, 10)
            for i in range(1, n_res):
                diff = residue_list[i]["angle"] - residue_list[i-1]["angle"]
                if diff < min_angular_sep:
                    residue_list[i]["angle"] = residue_list[i-1]["angle"] + min_angular_sep

        # Multi-Tier Radial Distribution
        base_rx = min(width, height) * 0.40
        base_ry = min(width, height) * 0.33
        tier_offsets = [0.0, 48.0, 85.0]

        nodes = []
        for idx, r in enumerate(residue_list):
            tier = idx % 2 if n_res <= 10 else idx % 3
            radial_offset = tier_offsets[tier]

            ang = r["angle"]
            rx = base_rx + radial_offset * 1.1
            ry = base_ry + radial_offset * 0.9

            px = cx + rx * math.cos(ang)
            py = cy + ry * math.sin(ang)

            px = max(60.0, min(width - 60.0, px))
            py = max(34.0, min(height - 52.0, py))

            nodes.append({
                "res": r,
                "x": px,
                "y": py,
                "angle": ang
            })

        # Step 3: Full 2D AABB Box Collision Relaxation (Both X and Y axes)
        box_w = 94.0
        box_h = 32.0

        for _ in range(30):
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    dx = nodes[i]["x"] - nodes[j]["x"]
                    dy = nodes[i]["y"] - nodes[j]["y"]
                    overlap_x = box_w - abs(dx)
                    overlap_y = box_h - abs(dy)

                    if overlap_x > 0 and overlap_y > 0:
                        if overlap_x < overlap_y:
                            shift = overlap_x / 2.0 + 1.5
                            sign = 1.0 if dx >= 0 else -1.0
                            nodes[i]["x"] += shift * sign
                            nodes[j]["x"] -= shift * sign
                        else:
                            shift = overlap_y / 2.0 + 1.5
                            sign = 1.0 if dy >= 0 else -1.0
                            nodes[i]["y"] += shift * sign
                            nodes[j]["y"] -= shift * sign

                        nodes[i]["x"] = max(58.0, min(width - 58.0, nodes[i]["x"]))
                        nodes[i]["y"] = max(32.0, min(height - 50.0, nodes[i]["y"]))
                        nodes[j]["x"] = max(58.0, min(width - 58.0, nodes[j]["x"]))
                        nodes[j]["y"] = max(32.0, min(height - 50.0, nodes[j]["y"]))

        # Step 4: Render annotations (lines, distance pills, badges)
        contact_elements_svg = []
        badge_elements_svg = []
        placed_pills = []

        for node in nodes:
            r = node["res"]
            res_x = node["x"]
            res_y = node["y"]
            res_label = r["label"]
            primary_style = interaction_styles[r["primary_type"]]
            types_class_str = " ".join([f"itype-{t}" for t in r["types_set"]])
            types_data_str = " ".join(r["types_set"])

            for c in r["contacts"]:
                ctype = c["type"]
                dist = c["distance"]
                target_pt = c["target_pt"]
                cstyle = interaction_styles[ctype]

                vx = res_x - target_pt[0]
                vy = res_y - target_pt[1]
                line_len = math.hypot(vx, vy)
                if line_len < 1.0:
                    line_len = 1.0

                t_start = min(0.12, 10.0 / line_len)
                t_end = max(0.85, (line_len - 44.0) / line_len)

                lx1 = target_pt[0] + vx * t_start
                ly1 = target_pt[1] + vy * t_start
                lx2 = target_pt[0] + vx * t_end
                ly2 = target_pt[1] + vy * t_end

                line_svg = (
                    f'<line class="interaction-line" x1="{lx1:.1f}" y1="{ly1:.1f}" x2="{lx2:.1f}" y2="{ly2:.1f}" '
                    f'stroke="{cstyle["color"]}" stroke-width="{cstyle["line_width"]}" '
                    f'stroke-dasharray="{cstyle["line_dash"]}" opacity="0.9"/>'
                )

                pill_svg = ""
                if line_len >= 52.0 and dist is not None:
                    t_pill = 0.44
                    px = target_pt[0] + vx * t_pill
                    py = target_pt[1] + vy * t_pill

                    d_res = math.hypot(px - res_x, py - res_y)
                    if d_res < 46.0:
                        px = res_x - (vx / line_len) * 46.0
                        py = res_y - (vy / line_len) * 46.0

                    d_target = math.hypot(px - target_pt[0], py - target_pt[1])
                    if d_target < 22.0:
                        px = target_pt[0] + (vx / line_len) * 22.0
                        py = target_pt[1] + (vy / line_len) * 22.0

                    for prev_p in placed_pills:
                        if abs(px - prev_p["x"]) < 36.0 and abs(py - prev_p["y"]) < 18.0:
                            norm_x = -vy / line_len
                            norm_y = vx / line_len
                            px += norm_x * 12.0
                            py += norm_y * 12.0
                            break

                    placed_pills.append({"x": px, "y": py})

                    pill_svg = (
                        f'<g class="interaction-distance-pill">'
                        f'<rect x="{px - 17.0:.1f}" y="{py - 7.5:.1f}" width="34" height="15" rx="3.5" '
                        f'fill="{cstyle["pill_bg"]}" stroke="{cstyle["pill_border"]}" stroke-width="1"/>'
                        f'<text x="{px:.1f}" y="{py + 3.5:.1f}" fill="{cstyle["pill_text"]}" '
                        f'font-family="system-ui, monospace" font-size="8.5" font-weight="bold" text-anchor="middle">{dist}Å</text>'
                        f'</g>'
                    )

                contact_elements_svg.append(
                    f'<g class="interaction-contact itype-{ctype}" data-type="{ctype}">\n'
                    f'  {line_svg}\n'
                    f'  {pill_svg}\n'
                    f'</g>'
                )

            dots_svg = []
            priority_order = ["salt_bridge", "hbond", "pi_stack", "pi_cation", "halogen", "hydrophobic"]
            sorted_types = sorted(list(r["types_set"]), key=lambda t: priority_order.index(t) if t in priority_order else 99)
            dot_start_x = -32.0
            for d_idx, dtype in enumerate(sorted_types[:3]):
                d_color = interaction_styles[dtype]["color"]
                dx_pos = dot_start_x + (d_idx * 7.5)
                dots_svg.append(f'<circle cx="{dx_pos:.1f}" cy="0" r="3.2" fill="{d_color}"/>')

            text_offset_x = 4.0 if len(sorted_types) <= 1 else 7.0

            badge_svg = (
                f'<g class="interaction-badge-node {types_class_str}" data-types="{types_data_str}" transform="translate({res_x:.1f},{res_y:.1f})">\n'
                f'  <rect x="-42" y="-12" width="84" height="24" rx="6" '
                f'fill="{primary_style["pill_bg"]}" stroke="{primary_style["badge_border"]}" stroke-width="1.4"/>\n'
                f'  {"".join(dots_svg)}\n'
                f'  <text x="{text_offset_x}" y="3.8" fill="{primary_style["pill_text"]}" '
                f'font-family="system-ui, -apple-system, sans-serif" font-size="9.5" font-weight="bold" text-anchor="middle">{res_label}</text>\n'
                f'</g>'
            )
            badge_elements_svg.append(badge_svg)

        legend_y = height - 20
        legend_svg = (
            f'<g id="bindora-diagram-legend" transform="translate(15, {legend_y})">'
            f'<rect x="-5" y="-13" width="{width - 25}" height="26" rx="5" '
            f'fill="var(--bd-legend-bg, #0b1329)" stroke="var(--bd-legend-border, #1e293b)" stroke-width="1"/>'
            f'<line x1="10" y1="0" x2="28" y2="0" stroke="#facc15" stroke-width="2" stroke-dasharray="4,3"/>'
            f'<text x="34" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9.5">Hydrogen Bond</text>'
            f'<line x1="95" y1="0" x2="113" y2="0" stroke="#ec4899" stroke-width="2" stroke-dasharray="4,3"/>'
            f'<text x="119" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9.5">Salt Bridge</text>'
            f'<line x1="195" y1="0" x2="213" y2="0" stroke="#10b981" stroke-width="2" stroke-dasharray="4,3"/>'
            f'<text x="219" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9.5">π-π Stack</text>'
            f'<line x1="288" y1="0" x2="306" y2="0" stroke="#f97316" stroke-width="1.8" stroke-dasharray="4,3"/>'
            f'<text x="312" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9.5">π-Cation</text>'
            f'<line x1="375" y1="0" x2="393" y2="0" stroke="#a855f7" stroke-width="1.8" stroke-dasharray="4,3"/>'
            f'<text x="399" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9.5">Halogen</text>'
            f'<line x1="460" y1="0" x2="478" y2="0" stroke="#94a3b8" stroke-width="1.5" stroke-dasharray="2,3"/>'
            f'<text x="484" y="3.5" fill="#cbd5e1" font-family="system-ui, sans-serif" font-size="9.5">Hydrophobic (Slate)</text>'
            f'<text x="{width - 40}" y="3.5" fill="var(--bd-legend-muted, #64748b)" font-family="system-ui, sans-serif" font-size="9.5" text-anchor="end">LigPlot-style 2D Map</text>'
            f'</g>'
        )

        style_defs = (
            '<defs>\n'
            '<style>\n'
            '  .bindora-svg-root {\n'
            '    --bd-bg: #091428;\n'
            '    --bd-bond: #e2e8f0;\n'
            '    --bd-text: #f1f5f9;\n'
            '    --bd-subtext: #cbd5e1;\n'
            '    --bd-legend-muted: #64748b;\n'
            '    --bd-legend-bg: #0b1329;\n'
            '    --bd-legend-border: #1e293b;\n'
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
            f'style="background-color: var(--bd-bg, #091428); border-radius: 8px; display: block; overflow: hidden;">\n'
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
