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
        width: int = 650,
        height: int = 500
    ) -> str:
        """
        Render a 2D interaction diagram showing the central ligand surrounded by
        contact residues, hydrogen bonds (dashed lines with distance), and hydrophobic contacts.
        Supports both Dark Mode and Light Mode seamlessly via CSS variables,
        with automated collision avoidance for residue badges and no molecule clipping.
        """
        if isinstance(smiles_or_mol, str):
            mol = Chem.MolFromSmiles(smiles_or_mol)
        else:
            mol = smiles_or_mol

        if not mol:
            # Fallback simple SVG
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
        cy = (height - 30) / 2.0  # Slightly offset upward to account for bottom legend bar

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
        hydrophobics: List[Dict[str, Any]] = interactions.get("hydrophobic_contacts", [])

        # Collect contacts for radial positioning
        raw_items = []

        # 1. Collect H-Bonds (limit to top 8)
        for idx, hb in enumerate(hbonds[:8]):
            res_name = hb.get('res_name', 'RES')
            res_num = hb.get('res_num', '')
            chain = hb.get('chain', 'A')
            res_label = f"{res_name} {res_num}:{chain}"
            dist = hb.get("distance", 3.0)

            # Match target atom coord
            target_pt = atom_coords[idx % len(atom_coords)] if atom_coords else (cx, cy)
            dx = target_pt[0] - cx
            dy = target_pt[1] - cy
            angle = math.atan2(dy, dx)

            raw_items.append({
                "type": "hbond",
                "label": res_label,
                "distance": dist,
                "target_pt": target_pt,
                "angle": angle
            })

        # 2. Collect Hydrophobic Contacts (limit to top 8)
        for idx, hp in enumerate(hydrophobics[:8]):
            res_name = hp.get('res_name', 'RES')
            res_num = hp.get('res_num', '')
            chain = hp.get('chain', 'A')
            res_label = f"{res_name} {res_num}:{chain}"

            target_pt = atom_coords[(idx * 2 + 1) % len(atom_coords)] if atom_coords else (cx, cy)
            dx = target_pt[0] - cx
            dy = target_pt[1] - cy
            angle = math.atan2(dy, dx) + (idx * 0.15)  # Slight offset to reduce initial collision

            raw_items.append({
                "type": "hydrophobic",
                "label": res_label,
                "distance": None,
                "target_pt": target_pt,
                "angle": angle
            })

        # Sort all items by angle around the central ligand
        raw_items.sort(key=lambda item: item["angle"])
        n_items = len(raw_items)

        # Enforce minimum angular separation
        if n_items > 1:
            min_sep = (2.0 * math.pi) / max(n_items + 1, 8)
            for i in range(1, n_items):
                diff = raw_items[i]["angle"] - raw_items[i-1]["angle"]
                if diff < min_sep:
                    raw_items[i]["angle"] = raw_items[i-1]["angle"] + min_sep

        # Compute initial radial positions on an ellipse matching aspect ratio
        rx = min(width, height) * 0.44
        ry = min(width, height) * 0.36
        placed_nodes = []

        for item in raw_items:
            ang = item["angle"]
            px = cx + rx * math.cos(ang)
            py = cy + ry * math.sin(ang)

            # Clamp within canvas boundaries with safe margin
            px = max(65, min(width - 65, px))
            py = max(38, min(height - 54, py))

            placed_nodes.append({
                "item": item,
                "x": px,
                "y": py
            })

        # Collision relaxation passes to guarantee badges never overlap
        for _ in range(6):
            for i in range(len(placed_nodes)):
                for j in range(i + 1, len(placed_nodes)):
                    dx = placed_nodes[i]["x"] - placed_nodes[j]["x"]
                    dy = placed_nodes[i]["y"] - placed_nodes[j]["y"]
                    if abs(dx) < 86 and abs(dy) < 28:
                        overlap_y = 28 - abs(dy)
                        shift = overlap_y / 2.0 + 1.5
                        if dy >= 0:
                            placed_nodes[i]["y"] = min(height - 54, placed_nodes[i]["y"] + shift)
                            placed_nodes[j]["y"] = max(38, placed_nodes[j]["y"] - shift)
                        else:
                            placed_nodes[i]["y"] = max(38, placed_nodes[i]["y"] - shift)
                            placed_nodes[j]["y"] = min(height - 54, placed_nodes[j]["y"] + shift)

        # Render annotations (lines, distance pills, badges)
        annotations_svg = []
        for node in placed_nodes:
            item = node["item"]
            res_x = node["x"]
            res_y = node["y"]
            target_pt = item["target_pt"]
            res_label = item["label"]

            if item["type"] == "hbond":
                dist = item["distance"]
                mid_x = (target_pt[0] + res_x) / 2.0
                mid_y = (target_pt[1] + res_y) / 2.0

                # Cyan dashed vector
                annotations_svg.append(
                    f'<line x1="{target_pt[0]:.1f}" y1="{target_pt[1]:.1f}" '
                    f'x2="{res_x:.1f}" y2="{res_y:.1f}" '
                    f'stroke="var(--bd-hb-line, #06b6d4)" stroke-width="2" stroke-dasharray="5,4" opacity="0.85"/>'
                )
                # Distance pill
                annotations_svg.append(
                    f'<rect x="{mid_x - 18:.1f}" y="{mid_y - 8:.1f}" width="36" height="16" rx="4" '
                    f'fill="var(--bd-tag-bg, #042f2e)" stroke="var(--bd-tag-border, #0891b2)" stroke-width="1"/>'
                    f'<text x="{mid_x:.1f}" y="{mid_y + 3.5:.1f}" fill="var(--bd-tag-text, #67e8f9)" '
                    f'font-family="system-ui, monospace" font-size="9" font-weight="bold" text-anchor="middle">{dist}Å</text>'
                )
                # Residue badge
                annotations_svg.append(
                    f'<g transform="translate({res_x:.1f},{res_y:.1f})">'
                    f'<rect x="-42" y="-12" width="84" height="24" rx="6" '
                    f'fill="var(--bd-hb-bg, #042f2e)" stroke="var(--bd-hb-border, #0d9488)" stroke-width="1.5"/>'
                    f'<circle cx="-32" cy="0" r="4" fill="#10b981"/>'
                    f'<text x="5" y="3.8" fill="var(--bd-hb-text, #e0f2fe)" '
                    f'font-family="system-ui, -apple-system, sans-serif" font-size="10" font-weight="bold" text-anchor="middle">{res_label}</text>'
                    f'</g>'
                )
            else:
                # Hydrophobic contact arc / spoke ray
                annotations_svg.append(
                    f'<line x1="{target_pt[0]:.1f}" y1="{target_pt[1]:.1f}" '
                    f'x2="{res_x:.1f}" y2="{res_y:.1f}" '
                    f'stroke="var(--bd-hp-line, #f59e0b)" stroke-width="1.5" stroke-dasharray="2,3" opacity="0.75"/>'
                )
                # Hydrophobic badge
                annotations_svg.append(
                    f'<g transform="translate({res_x:.1f},{res_y:.1f})">'
                    f'<rect x="-40" y="-11" width="80" height="22" rx="5" '
                    f'fill="var(--bd-hp-bg, #2d1d03)" stroke="var(--bd-hp-border, #b45309)" stroke-width="1.2"/>'
                    f'<circle cx="-30" cy="0" r="3.5" fill="#f59e0b"/>'
                    f'<text x="4" y="3.5" fill="var(--bd-hp-text, #fef3c7)" '
                    f'font-family="system-ui, -apple-system, sans-serif" font-size="9.5" font-weight="600" text-anchor="middle">{res_label}</text>'
                    f'</g>'
                )

        # 3. Legend bar at bottom
        legend_y = height - 20
        legend_svg = (
            f'<g transform="translate(15, {legend_y})">'
            f'<rect x="-5" y="-14" width="{width - 30}" height="28" rx="6" '
            f'fill="var(--bd-legend-bg, #0b1329)" stroke="var(--bd-legend-border, #1e293b)" stroke-width="1"/>'
            f'<line x1="15" y1="0" x2="35" y2="0" stroke="var(--bd-hb-line, #06b6d4)" stroke-width="2" stroke-dasharray="4,3"/>'
            f'<text x="42" y="3.5" fill="var(--bd-subtext, #cbd5e1)" font-family="system-ui, sans-serif" font-size="10">Hydrogen Bond (H-Bond)</text>'
            f'<line x1="210" y1="0" x2="230" y2="0" stroke="var(--bd-hp-line, #f59e0b)" stroke-width="2" stroke-dasharray="2,3"/>'
            f'<text x="237" y="3.5" fill="var(--bd-subtext, #cbd5e1)" font-family="system-ui, sans-serif" font-size="10">Hydrophobic Contact</text>'
            f'<text x="{width - 45}" y="3.5" fill="var(--bd-legend-muted, #64748b)" font-family="system-ui, sans-serif" font-size="10" text-anchor="end">LigPlot-style 2D Map</text>'
            f'</g>'
        )

        # Embedded CSS with dual-theme variable definitions
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
            '    --bd-hb-bg: #042f2e;\n'
            '    --bd-hb-border: #0d9488;\n'
            '    --bd-hb-text: #e0f2fe;\n'
            '    --bd-hb-line: #06b6d4;\n'
            '    --bd-hp-bg: #2d1d03;\n'
            '    --bd-hp-border: #b45309;\n'
            '    --bd-hp-text: #fef3c7;\n'
            '    --bd-hp-line: #f59e0b;\n'
            '    --bd-tag-bg: #042f2e;\n'
            '    --bd-tag-border: #0891b2;\n'
            '    --bd-tag-text: #67e8f9;\n'
            '  }\n'
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
            '    --bd-hb-bg: #f0fdfa;\n'
            '    --bd-hb-border: #0d9488;\n'
            '    --bd-hb-text: #115e59;\n'
            '    --bd-hb-line: #0891b2;\n'
            '    --bd-hp-bg: #fffbeb;\n'
            '    --bd-hp-border: #d97706;\n'
            '    --bd-hp-text: #78350f;\n'
            '    --bd-hp-line: #d97706;\n'
            '    --bd-tag-bg: #e0f2fe;\n'
            '    --bd-tag-border: #0284c7;\n'
            '    --bd-tag-text: #0369a1;\n'
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
            '      --bd-hb-bg: #f0fdfa !important;\n'
            '      --bd-hb-border: #0d9488 !important;\n'
            '      --bd-hb-text: #115e59 !important;\n'
            '      --bd-hb-line: #0891b2 !important;\n'
            '      --bd-hp-bg: #fffbeb !important;\n'
            '      --bd-hp-border: #d97706 !important;\n'
            '      --bd-hp-text: #78350f !important;\n'
            '      --bd-hp-line: #d97706 !important;\n'
            '      --bd-tag-bg: #e0f2fe !important;\n'
            '      --bd-tag-border: #0284c7 !important;\n'
            '      --bd-tag-text: #0369a1 !important;\n'
            '    }\n'
            '  }\n'
            '</style>\n'
            '</defs>\n'
        )

        # Extract molecule paths from RDKit output
        header_end = base_svg.find('<!-- END OF HEADER -->')
        svg_end = base_svg.rfind('</svg>')
        if header_end != -1 and svg_end != -1:
            molecule_paths = base_svg[header_end + len('<!-- END OF HEADER -->'):svg_end].strip()
        else:
            molecule_paths = base_svg

        # Replace RDKit default black stroke with adaptive theme variable
        molecule_paths = molecule_paths.replace('stroke:#000000', 'stroke:var(--bd-bond, #e2e8f0)')

        all_annotations = "\n".join(annotations_svg)

        final_svg = (
            f'<svg class="bindora-svg-root" id="bindora-2d-interaction-svg" '
            f'xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'viewBox="0 0 {width} {height}" width="100%" height="100%" '
            f'style="background-color: var(--bd-bg, #091428); border-radius: 8px; display: block; overflow: hidden;">\n'
            f'{style_defs}\n'
            f'  <!-- Zoomable and Pannable Diagram Content -->\n'
            f'  <g id="bindora-diagram-content" transform="matrix(1 0 0 1 0 0)">\n'
            f'    <g id="bindora-ligand-layer">\n{molecule_paths}\n    </g>\n'
            f'    <g id="bindora-annotations-layer">\n{all_annotations}\n    </g>\n'
            f'  </g>\n'
            f'  <!-- Pinned Legend Bar (Always in fixed position) -->\n'
            f'  <g id="bindora-diagram-legend">\n{legend_svg}\n  </g>\n'
            f'</svg>'
        )

        return final_svg
