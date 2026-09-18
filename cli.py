#!/usr/bin/env python3
"""
cli.py - Bindora Dock Command Line Interface (CLI)
===================================================
Run molecular docking, pocket detection, ADME profiling, and interaction analysis
directly from the terminal on your own proteins and ligands without needing a browser.

Examples:
  # 1. Dock by PDB ID and SMILES:
  python cli.py --receptor 1CX2 --ligand "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin

  # 2. Dock local PDB file with local SDF ligand:
  python cli.py --receptor ./my_target.pdb --ligand ./my_drug.sdf --out pose.pdbqt

  # 3. Fetch drug by generic name from PubChem and dock:
  python cli.py --receptor 1IEP --drug "Imatinib" --exhaustiveness 8

  # 4. Save results as JSON report:
  python cli.py --receptor 1CX2 --ligand "CC(=O)Oc1ccccc1C(=O)O" --json result.json
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Add project root to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine
from backend.services.adme import ADMEProfiler
from backend.services.bioactivity import BioactivityService

def main():
    parser = argparse.ArgumentParser(
        description="Bindora Dock CLI - Research Molecular Docking Engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cli.py --receptor 1CX2 --ligand "CC(=O)Oc1ccccc1C(=O)O" --name Aspirin
  python cli.py --receptor target.pdb --ligand drug.sdf --exhaustiveness 12 --out docked.pdbqt
  python cli.py --receptor 1IEP --drug Imatinib --adme
        """
    )

    parser.add_argument("--receptor", required=True, help="RCSB PDB ID (e.g. 1CX2) or path to .pdb/.cif file")
    parser.add_argument("--ligand", help="SMILES string, or path to .sdf/.mol/.pdbqt file")
    parser.add_argument("--drug", help="Search PubChem for drug name (e.g. Aspirin, Imatinib) to fetch SMILES")
    parser.add_argument("--name", default="Ligand", help="Display name for the ligand")
    parser.add_argument("--chain", default=None, help="Target protein chain (default: auto-detect)")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina search exhaustiveness (default: 8)")
    parser.add_argument("--modes", type=int, default=9, help="Number of binding poses to generate (default: 9)")
    parser.add_argument("--cpu", type=int, default=None, help="Number of CPU cores for AutoDock Vina (default: all available)")
    parser.add_argument("--out", help="Output file to save the top docked pose (PDBQT format)")
    parser.add_argument("--json", help="Save complete computational run metrics to a JSON file")
    parser.add_argument("--adme", action="store_true", help="Print RDKit ADME & Drug-likeness profile")
    parser.add_argument("--contacts", action="store_true", help="Print active site residue contacts & H-bonds")

    args = parser.parse_args()

    print("=" * 65)
    print("  BINDORA DOCK v2.0 - COMMAND LINE RESEARCH INTERFACE")
    print("=" * 65)

    # 1. Resolve Receptor
    receptor_input = args.receptor.strip()
    pdb_content = ""
    rec_display = receptor_input

    if Path(receptor_input).is_file():
        print(f"[*] Reading local receptor file: {receptor_input}")
        pdb_content = Path(receptor_input).read_text(encoding="utf-8", errors="replace")
    elif len(receptor_input) == 4 and receptor_input.isalnum():
        print(f"[*] Ingesting RCSB PDB structure: {receptor_input.upper()}...")
        fetched = StructureFetcher.fetch_rcsb_pdb(receptor_input)
        if not fetched or "pdb_content" not in fetched:
            print(f"[!] Error: Could not download PDB '{receptor_input}' from RCSB.", file=sys.stderr)
            sys.exit(1)
        pdb_content = fetched["pdb_content"]
        rec_display = f"{receptor_input.upper()} - {fetched.get('title', '')[:40]}"
    else:
        print(f"[!] Error: Receptor must be an existing file or a 4-letter PDB ID.", file=sys.stderr)
        sys.exit(1)

    # 2. Prepare Receptor
    print("[*] Preparing receptor (water/ion stripping, physiological protonation, pocket search)...")
    try:
        rec_prep = DockingEngine.prepare_receptor(pdb_content, target_chain=args.chain)
    except Exception as e:
        print(f"[!] Receptor preparation failed: {e}", file=sys.stderr)
        sys.exit(1)

    pocket = rec_prep.get("detected_pocket") or rec_prep.get("blind_docking_box")
    print(f"    -> Atoms: {rec_prep.get('atom_count')}, Residues: {rec_prep.get('residue_count')}")
    print(f"    -> Pocket: {pocket.get('description', 'Auto-detected')}")
    print(f"    -> Grid Center: ({pocket['center']['x']}, {pocket['center']['y']}, {pocket['center']['z']})")
    print(f"    -> Grid Box Size: ({pocket['size']['x']} x {pocket['size']['y']} x {pocket['size']['z']}) A")

    # 3. Resolve Ligand
    ligand_input = args.ligand
    is_sdf = False

    if args.drug and not ligand_input:
        print(f"[*] Searching PubChem for drug: '{args.drug}'...")
        pub_res = StructureFetcher.search_pubchem(args.drug)
        if not pub_res or "smiles" not in pub_res:
            print(f"[!] Error: Could not find compound '{args.drug}' in PubChem.", file=sys.stderr)
            sys.exit(1)
        ligand_input = pub_res["smiles"]
        args.name = pub_res.get("iupac_name") or args.drug
        print(f"    -> Found CID: {pub_res.get('cid')}, SMILES: {ligand_input}")

    if not ligand_input:
        print("[!] Error: Either --ligand (SMILES or file) or --drug (name) is required.", file=sys.stderr)
        sys.exit(1)

    if Path(ligand_input).is_file():
        ext = Path(ligand_input).suffix.lower()
        print(f"[*] Reading local ligand file: {ligand_input}")
        is_sdf = ext in (".sdf", ".mol")
        lig_raw = Path(ligand_input).read_text(encoding="utf-8", errors="replace")
    else:
        lig_raw = ligand_input

    # 4. Prepare Ligand (3D Conformer & PDBQT)
    print(f"[*] Preparing ligand '{args.name}' (MMFF94 3D optimization, Meeko flexible torsion tree)...")
    try:
        lig_prep = DockingEngine.prepare_ligand(lig_raw, is_sdf=is_sdf)
    except Exception as e:
        print(f"[!] Ligand preparation failed: {e}", file=sys.stderr)
        sys.exit(1)

    smiles = lig_prep.get("canonical_smiles") or lig_raw
    print(f"    -> Heavy Atoms: {lig_prep.get('heavy_atom_count')}, Rotatable Torsions: {lig_prep.get('rotatable_bonds')}")

    # 5. Run Docking
    print(f"[*] Executing AutoDock Vina (exhaustiveness={args.exhaustiveness}, modes={args.modes})...")
    try:
        poses = DockingEngine.run_docking(
            rec_prep["pdbqt_text"],
            lig_prep["pdbqt_text"],
            pocket["center"],
            pocket["size"],
            exhaustiveness=args.exhaustiveness,
            num_modes=args.modes,
            cpu=args.cpu
        )
    except Exception as e:
        print(f"[!] Docking execution failed: {e}", file=sys.stderr)
        sys.exit(1)

    if not poses:
        print("[!] No docking poses returned by Vina engine.", file=sys.stderr)
        sys.exit(1)

    top_pose = poses[0]
    affinity = top_pose["affinity_kcal"]

    # 6. Analyze Interactions & Thermodynamics
    mw = float(lig_prep.get("adme", {}).get("physicochemical", {}).get("molecular_weight", {}).get("value", 300.0))
    ha = int(lig_prep.get("heavy_atom_count", 20))
    thermo = BioactivityService.calculate_thermodynamics(affinity, ha, mw)
    contacts = DockingEngine.analyze_interactions(rec_prep["cleaned_pdb"], top_pose["pdbqt_content"])

    # 7. Print Results Table
    print("\n" + "=" * 65)
    print("  DOCKING RESULTS SUMMARY")
    print("=" * 65)
    print(f"Target Receptor:      {rec_display}")
    print(f"Docked Compound:      {args.name}")
    print(f"Best Affinity (dG):   {affinity:.2f} kcal/mol")
    print(f"Theoretical Kd:       {thermo['theoretical_kd_nm']} nM ({thermo['theoretical_kd_um']} uM)")
    print(f"Ligand Efficiency:    {thermo['ligand_efficiency']['value']} kcal/mol/heavy atom")
    print(f"Potency Tier:         {thermo['potency_class']}")
    print(f"Hydrogen Bonds:       {contacts['total_hbond_count']}")
    print(f"Hydrophobic Contacts: {contacts['total_hydrophobic_count']}")

    print("\n[Binding Modes Table]")
    print(f"Mode | Affinity (kcal/mol) | RMSD Lower Bound | RMSD Upper Bound")
    print("-" * 60)
    for p in poses:
        print(f"  {p['mode']:2d} | {p['affinity_kcal']:19.2f} | {p.get('rmsd_lower_bound', 0.0):16.2f} | {p.get('rmsd_upper_bound', 0.0):16.2f}")

    # Optional Contacts Print
    if args.contacts or True:
        print("\n[Interacting Active Site Residues]")
        hb_list = [f"{h['protein_residue']} ({h['distance_angstroms']} A)" for h in contacts.get("hydrogen_bonds", [])]
        print(f"  H-Bonds:     {', '.join(hb_list) if hb_list else 'None within 3.5 A'}")
        hp_list = contacts.get("interacting_residues", [])[:8]
        print(f"  Hydrophobic: {', '.join(hp_list) if hp_list else 'None'}")

    # Optional ADME Print
    if args.adme:
        adme = ADMEProfiler.calculate_adme(smiles)
        lip = adme.get("drug_likeness", {}).get("lipinski", {})
        print("\n[RDKit Physicochemical & ADME Profile]")
        print(f"  Molecular Weight: {adme['physicochemical']['molecular_weight']['value']} Da")
        print(f"  LogP:             {adme['physicochemical']['logp']['value']}")
        print(f"  H-Bond Donors:    {adme['physicochemical']['hbd']['value']}")
        print(f"  H-Bond Acceptors: {adme['physicochemical']['hba']['value']}")
        print(f"  Lipinski Rule:    {lip.get('status')} ({lip.get('violations')} violations)")
        print(f"  GI Absorption:    {adme.get('pharmacokinetics', {}).get('gi_absorption', 'Unknown')}")

    # Save Output Pose if requested
    if args.out:
        out_path = Path(args.out)
        out_path.write_text(top_pose["pdbqt_content"], encoding="utf-8")
        print(f"\n[+] Saved top binding pose to: {out_path.resolve()}")

    # Save Output JSON if requested
    if args.json:
        json_path = Path(args.json)
        full_result = {
            "receptor": receptor_input,
            "ligand_name": args.name,
            "smiles": smiles,
            "best_affinity_kcal": affinity,
            "thermodynamics": thermo,
            "interactions": contacts,
            "poses": poses
        }
        json_path.write_text(json.dumps(full_result, indent=2), encoding="utf-8")
        print(f"[+] Saved complete JSON metrics to: {json_path.resolve()}")

    print("=" * 65 + "\n")

if __name__ == "__main__":
    main()