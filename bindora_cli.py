#!/usr/bin/env python3
"""
Bindora Dock v2.0 — Research-Grade Computational Molecular Docking CLI
Official Interactive & Scriptable Command-Line Interface.
Features Arrow-Key Navigation (TUI), direct CLI arguments, and publication-grade output.
"""

import os
import sys
import time
import json
import argparse
from pathlib import Path

# 1. Initialize Colorama to ensure Windows CMD parses ANSI colors natively (No raw ←[96m codes)
try:
    import colorama
    colorama.init(autoreset=False)
except Exception:
    pass

# Ensure UTF-8 output encoding on Windows terminals (cp1252 fix)
if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
if sys.stderr and hasattr(sys.stderr, "reconfigure"):
    try:
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

# Color codes for terminal
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    BG_CYAN = '\033[46m\033[30m'
    BG_GREEN = '\033[42m\033[30m'
    BG_BLUE = '\033[44m\033[37m'
    RESET = '\033[0m'

BANNER = rf"""{Colors.CYAN}{Colors.BOLD}
================================================================================
   ____  _           _                   ____             _      {Colors.YELLOW}[v2.0]{Colors.CYAN}
  | __ )(_)_ __   __| | ___  _ __ __ _  |  _ \  ___   ___| | __  
  |  _ \| | '_ \ / _` |/ _ \| '__/ _` | | | | |/ _ \ / __| |/ /  
  | |_) | | | | | (_| | (_) | | | (_| | | |_| | (_) | (__|   <   
  |____/|_|_| |_|\__,_|\___/|_|  \__,_| |____/ \___/ \___|_|\_\  
================================================================================{Colors.RESET}
  {Colors.BOLD}RESEARCH-GRADE COMPUTATIONAL DRUG-RECEPTOR DOCKING PLATFORM{Colors.RESET}
  {Colors.DIM}Scripps AutoDock Vina 1.2.7 | RDKit | Meeko | CASF-2016 Core Verified{Colors.RESET}
"""

def print_banner():
    print(BANNER)

def print_section(title: str):
    print(f"\n{Colors.BOLD}{Colors.CYAN}═══ {title} ═══{Colors.RESET}\n")

def print_success(msg: str):
    print(f"{Colors.GREEN}✓ {msg}{Colors.RESET}")

def print_info(msg: str):
    print(f"{Colors.CYAN}ℹ {msg}{Colors.RESET}")

def print_warning(msg: str):
    print(f"{Colors.YELLOW}[!] {msg}{Colors.RESET}")

def print_error(msg: str):
    print(f"{Colors.RED}✖ {msg}{Colors.RESET}")


# ---------------------------------------------------------------------------
# Cross-Platform Arrow Key Input Reader
# ---------------------------------------------------------------------------
def read_key():
    """Reads a single keypress, returning 'UP', 'DOWN', 'ENTER', 'QUIT', digits, etc."""
    if os.name == 'nt':
        import msvcrt
        ch = msvcrt.getch()
        if ch in (b'\x00', b'\xe0'):
            ext = msvcrt.getch()
            if ext == b'H': return 'UP'
            if ext == b'P': return 'DOWN'
            if ext == b'K': return 'LEFT'
            if ext == b'M': return 'RIGHT'
        elif ch in (b'\r', b'\n'):
            return 'ENTER'
        elif ch in (b'\x1b', b'q', b'Q'):
            return 'QUIT'
        elif ch == b'\x03':
            raise KeyboardInterrupt
        elif ch in (b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9'):
            return ch.decode('ascii')
        return None
    else:
        import select, termios, tty
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)
            if ch == '\x1b':
                r, _, _ = select.select([sys.stdin], [], [], 0.08)
                if r:
                    ch2 = sys.stdin.read(1)
                    if ch2 == '[':
                        ch3 = sys.stdin.read(1)
                        if ch3 == 'A': return 'UP'
                        if ch3 == 'B': return 'DOWN'
                        if ch3 == 'D': return 'LEFT'
                        if ch3 == 'C': return 'RIGHT'
                return 'QUIT'
            elif ch in ('\r', '\n'):
                return 'ENTER'
            elif ch in ('q', 'Q'):
                return 'QUIT'
            elif ch in '123456789':
                return ch
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return None


# ---------------------------------------------------------------------------
# Interactive Arrow-Key Menu (TUI)
# ---------------------------------------------------------------------------
def interactive_menu(title: str, options: list, default_idx: int = 0) -> int:
    """
    Renders a live interactive menu navigated via Up/Down arrow keys.
    Pressing Enter selects the option.
    Pressing numbers (1-9) jumps directly.
    Pressing Q returns -1.
    """
    idx = default_idx
    n = len(options)

    # If terminal is non-interactive, fall back to standard text input
    if not sys.stdin.isatty():
        print(f"\n{title}")
        for i, opt in enumerate(options, 1):
            print(f"  {i}. {opt}")
        choice = input(f"Enter choice (1-{n}) [default: 1]: ").strip() or "1"
        try:
            return max(0, min(n - 1, int(choice) - 1))
        except ValueError:
            return 0

    while True:
        os.system('cls' if os.name == 'nt' else 'clear')
        print_banner()
        print(f"  {Colors.BOLD}{Colors.YELLOW}▶ {title}{Colors.RESET}")
        print(f"  {Colors.DIM}[ Use ↑ / ↓ Arrow Keys to Move | Enter to Select | Q to Back/Quit ]{Colors.RESET}\n")

        for i, opt in enumerate(options):
            if i == idx:
                # Active selection line
                print(f"  {Colors.BG_CYAN}{Colors.BOLD} ❯  {opt} {Colors.RESET}")
            else:
                print(f"     {Colors.DIM}{opt}{Colors.RESET}")

        print(f"\n  {Colors.DIM}────────────────────────────────────────────────────────────────────────{Colors.RESET}")

        try:
            key = read_key()
        except KeyboardInterrupt:
            print("\nExiting...")
            sys.exit(0)

        if key == 'UP':
            idx = (idx - 1) % n
        elif key == 'DOWN':
            idx = (idx + 1) % n
        elif key == 'ENTER':
            return idx
        elif key == 'QUIT':
            return -1
        elif key and key.isdigit():
            k_num = int(key) - 1
            if 0 <= k_num < n:
                return k_num


# ---------------------------------------------------------------------------
# CLI Command: Custom Docking
# ---------------------------------------------------------------------------
def cmd_dock(args):
    from backend.services.fetcher import StructureFetcher
    from backend.services.docking import DockingEngine
    from backend.services.bioactivity import BioactivityService

    print_section("STEP 1: Receptor Ingestion & Pocket Preparation")
    receptor_input = args.receptor
    pdb_content = None
    pdb_id = "CUSTOM"

    if os.path.exists(receptor_input):
        print_info(f"Loading local receptor PDB file: {receptor_input}")
        with open(receptor_input, "r", encoding="utf-8", errors="replace") as f:
            pdb_content = f.read()
        pdb_id = Path(receptor_input).stem.upper()
    elif len(receptor_input) == 4 and receptor_input.isalnum():
        print_info(f"Fetching crystallographic structure from RCSB PDB: {receptor_input.upper()}...")
        rec_meta = StructureFetcher.fetch_rcsb_pdb(receptor_input.upper())
        if not rec_meta or "pdb_content" not in rec_meta:
            print_error(f"Failed to fetch PDB ID '{receptor_input}' from RCSB.")
            return False
        pdb_content = rec_meta["pdb_content"]
        pdb_id = receptor_input.upper()
        print_success(f"Fetched PDB '{pdb_id}': {rec_meta.get('title', '')[:65]}...")
    else:
        print_error(f"Invalid receptor input '{receptor_input}'. Provide a valid 4-letter PDB ID or local .pdb path.")
        return False

    print_info("Processing receptor: stripping crystallographic waters, assigning physiological pH 7.4...")
    try:
        rec_prep = DockingEngine.prepare_receptor(pdb_content, target_chain=args.chain)
    except Exception as e:
        print_error(f"Receptor preparation failed: {e}")
        return False

    print_success(f"Receptor prepared ({rec_prep.get('atom_count', 0)} clean heavy atoms)")

    # Pocket Selection
    if args.center and args.size:
        pocket = {
            "center": {"x": args.center[0], "y": args.center[1], "z": args.center[2]},
            "size": {"x": args.size[0], "y": args.size[1], "z": args.size[2]},
            "description": "User-defined custom bounding box"
        }
    elif rec_prep.get("detected_pocket"):
        pocket = rec_prep["detected_pocket"]
        print_success(f"Pocket: {pocket.get('description', 'Auto-detected binding pocket')}")
    else:
        pocket = rec_prep.get("blind_docking_box")
        print_warning("No co-crystallized ligand found. Using Blind Docking search box.")

    c = pocket["center"]
    s = pocket["size"]
    print_info(f"Search Space Grid: Center=({c['x']:.1f}, {c['y']:.1f}, {c['z']:.1f}), Size=({s['x']:.1f}, {s['y']:.1f}, {s['z']:.1f}) Å")

    print_section("STEP 2: Ligand Chemical Ingestion & Conformation")
    lig_input = args.ligand
    smiles = None
    lig_name = "CustomLigand"

    if os.path.exists(lig_input):
        print_info(f"Loading local ligand file: {lig_input}")
        with open(lig_input, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
        is_sdf = lig_input.lower().endswith(".sdf")
        lig_prep = DockingEngine.prepare_ligand(content, is_sdf=is_sdf)
        smiles = lig_prep.get("canonical_smiles")
        lig_name = Path(lig_input).stem
    elif "C" in lig_input or "=" in lig_input or "#" in lig_input:
        smiles = lig_input
        print_info(f"Processing SMILES string: {smiles}")
        lig_prep = DockingEngine.prepare_ligand(smiles)
    else:
        print_info(f"Searching PubChem for compound name: '{lig_input}'...")
        pub_res = StructureFetcher.search_pubchem(lig_input)
        if not pub_res or "smiles" not in pub_res:
            print_error(f"Could not find compound '{lig_input}' on PubChem.")
            return False
        smiles = pub_res["smiles"]
        lig_name = pub_res.get("name", lig_input)
        print_success(f"Found on PubChem: {lig_name} (CID: {pub_res.get('cid')}, Formula: {pub_res.get('formula')})")
        lig_prep = DockingEngine.prepare_ligand(smiles)

    print_success(f"Ligand 3D conformation minimized with MMFF94 forcefield ({lig_prep.get('heavy_atom_count', 0)} heavy atoms, {lig_prep.get('rotatable_bonds', 0)} rotatable bonds)")

    print_section("STEP 3: Scripps AutoDock Vina v1.2.7 Docking Execution")
    from backend.services.hardware_profiler import HardwareProfiler
    hw = HardwareProfiler.get_hardware_profile()
    power_mode = getattr(args, "power_mode", "smart")

    exh, exh_reason = HardwareProfiler.calculate_adaptive_exhaustiveness(
        requested_exhaustiveness=args.exhaustiveness,
        rotatable_bonds=lig_prep.get("rotatable_bonds", 0),
        heavy_atoms=lig_prep.get("heavy_atom_count", 0),
        power_mode=power_mode
    )
    modes = args.modes
    print_info(f"Hardware Profile: {hw['icon']} {hw['tier_name']} ({hw['badge']})")
    print_info(f"Conformational Strategy: {exh_reason}")
    print_info(f"Starting Monte Carlo search (Exhaustiveness = {exh}, Modes = {modes}, Seed = {args.seed}, Power = {power_mode})...")

    t_dock0 = time.time()
    try:
        poses = DockingEngine.run_docking(
            receptor_pdbqt=rec_prep["pdbqt_text"],
            ligand_pdbqt=lig_prep["pdbqt_text"],
            center=pocket["center"],
            size=pocket["size"],
            exhaustiveness=exh,
            num_modes=modes,
            seed=args.seed,
            cpu=getattr(args, "cpu", None),
            power_mode=power_mode
        )
    except Exception as e:
        print_error(f"AutoDock Vina execution failed: {e}")
        return False

    dock_time = round(time.time() - t_dock0, 2)
    if not poses:
        print_error("AutoDock Vina completed but returned no binding poses.")
        return False

    print_success(f"Docking calculation completed in {dock_time}s! Produced {len(poses)} poses.")

    # Output directory
    out_dir = Path(args.out or f"docking_{pdb_id}_{lig_name}")
    out_dir.mkdir(parents=True, exist_ok=True)

    top_pose = poses[0]
    affinity = top_pose["affinity_kcal"]
    heavy_atoms = lig_prep.get("heavy_atom_count", 15)
    mw = lig_prep.get("molecular_weight", 300.0)

    thermo = BioactivityService.calculate_thermodynamics(affinity, heavy_atoms, mw)

    print_section("STEP 4: Binding Energetics & Thermodynamics")
    print(f"  {Colors.BOLD}Top Pose Binding Affinity (ΔG):{Colors.RESET} {Colors.GREEN}{affinity:.2f} kcal/mol{Colors.RESET}")
    print(f"  {Colors.BOLD}Predicted Dissociation Constant (Kd):{Colors.RESET} {thermo.get('theoretical_kd_nm', 'N/A')} nM ({thermo.get('theoretical_kd_um', 'N/A')} µM)")
    print(f"  {Colors.BOLD}pKd (-log10 Kd):{Colors.RESET} {thermo.get('pkd', 'N/A')}")
    print(f"  {Colors.BOLD}Ligand Efficiency (LE):{Colors.RESET} {thermo.get('ligand_efficiency', {}).get('value')} kcal/mol/heavy atom")
    print(f"  {Colors.BOLD}Potency Tier:{Colors.RESET} {thermo.get('potency_class', 'N/A')}")

    # Analyze Contacts
    contacts = DockingEngine.analyze_interactions(rec_prep["cleaned_pdb"], top_pose["pdbqt_content"])
    print_section("STEP 5: Intermolecular Active-Site Contacts")
    print(f"  Hydrogen Bonds:     {contacts.get('total_hbond_count', 0)}")
    print(f"  Hydrophobic Bonds:  {contacts.get('total_hydrophobic_count', 0)}")
    print(f"  Key Residues:       {', '.join(contacts.get('interacting_residues', [])[:8])}")

    # Save output files
    pose_pdbqt_file = out_dir / f"{lig_name}_top_pose.pdbqt"
    pose_pdbqt_file.write_text(top_pose["pdbqt_content"], encoding="utf-8")
    
    receptor_clean_file = out_dir / f"{pdb_id}_clean_receptor.pdb"
    receptor_clean_file.write_text(rec_prep["cleaned_pdb"], encoding="utf-8")

    # Summary JSON report
    report_data = {
        "pdb_id": pdb_id,
        "ligand_name": lig_name,
        "smiles": smiles,
        "binding_affinity_kcal": affinity,
        "exhaustiveness": exh,
        "docking_time_seconds": dock_time,
        "thermodynamics": thermo,
        "interactions": contacts,
        "poses": [{"mode": p["mode"], "affinity_kcal": p["affinity_kcal"]} for p in poses]
    }
    report_file = out_dir / "docking_summary.json"
    with open(report_file, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=2)

    print_section("STEP 6: Output Files & 3D Visualization Instructions")
    print_success(f"Top Docked Pose:       {pose_pdbqt_file}")
    print_success(f"Cleaned Receptor:      {receptor_clean_file}")
    print_success(f"Full Scientific JSON:  {report_file}")
    
    print(f"""
{Colors.CYAN}{Colors.BOLD}HOW TO VISUALIZE YOUR 3D DOCKED COMPLEX:{Colors.RESET}
  1. {Colors.BOLD}In PyMOL:{Colors.RESET}
     Open terminal and type:
     {Colors.YELLOW}pymol {receptor_clean_file} {pose_pdbqt_file}{Colors.RESET}

  2. {Colors.BOLD}In UCSF Chimera / ChimeraX:{Colors.RESET}
     {Colors.YELLOW}chimerax {receptor_clean_file} {pose_pdbqt_file}{Colors.RESET}

  3. {Colors.BOLD}In Bindora Web Studio (Browser 3D Viewer):{Colors.RESET}
     Run {Colors.YELLOW}python backend/app.py{Colors.RESET} and open {Colors.YELLOW}http://localhost:5000{Colors.RESET}
""")
    return True


# ---------------------------------------------------------------------------
# CLI Command: ADME Profiler
# ---------------------------------------------------------------------------
def cmd_adme(args):
    from backend.services.adme import ADMEProfiler

    print_section(f"ADME & Pharmacokinetics Profiler: {args.smiles}")
    res = ADMEProfiler.calculate_adme(args.smiles)
    if "error" in res:
        print_error(f"ADME calculation failed: {res['error']}")
        return False

    pc = res.get("physicochemical", {})
    dl = res.get("drug_likeness", {})
    lip = dl.get("lipinski", {})
    pk = res.get("pharmacokinetics", {})
    medchem = res.get("medicinal_chemistry", {})

    mw_val = pc.get('molecular_weight', {}).get('value', 'N/A')
    logp_val = pc.get('logp', {}).get('value', 'N/A')
    hbd_val = pc.get('hbd', {}).get('value', 'N/A')
    hba_val = pc.get('hba', {}).get('value', 'N/A')
    rotb_val = pc.get('rotatable_bonds', {}).get('value', 'N/A')
    tpsa_val = pc.get('tpsa', {}).get('value', 'N/A')
    viol_count = len(lip.get('violations', [])) if isinstance(lip.get('violations'), list) else 0

    print(f"  {Colors.BOLD}Molecular Weight:{Colors.RESET}      {mw_val} Da")
    print(f"  {Colors.BOLD}LogP (Lipophilicity):{Colors.RESET}  {logp_val}")
    print(f"  {Colors.BOLD}H-Bond Donors:{Colors.RESET}         {hbd_val}")
    print(f"  {Colors.BOLD}H-Bond Acceptors:{Colors.RESET}      {hba_val}")
    print(f"  {Colors.BOLD}Rotatable Bonds:{Colors.RESET}       {rotb_val}")
    print(f"  {Colors.BOLD}TPSA:{Colors.RESET}                  {tpsa_val} Å²")
    status_col = Colors.GREEN if lip.get('status') == 'Pass' else Colors.RED
    gi_val = pk.get('gi_absorption', {}).get('level', 'N/A') if isinstance(pk.get('gi_absorption'), dict) else pk.get('gi_absorption', 'N/A')
    bbb_val = pk.get('bbb_permeation', {}).get('status', 'N/A') if isinstance(pk.get('bbb_permeation'), dict) else pk.get('bbb_permeant', 'N/A')
    ppb_val = pk.get('plasma_protein_binding', {}).get('tier', 'N/A') if isinstance(pk.get('plasma_protein_binding'), dict) else 'N/A'

    print(f"  {Colors.BOLD}Lipinski Rule of 5:{Colors.RESET}    {status_col}{lip.get('status', 'N/A')} ({viol_count} violations){Colors.RESET}")
    print(f"  {Colors.BOLD}GI Absorption:{Colors.RESET}         {gi_val}")
    print(f"  {Colors.BOLD}BBB Permeation:{Colors.RESET}        {bbb_val}")
    print(f"  {Colors.BOLD}Plasma Protein Binding:{Colors.RESET} {ppb_val}")
    print(f"  {Colors.BOLD}PAINS Filter:{Colors.RESET}          {medchem.get('pains_alerts', {}).get('count', 0)} alerts")
    return True


# ---------------------------------------------------------------------------
# Interactive Arrow-Key Guided Wizard
# ---------------------------------------------------------------------------
def run_interactive_wizard():
    while True:
        main_options = [
            "Run Custom Molecular Docking (Receptor PDB + Ligand)",
            "Calculate ADME & Pharmacokinetics Profile (from SMILES)",
            "Run CASF-2016 Benchmark Suite (Field-Standard 285 complexes)",
            "Run DUD-E Virtual Screening Benchmark (ROC-AUC / Enrichment)",
            "Launch Web Studio UI (localhost:5000 in Browser)",
            "Help & Command Line Cheatsheet",
            "Exit"
        ]

        choice = interactive_menu("MAIN MENU — Select Action", main_options)

        if choice == 0:
            # Submenu: Docking setup
            dock_options = [
                "Quick Demo: Cyclooxygenase-2 (1CX2) vs Aspirin",
                "Quick Demo: Abl1 Kinase (1IEP) vs Imatinib",
                "Quick Demo: CDK2 Kinase (1AQ1) vs Staurosporine",
                "Custom RCSB PDB ID + Drug Name (Auto-downloaded)",
                "Local Protein .PDB File + Chemical SMILES",
                "Local Protein .PDB File + Multi-conformer .SDF File",
                "<-- Back to Main Menu"
            ]

            dock_choice = interactive_menu("DOCKING WIZARD — Choose Input Type", dock_options)

            if dock_choice == 6 or dock_choice == -1:
                continue

            rec = "1CX2"
            lig = "aspirin"
            if dock_choice == 0:
                rec, lig = "1CX2", "aspirin"
            elif dock_choice == 1:
                rec, lig = "1IEP", "imatinib"
            elif dock_choice == 2:
                rec, lig = "1AQ1", "staurosporine"
            elif dock_choice == 3:
                os.system('cls' if os.name == 'nt' else 'clear')
                print_banner()
                print_section("Custom RCSB PDB & PubChem Ingestion")
                rec = input(f"{Colors.BOLD}Enter 4-letter RCSB PDB ID (e.g. 1CX2, 1IEP): {Colors.RESET}").strip() or "1CX2"
                lig = input(f"{Colors.BOLD}Enter Drug Name or SMILES (e.g. aspirin, gefitinib): {Colors.RESET}").strip() or "aspirin"
            elif dock_choice == 4:
                os.system('cls' if os.name == 'nt' else 'clear')
                print_banner()
                print_section("Local PDB File & SMILES")
                rec = input(f"{Colors.BOLD}Enter path to local .pdb file (e.g. receptor.pdb): {Colors.RESET}").strip()
                lig = input(f"{Colors.BOLD}Enter SMILES string (e.g. CC(=O)Oc1ccccc1C(=O)O): {Colors.RESET}").strip()
            elif dock_choice == 5:
                os.system('cls' if os.name == 'nt' else 'clear')
                print_banner()
                print_section("Local PDB & SDF File")
                rec = input(f"{Colors.BOLD}Enter path to local .pdb file: {Colors.RESET}").strip()
                lig = input(f"{Colors.BOLD}Enter path to local .sdf file: {Colors.RESET}").strip()

            # Submenu: Exhaustiveness Selection
            exh_options = [
                "●  Standard Research (Exhaustiveness 8) — Recommended (~10-15s)",
                "●  Rapid Screening (Exhaustiveness 4)   — Fast (~3-5s)",
                "●  Publication Grade (Exhaustiveness 16) — Deep Sampling (~30s)",
                "●  Maximum Thoroughness (Exhaustiveness 32) — Exhaustive (~2m)"
            ]
            exh_choice = interactive_menu("SEARCH DEPTH — Select Vina Exhaustiveness", exh_options, default_idx=0)
            exh_map = [8, 4, 16, 32]
            exh_val = exh_map[exh_choice] if exh_choice != -1 else 8

            os.system('cls' if os.name == 'nt' else 'clear')
            print_banner()

            class DockArgs:
                receptor = rec
                ligand = lig
                exhaustiveness = exh_val
                modes = 9
                seed = 42
                chain = None
                center = None
                size = None
                out = None

            cmd_dock(DockArgs)
            input(f"\n{Colors.YELLOW}Press Enter to return to main menu...{Colors.RESET}")

        elif choice == 1:
            os.system('cls' if os.name == 'nt' else 'clear')
            print_banner()
            print_section("ADME & Pharmacokinetics Analysis")
            smi = input(f"{Colors.BOLD}Enter SMILES string [default: Aspirin]: {Colors.RESET}").strip() or "CC(=O)Oc1ccccc1C(=O)O"
            class AdmeArgs:
                smiles = smi
            cmd_adme(AdmeArgs)
            input(f"\n{Colors.YELLOW}Press Enter to return to main menu...{Colors.RESET}")

        elif choice == 2:
            bm_options = [
                "●  Quick Smoke Test (5 Complexes, ~5 min)",
                "●  Medium Validation (20 Complexes, ~20 min)",
                "●  Full CASF-2016 Core Set (285 Complexes, Overnight / Resumable)",
                "<-- Back to Main Menu"
            ]
            bm_choice = interactive_menu("CASF-2016 BENCHMARK — Select Evaluation Scale", bm_options)
            if bm_choice == 3 or bm_choice == -1:
                continue
            nums = [5, 20, 285]
            num_val = nums[bm_choice]

            os.system('cls' if os.name == 'nt' else 'clear')
            print_banner()
            print_section(f"Starting CASF-2016 Benchmark on {num_val} Complexes")
            cmd = f"python tests/benchmark_casf2016.py --num {num_val} --exhaustiveness 4 --resume"
            print_info(f"Executing: {cmd}")
            os.system(cmd)
            input(f"\n{Colors.YELLOW}Press Enter to return to main menu...{Colors.RESET}")

        elif choice == 3:
            dud_options = [
                "●  Acetylcholinesterase (ache)",
                "●  HIV-1 Protease (hivpr)",
                "●  Proto-oncogene Tyrosine Kinase Src (src)",
                "●  Vascular Endothelial Growth Factor Receptor 2 (vegfr2)",
                "●  Custom Target Name",
                "<-- Back to Main Menu"
            ]
            dud_choice = interactive_menu("DUD-E VIRTUAL SCREENING — Select Target", dud_options)
            if dud_choice == 5 or dud_choice == -1:
                continue
            dud_targets = ["ache", "hivpr", "src", "vegfr2"]
            if dud_choice < 4:
                tgt = dud_targets[dud_choice]
            else:
                tgt = input(f"{Colors.BOLD}Enter DUD-E target name (e.g. cdk2, egfr): {Colors.RESET}").strip() or "ache"

            os.system('cls' if os.name == 'nt' else 'clear')
            print_banner()
            print_section(f"Starting DUD-E Screening for Target: {tgt}")
            cmd = f"python tests/benchmark_screening.py --target {tgt} --exhaustiveness 4 --max-actives 20 --max-decoys 100 --resume"
            print_info(f"Executing: {cmd}")
            os.system(cmd)
            input(f"\n{Colors.YELLOW}Press Enter to return to main menu...{Colors.RESET}")

        elif choice == 4:
            import urllib.request
            import webbrowser
            import time
            os.system('cls' if os.name == 'nt' else 'clear')
            print_banner()
            print_section("Bindora Web Studio UI (3D Mol* Viewer & Interactive Analytics)")

            is_running = False
            try:
                with urllib.request.urlopen("http://127.0.0.1:5000/api/health", timeout=1) as response:
                    if response.status == 200:
                        is_running = True
            except Exception:
                is_running = False

            if is_running:
                print_success("Bindora Web Studio server is already running on http://localhost:5000")
            else:
                print_info("Starting Bindora Web Studio server in the background on http://localhost:5000 ...")
                if os.name == 'nt':
                    subprocess.Popen(
                        [sys.executable, "backend/app.py"],
                        creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, 'CREATE_NO_WINDOW') else 0,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL
                    )
                else:
                    subprocess.Popen(
                        [sys.executable, "backend/app.py"],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                        start_new_session=True
                    )
                # Wait up to 3 seconds for port 5000 to become responsive
                for _ in range(6):
                    time.sleep(0.5)
                    try:
                        with urllib.request.urlopen("http://127.0.0.1:5000/api/health", timeout=1) as response:
                            if response.status == 200:
                                is_running = True
                                break
                    except Exception:
                        pass

            print_info("Opening http://localhost:5000 in your web browser...")
            webbrowser.open("http://localhost:5000")
            print_success("Web Studio launched! You can interact with 3D poses and charts in your browser.")
            input(f"\n{Colors.YELLOW}Press Enter to return to CLI main menu...{Colors.RESET}")

        elif choice == 5:
            os.system('cls' if os.name == 'nt' else 'clear')
            print_banner()
            print_section("Command-Line Cheat Sheet & Scripting Reference")
            print(f"""
{Colors.BOLD}1. Direct Molecular Docking:{Colors.RESET}
   python bindora_cli.py dock --receptor 1CX2 --ligand aspirin --exhaustiveness 8
   python bindora_cli.py dock --receptor my_prot.pdb --ligand "CC(=O)Oc1ccccc1C(=O)O" --out my_results
   python bindora_cli.py dock --receptor my_prot.pdb --ligand candidate.sdf --exhaustiveness 16

{Colors.BOLD}2. ADME Pharmacokinetics:{Colors.RESET}
   python bindora_cli.py adme --smiles "CC(=O)Oc1ccccc1C(=O)O"

{Colors.BOLD}3. CASF-2016 Core Set Benchmark:{Colors.RESET}
   python bindora_cli.py benchmark --type casf2016 --num 10 --exhaustiveness 4

{Colors.BOLD}4. DUD-E Virtual Screening:{Colors.RESET}
   python bindora_cli.py screen --target ache --max-actives 20 --max-decoys 100

{Colors.BOLD}5. 3D Visualizer Commands:{Colors.RESET}
   PyMOL:     pymol my_results/protein_clean.pdb my_results/ligand_top_pose.pdbqt
   ChimeraX:  chimerax my_results/protein_clean.pdb my_results/ligand_top_pose.pdbqt
""")
            input(f"\n{Colors.YELLOW}Press Enter to return to main menu...{Colors.RESET}")

        elif choice == 6 or choice == -1:
            os.system('cls' if os.name == 'nt' else 'clear')
            print(f"\n{Colors.CYAN}Thank you for using Bindora Dock v2.0! Happy Research.{Colors.RESET}\n")
            sys.exit(0)


# ---------------------------------------------------------------------------
# Main CLI Parser (Preserves Direct Scripting Capability)
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) == 1:
        run_interactive_wizard()
        return

    parser = argparse.ArgumentParser(
        description="Bindora Dock v2.0 — Research-Grade Molecular Docking CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive arrow-key guided mode
  python bindora_cli.py

  # Dock a drug against a crystal structure from RCSB
  python bindora_cli.py dock --receptor 1CX2 --ligand aspirin --exhaustiveness 8

  # Dock a local PDB against a SMILES string and save to folder
  python bindora_cli.py dock --receptor receptor.pdb --ligand "CC(=O)Oc1ccccc1C(=O)O" --out my_results

  # Calculate ADME properties
  python bindora_cli.py adme --smiles "CC(=O)Oc1ccccc1C(=O)O"

  # Run CASF-2016 benchmark on top 10 complexes
  python bindora_cli.py benchmark --type casf2016 --num 10
        """
    )

    subparsers = parser.add_subparsers(dest="subcommand", help="Available commands")

    # Command: dock
    dock_parser = subparsers.add_parser("dock", help="Run molecular docking on a receptor and ligand")
    dock_parser.add_argument("--receptor", required=True, help="4-letter RCSB PDB ID (e.g. 1CX2) or path to local .pdb file")
    dock_parser.add_argument("--ligand", required=True, help="Drug name (e.g. aspirin), SMILES string, or path to .sdf file")
    dock_parser.add_argument("--chain", default=None, help="Receptor chain ID (optional, default: auto)")
    dock_parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina search exhaustiveness (default: 8)")
    dock_parser.add_argument("--modes", type=int, default=9, help="Number of binding poses to produce (default: 9)")
    dock_parser.add_argument("--seed", type=int, default=42, help="Random seed for full reproducibility (default: 42)")
    dock_parser.add_argument("--cpu", type=int, default=None, help="Number of CPU cores for AutoDock Vina (default: all available)")
    dock_parser.add_argument("--power-mode", choices=["smart", "performance", "eco"], default="smart", help="Adaptive compute mode (default: smart)")
    dock_parser.add_argument("--center", nargs=3, type=float, default=None, help="Search grid center: X Y Z (Å)")
    dock_parser.add_argument("--size", nargs=3, type=float, default=None, help="Search grid size: X Y Z (Å)")
    dock_parser.add_argument("--out", default=None, help="Directory to save docked pose PDBQT, cleaned PDB, and JSON report")

    # Command: adme
    adme_parser = subparsers.add_parser("adme", help="Calculate RDKit ADME & physicochemical properties")
    adme_parser.add_argument("--smiles", required=True, help="SMILES string of compound")

    # Command: benchmark
    bm_parser = subparsers.add_parser("benchmark", help="Run standardized docking accuracy benchmarks")
    bm_parser.add_argument("--type", choices=["casf2016", "accuracy"], default="casf2016", help="Benchmark suite")
    bm_parser.add_argument("--num", type=int, default=10, help="Number of complexes to evaluate")
    bm_parser.add_argument("--exhaustiveness", type=int, default=4, help="Exhaustiveness setting")
    bm_parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")

    # Command: screen
    screen_parser = subparsers.add_parser("screen", help="Run DUD-E virtual screening enrichment evaluation")
    screen_parser.add_argument("--target", required=True, help="DUD-E target name (e.g. ache, hivpr, src)")
    screen_parser.add_argument("--exhaustiveness", type=int, default=4, help="Exhaustiveness setting")
    screen_parser.add_argument("--max-actives", type=int, default=None, help="Limit number of active molecules")
    screen_parser.add_argument("--max-decoys", type=int, default=None, help="Limit number of decoys")
    screen_parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")

    # Command: interactive
    subparsers.add_parser("interactive", help="Start the interactive arrow-key wizard")

    args = parser.parse_args()

    if args.subcommand == "dock":
        cmd_dock(args)
    elif args.subcommand == "adme":
        cmd_adme(args)
    elif args.subcommand == "benchmark":
        cmd_benchmark(args)
    elif args.subcommand == "screen":
        cmd_screen(args)
    elif args.subcommand == "interactive" or not args.subcommand:
        run_interactive_wizard()

if __name__ == "__main__":
    main()
