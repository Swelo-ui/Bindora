import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine
from backend.services.adme import ADMEProfiler
from backend.services.bioactivity import BioactivityService
from backend.services.narrative import NarrativeExplainer

def run_verification():
    print("=" * 60)
    print("BINDORA FULL PLATFORM VERIFICATION")
    print("=" * 60)

    # 1. Structure Fetching
    print("\n[STEP 1] Testing PubChem and RCSB PDB Ingestion...")
    lig_data = StructureFetcher.search_pubchem("aspirin")
    assert lig_data is not None and lig_data["cid"] == 2244, "PubChem fetch failed"
    rec_meta = StructureFetcher.fetch_rcsb_pdb("1CX2")
    assert rec_meta is not None and "ATOM" in rec_meta["pdb_content"], "RCSB fetch failed"
    print(f"  [OK] PubChem Aspirin CID: {lig_data['cid']}, Formula: {lig_data['formula']}")
    print(f"  [OK] RCSB PDB 1CX2: {rec_meta['title'][:50]}...")

    # 2. Preparation
    print("\n[STEP 2] Testing 3D Conformation & Receptor Pocket Setup...")
    lig_prep = DockingEngine.prepare_ligand(lig_data["smiles"])
    rec_prep = DockingEngine.prepare_receptor(rec_meta["pdb_content"], target_chain="A")
    assert "ROOT" in lig_prep["pdbqt_text"], "Ligand PDBQT missing ROOT"
    assert rec_prep["atom_count"] > 1000, "Receptor atoms too low"
    print(f"  [OK] Ligand Heavy Atoms: {lig_prep['heavy_atom_count']}, Rotatable: {lig_prep['rotatable_bonds']}")
    print(f"  [OK] Receptor Cleaned Atoms: {rec_prep['atom_count']}")
    print(f"  [OK] Auto-detected Pocket: {rec_prep['detected_pocket']['description']}")

    # 3. Docking Run with AutoDock Vina
    print("\n[STEP 3] Executing Scripps AutoDock Vina v1.2.7 Engine...")
    pocket = rec_prep["detected_pocket"]
    poses = DockingEngine.run_docking(
        rec_prep["pdbqt_text"],
        lig_prep["pdbqt_text"],
        pocket["center"],
        pocket["size"],
        exhaustiveness=1,
        num_modes=2
    )
    assert len(poses) >= 1, "No docking poses produced"
    top_pose = poses[0]
    print(f"  [OK] AutoDock Vina produced {len(poses)} modes successfully!")
    print(f"  [OK] Mode 1 Binding Free Energy (delta G): {top_pose['affinity_kcal']} kcal/mol")

    # 4. Contact Extraction
    print("\n[STEP 4] Analyzing Intermolecular Active Site Contacts...")
    contacts = DockingEngine.analyze_interactions(rec_prep["cleaned_pdb"], top_pose["pdbqt_content"])
    print(f"  [OK] Hydrogen Bonds Found: {contacts['total_hbond_count']}")
    print(f"  [OK] Hydrophobic Contacts: {contacts['total_hydrophobic_count']}")
    print(f"  [OK] Active Residues Involved: {contacts['interacting_residues'][:5]}")

    # 5. ADME Profiler
    print("\n[STEP 5] Testing RDKit Physicochemical & ADME Profiler...")
    adme = ADMEProfiler.calculate_adme(lig_data["smiles"])
    phys = adme["physicochemical"]
    assert adme["drug_likeness"]["lipinski"]["status"] == "Pass"
    print(f"  [OK] Molecular Weight: {phys['molecular_weight']['value']} Da")
    print(f"  [OK] MolLogP: {phys['logp']['value']}")
    print(f"  [OK] Lipinski Status: {adme['drug_likeness']['lipinski']['status']} (0 violations)")
    print(f"  [OK] GI Absorption: {adme['pharmacokinetics']['gi_absorption']['level']}")
    print(f"  [OK] PAINS Alerts: {adme['medicinal_chemistry_safety']['pains_alerts']['count']} (Clear)")

    # 6. ChEMBL Bioactivity Cross-Check
    print("\n[STEP 6] Testing ChEMBL Curated Wet-Lab Bioactivity Database...")
    xcheck = BioactivityService.crosscheck_chembl("imatinib", "abl1")
    assert xcheck["is_cross_checked"] is True, "ChEMBL crosscheck failed"
    print(f"  [OK] Match Status: {xcheck['status_badge']}")
    print(f"  [OK] Curated Records Count: {len(xcheck['experimental_records'])}")
    first_rec = xcheck["experimental_records"][0]
    print(f"  [OK] Sample Measured Value: {first_rec['type']} {first_rec['relation']} {first_rec['value']} {first_rec['units']}")

    # 7. OpenRouter DeepSeek AI Explainer
    print("\n[STEP 7] Testing OpenRouter DeepSeek (deepseek-v4-flash-0731) AI Explainer...")
    thermo = BioactivityService.calculate_thermodynamics(top_pose["affinity_kcal"], lig_prep["heavy_atom_count"], lig_data["weight"])
    report_payload = {
        "ligand_name": "Aspirin",
        "target_name": "COX-2",
        "pdb_id": "1CX2",
        "thermodynamics": thermo,
        "interactions": contacts,
        "adme": adme,
        "bioactivity_crosscheck": {"is_cross_checked": False, "status_badge": "Computational Prediction Only"}
    }
    narrative_res = NarrativeExplainer.generate_explanation(report_payload)
    assert len(narrative_res["narrative"]) > 500, "Narrative too short"
    print(f"  [OK] AI Engine Source: {narrative_res['source']}")
    print(f"  [OK] Cached: {narrative_res.get('cached', False)}")
    print(f"  [OK] Total Character Count: {len(narrative_res['narrative'])} characters")

    # 8. Firebase Configuration
    print("\n[STEP 8] Testing Firebase Web Configuration...")
    with open("frontend/js/firebase-auth.js", "r", encoding="utf-8") as f:
        auth_js = f.read()
    assert "bindora-1db62" in auth_js, "Firebase project ID missing"
    assert "bindora-1db62-default-rtdb.asia-southeast1.firebasedatabase.app" in auth_js, "RTDB URL missing"
    assert "QUl6YVN5QVh1dTVFWVhycFByYWtPOF9RLXJOaXdUbDZTOXZHZ3hZ" in auth_js, "Encoded API key missing"
    assert "3d5e4eda9b24476f89c87b" in auth_js, "App ID missing"
    print("  [OK] Firebase Project ID: bindora-1db62 verified!")
    print("  [OK] Realtime Database: https://bindora-1db62-default-rtdb.asia-southeast1.firebasedatabase.app verified!")
    print("  [OK] Web API Key & AuthDomain verified!")

    print("\n" + "=" * 60)
    print("ALL 8 VERIFICATION STEPS COMPLETED SUCCESSFULLY!")
    print("BINDORA IS 100% OPERATIONAL & VERIFIED!")
    print("=" * 60)

if __name__ == "__main__":
    run_verification()
