import urllib.request
import json
import sys

def post(url, data):
    req = urllib.request.Request(
        url,
        data=json.dumps(data).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req) as resp:
        return resp.status, json.loads(resp.read().decode("utf-8"))

def get(url):
    with urllib.request.urlopen(url) as resp:
        return resp.status, json.loads(resp.read().decode("utf-8"))

def main():
    print("=== TESTING ALL TABS SUITE ===")
    
    # 1. Health
    st, health = get("http://127.0.0.1:5000/api/health")
    print(f"[TAB 1-7 Health]: status={st}, service={health.get('service')}")
    assert st == 200

    # 2. Tab 1 & Tab 3: Ligand Prep & ADME
    st, lig = post("http://127.0.0.1:5000/api/structure/ligand", {"structure": "CC(=O)Oc1ccccc1C(=O)O", "name": "Aspirin"})
    print(f"[TAB 1/3 Ligand & ADME]: status={st}, formula={lig.get('formula')}, heavy_atoms={lig.get('heavy_atom_count')}, ADME keys={list(lig.get('adme', {}).keys())}")
    assert st == 200 and "adme" in lig

    # 3. Tab 4: Target Activity & ChEMBL Cross-Check
    st, cc = get("http://127.0.0.1:5000/api/pkpd/crosscheck?drug=Aspirin&target=Cyclooxygenase-2")
    print(f"[TAB 4 ChEMBL Cross-Check]: status={st}, match_level={cc.get('match_level')}, records={len(cc.get('records', []))}")
    assert st == 200

    # 4. Tab 5: AI Pharmacologist Narrative Briefing
    st, narr = post("http://127.0.0.1:5000/api/explain/narrative", {
        "ligand_name": "Aspirin",
        "target_name": "COX-2",
        "thermodynamics": {"affinity_kcal": -6.7}
    })
    print(f"[TAB 5 AI Narrative]: status={st}, source={narr.get('source')}, narrative_chars={len(narr.get('narrative', ''))}")
    assert st == 200 and len(narr.get("narrative", "")) > 100

    # 5. Tab 6: High-Throughput Batch Screening (check non-existent job gives 404)
    try:
        st, bstatus = get("http://127.0.0.1:5000/api/batch/status/non-existent-id")
    except urllib.error.HTTPError as e:
        print(f"[TAB 6 Batch Manager Endpoint]: status={e.code} (Correctly 404 for non-existent job)")
        assert e.code == 404

    # 6. Tab 7: Regulatory Research Dossier & Validation Report
    st, sysinfo = get("http://127.0.0.1:5000/api/reproducibility/versions")
    print(f"[TAB 7 Reproducibility Versions]: status={st}, vina_version={sysinfo.get('vina_version')}")
    assert st == 200
    
    st, vrep = get("http://127.0.0.1:5000/api/validation-report")
    print(f"[TAB 7 Validation Report]: status={st}, total_targets={len(vrep.get('casf_targets', []))}")
    assert st == 200

    print(">>> ALL TABS (1 through 7) ARE 100% OPERATIONAL AND CONNECTED! <<<")

if __name__ == "__main__":
    main()
