import os
import sys
import json
import traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

from backend.config import (
    BASE_DIR, FRONTEND_DIR, BENCHMARKS_DIR, VINA_EXE, HOST, PORT, DEBUG
)
from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine
from backend.services.adme import ADMEProfiler
from backend.services.bioactivity import BioactivityService
from backend.services.narrative import NarrativeExplainer
from backend.services.batch import BatchScreeningService
from backend.utils.vina_setup import ensure_vina

app = Flask(__name__, static_folder=str(FRONTEND_DIR), static_url_path="")
CORS(app)

# Ensure Vina binary is ready on startup
try:
    ensure_vina()
except Exception as e:
    print(f"[STARTUP WARNING] Vina initialization warning: {e}", file=sys.stderr)

# ----------------- Static Frontend Routes -----------------

@app.route("/")
def index():
    return send_from_directory(FRONTEND_DIR, "index.html")

@app.route("/<path:path>")
def static_proxy(path):
    file_path = FRONTEND_DIR / path
    if file_path.exists():
        return send_from_directory(FRONTEND_DIR, path)
    return send_from_directory(FRONTEND_DIR, "index.html")

# ----------------- API Endpoints -----------------

@app.route("/api/health", methods=["GET"])
def health():
    vina_ok = VINA_EXE.exists() and VINA_EXE.stat().st_size > 100000
    return jsonify({
        "status": "healthy",
        "service": "Bindora 3D Drug-Receptor & PK/PD Analyzer",
        "vina_available": vina_ok,
        "vina_path": str(VINA_EXE) if vina_ok else None
    })

@app.route("/api/benchmarks", methods=["GET"])
def get_benchmarks():
    bm_file = BENCHMARKS_DIR / "benchmarks.json"
    if bm_file.exists():
        with open(bm_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        return jsonify({"benchmarks": data})
    return jsonify({"benchmarks": []})

@app.route("/api/search/pubchem", methods=["GET"])
def search_pubchem():
    query = request.args.get("query", "").strip()
    if not query:
        return jsonify({"error": "Query parameter 'query' is required"}), 400
    
    result = StructureFetcher.search_pubchem(query)
    if not result:
        return jsonify({"error": f"No compound found in PubChem matching '{query}'"}), 404
    
    # Compute ADME immediately for fast UI response
    smiles = result.get("smiles")
    adme_data = ADMEProfiler.calculate_adme(smiles) if smiles else {}
    result["adme"] = adme_data
    
    return jsonify(result)

@app.route("/api/search/pubchem/cid/<cid>", methods=["GET"])
def search_pubchem_cid(cid):
    result = StructureFetcher.fetch_pubchem_by_cid(cid)
    if not result:
        return jsonify({"error": f"No compound found in PubChem with CID '{cid}'"}), 404
    
    smiles = result.get("smiles")
    adme_data = ADMEProfiler.calculate_adme(smiles) if smiles else {}
    result["adme"] = adme_data
    return jsonify(result)

@app.route("/api/search/rcsb", methods=["GET"])
def search_rcsb():
    query = request.args.get("query", "").strip()
    if not query:
        return jsonify({"error": "Query parameter 'query' is required"}), 400
    
    # If 4 characters, treat as direct PDB ID
    if len(query) == 4 and query.isalnum():
        result = StructureFetcher.fetch_rcsb_pdb(query)
        if result:
            # Also fetch UniProt annotation if available
            uniprot_info = StructureFetcher.search_uniprot(result.get("title", query))
            result["uniprot"] = uniprot_info
            return jsonify({"direct": True, "entry": result})
        return jsonify({"error": f"PDB structure '{query}' not found"}), 404
    
    # Otherwise perform full text keyword search
    results = StructureFetcher.search_rcsb_keywords(query, max_results=6)
    return jsonify({"direct": False, "results": results})

@app.route("/api/structure/ligand", methods=["POST"])
def prepare_ligand():
    data = request.get_json() or {}
    smiles_or_sdf = data.get("structure", "").strip()
    is_sdf = data.get("is_sdf", False)
    
    if not smiles_or_sdf:
        return jsonify({"error": "Structure content or SMILES string is required"}), 400

    try:
        prep_result = DockingEngine.prepare_ligand(smiles_or_sdf, is_sdf=is_sdf)
        smiles = prep_result.get("canonical_smiles")
        adme = ADMEProfiler.calculate_adme(smiles)
        prep_result["adme"] = adme
        return jsonify(prep_result)
    except Exception as e:
        return jsonify({"error": f"Ligand preparation failed: {str(e)}"}), 400

@app.route("/api/structure/receptor", methods=["POST"])
def prepare_receptor():
    data = request.get_json() or {}
    pdb_content = data.get("pdb_content", "")
    pdb_id = data.get("pdb_id", "")
    target_chain = data.get("target_chain")

    if not pdb_content and pdb_id:
        fetched = StructureFetcher.fetch_rcsb_pdb(pdb_id)
        if fetched:
            pdb_content = fetched.get("pdb_content", "")

    if not pdb_content:
        return jsonify({"error": "Receptor PDB content or valid PDB ID is required"}), 400

    try:
        rec_result = DockingEngine.prepare_receptor(pdb_content, target_chain=target_chain)
        return jsonify(rec_result)
    except Exception as e:
        return jsonify({"error": f"Receptor preparation failed: {str(e)}"}), 400

@app.route("/api/docking/run", methods=["POST"])
def run_docking():
    data = request.get_json() or {}
    receptor_pdbqt = data.get("receptor_pdbqt")
    ligand_pdbqt = data.get("ligand_pdbqt")
    receptor_pdb = data.get("receptor_pdb")
    center = data.get("center")
    size = data.get("size")
    exhaustiveness = int(data.get("exhaustiveness", 8))
    num_modes = int(data.get("num_modes", 9))
    replicates = int(data.get("replicates", 1))
    heavy_atoms = int(data.get("heavy_atoms", 20))
    mw = float(data.get("molecular_weight", 300.0))

    if not receptor_pdbqt or not ligand_pdbqt or not center or not size:
        return jsonify({"error": "Missing required docking parameters (receptor, ligand, center, size)"}), 400

    try:
        # Run AutoDock Vina
        poses = DockingEngine.run_docking(
            receptor_pdbqt,
            ligand_pdbqt,
            center,
            size,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            replicates=replicates
        )

        if not poses:
            return jsonify({"error": "AutoDock Vina finished but returned no binding poses."}), 500

        best_pose = poses[0]
        affinity = best_pose["affinity_kcal"]
        replicate_stats = best_pose.get("replicate_stats")

        # Calculate thermodynamics and Ligand Efficiency
        thermo = BioactivityService.calculate_thermodynamics(affinity, heavy_atoms, mw)

        # Calculate atomic interactions (H-bonds, hydrophobic)
        contacts = {}
        if receptor_pdb:
            contacts = DockingEngine.analyze_interactions(receptor_pdb, best_pose["pdbqt_content"])

        return jsonify({
            "poses": poses,
            "top_pose": best_pose,
            "thermodynamics": thermo,
            "interactions": contacts,
            "replicate_stats": replicate_stats
        })
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Docking execution failed: {str(e)}"}), 500

@app.route("/api/docking/redock-validate", methods=["POST"])
def redock_validate():
    data = request.get_json() or {}
    receptor_pdbqt = data.get("receptor_pdbqt")
    native_ligand_pdb = data.get("native_ligand_pdb")
    center = data.get("center")
    size = data.get("size")
    exhaustiveness = int(data.get("exhaustiveness", 8))

    if not receptor_pdbqt or not native_ligand_pdb or not center or not size:
        return jsonify({"error": "Missing required parameters for redocking validation (receptor_pdbqt, native_ligand_pdb, center, size)"}), 400

    try:
        validation_result = DockingEngine.run_redocking_validation(
            receptor_pdbqt,
            native_ligand_pdb,
            center,
            size,
            exhaustiveness=exhaustiveness
        )
        return jsonify(validation_result)
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Redocking validation failed: {str(e)}"}), 500

@app.route("/api/reproducibility/versions", methods=["GET"])
def get_reproducibility_versions():
    import gemmi
    import rdkit
    import datetime
    return jsonify({
        "autodock_vina": "AutoDock Vina 1.2.5 (Scripps CCSB)",
        "rdkit": rdkit.__version__,
        "gemmi": gemmi.__version__,
        "meeko": "0.5.x (MoleculePreparation / PDBQTWriterLegacy)",
        "python": sys.version.split()[0],
        "platform": sys.platform,
        "utc_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "chembl_rest": "EMBL-EBI ChEMBL REST API v33",
        "rcsb_pdb_rest": "RCSB PDB REST API v1",
        "scoring_function": "AutoDock Vina Iterated Local Search & Monte Carlo"
    })

@app.route("/api/docking/interactions", methods=["POST"])
def analyze_interactions():
    data = request.get_json() or {}
    receptor_pdb = data.get("receptor_pdb", "")
    pose_pdbqt = data.get("pose_pdbqt", "")
    
    if not receptor_pdb or not pose_pdbqt:
        return jsonify({"error": "Both receptor_pdb and pose_pdbqt are required"}), 400

    try:
        contacts = DockingEngine.analyze_interactions(receptor_pdb, pose_pdbqt)
        return jsonify(contacts)
    except Exception as e:
        return jsonify({"error": f"Interaction analysis failed: {str(e)}"}), 400

@app.route("/api/pkpd/adme", methods=["POST"])
def calculate_adme():
    data = request.get_json() or {}
    smiles = data.get("smiles", "").strip()
    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    result = ADMEProfiler.calculate_adme(smiles)
    return jsonify(result)

@app.route("/api/pkpd/crosscheck", methods=["GET"])
def crosscheck_bioactivity():
    drug = request.args.get("drug", "").strip()
    target = request.args.get("target", "").strip()
    if not drug or not target:
        return jsonify({"error": "Parameters 'drug' and 'target' are required"}), 400

    result = BioactivityService.crosscheck_chembl(drug, target)
    return jsonify(result)

@app.route("/api/narrative/explain", methods=["POST"])
def explain_results():
    data = request.get_json() or {}
    api_key = data.get("api_key")
    provider = data.get("provider", "auto")
    
    explanation = NarrativeExplainer.generate_explanation(data, api_key=api_key, provider=provider)
    return jsonify(explanation)

@app.route("/api/batch/run", methods=["POST"])
def batch_docking():
    data = request.get_json() or {}
    receptor_pdbqt = data.get("receptor_pdbqt")
    receptor_pdb = data.get("receptor_pdb")
    center = data.get("center")
    size = data.get("size")
    ligands = data.get("ligands", [])
    exhaustiveness = int(data.get("exhaustiveness", 4))

    if not receptor_pdbqt or not receptor_pdb or not center or not size or not ligands:
        return jsonify({"error": "Missing required batch parameters"}), 400

    try:
        leaderboard = BatchScreeningService.run_batch(
            receptor_pdbqt,
            receptor_pdb,
            center,
            size,
            ligands,
            exhaustiveness=exhaustiveness
        )
        return jsonify({"leaderboard": leaderboard, "total_screened": len(leaderboard)})
    except Exception as e:
        return jsonify({"error": f"Batch docking failed: {str(e)}"}), 500

if __name__ == "__main__":
    print(f"Starting Bindora Server on http://{HOST}:{PORT}")
    app.run(host=HOST, port=PORT, debug=DEBUG)
