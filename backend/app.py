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
from backend.services.ensemble import EnsembleDockingService
from backend.services.pharmacophore import PharmacophoreService
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

@app.route("/api/validation-report", methods=["GET"])
@app.route("/api/validation-reports", methods=["GET"])
def get_validation_report():
    val_file = BENCHMARKS_DIR / "validation_report_v1.json"
    val_data = {}
    if val_file.exists():
        try:
            with open(val_file, "r", encoding="utf-8") as f:
                val_data = json.load(f)
        except Exception as e:
            val_data = {"error": str(e)}

    hsg_file = BENCHMARKS_DIR / "1hsg_benchmark_result.json"
    hsg_data = {}
    if hsg_file.exists():
        try:
            with open(hsg_file, "r", encoding="utf-8") as f:
                hsg_data = json.load(f)
        except Exception:
            pass

    aq1_file = BENCHMARKS_DIR / "1aq1_benchmark_result.json"
    aq1_data = {}
    if aq1_file.exists():
        try:
            with open(aq1_file, "r", encoding="utf-8") as f:
                aq1_data = json.load(f)
        except Exception:
            pass

    return jsonify({
        "status": "success",
        "parent_company": "NexPharmaTech",
        "support_email": "sharmaji.pharmatech.info@gmail.com",
        "validation_report": val_data,
        "flagship_targets": {
            "1HSG": hsg_data,
            "1AQ1": aq1_data
        }
    })

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

@app.route("/api/ligand/similar", methods=["GET"])
def get_similar_compounds():
    smiles = request.args.get("smiles", "").strip()
    if not smiles:
        return jsonify({"error": "Query parameter 'smiles' is required"}), 400
    try:
        threshold = int(request.args.get("threshold", 85))
    except (ValueError, TypeError):
        threshold = 85
    try:
        max_records = int(request.args.get("max", 5))
    except (ValueError, TypeError):
        max_records = 5

    results = StructureFetcher.search_similar_compounds(smiles, threshold=threshold, max_records=max_records)
    return jsonify({"similar_compounds": results, "total": len(results)})

@app.route("/api/search/rcsb", methods=["GET"])
def search_rcsb():
    query = request.args.get("query", "").strip()
    if not query:
        return jsonify({"error": "Query parameter 'query' is required"}), 400
    
    # If 4 characters, treat as direct PDB ID
    if len(query) == 4 and query.isalnum():
        result = StructureFetcher.fetch_rcsb_pdb(query)
        if result:
            # Also fetch UniProt annotation via PDB cross-reference if not already populated
            if "uniprot" not in result or not result["uniprot"]:
                uniprot_info = StructureFetcher.fetch_uniprot_by_pdb_id(query)
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
    try:
        exhaustiveness = int(data.get("exhaustiveness", 8))
        if exhaustiveness < 1 or exhaustiveness > 64:
            exhaustiveness = 8
    except (ValueError, TypeError):
        exhaustiveness = 8
    num_modes = int(data.get("num_modes", 9))
    replicates = int(data.get("replicates", 1))
    heavy_atoms = int(data.get("heavy_atoms", 20))
    mw = float(data.get("molecular_weight", 300.0))
    smiles = data.get("smiles", "")
    flexible_residues = data.get("flexible_residues")

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
            replicates=replicates,
            flexible_residues=flexible_residues,
            receptor_pdb=receptor_pdb
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
            if smiles:
                try:
                    from backend.services.interaction_diagram import InteractionDiagramGenerator
                    contacts["diagram_svg"] = InteractionDiagramGenerator.generate_diagram_svg(smiles, contacts)
                except Exception as ex:
                    print(f"[DIAGRAM ERROR] {ex}")

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

@app.route("/api/docking/interaction-diagram", methods=["POST"])
def get_interaction_diagram():
    data = request.get_json() or {}
    smiles = data.get("smiles", "")
    interactions = data.get("interactions", {})
    if not smiles or not interactions:
        return jsonify({"error": "Missing 'smiles' or 'interactions' in request body"}), 400
    try:
        from backend.services.interaction_diagram import InteractionDiagramGenerator
        svg = InteractionDiagramGenerator.generate_diagram_svg(smiles, interactions)
        return jsonify({"diagram_svg": svg})
    except Exception as e:
        return jsonify({"error": f"Failed to generate diagram: {str(e)}"}), 500

@app.route("/api/docking/redock-validate", methods=["POST"])
def redock_validate():
    data = request.get_json() or {}
    receptor_pdbqt = data.get("receptor_pdbqt")
    native_ligand_pdb = data.get("native_ligand_pdb")
    center = data.get("center")
    size = data.get("size")
    try:
        exhaustiveness = int(data.get("exhaustiveness", 8))
        if exhaustiveness < 1 or exhaustiveness > 64:
            exhaustiveness = 8
    except (ValueError, TypeError):
        exhaustiveness = 8

    if not receptor_pdbqt or not native_ligand_pdb or not center or not size:
        return jsonify({"error": "Missing required parameters for redocking validation (receptor_pdbqt, native_ligand_pdb, center, size)"}), 400

    seed = int(data.get("seed", 42))

    try:
        validation_result = DockingEngine.run_redocking_validation(
            receptor_pdbqt,
            native_ligand_pdb,
            center,
            size,
            exhaustiveness=exhaustiveness,
            seed=seed
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
        smiles = data.get("smiles", "").strip()
        if smiles:
            try:
                from backend.services.interaction_diagram import InteractionDiagramGenerator
                contacts["diagram_svg"] = InteractionDiagramGenerator.generate_diagram_svg(smiles, contacts)
            except Exception as se:
                contacts["diagram_svg_error"] = str(se)
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

@app.route("/api/ensemble/structures", methods=["GET"])
def get_ensemble_structures():
    accession = request.args.get("accession", "").strip()
    if not accession:
        return jsonify({"error": "Query parameter 'accession' is required"}), 400
    structures = EnsembleDockingService.fetch_ensemble_structures(accession, limit=8)
    return jsonify({"structures": structures, "total": len(structures)})

@app.route("/api/ensemble/run", methods=["POST"])
def run_ensemble():
    data = request.get_json() or {}
    pdb_ids = data.get("pdb_ids", [])
    ligand_smiles = data.get("ligand_smiles", "").strip()
    exhaustiveness = int(data.get("exhaustiveness", 4))

    if not pdb_ids or not ligand_smiles:
        return jsonify({"error": "Missing 'pdb_ids' or 'ligand_smiles' in request body"}), 400

    results = EnsembleDockingService.run_ensemble_docking(pdb_ids, ligand_smiles, exhaustiveness=exhaustiveness)
    return jsonify(results)

@app.route("/api/pharmacophore/actives", methods=["GET"])
def get_pharmacophore_actives():
    target = request.args.get("target", "").strip()
    chembl_id = request.args.get("chembl_id", "").strip() or None
    pdb_id = request.args.get("pdb_id", "").strip() or None
    uniprot_acc = request.args.get("uniprot_acc", "").strip() or None
    max_actives = int(request.args.get("max_actives", 10))

    if not target and not chembl_id and not pdb_id and not uniprot_acc:
        return jsonify({"error": "Target identifier (target, chembl_id, pdb_id, or uniprot_acc) is required"}), 400

    actives = PharmacophoreService.fetch_target_actives(
        target_name=target,
        chembl_target_id=chembl_id,
        pdb_id=pdb_id,
        uniprot_accession=uniprot_acc,
        max_actives=max_actives
    )
    profile = PharmacophoreService.build_consensus_profile(actives) if len(actives) >= 3 else None
    return jsonify({
        "actives_count": len(actives),
        "eligible": len(actives) >= 3,
        "actives": actives,
        "consensus_profile": profile
    })

@app.route("/api/pharmacophore/screen", methods=["POST"])
def screen_pharmacophore():
    data = request.get_json() or {}
    candidates = data.get("candidates", [])
    consensus_profile = data.get("consensus_profile")

    if not candidates or not consensus_profile:
        return jsonify({"error": "Missing 'candidates' or 'consensus_profile' in request body"}), 400

    results = []
    for cand in candidates:
        smi = cand.get("smiles", "")
        name = cand.get("name", "Candidate")
        match = PharmacophoreService.match_candidate(smi, consensus_profile)
        results.append({
            "name": name,
            "smiles": smi,
            **match
        })

    results.sort(key=lambda x: x.get("match_score_pct", 0.0), reverse=True)
    return jsonify({"screened": results, "total": len(results)})

if __name__ == "__main__":
    print(f"Starting Bindora Server on http://{HOST}:{PORT}")
    app.run(host=HOST, port=PORT, debug=DEBUG, use_reloader=False)
