import os
import sys
import json
import datetime
import traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from flask import Flask, request, jsonify, send_from_directory, Response
from flask_cors import CORS
from werkzeug.exceptions import HTTPException, BadRequest

from backend.config import (
    BASE_DIR, FRONTEND_DIR, BENCHMARKS_DIR, VINA_EXE, HOST, PORT, DEBUG,
    MAX_CONTENT_LENGTH, CORS_ORIGINS, DATABASE_URL
)
from backend.utils.validators import validate_smiles, validate_pdb_content, validate_grid_box
from backend.db.database import init_db, get_db
from backend.db.models import DockingSession
from backend.routes.session_routes import session_bp
from backend.services.fetcher import StructureFetcher
from backend.services.docking import DockingEngine
from backend.services.adme import ADMEProfiler
from backend.services.bioactivity import BioactivityService
from backend.services.narrative import NarrativeExplainer
from backend.services.batch import BatchScreeningService
from backend.services.batch_manager import BatchScreeningManager
from backend.services.ensemble import EnsembleDockingService
from backend.services.pharmacophore import PharmacophoreService
from backend.services.hardware_profiler import HardwareTelemetrySampler
from backend.utils.vina_setup import ensure_vina
from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__, static_folder=str(FRONTEND_DIR), static_url_path="")
_is_docking_active = False
hw_sampler = HardwareTelemetrySampler.get_instance()

# Security configuration
app.config["MAX_CONTENT_LENGTH"] = MAX_CONTENT_LENGTH
CORS(app, resources={r"/api/*": {"origins": "*"}}, supports_credentials=True)

# Standardized structured JSON error handlers for API endpoints
@app.errorhandler(HTTPException)
def handle_http_exception(e):
    if request.path.startswith("/api/"):
        return jsonify({
            "error": e.description,
            "status": e.code
        }), e.code
    return e.get_response()

@app.errorhandler(Exception)
def handle_generic_exception(e):
    if request.path.startswith("/api/"):
        traceback.print_exc()
        return jsonify({
            "error": f"Internal server error: {str(e)}",
            "status": 500
        }), 500
    raise e

def get_request_json() -> dict:
    """Safely parse request JSON body as a dictionary. Raises BadRequest (400) if malformed or non-dict."""
    if request.content_length and request.content_length > 0:
        try:
            data = request.get_json(force=False, silent=False)
        except Exception as e:
            raise BadRequest(f"Malformed JSON payload: {str(e)}")
    else:
        data = request.get_json(silent=True)
    
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise BadRequest("Request JSON payload must be a JSON object (dictionary)")
    return data


# Ensure Vina binary is ready on startup
try:
    ensure_vina()
except Exception as e:
    print(f"[STARTUP WARNING] Vina initialization warning: {e}", file=sys.stderr)

# Initialize database
try:
    init_db(DATABASE_URL)
    print(f"[DATABASE] Ready at {DATABASE_URL}")
except Exception as e:
    print(f"[STARTUP WARNING] Database initialization warning: {e}", file=sys.stderr)

# Register blueprints
app.register_blueprint(session_bp)

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
        "vina_path": str(VINA_EXE) if vina_ok else None,
        "cpu_count": os.cpu_count() or 1
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

    m17_file = BENCHMARKS_DIR / "1m17_benchmark_result.json"
    m17_data = {}
    if m17_file.exists():
        try:
            with open(m17_file, "r", encoding="utf-8") as f:
                m17_data = json.load(f)
        except Exception:
            pass

    t46_file = BENCHMARKS_DIR / "1t46_benchmark_result.json"
    t46_data = {}
    if t46_file.exists():
        try:
            with open(t46_file, "r", encoding="utf-8") as f:
                t46_data = json.load(f)
        except Exception:
            pass

    iep_file = BENCHMARKS_DIR / "1iep_benchmark_result.json"
    iep_data = {}
    if iep_file.exists():
        try:
            with open(iep_file, "r", encoding="utf-8") as f:
                iep_data = json.load(f)
        except Exception:
            pass

    return jsonify({
        "status": "success",
        "parent_company": "NexPharmaTech",
        "support_email": "sharmaji.pharmatech.info@gmail.com",
        "validation_report": val_data,
        "flagship_targets": {
            "1HSG": hsg_data,
            "1AQ1": aq1_data,
            "1M17": m17_data,
            "1T46": t46_data,
            "1IEP": iep_data
        }
    })

@app.route("/api/search/pubchem", methods=["GET"])
def search_pubchem():
    query = request.args.get("query", "").strip()
    if not query:
        return jsonify({"error": "Query parameter 'query' is required"}), 400
    
    result = StructureFetcher.search_pubchem(query)
    if not result:
        return jsonify({"error": f"No compound found matching '{query}'. Please verify spelling or try pasting a SMILES string directly."}), 404
    
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
    data = get_request_json()
    smiles_or_sdf = data.get("structure", "").strip()
    is_sdf = data.get("is_sdf", False)
    
    if not smiles_or_sdf:
        return jsonify({"error": "Structure content or SMILES string is required"}), 400
    
    # Validate SMILES if not SDF
    if not is_sdf:
        validation = validate_smiles(smiles_or_sdf)
        if not validation.valid:
            return jsonify({"error": validation.error, "field": "structure"}), 400

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
    data = get_request_json()
    pdb_content = data.get("pdb_content", "")
    pdb_id = data.get("pdb_id", "")
    target_chain = data.get("target_chain")

    if not pdb_content and pdb_id:
        fetched = StructureFetcher.fetch_rcsb_pdb(pdb_id)
        if fetched:
            pdb_content = fetched.get("pdb_content", "")

    if not pdb_content:
        return jsonify({"error": "Receptor PDB content or valid PDB ID is required"}), 400
    
    # Validate PDB content
    validation = validate_pdb_content(pdb_content)
    if not validation.valid:
        return jsonify({"error": validation.error, "field": "pdb_content"}), 400

    try:
        retain_waters = bool(data.get("retain_structural_waters", False))
        rec_result = DockingEngine.prepare_receptor(
            pdb_content,
            target_chain=target_chain,
            retain_structural_waters=retain_waters
        )
        return jsonify(rec_result)
    except Exception as e:
        return jsonify({"error": f"Receptor preparation failed: {str(e)}"}), 400

@app.route("/api/docking/run", methods=["POST"])
def run_docking():
    data = get_request_json()
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
    seed = data.get("seed")
    if seed is not None:
        try:
            seed = int(seed)
        except (ValueError, TypeError):
            seed = None
    else:
        seed = None

    try:
        num_modes = int(data.get("num_modes", 9))
        if num_modes < 1 or num_modes > 50:
            num_modes = 9
    except (ValueError, TypeError):
        num_modes = 9

    try:
        replicates = int(data.get("replicates", 1))
        if replicates < 1 or replicates > 10:
            replicates = 1
    except (ValueError, TypeError):
        replicates = 1

    smiles = data.get("smiles", "")
    lig_mol = None
    if smiles:
        lig_mol = Chem.MolFromSmiles(smiles)
        if lig_mol is None:
            return jsonify({
                "error": f"Invalid ligand SMILES '{smiles}': unable to parse into valid chemical structure",
                "heavy_atoms": None,
                "mw": None
            }), 400

    # Derive dynamic heavy atoms and MW strictly from sanitized RDKit molecule (never hardcoded 20 / 300.0)
    if lig_mol is not None:
        heavy_atoms = int(lig_mol.GetNumHeavyAtoms())
        mw = round(float(Descriptors.MolWt(lig_mol)), 2)
    else:
        heavy_atoms = None
        mw = None

    flexible_residues = data.get("flexible_residues")
    cpu = data.get("cpu")
    if cpu is not None:
        try:
            cpu = int(cpu)
        except (ValueError, TypeError):
            cpu = None

    if not receptor_pdbqt or not ligand_pdbqt or not center or not size:
        return jsonify({"error": "Missing required docking parameters (receptor, ligand, center, size)"}), 400
    
    # Validate grid box parameters
    grid_validation = validate_grid_box(center, size)
    if not grid_validation.valid:
        return jsonify({"error": grid_validation.error, "field": "grid_box"}), 400

    global _is_docking_active
    _is_docking_active = True
    try:
        # Run AutoDock Vina with multi-core parallel optimization
        poses = DockingEngine.run_docking(
            receptor_pdbqt,
            ligand_pdbqt,
            center,
            size,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            replicates=replicates,
            seed=seed,
            flexible_residues=flexible_residues,
            receptor_pdb=receptor_pdb,
            cpu=cpu,
            use_gpu=data.get("use_gpu", False),
            ligand_smiles=smiles or data.get("ligand_smiles")
        )

        if not poses:
            return jsonify({"error": "AutoDock Vina finished but returned no binding poses."}), 500

        lig_mol = Chem.MolFromSmiles(smiles) if smiles else None

        # Pre-compute biophysical interactions for ALL poses so switching modes is instant
        for p in poses:
            try:
                p_contacts = DockingEngine.analyze_interactions(receptor_pdb, p["pdbqt_content"], ligand_mol=lig_mol)
                if smiles:
                    try:
                        from backend.services.interaction_diagram import InteractionDiagramGenerator
                        p_contacts["diagram_svg"] = InteractionDiagramGenerator.generate_diagram_svg(smiles, p_contacts)
                    except Exception:
                        pass
                p["interactions"] = p_contacts
            except Exception as pe:
                p["interactions"] = {}

        best_pose = poses[0]
        affinity = best_pose["affinity_kcal"]
        replicate_stats = best_pose.get("replicate_stats")

        # Calculate thermodynamics and Ligand Efficiency
        thermo = BioactivityService.calculate_thermodynamics(affinity, heavy_atoms, mw)
        contacts = best_pose.get("interactions", {})
        
        # Classify docking experiment (Native Redocking vs Cross-Docking vs Targeted Docking)
        native_ligand_info = data.get("native_ligand_info")
        ligand_name = data.get("ligand_name", "Investigational Ligand")
        experiment_classification = DockingEngine.classify_docking_experiment(
            docked_smiles=smiles,
            native_ligand_info=native_ligand_info,
            docked_ligand_name=ligand_name
        )

        # Record session to database
        session_id = None
        try:
            with get_db() as db:
                pdb_id = data.get("pdb_id", "")
                
                session = DockingSession(
                    pdb_id=pdb_id if pdb_id else None,
                    ligand_name=ligand_name,
                    ligand_smiles=smiles if smiles else None,
                    affinity_kcal=affinity,
                    rmsd_lb=best_pose.get("rmsd_lower_bound"),
                    rmsd_ub=best_pose.get("rmsd_upper_bound"),
                    execution_duration_s=best_pose.get("execution_duration_s", 0.0),
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes
                )
                db.add(session)
                db.commit()
                session_id = session.id
                print(f"[DATABASE] Recorded docking session: {session_id}")
        except Exception as db_err:
            print(f"[DATABASE WARNING] Failed to record session: {db_err}", file=sys.stderr)

        # Compute cryptographic provenance hash for auditable reproducibility
        prov_hash = None
        try:
            from backend.utils.report_emitter import compute_provenance_hash
            prov_hash = "SHA256:" + compute_provenance_hash({
                "affinity": affinity,
                "poses_count": len(poses),
                "seed": best_pose.get("seed_used"),
                "exhaustiveness": exhaustiveness,
                "pdb_id": data.get("pdb_id", "") or "PDB",
                "ligand": ligand_name
            })
        except Exception:
            prov_hash = None

        return jsonify({
            "poses": poses,
            "top_pose": best_pose,
            "seed_used": best_pose.get("seed_used"),
            "exhaustiveness_used": best_pose.get("exhaustiveness_used", exhaustiveness),
            "adaptive_sampling": best_pose.get("adaptive_sampling_note"),
            "thermodynamics": thermo,
            "interactions": contacts,
            "replicate_stats": replicate_stats,
            "execution_duration_s": best_pose.get("execution_duration_s", 0.0),
            "cpu_count": best_pose.get("cpu_count_used", os.cpu_count() or 1),
            "session_id": session_id,
            "experiment_classification": experiment_classification,
            "provenance_hash": prov_hash,
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat()
        })
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Docking execution failed: {str(e)}"}), 500
    finally:
        _is_docking_active = False

@app.route("/api/docking/classify-experiment", methods=["POST"])
def classify_experiment_endpoint():
    data = get_request_json()
    smiles = data.get("smiles", "")
    native_ligand_info = data.get("native_ligand_info")
    ligand_name = data.get("ligand_name", "")
    if not smiles and not native_ligand_info and not ligand_name:
        return jsonify({"error": "Missing ligand SMILES, name, or native_ligand_info"}), 400
    try:
        classification = DockingEngine.classify_docking_experiment(
            docked_smiles=smiles,
            native_ligand_info=native_ligand_info,
            docked_ligand_name=ligand_name
        )
        return jsonify(classification)
    except Exception as e:
        return jsonify({"error": f"Experiment classification failed: {str(e)}"}), 400

@app.route("/api/docking/analyze-interactions", methods=["POST"])
def analyze_interactions_endpoint():
    data = get_request_json()
    receptor_pdb = data.get("receptor_pdb", "")
    pose_pdbqt = data.get("pose_pdbqt", "")
    smiles = data.get("smiles", "")
    if not receptor_pdb or not pose_pdbqt:
        return jsonify({"error": "Missing 'receptor_pdb' or 'pose_pdbqt'"}), 400
    try:
        lig_mol = Chem.MolFromSmiles(smiles) if smiles else None
        contacts = DockingEngine.analyze_interactions(receptor_pdb, pose_pdbqt, ligand_mol=lig_mol)
        if smiles:
            try:
                from backend.services.interaction_diagram import InteractionDiagramGenerator
                contacts["diagram_svg"] = InteractionDiagramGenerator.generate_diagram_svg(smiles, contacts)
            except Exception:
                pass
        return jsonify(contacts)
    except Exception as e:
        return jsonify({"error": f"Interaction analysis failed: {str(e)}"}), 400

@app.route("/api/docking/interaction-diagram", methods=["POST"])
def get_interaction_diagram():
    data = get_request_json()
    smiles = data.get("smiles", "")
    interactions = data.get("interactions", {})
    if not smiles or not interactions:
        return jsonify({"error": "Missing 'smiles' or 'interactions' in request body"}), 400
    if not isinstance(interactions, dict):
        return jsonify({"error": "'interactions' must be a dictionary"}), 400
    try:
        from backend.services.interaction_diagram import InteractionDiagramGenerator
        svg = InteractionDiagramGenerator.generate_diagram_svg(smiles, interactions)
        return jsonify({"diagram_svg": svg})
    except Exception as e:
        return jsonify({"error": f"Failed to generate diagram: {str(e)}"}), 400

@app.route("/api/docking/refine", methods=["POST"])
def refine_docked_pose():
    data = get_request_json()
    receptor_pdb = data.get("receptor_pdb", "")
    docked_pdb = data.get("docked_pdb", "")
    smiles = data.get("smiles", "")

    if not receptor_pdb or not docked_pdb:
        return jsonify({"error": "Both receptor_pdb and docked_pdb are required"}), 400

    try:
        from backend.services.refinement import ComplexRefinementService
        refinement_result = ComplexRefinementService.refine_pose(
            receptor_pdb=receptor_pdb,
            docked_pdb_or_pdbqt=docked_pdb,
            smiles=smiles
        )
        return jsonify(refinement_result)
    except Exception as e:
        return jsonify({"error": f"Pose refinement failed: {str(e)}"}), 400

@app.route("/api/docking/redock-validate", methods=["POST"])
def redock_validate():
    data = get_request_json()
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

    grid_val = validate_grid_box(center, size)
    if not grid_val.valid:
        return jsonify({"error": grid_val.error, "field": "grid_box"}), 400

    raw_seed = data.get("seed")
    try:
        seed = int(raw_seed) if raw_seed is not None else None
    except (ValueError, TypeError):
        seed = None

    global _is_docking_active
    _is_docking_active = True
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
        return jsonify({"error": f"Redocking validation failed: {str(e)}"}), 400
    finally:
        _is_docking_active = False

@app.route("/api/reproducibility/versions", methods=["GET"])
def get_reproducibility_versions():
    import gemmi
    import rdkit
    import datetime
    try:
        import sklearn
        sklearn_ver = sklearn.__version__
    except Exception:
        sklearn_ver = "1.9+"
    try:
        import openmm
        openmm_ver = openmm.__version__
    except Exception:
        openmm_ver = "8.6+"
    return jsonify({
        "bindora_dock": "v2.2 (Research-Grade)",
        "autodock_vina": "AutoDock Vina 1.2.5 / 1.2.7 (Scripps CCSB)",
        "solvation_engine": f"OpenMM MM-GBSA v{openmm_ver} (OBC2 / AMBER99SB-ILDN / GAFF2)",
        "ml_engine": f"Supervised ML Ensemble v{sklearn_ver} (Wang et al. 2011)",
        "rdkit": rdkit.__version__,
        "gemmi": gemmi.__version__,
        "meeko": "0.5.x (3-Tier Flexible Torsion Engine)",
        "python": sys.version.split()[0],
        "platform": sys.platform,
        "utc_timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "chembl_rest": "EMBL-EBI ChEMBL REST API v33",
        "rcsb_pdb_rest": "RCSB PDB REST API v1",
        "scoring_function": "AutoDock Vina & Vinardo Iterated Local Search",
        "provenance_standard": "SHA-256 Deterministic Provenance & Zero Hardcoding"
    })

@app.route("/api/system/hardware", methods=["GET"])
def get_hardware_info():
    from backend.services.docking import _ACTIVE_SUBPROCESSES
    is_docking = bool(_is_docking_active or len(_ACTIVE_SUBPROCESSES) > 0)
    info = hw_sampler.get_telemetry(is_docking_active=is_docking)
    return jsonify(info)

@app.route("/api/docking/interactions", methods=["POST"])
def analyze_interactions():
    data = get_request_json()
    receptor_pdb = data.get("receptor_pdb", "")
    pose_pdbqt = data.get("pose_pdbqt", "")
    
    if not receptor_pdb or not pose_pdbqt:
        return jsonify({"error": "Both receptor_pdb and pose_pdbqt are required"}), 400

    try:
        smiles = data.get("smiles", "")
        if isinstance(smiles, str):
            smiles = smiles.strip()
        else:
            smiles = ""
        lig_mol = Chem.MolFromSmiles(smiles) if smiles else None
        contacts = DockingEngine.analyze_interactions(receptor_pdb, pose_pdbqt, ligand_mol=lig_mol)
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
@app.route("/api/adme/profile", methods=["POST"])
def calculate_adme():
    data = get_request_json()
    smiles = data.get("smiles", "")
    if isinstance(smiles, str):
        smiles = smiles.strip()
    else:
        smiles = ""
    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400
    
    # Validate SMILES
    validation = validate_smiles(smiles)
    if not validation.valid:
        return jsonify({"error": validation.error, "field": "smiles"}), 400

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

@app.route("/api/narrative/explain", methods=["POST", "GET"])
@app.route("/api/explain/narrative", methods=["POST", "GET"])
@app.route("/api/narrative", methods=["POST", "GET"])
def explain_results():
    if request.method == "GET":
        data = request.args.to_dict()
    else:
        data = get_request_json()
    api_key = data.get("api_key")
    provider = data.get("provider", "auto")
    
    explanation = NarrativeExplainer.generate_explanation(data, api_key=api_key, provider=provider)
    return jsonify(explanation)

@app.route("/api/batch/start", methods=["POST"])
def start_batch_docking():
    data = get_request_json()
    receptor_pdbqt = data.get("receptor_pdbqt")
    receptor_pdb = data.get("receptor_pdb")
    center = data.get("center")
    size = data.get("size")
    ligands = data.get("ligands", [])
    try:
        exhaustiveness = int(data.get("exhaustiveness", 4))
    except (ValueError, TypeError):
        exhaustiveness = 4

    if not receptor_pdbqt or not receptor_pdb or not center or not size or not ligands:
        return jsonify({"error": "Missing required batch parameters"}), 400

    manager = BatchScreeningManager.get_instance()
    job_id = manager.start_batch_job(
        receptor_pdbqt=receptor_pdbqt,
        receptor_pdb=receptor_pdb,
        pocket_center=center,
        pocket_size=size,
        ligand_list=ligands,
        exhaustiveness=exhaustiveness
    )
    return jsonify({"job_id": job_id, "total": len(ligands), "status": "queued"})

@app.route("/api/batch/status/<job_id>", methods=["GET"])
def get_batch_status(job_id):
    manager = BatchScreeningManager.get_instance()
    job = manager.get_job(job_id)
    if not job:
        return jsonify({"error": "Job not found"}), 404
    return jsonify(job)

@app.route("/api/batch/stream/<job_id>", methods=["GET"])
def stream_batch_events(job_id):
    manager = BatchScreeningManager.get_instance()
    return Response(manager.stream_job_events(job_id), mimetype="text/event-stream")

@app.route("/api/batch/cancel/<job_id>", methods=["POST"])
def cancel_batch_docking(job_id):
    manager = BatchScreeningManager.get_instance()
    ok = manager.cancel_job(job_id)
    return jsonify({"cancelled": ok, "job_id": job_id})

@app.route("/api/batch/run", methods=["POST"])
@app.route("/api/batch-screen", methods=["POST"])
def batch_docking():
    data = get_request_json()
    receptor_pdbqt = data.get("receptor_pdbqt")
    receptor_pdb = data.get("receptor_pdb")
    center = data.get("center")
    size = data.get("size")
    ligands = data.get("ligands", [])
    try:
        exhaustiveness = int(data.get("exhaustiveness", 4))
        if exhaustiveness < 1 or exhaustiveness > 64:
            exhaustiveness = 4
    except (ValueError, TypeError):
        exhaustiveness = 4

    if not receptor_pdbqt or not receptor_pdb or not center or not size or not ligands:
        return jsonify({"error": "Missing required batch parameters (receptor_pdbqt, receptor_pdb, center, size, ligands)"}), 400

    if not isinstance(ligands, list) or len(ligands) == 0:
        return jsonify({"error": "Parameter 'ligands' must be a non-empty list"}), 400

    grid_val = validate_grid_box(center, size)
    if not grid_val.valid:
        return jsonify({"error": grid_val.error, "field": "grid_box"}), 400

    global _is_docking_active
    _is_docking_active = True
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
        return jsonify({"error": f"Batch docking failed: {str(e)}"}), 400
    finally:
        _is_docking_active = False

@app.route("/api/ensemble/structures", methods=["GET"])
def get_ensemble_structures():
    accession = request.args.get("accession", "").strip()
    if not accession:
        return jsonify({"error": "Query parameter 'accession' is required"}), 400
    structures = EnsembleDockingService.fetch_ensemble_structures(accession, limit=8)
    return jsonify({"structures": structures, "total": len(structures)})

@app.route("/api/ensemble/run", methods=["POST"])
def run_ensemble():
    data = get_request_json()
    pdb_ids = data.get("pdb_ids", [])
    ligand_smiles = data.get("ligand_smiles", "")
    if isinstance(ligand_smiles, str):
        ligand_smiles = ligand_smiles.strip()
    else:
        ligand_smiles = ""
    try:
        exhaustiveness = int(data.get("exhaustiveness", 4))
    except (ValueError, TypeError):
        exhaustiveness = 4

    if not pdb_ids or not ligand_smiles:
        return jsonify({"error": "Missing 'pdb_ids' or 'ligand_smiles' in request body"}), 400

    global _is_docking_active
    _is_docking_active = True
    try:
        results = EnsembleDockingService.run_ensemble_docking(pdb_ids, ligand_smiles, exhaustiveness=exhaustiveness)
        return jsonify(results)
    finally:
        _is_docking_active = False

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
    data = get_request_json()
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

# ==========================================
# PHASE 2 API ENDPOINTS (Days 61-120)
# ==========================================

@app.route("/api/covalent/evaluate", methods=["POST"])
def evaluate_covalent():
    """Evaluate covalent binding feasibility for an electrophilic ligand against receptor nucleophiles."""
    data = get_request_json()
    docked_pose = data.get("docked_pose_pdbqt") or data.get("docked_pose_pdb", "")
    receptor = data.get("receptor_pdbqt") or data.get("receptor_pdb", "")
    smiles = data.get("smiles", "")
    if isinstance(smiles, str):
        smiles = smiles.strip() or None
    else:
        smiles = None
    pocket_center = data.get("pocket_center")

    if not docked_pose or not receptor:
        return jsonify({"error": "Missing 'docked_pose_pdbqt'/'docked_pose_pdb' or 'receptor_pdbqt'/'receptor_pdb'"}), 400

    from backend.services.covalent import CovalentDockingService
    result = CovalentDockingService.evaluate_covalent_geometry(
        docked_pose_pdb_or_pdbqt=docked_pose,
        receptor_pdb_or_pdbqt=receptor,
        smiles=smiles,
        pocket_center=pocket_center
    )
    return jsonify(result)

@app.route("/api/macrocycle/sample", methods=["POST"])
def sample_macrocycle():
    """Sample conformational ensemble for macrocycle using srETKDGv3 and MMFF94."""
    data = get_request_json()
    smiles = data.get("smiles", "")
    if isinstance(smiles, str):
        smiles = smiles.strip()
    else:
        smiles = ""
    try:
        num_confs = int(data.get("num_confs", 20))
    except (ValueError, TypeError):
        num_confs = 20
    try:
        energy_window = float(data.get("energy_window", 15.0))
    except (ValueError, TypeError):
        energy_window = 15.0
    try:
        rmsd_threshold = float(data.get("rmsd_threshold", 0.5))
    except (ValueError, TypeError):
        rmsd_threshold = 0.5

    if not smiles:
        return jsonify({"error": "Missing 'smiles' in request body"}), 400

    from backend.services.macrocycle import MacrocycleConformerEngine
    try:
        result = MacrocycleConformerEngine.sample_macrocycle_conformers(
            smiles,
            num_confs=num_confs,
            energy_window=energy_window,
            rmsd_threshold=rmsd_threshold
        )
        # Exclude non-serializable mol objects from direct json response
        return jsonify({
            "is_macrocycle": result["is_macrocycle"],
            "max_ring_size": result["max_ring_size"],
            "macrocycle_ring_sizes": result["macrocycle_ring_sizes"],
            "sampling_engine": result["sampling_engine"],
            "force_field": result["force_field"],
            "initial_conformers_sampled": result["initial_conformers_sampled"],
            "conformers_within_energy_window": result["conformers_within_energy_window"],
            "conformers_retained_after_rmsd_pruning": result["conformers_retained_after_rmsd_pruning"],
            "global_min_energy_kcal": result["global_min_energy_kcal"],
            "relative_energies_kcal": result["relative_energies_kcal"]
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 400

@app.route("/api/openmm/export", methods=["POST"])
def export_openmm_script():
    """Generate standalone OpenMM explicit-solvent MD simulation package & MM-PBSA script."""
    data = get_request_json()
    receptor_pdb = data.get("receptor_pdb", "")
    docked_pose = data.get("docked_pose_pdbqt") or data.get("docked_pose_pdb", "")
    output_dir = data.get("output_dir", "")
    smiles = data.get("smiles", "")
    if isinstance(smiles, str):
        smiles = smiles.strip() or None
    else:
        smiles = None
    job_name = data.get("job_name", "bindora_complex_md")
    try:
        sim_time_ns = float(data.get("sim_time_ns", 1.0))
    except (ValueError, TypeError):
        sim_time_ns = 1.0

    if not receptor_pdb or not docked_pose:
        return jsonify({"error": "Missing 'receptor_pdb' or 'docked_pose_pdbqt'/'docked_pose_pdb'"}), 400

    if not output_dir:
        import tempfile
        output_dir = tempfile.mkdtemp(prefix="bindora_openmm_")

    from backend.services.md_export import OpenMMExportService
    try:
        pkg = OpenMMExportService.generate_simulation_package(
            receptor_pdb=receptor_pdb,
            docked_pose_pdb_or_pdbqt=docked_pose,
            output_dir=output_dir,
            ligand_smiles=smiles,
            job_name=job_name,
            sim_time_ns=sim_time_ns
        )
        return jsonify(pkg)
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route("/api/benchmark/pdbbind", methods=["GET", "POST"])
def pdbbind_benchmark():
    """Query PDBbind core set reference complexes or evaluate docking predictions against benchmark."""
    from backend.services.pdbbind_validation import PDBbindValidationEngine
    if request.method == "GET":
        ref_data = PDBbindValidationEngine.get_reference_dataset()
        return jsonify({
            "dataset_name": "CASF-2016 / PDBbind Core Set Curated Benchmark",
            "complexes_count": len(ref_data),
            "complexes": ref_data
        })
    else:
        data = get_request_json()
        predictions = data.get("predictions", [])
        if not predictions:
            return jsonify({"error": "Missing 'predictions' array in request body"}), 400
        result = PDBbindValidationEngine.evaluate_benchmark(predictions)
        return jsonify(result)

@app.route("/api/docking/induced-fit", methods=["POST"])
@app.route("/api/induced-fit", methods=["POST"])
def run_induced_fit():
    """Execute Monte Carlo loop and backbone phi/psi induced-fit docking (IFD)."""
    data = get_request_json()
    receptor_pdb = data.get("receptor_pdb", "")
    ligand_sdf_or_pdbqt = data.get("ligand", "") or data.get("ligand_sdf_or_pdbqt", "") or data.get("ligand_smiles", "")
    pocket_center = data.get("center", {})
    pocket_size = data.get("size", {})
    try:
        exhaustiveness = int(data.get("exhaustiveness", 8))
        if exhaustiveness < 1 or exhaustiveness > 64:
            exhaustiveness = 8
    except (ValueError, TypeError):
        exhaustiveness = 8

    try:
        loop_radius = float(data.get("loop_radius", 8.5))
        if loop_radius <= 0 or loop_radius > 50.0:
            loop_radius = 8.5
    except (ValueError, TypeError):
        loop_radius = 8.5

    try:
        num_iterations = int(data.get("num_iterations", 5))
        if num_iterations < 1 or num_iterations > 50:
            num_iterations = 5
    except (ValueError, TypeError):
        num_iterations = 5

    raw_seed = data.get("seed")
    try:
        seed = int(raw_seed) if raw_seed is not None else None
    except (ValueError, TypeError):
        seed = None

    if not receptor_pdb or not ligand_sdf_or_pdbqt or not pocket_center or not pocket_size:
        return jsonify({"error": "Missing required fields (receptor_pdb, ligand, center, size)"}), 400

    grid_val = validate_grid_box(pocket_center, pocket_size)
    if not grid_val.valid:
        return jsonify({"error": grid_val.error, "field": "grid_box"}), 400

    from backend.services.induced_fit import InducedFitService
    try:
        result = InducedFitService.run_induced_fit_docking(
            receptor_pdb=receptor_pdb,
            ligand_input=ligand_sdf_or_pdbqt,
            pocket_center=pocket_center,
            pocket_size=pocket_size,
            loop_radius=loop_radius,
            num_conformations=num_iterations,
            exhaustiveness=exhaustiveness,
            seed=seed
        )
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": f"Induced-fit docking failed: {str(e)}"}), 400

@app.route("/api/covalent/build-adduct", methods=["POST"])
def build_covalent_adduct():
    """Build physical covalent adduct complex with bidirectional CONECT records."""
    data = get_request_json()
    receptor_pdb = data.get("receptor_pdb", "")
    ligand_pose_pdbqt = data.get("pose_pdbqt", "")
    warhead_type = data.get("warhead_type", None)
    target_residue = data.get("target_residue", None)

    if not receptor_pdb or not ligand_pose_pdbqt:
        return jsonify({"error": "Missing 'receptor_pdb' or 'pose_pdbqt'"}), 400

    from backend.services.covalent import CovalentDockingEngine
    try:
        result = CovalentDockingEngine.build_covalent_adduct_complex(
            receptor_pdb=receptor_pdb,
            ligand_pose_pdbqt=ligand_pose_pdbqt,
            warhead_type=warhead_type,
            target_residue=target_residue
        )
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    print(f"Starting Bindora Server on http://{HOST}:{PORT}")
    app.run(host=HOST, port=PORT, debug=DEBUG, use_reloader=False)
