import json
import time
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Dict, Any, Optional, List
from backend.config import (
    CACHE_DIR,
    PUBCHEM_BASE_URL,
    RCSB_DATA_URL,
    RCSB_FILE_URL,
    RCSB_SEARCH_URL,
    CHEMBL_BASE_URL,
    UNIPROT_BASE_URL,
)

HEADERS = {"User-Agent": "Bindora-Research-Tool/1.0 (academic; +https://github.com/Swelo-ui/Bindora)"}

def _get_json(url: str, timeout: int = 15) -> Optional[Dict[str, Any]]:
    """Helper to fetch and parse JSON with error handling."""
    try:
        req = urllib.request.Request(url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return None

def _get_text(url: str, timeout: int = 20) -> Optional[str]:
    """Helper to fetch raw text content."""
    try:
        req = urllib.request.Request(url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except Exception as e:
        return None

class StructureFetcher:
    """Service to search and retrieve chemical and macromolecular structures from public repositories."""

    @staticmethod
    def search_pubchem(name: str) -> Optional[Dict[str, Any]]:
        """Search compound by name and retrieve properties, SMILES, and 3D/2D coordinates."""
        cache_file = CACHE_DIR / f"pubchem_{urllib.parse.quote_plus(name.lower())}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                pass

        # 1. Fetch properties and CID
        prop_url = (
            f"{PUBCHEM_BASE_URL}/compound/name/{urllib.parse.quote(name)}/property/"
            f"MolecularFormula,MolecularWeight,CanonicalSMILES,ConnectivitySMILES,IUPACName,"
            f"XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
        )
        data = _get_json(prop_url)
        if not data or "PropertyTable" not in data or not data["PropertyTable"]["Properties"]:
            return None

        props = data["PropertyTable"]["Properties"][0]
        cid = props.get("CID")
        smiles = props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")

        # 2. Try fetching 3D SDF first, fallback to 2D SDF
        sdf_3d_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/SDF?record_type=3d"
        sdf_content = _get_text(sdf_3d_url)
        is_3d = True
        if not sdf_content or "CONECT" not in sdf_content:
            sdf_2d_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/SDF"
            sdf_content = _get_text(sdf_2d_url)
            is_3d = False

        result = {
            "source": "PubChem",
            "cid": cid,
            "name": name.capitalize(),
            "iupac_name": props.get("IUPACName", ""),
            "formula": props.get("MolecularFormula", ""),
            "weight": float(props.get("MolecularWeight", 0.0)),
            "smiles": smiles,
            "xlogp": float(props.get("XLogP", 0.0)) if props.get("XLogP") is not None else None,
            "tpsa": float(props.get("TPSA", 0.0)) if props.get("TPSA") is not None else None,
            "hbd": int(props.get("HBondDonorCount", 0)),
            "hba": int(props.get("HBondAcceptorCount", 0)),
            "rotb": int(props.get("RotatableBondCount", 0)),
            "sdf": sdf_content,
            "is_3d": is_3d,
            "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        }

        # Cache result
        try:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2)
        except Exception:
            pass

        return result

    @staticmethod
    def fetch_pubchem_by_cid(cid: str) -> Optional[Dict[str, Any]]:
        """Fetch compound metadata, title, and properties by PubChem CID."""
        cid = str(cid).strip()
        if not cid.isdigit():
            return None
        cache_file = CACHE_DIR / f"pubchem_cid_{cid}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                pass

        prop_url = (
            f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/"
            f"Title,MolecularFormula,MolecularWeight,CanonicalSMILES,ConnectivitySMILES,IUPACName,"
            f"XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
        )
        data = _get_json(prop_url)
        if not data or "PropertyTable" not in data or not data["PropertyTable"]["Properties"]:
            return None

        props = data["PropertyTable"]["Properties"][0]
        smiles = props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
        title = props.get("Title") or props.get("IUPACName") or f"Compound #{cid}"

        result = {
            "source": "PubChem",
            "cid": int(cid),
            "name": title,
            "iupac_name": props.get("IUPACName", ""),
            "formula": props.get("MolecularFormula", ""),
            "weight": float(props.get("MolecularWeight", 0.0)),
            "smiles": smiles,
            "xlogp": float(props.get("XLogP", 0.0)) if props.get("XLogP") is not None else None,
            "tpsa": float(props.get("TPSA", 0.0)) if props.get("TPSA") is not None else None,
            "hbd": int(props.get("HBondDonorCount", 0)),
            "hba": int(props.get("HBondAcceptorCount", 0)),
            "rotb": int(props.get("RotatableBondCount", 0)),
            "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        }
        try:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2)
        except Exception:
            pass
        return result

    @staticmethod
    def fetch_rcsb_pdb(pdb_id: str) -> Optional[Dict[str, Any]]:
        """Fetch macromolecule PDB file and metadata by 4-letter PDB code."""
        pdb_id = pdb_id.strip().upper()
        if len(pdb_id) != 4:
            return None

        cache_pdb = CACHE_DIR / f"{pdb_id}.pdb"
        cache_meta = CACHE_DIR / f"{pdb_id}_meta.json"

        metadata = {}
        if cache_meta.exists():
            try:
                with open(cache_meta, "r", encoding="utf-8") as f:
                    metadata = json.load(f)
            except Exception:
                pass

        if not metadata:
            meta_url = f"{RCSB_DATA_URL}/{pdb_id}"
            raw_meta = _get_json(meta_url)
            if raw_meta:
                title = raw_meta.get("struct", {}).get("title", "")
                exptl = raw_meta.get("exptl", [{}])[0].get("method", "Unknown")
                res = raw_meta.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
                dep_date = raw_meta.get("rcsb_accession_info", {}).get("deposit_date", "")
                organism = "Homo sapiens"
                try:
                    organism = raw_meta["rcsb_entry_container_identifiers"]["source_organism_scientific_names"][0]
                except Exception:
                    pass
                metadata = {
                    "pdb_id": pdb_id,
                    "title": title,
                    "method": exptl,
                    "resolution": res,
                    "deposit_date": dep_date,
                    "organism": organism,
                    "rcsb_url": f"https://www.rcsb.org/structure/{pdb_id}"
                }
                with open(cache_meta, "w", encoding="utf-8") as f:
                    json.dump(metadata, f, indent=2)

        pdb_text = None
        if cache_pdb.exists():
            try:
                with open(cache_pdb, "r", encoding="utf-8") as f:
                    pdb_text = f.read()
            except Exception:
                pass

        if not pdb_text:
            download_url = f"{RCSB_FILE_URL}/{pdb_id}.pdb"
            pdb_text = _get_text(download_url)
            if not pdb_text or "ATOM" not in pdb_text:
                # Fallback to .cif download and convert via gemmi
                cif_url = f"{RCSB_FILE_URL}/{pdb_id}.cif"
                cif_text = _get_text(cif_url)
                if cif_text and ("_atom_site" in cif_text or "data_" in cif_text):
                    try:
                        import gemmi
                        st = gemmi.read_structure_string(cif_text, format=gemmi.CoorFormat.Detect)
                        pdb_text = st.make_pdb_string()
                    except Exception as e:
                        print(f"[RCSB FETCH] CIF fallback conversion failed for {pdb_id}: {e}")

            if pdb_text and "ATOM" in pdb_text:
                with open(cache_pdb, "w", encoding="utf-8") as f:
                    f.write(pdb_text)

        if not pdb_text:
            return None

        metadata["pdb_content"] = pdb_text
        return metadata

    @staticmethod
    def search_rcsb_keywords(keyword: str, max_results: int = 6) -> List[Dict[str, Any]]:
        """Search RCSB PDB for target macromolecule entries matching keyword."""
        query_payload = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": keyword
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": max_results
                }
            }
        }
        try:
            req = urllib.request.Request(
                RCSB_SEARCH_URL,
                data=json.dumps(query_payload).encode("utf-8"),
                headers={"Content-Type": "application/json", "User-Agent": "Bindora/1.0"}
            )
            with urllib.request.urlopen(req, timeout=15) as resp:
                search_data = json.loads(resp.read().decode("utf-8"))
                result_ids = [item["identifier"] for item in search_data.get("result_set", [])]
        except Exception:
            return []

        results = []
        for pdb_id in result_ids:
            meta = StructureFetcher.fetch_rcsb_pdb(pdb_id)
            if meta:
                # Return summary without full pdb content
                summary = {k: v for k, v in meta.items() if k != "pdb_content"}
                results.append(summary)
        return results

    @staticmethod
    def search_uniprot(query: str) -> Optional[Dict[str, Any]]:
        """Fetch protein functional and pathway annotations from UniProt."""
        url = f"{UNIPROT_BASE_URL}/search?query={urllib.parse.quote(query)}&format=json&size=1"
        data = _get_json(url)
        if not data or "results" not in data or not data["results"]:
            return None

        entry = data["results"][0]
        accession = entry.get("primaryAccession", "")
        gene_name = ""
        try:
            gene_name = entry["genes"][0]["geneName"]["value"]
        except Exception:
            pass

        rec_name = ""
        try:
            rec_name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
        except Exception:
            rec_name = gene_name or query

        organism = entry.get("organism", {}).get("scientificName", "")
        
        # Extract function comment
        function_desc = ""
        catalytic_activity = ""
        subcellular_loc = []
        for comment in entry.get("comments", []):
            ctype = comment.get("commentType")
            if ctype == "FUNCTION" and not function_desc:
                function_desc = comment.get("texts", [{}])[0].get("value", "")
            elif ctype == "CATALYTIC ACTIVITY" and not catalytic_activity:
                catalytic_activity = comment.get("reaction", {}).get("name", "")
            elif ctype == "SUBCELLULAR LOCATION":
                for loc in comment.get("subcellularLocations", []):
                    subcellular_loc.append(loc.get("location", {}).get("value", ""))

        return {
            "accession": accession,
            "protein_name": rec_name,
            "gene_name": gene_name,
            "organism": organism,
            "function": function_desc,
            "catalytic_activity": catalytic_activity,
            "subcellular_location": subcellular_loc,
            "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}"
        }
