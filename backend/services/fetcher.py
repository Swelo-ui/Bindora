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
    ALPHAFOLD_BASE_URL,
    ALPHAFOLD_FILE_URL,
)

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
    "Accept": "application/json, text/plain, */*"
}

KNOWN_DRUG_CIDS: Dict[str, int] = {
    "indapamide": 3702,
    "imatinib": 5291,
    "aspirin": 2244,
    "paracetamol": 1983,
    "acetaminophen": 1983,
    "erlotinib": 176870,
    "gefitinib": 123631,
    "staurosporine": 44259,
    "indinavir": 5362440,
    "caffeine": 2519,
    "metformin": 4091,
    "losartan": 3961,
    "atorvastatin": 60823,
    "ibuprofen": 3672,
    "remdesivir": 121304016,
    "benzene": 241,
    "etoricoxib": 123616,
    "celecoxib": 2662,
    "dasatinib": 3062316,
    "sorafenib": 216239,
    "sunitinib": 5329102,
}

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

def _enrich_mol_properties_and_sdf(
    smiles: str,
    name: str,
    cid: Any = None,
    source: str = "Online Repository",
    raw_sdf: Optional[str] = None,
    iupac_name: str = ""
) -> Optional[Dict[str, Any]]:
    """Generate 3D coordinates and compute physicochemical descriptors locally with RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, AllChem
    except ImportError:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        if raw_sdf:
            try:
                mol = Chem.MolFromMolBlock(raw_sdf)
            except Exception:
                mol = None
        if not mol:
            return None

    sdf_content = raw_sdf
    is_3d = False
    if sdf_content and "CONECT" in sdf_content and len(sdf_content) > 100:
        is_3d = True
    else:
        try:
            mol_h = Chem.AddHs(mol)
            embed_res = AllChem.EmbedMolecule(mol_h, randomSeed=42)
            if embed_res == 0:
                AllChem.MMFFOptimizeMolecule(mol_h, maxIters=200)
                sdf_content = Chem.MolToMolBlock(mol_h)
                is_3d = True
            else:
                sdf_content = Chem.MolToMolBlock(mol)
        except Exception:
            sdf_content = Chem.MolToMolBlock(mol)

    mol_no_h = Chem.RemoveHs(mol)
    try:
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol_no_h)
    except Exception:
        formula = ""
    try:
        weight = round(float(Descriptors.MolWt(mol_no_h)), 2)
    except Exception:
        weight = 0.0
    try:
        xlogp = round(float(Descriptors.MolLogP(mol_no_h)), 2)
    except Exception:
        xlogp = None
    try:
        tpsa = round(float(Descriptors.TPSA(mol_no_h)), 2)
    except Exception:
        tpsa = None
    try:
        hbd = int(Lipinski.NumHDonors(mol_no_h))
    except Exception:
        hbd = 0
    try:
        hba = int(Lipinski.NumHAcceptors(mol_no_h))
    except Exception:
        hba = 0
    try:
        rotb = int(Lipinski.NumRotatableBonds(mol_no_h))
    except Exception:
        rotb = 0

    norm_name = name.strip()
    canon_smiles = Chem.MolToSmiles(mol_no_h, isomericSmiles=True)

    clean_cid = cid
    if clean_cid is None:
        clean_cid = KNOWN_DRUG_CIDS.get(norm_name.lower())
    if clean_cid is None:
        clean_cid = norm_name.capitalize()

    if str(clean_cid).isdigit():
        comp_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{clean_cid}"
    elif str(clean_cid).startswith("CHEMBL"):
        comp_url = f"https://www.ebi.ac.uk/chembl/compound_report_card/{clean_cid}/"
    else:
        comp_url = f"https://pubchem.ncbi.nlm.nih.gov/#query={urllib.parse.quote(norm_name)}"

    return {
        "source": source,
        "cid": clean_cid,
        "name": norm_name.capitalize(),
        "iupac_name": iupac_name or "",
        "formula": formula,
        "weight": weight,
        "smiles": canon_smiles,
        "xlogp": xlogp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rotb": rotb,
        "sdf": sdf_content,
        "is_3d": is_3d,
        "url": comp_url
    }

class StructureFetcher:
    """Service to search and retrieve chemical and macromolecular structures from public repositories."""

    @staticmethod
    def search_pubchem(name: str) -> Optional[Dict[str, Any]]:
        """Search compound by name with multi-tier fallback (PubChem -> NCI CIR -> ChEMBL -> Local RDKit)."""
        name = name.strip()
        if not name:
            return None

        # Direct CID lookup if user typed pure digits
        if name.isdigit():
            return StructureFetcher.fetch_pubchem_by_cid(name)

        cache_file = CACHE_DIR / f"pubchem_{urllib.parse.quote_plus(name.lower())}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    if cached and (cached.get("smiles") or cached.get("sdf")):
                        return cached
            except Exception:
                pass

        # Tier 1: PubChem PUG REST
        try:
            prop_url = (
                f"{PUBCHEM_BASE_URL}/compound/name/{urllib.parse.quote(name)}/property/"
                f"MolecularFormula,MolecularWeight,CanonicalSMILES,ConnectivitySMILES,IUPACName,"
                f"XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
            )
            data = _get_json(prop_url, timeout=5)
            if data and "PropertyTable" in data and data["PropertyTable"].get("Properties"):
                props = data["PropertyTable"]["Properties"][0]
                cid = props.get("CID")
                smiles = props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")

                sdf_3d_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/SDF?record_type=3d"
                sdf_content = _get_text(sdf_3d_url, timeout=6)
                is_3d = True
                if not sdf_content or "CONECT" not in sdf_content:
                    sdf_2d_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/SDF"
                    sdf_content = _get_text(sdf_2d_url, timeout=6)
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

                if not result["sdf"] and result["smiles"]:
                    enriched = _enrich_mol_properties_and_sdf(result["smiles"], name, cid=cid, source="PubChem")
                    if enriched:
                        result["sdf"] = enriched["sdf"]
                        result["is_3d"] = enriched["is_3d"]

                try:
                    with open(cache_file, "w", encoding="utf-8") as f:
                        json.dump(result, f, indent=2)
                except Exception:
                    pass
                return result
        except Exception:
            pass

        # Tier 2: NCI Chemical Identifier Resolver (Cactus CIR)
        try:
            cir_smi_url = f"https://cactus.nci.nih.gov/chemical/structure/{urllib.parse.quote(name)}/smiles"
            cir_smiles = _get_text(cir_smi_url, timeout=6)
            if cir_smiles and not cir_smiles.startswith("<!DOCTYPE") and "Page not found" not in cir_smiles:
                cir_smiles = cir_smiles.strip().splitlines()[0]
                cir_sdf_url = f"https://cactus.nci.nih.gov/chemical/structure/{urllib.parse.quote(name)}/sdf"
                cir_sdf = _get_text(cir_sdf_url, timeout=6)
                cir_iupac_url = f"https://cactus.nci.nih.gov/chemical/structure/{urllib.parse.quote(name)}/iupac_name"
                cir_iupac = _get_text(cir_iupac_url, timeout=4) or ""
                cir_iupac = cir_iupac.strip() if not cir_iupac.startswith("<!DOCTYPE") else ""

                cid = KNOWN_DRUG_CIDS.get(name.lower(), f"CIR-{name.capitalize()}")
                enriched = _enrich_mol_properties_and_sdf(
                    cir_smiles,
                    name,
                    cid=cid,
                    source="PubChem (NCI CIR fallback)",
                    raw_sdf=cir_sdf,
                    iupac_name=cir_iupac
                )
                if enriched:
                    try:
                        with open(cache_file, "w", encoding="utf-8") as f:
                            json.dump(enriched, f, indent=2)
                    except Exception:
                        pass
                    return enriched
        except Exception:
            pass

        # Tier 3: ChEMBL REST API
        try:
            chembl_url = f"{CHEMBL_BASE_URL}/molecule/search?q={urllib.parse.quote(name)}&format=json"
            chembl_data = _get_json(chembl_url, timeout=6)
            if chembl_data and chembl_data.get("molecules"):
                best_mol = None
                for m in chembl_data["molecules"]:
                    struct = m.get("molecule_structures")
                    if struct and struct.get("canonical_smiles"):
                        pref = (m.get("pref_name") or "").lower()
                        if pref == name.lower():
                            best_mol = m
                            break
                        if best_mol is None:
                            best_mol = m
                if best_mol:
                    chembl_smiles = best_mol["molecule_structures"]["canonical_smiles"]
                    chembl_id = best_mol.get("molecule_chembl_id")
                    cid = KNOWN_DRUG_CIDS.get(name.lower(), chembl_id)
                    enriched = _enrich_mol_properties_and_sdf(
                        chembl_smiles,
                        best_mol.get("pref_name") or name,
                        cid=cid,
                        source="PubChem (ChEMBL fallback)"
                    )
                    if enriched:
                        try:
                            with open(cache_file, "w", encoding="utf-8") as f:
                                json.dump(enriched, f, indent=2)
                        except Exception:
                            pass
                        return enriched
        except Exception:
            pass

        # Tier 4: PubChem Autocomplete suggestion resolution
        try:
            auto_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{urllib.parse.quote(name)}/json?limit=3"
            auto_data = _get_json(auto_url, timeout=5)
            terms = auto_data.get("dictionary_terms", {}).get("compound", []) if auto_data else []
            for term in terms:
                if term.lower() != name.lower():
                    cir_url = f"https://cactus.nci.nih.gov/chemical/structure/{urllib.parse.quote(term)}/smiles"
                    cir_smi = _get_text(cir_url, timeout=5)
                    if cir_smi and not cir_smi.startswith("<!DOCTYPE") and "Page not found" not in cir_smi:
                        cir_smi = cir_smi.strip().splitlines()[0]
                        cid = KNOWN_DRUG_CIDS.get(term.lower(), f"CIR-{term.capitalize()}")
                        enriched = _enrich_mol_properties_and_sdf(
                            cir_smi,
                            term,
                            cid=cid,
                            source="PubChem (Autocomplete fallback)"
                        )
                        if enriched:
                            try:
                                with open(cache_file, "w", encoding="utf-8") as f:
                                    json.dump(enriched, f, indent=2)
                            except Exception:
                                pass
                            return enriched
        except Exception:
            pass

        return None

    @staticmethod
    def fetch_pubchem_by_cid(cid: str) -> Optional[Dict[str, Any]]:
        """Fetch compound metadata, title, and properties by PubChem CID with PUG View & ChEMBL fallbacks."""
        cid = str(cid).strip()
        if not cid:
            return None

        cache_file = CACHE_DIR / f"pubchem_cid_{cid}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    if cached and (cached.get("smiles") or cached.get("sdf")):
                        return cached
            except Exception:
                pass

        # Tier 1: PubChem PUG REST
        if cid.isdigit():
            try:
                prop_url = (
                    f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/"
                    f"Title,MolecularFormula,MolecularWeight,CanonicalSMILES,ConnectivitySMILES,IUPACName,"
                    f"XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
                )
                data = _get_json(prop_url, timeout=5)
                if data and "PropertyTable" in data and data["PropertyTable"].get("Properties"):
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
            except Exception:
                pass

        # Tier 2: PubChem PUG View (works during PUG REST 503 ServerBusy)
        if cid.isdigit():
            try:
                view_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
                view_data = _get_json(view_url, timeout=6)
                if view_data and "Record" in view_data:
                    record = view_data["Record"]
                    rec_title = record.get("RecordTitle") or f"Compound #{cid}"

                    def _extract_smiles_from_sections(sections):
                        for s in sections:
                            h = s.get("TOCHeading", "")
                            if h in ("Canonical SMILES", "Isomeric SMILES", "SMILES"):
                                for info in s.get("Information", []):
                                    val_str = info.get("Value", {}).get("StringWithMarkup", [{}])[0].get("String")
                                    if val_str:
                                        return val_str
                            sub = s.get("Section", [])
                            if sub:
                                found = _extract_smiles_from_sections(sub)
                                if found:
                                    return found
                        return None

                    smiles_found = _extract_smiles_from_sections(record.get("Section", []))
                    if smiles_found:
                        enriched = _enrich_mol_properties_and_sdf(
                            smiles_found,
                            rec_title,
                            cid=int(cid),
                            source="PubChem (PUG View)"
                        )
                        if enriched:
                            try:
                                with open(cache_file, "w", encoding="utf-8") as f:
                                    json.dump(enriched, f, indent=2)
                            except Exception:
                                pass
                            return enriched
            except Exception:
                pass

        # Tier 3: ChEMBL if CID is a CHEMBL ID
        if cid.startswith("CHEMBL"):
            try:
                ch_url = f"{CHEMBL_BASE_URL}/molecule/{cid}?format=json"
                ch_data = _get_json(ch_url, timeout=6)
                if ch_data and ch_data.get("molecule_structures"):
                    smi = ch_data["molecule_structures"].get("canonical_smiles")
                    if smi:
                        enriched = _enrich_mol_properties_and_sdf(
                            smi,
                            ch_data.get("pref_name") or cid,
                            cid=cid,
                            source="PubChem (ChEMBL)"
                        )
                        if enriched:
                            try:
                                with open(cache_file, "w", encoding="utf-8") as f:
                                    json.dump(enriched, f, indent=2)
                            except Exception:
                                pass
                            return enriched
            except Exception:
                pass

        return None

    @staticmethod
    def search_similar_compounds(smiles: str, threshold: int = 85, max_records: int = 5) -> List[Dict[str, Any]]:
        """Search structurally similar compounds with PubChem fastsimilarity & ChEMBL fallback."""
        if not smiles or not smiles.strip():
            return []

        cache_key = f"pubchem_sim_{urllib.parse.quote_plus(smiles.strip())}_{threshold}_{max_records}.json"
        cache_file = CACHE_DIR / cache_key
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                pass

        results = []
        # Tier 1: PubChem fastsimilarity_2d
        try:
            url = f"{PUBCHEM_BASE_URL}/compound/fastsimilarity_2d/smiles/{urllib.parse.quote(smiles.strip())}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}"
            data = _get_json(url, timeout=8)
            cids = data.get("IdentifierList", {}).get("CID", []) if data else []
            if cids:
                cids_str = ",".join(map(str, cids[:max_records]))
                prop_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cids_str}/property/Title,MolecularWeight,ConnectivitySMILES,CanonicalSMILES/JSON"
                prop_data = _get_json(prop_url, timeout=8)
                if prop_data and "PropertyTable" in prop_data:
                    for item in prop_data["PropertyTable"].get("Properties", []):
                        item_smi = item.get("ConnectivitySMILES") or item.get("CanonicalSMILES") or ""
                        comp_name = item.get("Title") or f"Compound #{item.get('CID')}"
                        comp_weight = float(item.get("MolecularWeight", 0.0))
                        results.append({
                            "cid": item.get("CID"),
                            "name": comp_name,
                            "title": comp_name,
                            "weight": comp_weight,
                            "molecular_weight": comp_weight,
                            "smiles": item_smi,
                            "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{item.get('CID')}"
                        })
        except Exception:
            pass

        # Tier 2: ChEMBL similarity fallback
        if not results:
            try:
                ch_sim = min(max(int(threshold), 60), 95)
                ch_url = f"{CHEMBL_BASE_URL}/similarity/{urllib.parse.quote(smiles.strip())}/{ch_sim}?format=json&limit={max_records}"
                ch_data = _get_json(ch_url, timeout=8)
                if ch_data and ch_data.get("molecules"):
                    for m in ch_data["molecules"][:max_records]:
                        struct = m.get("molecule_structures") or {}
                        m_smi = struct.get("canonical_smiles") or ""
                        m_id = m.get("molecule_chembl_id")
                        m_name = m.get("pref_name") or f"ChEMBL {m_id}"
                        props = m.get("molecule_properties") or {}
                        try:
                            m_wt = float(props.get("full_mwt") or 0.0)
                        except (ValueError, TypeError):
                            m_wt = 0.0
                        if m_smi:
                            results.append({
                                "cid": m_id,
                                "name": m_name,
                                "title": m_name,
                                "weight": m_wt,
                                "molecular_weight": m_wt,
                                "smiles": m_smi,
                                "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{m_id}/"
                            })
            except Exception:
                pass

        if results:
            try:
                with open(cache_file, "w", encoding="utf-8") as f:
                    json.dump(results, f, indent=2)
            except Exception:
                pass

        return results

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
        
        # Extract detected bound ligands from HETATM records
        lig_set = set()
        solvent_ions = {"HOH", "WAT", "DOD", "TIP", "NA", "CL", "MG", "ZN", "CA", "MN", "K", "FE", "SO4", "PO4", "ACT", "EDO", "DMS", "GOL", "PEG", "MPD", "BME", "TRS", "FMT"}
        for line in pdb_text.splitlines():
            if line.startswith("HETATM"):
                rname = line[17:20].strip()
                if rname and rname not in solvent_ions and len(rname) >= 2:
                    lig_set.add(rname)
        metadata["ligand_names"] = sorted(list(lig_set))
        metadata["has_ligand"] = bool(lig_set)

        if "uniprot" not in metadata or not metadata["uniprot"]:
            metadata["uniprot"] = StructureFetcher.fetch_uniprot_by_pdb_id(pdb_id)
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
    def search_rcsb_by_uniprot(uniprot_acc: str, max_results: int = 4) -> List[Dict[str, Any]]:
        """Search RCSB PDB for experimental structures matching a specific UniProt accession."""
        uniprot_acc = uniprot_acc.strip().upper()
        if not uniprot_acc:
            return []

        cache_file = CACHE_DIR / f"rcsb_uniprot_{uniprot_acc}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    if cached and isinstance(cached, list):
                        return cached[:max_results]
            except Exception:
                pass

        query_payload = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_acc
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": max_results * 2
                },
                "sort": [
                    {
                        "sort_by": "rcsb_entry_info.resolution_combined",
                        "direction": "asc"
                    }
                ]
            }
        }
        try:
            req = urllib.request.Request(
                RCSB_SEARCH_URL,
                data=json.dumps(query_payload).encode("utf-8"),
                headers={"Content-Type": "application/json", "User-Agent": "Bindora/1.0"}
            )
            with urllib.request.urlopen(req, timeout=12) as resp:
                search_data = json.loads(resp.read().decode("utf-8"))
                result_ids = [item["identifier"] for item in search_data.get("result_set", [])]
        except Exception:
            result_ids = []

        results = []
        for pdb_id in result_ids:
            meta = StructureFetcher.fetch_rcsb_pdb(pdb_id)
            if meta:
                summary = {k: v for k, v in meta.items() if k != "pdb_content"}
                summary["has_ligand"] = bool(meta.get("ligand_names"))
                results.append(summary)
            if len(results) >= max_results:
                break

        if results:
            try:
                with open(cache_file, "w", encoding="utf-8") as f:
                    json.dump(results, f, indent=2)
            except Exception:
                pass
        return results

    @staticmethod
    def _parse_uniprot_entry(entry: Dict[str, Any]) -> Dict[str, Any]:
        """Normalize raw UniProtKB JSON entry into Bindora standard annotation dict."""
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
            rec_name = gene_name or accession

        organism = entry.get("organism", {}).get("scientificName", "")

        function_desc = ""
        catalytic_activity = ""
        pathway_desc = ""
        tissue_spec = ""
        subcellular_loc = []
        for comment in entry.get("comments", []):
            ctype = comment.get("commentType")
            if ctype == "FUNCTION" and not function_desc:
                function_desc = comment.get("texts", [{}])[0].get("value", "")
            elif ctype == "CATALYTIC ACTIVITY" and not catalytic_activity:
                catalytic_activity = comment.get("reaction", {}).get("name", "")
            elif ctype == "PATHWAY" and not pathway_desc:
                pathway_desc = comment.get("texts", [{}])[0].get("value", "")
            elif ctype == "TISSUE SPECIFICITY" and not tissue_spec:
                tissue_spec = comment.get("texts", [{}])[0].get("value", "")
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
            "pathway": pathway_desc,
            "tissue_specificity": tissue_spec,
            "subcellular_location": subcellular_loc,
            "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}"
        }

    @staticmethod
    def fetch_uniprot_by_pdb_id(pdb_id: str) -> Optional[Dict[str, Any]]:
        """Fetch protein functional and pathway annotations from UniProt using PDB cross-reference."""
        pdb_id = pdb_id.strip().upper()
        if len(pdb_id) != 4:
            return None
        cache_file = CACHE_DIR / f"uniprot_pdb_{pdb_id}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                pass

        url = f"{UNIPROT_BASE_URL}/search?query=xref:pdb-{pdb_id}&format=json&size=1"
        data = _get_json(url)
        if not data or "results" not in data or not data["results"]:
            url_lower = f"{UNIPROT_BASE_URL}/search?query=xref:pdb-{pdb_id.lower()}&format=json&size=1"
            data = _get_json(url_lower)

        if not data or "results" not in data or not data["results"]:
            return None

        result = StructureFetcher._parse_uniprot_entry(data["results"][0])
        try:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2)
        except Exception:
            pass
        return result

    @staticmethod
    def search_uniprot(query: str) -> Optional[Dict[str, Any]]:
        """Fetch protein functional and pathway annotations from UniProt, prioritizing reviewed Swiss-Prot records."""
        q_clean = query.strip()
        if not q_clean:
            return None

        # Direct 4-character PDB code
        if len(q_clean) == 4 and q_clean.isalnum():
            pdb_res = StructureFetcher.fetch_uniprot_by_pdb_id(q_clean)
            if pdb_res:
                return pdb_res

        # Direct UniProt Accession (e.g. P00533, P24941)
        if (len(q_clean) in (6, 10)) and (q_clean[0].upper() in "OPQ" or q_clean.isalnum()):
            direct_data = _get_json(f"{UNIPROT_BASE_URL}/{q_clean.upper()}.json", timeout=6)
            if direct_data and "primaryAccession" in direct_data:
                return StructureFetcher._parse_uniprot_entry(direct_data)

        # Keyword Search: prioritize reviewed (Swiss-Prot) and organism if specified
        is_human = any(w in q_clean.lower() for w in ("human", "homo", "sapiens", "hsa"))
        stopwords = {"human", "homo", "sapiens", "protein", "receptor", "kinase", "inhibitor", "structure", "complex"}
        tokens = [w for w in q_clean.split() if w.lower() not in stopwords]
        core_query = " ".join(tokens) if tokens else q_clean

        candidate_queries = []
        if is_human:
            candidate_queries.append(f"(gene:{urllib.parse.quote(core_query)} OR {urllib.parse.quote(core_query)}) AND reviewed:true AND organism_id:9606")
            candidate_queries.append(f"({urllib.parse.quote(q_clean)}) AND reviewed:true AND organism_id:9606")
        candidate_queries.append(f"(gene:{urllib.parse.quote(core_query)} OR {urllib.parse.quote(core_query)}) AND reviewed:true")
        candidate_queries.append(f"({urllib.parse.quote(q_clean)}) AND reviewed:true")
        candidate_queries.append(f"{urllib.parse.quote(q_clean)}")

        data = None
        for q_str in candidate_queries:
            url = f"{UNIPROT_BASE_URL}/search?query={q_str}&format=json&size=5"
            data = _get_json(url, timeout=7)
            if data and "results" in data and data["results"]:
                break

        if not data or "results" not in data or not data["results"]:
            return None

        results = data["results"]
        best_entry = results[0]
        # Prefer exact gene or organism match if multiple
        for item in results:
            genes = [g.get("geneName", {}).get("value", "").lower() for g in item.get("genes", [])]
            if core_query.lower() in genes:
                best_entry = item
                break

        return StructureFetcher._parse_uniprot_entry(best_entry)

    @staticmethod
    def fetch_alphafold_structure(uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Fetch AlphaFold predicted structure (.pdb) and confidence metrics from EBI AlphaFold Database.
        """
        uniprot_id = uniprot_id.strip().upper()
        if not uniprot_id:
            return None

        clean_acc = uniprot_id.replace("AF-", "").split("-")[0]
        cache_pdb = CACHE_DIR / f"AF-{clean_acc}_predicted.pdb"
        cache_meta = CACHE_DIR / f"AF-{clean_acc}_meta.json"

        if cache_meta.exists() and cache_pdb.exists():
            try:
                with open(cache_meta, "r", encoding="utf-8") as f:
                    meta = json.load(f)
                with open(cache_pdb, "r", encoding="utf-8") as f:
                    pdb_text = f.read()
                if pdb_text and len(pdb_text) > 100:
                    meta["pdb_content"] = pdb_text
                    return meta
            except Exception:
                pass

        api_url = f"{ALPHAFOLD_BASE_URL}/{clean_acc}"
        data = _get_json(api_url, timeout=12)
        if not data or not isinstance(data, list) or len(data) == 0:
            return None

        entry = data[0]
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            return None

        pdb_content = _get_text(pdb_url, timeout=20)
        if not pdb_content or len(pdb_content) < 100:
            return None

        try:
            cache_pdb.write_text(pdb_content, encoding="utf-8")
        except Exception:
            pass

        meta = {
            "pdb_id": f"AF-{clean_acc}",
            "uniprot_id": clean_acc,
            "title": f"AlphaFold Predicted: {entry.get('uniprotDescription', clean_acc)} ({entry.get('organismScientificName', 'Homo sapiens')})",
            "gene": entry.get("gene", ""),
            "organism": entry.get("organismScientificName", "Homo sapiens"),
            "sequence_length": entry.get("sequenceEnd", 0) - entry.get("sequenceStart", 1) + 1,
            "plddt": round(float(entry.get("globalMetricValue", 0)), 1),
            "fraction_confident": round(float(entry.get("fractionPlddtConfident", 0)) * 100, 1),
            "fraction_very_high": round(float(entry.get("fractionPlddtVeryHigh", 0)) * 100, 1),
            "source": "AlphaFold DB",
            "is_alphafold": True,
            "pdb_url": pdb_url,
            "pae_image_url": entry.get("paeImageUrl", ""),
        }

        try:
            with open(cache_meta, "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)
        except Exception:
            pass

        meta["pdb_content"] = pdb_content
        return meta

    @staticmethod
    def search_unified_target(query: str, max_rcsb_results: int = 4) -> Dict[str, Any]:
        """
        Unified target resolver supporting:
        1. 4-letter PDB ID (e.g. 1M17) -> Returns crystal PDB + linked AlphaFold model
        2. UniProt Accession (e.g. P00533) -> Returns AlphaFold model + top crystal PDBs
        3. Protein/gene name (e.g. 'Human EGFR', 'CDK2') -> Fetches both with intelligent recommendation
        """
        query_clean = query.strip()
        if not query_clean:
            return {"error": "Empty query"}

        target_name = ""
        gene_name = ""
        organism = "Homo sapiens"
        uniprot_acc = ""
        rcsb_entries: List[Dict[str, Any]] = []
        af_model: Optional[Dict[str, Any]] = None

        # Case 1: 4-character PDB code
        if len(query_clean) == 4 and query_clean.isalnum():
            pdb_id = query_clean.upper()
            direct_pdb = StructureFetcher.fetch_rcsb_pdb(pdb_id)
            if direct_pdb:
                direct_summary = {k: v for k, v in direct_pdb.items() if k != "pdb_content"}
                direct_summary["has_ligand"] = bool(direct_pdb.get("ligand_names"))
                rcsb_entries.append(direct_summary)
                target_name = direct_pdb.get("title", pdb_id)

                uniprot_data = StructureFetcher.fetch_uniprot_by_pdb_id(pdb_id)
                if uniprot_data:
                    uniprot_acc = uniprot_data.get("accession", "")
                    gene_name = uniprot_data.get("gene_name", "")
                    organism = uniprot_data.get("organism", organism)

        # Case 2: UniProt Accession ID (e.g., P00533, P24941, Q9ULV8)
        elif len(query_clean) in (6, 10) and (query_clean[0].upper() in "OPQ" or query_clean.isalnum()):
            uniprot_info = StructureFetcher.search_uniprot(query_clean)
            if uniprot_info and uniprot_info.get("accession"):
                uniprot_acc = uniprot_info.get("accession")
                target_name = uniprot_info.get("protein_name", "")
                gene_name = uniprot_info.get("gene_name", "")
                organism = uniprot_info.get("organism", organism)

        # Case 3: Keyword / Protein name search
        if not uniprot_acc:
            uniprot_info = StructureFetcher.search_uniprot(query_clean)
            if uniprot_info and uniprot_info.get("accession"):
                uniprot_acc = uniprot_info.get("accession")
                target_name = uniprot_info.get("protein_name", query_clean)
                gene_name = uniprot_info.get("gene_name", "")
                organism = uniprot_info.get("organism", organism)

        # If we have a UniProt accession, fetch AlphaFold model
        if uniprot_acc:
            af_raw = StructureFetcher.fetch_alphafold_structure(uniprot_acc)
            if af_raw:
                af_model = {k: v for k, v in af_raw.items() if k != "pdb_content"}

        # If RCSB entries not already populated, search RCSB by UniProt accession first, then keywords
        if not rcsb_entries:
            if uniprot_acc:
                by_up = StructureFetcher.search_rcsb_by_uniprot(uniprot_acc, max_results=max_rcsb_results)
                if by_up:
                    rcsb_entries.extend(by_up)

            if not rcsb_entries:
                search_terms = [query_clean]
                if gene_name and gene_name.lower() not in query_clean.lower():
                    search_terms.insert(0, f"{gene_name} {organism}")

                for term in search_terms:
                    found_rcsb = StructureFetcher.search_rcsb_keywords(term, max_results=max_rcsb_results)
                    if found_rcsb:
                        for item in found_rcsb:
                            item["has_ligand"] = bool(item.get("ligand_names"))
                            rcsb_entries.append(item)
                        break

        best_rcsb = None
        for entry in rcsb_entries:
            res = entry.get("resolution")
            if res and float(res) <= 2.8 and entry.get("has_ligand"):
                best_rcsb = entry
                break
        if not best_rcsb and rcsb_entries:
            best_rcsb = rcsb_entries[0]

        if best_rcsb and best_rcsb.get("has_ligand") and best_rcsb.get("resolution") and float(best_rcsb.get("resolution")) <= 2.8:
            recommendation = {
                "source": "rcsb",
                "recommended_id": best_rcsb.get("pdb_id"),
                "badge": "Recommended for Pocket Docking",
                "reason": f"High-resolution ({best_rcsb.get('resolution')} Å) co-crystal with bound active pocket."
            }
        elif af_model and af_model.get("plddt", 0) >= 70:
            recommendation = {
                "source": "alphafold",
                "recommended_id": af_model.get("pdb_id"),
                "badge": "Recommended (Full Length)",
                "reason": f"Complete full-length predicted chain ({af_model.get('sequence_length')} aa, pLDDT: {af_model.get('plddt')}) with zero missing loops."
            }
        elif best_rcsb:
            recommendation = {
                "source": "rcsb",
                "recommended_id": best_rcsb.get("pdb_id"),
                "badge": "Experimental Structure",
                "reason": f"Experimental crystal structure ({best_rcsb.get('resolution') or 'X-ray'} Å)."
            }
        else:
            recommendation = {
                "source": "none",
                "recommended_id": "",
                "badge": "",
                "reason": "No high-confidence structure found."
            }

        return {
            "query": query_clean,
            "target_name": target_name or query_clean,
            "gene_name": gene_name,
            "organism": organism,
            "uniprot_id": uniprot_acc,
            "alphafold": af_model,
            "rcsb_entries": rcsb_entries[:max_rcsb_results],
            "recommendation": recommendation
        }
