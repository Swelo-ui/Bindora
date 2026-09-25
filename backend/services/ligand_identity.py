"""
backend/services/ligand_identity.py
=====================================
Bindora Dock - Canonical Ligand Identity Resolution Service

Purpose:
    Determine whether two ligands are chemically identical using structure-first
    molecular graph comparison. This service is the authoritative source for
    docking experiment classification (Native Redocking vs Cross-Docking).

Priority order (per scientific spec):
    Tier 1: InChIKey comparison (connectivity layer, ignoring stereochemistry)
    Tier 2: Canonical SMILES comparison (isomericSmiles=False for atom-graph equality)
    Tier 3: Molecular formula + heavy atom count agreement
    Tier 4: Synonym / known-alias lookup against RCSB Chemical Component Dictionary
    Tier 5: Name substring heuristic (LAST RESORT ONLY, must set confidence=low)

Do NOT classify as identical solely from name similarity when structure fails.
Do NOT hard-code compound-specific lookup tables.
"""

import json
import hashlib
import urllib.request
import urllib.parse
from typing import Any, Dict, List, Optional
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import inchi as rdkit_inchi

# Cache directory (filled in at runtime from config)
try:
    from backend.config import CACHE_DIR
except Exception:
    CACHE_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "cache"

RCSB_CHEMCOMP_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{code}"
_HEADERS = {"User-Agent": "Bindora-Research-Tool/1.0 (academic)"}

# ─────────────────────────────────────────────────────────────────────────────
# InChIKey helpers
# ─────────────────────────────────────────────────────────────────────────────

def _inchikey_from_smiles(smiles: str) -> Optional[str]:
    """Return InChIKey for a SMILES string, or None if conversion fails."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return rdkit_inchi.MolToInchiKey(mol)
    except Exception:
        pass
    return None


def _inchikey_from_mol(mol: Chem.Mol) -> Optional[str]:
    """Return InChIKey for an RDKit Mol, or None."""
    try:
        return rdkit_inchi.MolToInchiKey(mol)
    except Exception:
        return None


def _canonical_smiles(smiles: str) -> Optional[str]:
    """Return RDKit canonical SMILES (isomericSmiles=False) or None."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        pass
    return None


# ─────────────────────────────────────────────────────────────────────────────
# RCSB Chemical Component lookup
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_rcsb_chemcomp(pdb_code: str) -> Optional[Dict[str, Any]]:
    """
    Fetch canonical SMILES, InChIKey, synonyms for a PDB 3-letter ligand code
    from the RCSB Chemical Component Dictionary REST API.

    Result is cached locally to avoid repeated network calls.
    Returns None on any network / parse error.
    """
    code = pdb_code.strip().upper()
    cache_path = CACHE_DIR / f"rcsb_chemcomp_{code}.json"
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass

    # Cache hit
    if cache_path.exists():
        try:
            with open(cache_path, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            pass

    url = RCSB_CHEMCOMP_URL.format(code=urllib.parse.quote(code))
    try:
        req = urllib.request.Request(url, headers=_HEADERS)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception:
        return None

    # Extract what we need
    desc = data.get("rcsb_chem_comp_descriptor", {})
    smiles_canonical = desc.get("SMILES") or desc.get("SMILES_stereo") or ""
    inchikey = desc.get("InChIKey", "")
    inchi = desc.get("InChI", "")

    synonyms: List[str] = []
    for s in data.get("rcsb_chem_comp_synonyms", []):
        name = s.get("name", "")
        if name:
            try:
                synonyms.append(name.encode("ascii", "ignore").decode("ascii").strip())
            except Exception:
                pass

    # Also add pdbx identifiers
    for ident in data.get("pdbx_chem_comp_identifier", []):
        identifier = ident.get("identifier", "")
        if identifier and ident.get("type") in ("SYSTEMATIC NAME", "SYNONYM"):
            try:
                synonyms.append(identifier.encode("ascii", "ignore").decode("ascii").strip())
            except Exception:
                pass

    chem_comp = data.get("chem_comp", {})
    name = chem_comp.get("name", "")
    formula = chem_comp.get("formula", "")

    result = {
        "pdb_code": code,
        "canonical_name": name,
        "formula": formula,
        "canonical_smiles": smiles_canonical,
        "inchikey": inchikey,
        "inchi": inchi,
        "synonyms": synonyms,
    }

    # Write cache
    try:
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=2, ensure_ascii=False)
    except Exception:
        pass

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Core identity comparison
# ─────────────────────────────────────────────────────────────────────────────

def compare_ligand_identity(
    docked_smiles: str,
    native_smiles: Optional[str] = None,
    native_pdb_code: Optional[str] = None,
    native_name: Optional[str] = None,
    docked_name: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compare two ligands and determine if they are the same chemical compound.

    Decision priority:
      1. InChIKey (connectivity layer) comparison between docked SMILES and:
         a. native_smiles (already resolved)
         b. RCSB chemcomp canonical SMILES for native_pdb_code (if code given)
      2. Canonical SMILES (isomericSmiles=False) equality
      3. Molecular formula + heavy-atom count agreement
      4. Synonym lookup: does docked_name appear in RCSB synonym list for native_pdb_code?
      5. Last resort: docked_name / native_name substring match (flags confidence=low)

    Returns dict with:
      is_same_molecule (bool)
      confidence       (str): "structure-confirmed" | "structure-probable" |
                              "formula-only" | "synonym-matched" | "name-heuristic" | "undetermined"
      method           (str): human-readable description of the tier that matched
      canonical_name   (str): resolved name for the native ligand
      canonical_smiles (str): RCSB canonical SMILES for the native ligand (if available)
      inchikey_native  (str): InChIKey of native ligand
      inchikey_docked  (str): InChIKey of docked ligand
      synonyms         (list[str]): known synonyms from RCSB
    """
    result: Dict[str, Any] = {
        "is_same_molecule": False,
        "confidence": "undetermined",
        "method": "No comparison performed",
        "canonical_name": native_name or native_pdb_code or "Unknown",
        "canonical_smiles": native_smiles or "",
        "inchikey_native": "",
        "inchikey_docked": "",
        "synonyms": [],
    }

    if not docked_smiles or not docked_smiles.strip():
        result["method"] = "Docked SMILES missing — cannot compare"
        return result

    # Get InChIKey for docked molecule
    ikey_docked = _inchikey_from_smiles(docked_smiles)
    result["inchikey_docked"] = ikey_docked or ""

    # ── Tier 0: Try to resolve native ligand from RCSB chemcomp ────────────
    rcsb_data: Optional[Dict[str, Any]] = None
    rcsb_smiles: Optional[str] = None
    rcsb_ikey: Optional[str] = None

    if native_pdb_code:
        rcsb_data = _fetch_rcsb_chemcomp(native_pdb_code)
        if rcsb_data:
            result["canonical_name"] = rcsb_data.get("canonical_name") or native_name or native_pdb_code
            result["canonical_smiles"] = rcsb_data.get("canonical_smiles") or native_smiles or ""
            result["synonyms"] = rcsb_data.get("synonyms", [])
            result["inchikey_native"] = rcsb_data.get("inchikey", "")
            rcsb_smiles = rcsb_data.get("canonical_smiles")
            rcsb_ikey = rcsb_data.get("inchikey")

    # Also compute InChIKey from native_smiles (from PDB parsing) if available
    ikey_native_from_smiles: Optional[str] = None
    if native_smiles:
        ikey_native_from_smiles = _inchikey_from_smiles(native_smiles)

    # ── Tier 1: InChIKey comparison (connectivity layer only) ──────────────
    if ikey_docked:
        # 1a: vs RCSB authoritative InChIKey
        if rcsb_ikey and rcsb_ikey.strip():
            # Compare first 14 chars (connectivity layer, ignoring stereo)
            if ikey_docked[:14] == rcsb_ikey[:14]:
                result["is_same_molecule"] = True
                result["confidence"] = "structure-confirmed"
                result["method"] = f"InChIKey connectivity match (RCSB CCD): docked={ikey_docked} == native={rcsb_ikey}"
                result["inchikey_native"] = rcsb_ikey
                return result

        # 1b: vs PDB-parsed native SMILES InChIKey
        if ikey_native_from_smiles:
            result["inchikey_native"] = ikey_native_from_smiles
            if ikey_docked[:14] == ikey_native_from_smiles[:14]:
                result["is_same_molecule"] = True
                result["confidence"] = "structure-confirmed"
                result["method"] = f"InChIKey connectivity match (PDB-parsed): docked={ikey_docked} == native={ikey_native_from_smiles}"
                return result

    # ── Tier 2: Canonical SMILES comparison ────────────────────────────────
    can_docked = _canonical_smiles(docked_smiles)
    if can_docked:
        # 2a: vs RCSB canonical SMILES
        if rcsb_smiles:
            can_rcsb = _canonical_smiles(rcsb_smiles)
            if can_rcsb and can_docked == can_rcsb:
                result["is_same_molecule"] = True
                result["confidence"] = "structure-confirmed"
                result["method"] = "Canonical SMILES equality (RCSB CCD)"
                return result

        # 2b: vs PDB-parsed native SMILES
        if native_smiles:
            can_native = _canonical_smiles(native_smiles)
            if can_native and can_docked == can_native:
                result["is_same_molecule"] = True
                result["confidence"] = "structure-confirmed"
                result["method"] = "Canonical SMILES equality (PDB-parsed)"
                return result

    # ── Tier 3: Molecular formula + heavy atom count ───────────────────────
    try:
        mol_docked = Chem.MolFromSmiles(docked_smiles)
        if mol_docked and rcsb_data and rcsb_data.get("formula"):
            formula_native = rcsb_data["formula"].replace(" ", "")
            # Get formula for docked molecule using RDKit
            from rdkit.Chem import rdMolDescriptors
            formula_docked = rdMolDescriptors.CalcMolFormula(mol_docked).replace(" ", "")
            if formula_native and formula_docked and formula_native == formula_docked:
                heavy_docked = mol_docked.GetNumHeavyAtoms()
                mol_native_rcsb = Chem.MolFromSmiles(rcsb_smiles) if rcsb_smiles else None
                heavy_native = mol_native_rcsb.GetNumHeavyAtoms() if mol_native_rcsb else 0
                if heavy_native > 0 and heavy_docked == heavy_native:
                    result["is_same_molecule"] = True
                    result["confidence"] = "formula-only"
                    result["method"] = f"Formula + heavy-atom count agreement: {formula_docked}"
                    return result
    except Exception:
        pass

    # ── Tier 4: RCSB synonym lookup ────────────────────────────────────────
    if docked_name and result["synonyms"]:
        docked_norm = docked_name.strip().lower()
        for syn in result["synonyms"]:
            syn_norm = syn.strip().lower()
            if syn_norm and (syn_norm == docked_norm or docked_norm in syn_norm or syn_norm in docked_norm):
                result["is_same_molecule"] = True
                result["confidence"] = "synonym-matched"
                result["method"] = f"RCSB synonym match: '{docked_name}' matches '{syn}'"
                return result

    # ── Tier 5: Name substring heuristic (last resort) ─────────────────────
    if docked_name and native_name:
        d_lower = docked_name.strip().lower()
        n_lower = native_name.strip().lower()
        if len(d_lower) >= 4 and len(n_lower) >= 4:
            if d_lower == n_lower or d_lower in n_lower or n_lower in d_lower:
                result["is_same_molecule"] = True
                result["confidence"] = "name-heuristic"
                result["method"] = f"Name heuristic (LOW CONFIDENCE): '{docked_name}' ~ '{native_name}'"
                return result

    result["method"] = "No structural or synonym match found — classified as different molecules"
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def resolve_pdb_ligand_identity(pdb_code: str) -> Dict[str, Any]:
    """
    Resolve the canonical chemical identity of a PDB 3-letter ligand code.

    Returns:
        {
          "pdb_ligand_code": str,
          "canonical_name": str,
          "synonyms": list[str],
          "canonical_smiles": str,
          "inchi_key": str,
          "formula": str,
          "identity_confidence": str,  # "structure-confirmed" | "name-only" | "unknown"
        }
    """
    code = pdb_code.strip().upper()
    rcsb = _fetch_rcsb_chemcomp(code)
    if not rcsb:
        return {
            "pdb_ligand_code": code,
            "canonical_name": code,
            "synonyms": [],
            "canonical_smiles": "",
            "inchi_key": "",
            "formula": "",
            "identity_confidence": "unknown",
        }

    ikey = rcsb.get("inchikey", "")
    confidence = "structure-confirmed" if (rcsb.get("canonical_smiles") or ikey) else "name-only"

    return {
        "pdb_ligand_code": code,
        "canonical_name": rcsb.get("canonical_name", code),
        "synonyms": rcsb.get("synonyms", []),
        "canonical_smiles": rcsb.get("canonical_smiles", ""),
        "inchi_key": ikey,
        "formula": rcsb.get("formula", ""),
        "identity_confidence": confidence,
    }
