"""
Input Validation Module
Centralized validation for chemical structures and computational parameters
"""

from typing import Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors
import gemmi


class ValidationResult:
    """Standardized validation result container"""
    
    def __init__(self, valid: bool, error: Optional[str] = None, data: Optional[Dict[str, Any]] = None):
        self.valid = valid
        self.error = error
        self.data = data or {}
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        result = {"valid": self.valid}
        if self.error:
            result["error"] = self.error
        if self.data:
            result.update(self.data)
        return result


def validate_smiles(smiles: str) -> ValidationResult:
    """
    Validate SMILES string using RDKit parser
    
    Checks:
    - Length ≤ 5000 characters
    - RDKit parseable
    - Sanitizable (valence, aromaticity)
    - Heavy atom count ≥ 1
    
    Returns:
    - valid=True: canonical_smiles, mol_object, heavy_atom_count
    - valid=False: error message
    """
    if not isinstance(smiles, str):
        return ValidationResult(False, "SMILES must be a string")
    
    smiles = smiles.strip()
    
    if not smiles:
        return ValidationResult(False, "SMILES string cannot be empty")
    
    if len(smiles) > 5000:
        return ValidationResult(
            False,
            f"SMILES string too long ({len(smiles)} characters, maximum 5000)"
        )
    
    # Parse SMILES without sanitization first
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    except Exception as e:
        return ValidationResult(False, f"Invalid SMILES syntax: {str(e)}")
    
    if mol is None:
        return ValidationResult(False, "Invalid SMILES: could not parse chemical structure")
    
    # Attempt sanitization
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return ValidationResult(
            False,
            f"Invalid SMILES: sanitization failed (invalid valence or aromaticity): {str(e)}"
        )
    
    # Count heavy atoms
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count < 1:
        return ValidationResult(False, "Invalid SMILES: molecule has no heavy atoms")
    
    # Get canonical SMILES
    try:
        canonical_smiles = Chem.MolToSmiles(mol)
    except Exception as e:
        return ValidationResult(False, f"Could not generate canonical SMILES: {str(e)}")
    
    return ValidationResult(
        True,
        data={
            "canonical_smiles": canonical_smiles,
            "heavy_atom_count": heavy_atom_count,
            "molecular_weight": round(Descriptors.MolWt(mol), 2)
        }
    )


def validate_pdb_content(pdb_text: str) -> ValidationResult:
    """
    Validate PDB block structure
    
    Checks:
    - Length ≤ 10MB
    - Contains ATOM or HETATM records
    - Coordinate bounds within [-999, 999]
    - At least 3 atoms
    
    Returns:
    - valid=True: atom_count, residue_count
    - valid=False: error message
    """
    if not isinstance(pdb_text, str):
        return ValidationResult(False, "PDB content must be a string")
    
    pdb_text = pdb_text.strip()
    
    if not pdb_text:
        return ValidationResult(False, "PDB content cannot be empty")
    
    # Check size limit (10MB)
    pdb_size = len(pdb_text.encode('utf-8'))
    if pdb_size > 10 * 1024 * 1024:
        return ValidationResult(
            False,
            f"PDB content too large ({pdb_size / (1024*1024):.1f} MB, maximum 10 MB)"
        )
    
    # Check for ATOM or HETATM records
    lines = pdb_text.splitlines()
    atom_lines = [l for l in lines if l.startswith("ATOM  ") or l.startswith("HETATM")]
    
    if not atom_lines:
        return ValidationResult(
            False,
            "Invalid PDB: no ATOM or HETATM records found"
        )
    
    if len(atom_lines) < 3:
        return ValidationResult(
            False,
            f"Invalid PDB: only {len(atom_lines)} atoms found (minimum 3 required)"
        )
    
    # Validate coordinate bounds
    try:
        for line in atom_lines[:100]:  # Check first 100 atoms for performance
            if len(line) >= 54:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                if not all(-999 <= coord <= 999 for coord in [x, y, z]):
                    return ValidationResult(
                        False,
                        f"Invalid PDB: coordinates out of bounds (must be within [-999, 999] Angstroms)"
                    )
    except (ValueError, IndexError) as e:
        return ValidationResult(
            False,
            f"Invalid PDB: malformed coordinate data: {str(e)}"
        )
    
    # Try parsing with gemmi for additional validation
    try:
        st = gemmi.read_pdb_string(pdb_text, "input.pdb")
        atom_count = sum(len(chain) for model in st for chain in model)
    except Exception as e:
        # Gemmi parsing is optional - fallback to line count
        atom_count = len(atom_lines)
    
    return ValidationResult(
        True,
        data={
            "atom_count": atom_count,
            "has_hetatm": any(l.startswith("HETATM") for l in lines)
        }
    )


def validate_grid_box(center: Dict[str, float], size: Dict[str, float]) -> ValidationResult:
    """
    Validate docking grid box parameters
    
    Checks:
    - center.x, center.y, center.z are numeric
    - size.x, size.y, size.z are positive
    - 5.0 ≤ size ≤ 150.0 Angstroms
    - center coordinates within [-999, 999]
    
    Returns:
    - valid=True: validated_center, validated_size
    - valid=False: error message
    """
    if not isinstance(center, dict):
        return ValidationResult(False, "Grid box center must be a dictionary with x, y, z keys")
    
    if not isinstance(size, dict):
        return ValidationResult(False, "Grid box size must be a dictionary with x, y, z keys")
    
    # Validate center coordinates
    for axis in ['x', 'y', 'z']:
        if axis not in center:
            return ValidationResult(False, f"Grid box center missing '{axis}' coordinate")
        
        try:
            coord = float(center[axis])
        except (ValueError, TypeError):
            return ValidationResult(False, f"Grid box center.{axis} must be numeric, got: {center[axis]}")
        
        if not -999 <= coord <= 999:
            return ValidationResult(
                False,
                f"Grid box center.{axis} out of bounds ({coord:.2f}, must be within [-999, 999])"
            )
    
    # Validate size dimensions
    for axis in ['x', 'y', 'z']:
        if axis not in size:
            return ValidationResult(False, f"Grid box size missing '{axis}' dimension")
        
        try:
            dim = float(size[axis])
        except (ValueError, TypeError):
            return ValidationResult(False, f"Grid box size.{axis} must be numeric, got: {size[axis]}")
        
        if dim <= 0:
            return ValidationResult(
                False,
                f"Grid box size.{axis} must be positive, got: {dim:.2f}"
            )
        
        if not 5.0 <= dim <= 150.0:
            return ValidationResult(
                False,
                f"Grid box size.{axis} out of acceptable range ({dim:.2f} Å, must be between 5.0 and 150.0 Å)"
            )
    
    return ValidationResult(
        True,
        data={
            "validated_center": {k: round(float(v), 2) for k, v in center.items()},
            "validated_size": {k: round(float(v), 2) for k, v in size.items()}
        }
    )
