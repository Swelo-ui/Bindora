import pytest
import gemmi
from pathlib import Path
from backend.services.docking import DockingEngine

SAMPLE_PDBQT_LIGAND = """REMARK  Name = Custom Ligand
ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00    +0.000 C
ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00    +0.000 C
ENDROOT
TORSDOF 0
"""

def test_cif_receptor_preparation():
    # Convert real 1IEP PDB to genuine mmCIF
    cache_file = Path("H:/AnuDock/data/cache/1IEP.pdb")
    if cache_file.exists():
        st = gemmi.read_structure(str(cache_file))
        cif_content = st.make_mmcif_document().as_string()
        result = DockingEngine.prepare_receptor(cif_content)
        assert result["atom_count"] > 0
        assert "cleaned_pdb" in result
        assert "pdbqt_text" in result
        assert "detected_pocket" in result
        assert "ATOM" in result["cleaned_pdb"]

def test_pdbqt_ligand_preparation():
    result = DockingEngine.prepare_ligand(SAMPLE_PDBQT_LIGAND)
    assert "pdbqt_text" in result
    assert "pdb_block" in result
    assert "canonical_smiles" in result
    assert result["heavy_atom_count"] >= 2
