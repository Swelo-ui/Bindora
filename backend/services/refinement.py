import io
import math
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

class ComplexRefinementService:
    """
    Post-docking pose energy minimization for binding-site relaxation.

    NOTE: This service performs structural energy minimization (force-field relaxation)
    to refine docked pose geometry. It does NOT compute MM-GBSA binding free energy.
    True MM-GBSA requires three separate simulations: complex, receptor-alone, and
    ligand-alone — each with full implicit solvation — which is outside the scope of
    this lightweight refinement step.

    Output fields are named to reflect what is actually computed:
      - complex_relaxation_delta_kcal : potential energy change after OpenMM minimization
      - ligand_strain_relaxation_kcal : MMFF94 intramolecular strain released in ligand
    Neither of these is a binding free energy estimate.
    """

    @classmethod
    def refine_pose(
        cls,
        receptor_pdb: str,
        docked_pdb_or_pdbqt: str,
        smiles: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Energy minimize docked pose in the receptor pocket.
        Uses OpenMM implicit solvent (GBn2) if available, otherwise RDKit MMFF94
        ligand-only relaxation. Returns geometry/strain metrics — NOT binding free energy.
        """
        # Try OpenMM GBn2 implicit solvent first
        openmm_res = cls._refine_openmm(receptor_pdb, docked_pdb_or_pdbqt)
        if openmm_res and "error" not in openmm_res:
            return openmm_res

        # Fall back to RDKit MMFF94 force field minimization
        return cls._refine_rdkit(receptor_pdb, docked_pdb_or_pdbqt, smiles=smiles)

    @classmethod
    def _refine_openmm(cls, receptor_pdb: str, docked_pdb: str) -> Optional[Dict[str, Any]]:
        """
        Run short OpenMM GBn2 implicit solvent minimization on the complex.
        Reports potential energy change (relaxation delta), NOT binding free energy.
        """
        try:
            import openmm as mm
            from openmm import app
            from openmm import unit

            # Combine receptor and docked ligand into single PDB
            combined_lines = []
            for l in receptor_pdb.splitlines():
                if l.startswith(("ATOM  ", "HETATM")):
                    combined_lines.append(l)
            for l in docked_pdb.splitlines():
                if l.startswith(("ATOM  ", "HETATM")):
                    combined_lines.append(l)

            pdb_io = io.StringIO("\n".join(combined_lines) + "\nEND\n")
            pdb = app.PDBFile(pdb_io)

            forcefield = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds
            )

            integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
            simulation = app.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)

            # Initial potential energy
            state_initial = simulation.context.getState(getEnergy=True)
            e_init = state_initial.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

            # Minimize 150 steps
            simulation.minimizeEnergy(maxIterations=150)

            state_min = simulation.context.getState(getEnergy=True, getPositions=True)
            e_min = state_min.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            delta_e = e_min - e_init

            return {
                "method": "OpenMM 8.x GBn2 Implicit Solvent — Complex Geometry Relaxation",
                "method_note": (
                    "Reports complex potential energy change after short minimization. "
                    "This is NOT a binding free energy (ΔG_bind). "
                    "True MM-GBSA requires separate receptor-alone and ligand-alone simulations."
                ),
                "initial_energy_kcal": round(e_init, 2),
                "minimized_energy_kcal": round(e_min, 2),
                "complex_relaxation_delta_kcal": round(delta_e, 2),
                "status": "Minimized & Solvated",
                "relaxation_steps": 150
            }
        except Exception as e:
            return {"error": str(e)}

    @classmethod
    def _refine_rdkit(
        cls,
        receptor_pdb: str,
        docked_pdb_or_pdbqt: str,
        smiles: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        RDKit MMFF94 ligand-only intramolecular strain relaxation.
        Receptor is NOT included. This measures how much intramolecular strain
        the ligand carries in its docked conformation — NOT binding free energy.
        """
        try:
            # Parse ligand pose into RDKit Mol
            mol = None
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
            if not mol:
                mol = Chem.MolFromPDBBlock(docked_pdb_or_pdbqt, removeHs=False)

            if not mol:
                return {
                    "method": "RDKit MMFF94 Ligand Strain Relaxation",
                    "status": "Could not parse docked pose for minimization",
                    "ligand_strain_relaxation_kcal": None,
                    "method_note": (
                        "Ligand-only MMFF94 strain relaxation. Receptor NOT included. "
                        "This is NOT a binding free energy estimate."
                    )
                }

            mol_h = Chem.AddHs(mol)
            # Embed if coordinates are missing
            if mol_h.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol_h, randomSeed=42)

            props = AllChem.MMFFGetMoleculeProperties(mol_h)
            if not props:
                return {
                    "method": "RDKit UFF Fallback",
                    "status": "MMFF94 parameterization unavailable for structure",
                    "ligand_strain_relaxation_kcal": None,
                    "method_note": (
                        "Ligand-only strain relaxation. Receptor NOT included. "
                        "This is NOT a binding free energy estimate."
                    )
                }

            ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
            e_init = ff.CalcEnergy() if ff else 0.0

            if ff:
                ff.Minimize(maxIts=200)
                e_min = ff.CalcEnergy()
            else:
                e_min = e_init

            strain_delta = round(abs(e_init - e_min), 2)
            return {
                "method": "RDKit MMFF94 Ligand Strain Relaxation",
                "method_note": (
                    "Ligand-only MMFF94 intramolecular strain relaxation. "
                    "Receptor NOT included. This is NOT a binding free energy estimate. "
                    "Lower strain delta = docked conformation closer to low-energy geometry."
                ),
                "initial_strain_kcal": round(e_init, 2),
                "minimized_strain_kcal": round(e_min, 2),
                "ligand_strain_relaxation_kcal": strain_delta,
                "status": "Pose Relaxed & Minimized",
                "relaxation_steps": 200
            }
        except Exception as e:
            return {
                "method": "RDKit Minimization",
                "error": str(e),
                "status": "Minimization skipped"
            }
