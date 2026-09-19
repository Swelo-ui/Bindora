// 3D Molecular Viewer controller using 3Dmol.js
// High-performance WebGL implementation tuned for low-memory systems (2GB RAM safe)

class MolecularViewer {
  constructor(elementId) {
    this.elementId = elementId;
    this.viewer = null;
    this.receptorModel = null;
    this.ligandModel = null;
    this.crystModel = null;
    this.surfaceObj = null;
    this.interactionShapes = [];
    this.gridBoxShapes = [];
    this.measureMode = false;
    this.measureAtoms = [];
    this.measureShapes = [];
    this.measureMarkers = [];
    this.measureLabels = [];   // track only distance labels for selective removal
    this.isSpinning = false;
    this.flexHighlights = new Map();
    this.lastInteractingResidues = new Set();
    
    // Default display settings
    this.settings = {
      proteinStyle: 'cartoon',
      proteinColor: 'spectrum',
      ligandStyle: 'ballstick',
      ligandColor: 'cyanCarbon',
      showSurface: false,
      surfaceType: 'VDW',
      surfaceOpacity: 0.55,
      showGridBox: false,
      showHbonds: true,
      showSaltBridges: true,
      showPiStacking: true,
      showPiCation: true,
      showHalogenBonds: true,
      showHydrophobic: false,
      showResidueLabels: true
    };

    // Surface tracking — prevents accumulation bug
    this._surfaceId = null;
    this._surfaceUpdating = false;
    this.lastInteractions = null;

    this.init();
  }

  init() {
    const el = document.getElementById(this.elementId);
    if (!el) return;

    const M = window["$" + "3Dmol"];
    if (M) {
      this.viewer = M.createViewer(el, {
        backgroundColor: '0x090a0f',
        disableFog: true,
        antialias: true
      });

      // Completely disable depth fogging so the receptor protein never fades into black fog
      try {
        if (this.viewer.getConfig()) {
          this.viewer.getConfig().disableFog = true;
        }
        if (typeof this.viewer.enableFog === 'function') {
          this.viewer.enableFog(false);
        }
        // Intercept setSlabAndFog to guarantee ample camera depth and zero black fog at all times
        if (typeof this.viewer.setSlabAndFog === 'function') {
          const origSetSlabAndFog = this.viewer.setSlabAndFog.bind(this.viewer);
          this.viewer.setSlabAndFog = () => {
            this.viewer.slabNear = Math.min(this.viewer.slabNear || -50, -300);
            this.viewer.slabFar = Math.max(this.viewer.slabFar || 50, 300);
            if (this.viewer.getConfig()) this.viewer.getConfig().disableFog = true;
            origSetSlabAndFog();
            if (this.viewer.scene && this.viewer.scene.fog) {
              this.viewer.scene.fog.near = this.viewer.scene.fog.far;
            }
          };
        }
      } catch (fe) {
        console.warn('[3Dmol] Fog configuration notice:', fe);
      }

      console.log('[3Dmol] Viewer initialized successfully with zero-fog high-depth camera.');
    } else {
      console.warn('[3Dmol] 3Dmol.js library not loaded yet.');
    }
  }

  clear() {
    if (!this.viewer) return;
    // Remove surface first to avoid accumulation
    this._removeSurface();
    this.viewer.clear();
    this.receptorModel = null;
    this.ligandModel = null;
    this.surfaceObj = null;
    this._surfaceId = null;
    this._surfaceUpdating = false;
    this.lastInteractions = null;
    this.clearInteractions();
    this.clearMeasurements();
    this.clearGridBox();
    this.clearRedockOverlay();
  }

  clearInteractions() {
    if (!this.viewer) return;
    this.interactionShapes.forEach(s => {
      try { this.viewer.removeShape(s); } catch (e) {}
    });
    this.interactionShapes = [];
    this.viewer.removeAllLabels();
    this.clearFlexibleHighlights();
  }

  clearMeasurements() {
    if (!this.viewer) return;
    this.measureShapes.forEach(s => {
      try { this.viewer.removeShape(s); } catch (e) {}
    });
    this.measureShapes = [];
    this.measureMarkers.forEach(m => {
      try { this.viewer.removeShape(m); } catch (e) {}
    });
    this.measureMarkers = [];
    // Remove only distance labels we added — not residue labels
    this.measureLabels.forEach(lbl => {
      try { this.viewer.removeLabel(lbl); } catch (e) {}
    });
    this.measureLabels = [];
    this.measureAtoms = [];
    this.viewer.render();
  }

  clearGridBox() {
    if (!this.viewer) return;
    this.gridBoxShapes.forEach(s => {
      try { this.viewer.removeShape(s); } catch (e) {}
    });
    this.gridBoxShapes = [];
    this.viewer.render();
  }

  loadReceptor(pdbContent, format = 'pdb') {
    if (!this.viewer || !pdbContent) return;
    
    if (this.receptorModel) {
      this.viewer.removeModel(this.receptorModel);
      this.receptorModel = null;
    }

    try {
      this.receptorModel = this.viewer.addModel(pdbContent, format);
      this.applyProteinStyle();
      this.viewer.zoomTo();
      this.viewer.render();

      if (this.settings.showSurface) {
        this.updateSurface();
      }
    } catch (e) {
      console.error('[3Dmol] Failed to load receptor:', e);
    }
  }

  applyProteinStyle() {
    if (!this.receptorModel) return;

    const style = this.settings.proteinStyle || 'cartoon';
    const colorType = this.settings.proteinColor || 'spectrum';

    // Remove any previously-added pocket surface when switching styles
    if (this._pocketSurface) {
      try { this.viewer.removeSurface(this._pocketSurface); } catch (e) {}
      this._pocketSurface = null;
    }

    let colorScheme = {};
    if (colorType === 'spectrum') {
      colorScheme = { color: 'spectrum' };
    } else if (colorType === 'chain') {
      colorScheme = { colorscheme: 'chain' };
    } else if (colorType === 'secondary') {
      colorScheme = { colorscheme: 'ssJmol' };
    } else {
      colorScheme = { color: 'lightgray' };
    }

    if (style === 'cartoon') {
      this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.98, thickness: 0.48 } });

    } else if (style === 'stick') {
      this.receptorModel.setStyle({}, { stick: { radius: 0.18, colorscheme: 'Jmol' } });

    } else if (style === 'sphere') {
      this.receptorModel.setStyle({}, { sphere: { radius: 0.75, colorscheme: 'Jmol' } });

    } else if (style === 'ribbon') {
      this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.98, style: 'trace', thickness: 0.6 } });

    } else if (style === 'tube') {
      // FIX: 'tube' is not a valid 3Dmol cartoon style — use 'trace' with thick width instead
      this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, style: 'trace', thickness: 0.85, opacity: 0.96 } });

    } else if (style === 'pocket_surface') {
      if (this.ligandModel) {
        // Faded cartoon ribbon base + translucent VDW pocket surface around ligand (5Å)
        this.receptorModel.setStyle({}, { cartoon: { color: 'lightgray', opacity: 0.38, thickness: 0.42 } });
        try {
          this._pocketSurface = this.viewer.addSurface(
            $3Dmol.SurfaceType.VDW,
            { opacity: 0.78, colorscheme: 'whiteCarbon' },
            { model: this.receptorModel },
            { model: this.ligandModel, within: { distance: 5.5, sel: {} } }
          );
        } catch (e) {
          // FIX: fallback to cartoon if surface fails
          console.warn('[3Dmol] Pocket surface failed, falling back to cartoon:', e);
          this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.98, thickness: 0.48 } });
        }
      } else {
        // FIX: No ligand loaded — show normal cartoon instead of near-invisible ghost
        this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.98, thickness: 0.48 } });
      }

    } else if (style === 'pocket_sticks') {
      if (this.ligandModel) {
        const ligandAtoms = this.ligandModel.selectedAtoms({});
        if (ligandAtoms && ligandAtoms.length > 0) {
          // FIX: Use visible ghost ribbon (0.28 opacity) — was 0.12 which is invisible on dark bg
          this.receptorModel.setStyle({}, { cartoon: { color: '#64748b', opacity: 0.28, thickness: 0.36 } });
          const receptorAtoms = this.receptorModel.selectedAtoms({});
          const nearResidues = new Set();
          const thresh2 = 5.5 * 5.5;
          receptorAtoms.forEach(ra => {
            for (const la of ligandAtoms) {
              const dx = ra.x - la.x, dy = ra.y - la.y, dz = ra.z - la.z;
              if (dx * dx + dy * dy + dz * dz <= thresh2) {
                nearResidues.add(`${ra.chain}:${ra.resi}`);
                break;
              }
            }
          });
          nearResidues.forEach(key => {
            const [chain, resi] = key.split(':');
            // Pocket residues: amino-colored sticks only (no cartoon override here)
            this.receptorModel.setStyle(
              { chain, resi: parseInt(resi) },
              { stick: { radius: 0.24, colorscheme: 'amino' } }
            );
          });
        } else {
          // Ligand model exists but no atoms — fallback
          this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.98, thickness: 0.48 } });
        }
      } else {
        // FIX: No ligand loaded — show normal cartoon as fallback
        this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.98, thickness: 0.48 } });
      }

    } else if (style === 'bfactor') {
      // FIX: Use wide range (0–100) so color variation always visible regardless of PDB source
      // rwb gradient: low B-factor = red (warm, disordered), high B-factor = blue (cool, rigid)
      // This matches ChimeraX default B-factor coloring convention
      this.receptorModel.setStyle({}, {
        cartoon: {
          colorscheme: { gradient: 'rwb', prop: 'b', min: 0, max: 100 },
          opacity: 0.97,
          thickness: 0.48
        }
      });

    } else if (style === 'electrostatic') {
      // Approximate electrostatic ribbon: Asp/Glu=red, Arg/Lys/His=blue, polar=green, neutral=white
      this.receptorModel.setStyle({}, { cartoon: { color: '#e2e8f0', opacity: 0.92, thickness: 0.46 } });
      // Acidic — red (negative charge)
      this.receptorModel.setStyle({ resn: 'ASP' }, { cartoon: { color: '#dc2626', opacity: 0.97, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'GLU' }, { cartoon: { color: '#ef4444', opacity: 0.97, thickness: 0.46 } });
      // Basic — blue (positive charge)
      this.receptorModel.setStyle({ resn: 'ARG' }, { cartoon: { color: '#1d4ed8', opacity: 0.97, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'LYS' }, { cartoon: { color: '#2563eb', opacity: 0.97, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'HIS' }, { cartoon: { color: '#60a5fa', opacity: 0.95, thickness: 0.46 } });
      // Polar uncharged — green
      this.receptorModel.setStyle({ resn: 'SER' }, { cartoon: { color: '#22c55e', opacity: 0.88, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'THR' }, { cartoon: { color: '#22c55e', opacity: 0.88, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'ASN' }, { cartoon: { color: '#4ade80', opacity: 0.88, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'GLN' }, { cartoon: { color: '#4ade80', opacity: 0.88, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'TYR' }, { cartoon: { color: '#86efac', opacity: 0.85, thickness: 0.46 } });
      this.receptorModel.setStyle({ resn: 'TRP' }, { cartoon: { color: '#86efac', opacity: 0.85, thickness: 0.46 } });

    } else if (style === 'hidden') {
      this.receptorModel.setStyle({}, {});
    }
    
    this.viewer.render();
  }



  loadLigand(pdbOrPdbqtContent, format = 'pdb') {
    if (!this.viewer || !pdbOrPdbqtContent) return;

    this.clearRedockOverlay();

    if (this.ligandModel) {
      this.viewer.removeModel(this.ligandModel);
      this.ligandModel = null;
    }

    try {
      this.ligandModel = this.viewer.addModel(pdbOrPdbqtContent, format);
      this.applyLigandStyle();

      this.viewer.zoomTo({ model: this.ligandModel }, 600);
      if (typeof this.viewer.setSlab === 'function') {
        this.viewer.setSlab(-300, 300);
      }
      this.viewer.render();
      setTimeout(() => {
        if (this.viewer && typeof this.viewer.setSlab === 'function') {
          this.viewer.setSlab(-300, 300);
          this.viewer.render();
        }
      }, 650);

      if (this.settings.showSurface) {
        this.updateSurface();
      }
    } catch (e) {
      console.error('[3Dmol] Failed to load ligand:', e);
    }
  }

  loadRedockOverlay(crystPdb, dockedPdb) {
    if (!this.viewer) return;

    if (this.ligandModel) {
      this.viewer.removeModel(this.ligandModel);
      this.ligandModel = null;
    }
    if (this.crystModel) {
      this.viewer.removeModel(this.crystModel);
      this.crystModel = null;
    }

    try {
      // 1. Crystal Reference Ligand (distinctive gold/amber stick representation)
      if (crystPdb) {
        this.crystModel = this.viewer.addModel(crystPdb, 'pdb');
        this.crystModel.setStyle({}, {
          stick: { radius: 0.22, color: 'goldenrod', opacity: 0.85 }
        });
      }

      // 2. Vina Redocked Pose (cyan ball-and-stick)
      if (dockedPdb) {
        this.ligandModel = this.viewer.addModel(dockedPdb, 'pdb');
        this.applyLigandStyle();
      }

      const focusModel = this.ligandModel || this.crystModel;
      if (focusModel) {
        this.viewer.zoomTo({ model: focusModel }, 600);
      }
      this.viewer.render();
    } catch (e) {
      console.error('[3Dmol] Failed to load redock overlay:', e);
    }
  }

  clearRedockOverlay() {
    if (this.crystModel && this.viewer) {
      try {
        this.viewer.removeModel(this.crystModel);
      } catch (e) {}
      this.crystModel = null;
      this.viewer.render();
    }
  }

  applyLigandStyle() {
    if (!this.ligandModel) return;

    const style = this.settings.ligandStyle || 'ballstick';
    const colorScheme = this.settings.ligandColor || 'cyanCarbon';

    if (style === 'ballstick') {
      this.ligandModel.setStyle({}, {
        stick: { radius: 0.25, colorscheme: colorScheme },
        sphere: { radius: 0.35, colorscheme: colorScheme }
      });
    } else if (style === 'stick') {
      this.ligandModel.setStyle({}, {
        stick: { radius: 0.32, colorscheme: colorScheme }
      });
    } else if (style === 'sphere') {
      this.ligandModel.setStyle({}, {
        sphere: { radius: 0.75, colorscheme: colorScheme }
      });
    } else if (style === 'line') {
      this.ligandModel.setStyle({}, {
        line: { linewidth: 2, colorscheme: colorScheme }
      });
    } else if (style === 'hidden') {
      this.ligandModel.setStyle({}, {});
    }

    this.viewer.render();
  }

  renderInteractions(interactions) {
    if (!this.viewer) return;
    if (interactions) {
      this.lastInteractions = interactions;
    } else {
      interactions = this.lastInteractions;
    }
    if (!interactions) return;
    this.clearInteractions();

    const {
      hydrogen_bonds = [],
      salt_bridges = [],
      pi_stacking = [],
      pi_cation = [],
      halogen_bonds = [],
      hydrophobic_contacts = [],
      interacting_residues = []
    } = interactions;

    this.lastInteractingResidues = new Set();
    if (!this.residueLabels) this.residueLabels = new Map();
    this.residueLabels.clear();

    if (this.receptorModel && interacting_residues.length > 0) {
      interacting_residues.forEach(resStr => {
        const parts = resStr.split(' ');
        const resName = parts[0];
        const numChain = parts[1] || '';
        const [resNumStr, chain] = numChain.split(':');
        const resNum = parseInt(resNumStr);
        const key = `${chain || 'A'}:${resNum}`;
        this.lastInteractingResidues.add(key);

        const sel = { resi: resNum };
        if (chain) sel.chain = chain;

        this.receptorModel.setStyle(sel, {
          cartoon: { color: '#38bdf8', opacity: 0.9 },
          stick: { radius: 0.22, colorscheme: 'amino' }
        });

        if (this.settings.showResidueLabels) {
          const atoms = this.receptorModel.selectedAtoms(sel);
          if (atoms && atoms.length > 0) {
            const ca = atoms.find(a => a.atom === 'CA') || atoms[0];
            const sidechainAtoms = atoms.filter(a => a.atom !== 'CA' && a.atom !== 'C' && a.atom !== 'N' && a.atom !== 'O');
            let pos = { x: ca.x, y: ca.y, z: ca.z };
            if (sidechainAtoms.length > 0) {
              const avgX = sidechainAtoms.reduce((s, a) => s + a.x, 0) / sidechainAtoms.length;
              const avgY = sidechainAtoms.reduce((s, a) => s + a.y, 0) / sidechainAtoms.length;
              const avgZ = sidechainAtoms.reduce((s, a) => s + a.z, 0) / sidechainAtoms.length;
              const dx = avgX - ca.x;
              const dy = avgY - ca.y;
              const dz = avgZ - ca.z;
              const dist = Math.hypot(dx, dy, dz) || 1.0;
              // Shift 1.6 Angstroms toward sidechain centroid to clear cartoon ribbon and prevent label collision
              pos = {
                x: ca.x + (dx / dist) * 1.6,
                y: ca.y + (dy / dist) * 1.6,
                z: ca.z + (dz / dist) * 1.6
              };
            }
            const displayRes = (resName ? resName + ' ' : '') + resNum;
            const lbl = this.viewer.addLabel(displayRes, {
              position: pos,
              backgroundColor: 'rgba(11, 12, 17, 0.92)',
              fontColor: '#38bdf8',
              fontSize: 11,
              borderThickness: 1,
              borderColor: '#0284c7'
            });
            if (lbl) this.residueLabels.set(key, lbl);
          }
        }
      });
    }

    // 3D Distance Label Collision Avoidance Placer
    const placed3DLabels = [];
    const getOptimalLabelPos = (start, end, defaultT = 0.5) => {
      const vx = end.x - start.x;
      const vy = end.y - start.y;
      const vz = end.z - start.z;
      const tCandidates = [defaultT, 0.32, 0.68, 0.22, 0.78];

      for (const t of tCandidates) {
        const candidate = {
          x: start.x + vx * t,
          y: start.y + vy * t,
          z: start.z + vz * t
        };

        const collision = placed3DLabels.some(pl => {
          const d = Math.hypot(candidate.x - pl.x, candidate.y - pl.y, candidate.z - pl.z);
          return d < 1.8;
        });

        if (!collision) {
          placed3DLabels.push(candidate);
          return candidate;
        }
      }

      let px = -vy, py = vx, pz = 0;
      let pLen = Math.hypot(px, py, pz);
      if (pLen < 0.01) {
        px = 0; py = -vz; pz = vy;
        pLen = Math.hypot(px, py, pz) || 1.0;
      }
      const shift = 1.4;
      const shifted = {
        x: start.x + vx * defaultT + (px / pLen) * shift,
        y: start.y + vy * defaultT + (py / pLen) * shift,
        z: start.z + vz * defaultT + (pz / pLen) * shift
      };
      placed3DLabels.push(shifted);
      return shifted;
    };

    // 1. Directional Hydrogen Bonds (Yellow)
    if (this.settings.showHbonds && hydrogen_bonds.length > 0) {
      hydrogen_bonds.forEach(hb => {
        if (!hb.start_coord || !hb.end_coord) return;
        const [lx, ly, lz] = hb.start_coord;
        const [rx, ry, rz] = hb.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.08,
          dashed: true,
          color: '#facc15',
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);

        const lblPos = getOptimalLabelPos({ x: lx, y: ly, z: lz }, { x: rx, y: ry, z: rz }, 0.44);
        this.viewer.addLabel(hb.distance + ' A', {
          position: lblPos,
          backgroundColor: 'rgba(11, 12, 17, 0.92)',
          fontColor: '#facc15',
          fontSize: 10,
          borderThickness: 1,
          borderColor: '#ca8a04'
        });
      });
    }

    // 2. Salt Bridges (Magenta)
    if (this.settings.showSaltBridges && salt_bridges.length > 0) {
      salt_bridges.forEach(sb => {
        if (!sb.start_coord || !sb.end_coord) return;
        const [lx, ly, lz] = sb.start_coord;
        const [rx, ry, rz] = sb.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.09,
          dashed: true,
          color: '#ec4899',
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);

        const lblPos = getOptimalLabelPos({ x: lx, y: ly, z: lz }, { x: rx, y: ry, z: rz }, 0.58);
        this.viewer.addLabel(`Salt ${sb.distance} A`, {
          position: lblPos,
          backgroundColor: 'rgba(11, 12, 17, 0.92)',
          fontColor: '#ec4899',
          fontSize: 10,
          borderThickness: 1,
          borderColor: '#db2777'
        });
      });
    }

    // 3. π-π Stacking (Emerald Green)
    if (this.settings.showPiStacking && pi_stacking.length > 0) {
      pi_stacking.forEach(ps => {
        if (!ps.start_coord || !ps.end_coord) return;
        const [lx, ly, lz] = ps.start_coord;
        const [rx, ry, rz] = ps.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.09,
          dashed: true,
          color: '#10b981',
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);

        const lblPos = getOptimalLabelPos({ x: lx, y: ly, z: lz }, { x: rx, y: ry, z: rz }, 0.50);
        this.viewer.addLabel(`π-π ${ps.distance} A`, {
          position: lblPos,
          backgroundColor: 'rgba(11, 12, 17, 0.92)',
          fontColor: '#10b981',
          fontSize: 10,
          borderThickness: 1,
          borderColor: '#059669'
        });
      });
    }

    // 4. π-Cation Interactions (Orange)
    if (this.settings.showPiCation && pi_cation.length > 0) {
      pi_cation.forEach(pc => {
        if (!pc.start_coord || !pc.end_coord) return;
        const [lx, ly, lz] = pc.start_coord;
        const [rx, ry, rz] = pc.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.08,
          dashed: true,
          color: '#f97316',
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);

        const lblPos = getOptimalLabelPos({ x: lx, y: ly, z: lz }, { x: rx, y: ry, z: rz }, 0.50);
        this.viewer.addLabel(`π-Cat ${pc.distance} A`, {
          position: lblPos,
          backgroundColor: 'rgba(11, 12, 17, 0.92)',
          fontColor: '#f97316',
          fontSize: 10,
          borderThickness: 1,
          borderColor: '#ea580c'
        });
      });
    }

    // 5. Halogen Bonds (Purple)
    if (this.settings.showHalogenBonds && halogen_bonds.length > 0) {
      halogen_bonds.forEach(hal => {
        if (!hal.start_coord || !hal.end_coord) return;
        const [lx, ly, lz] = hal.start_coord;
        const [rx, ry, rz] = hal.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.08,
          dashed: true,
          color: '#a855f7',
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);

        const lblPos = getOptimalLabelPos({ x: lx, y: ly, z: lz }, { x: rx, y: ry, z: rz }, 0.50);
        this.viewer.addLabel(`Hal ${hal.distance} A`, {
          position: lblPos,
          backgroundColor: 'rgba(11, 12, 17, 0.92)',
          fontColor: '#a855f7',
          fontSize: 10,
          borderThickness: 1,
          borderColor: '#9333ea'
        });
      });
    }

    // 6. Hydrophobic Contacts (Sky Blue)
    if (this.settings.showHydrophobic && hydrophobic_contacts.length > 0) {
      hydrophobic_contacts.forEach(hp => {
        if (!hp.start_coord || !hp.end_coord) return;
        const [lx, ly, lz] = hp.start_coord;
        const [rx, ry, rz] = hp.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.06,
          dashed: true,
          color: '#94a3b8',
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);
      });
    }

    this.viewer.render();
  }

  toggleHBonds(show) {
    this.settings.showHbonds = show;
    if (this.lastInteractions) this.renderInteractions(this.lastInteractions);
  }

  toggleSaltBridges(show) {
    this.settings.showSaltBridges = show;
    if (this.lastInteractions) this.renderInteractions(this.lastInteractions);
  }

  togglePiStacking(show) {
    this.settings.showPiStacking = show;
    if (this.lastInteractions) this.renderInteractions(this.lastInteractions);
  }

  togglePiCation(show) {
    this.settings.showPiCation = show;
    if (this.lastInteractions) this.renderInteractions(this.lastInteractions);
  }

  toggleHalogenBonds(show) {
    this.settings.showHalogenBonds = show;
    if (this.lastInteractions) this.renderInteractions(this.lastInteractions);
  }

  toggleHydrophobic(show) {
    this.settings.showHydrophobic = show;
    if (this.lastInteractions) this.renderInteractions(this.lastInteractions);
  }

  // Internal helper: safely remove ALL surfaces (prevents accumulation)
  _removeSurface() {
    if (!this.viewer) return;
    // Try removeAllSurfaces first (safest, removes every surface layer)
    try {
      if (typeof this.viewer.removeAllSurfaces === 'function') {
        this.viewer.removeAllSurfaces();
        this._surfaceId = null;
        this.surfaceObj = null;
        return;
      }
    } catch (e) {}
    // Fallback: remove by stored ID
    if (this._surfaceId !== null && this._surfaceId !== undefined) {
      try { this.viewer.removeSurface(this._surfaceId); } catch (e) {}
      this._surfaceId = null;
    }
    if (this.surfaceObj !== null && this.surfaceObj !== undefined) {
      try { this.viewer.removeSurface(this.surfaceObj); } catch (e) {}
      this.surfaceObj = null;
    }
  }

  toggleSurface(show) {
    this.settings.showSurface = show;
    this._surfaceUpdating = false; // Reset lock so new state takes effect
    this.updateSurface();
  }

  updateSurface() {
    if (!this.viewer) return;
    // Prevent concurrent surface additions (race condition guard)
    if (this._surfaceUpdating) return;

    // Always remove all existing surfaces before doing anything
    this._removeSurface();

    if (!this.settings.showSurface) {
      this.viewer.render();
      return;
    }

    const M = window["$" + "3Dmol"];
    const SType = M ? M.SurfaceType.VDW : 1;
    this._surfaceUpdating = true;

    try {
      let surfResult;
      if (this.receptorModel && this.ligandModel) {
        surfResult = this.viewer.addSurface(SType, {
          opacity: this.settings.surfaceOpacity || 0.5,
          color: '#38bdf8'
        }, { model: this.receptorModel, within: { distance: 6.0, sel: { model: this.ligandModel } } });
      } else if (this.ligandModel) {
        surfResult = this.viewer.addSurface(SType, {
          opacity: 0.55,
          color: '#06b6d4'
        }, { model: this.ligandModel });
      } else if (this.receptorModel) {
        surfResult = this.viewer.addSurface(SType, {
          opacity: 0.4,
          color: '#64748b'
        }, { model: this.receptorModel });
      }

      // Handle both Promise (newer 3Dmol) and numeric ID (older 3Dmol)
      if (surfResult && typeof surfResult.then === 'function') {
        surfResult.then(id => {
          this._surfaceId = id;
          this.surfaceObj = id;
          this._surfaceUpdating = false;
        }).catch(err => {
          console.error('[3Dmol] Surface async error:', err);
          this._surfaceUpdating = false;
        });
      } else if (surfResult !== undefined && surfResult !== null) {
        this._surfaceId = surfResult;
        this.surfaceObj = surfResult;
        this._surfaceUpdating = false;
      } else {
        this._surfaceUpdating = false;
      }
    } catch (err) {
      console.error('[3Dmol] Surface generation error:', err);
      this._surfaceUpdating = false;
    }

    this.viewer.render();
  }

  renderGridBox(center, size, visible = true) {
    this.clearGridBox();
    if (!this.viewer || !visible || !center || !size) return;

    const cx = parseFloat(center.x) || 0;
    const cy = parseFloat(center.y) || 0;
    const cz = parseFloat(center.z) || 0;

    const sx = Math.max(1, parseFloat(size.x) || 20);
    const sy = Math.max(1, parseFloat(size.y) || 20);
    const sz = Math.max(1, parseFloat(size.z) || 20);

    const hx = sx / 2;
    const hy = sy / 2;
    const hz = sz / 2;

    const p = [
      { x: cx - hx, y: cy - hy, z: cz - hz },
      { x: cx + hx, y: cy - hy, z: cz - hz },
      { x: cx + hx, y: cy + hy, z: cz - hz },
      { x: cx - hx, y: cy + hy, z: cz - hz },
      { x: cx - hx, y: cy - hy, z: cz + hz },
      { x: cx + hx, y: cy - hy, z: cz + hz },
      { x: cx + hx, y: cy + hy, z: cz + hz },
      { x: cx - hx, y: cy + hy, z: cz + hz }
    ];

    const edges = [
      [0,1],[1,2],[2,3],[3,0],
      [4,5],[5,6],[6,7],[7,4],
      [0,4],[1,5],[2,6],[3,7]
    ];

    edges.forEach(([a, b]) => {
      const cyl = this.viewer.addCylinder({
        start: p[a],
        end: p[b],
        radius: 0.12,
        color: '#06b6d4'
      });
      this.gridBoxShapes.push(cyl);
    });

    const centerSphere = this.viewer.addSphere({
      center: { x: cx, y: cy, z: cz },
      radius: 0.45,
      color: '#f59e0b'
    });
    this.gridBoxShapes.push(centerSphere);

    this.viewer.render();
  }

  enableMeasurementMode(enable) {
    this.measureMode = enable;
    if (!this.viewer) return;

    const hud = document.getElementById('measure-hud');
    const hudText = document.getElementById('measure-hud-text');

    if (!enable) {
      if (hud) hud.classList.add('hidden');

      // Disable click handlers on all models
      if (this.receptorModel) {
        try { this.receptorModel.setClickable({}, false, null); } catch (e) {}
      }
      if (this.ligandModel) {
        try { this.ligandModel.setClickable({}, false, null); } catch (e) {}
      }
      if (this.crystModel) {
        try { this.crystModel.setClickable({}, false, null); } catch (e) {}
      }

      this.measureMarkers.forEach(m => {
        try { this.viewer.removeShape(m); } catch (e) {}
      });
      this.measureMarkers = [];
      this.measureAtoms = [];
      this.viewer.render();
      return;
    }

    if (hud) {
      hud.classList.remove('hidden');
      if (hudText) hudText.textContent = 'Measure Mode: Click any atom to begin measurement...';
    }

    // The atom click handler — shared across all models
    const atomClickHandler = (atom) => {
      if (!this.measureMode || !atom) return;
      // Ensure atom has valid coordinates
      const ax = atom.x !== undefined ? atom.x : (atom.xyz ? atom.xyz[0] : null);
      const ay = atom.y !== undefined ? atom.y : (atom.xyz ? atom.xyz[1] : null);
      const az = atom.z !== undefined ? atom.z : (atom.xyz ? atom.xyz[2] : null);
      if (ax === null || ay === null || az === null) return;

      const atomWithCoords = { ...atom, x: ax, y: ay, z: az };
      this.measureAtoms.push(atomWithCoords);

      if (this.measureAtoms.length === 1) {
        // First atom: show a red marker sphere
        const marker = this.viewer.addSphere({
          center: { x: ax, y: ay, z: az },
          radius: 0.45,
          color: '#f43f5e'
        });
        this.measureMarkers.push(marker);
        this.viewer.render();

        const atomName = atom.atom || atom.elem || 'Atom';
        const resName = (atom.resn || '') + (atom.resi !== undefined ? atom.resi : '');
        if (hudText) hudText.textContent = `Point 1: [${atomName} ${resName}] selected. Now click 2nd atom...`;

      } else if (this.measureAtoms.length >= 2) {
        const a1 = this.measureAtoms[0];
        const a2 = this.measureAtoms[1];
        const dx = a1.x - a2.x;
        const dy = a1.y - a2.y;
        const dz = a1.z - a2.z;
        const dist = Math.sqrt(dx*dx + dy*dy + dz*dz).toFixed(2);

        // Draw measurement cylinder
        const cyl = this.viewer.addCylinder({
          start:  { x: a1.x, y: a1.y, z: a1.z },
          end:    { x: a2.x, y: a2.y, z: a2.z },
          radius: 0.12,
          color: '#f43f5e',
          dashed: false,
          fromCap: 2,
          toCap: 2
        });
        this.measureShapes.push(cyl);

        // Label at midpoint — show distance in Å
        const mid = {
          x: (a1.x + a2.x) / 2,
          y: (a1.y + a2.y) / 2 + 0.5,
          z: (a1.z + a2.z) / 2
        };
        const lbl = this.viewer.addLabel(dist + ' Å', {
          position: mid,
          backgroundColor: 'rgba(244, 63, 94, 0.95)',
          fontColor: '#ffffff',
          fontSize: 13,
          borderThickness: 1,
          borderColor: '#ffffff',
          showBackground: true,
          inFront: true
        });
        if (lbl) this.measureLabels.push(lbl);

        // Clear first-atom marker
        this.measureMarkers.forEach(m => {
          try { this.viewer.removeShape(m); } catch (e) {}
        });
        this.measureMarkers = [];
        this.measureAtoms = [];

        const desc1 = (a1.atom || a1.elem || 'A') + (a1.resn ? ' ' + a1.resn + (a1.resi || '') : '');
        const desc2 = (a2.atom || a2.elem || 'A') + (a2.resn ? ' ' + a2.resn + (a2.resi || '') : '');
        if (hudText) hudText.textContent = `✓ ${dist} Å [${desc1} ↔ ${desc2}] — Click another pair to measure.`;
        this.viewer.render();
      }
    };

    // Apply clickable to ALL models with a generous click-sphere radius (0.8 Å)
    if (this.receptorModel) {
      this.receptorModel.setClickable({}, true, atomClickHandler);
    }
    if (this.ligandModel) {
      this.ligandModel.setClickable({}, true, atomClickHandler);
    }
    if (this.crystModel) {
      this.crystModel.setClickable({}, true, atomClickHandler);
    }

    this.viewer.render();
  }

  toggleSpin() {
    if (!this.viewer) return;
    this.isSpinning = !this.isSpinning;
    this.viewer.spin(this.isSpinning ? 'y' : false, 1);
    const spinBtn = document.getElementById('btn-toggle-spin');
    if (spinBtn) {
      if (this.isSpinning) {
        spinBtn.classList.add('text-cyan-400', 'bg-cyan-950/60');
        spinBtn.classList.remove('text-slate-400');
      } else {
        spinBtn.classList.remove('text-cyan-400', 'bg-cyan-950/60');
        spinBtn.classList.add('text-slate-400');
      }
    }
  }

  downloadScreenshot() {
    if (!this.viewer) return;
    try {
      const uri = this.viewer.pngURI();
      const a = document.createElement('a');
      a.href = uri;
      a.download = 'bindora_3d_structure_' + Date.now() + '.png';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    } catch (e) {
      console.error('[3Dmol] Screenshot failed:', e);
    }
  }

  focusResidue(chain, resNum, enable = true, resName = '') {
    if (!this.viewer || !this.receptorModel) return;
    try {
      const num = parseInt(resNum);
      const sel = { resi: num };
      if (chain) sel.chain = chain;
      const key = `${chain || 'A'}:${num}`;

      if (!this.flexHighlights) {
        this.flexHighlights = new Map();
      }

      // If disabling / unselecting: Cleanly remove wireframe surface, shapes, labels and restore original residue style
      if (!enable) {
        if (this.flexHighlights.has(key)) {
          const item = this.flexHighlights.get(key);
          if (item.surfaceId !== undefined && item.surfaceId !== null) {
            try { this.viewer.removeSurface(item.surfaceId); } catch (e) {}
          }
          if (item.shapes) {
            item.shapes.forEach(s => {
              try { this.viewer.removeShape(s); } catch (e) {}
            });
          }
          if (item.labels) {
            item.labels.forEach(l => {
              try { this.viewer.removeLabel(l); } catch (e) {}
            });
          }
          this.flexHighlights.delete(key);
        }

        // Restore original style cleanly (no lingering color)
        const isInteracting = this.lastInteractingResidues && this.lastInteractingResidues.has(key);
        if (isInteracting) {
          this.receptorModel.setStyle(sel, {
            cartoon: { color: '#38bdf8', opacity: 0.9 },
            stick: { radius: 0.22, colorscheme: 'amino' }
          });
        } else {
          const currentStyle = this.settings.proteinStyle || 'cartoon';
          if (currentStyle === 'cartoon') {
            this.receptorModel.setStyle(sel, {
              cartoon: { colorscheme: 'chain', opacity: 0.98, thickness: 0.48 }
            });
          } else {
            this.applyReceptorStyle();
          }
        }

        // Restore base residue label if active (ensuring no missing labels on deselect)
        if (this.settings.showResidueLabels && isInteracting) {
          const atoms = this.receptorModel.selectedAtoms(sel);
          if (atoms && atoms.length > 0) {
            const ca = atoms.find(a => a.atom === 'CA') || atoms[0];
            const displayRes = (resName ? resName + ' ' : '') + num;
            const baseLbl = this.viewer.addLabel(displayRes, {
              position: { x: ca.x, y: ca.y, z: ca.z },
              backgroundColor: 'rgba(15, 23, 42, 0.8)',
              fontColor: '#38bdf8',
              fontSize: 11,
              borderThickness: 1,
              borderColor: '#0284c7'
            });
            if (baseLbl) {
              if (!this.residueLabels) this.residueLabels = new Map();
              this.residueLabels.set(key, baseLbl);
            }
          }
        }

        this.viewer.render();
        return;
      }

      // If enabling / selecting: Professional Crystallographic Focus with Wireframe Structure
      this.viewer.zoomTo(sel, 600);

      // Clean existing shapes/labels/surface if re-selecting
      if (this.flexHighlights.has(key)) {
        const old = this.flexHighlights.get(key);
        if (old.surfaceId !== undefined && old.surfaceId !== null) {
          try { this.viewer.removeSurface(old.surfaceId); } catch (e) {}
        }
        if (old.shapes) old.shapes.forEach(s => { try { this.viewer.removeShape(s); } catch (e) {} });
        if (old.labels) old.labels.forEach(l => { try { this.viewer.removeLabel(l); } catch (e) {} });
      }

      // Suppress base interaction label to avoid double-label collision / overlapping
      if (this.residueLabels && this.residueLabels.has(key)) {
        try { this.viewer.removeLabel(this.residueLabels.get(key)); } catch (e) {}
        this.residueLabels.delete(key);
      }

      // Render sidechain with high-contrast elemental sticks (cyan carbons, red oxygen, blue nitrogen, yellow sulfur)
      this.receptorModel.setStyle(sel, {
        cartoon: { color: '#38bdf8', opacity: 0.95, thickness: 0.45 },
        stick: { radius: 0.28, colorscheme: 'cyanCarbon' }
      });

      // Add research-grade glowing wireframe VDW surface around the side chain
      const M = window["$" + "3Dmol"];
      const SType = M ? M.SurfaceType.VDW : 1;
      let surfPromise = null;
      try {
        const surfSel = { model: this.receptorModel, resi: num };
        if (chain) surfSel.chain = chain;
        surfPromise = this.viewer.addSurface(SType, {
          wireframe: true,
          opacity: 0.8,
          color: '#00f0ff'
        }, surfSel);
      } catch (e) {
        console.warn('[3Dmol] Residue wireframe surface warning:', e);
      }

      const shapes = [];
      const labels = [];
      const atoms = this.receptorModel.selectedAtoms(sel);

      const highlightItem = { shapes, labels, sel, surfaceId: null, resName };

      if (surfPromise && typeof surfPromise.then === 'function') {
        surfPromise.then(id => {
          if (this.flexHighlights && this.flexHighlights.has(key)) {
            this.flexHighlights.get(key).surfaceId = id;
          } else {
            try { this.viewer.removeSurface(id); } catch (e) {}
          }
        }).catch(err => {
          console.warn('[3Dmol] Residue wireframe async error:', err);
        });
      } else if (surfPromise !== undefined && surfPromise !== null) {
        highlightItem.surfaceId = surfPromise;
      }

      if (atoms && atoms.length > 0) {
        // Clean single academic 3D label without collision
        const ca = atoms.find(a => a.atom === 'CA') || atoms[0];
        const displayLabel = resName ? `${resName} ${num}:${chain || 'A'}` : `${num}:${chain || 'A'}`;
        const lbl = this.viewer.addLabel(`[Flex] ${displayLabel}`, {
          position: { x: ca.x, y: ca.y + 1.2, z: ca.z },
          backgroundColor: 'rgba(15, 23, 42, 0.95)',
          fontColor: '#38bdf8',
          fontSize: 11,
          borderThickness: 1.5,
          borderColor: '#0284c7',
          inFront: true
        });
        if (lbl) labels.push(lbl);
      }

      this.flexHighlights.set(key, highlightItem);
      this.viewer.render();
    } catch (e) {
      console.warn('[3Dmol] focusResidue notice:', e);
    }
  }

  clearFlexibleHighlights() {
    if (!this.viewer) return;
    if (this.flexHighlights && this.flexHighlights.size > 0) {
      this.flexHighlights.forEach((item, key) => {
        if (item.surfaceId !== undefined && item.surfaceId !== null) {
          try { this.viewer.removeSurface(item.surfaceId); } catch (e) {}
        }
        if (item.shapes) {
          item.shapes.forEach(s => {
            try { this.viewer.removeShape(s); } catch (e) {}
          });
        }
        if (item.labels) {
          item.labels.forEach(l => {
            try { this.viewer.removeLabel(l); } catch (e) {}
          });
        }
        if (this.receptorModel && item.sel) {
          const isInteracting = this.lastInteractingResidues && this.lastInteractingResidues.has(key);
          if (isInteracting) {
            this.receptorModel.setStyle(item.sel, {
              cartoon: { color: '#38bdf8', opacity: 0.9 },
              stick: { radius: 0.22, colorscheme: 'amino' }
            });
            if (this.settings.showResidueLabels) {
              const atoms = this.receptorModel.selectedAtoms(item.sel);
              if (atoms && atoms.length > 0) {
                const ca = atoms.find(a => a.atom === 'CA') || atoms[0];
                const parts = key.split(':');
                const num = parts[1] || '';
                const displayRes = (item.resName ? item.resName + ' ' : '') + num;
                const baseLbl = this.viewer.addLabel(displayRes, {
                  position: { x: ca.x, y: ca.y, z: ca.z },
                  backgroundColor: 'rgba(15, 23, 42, 0.8)',
                  fontColor: '#38bdf8',
                  fontSize: 11,
                  borderThickness: 1,
                  borderColor: '#0284c7'
                });
                if (baseLbl) {
                  if (!this.residueLabels) this.residueLabels = new Map();
                  this.residueLabels.set(key, baseLbl);
                }
              }
            }
          } else {
            this.receptorModel.setStyle(item.sel, {
              cartoon: { colorscheme: 'chain', opacity: 0.85, thickness: 0.4 }
            });
          }
        }
      });
      this.flexHighlights.clear();
      this.viewer.render();
    }
  }

  resetCamera() {
    if (!this.viewer) return;
    if (this.ligandModel && this.receptorModel) {
      this.viewer.zoomTo();
    } else if (this.ligandModel) {
      this.viewer.zoomTo({ model: this.ligandModel }, 500);
    } else if (this.receptorModel) {
      this.viewer.zoomTo({ model: this.receptorModel }, 500);
    }
    this.viewer.render();
  }
}

window.MolecularViewer = MolecularViewer;
