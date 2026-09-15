// 3D Molecular Viewer controller using 3Dmol.js
// High-performance WebGL implementation tuned for low-memory systems (2GB RAM safe)

class MolecularViewer {
  constructor(elementId) {
    this.elementId = elementId;
    this.viewer = null;
    this.receptorModel = null;
    this.ligandModel = null;
    this.surfaceObj = null;
    this.interactionShapes = [];
    this.gridBoxShapes = [];
    this.measureMode = false;
    this.measureAtoms = [];
    this.measureShapes = [];
    this.measureMarkers = [];
    this.isSpinning = false;
    
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
      showHydrophobic: false,
      showResidueLabels: true
    };

    // Surface tracking — prevents accumulation bug
    this._surfaceId = null;
    this._surfaceUpdating = false;

    this.init();
  }

  init() {
    const el = document.getElementById(this.elementId);
    if (!el) return;

    const M = window["$" + "3Dmol"];
    if (M) {
      this.viewer = M.createViewer(el, {
        backgroundColor: '0x0f172a',
        antialias: true
      });
      console.log('[3Dmol] Viewer initialized successfully.');
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
    this.clearInteractions();
    this.clearMeasurements();
    this.clearGridBox();
  }

  clearInteractions() {
    if (!this.viewer) return;
    this.interactionShapes.forEach(s => {
      try { this.viewer.removeShape(s); } catch (e) {}
    });
    this.interactionShapes = [];
    this.viewer.removeAllLabels();
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
    this.measureAtoms = [];
    this.viewer.removeAllLabels();
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
      this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.85, thickness: 0.4 } });
    } else if (style === 'stick') {
      this.receptorModel.setStyle({}, { stick: { radius: 0.18, colorscheme: 'Jmol' } });
    } else if (style === 'sphere') {
      this.receptorModel.setStyle({}, { sphere: { radius: 0.75, colorscheme: 'Jmol' } });
    } else if (style === 'ribbon') {
      this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.9, style: 'trace', thickness: 0.6 } });
    } else if (style === 'hidden') {
      this.receptorModel.setStyle({}, {});
    }
    
    this.viewer.render();
  }

  loadLigand(pdbOrPdbqtContent, format = 'pdb') {
    if (!this.viewer || !pdbOrPdbqtContent) return;

    if (this.ligandModel) {
      this.viewer.removeModel(this.ligandModel);
      this.ligandModel = null;
    }

    try {
      this.ligandModel = this.viewer.addModel(pdbOrPdbqtContent, format);
      this.applyLigandStyle();

      this.viewer.zoomTo({ model: this.ligandModel }, 600);
      this.viewer.render();

      if (this.settings.showSurface) {
        this.updateSurface();
      }
    } catch (e) {
      console.error('[3Dmol] Failed to load ligand:', e);
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
    if (!this.viewer || !interactions) return;
    this.clearInteractions();

    const { hydrogen_bonds = [], hydrophobic_contacts = [], interacting_residues = [] } = interactions;

    if (this.receptorModel && interacting_residues.length > 0) {
      interacting_residues.forEach(resStr => {
        const parts = resStr.split(' ');
        const resName = parts[0];
        const numChain = parts[1] || '';
        const [resNumStr, chain] = numChain.split(':');
        const resNum = parseInt(resNumStr);

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
            this.viewer.addLabel(resName + resNum, {
              position: { x: ca.x, y: ca.y, z: ca.z },
              backgroundColor: 'rgba(15, 23, 42, 0.8)',
              fontColor: '#38bdf8',
              fontSize: 11,
              borderThickness: 1,
              borderColor: '#0284c7'
            });
          }
        }
      });
    }

    if (this.settings.showHbonds && hydrogen_bonds.length > 0) {
      hydrogen_bonds.forEach(hb => {
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

        const mid = { x: (lx + rx) / 2, y: (ly + ry) / 2, z: (lz + rz) / 2 };
        this.viewer.addLabel(hb.distance + ' A', {
          position: mid,
          backgroundColor: 'rgba(17, 24, 39, 0.85)',
          fontColor: '#facc15',
          fontSize: 10
        });
      });
    }

    if (this.settings.showHydrophobic && hydrophobic_contacts.length > 0) {
      hydrophobic_contacts.forEach(hp => {
        const [lx, ly, lz] = hp.start_coord;
        const [rx, ry, rz] = hp.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.06,
          dashed: true,
          color: '#38bdf8',
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);
      });
    }

    this.viewer.render();
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
      this.viewer.setClickable({}, false, null);
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

    this.viewer.setClickable({}, true, (atom) => {
      if (!this.measureMode || !atom) return;

      this.measureAtoms.push(atom);

      if (this.measureAtoms.length === 1) {
        const marker = this.viewer.addSphere({
          center: { x: atom.x, y: atom.y, z: atom.z },
          radius: 0.45,
          color: '#f43f5e'
        });
        this.measureMarkers.push(marker);
        this.viewer.render();

        const atomDesc = (atom.atom || atom.elem || 'Atom') + ' (' + (atom.resn || '') + (atom.resi || '') + ')';
        if (hudText) hudText.textContent = 'Point 1: [' + atomDesc + '] selected. Click 2nd atom...';
      } else if (this.measureAtoms.length >= 2) {
        const a1 = this.measureAtoms[0];
        const a2 = this.measureAtoms[1];
        const dx = a1.x - a2.x;
        const dy = a1.y - a2.y;
        const dz = a1.z - a2.z;
        const dist = Math.sqrt(dx*dx + dy*dy + dz*dz).toFixed(2);

        const cyl = this.viewer.addCylinder({
          start: { x: a1.x, y: a1.y, z: a1.z },
          end: { x: a2.x, y: a2.y, z: a2.z },
          radius: 0.12,
          color: '#f43f5e'
        });
        this.measureShapes.push(cyl);

        const mid = { x: (a1.x + a2.x) / 2, y: (a1.y + a2.y) / 2, z: (a1.z + a2.z) / 2 };
        this.viewer.addLabel(dist + ' A', {
          position: mid,
          backgroundColor: 'rgba(244, 63, 94, 0.95)',
          fontColor: '#ffffff',
          fontSize: 12,
          borderThickness: 1,
          borderColor: '#ffffff'
        });

        this.measureMarkers.forEach(m => {
          try { this.viewer.removeShape(m); } catch (e) {}
        });
        this.measureMarkers = [];
        this.measureAtoms = [];

        const desc1 = a1.atom || a1.elem || 'A1';
        const desc2 = a2.atom || a2.elem || 'A2';
        if (hudText) hudText.textContent = 'Measured: ' + dist + ' A (' + desc1 + ' - ' + desc2 + '). Click another pair to measure.';
        this.viewer.render();
      }
    });
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
