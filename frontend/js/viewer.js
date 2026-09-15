// 3D Molecular Viewer controller using 3Dmol.js
// Ultra-lightweight WebGL implementation tuned for low-memory systems (2GB RAM safe)

class MolecularViewer {
  constructor(elementId) {
    this.elementId = elementId;
    this.viewer = null;
    this.receptorModel = null;
    this.ligandModel = null;
    this.surfaceObj = null;
    this.interactionShapes = [];
    this.measureMode = false;
    this.measureAtoms = [];
    this.measureShapes = [];
    
    // Default display settings
    this.settings = {
      proteinStyle: "cartoon",
      proteinColor: "spectrum", // 'spectrum', 'chain', 'white'
      ligandStyle: "stick",
      showSurface: false,
      surfaceType: "SES", // Solvent Excluded Surface
      surfaceOpacity: 0.5,
      showHbonds: true,
      showHydrophobic: false,
      showResidueLabels: true
    };

    this.init();
  }

  init() {
    const el = document.getElementById(this.elementId);
    if (!el) return;

    if (window.$3Dmol) {
      this.viewer = $3Dmol.createViewer(el, {
        backgroundColor: "0x0f172a",
        antialias: true
      });
      console.log("[3Dmol] Viewer initialized successfully.");
    } else {
      console.warn("[3Dmol] 3Dmol.js library not loaded yet.");
    }
  }

  clear() {
    if (!this.viewer) return;
    this.viewer.clear();
    this.receptorModel = null;
    this.ligandModel = null;
    this.surfaceObj = null;
    this.clearInteractions();
    this.clearMeasurements();
  }

  clearInteractions() {
    if (!this.viewer) return;
    this.interactionShapes.forEach(s => this.viewer.removeShape(s));
    this.interactionShapes = [];
    this.viewer.removeAllLabels();
  }

  clearMeasurements() {
    if (!this.viewer) return;
    this.measureShapes.forEach(s => this.viewer.removeShape(s));
    this.measureShapes = [];
    this.measureAtoms = [];
  }

  loadReceptor(pdbContent, format = "pdb") {
    if (!this.viewer || !pdbContent) return;
    
    // Remove previous receptor
    if (this.receptorModel) {
      this.viewer.removeModel(this.receptorModel);
    }

    this.receptorModel = this.viewer.addModel(pdbContent, format);
    this.applyProteinStyle();
    this.viewer.zoomTo();
    this.viewer.render();
  }

  applyProteinStyle() {
    if (!this.receptorModel) return;

    let colorScheme = {};
    if (this.settings.proteinColor === "spectrum") {
      colorScheme = { color: "spectrum" };
    } else if (this.settings.proteinColor === "chain") {
      colorScheme = { colorscheme: "chain" };
    } else {
      colorScheme = { color: "lightgray" };
    }

    if (this.settings.proteinStyle === "cartoon") {
      this.receptorModel.setStyle({}, { cartoon: { ...colorScheme, opacity: 0.85, thickness: 0.4 } });
    } else if (this.settings.proteinStyle === "stick") {
      this.receptorModel.setStyle({}, { stick: { radius: 0.2, colorscheme: "Jmol" } });
    } else if (this.settings.proteinStyle === "sphere") {
      this.receptorModel.setStyle({}, { sphere: { radius: 0.8, colorscheme: "Jmol" } });
    } else if (this.settings.proteinStyle === "ribbon") {
      this.receptorModel.setStyle({}, { ribbon: { ...colorScheme, opacity: 0.9 } });
    }
    
    this.viewer.render();
  }

  loadLigand(pdbOrPdbqtContent, format = "pdb") {
    if (!this.viewer || !pdbOrPdbqtContent) return;

    if (this.ligandModel) {
      this.viewer.removeModel(this.ligandModel);
    }

    this.ligandModel = this.viewer.addModel(pdbOrPdbqtContent, format);
    
    // Style ligand as distinct, high-contrast stick or ball & stick
    this.ligandModel.setStyle({}, {
      stick: { radius: 0.28, colorscheme: "cyanCarbon" },
      sphere: { radius: 0.35, colorscheme: "cyanCarbon" }
    });

    // Center view smoothly on the docked ligand
    this.viewer.zoomTo({ model: this.ligandModel }, 600);
    this.viewer.render();

    if (this.settings.showSurface) {
      this.updateSurface();
    }
  }

  renderInteractions(interactions) {
    if (!this.viewer || !interactions) return;
    this.clearInteractions();

    const { hydrogen_bonds = [], hydrophobic_contacts = [], interacting_residues = [] } = interactions;

    // 1. Highlight interacting amino acid residues on receptor
    if (this.receptorModel && interacting_residues.length > 0) {
      interacting_residues.forEach(resStr => {
        // resStr format e.g. "ARG 120:A"
        const parts = resStr.split(" ");
        const resName = parts[0];
        const numChain = parts[1] || "";
        const [resNumStr, chain] = numChain.split(":");
        const resNum = parseInt(resNumStr);

        const sel = { resi: resNum };
        if (chain) sel.chain = chain;

        // Display sidechain stick
        this.receptorModel.setStyle(sel, {
          cartoon: { color: "#38bdf8", opacity: 0.9 },
          stick: { radius: 0.22, colorscheme: "amino" }
        });

        // Add residue label
        if (this.settings.showResidueLabels) {
          const atoms = this.receptorModel.selectedAtoms(sel);
          if (atoms && atoms.length > 0) {
            const ca = atoms.find(a => a.atom === "CA") || atoms[0];
            this.viewer.addLabel(`${resName}${resNum}`, {
              position: { x: ca.x, y: ca.y, z: ca.z },
              backgroundColor: "rgba(15, 23, 42, 0.75)",
              fontColor: "#38bdf8",
              fontSize: 11,
              borderThickness: 1,
              borderColor: "#0284c7"
            });
          }
        }
      });
    }

    // 2. Render Hydrogen Bonds as yellow dashed cylinders
    if (this.settings.showHbonds && hydrogen_bonds.length > 0) {
      hydrogen_bonds.forEach(hb => {
        const [lx, ly, lz] = hb.start_coord;
        const [rx, ry, rz] = hb.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.08,
          dashed: true,
          color: "#facc15", // Bright yellow
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);

        // Distance label at midpoint
        const mid = { x: (lx + rx) / 2, y: (ly + ry) / 2, z: (lz + rz) / 2 };
        this.viewer.addLabel(`${hb.distance} Å`, {
          position: mid,
          backgroundColor: "rgba(17, 24, 39, 0.8)",
          fontColor: "#facc15",
          fontSize: 10
        });
      });
    }

    // 3. Render Hydrophobic contacts if enabled (cyan dashed)
    if (this.settings.showHydrophobic && hydrophobic_contacts.length > 0) {
      hydrophobic_contacts.forEach(hp => {
        const [lx, ly, lz] = hp.start_coord;
        const [rx, ry, rz] = hp.end_coord;

        const cyl = this.viewer.addCylinder({
          start: { x: lx, y: ly, z: lz },
          end: { x: rx, y: ry, z: rz },
          radius: 0.06,
          dashed: true,
          color: "#38bdf8",
          fromCap: 1,
          toCap: 1
        });
        this.interactionShapes.push(cyl);
      });
    }

    this.viewer.render();
  }

  toggleSurface(show) {
    this.settings.showSurface = show;
    this.updateSurface();
  }

  updateSurface() {
    if (!this.viewer) return;

    if (this.surfaceObj) {
      this.viewer.removeSurface(this.surfaceObj);
      this.surfaceObj = null;
    }

    if (this.settings.showSurface) {
      // Create pocket surface around ligand
      const sel = this.ligandModel ? { model: this.receptorModel, within: { distance: 6.0, sel: { model: this.ligandModel } } } : { model: this.receptorModel };
      this.surfaceObj = this.viewer.addSurface($3Dmol.SurfaceType.MS, {
        opacity: this.settings.surfaceOpacity,
        color: "#64748b"
      }, sel);
    }
    this.viewer.render();
  }

  enableMeasurementMode(enable) {
    this.measureMode = enable;
    this.clearMeasurements();
    if (!this.viewer) return;

    if (enable) {
      this.viewer.setClickable({}, true, (atom, viewer, event, container) => {
        if (!this.measureMode) return;
        this.measureAtoms.push(atom);
        
        if (this.measureAtoms.length === 1) {
          console.log("Atom 1 selected for measurement:", atom.atom, atom.resi);
        } else if (this.measureAtoms.length === 2) {
          const a1 = this.measureAtoms[0];
          const a2 = this.measureAtoms[1];
          const dx = a1.x - a2.x;
          const dy = a1.y - a2.y;
          const dz = a1.z - a2.z;
          const dist = Math.sqrt(dx*dx + dy*dy + dz*dz).toFixed(2);

          const cyl = this.viewer.addCylinder({
            start: { x: a1.x, y: a1.y, z: a1.z },
            end: { x: a2.x, y: a2.y, z: a2.z },
            radius: 0.1,
            color: "#f43f5e" // Rose
          });
          this.measureShapes.push(cyl);

          const mid = { x: (a1.x + a2.x)/2, y: (a1.y + a2.y)/2, z: (a1.z + a2.z)/2 };
          this.viewer.addLabel(`${dist} Å`, {
            position: mid,
            backgroundColor: "rgba(244, 63, 94, 0.9)",
            fontColor: "#ffffff",
            fontSize: 12
          });

          this.measureAtoms = [];
          this.viewer.render();
        }
      });
    } else {
      this.viewer.setClickable({}, false, null);
    }
  }

  resetCamera() {
    if (!this.viewer) return;
    if (this.ligandModel) {
      this.viewer.zoomTo({ model: this.ligandModel }, 600);
    } else if (this.receptorModel) {
      this.viewer.zoomTo({ model: this.receptorModel }, 600);
    }
  }
}

window.MolecularViewer = MolecularViewer;
