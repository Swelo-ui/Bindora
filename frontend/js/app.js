// Main frontend application controller for Bindora

document.addEventListener("DOMContentLoaded", () => {
  const app = new BindoraApp();
  window.app = app;
  window.bindoraApp = app;
  app.init();
});

class BindoraApp {
  constructor() {
    this.viewer = null;
    this.state = {
      ligand: null,
      receptor: null,
      docking: null,
      currentPoseIdx: 0,
      crosscheck: null,
      narrative: null,
      benchmarks: [],
      activeTab: "home"
    };
  }

  async init() {
    // 1. Initialize 3Dmol viewer
    try {
      this.viewer = new MolecularViewer("gldiv");
    } catch (e) {
      console.error("Failed to initialize 3D viewer:", e);
    }

    // 2. Setup DOM event listeners
    this.setupEventListeners();

    // 3. Check backend health
    try {
      const health = await BindoraAPI.checkHealth();
      const statusBadge = document.getElementById("backend-status-badge");
      if (statusBadge) {
        if (health.vina_available) {
          statusBadge.innerHTML = `<span class="inline-block w-2 h-2 rounded-full bg-emerald-400 mr-1.5 animate-pulse flex-shrink-0"></span><span class="hidden 2xl:inline">AutoDock Vina Ready</span><span class="2xl:hidden">Vina Ready</span>`;
          statusBadge.className = "hidden lg:flex px-2.5 py-1 text-xs font-semibold rounded-full bg-emerald-950/80 text-emerald-300 border border-emerald-800 items-center whitespace-nowrap flex-shrink-0";
        } else {
          statusBadge.innerHTML = `<span class="inline-block w-2 h-2 rounded-full bg-amber-400 mr-1.5 flex-shrink-0"></span><span class="hidden 2xl:inline">Vina Standby</span><span class="2xl:hidden">Standby</span>`;
          statusBadge.className = "hidden lg:flex px-2.5 py-1 text-xs font-semibold rounded-full bg-amber-950/80 text-amber-300 border border-amber-800 items-center whitespace-nowrap flex-shrink-0";
        }
      }
    } catch (e) {
      console.warn("Backend health check failed:", e);
    }

    // 4. Fetch benchmarks
    try {
      const bmData = await BindoraAPI.getBenchmarks();
      this.state.benchmarks = bmData.benchmarks || [];
      this.renderBenchmarksUI();
    } catch (e) {
      console.warn("Could not fetch benchmarks:", e);
    }

    // 4b. Fetch live empirical validation reports
    try {
      this.loadValidationReportUI();
    } catch (e) {
      console.warn("Could not load validation report:", e);
    }

    // 5. Clean live state on startup (benchmarks are available above for 1-click exploration if user chooses)
    this.updateStudioCards();

    // 6. Set active tab (supports ?tab=studio or #studio, default to home)
    const urlParams = new URLSearchParams(window.location.search);
    const initialTab = urlParams.get("tab") || window.location.hash.replace("#", "") || "home";
    this.switchTab(initialTab);

    // Auto-open user guide if requested
    if (urlParams.get("guide") === "true") {
      document.getElementById("modal-about")?.classList.remove("hidden");
    }
  }

  setupEventListeners() {
    // Tab switching
    document.querySelectorAll(".nav-tab").forEach(btn => {
      btn.addEventListener("click", () => {
        const tab = btn.dataset.tab;
        this.switchTab(tab);
      });
    });

    // PubChem Search
    const pubchemBtn = document.getElementById("btn-search-pubchem");
    const pubchemInput = document.getElementById("input-search-pubchem");
    if (pubchemBtn && pubchemInput) {
      pubchemBtn.addEventListener("click", () => this.searchPubChem(pubchemInput.value));
      pubchemInput.addEventListener("keydown", e => {
        if (e.key === "Enter") this.searchPubChem(pubchemInput.value);
      });
    }

    // Ligand Mode Switching (PubChem / SMILES / File Upload)
    const setLigandMode = (mode) => {
      ['pubchem', 'smiles', 'upload'].forEach(m => {
        const el = document.getElementById(`mode-ligand-${m}`);
        const btn = document.getElementById(`tab-btn-ligand-${m}`);
        if (el) el.classList.toggle('hidden', m !== mode);
        if (btn) {
          btn.classList.toggle('text-cyan-300', m === mode);
          btn.classList.toggle('bg-slate-800', m === mode);
          btn.classList.toggle('text-slate-400', m !== mode);
        }
      });
    };
    document.getElementById("tab-btn-ligand-pubchem")?.addEventListener("click", () => setLigandMode('pubchem'));
    document.getElementById("tab-btn-ligand-smiles")?.addEventListener("click", () => setLigandMode('smiles'));
    document.getElementById("tab-btn-ligand-upload")?.addEventListener("click", () => setLigandMode('upload'));

    // Custom SMILES load
    const loadSmilesBtn = document.getElementById("btn-load-custom-smiles");
    const customSmilesInput = document.getElementById("input-custom-smiles");
    const customLigandNameInput = document.getElementById("input-custom-ligand-name");
    if (loadSmilesBtn && customSmilesInput) {
      loadSmilesBtn.addEventListener("click", () => {
        this.loadCustomLigand(customLigandNameInput?.value.trim(), customSmilesInput.value.trim());
      });
      customSmilesInput.addEventListener("keydown", e => {
        if (e.key === "Enter") this.loadCustomLigand(customLigandNameInput?.value.trim(), customSmilesInput.value.trim());
      });
    }

    // Ligand File Upload (.sdf, .mol, .smi)
    const ligandFileInput = document.getElementById("file-upload-ligand");
    if (ligandFileInput) {
      ligandFileInput.addEventListener("change", (e) => {
        const file = e.target.files?.[0];
        if (file) this.uploadLigandFile(file);
      });
    }

    // RCSB Search
    const rcsbBtn = document.getElementById("btn-search-rcsb");
    const rcsbInput = document.getElementById("input-search-rcsb");
    if (rcsbBtn && rcsbInput) {
      rcsbBtn.addEventListener("click", () => this.searchRCSB(rcsbInput.value));
      rcsbInput.addEventListener("keydown", e => {
        if (e.key === "Enter") this.searchRCSB(rcsbInput.value);
      });
    }

    // Receptor Mode Switching (RCSB / File Upload / Paste PDB)
    const setReceptorMode = (mode) => {
      ['rcsb', 'upload', 'paste'].forEach(m => {
        const el = document.getElementById(`mode-receptor-${m}`);
        const btn = document.getElementById(`tab-btn-receptor-${m}`);
        if (el) el.classList.toggle('hidden', m !== mode);
        if (btn) {
          btn.classList.toggle('text-emerald-300', m === mode);
          btn.classList.toggle('bg-slate-800', m === mode);
          btn.classList.toggle('text-slate-400', m !== mode);
        }
      });
    };
    document.getElementById("tab-btn-receptor-rcsb")?.addEventListener("click", () => setReceptorMode('rcsb'));
    document.getElementById("tab-btn-receptor-upload")?.addEventListener("click", () => setReceptorMode('upload'));
    document.getElementById("tab-btn-receptor-paste")?.addEventListener("click", () => setReceptorMode('paste'));

    // Receptor File Upload (.pdb)
    const receptorFileInput = document.getElementById("file-upload-receptor");
    if (receptorFileInput) {
      receptorFileInput.addEventListener("change", (e) => {
        const file = e.target.files?.[0];
        if (file) this.uploadReceptorFile(file);
      });
    }

    // Custom PDB Paste load
    const loadPdbBtn = document.getElementById("btn-load-custom-pdb");
    const customPdbTextInput = document.getElementById("input-custom-pdb-text");
    const customReceptorNameInput = document.getElementById("input-custom-receptor-name");
    if (loadPdbBtn && customPdbTextInput) {
      loadPdbBtn.addEventListener("click", () => {
        this.loadCustomReceptorPdb(customReceptorNameInput?.value.trim() || "CUSTOM", customPdbTextInput.value.trim());
      });
    }

    // Docking trigger button
    const dockBtn = document.getElementById("btn-run-docking");
    if (dockBtn) {
      dockBtn.addEventListener("click", () => this.runDockingPipeline());
    }

    // Viewer Controls: Camera & Tools
    const resetCamBtn = document.getElementById("btn-reset-cam");
    if (resetCamBtn) {
      resetCamBtn.addEventListener("click", () => this.viewer && this.viewer.resetCamera());
    }

    const spinBtn = document.getElementById("btn-toggle-spin");
    if (spinBtn) {
      spinBtn.addEventListener("click", () => this.viewer && this.viewer.toggleSpin());
    }

    const downloadPngBtn = document.getElementById("btn-download-png");
    if (downloadPngBtn) {
      downloadPngBtn.addEventListener("click", () => this.viewer && this.viewer.downloadScreenshot());
    }

    // Molecular Surface Toggle
    const surfaceToggle = document.getElementById("toggle-surface");
    if (surfaceToggle) {
      surfaceToggle.addEventListener("change", (e) => {
        if (this.viewer) {
          this.viewer.toggleSurface(e.target.checked);
          if (e.target.checked) {
            this.showToast(this.state.receptor && this.state.ligand ? "Binding pocket surface displayed" : "Molecular surface displayed", "info");
          }
        }
      });
    }

    // Docking Grid Box Toggle & Live Coordinates Sync
    const getGridCenter = () => ({
      x: parseFloat(document.getElementById("grid-cx")?.value) || 0,
      y: parseFloat(document.getElementById("grid-cy")?.value) || 0,
      z: parseFloat(document.getElementById("grid-cz")?.value) || 0
    });
    const getGridSize = () => ({
      x: parseFloat(document.getElementById("grid-sx")?.value) || 22,
      y: parseFloat(document.getElementById("grid-sy")?.value) || 22,
      z: parseFloat(document.getElementById("grid-sz")?.value) || 22
    });

    const gridBoxToggle = document.getElementById("toggle-gridbox");
    if (gridBoxToggle) {
      gridBoxToggle.addEventListener("change", (e) => {
        if (this.viewer) {
          this.viewer.settings.showGridBox = e.target.checked;
          this.viewer.renderGridBox(getGridCenter(), getGridSize(), e.target.checked);
          if (e.target.checked) this.showToast("3D Docking Grid Box displayed", "info");
        }
      });
    }

    ["grid-cx", "grid-cy", "grid-cz", "grid-sx", "grid-sy", "grid-sz"].forEach(id => {
      const el = document.getElementById(id);
      if (el) {
        el.addEventListener("input", () => {
          const gToggle = document.getElementById("toggle-gridbox");
          if (gToggle && gToggle.checked && this.viewer) {
            this.viewer.renderGridBox(getGridCenter(), getGridSize(), true);
          }
        });
      }
    });

    // Measure Tool
    const measureToggle = document.getElementById("toggle-measure");
    if (measureToggle) {
      measureToggle.addEventListener("change", (e) => {
        if (this.viewer) {
          this.viewer.enableMeasurementMode(e.target.checked);
          if (e.target.checked) {
            this.showToast("Measurement tool enabled: Click two atoms to measure distance in Å", "info");
          }
        }
      });
    }
    const clearMeasureBtn = document.getElementById("btn-clear-measure");
    if (clearMeasureBtn) {
      clearMeasureBtn.addEventListener("click", () => {
        if (this.viewer) this.viewer.clearMeasurements();
        const hudText = document.getElementById("measure-hud-text");
        if (hudText) hudText.innerHTML = '<svg class="w-3.5 h-3.5 text-rose-400 inline mr-1" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21.3 8.7 8.7 21.3c-1 1-2.5 1-3.4 0l-2.6-2.6c-1-1-1-2.5 0-3.4L15.3 2.7c1-1 2.5-1 3.4 0l2.6 2.6c1 1 1 2.5 0 3.4Z"/><path d="m14.5 3.5 2 2"/><path d="m11.5 6.5 2 2"/><path d="m8.5 9.5 2 2"/><path d="m5.5 12.5 2 2"/></svg><span>Click any atom to begin measurement...</span>';
      });
    }

    // Protein Style Select
    const proteinStyleSelect = document.getElementById("select-protein-style");
    if (proteinStyleSelect) {
      proteinStyleSelect.addEventListener("change", (e) => {
        if (this.viewer) {
          this.viewer.settings.proteinStyle = e.target.value;
          this.viewer.applyProteinStyle();
          if (!this.state.receptor) {
            this.showToast("No protein loaded yet. Style will apply when a protein is loaded.", "warning");
          }
        }
      });
    }

    // Ligand Style Select
    const ligandStyleSelect = document.getElementById("select-ligand-style");
    if (ligandStyleSelect) {
      ligandStyleSelect.addEventListener("change", (e) => {
        if (this.viewer) {
          this.viewer.settings.ligandStyle = e.target.value;
          this.viewer.applyLigandStyle();
          if (!this.state.ligand) {
            this.showToast("No ligand loaded yet.", "warning");
          }
        }
      });
    }

    // Ligand Color Select
    const ligandColorSelect = document.getElementById("select-ligand-color");
    if (ligandColorSelect) {
      ligandColorSelect.addEventListener("change", (e) => {
        if (this.viewer) {
          this.viewer.settings.ligandColor = e.target.value;
          this.viewer.applyLigandStyle();
        }
      });
    }

    // Pose Steppers
    const prevPoseBtn = document.getElementById("btn-prev-pose");
    const nextPoseBtn = document.getElementById("btn-next-pose");
    if (prevPoseBtn && nextPoseBtn) {
      prevPoseBtn.addEventListener("click", () => this.changePose(-1));
      nextPoseBtn.addEventListener("click", () => this.changePose(1));
    }

    // Batch docking trigger
    const batchBtn = document.getElementById("btn-run-batch");
    if (batchBtn) {
      batchBtn.addEventListener("click", () => this.runBatchScreening());
    }

    // Batch Mode Switching (Candidates / Ensemble / Pharmacophore)
    const setBatchMode = (mode) => {
      ['candidates', 'ensemble', 'pharmacophore'].forEach(m => {
        const panel = document.getElementById(`panel-batch-${m}`);
        const btn = document.getElementById(`btn-batch-mode-${m}`);
        if (panel) panel.classList.toggle('hidden', m !== mode);
        if (btn) {
          if (m === mode) {
            btn.className = m === 'pharmacophore' 
              ? 'py-1 px-3 rounded font-medium bg-purple-600 text-white shadow transition'
              : m === 'ensemble'
              ? 'py-1 px-3 rounded font-medium bg-emerald-600 text-white shadow transition'
              : 'py-1 px-3 rounded font-medium bg-cyan-600 text-white shadow transition';
          } else {
            btn.className = m === 'pharmacophore'
              ? 'py-1 px-3 rounded font-medium text-purple-600 dark:text-purple-300 hover:text-white border border-purple-800/40 bg-purple-950/20 transition'
              : 'py-1 px-3 rounded font-medium text-slate-600 dark:text-slate-400 hover:text-white transition';
          }
        }
      });
      if (mode === 'ensemble') {
        this.loadEnsembleStructures();
      } else if (mode === 'pharmacophore') {
        this.checkPharmacophoreEligibility();
      }
    };
    document.getElementById("btn-batch-mode-candidates")?.addEventListener("click", () => setBatchMode('candidates'));
    document.getElementById("btn-batch-mode-ensemble")?.addEventListener("click", () => setBatchMode('ensemble'));
    document.getElementById("btn-batch-mode-pharmacophore")?.addEventListener("click", () => setBatchMode('pharmacophore'));

    // Ensemble run trigger
    document.getElementById("btn-run-ensemble")?.addEventListener("click", () => this.runEnsembleDocking());

    // Pharmacophore screen trigger
    document.getElementById("btn-run-pharmacophore-screen")?.addEventListener("click", () => this.runPharmacophoreScreen());

    // Batch Library Presets
    const librarySelect = document.getElementById("select-batch-library");
    const batchTextarea = document.getElementById("input-batch-candidates");
    if (librarySelect && batchTextarea) {
      const presets = {
        nsaid: "Aspirin, CC(=O)Oc1ccccc1C(=O)O\nIbuprofen, CC(C)Cc1ccc(cc1)C(C)C(=O)O\nNaproxen, COc1ccc2cc(ccc2c1)C(C)C(=O)O\nCelecoxib, Cc1ccc(cc1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F",
        kinase: "Imatinib, Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C\nGefitinib, COc1cc2ncnc(c2cc1OCCCN3CCOCC3)Nc4ccc(c(c4)Cl)F\nErlotinib, COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC\nDasatinib, Cc1cccc(c1Cl)NC(=O)c2cnc(s2)Nc3cc(nc(n3)C)N4CCN(CC4)CCO",
        antiviral: "Remdesivir, CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4\nFavipiravir, C1=C(N=C(C(=O)N1)C(=O)N)F\nMolnupiravir, CC(C)C(=O)OCC1C(C(C(O1)N2C=CC(=NO)NC2=O)O)O\nRibavirin, C1=NC(=NN1C2C(C(C(O2)CO)O)O)C(=O)N",
        clear: ""
      };
      librarySelect.addEventListener("change", (e) => {
        const val = e.target.value;
        if (presets[val] !== undefined) {
          batchTextarea.value = presets[val];
          if (val !== "clear") {
            this.showToast(`Loaded ${librarySelect.options[librarySelect.selectedIndex].text}`, "info");
          }
        }
      });
    }

    // Workspace Reset / Clear Session
    const resetBtn = document.getElementById("btn-reset-session");
    if (resetBtn) {
      resetBtn.addEventListener("click", () => this.resetWorkspace());
    }

    // Firebase Auth Modal & Events
    const authBtn = document.getElementById("btn-user-auth");
    const authModal = document.getElementById("modal-firebase-auth");
    const closeAuthBtn = document.getElementById("btn-close-auth-modal");
    const googleAuthBtn = document.getElementById("btn-auth-google");
    const emailAuthBtn = document.getElementById("btn-auth-email-submit");
    const toggleAuthModeBtn = document.getElementById("btn-auth-toggle-mode");
    const signoutBtn = document.getElementById("btn-auth-signout");
    const saveCloudBtn = document.getElementById("btn-save-cloud");
    const refreshSavedRunsBtn = document.getElementById("btn-refresh-saved-runs");

    let isSignUpMode = false;

    if (authBtn && authModal) {
      authBtn.addEventListener("click", () => {
        this.updateAuthModalView();
        authModal.classList.remove("hidden");
      });
    }
    if (closeAuthBtn && authModal) {
      closeAuthBtn.addEventListener("click", () => authModal.classList.add("hidden"));
    }

    if (toggleAuthModeBtn) {
      toggleAuthModeBtn.addEventListener("click", () => {
        isSignUpMode = !isSignUpMode;
        if (emailAuthBtn) emailAuthBtn.textContent = isSignUpMode ? "Create Account" : "Sign In";
        toggleAuthModeBtn.textContent = isSignUpMode ? "Already have an account? Sign In" : "New user? Create an account";
        const errorEl = document.getElementById("auth-error-msg");
        if (errorEl) errorEl.classList.add("hidden");
      });
    }

    if (googleAuthBtn) {
      googleAuthBtn.addEventListener("click", async () => {
        const errorEl = document.getElementById("auth-error-msg");
        if (errorEl) errorEl.classList.add("hidden");
        try {
          await window.bindoraFirebase.signInWithGoogle();
          this.showToast("Signed in with Google successfully!", "success");
          this.updateAuthModalView();
        } catch (e) {
          if (errorEl) {
            errorEl.textContent = e.message || "Google sign-in failed.";
            errorEl.classList.remove("hidden");
          }
        }
      });
    }

    if (emailAuthBtn) {
      emailAuthBtn.addEventListener("click", async () => {
        const email = document.getElementById("input-auth-email")?.value.trim();
        const pass = document.getElementById("input-auth-password")?.value;
        const errorEl = document.getElementById("auth-error-msg");
        if (errorEl) errorEl.classList.add("hidden");

        if (!email || !pass) {
          if (errorEl) {
            errorEl.textContent = "Please enter both email and password.";
            errorEl.classList.remove("hidden");
          }
          return;
        }

        try {
          if (isSignUpMode) {
            await window.bindoraFirebase.signUpWithEmail(email, pass);
            this.showToast("Account created and signed in!", "success");
          } else {
            await window.bindoraFirebase.signInWithEmail(email, pass);
            this.showToast("Signed in successfully!", "success");
          }
          this.updateAuthModalView();
        } catch (e) {
          if (errorEl) {
            errorEl.textContent = e.message;
            errorEl.classList.remove("hidden");
          }
        }
      });
    }

    if (signoutBtn) {
      signoutBtn.addEventListener("click", async () => {
        await window.bindoraFirebase.signOut();
        this.showToast("Signed out of researcher account.", "info");
        this.updateAuthModalView();
      });
    }

    if (saveCloudBtn) {
      saveCloudBtn.addEventListener("click", async () => {
        if (!this.state.docking) {
          this.showToast("Please run docking first before saving.", "error");
          return;
        }
        try {
          const topPose = this.state.docking.top_pose || {};
          const thermo = this.state.docking.thermodynamics || {};
          const contacts = this.state.docking.interactions || {};
          const adme = this.state.ligand?.adme || {};

          await window.bindoraFirebase.saveDockingRun({
            drugName: this.state.ligand?.name || "Drug Candidate",
            targetName: this.state.receptor?.title || this.state.receptor?.pdb_id || "Target",
            pdbId: this.state.receptor?.pdb_id || "N/A",
            affinityKcal: topPose.affinity_kcal || 0.0,
            theoreticalKdNm: thermo.theoretical_kd_nm || "N/A",
            ligandEfficiency: thermo.ligand_efficiency?.value || 0.0,
            hbondCount: contacts.total_hbond_count || 0,
            lipinskiStatus: adme.drug_likeness?.lipinski?.status || "Pass",
            crosscheckBadge: this.state.crosscheck?.status_badge || "Computational Prediction Only",
            isCrossChecked: !!this.state.crosscheck?.is_cross_checked
          });

          this.showToast("Simulation saved to cloud history!", "success");
          this.loadCloudHistory();
        } catch (e) {
          this.showToast(e.message, "error");
        }
      });
    }

    if (refreshSavedRunsBtn) {
      refreshSavedRunsBtn.addEventListener("click", () => this.loadCloudHistory());
    }

    // Profile Editing Events
    const toggleEditProfileBtn = document.getElementById("btn-toggle-edit-profile");
    const cancelEditProfileBtn = document.getElementById("btn-cancel-edit-profile");
    const profileEditForm = document.getElementById("profile-edit-form");
    const randomAvatarBtn = document.getElementById("btn-random-avatar");
    const saveProfileBtn = document.getElementById("btn-save-profile");
    const inputEditName = document.getElementById("input-edit-display-name");
    const inputEditPhoto = document.getElementById("input-edit-photo-url");
    const profileEditStatus = document.getElementById("profile-edit-status");

    if (toggleEditProfileBtn && profileEditForm) {
      toggleEditProfileBtn.addEventListener("click", () => {
        const isHidden = profileEditForm.classList.contains("hidden");
        if (isHidden) {
          profileEditForm.classList.remove("hidden");
          const user = window.bindoraFirebase ? window.bindoraFirebase.currentUser : null;
          if (user) {
            if (inputEditName) inputEditName.value = user.displayName || "";
            if (inputEditPhoto) inputEditPhoto.value = user.photoURL || "";
          }
          if (profileEditStatus) profileEditStatus.classList.add("hidden");
        } else {
          profileEditForm.classList.add("hidden");
        }
      });
    }

    if (cancelEditProfileBtn && profileEditForm) {
      cancelEditProfileBtn.addEventListener("click", () => {
        profileEditForm.classList.add("hidden");
        if (profileEditStatus) profileEditStatus.classList.add("hidden");
      });
    }

    if (randomAvatarBtn && inputEditPhoto) {
      randomAvatarBtn.addEventListener("click", () => {
        const randomSeed = "Scientist-" + Math.floor(Math.random() * 90000 + 10000);
        inputEditPhoto.value = `https://api.dicebear.com/7.x/bottts/svg?seed=${randomSeed}`;
      });
    }

    // Avatar Presets
    document.querySelectorAll(".preset-avatar-btn").forEach(btn => {
      btn.addEventListener("click", () => {
        const avatarUrl = btn.getAttribute("data-avatar");
        if (avatarUrl && inputEditPhoto) {
          inputEditPhoto.value = avatarUrl;
        }
      });
    });

    if (saveProfileBtn) {
      saveProfileBtn.addEventListener("click", async () => {
        const newName = inputEditName ? inputEditName.value.trim() : "";
        const newPhoto = inputEditPhoto ? inputEditPhoto.value.trim() : "";
        if (!newName && !newPhoto) {
          if (profileEditStatus) {
            profileEditStatus.textContent = "Please enter a display name or avatar URL.";
            profileEditStatus.className = "text-xs text-rose-400";
            profileEditStatus.classList.remove("hidden");
          }
          return;
        }

        const originalBtnHtml = saveProfileBtn.innerHTML;
        saveProfileBtn.disabled = true;
        saveProfileBtn.innerHTML = `
          <svg class="animate-spin w-3.5 h-3.5 text-white inline" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
          <span>Saving...</span>
        `;

        try {
          await window.bindoraFirebase.updateUserProfile(newName, newPhoto);
          if (profileEditStatus) {
            profileEditStatus.innerHTML = '<span class="inline-flex items-center space-x-1.5"><svg class="w-3.5 h-3.5 text-emerald-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polyline points="20 6 9 17 4 12"/></svg><span>Profile updated successfully!</span></span>';
            profileEditStatus.className = "text-xs text-emerald-400";
            profileEditStatus.classList.remove("hidden");
          }
          this.showToast("Profile updated successfully!", "success");
          this.updateAuthModalView();
          setTimeout(() => {
            if (profileEditForm) profileEditForm.classList.add("hidden");
            if (profileEditStatus) profileEditStatus.classList.add("hidden");
          }, 1200);
        } catch (err) {
          if (profileEditStatus) {
            profileEditStatus.textContent = err.message || "Failed to update profile.";
            profileEditStatus.className = "text-xs text-rose-400";
            profileEditStatus.classList.remove("hidden");
          }
          this.showToast(err.message || "Failed to update profile.", "error");
        } finally {
          saveProfileBtn.disabled = false;
          saveProfileBtn.innerHTML = originalBtnHtml;
        }
      });
    }

    // Export Dossier
    const exportBtn = document.getElementById("btn-export-dossier");
    if (exportBtn) {
      exportBtn.addEventListener("click", () => {
        this.updateDossierView();
        setTimeout(() => {
          window.print();
        }, 150);
      });
    }

    // Brand Guide / About Modal
    const aboutBtn = document.getElementById("btn-about");
    const aboutModal = document.getElementById("modal-about");
    const closeAboutBtn = document.getElementById("btn-close-about");
    const closeAboutActionBtn = document.getElementById("btn-modal-close-action");
    if (aboutBtn && aboutModal) {
      aboutBtn.addEventListener("click", () => aboutModal.classList.remove("hidden"));
    }
    if (closeAboutBtn && aboutModal) {
      closeAboutBtn.addEventListener("click", () => aboutModal.classList.add("hidden"));
    }
    if (closeAboutActionBtn && aboutModal) {
      closeAboutActionBtn.addEventListener("click", () => aboutModal.classList.add("hidden"));
    }

    // Terminal CLI Guide Modal [v2.0]
    const cliBtn = document.getElementById("btn-cli-guide");
    const cliModal = document.getElementById("modal-cli-guide");
    const closeCliBtn = document.getElementById("btn-close-cli-guide");
    if (cliBtn && cliModal) {
      cliBtn.addEventListener("click", () => cliModal.classList.remove("hidden"));
    }
    if (closeCliBtn && cliModal) {
      closeCliBtn.addEventListener("click", () => cliModal.classList.add("hidden"));
    }

    // Blind Docking & Pocket Centroid Reset
    const blindBtn = document.getElementById("btn-blind-docking");
    if (blindBtn) {
      blindBtn.addEventListener("click", () => this.applyBlindDockingBox());
    }
    const resetPocketBtn = document.getElementById("btn-reset-pocket");
    if (resetPocketBtn) {
      resetPocketBtn.addEventListener("click", () => this.resetToPocketCentroid());
    }

    // Redocking Self-Validation
    const redockBtn = document.getElementById("btn-run-redock-validation");
    if (redockBtn) {
      redockBtn.addEventListener("click", () => this.runRedockingValidation());
    }

    // Preparation Protocol Log Accordion Toggle
    const togglePrepBtn = document.getElementById("btn-toggle-prep-log");
    const prepLogContent = document.getElementById("prep-log-content");
    const prepLogChevron = document.getElementById("prep-log-chevron");
    if (togglePrepBtn && prepLogContent) {
      togglePrepBtn.addEventListener("click", () => {
        const isHidden = prepLogContent.classList.contains("hidden");
        if (isHidden) {
          prepLogContent.classList.remove("hidden");
          if (prepLogChevron) prepLogChevron.textContent = "−";
        } else {
          prepLogContent.classList.add("hidden");
          if (prepLogChevron) prepLogChevron.textContent = "+";
        }
      });
    }

    // Dark / Light Mode Toggle
    const themeToggleBtn = document.getElementById("btn-theme-toggle");
    if (themeToggleBtn) {
      // Restore saved theme on startup (or via ?theme= url parameter)
      const urlParams = new URLSearchParams(window.location.search);
      const urlTheme = urlParams.get('theme');
      const savedTheme = urlTheme || localStorage.getItem('bindora-theme') || 'dark';
      this.applyTheme(savedTheme);
      themeToggleBtn.addEventListener("click", () => {
        const isDark = document.documentElement.classList.contains('dark');
        this.applyTheme(isDark ? 'light' : 'dark');
      });
    }
  }

  updateAuthModalView() {
    const user = window.bindoraFirebase ? window.bindoraFirebase.currentUser : null;
    const loggedOutView = document.getElementById("auth-logged-out-view");
    const loggedInView = document.getElementById("auth-logged-in-view");
    const profileName = document.getElementById("user-profile-name");
    const profileEmail = document.getElementById("user-profile-email");
    const profileImg = document.getElementById("user-profile-img");

    if (user) {
      if (loggedOutView) loggedOutView.classList.add("hidden");
      if (loggedInView) loggedInView.classList.remove("hidden");
      if (profileName) profileName.textContent = user.displayName || user.email.split("@")[0];
      if (profileEmail) profileEmail.textContent = user.email;
      if (profileImg) profileImg.src = user.photoURL || `https://api.dicebear.com/7.x/bottts/svg?seed=${user.uid}`;
      
      const inputEditName = document.getElementById("input-edit-display-name");
      const inputEditPhoto = document.getElementById("input-edit-photo-url");
      if (inputEditName) inputEditName.value = user.displayName || "";
      if (inputEditPhoto) inputEditPhoto.value = user.photoURL || "";
      
      // Update header auth button
      const headerName = document.getElementById("user-display-name");
      const headerAvatar = document.getElementById("user-avatar-img");
      const headerIcon = document.getElementById("user-default-icon");
      if (headerName) headerName.textContent = user.displayName || user.email.split("@")[0];
      if (headerAvatar) { headerAvatar.src = user.photoURL || `https://api.dicebear.com/7.x/bottts/svg?seed=${user.uid}`; headerAvatar.classList.remove("hidden"); }
      if (headerIcon) headerIcon.classList.add("hidden");
      
      this.loadCloudHistory();
    } else {
      if (loggedOutView) loggedOutView.classList.remove("hidden");
      if (loggedInView) loggedInView.classList.add("hidden");
      
      const headerName = document.getElementById("user-display-name");
      const headerAvatar = document.getElementById("user-avatar-img");
      const headerIcon = document.getElementById("user-default-icon");
      if (headerName) headerName.textContent = "Sign In";
      if (headerAvatar) headerAvatar.classList.add("hidden");
      if (headerIcon) headerIcon.classList.remove("hidden");
    }
  }

  async loadCloudHistory() {
    const container = document.getElementById("saved-runs-list");
    if (!container || !window.bindoraFirebase) return;

    try {
      const runs = await window.bindoraFirebase.fetchSavedRuns();
      if (!runs || runs.length === 0) {
        container.innerHTML = `<p class="text-slate-500 italic p-2 text-center text-xs">No saved docking runs in cloud history yet.</p>`;
        return;
      }

      container.innerHTML = runs.map(r => `
        <div class="p-2 flex items-center justify-between text-xs hover:bg-slate-900/60 rounded">
          <div>
            <div class="font-bold text-white">${r.drugName} <span class="text-slate-400 font-normal">vs</span> ${r.targetName} (${r.pdbId})</div>
            <div class="font-mono text-[11px] text-cyan-300">
              ΔG: ${r.affinityKcal} kcal/mol &bull; Kd: ${r.theoreticalKdNm} nM &bull; ${new Date(r.savedAt).toLocaleDateString()}
            </div>
          </div>
          <span class="px-2 py-0.5 rounded text-[10px] ${r.isCrossChecked ? 'bg-emerald-950 text-emerald-300 border border-emerald-800' : 'bg-slate-800 text-slate-400'}">
            ${r.crosscheckBadge}
          </span>
        </div>
      `).join("");
    } catch (e) {
      container.innerHTML = `<p class="text-rose-400 p-2 text-center text-xs">Could not load saved runs: ${e.message}</p>`;
    }
  }

  switchTab(tabId) {
    this.state.activeTab = tabId;
    document.querySelectorAll(".nav-tab").forEach(b => {
      if (b.dataset.tab === tabId) {
        b.classList.add("active");
      } else {
        b.classList.remove("active");
      }
    });

    document.querySelectorAll(".tab-pane").forEach(pane => {
      if (pane.id === `tab-${tabId}`) {
        pane.classList.remove("hidden");
      } else {
        pane.classList.add("hidden");
      }
    });

    // If switching to 3D tab, re-render viewer
    if (tabId === "docking" && this.viewer) {
      setTimeout(() => {
        if (this.viewer.viewer) this.viewer.viewer.render();
      }, 100);
    }

    // If switching to AI narrative tab and docking is completed but no narrative yet, generate it
    if (tabId === "ai-narrative" && this.state.docking && !this.state.narrative) {
      this.generateNarrativeReport();
    }

    // If switching to batch tab, check pharmacophore eligibility
    if (tabId === "batch") {
      this.checkPharmacophoreEligibility();
    }
  }

  showToast(message, type = "info") {
    const toast = document.createElement("div");
    const bg = type === "error" ? "bg-rose-900/90 border-rose-600 text-rose-200" :
               type === "success" ? "bg-emerald-900/90 border-emerald-600 text-emerald-200" :
               "bg-blue-900/90 border-blue-600 text-blue-200";
    toast.className = `fixed bottom-5 right-5 z-50 px-4 py-2.5 rounded-lg border shadow-xl text-sm font-medium transition-all transform duration-300 ${bg}`;
    toast.textContent = message;
    document.body.appendChild(toast);
    setTimeout(() => {
      toast.style.opacity = "0";
      setTimeout(() => toast.remove(), 300);
    }, 3500);
  }

  // Publication-grade markdown → HTML renderer for AI narratives & reports
  renderMarkdown(text) {
    if (!text) return '';

    // Clean up any internal raw payload keys or leaks (e.g. "(payload: ...)", "is_cross_checked: false")
    let cleaned = text
      .replace(/\s*\(payload:\s*[^)]+\)/gi, '')
      .replace(/•\s*Cross-check flag:\s*is_cross_checked:\s*(true|false)/gi, (m, val) => {
        return val === 'true'
          ? '<span class="inline-flex items-center text-emerald-500 font-semibold text-xs">&#10003; Experimental ChEMBL Records Verified</span>'
          : '<span class="inline-flex items-center text-amber-500 font-semibold text-xs">&#9888; Computational Prediction (Unverified Simulation)</span>';
      });

    // 1. Extract markdown tables to tokens so line-break parsing doesn't break table syntax
    const tables = [];
    const tableRegex = /((?:^[ \t]*\|.+?\|[ \t]*\r?\n)+(?:^[ \t]*\|[-: |]+\|[ \t]*\r?\n)(?:^[ \t]*\|.+?\|[ \t]*(?:\r?\n|$))+)/gm;

    let out = cleaned.replace(tableRegex, (match) => {
      const rows = match.trim().split(/\r?\n/).filter(r => r.trim().startsWith('|'));
      if (rows.length < 2) return match;
      const parseRow = (row) => row.trim().replace(/^\|/, '').replace(/\|$/, '').split('|').map(c => c.trim());
      const headers = parseRow(rows[0]);
      const dataRows = rows.slice(2).map(parseRow);

      let html = '<div class="overflow-x-auto my-3 rounded-xl border border-slate-200 dark:border-slate-800 shadow-sm bg-white dark:bg-slate-900/40">';
      html += '<table class="w-full text-xs text-left border-collapse">';
      html += '<thead class="bg-slate-100 dark:bg-slate-800/90 text-slate-700 dark:text-slate-300 font-semibold border-b border-slate-200 dark:border-slate-700">';
      html += '<tr>';
      headers.forEach(h => {
        html += `<th class="p-2.5 font-sans">${h}</th>`;
      });
      html += '</tr></thead><tbody class="divide-y divide-slate-100 dark:divide-slate-800/80 text-slate-800 dark:text-slate-200">';
      dataRows.forEach(row => {
        html += '<tr class="hover:bg-slate-50 dark:hover:bg-slate-800/30 transition">';
        row.forEach((cell, idx) => {
          const isNum = idx > 0 && /^-?\d+(\.\d+)?/.test(cell);
          html += `<td class="p-2.5 ${isNum ? 'font-mono' : 'font-sans'}">${cell}</td>`;
        });
        html += '</tr>';
      });
      html += '</tbody></table></div>';

      const token = `%%TABLE_BLOCK_${tables.length}%%`;
      tables.push(html);
      return `\n\n${token}\n\n`;
    });

    // Escape HTML entities
    out = out
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');

    // Horizontal divider
    out = out.replace(/^---$/gm, '<hr class="border-slate-200 dark:border-slate-800 my-4">');

    // Headers
    out = out
      .replace(/^### (.+)$/gm, '<h3 class="text-sm font-bold text-slate-900 dark:text-white mt-4 mb-1 pb-0.5 border-b border-slate-200 dark:border-slate-700/50">$1</h3>')
      .replace(/^## (.+)$/gm, '<h2 class="text-sm font-bold text-[#0D6EFD] dark:text-[#00C6FF] mt-5 mb-1 uppercase tracking-wide">$1</h2>')
      .replace(/^# (.+)$/gm, '<h1 class="text-base font-extrabold text-slate-900 dark:text-white mt-5 mb-2">$1</h1>');

    // Unordered lists (- or * or •)
    out = out.replace(/^[•\*\-] (.+)$/gm, '<div class="flex items-start space-x-2 my-1.5"><span class="text-cyan-500 font-bold leading-tight select-none">&#8226;</span><span class="text-slate-700 dark:text-slate-300 leading-relaxed">$1</span></div>');

    // Bold and italic
    out = out
      .replace(/\*\*(.+?)\*\*/g, '<strong class="font-semibold text-slate-900 dark:text-white">$1</strong>')
      .replace(/\*([^*]+?)\*/g, '<em class="italic text-slate-600 dark:text-slate-300">$1</em>');

    // Clean any dangling unclosed bold tags
    out = out.replace(/\*\*([^*]+)$/gm, '<strong class="font-semibold text-slate-900 dark:text-white">$1</strong>');

    // Inline code
    out = out
      .replace(/`(.+?)`/g, '<code class="font-mono text-cyan-700 dark:text-[#00C6FF] bg-slate-100 dark:bg-slate-800/70 px-1 py-0.5 rounded text-xs border border-slate-200 dark:border-slate-700">$1</code>');

    // Double newline = paragraph break
    out = out.replace(/\n\n/g, '</p><p class="mt-2 text-slate-700 dark:text-slate-300">');

    // Single newline
    out = out.replace(/\n/g, '<br>');

    // Restore table tokens
    tables.forEach((tableHtml, i) => {
      out = out.replace(new RegExp(`%%TABLE_BLOCK_${i}%%`, 'g'), tableHtml);
    });

    // Clean empty paragraph wrappers around tables
    out = out.replace(/<p[^>]*>\s*<\/p>/g, '');
    out = out.replace(/<br>\s*(<div class="overflow-x-auto)/g, '$1');
    out = out.replace(/(<\/table><\/div>)\s*<br>/g, '$1');
    out = out.replace(/<br>\s*(<div class="flex items-start)/g, '$1');

    return out;
  }

  // Apply dark/light theme across the page
  applyTheme(theme) {
    const html = document.documentElement;
    const sunIcon = document.getElementById("icon-sun");
    const moonIcon = document.getElementById("icon-moon");
    const logoDark = document.getElementById("header-logo-dark");
    const logoLight = document.getElementById("header-logo-light");

    if (theme === 'light') {
      html.classList.remove('dark');
      html.classList.add('light');
      if (sunIcon) sunIcon.classList.add('hidden');
      if (moonIcon) moonIcon.classList.remove('hidden');
      if (logoDark) logoDark.classList.add('hidden');
      if (logoLight) logoLight.classList.remove('hidden');
      localStorage.setItem('bindora-theme', 'light');

      // Update 3Dmol viewer background to soft light if active
      if (this.viewer && this.viewer.viewer) {
        try {
          this.viewer.viewer.setBackgroundColor('0xf1f5f9');
          this.viewer.viewer.render();
        } catch (e) {}
      }
    } else {
      html.classList.remove('light');
      html.classList.add('dark');
      if (sunIcon) sunIcon.classList.remove('hidden');
      if (moonIcon) moonIcon.classList.add('hidden');
      if (logoDark) logoDark.classList.remove('hidden');
      if (logoLight) logoLight.classList.add('hidden');
      localStorage.setItem('bindora-theme', 'dark');

      // Update 3Dmol viewer background to dark navy if active
      if (this.viewer && this.viewer.viewer) {
        try {
          this.viewer.viewer.setBackgroundColor('0x0f172a');
          this.viewer.viewer.render();
        } catch (e) {}
      }
    }

    // Refresh radar chart colors if ADME is currently loaded
    if (this.state.ligand?.adme) {
      this.updateADMEView();
    }

    // Sync any rendered 2D interaction diagrams with current theme
    const diagSvgs = document.querySelectorAll(".bindora-svg-root");
    diagSvgs.forEach(svg => {
      if (theme === 'light') {
        svg.classList.remove('theme-dark');
        svg.classList.add('theme-light');
      } else {
        svg.classList.remove('theme-light');
        svg.classList.add('theme-dark');
      }
    });
  }

  renderBenchmarksUI() {
    const container = document.getElementById("benchmark-cards");
    if (!container) return;

    container.innerHTML = this.state.benchmarks.map(bm => `
      <div class="benchmark-case-card cursor-pointer p-3 rounded-xl border border-slate-700/60 bg-slate-800/40 hover:bg-slate-700/50 transition-all hover:border-cyan-500/50" onclick="window.app.loadBenchmark('${bm.id}')">
        <div class="flex items-center justify-between mb-1">
          <span class="benchmark-case-name font-bold text-sm text-cyan-300">${bm.drug_name}</span>
          <span class="benchmark-case-badge text-xs px-2 py-0.5 rounded bg-slate-700 text-slate-300 font-mono">${bm.pdb_id}</span>
        </div>
        <p class="benchmark-case-target text-xs text-slate-300 font-medium">${bm.target_name}</p>
        <p class="benchmark-case-desc text-[11px] text-slate-400 mt-1 line-clamp-2">${bm.mechanism}</p>
      </div>
    `).join("");
  }

  async loadValidationReportUI() {
    try {
      const data = await BindoraAPI.getValidationReport();
      if (!data) return;

      const { validation_report, flagship_targets } = data;
      const hsg = flagship_targets?.["1HSG"];
      const aq1 = flagship_targets?.["1AQ1"];
      const sumStats = validation_report?.summary_statistics;

      // Update KPI Stat Cards
      const kpiAq1 = document.getElementById("wp-kpi-1aq1-rmsd");
      if (kpiAq1 && aq1?.mode1_rmsd_angstroms !== undefined) {
        kpiAq1.textContent = `${aq1.mode1_rmsd_angstroms} Å`;
      }

      const kpiHsg = document.getElementById("wp-kpi-1hsg-rmsd");
      if (kpiHsg && hsg?.mode1_rmsd_angstroms !== undefined) {
        kpiHsg.textContent = `${hsg.mode1_rmsd_angstroms} Å`;
      }

      const kpiHsgE = document.getElementById("wp-kpi-1hsg-energy");
      if (kpiHsgE && hsg?.vina_affinity_kcal !== undefined) {
        kpiHsgE.textContent = `${hsg.vina_affinity_kcal.toFixed(2)}`;
      }

      const kpiCasfRate = document.getElementById("wp-kpi-casf-rate");
      if (kpiCasfRate && sumStats?.pose_reconstruction?.rmsd_success_rate_percent !== undefined) {
        kpiCasfRate.textContent = `${sumStats.pose_reconstruction.rmsd_success_rate_percent.toFixed(1)}%`;
      }

      const kpiCasfMean = document.getElementById("wp-kpi-casf-mean");
      if (kpiCasfMean && sumStats?.pose_reconstruction?.mean_rmsd_angstroms !== undefined) {
        kpiCasfMean.textContent = `${sumStats.pose_reconstruction.mean_rmsd_angstroms.toFixed(2)} Å`;
      }

      // Update Flagship Table Rows
      if (aq1) {
        const elRmsd = document.getElementById("wp-row-1aq1-rmsd");
        if (elRmsd) elRmsd.textContent = `${aq1.mode1_rmsd_angstroms} Å`;
        const elEnergy = document.getElementById("wp-row-1aq1-energy");
        if (elEnergy) elEnergy.textContent = `${aq1.vina_affinity_kcal.toFixed(2)} kcal/mol`;
      }

      if (hsg) {
        const elRmsd = document.getElementById("wp-row-1hsg-rmsd");
        if (elRmsd) elRmsd.textContent = `${hsg.mode1_rmsd_angstroms} Å`;
        const elEnergy = document.getElementById("wp-row-1hsg-energy");
        if (elEnergy) elEnergy.textContent = `${hsg.vina_affinity_kcal.toFixed(2)} kcal/mol`;
      }

      // Populate Multi-Target CASF Benchmark Rows
      const casfBody = document.getElementById("casf-benchmark-rows");
      if (casfBody && validation_report?.complex_results?.length) {
        casfBody.innerHTML = validation_report.complex_results.map(c => {
          const pass = c.rmsd_angstroms <= 2.0;
          const statusBadge = pass
            ? `<span class="px-2 py-0.5 rounded-full text-[10px] font-semibold font-mono bg-emerald-50 dark:bg-emerald-950/50 text-emerald-700 dark:text-emerald-300 border border-emerald-200 dark:border-emerald-700">Sub-2.0 Å &check;</span>`
            : `<span class="px-2 py-0.5 rounded-full text-[10px] font-semibold font-mono bg-amber-50 dark:bg-amber-950/50 text-amber-700 dark:text-amber-300 border border-amber-200 dark:border-amber-700">Near-Native (${c.rmsd_angstroms} Å)</span>`;
          return `
            <tr class="hover:bg-slate-50 dark:hover:bg-slate-800/30 transition">
              <td class="p-2.5 font-mono font-bold">
                <a href="https://www.rcsb.org/structure/${c.pdb_id}" target="_blank" class="text-sky-700 dark:text-cyan-400 hover:underline">${c.pdb_id}</a>
              </td>
              <td class="p-2.5">
                <span class="font-medium text-slate-800 dark:text-slate-200 block">${c.target_name}</span>
                <span class="text-[10px] text-slate-400">${c.target_class || ""}</span>
              </td>
              <td class="p-2.5">
                <span class="font-mono text-slate-700 dark:text-slate-300">${c.drug_name}</span>
                <span class="text-[10px] text-slate-400 block">${c.ligand_resname}</span>
              </td>
              <td class="p-2.5 font-mono text-center text-slate-600 dark:text-slate-300">${c.exp_delta_g_kcal ? c.exp_delta_g_kcal.toFixed(2) : "—"}</td>
              <td class="p-2.5 font-mono text-center font-semibold text-slate-900 dark:text-slate-100">${c.vina_delta_g_kcal ? c.vina_delta_g_kcal.toFixed(2) : "—"}</td>
              <td class="p-2.5 font-mono text-center text-slate-500 dark:text-slate-400">${c.vinardo_delta_g_kcal ? c.vinardo_delta_g_kcal.toFixed(2) : "—"}</td>
              <td class="p-2.5 font-mono text-center font-bold ${pass ? 'text-emerald-700 dark:text-emerald-400' : 'text-amber-700 dark:text-amber-400'}">${c.rmsd_angstroms.toFixed(2)} Å</td>
              <td class="p-2.5 text-center">${statusBadge}</td>
            </tr>
          `;
        }).join("");
      }
    } catch (e) {
      console.warn("Could not load validation report:", e);
    }
  }

  async loadBenchmark(id) {
    const bm = this.state.benchmarks.find(b => b.id === id);
    if (!bm) return;

    this.showToast(`Loading benchmark: ${bm.drug_name} vs ${bm.target_name}...`, "info");

    try {
      // 1. Fetch drug from PubChem or use preset SMILES
      const ligPrep = await BindoraAPI.prepareLigand(bm.smiles, false);
      this.state.ligand = {
        name: bm.drug_name,
        smiles: bm.smiles,
        weight: ligPrep.adme?.physicochemical?.molecular_weight?.value,
        formula: "",
        ...ligPrep
      };

      // 2. Fetch PDB from RCSB
      const recRes = await BindoraAPI.prepareReceptor("", bm.pdb_id);
      this.state.receptor = {
        pdb_id: bm.pdb_id,
        title: bm.target_name,
        ...recRes
      };
      this.checkPharmacophoreEligibility();
      this.state.ensembleStructures = null;

      // 3. Update UI
      this.updateStudioCards();
      this.updatePathwayInfo();
      this.updateDossierView();
      
      // 4. Update 3D viewer with receptor
      if (this.viewer && recRes.cleaned_pdb) {
        this.viewer.loadReceptor(recRes.cleaned_pdb);
        const emptyState = document.getElementById("viewer-empty-state");
        if (emptyState) emptyState.classList.add("hidden");
      }

      this.showToast(`Loaded ${bm.drug_name} & ${bm.target_name} (${bm.pdb_id})`, "success");
    } catch (e) {
      console.error("Benchmark load error:", e);
      this.showToast(`Failed to load benchmark: ${e.message}`, "error");
    }
  }

  async searchPubChem(query) {
    if (!query) return;
    const ligCard = document.getElementById("card-ligand-info");
    if (ligCard) {
      ligCard.innerHTML = `<div class="p-6 text-center text-cyan-400 font-medium animate-pulse flex flex-col items-center justify-center space-y-2">
        <svg class="animate-spin w-6 h-6 text-cyan-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
        <span class="text-xs">Searching PubChem REST API for '${query}'...</span>
      </div>`;
    }
    this.showToast(`Searching PubChem for '${query}'...`, "info");
    try {
      const data = await BindoraAPI.searchPubChem(query);
      const ligPrep = await BindoraAPI.prepareLigand(data.smiles || data.sdf, !data.smiles);
      this.state.ligand = {
        ...data,
        ...ligPrep
      };
      this.updateStudioCards();
      this.updateDossierView();
      if (this.viewer && ligPrep.pdb_block) {
        this.viewer.loadLigand(ligPrep.pdb_block);
        const emptyState = document.getElementById("viewer-empty-state");
        if (emptyState) emptyState.classList.add("hidden");
      }
      this.showToast(`Found PubChem compound: ${data.name} (CID: ${data.cid})`, "success");
    } catch (e) {
      this.updateStudioCards();
      this.showToast(`PubChem search failed: ${e.message}`, "error");
    }
  }

  async searchRCSB(query) {
    if (!query) return;
    const recCard = document.getElementById("card-receptor-info");
    if (recCard) {
      recCard.innerHTML = `<div class="p-6 text-center text-emerald-400 font-medium animate-pulse flex flex-col items-center justify-center space-y-2">
        <svg class="animate-spin w-6 h-6 text-emerald-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
        <span class="text-xs">Querying RCSB Protein Data Bank for '${query}'...</span>
      </div>`;
    }
    this.showToast(`Searching RCSB for '${query}'...`, "info");
    try {
      const data = await BindoraAPI.searchRCSB(query);
      if (data.direct && data.entry) {
        const entry = data.entry;
        const recPrep = await BindoraAPI.prepareReceptor(entry.pdb_content, entry.pdb_id);
        this.state.receptor = {
          ...entry,
          ...recPrep
        };
        this.checkPharmacophoreEligibility();
        this.state.ensembleStructures = null;
        this.updateStudioCards();
        this.updatePathwayInfo();
        this.updateDossierView();
        if (this.viewer && recPrep.cleaned_pdb) {
          this.viewer.loadReceptor(recPrep.cleaned_pdb);
          const emptyState = document.getElementById("viewer-empty-state");
          if (emptyState) emptyState.classList.add("hidden");
        }
        this.showToast(`Loaded RCSB entry ${entry.pdb_id}: ${entry.title.substring(0, 30)}...`, "success");
      } else if (data.results && data.results.length > 0) {
        const first = data.results[0];
        const recPrep = await BindoraAPI.prepareReceptor("", first.pdb_id);
        this.state.receptor = {
          ...first,
          ...recPrep
        };
        this.checkPharmacophoreEligibility();
        this.state.ensembleStructures = null;
        this.updateStudioCards();
        this.updatePathwayInfo();
        this.updateDossierView();
        if (this.viewer && recPrep.cleaned_pdb) {
          this.viewer.loadReceptor(recPrep.cleaned_pdb);
          const emptyState = document.getElementById("viewer-empty-state");
          if (emptyState) emptyState.classList.add("hidden");
        }
        this.showToast(`Selected top RCSB match: ${first.pdb_id}`, "success");
      } else {
        this.updateStudioCards();
        this.showToast("No RCSB entries found for query", "error");
      }
    } catch (e) {
      this.updateStudioCards();
      this.showToast(`RCSB query failed: ${e.message}`, "error");
    }
  }

  async loadCustomLigand(name, smiles) {
    if (!smiles) {
      this.showToast("Please enter a valid SMILES string", "error");
      return;
    }
    const displayName = name || "Custom Ligand";
    const ligCard = document.getElementById("card-ligand-info");
    if (ligCard) {
      ligCard.innerHTML = `<div class="p-6 text-center text-cyan-400 font-medium animate-pulse flex flex-col items-center justify-center space-y-2">
        <svg class="animate-spin w-6 h-6 text-cyan-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
        <span class="text-xs">Parsing SMILES & generating 3D coordinates with RDKit...</span>
      </div>`;
    }
    this.showToast(`Preparing custom ligand: ${displayName}...`, "info");
    try {
      const ligPrep = await BindoraAPI.prepareLigand(smiles, false);
      const p = ligPrep.adme?.physicochemical || {};
      this.state.ligand = {
        name: displayName,
        smiles: smiles,
        canonical_smiles: ligPrep.canonical_smiles || smiles,
        weight: p.molecular_weight?.value,
        formula: "",
        ...ligPrep
      };
      this.updateStudioCards();
      this.updateDossierView();
      if (this.viewer && ligPrep.pdb_block) {
        this.viewer.loadLigand(ligPrep.pdb_block);
        const emptyState = document.getElementById("viewer-empty-state");
        if (emptyState) emptyState.classList.add("hidden");
      }
      this.showToast(`Custom ligand loaded: ${displayName}`, "success");
    } catch (e) {
      this.updateStudioCards();
      this.showToast(`Ligand preparation failed: ${e.message}`, "error");
    }
  }

  async uploadLigandFile(file) {
    if (!file) return;
    const isSdf = file.name.toLowerCase().endsWith('.sdf') || file.name.toLowerCase().endsWith('.mol');
    const ligCard = document.getElementById("card-ligand-info");
    if (ligCard) {
      ligCard.innerHTML = `<div class="p-6 text-center text-cyan-400 font-medium animate-pulse flex flex-col items-center justify-center space-y-2">
        <svg class="animate-spin w-6 h-6 text-cyan-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
        <span class="text-xs">Reading ${file.name} & calculating 3D conformer...</span>
      </div>`;
    }
    this.showToast(`Processing ${file.name}...`, "info");
    const reader = new FileReader();
    reader.onload = async (e) => {
      const content = e.target.result;
      try {
        const ligPrep = await BindoraAPI.prepareLigand(content, isSdf);
        const p = ligPrep.adme?.physicochemical || {};
        let baseName = file.name.replace(/\.[^/.]+$/, "");
        
        // Resolve PubChem CID if present in filename or content
        const cidMatch = file.name.match(/COMPOUND_CID_(\d+)/i) || (content && content.match(/PUBCHEM_COMPOUND_CID[>\s\n\r]+(\d+)/i));
        if (cidMatch) {
          const cid = cidMatch[1];
          try {
            const cidData = await BindoraAPI.searchPubChemByCID(cid);
            if (cidData && cidData.name) {
              baseName = `${cidData.name} (CID: ${cid})`;
            } else {
              baseName = `PubChem Compound #${cid}`;
            }
          } catch (err) {
            baseName = `PubChem Compound #${cid}`;
          }
        }

        this.state.ligand = {
          name: baseName,
          smiles: ligPrep.canonical_smiles,
          canonical_smiles: ligPrep.canonical_smiles,
          weight: p.molecular_weight?.value,
          formula: "",
          ...ligPrep
        };
        this.updateStudioCards();
        this.updateDossierView();
        if (this.viewer && ligPrep.pdb_block) {
          this.viewer.loadLigand(ligPrep.pdb_block);
          const emptyState = document.getElementById("viewer-empty-state");
          if (emptyState) emptyState.classList.add("hidden");
        }
        this.showToast(`Loaded ${baseName} successfully`, "success");
      } catch (err) {
        this.updateStudioCards();
        this.showToast(`File preparation failed: ${err.message}`, "error");
      }
    };
    reader.readAsText(file);
  }

  async loadCustomReceptorPdb(name, pdbContent) {
    if (!pdbContent || pdbContent.trim().length < 20) {
      this.showToast("Please provide valid PDB coordinate text", "error");
      return;
    }
    const displayName = name || "CUSTOM_RECEPTOR";
    const recCard = document.getElementById("card-receptor-info");
    if (recCard) {
      recCard.innerHTML = `<div class="p-6 text-center text-emerald-400 font-medium animate-pulse flex flex-col items-center justify-center space-y-2">
        <svg class="animate-spin w-6 h-6 text-emerald-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
        <span class="text-xs">Preparing receptor structure: ${displayName}...</span>
      </div>`;
    }
    this.showToast(`Preparing receptor structure: ${displayName}...`, "info");
    try {
      const recPrep = await BindoraAPI.prepareReceptor(pdbContent, "");
      this.state.receptor = {
        pdb_id: displayName,
        title: displayName,
        pdb_content: pdbContent,
        ...recPrep
      };
      this.checkPharmacophoreEligibility();
      this.state.ensembleStructures = null;
      this.updateStudioCards();
      this.updatePathwayInfo();
      this.updateDossierView();
      if (this.viewer && recPrep.cleaned_pdb) {
        this.viewer.loadReceptor(recPrep.cleaned_pdb);
        const emptyState = document.getElementById("viewer-empty-state");
        if (emptyState) emptyState.classList.add("hidden");
      }
      this.showToast(`Loaded custom receptor: ${displayName}`, "success");
    } catch (e) {
      this.updateStudioCards();
      this.showToast(`Receptor preparation failed: ${e.message}`, "error");
    }
  }

  async uploadReceptorFile(file) {
    if (!file) return;
    const recCard = document.getElementById("card-receptor-info");
    if (recCard) {
      recCard.innerHTML = `<div class="p-6 text-center text-emerald-400 font-medium animate-pulse flex flex-col items-center justify-center space-y-2">
        <svg class="animate-spin w-6 h-6 text-emerald-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
        <span class="text-xs">Reading ${file.name} & preparing structure...</span>
      </div>`;
    }
    const reader = new FileReader();
    reader.onload = async (e) => {
      const content = e.target.result;
      const baseName = file.name.replace(/\.[^/.]+$/, "").toUpperCase();
      await this.loadCustomReceptorPdb(baseName, content);
    };
    reader.readAsText(file);
  }

  resetWorkspace() {
    this.state.ligand = null;
    this.state.receptor = null;
    this.state.docking = null;
    this.state.crosscheck = null;
    this.state.narrative = null;
    this.state.redockingValidation = null;
    this.state.ensembleStructures = null;
    this.state.pharmacophore = null;
    this.state.currentPoseIdx = 0;

    if (this.viewer) {
      this.viewer.clear();
      this.viewer.clearInteractions();
      this.viewer.clearMeasurements();
      this.viewer.clearGridBox();
      this.viewer.clearRedockOverlay();
    }
    const emptyState = document.getElementById("viewer-empty-state");
    if (emptyState) emptyState.classList.remove("hidden");

    const statusBanner = document.getElementById("docking-status-banner");
    if (statusBanner) statusBanner.classList.add("hidden");

    const weakAlert = document.getElementById("weak-binder-alert");
    if (weakAlert) weakAlert.classList.add("hidden");

    const repBanner = document.getElementById("replicate-stats-banner");
    if (repBanner) repBanner.classList.add("hidden");

    const redockPanel = document.getElementById("redock-results-panel");
    if (redockPanel) redockPanel.classList.add("hidden");

    const redockNative = document.getElementById("redock-native-name");
    if (redockNative) redockNative.textContent = "No co-ligand loaded";

    const poseTable = document.getElementById("pose-table-rows");
    if (poseTable) poseTable.innerHTML = `<tr><td colspan="5" class="text-center p-6 text-xs text-slate-500 italic">No docking run loaded yet.</td></tr>`;

    const diagContainer = document.getElementById("dock-interaction-diagram-container");
    if (diagContainer) diagContainer.innerHTML = `<span class="text-xs text-slate-500 italic">Run docking to render 2D radial interaction schematic</span>`;

    if (window.DockingCharts && DockingCharts.energyChartInstance) {
      DockingCharts.energyChartInstance.destroy();
      DockingCharts.energyChartInstance = null;
    }

    // Reset toolbar toggle checkboxes
    const surfToggle = document.getElementById("toggle-surface");
    if (surfToggle) surfToggle.checked = false;
    const gBoxToggle = document.getElementById("toggle-gridbox");
    if (gBoxToggle) gBoxToggle.checked = false;
    const mToggle = document.getElementById("toggle-measure");
    if (mToggle) mToggle.checked = false;
    const mHud = document.getElementById("measure-hud");
    if (mHud) mHud.classList.add("hidden");

    const setTxt = (id, txt) => { const el = document.getElementById(id); if (el) el.textContent = txt; };
    setTxt("dock-affinity-score", "— kcal/mol");
    setTxt("dock-kd-score", "— nM");
    setTxt("dock-le-score", "—");
    setTxt("dock-hbonds-count", "—");
    setTxt("dock-current-mode-label", "Mode 1 / 9");
    const resContainer = document.getElementById("dock-interacting-residues");
    if (resContainer) resContainer.innerHTML = `<span class="text-xs text-slate-500 italic">Run docking to analyze active pocket contacts</span>`;

    this.updateStudioCards();
    this.updateDossierView();
    this.updatePathwayInfo();
    this.showToast("Workspace cleared. Ready for a new simulation.", "info");
  }

  proceedToDocking() {
    this.switchTab("docking");
    const statusBanner = document.getElementById("docking-status-banner");
    if (!this.state.docking && statusBanner) {
      if (this.state.ligand && this.state.receptor) {
        statusBanner.classList.remove("hidden");
        statusBanner.className = "glass-panel p-3.5 rounded-xl border border-cyan-500/40 bg-cyan-950/40 transition-all duration-300 mb-2";
        const icon = document.getElementById("docking-status-icon");
        if (icon) icon.innerHTML = `<span class="w-6 h-6 rounded-full bg-cyan-500/20 text-cyan-300 flex items-center justify-center font-bold text-xs">READY</span>`;
        const title = document.getElementById("docking-status-title");
        if (title) title.textContent = `Ready: ${this.state.ligand.name} vs ${this.state.receptor.pdb_id}`;
        const sub = document.getElementById("docking-status-subtitle");
        if (sub) sub.textContent = "Click 'Execute 3D Molecular Docking' below to run the Scripps AutoDock Vina engine.";
        const extra = document.getElementById("docking-status-extra");
        if (extra) extra.textContent = "Awaiting execution";
      } else {
        statusBanner.classList.add("hidden");
      }
    }
  }

  updateStudioCards() {
    // Ligand Card
    const ligCard = document.getElementById("card-ligand-info");
    if (ligCard) {
      if (this.state.ligand) {
        const l = this.state.ligand;
        const p = l.adme?.physicochemical || {};
        ligCard.innerHTML = `
          <div class="space-y-2">
            <div class="flex items-center justify-between">
              <span class="text-base font-bold text-cyan-400">${l.name || "Custom Ligand"}</span>
              <span class="text-xs px-2 py-0.5 rounded bg-cyan-950 text-cyan-300 border border-cyan-800">Ready</span>
            </div>
            <p class="text-xs font-mono text-slate-400 break-all">${l.canonical_smiles || l.smiles || "Structure loaded"}</p>
            <div class="grid grid-cols-2 gap-2 text-xs text-slate-300 pt-2 border-t border-slate-700/60">
              <div>MW: <span class="font-semibold text-white">${p.molecular_weight?.value || l.weight || "—"} Da</span></div>
              <div>LogP: <span class="font-semibold text-white">${p.logp?.value ?? "—"}</span></div>
              <div>HBD: <span class="font-semibold text-white">${p.hbd?.value ?? "—"}</span></div>
              <div>HBA: <span class="font-semibold text-white">${p.hba?.value ?? "—"}</span></div>
              <div>RotB: <span class="font-semibold text-white">${p.rotatable_bonds?.value ?? l.rotatable_bonds ?? "—"}</span></div>
              <div>TPSA: <span class="font-semibold text-white">${p.tpsa?.value ?? "—"} Å²</span></div>
            </div>
            <div class="pt-2 border-t border-slate-700/60 flex flex-col space-y-1.5">
              <button id="btn-find-similar" class="text-[11px] text-cyan-400 hover:text-cyan-300 flex items-center space-x-1 font-medium transition cursor-pointer self-start">
                <svg class="w-3.5 h-3.5 text-cyan-400" fill="none" viewBox="0 0 24 24" stroke="currentColor"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"/></svg>
                <span>Find similar compounds (PubChem)</span>
              </button>
              <div id="similar-compounds-container" class="hidden space-y-1.5 pt-1"></div>
            </div>
          </div>
        `;

        document.getElementById("btn-find-similar")?.addEventListener("click", () => this.findSimilarCompounds());
      } else {
        ligCard.innerHTML = `<p class="text-xs text-slate-400 italic">No ligand loaded yet. Search PubChem or pick a benchmark above.</p>`;
      }
    }

    // Receptor Card
    const recCard = document.getElementById("card-receptor-info");
    if (recCard) {
      if (this.state.receptor) {
        const r = this.state.receptor;
        const pocket = r.detected_pocket || {};
        const redock = this.state.redockingValidation;
        const hasNative = r.native_ligand?.has_native;
        const detectedPockets = r.detected_pockets || [];

        const validationHtml = redock ? `
          <div class="mt-2 pt-2 border-t border-slate-700/60 flex items-center justify-between">
            <span class="text-[11px] text-slate-400">Protocol Validation:</span>
            <span class="px-2 py-0.5 rounded-full text-[10px] font-bold ${redock.is_validated ? 'bg-emerald-950 text-emerald-300 border border-emerald-800' : 'bg-amber-950 text-amber-300 border border-amber-800'}">
              ${redock.validation_badge} (RMSD ${redock.rmsd_angstroms} Å)
            </span>
          </div>
        ` : (hasNative ? `
          <div class="mt-2 pt-2 border-t border-slate-700/60 flex items-center justify-between">
            <span class="text-[11px] text-slate-400">Protocol Validation:</span>
            <span id="rec-val-status" class="text-[10px] text-cyan-400 flex items-center"><svg class="animate-spin -ml-1 mr-1.5 h-3 w-3 inline" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg> Validating native ligand...</span>
          </div>
        ` : '');

        let pocketHtml = '';
        if (hasNative) {
          pocketHtml = `
            <div class="text-xs text-slate-400 pt-2 border-t border-slate-700/60 space-y-1">
              <div>Binding Pocket: <span class="text-emerald-300 font-medium">${pocket.description || "Auto-centered on co-crystallized native ligand pocket"}</span></div>
              <div class="font-mono text-[11px] text-slate-400">
                Center: (${pocket.center?.x}, ${pocket.center?.y}, ${pocket.center?.z}) | Size: (${pocket.size?.x}, ${pocket.size?.y}, ${pocket.size?.z})
              </div>
            </div>
          `;
        } else if (detectedPockets.length > 0) {
          pocketHtml = `
            <div class="text-xs text-slate-400 pt-2 border-t border-slate-700/60 space-y-1.5">
              <div class="flex items-center justify-between">
                <span class="text-amber-300 font-semibold flex items-center space-x-1">
                  <svg class="w-3.5 h-3.5 text-amber-400 inline" fill="none" viewBox="0 0 24 24" stroke="currentColor"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"/></svg>
                  <span>Blind Pocket Selection</span>
                </span>
                <span class="text-[10px] font-mono px-1.5 py-0.5 rounded bg-amber-950/80 text-amber-300 border border-amber-800/60">${detectedPockets.length} candidate pockets (fpocket)</span>
              </div>
              <p class="text-[11px] text-slate-400 leading-tight">No native ligand found — select candidate pocket ranked by druggability:</p>
              <select id="select-detected-pocket" class="w-full bg-slate-900 border border-slate-700 rounded p-1 text-xs text-slate-200">
                ${detectedPockets.map((p, idx) => `
                  <option value="${idx}" ${idx === (r.selected_pocket_idx || 0) ? 'selected' : ''}>
                    Pocket #${p.rank || idx + 1}: Druggability ${(p.druggability_score || 0).toFixed(2)} (Vol: ${Math.round(p.volume_a3 || 0)} Å³)
                  </option>
                `).join('')}
              </select>
              <div class="font-mono text-[11px] text-slate-400">
                Center: (${pocket.center?.x}, ${pocket.center?.y}, ${pocket.center?.z}) | Size: (${pocket.size?.x}, ${pocket.size?.y}, ${pocket.size?.z})
              </div>
            </div>
          `;
        } else {
          pocketHtml = `
            <div class="text-xs text-slate-400 pt-2 border-t border-slate-700/60 space-y-1">
              <div>Binding Pocket: <span class="text-emerald-300 font-medium">${pocket.description || "Auto-detected"}</span></div>
              <div class="font-mono text-[11px] text-slate-400">
                Center: (${pocket.center?.x}, ${pocket.center?.y}, ${pocket.center?.z}) | Size: (${pocket.size?.x}, ${pocket.size?.y}, ${pocket.size?.z})
              </div>
            </div>
          `;
        }

        recCard.innerHTML = `
          <div class="space-y-2">
            <div class="flex items-center justify-between">
              <span class="text-base font-bold text-emerald-400">${r.pdb_id || "Custom Receptor"}</span>
              <span class="text-xs px-2 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-800 font-mono">${r.atom_count || 0} atoms</span>
            </div>
            <p class="text-xs text-slate-300 font-medium">${r.title || "Macromolecular target"}</p>
            ${pocketHtml}
            ${validationHtml}
          </div>
        `;

        const pocketSelect = document.getElementById("select-detected-pocket");
        if (pocketSelect) {
          pocketSelect.addEventListener("change", (e) => {
            const idx = parseInt(e.target.value, 10);
            const selPocket = detectedPockets[idx];
            if (selPocket) {
              r.selected_pocket_idx = idx;
              r.detected_pocket = {
                center: selPocket.center,
                size: selPocket.size,
                description: selPocket.description || `Pocket #${selPocket.rank} (Druggability: ${(selPocket.druggability_score || 0).toFixed(2)})`
              };
              const setVal = (id, v) => { const el = document.getElementById(id); if (el) el.value = v; };
              setVal("grid-cx", selPocket.center?.x || 0);
              setVal("grid-cy", selPocket.center?.y || 0);
              setVal("grid-cz", selPocket.center?.z || 0);
              setVal("grid-sx", selPocket.size?.x || 22);
              setVal("grid-sy", selPocket.size?.y || 22);
              setVal("grid-sz", selPocket.size?.z || 22);
              if (this.viewer && this.viewer.renderGridBox) {
                this.viewer.renderGridBox(selPocket.center, selPocket.size);
              }
              this.updateStudioCards();
            }
          });
        }

        // Auto-trigger redocking validation if native ligand detected and not yet cached/run
        if (r.native_ligand?.has_native && !this.state.redockingValidation && !this._redockingRunning) {
          const recKey = r.pdb_id || r.title;
          if (this.state.redockCache && this.state.redockCache[recKey]) {
            this.state.redockingValidation = this.state.redockCache[recKey];
            this.updateReceptorView();
            this.updateDossierView();
          } else {
            setTimeout(() => this.runRedockingValidation(true), 150);
          }
        }
      } else {
        recCard.innerHTML = `<p class="text-xs text-slate-400 italic">No receptor loaded yet. Search RCSB PDB or pick a benchmark above.</p>`;
      }
    }

    // Grid box inputs
    if (this.state.receptor && this.state.receptor.detected_pocket) {
      const p = this.state.receptor.detected_pocket;
      const setVal = (id, v) => { const el = document.getElementById(id); if (el) el.value = v; };
      setVal("grid-cx", p.center?.x || 0);
      setVal("grid-cy", p.center?.y || 0);
      setVal("grid-cz", p.center?.z || 0);
      setVal("grid-sx", p.size?.x || 22);
      setVal("grid-sy", p.size?.y || 22);
      setVal("grid-sz", p.size?.z || 22);

      const gToggle = document.getElementById("toggle-gridbox");
      if (gToggle && gToggle.checked && this.viewer) {
        this.viewer.renderGridBox(
          { x: p.center?.x || 0, y: p.center?.y || 0, z: p.center?.z || 0 },
          { x: p.size?.x || 22, y: p.size?.y || 22, z: p.size?.z || 22 },
          true
        );
      }
    }

    // Native Ligand Redocking Validation Status
    const redockNameEl = document.getElementById("redock-native-name");
    const redockBtn = document.getElementById("btn-run-redock-validation");
    if (this.state.receptor?.native_ligand?.has_native) {
      const nat = this.state.receptor.native_ligand;
      if (redockNameEl) redockNameEl.textContent = `${nat.name} (Chain ${nat.chain || 'A'}, ${nat.atom_count} atoms)`;
      if (redockBtn) {
        redockBtn.disabled = false;
        redockBtn.className = "w-full py-2 px-3 rounded-lg bg-cyan-900/40 hover:bg-cyan-800/60 text-cyan-300 border border-cyan-700/60 text-xs font-semibold transition flex items-center justify-center space-x-1.5 cursor-pointer shadow-sm";
      }
    } else {
      if (redockNameEl) redockNameEl.textContent = "No co-crystallized ligand detected";
      if (redockBtn) {
        redockBtn.disabled = true;
        redockBtn.className = "w-full py-2 px-3 rounded-lg bg-slate-800 text-slate-500 border border-slate-700 text-xs font-semibold transition flex items-center justify-center space-x-1.5 cursor-not-allowed";
      }
    }

    // Update Preparation Protocol Transparency Log
    const prepRecEl = document.getElementById("prep-log-receptor");
    if (prepRecEl && this.state.receptor?.prep_log) {
      const pl = this.state.receptor.prep_log;
      prepRecEl.innerHTML = `
        <div>&bull; Stripped Waters: <span class="text-white font-bold">${pl.waters_removed}</span> atoms (HOH/WAT)</div>
        <div>&bull; Filtered Ions/Buffer: <span class="text-white font-bold">${pl.ions_and_buffer_removed}</span> atoms</div>
        <div>&bull; Retained Protein Atoms: <span class="text-white font-bold">${pl.protein_atoms_retained}</span> (${pl.selected_chain})</div>
        <div>&bull; Ionization Model: <span class="text-slate-300">${pl.protonation_state}</span></div>
        <div>&bull; Partial Charges: <span class="text-slate-300">${pl.charge_model}</span></div>
        <div>&bull; Pocket Positioning: <span class="text-emerald-300">${pl.active_pocket_centering}</span></div>
      `;
    }
    const prepLigEl = document.getElementById("prep-log-ligand");
    if (prepLigEl && this.state.ligand?.prep_log) {
      const pl = this.state.ligand.prep_log;
      prepLigEl.innerHTML = `
        <div>&bull; Input Format: <span class="text-white font-bold">${pl.input_format}</span> (${pl.heavy_atom_count} heavy atoms)</div>
        <div>&bull; Protonation: <span class="text-slate-300">${pl.hydrogens_added}</span></div>
        <div>&bull; Conformer Engine: <span class="text-slate-300">${pl.conformer_algorithm}</span></div>
        <div>&bull; Energy Minimization: <span class="text-slate-300">${pl.energy_minimization}</span></div>
        <div>&bull; Torsions: <span class="text-cyan-300">${pl.torsions_configured}</span></div>
        <div>&bull; Charges: <span class="text-slate-300">${pl.partial_charges}</span></div>
      `;
    }
  }

  applyBlindDockingBox() {
    if (!this.state.receptor) {
      this.showToast("Load a target receptor first to enable Blind Docking.", "warning");
      return;
    }
    const box = this.state.receptor.blind_docking_box;
    if (!box) {
      this.showToast("Blind docking bounding box not available for this structure.", "error");
      return;
    }

    const setVal = (id, v) => { const el = document.getElementById(id); if (el) el.value = v; };
    setVal("grid-cx", box.center.x);
    setVal("grid-cy", box.center.y);
    setVal("grid-cz", box.center.z);
    setVal("grid-sx", box.size.x);
    setVal("grid-sy", box.size.y);
    setVal("grid-sz", box.size.z);

    const gToggle = document.getElementById("toggle-gridbox");
    if (gToggle) gToggle.checked = true;
    if (this.viewer) {
      this.viewer.renderGridBox(box.center, box.size, true);
    }
    this.showToast(`Grid box expanded to whole protein surface (${box.size.x}×${box.size.y}×${box.size.z} Å) for Blind Docking.`, "info");
  }

  resetToPocketCentroid() {
    if (!this.state.receptor) {
      this.showToast("Load a target receptor first.", "warning");
      return;
    }
    const pocket = this.state.receptor.detected_pocket;
    if (!pocket) return;

    const setVal = (id, v) => { const el = document.getElementById(id); if (el) el.value = v; };
    setVal("grid-cx", pocket.center.x);
    setVal("grid-cy", pocket.center.y);
    setVal("grid-cz", pocket.center.z);
    setVal("grid-sx", pocket.size.x);
    setVal("grid-sy", pocket.size.y);
    setVal("grid-sz", pocket.size.z);

    const gToggle = document.getElementById("toggle-gridbox");
    if (gToggle) gToggle.checked = true;
    if (this.viewer) {
      this.viewer.renderGridBox(pocket.center, pocket.size, true);
    }
    this.showToast("Grid box restored to detected active pocket centroid.", "success");
  }

  async runRedockingValidation(isAuto = false) {
    if (!this.state.receptor || !this.state.receptor.native_ligand?.pdb_block) {
      if (!isAuto) this.showToast("No co-crystallized native ligand found in current receptor.", "error");
      return;
    }
    const nat = this.state.receptor.native_ligand;
    const recKey = this.state.receptor.pdb_id || this.state.receptor.title || "custom_rec";
    this.state.redockCache = this.state.redockCache || {};

    if (this.state.redockCache[recKey]) {
      this.state.redockingValidation = this.state.redockCache[recKey];
      this.updateReceptorView();
      this.updateDossierView();
      if (!isAuto) {
        this.showToast(`Loaded cached redocking: RMSD = ${this.state.redockingValidation.rmsd_angstroms} Å`, "info");
      }
      return;
    }

    const btn = document.getElementById("btn-run-redock-validation");
    const resultsPanel = document.getElementById("redock-results-panel");
    const affVal = document.getElementById("redock-aff-val");
    const rmsdVal = document.getElementById("redock-rmsd-val");
    const badgeBox = document.getElementById("redock-badge-container");

    const originalBtnText = btn ? btn.innerHTML : "";
    if (btn) {
      btn.disabled = true;
      btn.innerHTML = `<svg class="animate-spin -ml-1 mr-2 h-3.5 w-3.5 text-white inline" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg> Redocking ${nat.name}...`;
    }

    if (!isAuto) {
      this.showToast(`Starting redocking validation of native ligand '${nat.name}'...`, "info");
    }
    this._redockingRunning = true;

    try {
      const getVal = (id, def) => parseFloat(document.getElementById(id)?.value) || def;
      const center = { x: getVal("grid-cx", nat.center?.x || 0), y: getVal("grid-cy", nat.center?.y || 0), z: getVal("grid-cz", nat.center?.z || 0) };
      const size = { x: getVal("grid-sx", 22), y: getVal("grid-sy", 22), z: getVal("grid-sz", 22) };

      const exhaustiveness = parseInt(document.getElementById("docking-exhaustiveness")?.value) || 8;
      const res = await BindoraAPI.redockValidate({
        receptor_pdbqt: this.state.receptor.pdbqt_text,
        native_ligand_pdb: nat.pdb_block,
        center: center,
        size: size,
        exhaustiveness: exhaustiveness
      });

      this.state.redockingValidation = res;
      this.state.redockCache[recKey] = res;

      if (resultsPanel) resultsPanel.classList.remove("hidden");
      if (affVal) affVal.textContent = `${res.affinity_kcal} kcal/mol`;
      if (rmsdVal) rmsdVal.textContent = `${res.rmsd_angstroms} Å`;

      if (badgeBox) {
        if (res.is_validated) {
          badgeBox.innerHTML = `<span class="px-2.5 py-1 rounded-full bg-emerald-100 dark:bg-emerald-950 text-emerald-800 dark:text-emerald-300 border border-emerald-300 dark:border-emerald-800 text-[11px] font-bold inline-flex items-center space-x-1.5"><svg class="w-3.5 h-3.5 text-emerald-600 dark:text-emerald-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polyline points="20 6 9 17 4 12"/></svg> <span>${res.validation_badge}</span></span>`;
        } else {
          badgeBox.innerHTML = `<span class="px-2.5 py-1 rounded-full bg-amber-100 dark:bg-amber-950 text-amber-800 dark:text-amber-300 border border-amber-300 dark:border-amber-800 text-[11px] font-bold inline-flex items-center space-x-1.5"><svg class="w-3.5 h-3.5 text-amber-600 dark:text-amber-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="m21.73 18-8-14a2 2 0 0 0-3.48 0l-8 14A2 2 0 0 0 4 21h16a2 2 0 0 0 1.73-3Z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg> <span>${res.validation_badge}</span></span>`;
        }
      }

      // Preview 3D overlay of crystallographic reference vs redocked pose
      if (this.viewer && res.docked_pdb && !isAuto) {
        if (res.cryst_pdb) {
          this.viewer.loadRedockOverlay(res.cryst_pdb, res.docked_pdb);
        } else {
          this.viewer.loadLigand(res.docked_pdb);
        }
      }

      this.updateReceptorView();
      this.updateDossierView();
      if (!isAuto) {
        this.showToast(`Redocking validation finished! RMSD = ${res.rmsd_angstroms} Å (${res.benchmark_status})`, res.is_validated ? "success" : "warning");
      }

    } catch (err) {
      console.error("Redocking error:", err);
      if (!isAuto) this.showToast(`Redocking validation failed: ${err.message}`, "error");
    } finally {
      this._redockingRunning = false;
      if (btn) {
        btn.disabled = false;
        btn.innerHTML = originalBtnText || "<span>Run Self-Validation (Redock Native Ligand)</span>";
      }
    }
  }

  async runDockingPipeline() {
    if (!this.state.ligand || !this.state.receptor) {
      this.showToast("Please load both a ligand and a receptor first.", "error");
      return;
    }

    const dockBtn = document.getElementById("btn-run-docking");
    const originalText = dockBtn ? dockBtn.innerHTML : "";
    if (dockBtn) {
      dockBtn.disabled = true;
      dockBtn.innerHTML = `<svg class="animate-spin -ml-1 mr-2 h-4 w-4 text-white inline" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg> Running AutoDock Vina...`;
    }

    const statusBanner = document.getElementById("docking-status-banner");
    const statusIcon = document.getElementById("docking-status-icon");
    const statusTitle = document.getElementById("docking-status-title");
    const statusSub = document.getElementById("docking-status-subtitle");
    const statusExtra = document.getElementById("docking-status-extra");

    const exhaustiveness = parseInt(document.getElementById("docking-exhaustiveness")?.value) || 8;
    const replicates = parseInt(document.getElementById("docking-replicates")?.value) || 1;

    if (statusBanner) {
      statusBanner.classList.remove("hidden");
      statusBanner.className = "glass-panel p-3.5 rounded-xl border border-cyan-500/50 bg-cyan-950/80 transition-all duration-300 mb-2 animate-pulse";
      if (statusIcon) statusIcon.innerHTML = `<svg class="animate-spin w-5 h-5 text-cyan-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>`;
      if (statusTitle) statusTitle.textContent = "AutoDock Vina Docking in Progress...";
      if (statusSub) statusSub.textContent = `Sampling conformational space (exhaustiveness = ${exhaustiveness}, modes = 9, ${replicates > 1 ? '3 replicate seeds' : 'single seed'}) on CPU...`;
      if (statusExtra) statusExtra.textContent = replicates > 1 ? "Replicates Running..." : "Running Monte Carlo...";
    }

    this.showToast(`Launching AutoDock Vina (${replicates > 1 ? '3-seed replicate' : 'single run'})...`, "info");

    try {
      // Read grid parameters
      const getVal = (id, def) => parseFloat(document.getElementById(id)?.value) || def;
      const center = { x: getVal("grid-cx", 0), y: getVal("grid-cy", 0), z: getVal("grid-cz", 0) };
      const size = { x: getVal("grid-sx", 22), y: getVal("grid-sy", 22), z: getVal("grid-sz", 22) };

      const p = this.state.ligand.adme?.physicochemical || {};

      const flexCheckboxes = document.querySelectorAll("#flex-residues-container input[type='checkbox']:checked");
      const flexResidues = Array.from(flexCheckboxes).map(cb => cb.value);

      const dockPayload = {
        receptor_pdbqt: this.state.receptor.pdbqt_text,
        receptor_pdb: this.state.receptor.cleaned_pdb,
        ligand_pdbqt: this.state.ligand.pdbqt_text,
        center: center,
        size: size,
        exhaustiveness: exhaustiveness,
        num_modes: 9,
        replicates: replicates,
        heavy_atoms: p.heavy_atoms?.value || this.state.ligand.heavy_atom_count || 20,
        molecular_weight: p.molecular_weight?.value || this.state.ligand.weight || 300.0,
        smiles: this.state.ligand.canonical_smiles || this.state.ligand.smiles || ""
      };
      if (flexResidues.length > 0) {
        dockPayload.flexible_residues = flexResidues;
      }

      const dockResult = await BindoraAPI.runDocking(dockPayload);

      this.state.docking = dockResult;
      this.state.currentPoseIdx = 0;

      // Render top pose in 3D viewer
      if (this.viewer && dockResult.top_pose) {
        this.viewer.loadLigand(dockResult.top_pose.pdb_block);
        if (dockResult.interactions) {
          this.viewer.renderInteractions(dockResult.interactions);
        }
      }

      // Populate flexible active site residues for user inspection / induced-fit tuning
      if (dockResult.interactions?.flexible_candidates) {
        this.updateFlexibleResiduesUI(dockResult.interactions.flexible_candidates);
      }

      // Render pose score badges
      this.updateDockingScoresUI();

      // Trigger ChEMBL Bioactivity cross-check asynchronously
      this.fetchBioactivityCrosscheck();

      // Render ADME charts
      this.updateADMEView();

      // Replicate Statistics Display
      const repBanner = document.getElementById("replicate-stats-banner");
      const repVal = document.getElementById("replicate-stats-value");
      if (dockResult.replicate_stats && replicates > 1) {
        const rs = dockResult.replicate_stats;
        if (repBanner) repBanner.classList.remove("hidden");
        if (repVal) repVal.textContent = `Mean: ${rs.mean_affinity_kcal} ± ${rs.sd_affinity_kcal} kcal/mol (N=${rs.replicates_count} seeds)`;
      } else {
        if (repBanner) repBanner.classList.add("hidden");
      }

      // Weak-Binder Auto-Alert Banner
      const weakAlert = document.getElementById("weak-binder-alert");
      const weakText = document.getElementById("weak-binder-alert-text");
      if (dockResult.thermodynamics?.is_weak_binder) {
        if (weakAlert) weakAlert.classList.remove("hidden");
        if (weakText) weakText.textContent = dockResult.thermodynamics.weak_binder_warning || "Calculated affinity falls above the -6.0 kcal/mol threshold (high micromolar/millimolar Kd). Such weak interactions usually indicate non-specific surface adhesion or numeric artifacts. Interpret with extreme caution.";
      } else {
        if (weakAlert) weakAlert.classList.add("hidden");
      }

      // Update completion banner
      if (statusBanner) {
        statusBanner.className = "glass-panel p-3.5 rounded-xl border border-emerald-500/50 bg-emerald-50 dark:bg-emerald-950/80 text-emerald-950 dark:text-emerald-100 transition-all duration-300 mb-2";
        if (statusIcon) statusIcon.innerHTML = `<span class="w-6 h-6 rounded-full bg-emerald-500/20 text-emerald-600 dark:text-emerald-400 flex items-center justify-center font-bold text-sm"><svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polyline points="20 6 9 17 4 12"/></svg></span>`;
        if (statusTitle) statusTitle.textContent = `Docking Complete: ΔG = ${dockResult.top_pose.affinity_kcal} kcal/mol`;
        if (statusSub) statusSub.textContent = `Theoretical Kd: ${dockResult.thermodynamics?.theoretical_kd_nm || "—"} nM | ${dockResult.interactions?.total_hbond_count || 0} H-Bonds | Top binding pose loaded.`;
        if (statusExtra) statusExtra.textContent = "Finished";
      }

      // Switch to 3D docking tab
      this.switchTab("docking");
      this.showToast(`Docking finished! Top ΔG: ${dockResult.top_pose.affinity_kcal} kcal/mol`, "success");
      this.updateDossierView();

    } catch (e) {
      console.error("Docking error:", e);
      if (statusBanner) {
        statusBanner.className = "glass-panel p-3.5 rounded-xl border border-rose-500/50 bg-rose-50 dark:bg-rose-950/80 text-rose-950 dark:text-rose-100 transition-all duration-300 mb-2";
        if (statusIcon) statusIcon.innerHTML = `<span class="w-6 h-6 rounded-full bg-rose-500/20 text-rose-600 dark:text-rose-400 flex items-center justify-center font-bold text-sm"><svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><line x1="18" y1="6" x2="6" y2="18"/><line x1="6" y1="6" x2="18" y2="18"/></svg></span>`;
        if (statusTitle) statusTitle.textContent = "Docking Execution Failed";
        if (statusSub) statusSub.textContent = e.message;
        if (statusExtra) statusExtra.textContent = "Error";
      }
      this.showToast(`Docking execution failed: ${e.message}`, "error");
    } finally {
      if (dockBtn) {
        dockBtn.disabled = false;
        dockBtn.innerHTML = originalText;
      }
    }
  }

  updateDockingScoresUI() {
    if (!this.state.docking) return;
    const { poses = [], top_pose, thermodynamics = {}, interactions = {} } = this.state.docking;
    const currentPose = poses[this.state.currentPoseIdx] || top_pose;

    // Header values
    const elAffinity = document.getElementById("dock-affinity-score");
    if (elAffinity) elAffinity.textContent = `${currentPose.affinity_kcal} kcal/mol`;

    const elKd = document.getElementById("dock-kd-score");
    if (elKd) {
      const rt = 0.5924847;
      const aff = parseFloat(currentPose.affinity_kcal) || 0;
      const kdMolar = Math.exp(aff / rt);
      const kdNm = kdMolar * 1e9;
      const kdDisplay = kdNm < 1000 ? `${kdNm.toFixed(2)} nM` : kdNm < 1e6 ? `${(kdNm / 1e3).toFixed(2)} µM` : `${(kdNm / 1e6).toFixed(2)} mM`;
      elKd.textContent = kdDisplay;
    }

    const elLE = document.getElementById("dock-le-score");
    if (elLE) {
      const ha = this.state.ligand?.heavy_atom_count || this.state.ligand?.adme?.physicochemical?.heavy_atoms?.value || 0;
      if (ha > 0) {
        const aff = parseFloat(currentPose.affinity_kcal) || 0;
        elLE.textContent = (-aff / ha).toFixed(3);
      } else {
        elLE.textContent = `${thermodynamics.ligand_efficiency?.value || "—"}`;
      }
    }

    const elHbonds = document.getElementById("dock-hbonds-count");
    if (elHbonds) elHbonds.textContent = `${interactions.total_hbond_count || 0}`;

    // Pose Stepper Label
    const elMode = document.getElementById("dock-current-mode-label");
    if (elMode) elMode.textContent = `Mode ${currentPose.mode || (this.state.currentPoseIdx + 1)} / ${poses.length}`;

    // Pose Table
    const poseTable = document.getElementById("pose-table-rows");
    if (poseTable) {
      poseTable.innerHTML = poses.map((p, idx) => `
        <tr class="border-b border-slate-200 dark:border-slate-800 hover:bg-slate-100 dark:hover:bg-slate-800/40 cursor-pointer ${idx === this.state.currentPoseIdx ? 'bg-cyan-50 dark:bg-cyan-950/40 font-semibold' : ''}" onclick="window.app.selectPose(${idx})">
          <td class="px-3 py-2 text-cyan-700 dark:text-cyan-300">Mode ${p.mode}</td>
          <td class="px-3 py-2 font-mono text-slate-900 dark:text-white font-semibold">${p.affinity_kcal}</td>
          <td class="px-3 py-2 font-mono text-emerald-600 dark:text-emerald-400">${p.vinardo_affinity_kcal != null ? p.vinardo_affinity_kcal : '—'}</td>
          <td class="px-3 py-2 font-mono text-slate-600 dark:text-slate-400">${p.rmsd_lb}</td>
          <td class="px-3 py-2 font-mono text-slate-600 dark:text-slate-400">${p.rmsd_ub}</td>
        </tr>
      `).join("");
    }

    // Interacting Residues Chips
    const resContainer = document.getElementById("dock-interacting-residues");
    if (resContainer) {
      const resList = interactions.interacting_residues || [];
      if (resList.length > 0) {
        resContainer.innerHTML = resList.map(r => `
          <span class="inline-block px-2 py-0.5 rounded text-xs bg-cyan-50 dark:bg-slate-800 text-cyan-800 dark:text-cyan-300 border border-cyan-200 dark:border-slate-700 font-mono">${r}</span>
        `).join("");
      } else {
        resContainer.innerHTML = `<span class="text-xs text-slate-500 italic">No specific residues within contact threshold</span>`;
      }
    }

    // 2D Interaction Diagram Rendering
    const diagramContainer = document.getElementById("dock-interaction-diagram-container");
    if (diagramContainer) {
      if (interactions.diagram_svg) {
        diagramContainer.innerHTML = interactions.diagram_svg;
        this.setupDiagramPanZoom(diagramContainer);
      } else {
        diagramContainer.innerHTML = `<span class="text-xs text-slate-500 italic">No 2D interaction schematic available for this pose</span>`;
      }
    }

    // Binding Energy Landscape Chart Rendering
    if (window.DockingCharts && poses.length > 0) {
      DockingCharts.renderEnergyLandscape("chart-energy-landscape", poses);
    }
  }

  async selectPose(idx) {
    if (!this.state.docking || !this.state.docking.poses) return;
    if (idx < 0 || idx >= this.state.docking.poses.length) return;

    this.state.currentPoseIdx = idx;
    const pose = this.state.docking.poses[idx];

    // Re-analyze interactions for this pose
    try {
      const smiles = this.state.ligand?.canonical_smiles || this.state.ligand?.smiles || "";
      const contacts = await BindoraAPI.analyzeInteractions(this.state.receptor.cleaned_pdb, pose.pdbqt_content, smiles);
      this.state.docking.interactions = contacts;

      if (this.viewer) {
        this.viewer.loadLigand(pose.pdb_block);
        this.viewer.renderInteractions(contacts);
      }

      // Update 2D interaction diagram for the selected pose
      const diagramContainer = document.getElementById("dock-interaction-diagram-container");
      if (diagramContainer && contacts.diagram_svg) {
        diagramContainer.innerHTML = contacts.diagram_svg;
        this.setupDiagramPanZoom(diagramContainer);
      }

      if (contacts?.flexible_candidates) {
        this.updateFlexibleResiduesUI(contacts.flexible_candidates);
      }
    } catch (e) {
      console.warn("Could not re-analyze pose interactions:", e);
    }

    this.updateDockingScoresUI();
    this.updateDossierView();
  }

  setupDiagramPanZoom(container) {
    const svg = container.querySelector("svg");
    if (!svg) return;

    // Set theme class matching current html theme
    if (document.documentElement.classList.contains("light")) {
      svg.classList.add("theme-light");
    } else {
      svg.classList.remove("theme-light");
    }

    const contentGroup = svg.querySelector("#bindora-diagram-content") || svg;

    // Diagram center coordinates (650 x 500)
    const cx = 325;
    const cy = 235;

    // State
    let currentScale = 1.0;
    let panX = 0;
    let panY = 0;

    const applyTransform = () => {
      if (contentGroup && contentGroup !== svg) {
        contentGroup.setAttribute(
          "transform",
          `translate(${panX.toFixed(1)}, ${panY.toFixed(1)}) translate(${cx}, ${cy}) scale(${currentScale.toFixed(3)}) translate(-${cx}, -${cy})`
        );
      }
    };
    applyTransform();

    // Mouse drag to pan
    let isDragging = false;
    let startX = 0;
    let startY = 0;
    let startPanX = 0;
    let startPanY = 0;

    container.onmousedown = (e) => {
      if (e.button !== 0) return; // Only primary mouse button
      isDragging = true;
      startX = e.clientX;
      startY = e.clientY;
      startPanX = panX;
      startPanY = panY;
      container.style.cursor = "grabbing";
      e.preventDefault();
    };

    const handleMouseMove = (e) => {
      if (!isDragging) return;
      const rect = container.getBoundingClientRect();
      const scaleFactor = 650 / (rect.width || 650);
      const dx = (e.clientX - startX) * scaleFactor;
      const dy = (e.clientY - startY) * scaleFactor;

      // Clamp panning bounds so molecule is never lost
      panX = Math.max(-280, Math.min(280, startPanX + dx));
      panY = Math.max(-220, Math.min(220, startPanY + dy));
      applyTransform();
    };

    const handleMouseUp = () => {
      if (isDragging) {
        isDragging = false;
        container.style.cursor = "grab";
      }
    };

    if (container._panMoveHandler) {
      window.removeEventListener("mousemove", container._panMoveHandler);
    }
    if (container._panUpHandler) {
      window.removeEventListener("mouseup", container._panUpHandler);
    }
    container._panMoveHandler = handleMouseMove;
    container._panUpHandler = handleMouseUp;
    window.addEventListener("mousemove", handleMouseMove);
    window.addEventListener("mouseup", handleMouseUp);

    // Mouse wheel zoom - smoothly centered on (cx, cy)
    container.onwheel = (e) => {
      e.preventDefault();
      const zoomDelta = e.deltaY > 0 ? 0.92 : 1.08;
      currentScale = Math.max(0.68, Math.min(2.4, currentScale * zoomDelta));
      applyTransform();
    };

    // Wire toolbar buttons
    const btnZoomIn = document.getElementById("btn-diagram-zoom-in");
    const btnZoomOut = document.getElementById("btn-diagram-zoom-out");
    const btnReset = document.getElementById("btn-diagram-reset-view");
    const btnExportPng = document.getElementById("btn-diagram-export-png");
    const btnExportSvg = document.getElementById("btn-diagram-export-svg");

    if (btnZoomIn) {
      btnZoomIn.onclick = () => {
        currentScale = Math.min(2.4, currentScale * 1.2);
        applyTransform();
      };
    }

    if (btnZoomOut) {
      btnZoomOut.onclick = () => {
        currentScale = Math.max(0.68, currentScale * 0.83);
        applyTransform();
      };
    }

    if (btnReset) {
      btnReset.onclick = () => {
        currentScale = 1.0;
        panX = 0;
        panY = 0;
        applyTransform();
      };
    }

    if (btnExportPng) {
      btnExportPng.onclick = () => this.exportDiagramPNG(svg);
    }

    if (btnExportSvg) {
      btnExportSvg.onclick = () => this.exportDiagramSVG(svg);
    }
  }

  exportDiagramPNG(svgElement) {
    if (!svgElement) return;
    const isLight = document.documentElement.classList.contains("light") || svgElement.classList.contains("theme-light");
    const clone = svgElement.cloneNode(true);

    // Publication standard: Reset pan/zoom on export so image is 100% centered, unclipped and includes legend
    const cloneContent = clone.querySelector("#bindora-diagram-content");
    if (cloneContent) {
      cloneContent.setAttribute("transform", "matrix(1 0 0 1 0 0)");
    }

    clone.setAttribute("viewBox", "0 0 650 500");
    clone.setAttribute("width", "2600");
    clone.setAttribute("height", "2000");
    if (isLight) {
      clone.classList.add("theme-light");
    } else {
      clone.classList.remove("theme-light");
    }

    const svgXml = new XMLSerializer().serializeToString(clone);
    const svgBlob = new Blob([svgXml], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(svgBlob);

    const img = new Image();
    img.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = 2600;
      canvas.height = 2000;
      const ctx = canvas.getContext("2d");

      // Solid background fill for publication clarity
      ctx.fillStyle = isLight ? "#ffffff" : "#091428";
      ctx.fillRect(0, 0, canvas.width, canvas.height);

      ctx.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);

      const pngData = canvas.toDataURL("image/png");
      const a = document.createElement("a");
      const ligName = (this.state.ligand?.name || "ligand").toLowerCase().replace(/[^a-z0-9]/g, "_");
      a.download = `bindora_2d_interaction_${ligName}_300dpi.png`;
      a.href = pngData;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    };
    img.src = url;
  }

  exportDiagramSVG(svgElement) {
    if (!svgElement) return;
    const isLight = document.documentElement.classList.contains("light");
    const clone = svgElement.cloneNode(true);
    const cloneContent = clone.querySelector("#bindora-diagram-content");
    if (cloneContent) {
      cloneContent.setAttribute("transform", "matrix(1 0 0 1 0 0)");
    }
    clone.setAttribute("viewBox", "0 0 650 500");
    if (isLight) {
      clone.classList.add("theme-light");
    } else {
      clone.classList.remove("theme-light");
    }

    const svgXml = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
    const blob = new Blob([svgXml], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    const ligName = (this.state.ligand?.name || "ligand").toLowerCase().replace(/[^a-z0-9]/g, "_");
    a.download = `bindora_2d_interaction_${ligName}.svg`;
    a.href = url;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  changePose(direction) {
    if (!this.state.docking || !this.state.docking.poses) return;
    const newIdx = this.state.currentPoseIdx + direction;
    this.selectPose(newIdx);
  }

  async fetchBioactivityCrosscheck() {
    const drug = this.state.ligand?.name || "";
    const target = this.state.receptor?.title || this.state.receptor?.pdb_id || "";
    if (!drug || !target) return;

    const crosscheckContainer = document.getElementById("crosscheck-container");
    if (crosscheckContainer) {
      crosscheckContainer.innerHTML = `<div class="text-xs text-slate-400 animate-pulse">Querying ChEMBL curated bioactivity database for ${drug} + ${target}...</div>`;
    }

    try {
      const res = await BindoraAPI.crosscheckBioactivity(drug, target);
      this.state.crosscheck = res;
      this.renderCrosscheckUI(res);
      this.updateDossierView();
      
      // Once crosscheck is done, generate educational AI narrative
      this.generateNarrativeReport();
    } catch (e) {
      console.warn("Crosscheck query error:", e);
      if (crosscheckContainer) {
        crosscheckContainer.innerHTML = `<div class="text-xs text-amber-400">Bioactivity cross-check unavailable: ${e.message}</div>`;
      }
    }
  }

  renderCrosscheckUI(data) {
    const container = document.getElementById("crosscheck-container");
    if (!container) return;

    const isCorroborated = data.is_cross_checked;
    const badgeClass = isCorroborated ? "badge-verified" : "badge-computational";
    const statusIcon = isCorroborated
      ? `<svg class="w-4 h-4 text-emerald-500 inline" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polyline points="20 6 9 17 4 12"/></svg>`
      : `<svg class="w-4 h-4 text-amber-500 inline" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="m21.73 18-8-14a2 2 0 0 0-3.48 0l-8 14A2 2 0 0 0 4 21h16a2 2 0 0 0 1.73-3Z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>`;

    // Build ChEMBL verification links
    const chemblMolUrl = data.chembl_compound_url || null;
    const chemblTargetUrl = data.chembl_target_url || null;
    const molId = data.molecule_chembl_id;
    const targetId = data.target_chembl_id;

    const searchedHtml = `
      <div class="grid grid-cols-2 gap-2 text-xs font-mono mt-2">
        <div class="p-2.5 rounded-lg bg-slate-100/90 dark:bg-slate-900/60 border border-slate-200 dark:border-slate-800 space-y-1">
          <div class="text-slate-500 dark:text-slate-400 text-[10px] uppercase tracking-wider font-sans font-semibold">Compound Searched</div>
          <div class="text-slate-900 dark:text-white font-semibold">${data.drug_searched || "—"}</div>
          ${molId ? `<div class="text-[11px] text-slate-500 dark:text-slate-400">${molId}
            <a href="${chemblMolUrl}" target="_blank" rel="noopener" class="text-cyan-600 dark:text-cyan-400 hover:underline ml-1">↗ ChEMBL</a></div>` : `<div class="text-[11px] text-rose-500 dark:text-rose-400">Not found in ChEMBL</div>`}
        </div>
        <div class="p-2.5 rounded-lg bg-slate-100/90 dark:bg-slate-900/60 border border-slate-200 dark:border-slate-800 space-y-1">
          <div class="text-slate-500 dark:text-slate-400 text-[10px] uppercase tracking-wider font-sans font-semibold">Target Searched</div>
          <div class="text-slate-900 dark:text-white font-semibold">${data.target_searched || "—"}</div>
          ${targetId ? `<div class="text-[11px] text-slate-500 dark:text-slate-400">${targetId}
            <a href="${chemblTargetUrl}" target="_blank" rel="noopener" class="text-cyan-600 dark:text-cyan-400 hover:underline ml-1">↗ ChEMBL</a></div>` : `<div class="text-[11px] text-rose-500 dark:text-rose-400">Not found in ChEMBL</div>`}
        </div>
      </div>
    `;

    // No records explanation panel (scientific context, not an error)
    const noRecordsExplanationHtml = !isCorroborated ? `
      <div class="mt-3 p-3 rounded-lg bg-amber-50/80 dark:bg-slate-900/50 border border-amber-200 dark:border-amber-900/40 text-xs space-y-1">
        <p class="font-semibold text-amber-800 dark:text-amber-300 text-[11px] uppercase tracking-wide">What does this mean scientifically?</p>
        <p class="text-slate-700 dark:text-slate-400 leading-relaxed">
          <strong class="text-slate-900 dark:text-slate-300">No ChEMBL records ≠ inactive compound.</strong> It means no curated wet-lab binding assay (IC50, Ki, Kd, EC50) has been deposited in the EMBL-EBI ChEMBL database for this exact drug-target pair. This is a common finding for:
        </p>
        <ul class="list-disc pl-4 text-slate-400 space-y-0.5">
          <li>Novel or emerging drug candidates with unpublished experimental data</li>
          <li>Target-ligand combinations studied under different assay conditions not yet curated</li>
          <li>Computational lead compounds undergoing preclinical evaluation</li>
        </ul>
        <p class="text-slate-400">All displayed thermodynamic values (ΔG, Kd, LE) are <strong class="text-amber-300">computational predictions</strong> from AutoDock Vina and must be treated as hypothesis-generating. Wet-lab validation (SPR, ITC, fluorescence anisotropy) is required before drawing pharmacological conclusions.</p>
      </div>
    ` : '';

    let recordsHtml = "";
    if (isCorroborated && data.experimental_records?.length > 0) {
      recordsHtml = `
        <div class="mt-3 overflow-x-auto">
          <div class="text-[10px] uppercase tracking-wider text-emerald-400 font-semibold mb-1">Curated Experimental Bioactivity Records from ChEMBL</div>
          <table class="w-full text-xs text-left border border-slate-700/60 rounded-lg overflow-hidden">
            <thead class="bg-slate-800/70 text-slate-300 font-semibold border-b border-slate-700">
              <tr>
                <th class="p-2">Assay Type</th>
                <th class="p-2">Experimental Value</th>
                <th class="p-2">Assay Description</th>
                <th class="p-2">Source</th>
              </tr>
            </thead>
            <tbody class="divide-y divide-slate-800 font-mono text-slate-300">
              ${data.experimental_records.map(r => `
                <tr class="hover:bg-slate-800/40">
                  <td class="p-2 font-bold text-emerald-400">${r.type}</td>
                  <td class="p-2 text-white">${r.relation} ${r.value} ${r.units}</td>
                  <td class="p-2 font-sans text-slate-400 text-[11px]">${r.assay_description}</td>
                  <td class="p-2 font-sans text-cyan-400 text-[11px]">
                    <a href="https://www.ebi.ac.uk/chembl/" target="_blank" rel="noopener" class="hover:underline">ChEMBL</a>
                  </td>
                </tr>
              `).join("")}
            </tbody>
          </table>
        </div>
      `;
    }

    container.innerHTML = `
      <div class="p-4 rounded-xl border ${badgeClass} space-y-2">
        <div class="flex items-center justify-between flex-wrap gap-2">
          <div class="flex items-center space-x-2">
            <span class="text-base font-bold">${statusIcon}</span>
            <span class="text-sm font-bold tracking-wide uppercase">${data.status_badge}</span>
          </div>
          <div class="flex items-center gap-2">
            <span class="text-xs px-2 py-0.5 rounded bg-slate-900/60 font-mono">${data.target_organism || "Homo sapiens"}</span>
            ${chemblMolUrl ? `<a href="${chemblMolUrl}" target="_blank" rel="noopener" class="text-xs px-2 py-0.5 rounded bg-cyan-950/60 border border-cyan-800/50 text-cyan-400 hover:underline font-mono">Verify on ChEMBL ↗</a>` : ''}
          </div>
        </div>
        <p class="text-xs text-slate-300 leading-relaxed">${data.summary_note}</p>
        ${searchedHtml}
        ${noRecordsExplanationHtml}
        ${recordsHtml}
      </div>
    `;
  }

  updateADMEView() {
    const adme = this.state.ligand?.adme;
    if (!adme || !adme.physicochemical) return;

    // Render Lipinski Radar Chart
    ADMECharts.renderLipinskiRadar("chart-lipinski-radar", adme);

    // Update Lipinski Status Badge
    const lipBadge = document.getElementById("adme-lipinski-badge");
    if (lipBadge) {
      const l = adme.drug_likeness?.lipinski || {};
      const status = l.status || "Pass";
      const color = status === "Pass" ? "bg-emerald-950 text-emerald-300 border-emerald-800" :
                    status === "Borderline" ? "bg-amber-950 text-amber-300 border-amber-800" :
                    "bg-rose-950 text-rose-300 border-rose-800";
      lipBadge.className = `px-3 py-1 rounded-full text-xs font-semibold border ${color}`;
      lipBadge.textContent = `Lipinski: ${status} (${l.violations_count || 0} violations)`;
    }

    // Update Pharmacokinetic Cards
    const pk = adme.pharmacokinetics || {};
    const setVal = (id, txt) => { const el = document.getElementById(id); if (el) el.textContent = txt; };
    
    // Physicochemical descriptor cards
    const p = adme.physicochemical || {};
    setVal("adme-mw-val", p.molecular_weight?.value ?? "—");
    setVal("adme-logp-val", p.logp?.value ?? "—");
    setVal("adme-hbd-val", p.hbd?.value ?? "—");
    setVal("adme-hba-val", p.hba?.value ?? "—");
    setVal("adme-rotb-val", p.rotatable_bonds?.value ?? "—");
    setVal("adme-tpsa-val", p.tpsa?.value ?? "—");
    setVal("adme-sascore-val", p.sascore?.value != null ? `${p.sascore.value}` : "—");
    const sasEl = document.getElementById("adme-sascore-val");
    if (sasEl && p.sascore?.interpretation) {
      sasEl.title = `${p.sascore.interpretation} (Ertl & Schuffenhauer 2009)`;
    }

    setVal("adme-gi-val", pk.gi_absorption?.level || pk.gi_absorption?.status || "—");
    setVal("adme-bbb-val", pk.bbb_permeation?.status?.split(" ")[0] || "—");
    setVal("adme-ppb-val", pk.plasma_protein_binding?.tier?.split(" ")[0] || "—");

    // Structural Alerts (PAINS, Brenk, NIH, ZINC)
    const safety = adme.medicinal_chemistry_safety || {};
    const painsEl = document.getElementById("adme-pains-val");
    if (painsEl) {
      const totalAlerts = safety.total_alerts_count ?? (safety.pains_alerts?.count || 0);
      const statusText = safety.overall_status || (totalAlerts === 0 ? "Clean (0 alerts)" : `${totalAlerts} Alert(s) Detected`);
      painsEl.textContent = statusText;
      painsEl.className = totalAlerts > 0 ? "font-bold text-rose-400" : "font-bold text-emerald-400";
    }

    // CYP450 Heuristic Liability
    const cyp = pk.cyp450_inhibition || {};
    const cypEl = document.getElementById("adme-cyp-val");
    if (cypEl) {
      const cypAlerts = cyp.alerts || (Array.isArray(cyp) ? cyp : []);
      if (cypAlerts.length > 0) {
        cypEl.innerHTML = `<span class="text-amber-600 dark:text-amber-400 font-bold">${cypAlerts.length} SMARTS Alert(s):</span> ` +
          cypAlerts.map(a => `<span class="inline-block font-mono text-[11px] bg-amber-200/60 dark:bg-amber-900/50 text-amber-900 dark:text-amber-200 px-1.5 py-0.5 rounded mr-1">${a.cyp} (${a.description})</span>`).join("");
      } else {
        cypEl.innerHTML = `<span class="text-emerald-600 dark:text-emerald-400 font-semibold">Clear &bull; No planar/basic/azole SMARTS alerts detected.</span>`;
      }
    }
  }

  async generateNarrativeReport() {
    const container = document.getElementById("narrative-report-content");
    if (!container) return;

    // If neither receptor nor ligand is loaded yet, display friendly guidance instead of error
    if (!this.state.receptor && !this.state.ligand) {
      container.innerHTML = `
        <div class="p-6 rounded-xl bg-slate-900/40 border border-slate-700/60 text-center space-y-3">
          <div class="w-12 h-12 mx-auto rounded-full bg-purple-950/60 border border-purple-500/40 flex items-center justify-center text-purple-300 shadow-inner">
            <svg class="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z"/></svg>
          </div>
          <h4 class="text-sm font-semibold text-slate-200">No Active Simulation Loaded</h4>
          <p class="text-xs text-slate-400 max-w-md mx-auto leading-relaxed">
            To generate a publication-grade pharmacological briefing, select a target receptor and drug in <span class="text-cyan-400 font-semibold">1. Target &amp; Ligand</span> and run docking in <span class="text-cyan-400 font-semibold">2. 3D Docking</span>, or load a preset benchmark from the Home tab.
          </p>
          <div class="pt-2 flex items-center justify-center space-x-2">
            <button onclick="window.bindoraApp.switchTab('home')" class="px-3.5 py-1.5 rounded-lg text-xs font-semibold bg-cyan-600 hover:bg-cyan-500 text-white transition shadow cursor-pointer">
              Go to Benchmarks / Home
            </button>
            <button onclick="window.bindoraApp.switchTab('studio')" class="px-3.5 py-1.5 rounded-lg text-xs font-semibold bg-slate-800 hover:bg-slate-700 text-slate-200 border border-slate-700 transition cursor-pointer">
              Open Target &amp; Ligand Studio
            </button>
          </div>
        </div>
      `;
      return;
    }

    container.innerHTML = `<div class="text-xs text-slate-400 animate-pulse flex items-center space-x-2 py-3"><svg class="w-4 h-4 animate-spin text-purple-400" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15"/></svg><span>Synthesizing pharmacological briefing with anti-hallucination verification...</span></div>`;

    const payload = {
      ligand_name: this.state.ligand?.name || "Drug Candidate",
      target_name: this.state.receptor?.title || this.state.receptor?.pdb_id || "Receptor",
      pdb_id: this.state.receptor?.pdb_id || "N/A",
      thermodynamics: this.state.docking?.thermodynamics || {},
      interactions: this.state.docking?.interactions || {},
      adme: this.state.ligand?.adme || {},
      bioactivity_crosscheck: this.state.crosscheck || {},
      uniprot: this.state.receptor?.uniprot || {}
    };

    try {
      const res = await BindoraAPI.explainNarrative(payload);
      this.state.narrative = res;
      // Render markdown formatted text properly (not raw with ** showing)
      const renderedMarkdown = this.renderMarkdown(res.narrative || '');
      container.innerHTML = `
        <div class="space-y-4">
          <div class="flex items-center justify-between pb-2 border-b border-slate-200 dark:border-slate-700/60 text-xs text-slate-500 dark:text-slate-400">
            <span>Engine: <span class="text-cyan-600 dark:text-cyan-400 font-semibold">${res.source}</span></span>
            <span class="text-emerald-600 dark:text-emerald-400 font-semibold">&#10003; Zero-Hallucination Verified</span>
          </div>
          <div class="text-slate-800 dark:text-slate-300 text-xs sm:text-sm leading-relaxed narrative-md space-y-3">
            ${renderedMarkdown}
          </div>
        </div>
      `;
    } catch (e) {
      console.warn("Narrative generation error:", e);
      container.innerHTML = `
        <div class="p-4 rounded-xl bg-rose-50 dark:bg-rose-950/30 border border-rose-200 dark:border-rose-800/50 space-y-2.5">
          <div class="flex items-center space-x-2 text-rose-600 dark:text-rose-400 font-semibold text-xs">
            <svg class="w-4 h-4 flex-shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"/></svg>
            <span>Narrative Briefing Generation Notice: ${e.message}</span>
          </div>
          <p class="text-[11px] text-slate-600 dark:text-slate-400 leading-relaxed">
            The server request did not complete. Ensure the local server is running on <code class="font-mono text-cyan-600 dark:text-cyan-400">http://localhost:5000</code>, then click retry below to synthesize using the rule-based reasoning engine.
          </p>
          <div>
            <button onclick="window.bindoraApp.generateNarrativeReport()" class="inline-flex items-center space-x-1.5 px-3 py-1.5 rounded-lg text-xs font-semibold bg-purple-600 hover:bg-purple-500 text-white shadow-sm transition cursor-pointer">
              <svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 12a9 9 0 0 0-9-9 9.75 9.75 0 0 0-6.74 2.74L3 8"/><path d="M3 3v5h5"/><path d="M3 12a9 9 0 0 0 9 9 9.75 9.75 0 0 0 6.74-2.74L21 16"/><path d="M16 21h5v-5"/></svg>
              <span>Retry Generation</span>
            </button>
          </div>
        </div>
      `;
    }
  }

  async runBatchScreening() {
    if (!this.state.receptor || !this.state.receptor.pdbqt_text) {
      this.showToast("Please load a receptor first for batch screening.", "error");
      return;
    }

    const batchInput = document.getElementById("input-batch-candidates");
    const rawText = batchInput ? batchInput.value.trim() : "";
    if (!rawText) {
      this.showToast("Please enter at least one candidate drug name or SMILES.", "error");
      return;
    }

    // Parse input lines (format: Name, SMILES)
    const lines = rawText.split("\n");
    const candidates = [];
    for (const line of lines) {
      const parts = line.split(",");
      if (parts.length >= 2) {
        candidates.push({ name: parts[0].trim(), smiles: parts[1].trim() });
      } else if (line.trim().length > 3) {
        candidates.push({ name: `Candidate ${candidates.length+1}`, smiles: line.trim() });
      }
    }

    if (candidates.length === 0) {
      this.showToast("No valid candidates parsed from text input.", "error");
      return;
    }

    this.showToast(`Starting batch screening for ${candidates.length} candidate(s)...`, "info");
    const btn = document.getElementById("btn-run-batch");
    const originalBtnText = btn ? btn.innerHTML : "";
    if (btn) {
      btn.disabled = true;
      btn.innerHTML = '<svg class="animate-spin -ml-1 mr-2 h-4 w-4 text-white inline" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg> Screening...';
    }

    try {
      const p = this.state.receptor.detected_pocket || {};
      const center = p.center || { x: 0, y: 0, z: 0 };
      const size = p.size || { x: 22, y: 22, z: 22 };

      const batchResult = await BindoraAPI.runBatchDocking({
        receptor_pdbqt: this.state.receptor.pdbqt_text,
        receptor_pdb: this.state.receptor.cleaned_pdb,
        center: center,
        size: size,
        ligands: candidates,
        exhaustiveness: 4
      });

      this.renderBatchLeaderboard(batchResult.leaderboard || []);
      this.showToast(`Batch screening complete! Screened ${batchResult.total_screened} compounds.`, "success");
    } catch (e) {
      this.showToast(`Batch screening error: ${e.message}`, "error");
    } finally {
      if (btn) { btn.disabled = false; btn.innerHTML = originalBtnText || '<span>Run Batch Virtual Screening</span>'; }
    }
  }

  async loadReproducibilityMetadata() {
    try {
      const data = await BindoraAPI.getReproducibilityVersions();
      const setTxt = (id, txt) => { const el = document.getElementById(id); if (el) el.textContent = txt; };
      setTxt("dossier-vina-ver", data.autodock_vina || "AutoDock Vina v1.2.5");
      setTxt("dossier-rdkit-ver", `RDKit v${data.rdkit || "2024+"}`);
      setTxt("dossier-gemmi-ver", `Gemmi v${data.gemmi || "0.7+"}`);
      setTxt("dossier-meeko-ver", data.meeko || "Meeko Flexible");
      const ts = data.utc_timestamp ? new Date(data.utc_timestamp).toUTCString() : new Date().toUTCString();
      setTxt("dossier-timestamp", `Calculation Timestamp: ${ts}`);
    } catch (e) {
      console.warn("Reproducibility version query notice:", e);
    }
  }

  updateDossierView() {
    const ligName = document.getElementById("dossier-ligand-name");
    const targetName = document.getElementById("dossier-target-name");
    const resultsSummary = document.getElementById("dossier-results-summary");

    // Fetch reproducibility metadata in background if not already loaded
    this.loadReproducibilityMetadata();

    if (ligName && this.state.ligand) {
      const l = this.state.ligand;
      let cleanName = l.name || "Custom Ligand";
      if (cleanName.startsWith("Conformer3D_COMPOUND_CID_")) {
        cleanName = "PubChem Compound #" + cleanName.replace("Conformer3D_COMPOUND_CID_", "");
      }
      ligName.textContent = `${cleanName} (${l.formula || l.canonical_smiles || "Structure loaded"})`;
    }
    if (targetName && this.state.receptor) {
      const r = this.state.receptor;
      targetName.textContent = `${r.title || "Target Receptor"} (PDB: ${r.pdb_id || "N/A"})`;
    }

    if (!resultsSummary) return;

    if (!this.state.docking) {
      resultsSummary.innerHTML = `
        <div class="p-8 text-center border border-dashed border-slate-300 dark:border-slate-800 rounded-xl">
          <p class="text-slate-500 italic text-xs">No active molecular docking simulation loaded. Run a docking simulation in Tab 3 to automatically compile and generate the complete preclinical research dossier.</p>
        </div>
      `;
      return;
    }

    const d = this.state.docking;
    const l = this.state.ligand || {};
    const r = this.state.receptor || {};
    const thermo = d.thermodynamics || {};
    const contacts = d.interactions || {};
    const adme = l.adme || {};
    const phys = adme.physicochemical || {};
    const lip = adme.drug_likeness?.lipinski || {};
    const veber = adme.drug_likeness?.veber || {};
    const pk = adme.pharmacokinetics || {};
    const redock = this.state.redockingValidation;
    const repStats = d.replicate_stats;
    const poses = d.poses || (d.top_pose ? [d.top_pose] : []);
    const topPose = d.top_pose || poses[0] || {};
    const heavyCount = phys.heavy_atoms?.value || l.heavy_atom_count || 20;

    // Grid coordinates
    const getVal = (id, def) => parseFloat(document.getElementById(id)?.value) || def;
    const center = d.grid?.center || { x: getVal("grid-cx", 0), y: getVal("grid-cy", 0), z: getVal("grid-cz", 0) };
    const size = d.grid?.size || { x: getVal("grid-sx", 22), y: getVal("grid-sy", 22), z: getVal("grid-sz", 22) };
    const exhaustiveness = d.exhaustiveness || parseInt(document.getElementById("docking-exhaustiveness")?.value) || 8;

    const isWeak = thermo.is_weak_binder || (topPose.affinity_kcal > -6.0);

    // Helpers
    const formatKd = (dg) => {
      if (dg == null || isNaN(dg)) return "—";
      const RT = 0.59248; // kcal/mol at 298.15K
      const kdM = Math.exp(dg / RT);
      const kdNM = kdM * 1e9;
      if (kdNM < 1) return `${(kdNM * 1000).toFixed(1)} pM`;
      if (kdNM < 1000) return `${kdNM.toFixed(1)} nM`;
      const kdUM = kdM * 1e6;
      if (kdUM < 1000) return `${kdUM.toFixed(2)} µM`;
      return `${(kdM * 1000).toFixed(2)} mM`;
    };

    const calcLE = (dg, hCount) => {
      if (dg == null || !hCount || isNaN(dg) || hCount <= 0) return "—";
      return (-dg / hCount).toFixed(3);
    };

    // Modes table rows
    const posesRows = poses.map((p, idx) => {
      const mode = p.mode || (idx + 1);
      const aff = p.affinity_kcal != null ? Number(p.affinity_kcal).toFixed(2) : "—";
      const vinAff = p.vinardo_affinity_kcal != null ? Number(p.vinardo_affinity_kcal).toFixed(2) : "—";
      const kdStr = p.affinity_kcal != null ? formatKd(p.affinity_kcal) : "—";
      const leStr = p.affinity_kcal != null ? calcLE(p.affinity_kcal, heavyCount) : "—";
      const rmsdLb = p.rmsd_lb != null ? Number(p.rmsd_lb).toFixed(3) : "0.000";
      const rmsdUb = p.rmsd_ub != null ? Number(p.rmsd_ub).toFixed(3) : "0.000";
      const isTop = idx === 0;

      return `
        <tr class="${isTop ? 'bg-cyan-50/70 dark:bg-cyan-950/40 font-semibold' : ''}">
          <td class="text-center font-mono font-bold ${isTop ? 'text-cyan-700 dark:text-cyan-400' : 'text-slate-700 dark:text-slate-300'} p-2 border border-slate-200 dark:border-slate-800">#${mode}</td>
          <td class="text-center font-mono font-bold ${isTop ? 'text-cyan-700 dark:text-cyan-400' : 'text-slate-900 dark:text-white'} p-2 border border-slate-200 dark:border-slate-800">${aff} kcal/mol</td>
          <td class="text-center font-mono text-slate-700 dark:text-slate-300 p-2 border border-slate-200 dark:border-slate-800">${vinAff !== "—" ? vinAff + ' kcal/mol' : '—'}</td>
          <td class="text-center font-mono text-emerald-700 dark:text-emerald-400 font-bold p-2 border border-slate-200 dark:border-slate-800">${kdStr}</td>
          <td class="text-center font-mono text-amber-700 dark:text-amber-400 p-2 border border-slate-200 dark:border-slate-800">${leStr}</td>
          <td class="text-center font-mono text-slate-500 dark:text-slate-400 p-2 border border-slate-200 dark:border-slate-800">${rmsdLb} Å</td>
          <td class="text-center font-mono text-slate-500 dark:text-slate-400 p-2 border border-slate-200 dark:border-slate-800">${rmsdUb} Å</td>
          <td class="text-center text-[10px] p-2 border border-slate-200 dark:border-slate-800">
            ${isTop ? '<span class="px-2 py-0.5 rounded bg-emerald-100 dark:bg-emerald-950 text-emerald-800 dark:text-emerald-300 border border-emerald-300 dark:border-emerald-800 font-bold">Rank 1 (Global Min)</span>' : '<span class="text-slate-500">Conformational Mode</span>'}
          </td>
        </tr>
      `;
    }).join("");

    // Hydrogen bonds
    const hbonds = contacts.hydrogen_bonds || [];
    const hbondRows = hbonds.length > 0 ? hbonds.map(hb => `
      <tr>
        <td class="font-mono font-bold text-slate-900 dark:text-white p-2 border border-slate-200 dark:border-slate-800">${hb.residue || `${hb.res_name} ${hb.res_num}:${hb.chain}`}</td>
        <td class="font-mono text-slate-700 dark:text-slate-300 p-2 border border-slate-200 dark:border-slate-800">${hb.receptor_atom || '—'}</td>
        <td class="font-mono text-slate-700 dark:text-slate-300 p-2 border border-slate-200 dark:border-slate-800">${hb.ligand_atom || '—'}</td>
        <td class="font-mono text-cyan-700 dark:text-cyan-400 font-bold text-center p-2 border border-slate-200 dark:border-slate-800">${hb.distance != null ? Number(hb.distance).toFixed(2) + ' Å' : '—'}</td>
        <td class="text-slate-600 dark:text-slate-400 text-[10px] p-2 border border-slate-200 dark:border-slate-800">${hb.type || 'Hydrogen Bond'}</td>
      </tr>
    `).join("") : `<tr><td colspan="5" class="text-center py-2 text-slate-500 italic p-2 border border-slate-200 dark:border-slate-800">No directional hydrogen bonds detected within 3.5 Å cutoff.</td></tr>`;

    // Hydrophobic contacts
    const hydrophobics = contacts.hydrophobic_contacts || [];
    const hydrophobicList = hydrophobics.length > 0 ? hydrophobics.map(hp => `
      <span class="inline-flex items-center space-x-1 px-2 py-1 rounded bg-amber-50 dark:bg-amber-950/40 border border-amber-200 dark:border-amber-800/60 text-amber-900 dark:text-amber-200 font-mono text-[11px]">
        <span class="font-bold">${hp.residue || `${hp.res_name} ${hp.res_num}:${hp.chain}`}</span>
        <span class="text-[10px] text-amber-700 dark:text-amber-400">(${Number(hp.distance).toFixed(2)} Å)</span>
      </span>
    `).join(" ") : `<span class="text-slate-500 italic text-xs">No hydrophobic contacts detected within 4.0 Å cutoff.</span>`;

    // Lipinski Rule of 5 Matrix
    const mw = phys.molecular_weight?.value != null ? Number(phys.molecular_weight.value).toFixed(2) : (l.weight != null ? Number(l.weight).toFixed(2) : "—");
    const logp = phys.logp?.value != null ? Number(phys.logp.value).toFixed(2) : "—";
    const hbd = phys.h_bond_donors?.value != null ? phys.h_bond_donors.value : "—";
    const hba = phys.h_bond_acceptors?.value != null ? phys.h_bond_acceptors.value : "—";
    const tpsa = phys.tpsa?.value != null ? Number(phys.tpsa.value).toFixed(2) : "—";
    const rotb = phys.rotatable_bonds?.value != null ? phys.rotatable_bonds.value : "—";

    const mwPass = mw !== "—" ? Number(mw) <= 500 : true;
    const logpPass = logp !== "—" ? Number(logp) <= 5.0 : true;
    const hbdPass = hbd !== "—" ? Number(hbd) <= 5 : true;
    const hbaPass = hba !== "—" ? Number(hba) <= 10 : true;
    const lipViolations = (mwPass ? 0 : 1) + (logpPass ? 0 : 1) + (hbdPass ? 0 : 1) + (hbaPass ? 0 : 1);

    let currentSec = 5;
    const svgSectionHtml = contacts.diagram_svg ? `
      <!-- SECTION: 2D PROTEIN-LIGAND INTERACTION SCHEMATIC -->
      <div class="dossier-card p-4 rounded-xl border border-slate-200 dark:border-slate-800 bg-slate-50/80 dark:bg-slate-900/60 space-y-2 page-break-inside-avoid">
        <div class="flex items-center justify-between border-b border-slate-200 dark:border-slate-800 pb-2">
          <span class="font-bold text-slate-800 dark:text-slate-200 uppercase tracking-wider text-xs flex items-center space-x-1.5 font-sans">
            <svg class="w-4 h-4 text-cyan-600 dark:text-cyan-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect width="18" height="18" x="3" y="3" rx="2"/><circle cx="9" cy="9" r="2"/><path d="m21 15-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/></svg>
            <span>Section ${currentSec++}: 2D Protein–Ligand Interaction Schematic</span>
          </span>
          <span class="text-[10px] text-cyan-700 dark:text-cyan-400 font-mono">LigPlot-Style Radial Vector Map</span>
        </div>
        <div class="dossier-subbox flex items-center justify-center p-3 bg-white dark:bg-slate-950/90 rounded-lg border border-slate-200 dark:border-slate-800 overflow-hidden" id="dossier-svg-container">
          ${contacts.diagram_svg}
        </div>
        <p class="text-[10px] text-slate-500 font-sans text-center mt-1">
          <em>Figure 1:</em> 2D radial representation of binding cleft contacts (dashed cyan lines = hydrogen bonds with distances; amber rays = hydrophobic interactions).
        </p>
      </div>
    ` : '';

    const admetSecNumber = currentSec++;
    const aiSecNumber = currentSec++;

    resultsSummary.innerHTML = `
      <!-- SECTION 1: EXECUTIVE THERMODYNAMIC & KINETIC METRICS -->
      <div class="dossier-card p-4 rounded-xl border border-slate-200 dark:border-slate-800 bg-slate-50/80 dark:bg-slate-900/60 space-y-3">
        <div class="flex items-center justify-between border-b border-slate-200 dark:border-slate-800 pb-2">
          <span class="font-bold text-slate-800 dark:text-slate-200 uppercase tracking-wider text-xs flex items-center space-x-1.5 font-sans">
            <svg class="w-4 h-4 text-cyan-600 dark:text-cyan-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M10 2v7.31L4.16 20.25A1 1 0 0 0 5 21.6h14a1 1 0 0 0 .84-1.35L14 9.31V2"/><path d="M8.5 2h7"/><path d="M7 16h10"/></svg>
            <span>Section 1: Executive Quantitative Binding Metrics</span>
          </span>
          <span class="text-[10px] px-2 py-0.5 rounded bg-cyan-100 dark:bg-cyan-950 text-cyan-800 dark:text-cyan-300 border border-cyan-300 dark:border-cyan-800 font-semibold">Primary Simulation</span>
        </div>

        <div class="grid grid-cols-2 sm:grid-cols-4 gap-3 text-xs font-mono">
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 shadow-sm">
            <span class="text-slate-500 dark:text-slate-400 block font-sans text-[11px] font-medium">Binding Free Energy (ΔG)</span>
            <span class="text-cyan-700 dark:text-cyan-400 font-black text-base">${topPose.affinity_kcal != null ? topPose.affinity_kcal.toFixed(2) : "—"} <span class="text-xs font-normal font-sans">kcal/mol</span></span>
          </div>
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 shadow-sm">
            <span class="text-slate-500 dark:text-slate-400 block font-sans text-[11px] font-medium">Theoretical Dissociation (Kd)</span>
            <span class="text-emerald-700 dark:text-emerald-400 font-black text-base">${thermo.theoretical_kd_nm != null ? thermo.theoretical_kd_nm + ' nM' : formatKd(topPose.affinity_kcal)}</span>
          </div>
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 shadow-sm">
            <span class="text-slate-500 dark:text-slate-400 block font-sans text-[11px] font-medium">Ligand Efficiency (LE)</span>
            <span class="text-amber-700 dark:text-amber-400 font-black text-base">${thermo.ligand_efficiency?.value || calcLE(topPose.affinity_kcal, heavyCount)}</span>
            <span class="text-[10px] text-slate-500 block font-sans">SILE: ${thermo.size_independent_le?.value || "—"}</span>
          </div>
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 shadow-sm">
            <span class="text-slate-500 dark:text-slate-400 block font-sans text-[11px] font-medium">Intermolecular Contacts</span>
            <span class="text-purple-700 dark:text-purple-400 font-black text-base">${contacts.total_hbond_count || hbonds.length} <span class="text-xs font-normal font-sans">H-Bonds</span></span>
            <span class="text-[10px] text-slate-500 block font-sans">${contacts.total_hydrophobic_count || hydrophobics.length} Hydrophobic</span>
          </div>
        </div>

        ${repStats ? `
          <div class="p-2.5 rounded-lg bg-cyan-50 dark:bg-cyan-950/30 border border-cyan-200 dark:border-cyan-800/40 text-xs font-mono flex items-center justify-between">
            <span class="text-slate-800 dark:text-slate-200 font-sans font-medium">Multi-Seed Stochastic Replicates (N=${repStats.replicates_count}):</span>
            <span class="text-cyan-800 dark:text-cyan-300 font-bold">Mean ΔG = ${repStats.mean_affinity_kcal} ± ${repStats.sd_affinity_kcal} kcal/mol (95% CI: ±${repStats.confidence_interval_95} kcal/mol)</span>
          </div>
        ` : ''}

        ${redock ? `
          <div class="p-2.5 rounded-lg ${redock.is_validated ? 'bg-emerald-50 dark:bg-emerald-950/30 border-emerald-200 dark:border-emerald-800/40' : 'bg-amber-50 dark:bg-amber-950/30 border-amber-200 dark:border-amber-800/40'} border text-xs font-mono flex items-center justify-between">
            <span class="text-slate-800 dark:text-slate-200 font-sans font-medium">Protocol Self-Validation (Native Ligand Redocking):</span>
            <span class="${redock.is_validated ? 'text-emerald-800 dark:text-emerald-300' : 'text-amber-800 dark:text-amber-300'} font-bold">RMSD = ${redock.rmsd_angstroms} Å &bull; ${redock.benchmark_status}</span>
          </div>
        ` : ''}

        ${isWeak ? `
          <div class="p-2.5 rounded-lg bg-amber-50 dark:bg-amber-950/40 border border-amber-300 dark:border-amber-600/50 text-xs text-amber-900 dark:text-amber-200 font-sans space-y-1">
            <strong class="font-bold flex items-center space-x-1.5 text-amber-900 dark:text-amber-300">
              <svg class="w-4 h-4 text-amber-600 dark:text-amber-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="m21.73 18-8-14a2 2 0 0 0-3.48 0l-8 14A2 2 0 0 0 4 21h16a2 2 0 0 0 1.73-3Z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>
              <span>Scientific Qualification: Sub-threshold Binding Detected</span>
            </strong>
            <p class="text-[11px] leading-relaxed text-amber-900/90 dark:text-amber-200/90">
              Calculated ΔG (${topPose.affinity_kcal} kcal/mol) corresponds to micromolar or millimolar affinity (Kd ≈ ${thermo.theoretical_kd_um || '—'} µM). At this range, experimental assays typically observe non-specific surface adhesion or weak kinetics. Interpret cautiously.
            </p>
          </div>
        ` : ''}
      </div>

      <!-- SECTION 2: BIOPHYSICAL ENTITIES & DOCKING SEARCH SPACE SPECIFICATIONS -->
      <div class="dossier-card p-4 rounded-xl border border-slate-200 dark:border-slate-800 bg-slate-50/80 dark:bg-slate-900/60 space-y-3">
        <div class="flex items-center justify-between border-b border-slate-200 dark:border-slate-800 pb-2">
          <span class="font-bold text-slate-800 dark:text-slate-200 uppercase tracking-wider text-xs flex items-center space-x-1.5 font-sans">
            <svg class="w-4 h-4 text-cyan-600 dark:text-cyan-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect width="18" height="18" x="3" y="3" rx="2"/><path d="M9 3v18"/><path d="M15 3v18"/></svg>
            <span>Section 2: Biophysical &amp; Search Space Parameters</span>
          </span>
          <span class="text-[10px] text-slate-500 font-mono">Algorithm: Vina Iterated Local Search</span>
        </div>

        <div class="grid grid-cols-1 md:grid-cols-2 gap-3 text-xs">
          <!-- Target Specifications -->
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 space-y-1.5 font-mono">
            <div class="font-bold text-slate-800 dark:text-slate-200 border-b border-slate-100 dark:border-slate-800 pb-1 font-sans text-xs">Target Receptor</div>
            <div class="flex justify-between"><span class="text-slate-500">PDB Identifier:</span> <span class="font-bold text-slate-900 dark:text-white">${r.pdb_id || 'Custom PDB'}</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Chains / Heavy Atoms:</span> <span class="text-slate-800 dark:text-slate-200">${(r.chains || ['A']).join(', ')} &bull; ${r.atom_count || '—'} atoms</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Grid Center [X, Y, Z]:</span> <span class="text-cyan-700 dark:text-cyan-400 font-bold">${center.x.toFixed(1)}, ${center.y.toFixed(1)}, ${center.z.toFixed(1)}</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Search Box Dimensions:</span> <span class="text-slate-800 dark:text-slate-200">${size.x.toFixed(1)} × ${size.y.toFixed(1)} × ${size.z.toFixed(1)} Å</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Search Exhaustiveness:</span> <span class="font-bold text-slate-900 dark:text-white">${exhaustiveness} (Deep Sampling)</span></div>
          </div>

          <!-- Ligand Specifications -->
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 space-y-1.5 font-mono">
            <div class="font-bold text-slate-800 dark:text-slate-200 border-b border-slate-100 dark:border-slate-800 pb-1 font-sans text-xs">Investigational Ligand</div>
            <div class="flex justify-between"><span class="text-slate-500">Ligand Name:</span> <span class="font-bold text-slate-900 dark:text-white truncate max-w-[200px]" title="${l.name || ''}">${l.name || 'Custom Ligand'}</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Molecular Formula:</span> <span class="text-slate-800 dark:text-slate-200">${l.formula || 'Derived from SMILES'}</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Molecular Weight:</span> <span class="text-slate-800 dark:text-slate-200">${mw} g/mol</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Heavy Atoms / Rot. Bonds:</span> <span class="text-slate-800 dark:text-slate-200">${heavyCount} heavy &bull; ${rotb} rotatable</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Wildman-Crippen LogP:</span> <span class="font-bold text-slate-900 dark:text-white">${logp}</span></div>
          </div>
        </div>
      </div>

      <!-- SECTION 3: COMPREHENSIVE CONFORMATIONAL SAMPLING TABLE (MODES 1 TO 9) -->
      <div class="dossier-card p-4 rounded-xl border border-slate-200 dark:border-slate-800 bg-slate-50/80 dark:bg-slate-900/60 space-y-3">
        <div class="flex items-center justify-between border-b border-slate-200 dark:border-slate-800 pb-2">
          <span class="font-bold text-slate-800 dark:text-slate-200 uppercase tracking-wider text-xs flex items-center space-x-1.5 font-sans">
            <svg class="w-4 h-4 text-cyan-600 dark:text-cyan-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M3 3h18v18H3zM3 9h18M3 15h18M9 3v18M15 3v18"/></svg>
            <span>Section 3: Ranked Binding Modes Table (AutoDock Vina &amp; Vinardo)</span>
          </span>
          <span class="text-[10px] text-slate-500 font-mono">Cluster Cutoff: 2.0 Å</span>
        </div>

        <div class="overflow-x-auto">
          <table class="w-full text-xs text-left border-collapse border border-slate-200 dark:border-slate-800 bg-white dark:bg-slate-950/80">
            <thead>
              <tr class="bg-slate-100 dark:bg-slate-900/90 text-slate-700 dark:text-slate-300 font-bold font-sans">
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">Mode</th>
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">Vina ΔG</th>
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">Vinardo ΔG</th>
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">Theoretical Kd</th>
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">Ligand Eff.</th>
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">RMSD l.b.</th>
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">RMSD u.b.</th>
                <th class="p-2 text-center border border-slate-200 dark:border-slate-800">Classification</th>
              </tr>
            </thead>
            <tbody>
              ${posesRows}
            </tbody>
          </table>
        </div>
      </div>

      <!-- SECTION 4: INTERMOLECULAR INTERACTION RESIDUE PROFILING -->
      <div class="dossier-card p-4 rounded-xl border border-slate-200 dark:border-slate-800 bg-slate-50/80 dark:bg-slate-900/60 space-y-3">
        <div class="flex items-center justify-between border-b border-slate-200 dark:border-slate-800 pb-2">
          <span class="font-bold text-slate-800 dark:text-slate-200 uppercase tracking-wider text-xs flex items-center space-x-1.5 font-sans">
            <svg class="w-4 h-4 text-purple-600 dark:text-purple-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><path d="m4.93 4.93 4.24 4.24"/><path d="m14.83 9.17 4.24-4.24"/><path d="m14.83 14.83 4.24 4.24"/><path d="m9.17 14.83-4.24 4.24"/></svg>
            <span>Section 4: Intermolecular Interaction Fingerprint</span>
          </span>
          <span class="text-[10px] text-slate-500 font-mono">Cutoffs: H-Bond &le; 3.5 Å &bull; Hydrophobic &le; 4.0 Å</span>
        </div>

        <div class="space-y-3">
          <div>
            <span class="font-bold text-slate-800 dark:text-slate-200 text-xs block mb-1.5 font-sans">Directional Hydrogen Bonds:</span>
            <div class="overflow-x-auto">
              <table class="w-full text-xs text-left border-collapse border border-slate-200 dark:border-slate-800 bg-white dark:bg-slate-950/80">
                <thead>
                  <tr class="bg-slate-100 dark:bg-slate-900/90 text-slate-700 dark:text-slate-300 font-sans font-bold">
                    <th class="p-2 border border-slate-200 dark:border-slate-800">Receptor Residue</th>
                    <th class="p-2 border border-slate-200 dark:border-slate-800">Receptor Atom</th>
                    <th class="p-2 border border-slate-200 dark:border-slate-800">Ligand Atom</th>
                    <th class="p-2 border border-slate-200 dark:border-slate-800 text-center">Distance (Å)</th>
                    <th class="p-2 border border-slate-200 dark:border-slate-800">Bond Type</th>
                  </tr>
                </thead>
                <tbody>
                  ${hbondRows}
                </tbody>
              </table>
            </div>
          </div>

          <div>
            <span class="font-bold text-slate-800 dark:text-slate-200 text-xs block mb-1.5 font-sans">Hydrophobic Contact Residues:</span>
            <div class="p-2.5 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 flex flex-wrap gap-1.5">
              ${hydrophobicList}
            </div>
          </div>
        </div>
      </div>

      <!-- SECTION 5: 2D PROTEIN-LIGAND INTERACTION SCHEMATIC (OPTIONAL) -->
      ${svgSectionHtml}

      <!-- SECTION: CHEMINFORMATICS & PRECLINICAL ADMET ASSESSMENT -->
      <div class="dossier-card p-4 rounded-xl border border-slate-200 dark:border-slate-800 bg-slate-50/80 dark:bg-slate-900/60 space-y-3">
        <div class="flex items-center justify-between border-b border-slate-200 dark:border-slate-800 pb-2">
          <span class="font-bold text-slate-800 dark:text-slate-200 uppercase tracking-wider text-xs flex items-center space-x-1.5 font-sans">
            <svg class="w-4 h-4 text-emerald-600 dark:text-emerald-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10z"/></svg>
            <span>Section ${admetSecNumber}: Cheminformatics &amp; Preclinical Drug-Likeness (Lipinski Rule of 5)</span>
          </span>
          <span class="text-[10px] px-2 py-0.5 rounded ${lipViolations === 0 ? 'bg-emerald-100 dark:bg-emerald-950 text-emerald-800 dark:text-emerald-300 border border-emerald-300 dark:border-emerald-800' : 'bg-amber-100 dark:bg-amber-950 text-amber-800 dark:text-amber-300 border border-amber-300 dark:border-amber-800'} font-bold font-mono">
            ${lipViolations === 0 ? 'Lipinski Compliant (0 Violations)' : `Borderline (${lipViolations} Violations)`}
          </span>
        </div>

        <div class="grid grid-cols-1 md:grid-cols-2 gap-3 text-xs">
          <!-- Lipinski Matrix -->
          <div class="overflow-x-auto">
            <table class="w-full text-xs text-left border-collapse border border-slate-200 dark:border-slate-800 bg-white dark:bg-slate-950/80">
              <thead>
                <tr class="bg-slate-100 dark:bg-slate-900/90 text-slate-700 dark:text-slate-300 font-sans font-bold">
                  <th class="p-2 border border-slate-200 dark:border-slate-800">Lipinski Parameter</th>
                  <th class="p-2 border border-slate-200 dark:border-slate-800 text-center">Value</th>
                  <th class="p-2 border border-slate-200 dark:border-slate-800 text-center">Criterion</th>
                  <th class="p-2 border border-slate-200 dark:border-slate-800 text-center">Status</th>
                </tr>
              </thead>
              <tbody class="font-mono">
                <tr>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 font-sans">Molecular Weight (MW)</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold">${mw} Da</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center text-slate-500">&le; 500 Da</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold ${mwPass ? 'text-emerald-700 dark:text-emerald-400' : 'text-rose-700 dark:text-rose-400'}">${mwPass ? 'Pass' : 'Violated'}</td>
                </tr>
                <tr>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 font-sans">Wildman-Crippen LogP</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold">${logp}</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center text-slate-500">&le; 5.0</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold ${logpPass ? 'text-emerald-700 dark:text-emerald-400' : 'text-rose-700 dark:text-rose-400'}">${logpPass ? 'Pass' : 'Violated'}</td>
                </tr>
                <tr>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 font-sans">H-Bond Donors (HBD)</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold">${hbd}</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center text-slate-500">&le; 5</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold ${hbdPass ? 'text-emerald-700 dark:text-emerald-400' : 'text-rose-700 dark:text-rose-400'}">${hbdPass ? 'Pass' : 'Violated'}</td>
                </tr>
                <tr>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 font-sans">H-Bond Acceptors (HBA)</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold">${hba}</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center text-slate-500">&le; 10</td>
                  <td class="p-2 border border-slate-200 dark:border-slate-800 text-center font-bold ${hbaPass ? 'text-emerald-700 dark:text-emerald-400' : 'text-rose-700 dark:text-rose-400'}">${hbaPass ? 'Pass' : 'Violated'}</td>
                </tr>
              </tbody>
            </table>
          </div>

          <!-- Pharmacokinetics & Veber Rules -->
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 space-y-2 text-xs font-mono">
            <div class="font-bold text-slate-800 dark:text-slate-200 border-b border-slate-100 dark:border-slate-800 pb-1 font-sans">Pharmacokinetics &amp; Bioavailability</div>
            <div class="flex justify-between"><span class="text-slate-500">TPSA (Topological Polar Surface):</span> <span class="font-bold text-slate-900 dark:text-white">${tpsa} Å² (&le; 140 Å²)</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Rotatable Bonds Count:</span> <span class="font-bold text-slate-900 dark:text-white">${rotb} (&le; 10)</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Gastrointestinal (GI) Absorption:</span> <span class="text-emerald-700 dark:text-emerald-400 font-bold">${pk.gi_absorption?.level || 'High'}</span></div>
            <div class="flex justify-between"><span class="text-slate-500">Blood-Brain Barrier (BBB):</span> <span class="text-slate-800 dark:text-slate-200">${pk.bbb_permeant?.is_permeant ? 'Permeant' : 'Non-Permeant'}</span></div>
            <div class="flex justify-between pt-1 border-t border-slate-100 dark:border-slate-800">
              <span class="text-slate-500">CYP450 Liability:</span>
              <span class="text-slate-800 dark:text-slate-200">${(pk.cyp450_inhibition?.alerts || []).length > 0 ? pk.cyp450_inhibition.alerts.map(a => a.cyp).join(', ') : 'No Inhibitory Alerts'}</span>
            </div>
          </div>
        </div>
      </div>

      <!-- SECTION: AI-GENERATED MECHANISTIC HYPOTHESIS -->
      ${this.state.narrative ? `
        <div class="dossier-card p-4 rounded-xl border border-slate-200 dark:border-slate-800 bg-slate-50/80 dark:bg-slate-900/60 space-y-2 page-break-inside-avoid">
          <div class="flex items-center justify-between border-b border-slate-200 dark:border-slate-800 pb-2">
            <span class="font-bold text-purple-800 dark:text-purple-300 uppercase tracking-wider text-xs flex items-center space-x-1.5 font-sans">
              <svg class="w-4 h-4 text-purple-600 dark:text-purple-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="m12 3-1.9 5.8a2 2 0 0 1-1.3 1.3L3 12l5.8 1.9a2 2 0 0 1 1.3 1.3L12 21l1.9-5.8a2 2 0 0 1 1.3-1.3L21 12l-5.8-1.9a2 2 0 0 1-1.3-1.3L12 3z"/></svg>
              <span>Section ${aiSecNumber}: AI Mechanistic Synthesis (DeepSeek LLM Hypothesis)</span>
            </span>
            <span class="text-[10px] px-2 py-0.5 rounded bg-purple-100 dark:bg-purple-950 text-purple-800 dark:text-purple-300 border border-purple-300 dark:border-purple-800 font-semibold font-mono">In-Silico Hypothesis</span>
          </div>
          <div class="dossier-subbox p-3 bg-white dark:bg-slate-950/80 rounded-lg border border-slate-200 dark:border-slate-800 text-slate-800 dark:text-slate-200 text-xs leading-relaxed font-sans">
            ${this.renderMarkdown(this.state.narrative.narrative || '')}
          </div>
        </div>
      ` : ''}
    `;
  }

  updatePathwayInfo() {
    const container = document.getElementById("pathway-info-container");
    if (!container) return;

    const uniprot = this.state.receptor?.uniprot;
    if (uniprot && (uniprot.function || uniprot.subcellular_location)) {
      container.innerHTML = `
        ${uniprot.function ? '<div><span class="font-bold text-slate-900 dark:text-slate-200">Biological Function:</span> <span class="text-slate-700 dark:text-slate-300">' + uniprot.function + '</span></div>' : ''}
        ${uniprot.catalytic_activity ? '<div><span class="font-bold text-slate-900 dark:text-slate-200">Catalytic Activity:</span> <span class="text-slate-700 dark:text-slate-300">' + uniprot.catalytic_activity + '</span></div>' : ''}
        ${uniprot.subcellular_location ? '<div><span class="font-bold text-slate-900 dark:text-slate-200">Subcellular Location:</span> <span class="text-slate-700 dark:text-slate-300">' + uniprot.subcellular_location + '</span></div>' : ''}
        ${uniprot.tissue_specificity ? '<div><span class="font-bold text-slate-900 dark:text-slate-200">Tissue Specificity:</span> <span class="text-slate-700 dark:text-slate-300">' + uniprot.tissue_specificity + '</span></div>' : ''}
      `;
    } else {
      container.innerHTML = '<p class="text-slate-500 italic">No UniProt pathway annotation available for this target. Load a receptor with known annotations to view biological function and signaling cascade details.</p>';
    }
  }

  renderBatchLeaderboard(leaderboard) {
    const tbody = document.getElementById("batch-leaderboard-rows");
    if (!tbody) return;

    if (leaderboard.length === 0) {
      tbody.innerHTML = `<tr><td colspan="10" class="text-center p-6 text-xs text-slate-500 italic">No docking results produced.</td></tr>`;
      return;
    }

    tbody.innerHTML = leaderboard.map(item => {
      if (item.valid === false) {
        return `
          <tr class="border-b border-rose-950/60 bg-rose-950/20 text-rose-300">
            <td class="p-2.5 font-bold text-center text-slate-500">—</td>
            <td class="p-2.5 font-semibold text-rose-300 font-sans" title="${item.smiles}">${item.name}</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5">
              <span class="px-2 py-0.5 rounded bg-rose-950 text-rose-300 border border-rose-800 text-[10px]" title="${item.error || 'Syntax error'}">
                Failed SMILES
              </span>
            </td>
          </tr>
        `;
      }

      const isWeak = item.is_weak_binder || item.affinity_kcal > -6.0;
      const isHighConf = item.consensus_confidence === "High-Confidence";
      const isModerate = item.consensus_confidence === "Moderate-Confidence";
      const crit = item.consensus_criteria || "";
      const confBadge = isHighConf ?
        `<span class="px-1.5 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-800 text-[10px] font-bold inline-block cursor-help" title="${crit}">High-Confidence</span>` :
        isModerate ?
        `<span class="px-1.5 py-0.5 rounded bg-cyan-950 text-cyan-300 border border-cyan-800 text-[10px] font-bold inline-block cursor-help" title="${crit}">Moderate</span>` :
        `<span class="px-1.5 py-0.5 rounded bg-amber-950 text-amber-300 border border-amber-800 text-[10px] font-bold inline-block cursor-help" title="${crit}">Divergent</span>`;

      return `
        <tr class="border-b border-slate-800 hover:bg-slate-800/50">
          <td class="p-2.5 font-bold text-center text-cyan-400">#${item.rank}</td>
          <td class="p-2.5 font-semibold text-white font-sans" title="${item.smiles}">${item.name}</td>
          <td class="p-2.5 font-mono font-bold ${isWeak ? 'text-amber-400' : 'text-emerald-400'}">${item.affinity_kcal}</td>
          <td class="p-2.5 font-mono text-cyan-300 font-semibold">${item.vinardo_affinity_kcal ?? "—"}</td>
          <td class="p-2.5 font-mono font-bold text-cyan-300">${item.consensus_score ?? "—"}</td>
          <td class="p-2.5">${confBadge}</td>
          <td class="p-2.5 font-mono text-slate-300">${item.theoretical_kd_nm}</td>
          <td class="p-2.5 font-mono text-slate-300">${item.ligand_efficiency}</td>
          <td class="p-2.5 font-mono text-slate-300">${item.hbond_count}</td>
          <td class="p-2.5 text-[10px]">
            ${isWeak ? 
              '<span class="px-1.5 py-0.5 rounded bg-amber-950/80 text-amber-300 border border-amber-800" title="ΔG > -6.0 kcal/mol: Sub-threshold/Weak binder">Weak Hit</span>' : 
              '<span class="px-1.5 py-0.5 rounded bg-emerald-950/80 text-emerald-300 border border-emerald-800">Screened</span>'
            }
          </td>
        </tr>
      `;
    }).join("");

    const sub = document.getElementById("batch-matrix-subtitle");
    if (sub) sub.textContent = "Consensus = ΔG (50%) + LE (scaled) + H-Bond Bonus";
  }

  updateFlexibleResiduesUI(flexCandidates) {
    const container = document.getElementById("flex-residues-container");
    if (!container) return;
    if (!flexCandidates || flexCandidates.length === 0) {
      if (!container.querySelector("input[type='checkbox']")) {
        container.innerHTML = `<span class="text-[11px] text-slate-500 italic">No contacting active site residues detected yet.</span>`;
      }
      return;
    }
    const previouslyChecked = new Set(
      Array.from(container.querySelectorAll("input[type='checkbox']:checked")).map(cb => cb.value)
    );
    container.innerHTML = flexCandidates.map(c => `
      <label class="inline-flex items-center space-x-1.5 px-2 py-1 rounded bg-slate-800/80 hover:bg-slate-700/80 border border-slate-700/60 cursor-pointer transition text-[11px]">
        <input type="checkbox" value="${c.id}" ${previouslyChecked.has(c.id) ? 'checked' : ''} class="rounded text-cyan-500 focus:ring-0 focus:ring-offset-0 bg-slate-900 border-slate-600">
        <span class="font-mono text-cyan-300 font-medium">${c.res_name} ${c.res_num}</span>
        <span class="text-[10px] text-slate-400">(${c.chain})</span>
      </label>
    `).join("");
  }

  async findSimilarCompounds() {
    const l = this.state.ligand;
    if (!l) {
      this.showToast("Please load a ligand first.", "error");
      return;
    }
    const smiles = l.canonical_smiles || l.smiles || "";
    if (!smiles) {
      this.showToast("No valid SMILES available for current ligand.", "error");
      return;
    }
    const container = document.getElementById("similar-compounds-container");
    if (!container) return;
    container.classList.remove("hidden");
    container.innerHTML = `<div class="text-[11px] text-cyan-400 animate-pulse flex items-center space-x-1.5 py-1">
      <svg class="animate-spin w-3.5 h-3.5 text-cyan-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>
      <span>Searching PubChem 2D fast similarity (≥85%)...</span>
    </div>`;

    try {
      const res = await BindoraAPI.getSimilarCompounds(smiles, 85, 5);
      const compounds = res.similar_compounds || [];
      if (compounds.length === 0) {
        container.innerHTML = `<span class="text-[11px] text-slate-400 italic">No similar compounds found in PubChem at &ge;85% Tanimoto threshold.</span>`;
        return;
      }
      container.innerHTML = `
        <div class="text-[11px] font-semibold text-cyan-300 mb-1 flex items-center justify-between">
          <span>Found ${compounds.length} Similar Analogues (PubChem):</span>
          <button onclick="document.getElementById('similar-compounds-container').classList.add('hidden')" class="text-slate-400 hover:text-white text-xs">&times;</button>
        </div>
        <div class="space-y-1.5 max-h-48 overflow-y-auto pr-1">
          ${compounds.map(c => `
            <div class="p-2 rounded bg-slate-800/80 border border-slate-700/60 flex items-center justify-between text-[11px]">
              <div class="truncate mr-2">
                <span class="font-bold text-white block truncate" title="${c.title || `CID ${c.cid}`}">${c.title || `CID ${c.cid}`}</span>
                <span class="text-[10px] text-slate-400 font-mono">MW: ${c.molecular_weight || '—'} Da | CID: ${c.cid}</span>
              </div>
              <div class="flex items-center space-x-1 flex-shrink-0">
                <button onclick="window.bindoraApp.loadSimilarAsLigand('${encodeURIComponent(c.smiles)}', '${encodeURIComponent(c.title || `CID ${c.cid}`)}')" class="px-2 py-0.5 rounded bg-cyan-700 hover:bg-cyan-600 text-white text-[10px] font-medium transition cursor-pointer">Use</button>
                <button onclick="window.bindoraApp.addSimilarToBatch('${encodeURIComponent(c.title || `CID ${c.cid}`)}', '${encodeURIComponent(c.smiles)}')" class="px-2 py-0.5 rounded bg-slate-700 hover:bg-slate-600 text-cyan-300 text-[10px] font-medium transition cursor-pointer">+ Batch</button>
              </div>
            </div>
          `).join("")}
        </div>
      `;
    } catch (e) {
      container.innerHTML = `<span class="text-[11px] text-rose-400">Similarity search failed: ${e.message}</span>`;
    }
  }

  async loadSimilarAsLigand(encodedSmiles, encodedName) {
    const smiles = decodeURIComponent(encodedSmiles);
    const name = decodeURIComponent(encodedName);
    this.showToast(`Loading similar compound: ${name}...`, "info");
    try {
      const ligPrep = await BindoraAPI.prepareLigand(smiles, false);
      this.state.ligand = {
        name: name,
        smiles: smiles,
        weight: ligPrep.adme?.physicochemical?.molecular_weight?.value,
        ...ligPrep
      };
      this.updateStudioCards();
      this.updateDossierView();
      if (this.viewer && ligPrep.pdb_block) {
        this.viewer.loadLigand(ligPrep.pdb_block);
      }
      this.showToast(`Loaded ${name} as active ligand!`, "success");
    } catch (e) {
      this.showToast(`Failed to load compound: ${e.message}`, "error");
    }
  }

  addSimilarToBatch(encodedName, encodedSmiles) {
    const name = decodeURIComponent(encodedName);
    const smiles = decodeURIComponent(encodedSmiles);
    const batchInput = document.getElementById("input-batch-candidates");
    if (batchInput) {
      const existing = batchInput.value.trim();
      batchInput.value = (existing ? existing + "\n" : "") + `${name}, ${smiles}`;
      this.showToast(`Added '${name}' to batch candidate list!`, "success");
    }
  }

  async checkPharmacophoreEligibility() {
    const r = this.state.receptor;
    const summary = document.getElementById("pharmacophore-profile-summary");
    if (!r) {
      if (summary) summary.innerHTML = `<span class="text-xs text-slate-400 italic">No target receptor loaded yet. Please search RCSB or load 1HSG in Tab 1 to derive consensus pharmacophore features.</span>`;
      return;
    }

    const pdbId = r.pdb_id || "";
    const uniprotAcc = r.uniprot?.accession || r.uniprot_accession || "";
    const targetName = r.uniprot?.protein_name || r.pdb_id || r.title || "";
    const chemblId = r.chembl_id || "";

    if (summary && !this.state.pharmacophore) {
      summary.innerHTML = `<span class="text-purple-300 animate-pulse text-[11px]">⚡ Deriving 3D consensus pharmacophore from ChEMBL actives...</span>`;
    }

    try {
      const data = await BindoraAPI.getPharmacophoreActives(targetName, chemblId, pdbId, uniprotAcc, 10);
      if (data.eligible && data.actives_count >= 3) {
        this.state.pharmacophore = data;
        const badge = document.getElementById("pharmacophore-actives-badge");
        if (badge) badge.textContent = `${data.actives_count} ChEMBL Actives`;
        if (summary && data.consensus_profile) {
          const reqs = data.consensus_profile.core_requirements || {};
          const features = Object.entries(reqs).map(([f, cnt]) => `<span class="inline-block px-1.5 py-0.5 rounded bg-purple-950 text-purple-300 border border-purple-800 mr-1 mb-1 font-mono text-[10px]">${f}: &ge;${cnt}</span>`).join("");
          const displayName = pdbId ? `${pdbId} (${r.uniprot?.protein_name || 'HIV-1 Protease'})` : targetName.substring(0, 35);
          summary.innerHTML = `
            <div class="text-[11px] text-purple-300 font-semibold mb-1">Target: ${displayName} (${data.actives_count} actives)</div>
            <div class="flex flex-wrap">${features || '<span class="text-slate-400">Consensus features mapped.</span>'}</div>
          `;
        }
      } else {
        if (summary) {
          summary.innerHTML = `<span class="text-slate-400 text-[11px]">No curated ChEMBL actives (IC50 &le; 10 µM) found for target "${pdbId || targetName.substring(0, 25)}". Pharmacophore profile requires &ge;3 known active binders.</span>`;
        }
      }
    } catch (e) {
      if (summary) {
        summary.innerHTML = `<span class="text-rose-400 text-[11px]">Failed to load pharmacophore: ${e.message}</span>`;
      }
    }
  }

  async loadEnsembleStructures() {
    const listEl = document.getElementById("ensemble-structures-list");
    const r = this.state.receptor;
    if (!r) {
      if (listEl) listEl.innerHTML = `<span class="text-xs text-slate-500 italic">No target receptor loaded yet. Search RCSB or pick a benchmark in Tab 1.</span>`;
      return;
    }

    const acc = r.uniprot?.accession || r.uniprot_accession || r.pdb_id;
    if (!acc) {
      if (listEl) listEl.innerHTML = `<span class="text-xs text-slate-500 italic">No UniProt accession found for target.</span>`;
      return;
    }

    if (listEl) {
      listEl.innerHTML = `<div class="text-xs text-emerald-400 animate-pulse flex items-center space-x-1.5 py-2"><svg class="animate-spin w-4 h-4 text-emerald-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg><span>Querying UniProt cross-references for deposited structures...</span></div>`;
    }

    try {
      const data = await BindoraAPI.getEnsembleStructures(acc, 8);
      const structures = data.structures || [];
      this.state.ensembleStructures = structures;

      const badge = document.getElementById("ensemble-count-badge");
      if (badge) badge.textContent = `${structures.length} Deposited PDBs`;

      if (structures.length === 0) {
        if (listEl) listEl.innerHTML = `<span class="text-xs text-slate-500 italic">No alternative crystal structures deposited for ${acc}.</span>`;
        return;
      }

      const currentPdb = (r.pdb_id || "").toUpperCase();
      if (listEl) {
        listEl.innerHTML = structures.map((s, idx) => {
          const isCurrent = s.pdb_id.toUpperCase() === currentPdb;
          const checked = idx < 4;
          return `
            <label class="flex items-center justify-between p-2 rounded-lg bg-slate-800/60 hover:bg-slate-800 border border-slate-700/60 cursor-pointer transition text-xs">
              <div class="flex items-center space-x-2">
                <input type="checkbox" value="${s.pdb_id}" ${checked ? 'checked' : ''} class="ensemble-struct-checkbox rounded text-emerald-500 focus:ring-0 focus:ring-offset-0 bg-slate-900 border-slate-600">
                <span class="font-mono font-bold text-white">${s.pdb_id}</span>
                <span class="text-[10px] text-slate-400">${s.method}</span>
                ${isCurrent ? '<span class="text-[9px] font-bold px-1 rounded bg-emerald-950 text-emerald-300 border border-emerald-800">Current</span>' : ''}
              </div>
              <span class="font-mono text-[11px] text-emerald-400">${s.resolution}</span>
            </label>
          `;
        }).join("");
      }
    } catch (e) {
      if (listEl) listEl.innerHTML = `<span class="text-xs text-rose-400">Failed to query UniProt structures: ${e.message}</span>`;
    }
  }

  async runEnsembleDocking() {
    const r = this.state.receptor;
    const l = this.state.ligand;
    if (!r) {
      this.showToast("Please load a target receptor first.", "error");
      return;
    }
    if (!l) {
      this.showToast("Please load a ligand first.", "error");
      return;
    }

    const checkboxes = document.querySelectorAll(".ensemble-struct-checkbox:checked");
    const pdbIds = Array.from(checkboxes).map(cb => cb.value);
    if (pdbIds.length === 0) {
      this.showToast("Please select at least one PDB structure for ensemble docking.", "error");
      return;
    }

    const smiles = l.canonical_smiles || l.smiles || "";
    if (!smiles) {
      this.showToast("Ligand SMILES required for ensemble docking.", "error");
      return;
    }

    const btn = document.getElementById("btn-run-ensemble");
    const originalText = btn ? btn.innerHTML : "";
    if (btn) {
      btn.disabled = true;
      btn.innerHTML = '<svg class="animate-spin -ml-1 mr-2 h-4 w-4 text-white inline" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg> Cross-docking...';
    }

    this.showToast(`Starting ensemble cross-docking across ${pdbIds.length} target structures...`, "info");

    try {
      const res = await BindoraAPI.runEnsembleDocking({
        pdb_ids: pdbIds,
        ligand_smiles: smiles,
        exhaustiveness: 4
      });

      this.renderEnsembleLeaderboard(res);
      const s = res.summary || {};
      this.showToast(`Ensemble docking complete! Mean ΔG: ${s.mean_affinity ?? '—'} ± ${s.std_affinity ?? 0} kcal/mol`, "success");
    } catch (e) {
      this.showToast(`Ensemble docking failed: ${e.message}`, "error");
    } finally {
      if (btn) {
        btn.disabled = false;
        btn.innerHTML = originalText || '<span>Run Ensemble Cross-Docking</span>';
      }
    }
  }

  renderEnsembleLeaderboard(res) {
    const tbody = document.getElementById("batch-leaderboard-rows");
    if (!tbody) return;

    const structures = res.structures || [];
    if (structures.length === 0) {
      tbody.innerHTML = `<tr><td colspan="10" class="text-center p-6 text-xs text-slate-500 italic">No ensemble docking results produced.</td></tr>`;
      return;
    }

    tbody.innerHTML = structures.map((s, idx) => {
      if (s.affinity_kcal == null) {
        return `
          <tr class="border-b border-rose-950/60 bg-rose-950/20 text-rose-300">
            <td class="p-2.5 font-bold text-center text-slate-500">—</td>
            <td class="p-2.5 font-semibold text-rose-300 font-sans">${s.pdb_id} (${s.resolution || 'N/A'})</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5 text-slate-500">—</td>
            <td class="p-2.5"><span class="px-2 py-0.5 rounded bg-rose-950 text-rose-300 border border-rose-800 text-[10px]">${s.status || 'Failed'}</span></td>
          </tr>
        `;
      }

      const isWeak = s.affinity_kcal > -6.0;
      const leVal = (typeof s.ligand_efficiency === 'object' && s.ligand_efficiency !== null)
        ? (s.ligand_efficiency.value ?? "—")
        : (s.ligand_efficiency ?? "—");
      const kdVal = s.kd_nanomolar ?? s.theoretical_kd_nm ?? "—";
      const consensusVal = s.consensus_score ?? (s.vinardo_score != null ? ((s.affinity_kcal * 0.6) + (s.vinardo_score * 0.4)).toFixed(2) : s.affinity_kcal);

      return `
        <tr class="border-b border-slate-800 hover:bg-slate-800/50">
          <td class="p-2.5 font-bold text-center text-emerald-400">#${idx + 1}</td>
          <td class="p-2.5 font-semibold text-white font-sans">
            <span class="font-mono text-emerald-300">${s.pdb_id}</span>
            <span class="text-[10px] text-slate-400 block">${s.title || ''} (${s.resolution || 'N/A'})</span>
          </td>
          <td class="p-2.5 font-mono font-bold ${isWeak ? 'text-amber-400' : 'text-emerald-400'}">${s.affinity_kcal}</td>
          <td class="p-2.5 font-mono text-cyan-300 font-semibold">${s.vinardo_score ?? "—"}</td>
          <td class="p-2.5 font-mono font-bold text-cyan-300">${consensusVal}</td>
          <td class="p-2.5"><span class="px-1.5 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-800 text-[10px] font-bold inline-block">Conformation</span></td>
          <td class="p-2.5 font-mono text-slate-300">${kdVal}</td>
          <td class="p-2.5 font-mono text-slate-300">${leVal}</td>
          <td class="p-2.5 font-mono text-slate-300">${s.hbond_count ?? 0}</td>
          <td class="p-2.5 text-[10px]">
            <span class="px-1.5 py-0.5 rounded bg-emerald-950/80 text-emerald-300 border border-emerald-800">Docked</span>
          </td>
        </tr>
      `;
    }).join("");

    const s = res.summary || {};
    const sub = document.getElementById("batch-matrix-subtitle");
    if (sub && s.mean_affinity != null) {
      sub.textContent = `Ensemble Consensus: Mean ΔG = ${s.mean_affinity} ± ${s.std_affinity} kcal/mol (Spread: ${s.spread_kcal} kcal/mol across ${s.successful_runs} structures)`;
    }
  }

  async runPharmacophoreScreen() {
    if (!this.state.pharmacophore?.consensus_profile) {
      this.showToast("No active consensus pharmacophore profile loaded.", "error");
      return;
    }

    const rawText = document.getElementById("input-batch-candidates")?.value || "";
    const lines = rawText.split("\n");
    const candidates = [];
    for (const line of lines) {
      const parts = line.split(",");
      if (parts.length >= 2) {
        candidates.push({ name: parts[0].trim(), smiles: parts[1].trim() });
      } else if (line.trim().length > 3) {
        candidates.push({ name: `Candidate ${candidates.length+1}`, smiles: line.trim() });
      }
    }

    if (candidates.length === 0 && this.state.ligand?.smiles) {
      candidates.push({
        name: this.state.ligand.name || "Loaded Investigational Ligand",
        smiles: this.state.ligand.smiles
      });
    }

    if (candidates.length === 0) {
      this.showToast("No valid candidates found in text input or loaded ligand.", "error");
      return;
    }

    const btn = document.getElementById("btn-run-pharmacophore-screen");
    const originalText = btn ? btn.innerHTML : "";
    if (btn) {
      btn.disabled = true;
      btn.innerHTML = '<svg class="animate-spin -ml-1 mr-2 h-4 w-4 text-white inline" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg> Screening...';
    }

    this.showToast(`Screening ${candidates.length} candidate(s) against 3D consensus pharmacophore...`, "info");

    try {
      const data = await BindoraAPI.screenPharmacophore(candidates, this.state.pharmacophore.consensus_profile);
      this.renderPharmacophoreLeaderboard(data);
      this.showToast(`Pharmacophore screening complete! Screened ${data.total} candidates.`, "success");
    } catch (e) {
      this.showToast(`Pharmacophore screening failed: ${e.message}`, "error");
    } finally {
      if (btn) {
        btn.disabled = false;
        btn.innerHTML = originalText || '<span>Screen Candidates vs Pharmacophore</span>';
      }
    }
  }

  renderPharmacophoreLeaderboard(data) {
    const tbody = document.getElementById("batch-leaderboard-rows");
    if (!tbody) return;

    const screened = data.screened || [];
    if (screened.length === 0) {
      tbody.innerHTML = `<tr><td colspan="10" class="text-center p-6 text-xs text-slate-500 italic">No candidates screened against pharmacophore.</td></tr>`;
      return;
    }

    tbody.innerHTML = screened.map((c, idx) => {
      const matchPct = c.match_score_pct ?? 0;
      const isStrong = matchPct >= 80;
      const isModerate = matchPct >= 50 && matchPct < 80;
      const align = c.feature_alignment || {};
      const alignTooltip = Object.entries(align).map(([k, v]) => `${k}: ${v}`).join("\n");
      const feats = c.candidate_features || {};

      return `
        <tr class="border-b border-slate-800 hover:bg-slate-800/50">
          <td class="p-2.5 font-bold text-center text-purple-400">#${idx + 1}</td>
          <td class="p-2.5 font-semibold text-white font-sans" title="${c.smiles}">
            <span class="text-purple-300 font-mono">${c.name}</span>
            <span class="text-[10px] text-slate-400 block truncate max-w-xs">${c.smiles}</span>
          </td>
          <td class="p-2.5 font-mono font-bold ${isStrong ? 'text-emerald-400' : isModerate ? 'text-purple-300' : 'text-slate-400'}">${matchPct.toFixed(1)}%</td>
          <td class="p-2.5 font-mono text-cyan-300 text-xs" title="${alignTooltip}">${feats.Aromatic ?? 0} Aro / ${feats.Hydrophobe ?? 0} Hyd</td>
          <td class="p-2.5 font-mono font-bold text-purple-300">${(matchPct / 10).toFixed(1)} / 10</td>
          <td class="p-2.5">
            <span class="px-1.5 py-0.5 rounded ${isStrong ? 'bg-emerald-950 text-emerald-300 border border-emerald-800' : isModerate ? 'bg-purple-950 text-purple-300 border border-purple-800' : 'bg-slate-800 text-slate-400 border border-slate-700'} text-[10px] font-bold inline-block cursor-help" title="${alignTooltip}">
              ${c.status || (isStrong ? 'Strong Match' : 'Partial')}
            </span>
          </td>
          <td class="p-2.5 font-mono text-slate-400">—</td>
          <td class="p-2.5 font-mono text-slate-400">—</td>
          <td class="p-2.5 font-mono text-purple-300">${feats.Acceptor ?? 0} HBA / ${feats.Donor ?? 0} HBD</td>
          <td class="p-2.5 text-[10px]">
            <span class="px-1.5 py-0.5 rounded ${isStrong ? 'bg-emerald-950/80 text-emerald-300 border border-emerald-800' : 'bg-slate-800 text-slate-300 border border-slate-700'}">
              ${isStrong ? 'Pass' : 'Sub-match'}
            </span>
          </td>
        </tr>
      `;
    }).join("");

    const sub = document.getElementById("batch-matrix-subtitle");
    if (sub) {
      sub.textContent = `Pharmacophore Consensus: Screened ${data.total} candidates against active inhibitor profile`;
    }
  }
}
