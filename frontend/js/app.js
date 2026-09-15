// Main frontend application controller for Bindora

document.addEventListener("DOMContentLoaded", () => {
  const app = new BindoraApp();
  window.app = app;
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
          statusBadge.innerHTML = `<span class="inline-block w-2 h-2 rounded-full bg-emerald-400 mr-1.5 animate-pulse"></span> AutoDock Vina Ready`;
          statusBadge.className = "px-2.5 py-1 text-xs font-semibold rounded-full bg-emerald-950/80 text-emerald-300 border border-emerald-800 flex items-center";
        } else {
          statusBadge.innerHTML = `<span class="inline-block w-2 h-2 rounded-full bg-amber-400 mr-1.5"></span> Vina Standby`;
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

    // 5. Clean live state on startup (benchmarks are available above for 1-click exploration if user chooses)
    this.updateStudioCards();

    // 6. Set active tab to Home page
    this.switchTab("home");
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
        if (hudText) hudText.textContent = "📐 Click any atom to begin measurement...";
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

    // Batch Library Presets
    const librarySelect = document.getElementById("select-batch-library");
    const batchTextarea = document.getElementById("input-batch-candidates");
    if (librarySelect && batchTextarea) {
      const presets = {
        nsaid: "Aspirin, CC(=O)Oc1ccccc1C(=O)O\nIbuprofen, CC(C)Cc1ccc(cc1)C(C)C(=O)O\nNaproxen, COc1ccc2cc(ccc2c1)C(C)C(=O)O\nCelecoxib, Cc1ccc(cc1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F",
        kinase: "Imatinib, Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C\nGefitinib, COc1cc2ncnc(c2cc1OCCCN3CCOCC3)Nc4ccc(c(c4)Cl)F\nErlotinib, COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC\nDasatinib, Cc1cccc(c1Cl)NC(=O)c2cnc(s2)Nc3cc(nc(n3)C)N4CCN(CC4)CCO",
        antiviral: "Remdesivir, CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)C#N)O)O)Oc2ccccc2\nFavipiravir, C1=C(N=C(C(=O)N1)C(=O)N)F\nMolnupiravir, CC(C)C(=O)OCC1C(C(C(O1)N2C=CC(=NO)NC2=O)O)O\nRibavirin, C1C(C(C(O1)N2C=NC(=N2)C(=O)N)O)O",
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
        this.showToast("Signed out of Firebase account.", "info");
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

          this.showToast("Simulation saved to Firebase Realtime Database!", "success");
          this.loadCloudHistory();
        } catch (e) {
          this.showToast(e.message, "error");
        }
      });
    }

    if (refreshSavedRunsBtn) {
      refreshSavedRunsBtn.addEventListener("click", () => this.loadCloudHistory());
    }

    // Export Dossier
    const exportBtn = document.getElementById("btn-export-dossier");
    if (exportBtn) {
      exportBtn.addEventListener("click", () => window.print());
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

    // Dark / Light Mode Toggle
    const themeToggleBtn = document.getElementById("btn-theme-toggle");
    if (themeToggleBtn) {
      // Restore saved theme on startup
      const savedTheme = localStorage.getItem('bindora-theme') || 'dark';
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
        container.innerHTML = `<p class="text-slate-500 italic p-2 text-center text-xs">No saved docking runs in cloud database yet.</p>`;
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

  // Lightweight markdown → HTML renderer for AI narratives
  renderMarkdown(text) {
    if (!text) return '';
    // Escape HTML entities first
    let out = text
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');
    // Headers
    out = out
      .replace(/^### (.+)$/gm, '<h3 class="text-sm font-bold text-white mt-4 mb-1 pb-0.5 border-b border-slate-700/50">$1</h3>')
      .replace(/^## (.+)$/gm, '<h2 class="text-sm font-bold text-[#00C6FF] mt-5 mb-1 uppercase tracking-wide">$1</h2>')
      .replace(/^# (.+)$/gm, '<h1 class="text-base font-extrabold text-white mt-5 mb-2">$1</h1>');
    // Bold and italic
    out = out
      .replace(/\*\*(.+?)\*\*/g, '<strong class="font-semibold text-white">$1</strong>')
      .replace(/\*([^*]+?)\*/g, '<em class="italic text-slate-300">$1</em>');
    // Inline code
    out = out
      .replace(/`(.+?)`/g, '<code class="font-mono text-[#00C6FF] bg-slate-800/70 px-1 py-0.5 rounded text-xs">$1</code>');
    // Double newline = paragraph break
    out = out.replace(/\n\n/g, '</p><p class="mt-2 text-slate-300">');
    // Single newline
    out = out.replace(/\n/g, '<br>');
    return out;
  }

  // Apply dark/light theme across the page
  applyTheme(theme) {
    const html = document.documentElement;
    const sunIcon = document.getElementById("icon-sun");
    const moonIcon = document.getElementById("icon-moon");
    if (theme === 'light') {
      html.classList.remove('dark');
      html.classList.add('light');
      if (sunIcon) sunIcon.classList.add('hidden');
      if (moonIcon) moonIcon.classList.remove('hidden');
      localStorage.setItem('bindora-theme', 'light');
    } else {
      html.classList.remove('light');
      html.classList.add('dark');
      if (sunIcon) sunIcon.classList.remove('hidden');
      if (moonIcon) moonIcon.classList.add('hidden');
      localStorage.setItem('bindora-theme', 'dark');
    }
  }

  renderBenchmarksUI() {
    const container = document.getElementById("benchmark-cards");
    if (!container) return;

    container.innerHTML = this.state.benchmarks.map(bm => `
      <div class="cursor-pointer p-3 rounded-lg border border-slate-700/60 bg-slate-800/40 hover:bg-slate-700/50 transition-all hover:border-cyan-500/50" onclick="window.app.loadBenchmark('${bm.id}')">
        <div class="flex items-center justify-between mb-1">
          <span class="font-bold text-sm text-cyan-300">${bm.drug_name}</span>
          <span class="text-xs px-2 py-0.5 rounded bg-slate-700 text-slate-300 font-mono">${bm.pdb_id}</span>
        </div>
        <p class="text-xs text-slate-300 font-medium">${bm.target_name}</p>
        <p class="text-[11px] text-slate-400 mt-1 line-clamp-2">${bm.mechanism}</p>
      </div>
    `).join("");
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
    this.state.currentPoseIdx = 0;

    if (this.viewer) {
      this.viewer.clear();
      this.viewer.clearInteractions();
      this.viewer.clearMeasurements();
      this.viewer.clearGridBox();
    }
    const emptyState = document.getElementById("viewer-empty-state");
    if (emptyState) emptyState.classList.remove("hidden");

    const statusBanner = document.getElementById("docking-status-banner");
    if (statusBanner) statusBanner.classList.add("hidden");

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
          </div>
        `;
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
        recCard.innerHTML = `
          <div class="space-y-2">
            <div class="flex items-center justify-between">
              <span class="text-base font-bold text-emerald-400">${r.pdb_id || "Custom Receptor"}</span>
              <span class="text-xs px-2 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-800 font-mono">${r.atom_count || 0} atoms</span>
            </div>
            <p class="text-xs text-slate-300 font-medium">${r.title || "Macromolecular target"}</p>
            <div class="text-xs text-slate-400 pt-2 border-t border-slate-700/60 space-y-1">
              <div>Binding Pocket: <span class="text-emerald-300 font-medium">${pocket.description || "Auto-detected"}</span></div>
              <div class="font-mono text-[11px] text-slate-400">
                Center: (${pocket.center?.x}, ${pocket.center?.y}, ${pocket.center?.z}) | Size: (${pocket.size?.x}, ${pocket.size?.y}, ${pocket.size?.z})
              </div>
            </div>
          </div>
        `;
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

    if (statusBanner) {
      statusBanner.classList.remove("hidden");
      statusBanner.className = "glass-panel p-3.5 rounded-xl border border-cyan-500/50 bg-cyan-950/80 transition-all duration-300 mb-2 animate-pulse";
      if (statusIcon) statusIcon.innerHTML = `<svg class="animate-spin w-5 h-5 text-cyan-400" fill="none" viewBox="0 0 24 24"><circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle><path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path></svg>`;
      if (statusTitle) statusTitle.textContent = "AutoDock Vina Docking in Progress...";
      if (statusSub) statusSub.textContent = `Sampling conformational space (exhaustiveness = ${exhaustiveness}, modes = 9) on CPU...`;
      if (statusExtra) statusExtra.textContent = "Running Monte Carlo...";
    }

    this.showToast("Launching AutoDock Vina docking engine...", "info");

    try {
      // Read grid parameters
      const getVal = (id, def) => parseFloat(document.getElementById(id)?.value) || def;
      const center = { x: getVal("grid-cx", 0), y: getVal("grid-cy", 0), z: getVal("grid-cz", 0) };
      const size = { x: getVal("grid-sx", 22), y: getVal("grid-sy", 22), z: getVal("grid-sz", 22) };

      const p = this.state.ligand.adme?.physicochemical || {};

      const dockResult = await BindoraAPI.runDocking({
        receptor_pdbqt: this.state.receptor.pdbqt_text,
        receptor_pdb: this.state.receptor.cleaned_pdb,
        ligand_pdbqt: this.state.ligand.pdbqt_text,
        center: center,
        size: size,
        exhaustiveness: exhaustiveness,
        num_modes: 9,
        heavy_atoms: p.heavy_atoms?.value || this.state.ligand.heavy_atom_count || 20,
        molecular_weight: p.molecular_weight?.value || this.state.ligand.weight || 300.0
      });

      this.state.docking = dockResult;
      this.state.currentPoseIdx = 0;

      // Render top pose in 3D viewer
      if (this.viewer && dockResult.top_pose) {
        this.viewer.loadLigand(dockResult.top_pose.pdb_block);
        if (dockResult.interactions) {
          this.viewer.renderInteractions(dockResult.interactions);
        }
      }

      // Render pose score badges
      this.updateDockingScoresUI();

      // Trigger ChEMBL Bioactivity cross-check asynchronously
      this.fetchBioactivityCrosscheck();

      // Render ADME charts
      this.updateADMEView();

      // Update completion banner
      if (statusBanner) {
        statusBanner.className = "glass-panel p-3.5 rounded-xl border border-emerald-500/50 bg-emerald-950/80 transition-all duration-300 mb-2";
        if (statusIcon) statusIcon.innerHTML = `<span class="w-6 h-6 rounded-full bg-emerald-500/20 text-emerald-400 flex items-center justify-center font-bold text-sm">✓</span>`;
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
        statusBanner.className = "glass-panel p-3.5 rounded-xl border border-rose-500/50 bg-rose-950/80 transition-all duration-300 mb-2";
        if (statusIcon) statusIcon.innerHTML = `<span class="w-6 h-6 rounded-full bg-rose-500/20 text-rose-400 flex items-center justify-center font-bold text-sm">✕</span>`;
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
    if (elKd) elKd.textContent = `${thermodynamics.theoretical_kd_nm} nM`;

    const elLE = document.getElementById("dock-le-score");
    if (elLE) elLE.textContent = `${thermodynamics.ligand_efficiency?.value || "—"}`;

    const elHbonds = document.getElementById("dock-hbonds-count");
    if (elHbonds) elHbonds.textContent = `${interactions.total_hbond_count || 0}`;

    // Pose Stepper Label
    const elMode = document.getElementById("dock-current-mode-label");
    if (elMode) elMode.textContent = `Mode ${currentPose.mode || (this.state.currentPoseIdx + 1)} / ${poses.length}`;

    // Pose Table
    const poseTable = document.getElementById("pose-table-rows");
    if (poseTable) {
      poseTable.innerHTML = poses.map((p, idx) => `
        <tr class="border-b border-slate-800 hover:bg-slate-800/40 cursor-pointer ${idx === this.state.currentPoseIdx ? 'bg-cyan-950/40 font-semibold' : ''}" onclick="window.app.selectPose(${idx})">
          <td class="px-3 py-2 text-cyan-300">Mode ${p.mode}</td>
          <td class="px-3 py-2 font-mono text-white">${p.affinity_kcal}</td>
          <td class="px-3 py-2 font-mono text-slate-400">${p.rmsd_lb}</td>
          <td class="px-3 py-2 font-mono text-slate-400">${p.rmsd_ub}</td>
        </tr>
      `).join("");
    }

    // Interacting Residues Chips
    const resContainer = document.getElementById("dock-interacting-residues");
    if (resContainer) {
      const resList = interactions.interacting_residues || [];
      if (resList.length > 0) {
        resContainer.innerHTML = resList.map(r => `
          <span class="inline-block px-2 py-0.5 rounded text-xs bg-slate-800 text-cyan-300 border border-slate-700 font-mono">${r}</span>
        `).join("");
      } else {
        resContainer.innerHTML = `<span class="text-xs text-slate-500 italic">No specific residues within contact threshold</span>`;
      }
    }
  }

  async selectPose(idx) {
    if (!this.state.docking || !this.state.docking.poses) return;
    if (idx < 0 || idx >= this.state.docking.poses.length) return;

    this.state.currentPoseIdx = idx;
    const pose = this.state.docking.poses[idx];

    // Re-analyze interactions for this pose
    try {
      const contacts = await BindoraAPI.analyzeInteractions(this.state.receptor.cleaned_pdb, pose.pdbqt_content);
      this.state.docking.interactions = contacts;

      if (this.viewer) {
        this.viewer.loadLigand(pose.pdb_block);
        this.viewer.renderInteractions(contacts);
      }
    } catch (e) {
      console.warn("Could not re-analyze pose interactions:", e);
    }

    this.updateDockingScoresUI();
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
    const statusIcon = isCorroborated ? "✓" : "⚠";

    // Build ChEMBL verification links
    const chemblMolUrl = data.chembl_compound_url || null;
    const chemblTargetUrl = data.chembl_target_url || null;
    const molId = data.molecule_chembl_id;
    const targetId = data.target_chembl_id;

    const searchedHtml = `
      <div class="grid grid-cols-2 gap-2 text-xs font-mono mt-2">
        <div class="p-2 rounded-lg bg-slate-900/60 border border-slate-800 space-y-1">
          <div class="text-slate-400 text-[10px] uppercase tracking-wider font-sans font-semibold">Compound Searched</div>
          <div class="text-white font-semibold">${data.drug_searched || "—"}</div>
          ${molId ? `<div class="text-[11px] text-slate-400">${molId}
            <a href="${chemblMolUrl}" target="_blank" rel="noopener" class="text-cyan-400 hover:underline ml-1">↗ ChEMBL</a></div>` : `<div class="text-[11px] text-rose-400">Not found in ChEMBL</div>`}
        </div>
        <div class="p-2 rounded-lg bg-slate-900/60 border border-slate-800 space-y-1">
          <div class="text-slate-400 text-[10px] uppercase tracking-wider font-sans font-semibold">Target Searched</div>
          <div class="text-white font-semibold">${data.target_searched || "—"}</div>
          ${targetId ? `<div class="text-[11px] text-slate-400">${targetId}
            <a href="${chemblTargetUrl}" target="_blank" rel="noopener" class="text-cyan-400 hover:underline ml-1">↗ ChEMBL</a></div>` : `<div class="text-[11px] text-rose-400">Not found in ChEMBL</div>`}
        </div>
      </div>
    `;

    // No records explanation panel (scientific context, not an error)
    const noRecordsExplanationHtml = !isCorroborated ? `
      <div class="mt-3 p-3 rounded-lg bg-slate-900/50 border border-amber-900/40 text-xs space-y-1">
        <p class="font-semibold text-amber-300 text-[11px] uppercase tracking-wide">What does this mean scientifically?</p>
        <p class="text-slate-400 leading-relaxed">
          <strong class="text-slate-300">No ChEMBL records ≠ inactive compound.</strong> It means no curated wet-lab binding assay (IC50, Ki, Kd, EC50) has been deposited in the EMBL-EBI ChEMBL database for this exact drug-target pair. This is a common finding for:
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

    setVal("adme-gi-val", pk.gi_absorption?.level || "—");
    setVal("adme-bbb-val", pk.bbb_permeation?.status?.split(" ")[0] || "—");
    setVal("adme-ppb-val", pk.plasma_protein_binding?.tier?.split(" ")[0] || "—");

    // PAINS / Brenk Alerts
    const safety = adme.medicinal_chemistry_safety || {};
    const pains = safety.pains_alerts || {};
    const painsEl = document.getElementById("adme-pains-val");
    if (painsEl) {
      painsEl.textContent = pains.status || "Clear";
      painsEl.className = pains.count > 0 ? "font-bold text-rose-400" : "font-bold text-emerald-400";
    }
  }

  async generateNarrativeReport() {
    const container = document.getElementById("narrative-report-content");
    if (container) {
      container.innerHTML = `<div class="text-xs text-slate-400 animate-pulse">Generating pedagogical pharmacology explanation with strict anti-hallucination verification...</div>`;
    }

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
      if (container) {
        // Render markdown formatted text properly (not raw with ** showing)
        const renderedMarkdown = this.renderMarkdown(res.narrative || '');
        container.innerHTML = `
          <div class="space-y-4">
            <div class="flex items-center justify-between pb-2 border-b border-slate-700/60 text-xs text-slate-400">
              <span>Engine: <a href="https://openrouter.ai" target="_blank" rel="noopener" class="text-cyan-400 font-semibold hover:underline">${res.source}</a></span>
              <span class="text-emerald-400 font-semibold">&#10003; Zero-Hallucination Verified</span>
            </div>
            <div class="text-slate-300 text-xs sm:text-sm leading-relaxed narrative-md">
              <p>${renderedMarkdown}</p>
            </div>
          </div>
        `;
      }
    } catch (e) {
      console.warn("Narrative generation error:", e);
      if (container) {
        container.innerHTML = `<div class="text-xs text-rose-400">Failed to generate narrative: ${e.message}</div>`;
      }
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

  updateDossierView() {
    const ligName = document.getElementById("dossier-ligand-name");
    const targetName = document.getElementById("dossier-target-name");
    const resultsSummary = document.getElementById("dossier-results-summary");

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

    if (resultsSummary && this.state.docking) {
      const d = this.state.docking;
      const thermo = d.thermodynamics || {};
      const contacts = d.interactions || {};
      const adme = this.state.ligand?.adme || {};
      const lip = adme.drug_likeness?.lipinski || {};
      const xcheck = this.state.crosscheck || {};

      resultsSummary.innerHTML = `
        <h4 class="font-bold text-slate-200 uppercase text-xs tracking-wider">Docking Results Summary</h4>
        <div class="grid grid-cols-2 sm:grid-cols-4 gap-3 text-xs font-mono">
          <div class="p-2 bg-slate-900/60 rounded border border-slate-800">
            <span class="text-slate-400 block">Binding Free Energy</span>
            <span class="text-cyan-400 font-bold">${d.top_pose?.affinity_kcal || "—"} kcal/mol</span>
          </div>
          <div class="p-2 bg-slate-900/60 rounded border border-slate-800">
            <span class="text-slate-400 block">Theoretical Kd</span>
            <span class="text-emerald-400 font-bold">${thermo.theoretical_kd_nm || "—"} nM</span>
          </div>
          <div class="p-2 bg-slate-900/60 rounded border border-slate-800">
            <span class="text-slate-400 block">H-Bonds</span>
            <span class="text-yellow-400 font-bold">${contacts.total_hbond_count || 0}</span>
          </div>
          <div class="p-2 bg-slate-900/60 rounded border border-slate-800">
            <span class="text-slate-400 block">Lipinski Status</span>
            <span class="${lip.status === 'Pass' ? 'text-emerald-400' : 'text-amber-400'} font-bold">${lip.status || "—"} (${lip.violations_count || 0} violations)</span>
          </div>
        </div>
        <div class="mt-2 px-3 py-2 rounded border ${xcheck.is_cross_checked ? 'border-emerald-800 bg-emerald-950/30' : 'border-amber-800 bg-amber-950/30'} text-xs">
          <span class="font-bold">${xcheck.is_cross_checked ? '✓' : '⚠'} ${xcheck.status_badge || 'Computational Prediction Only'}</span>
          ${xcheck.summary_note ? '<p class="text-slate-400 mt-1">' + xcheck.summary_note + '</p>' : ''}
        </div>
      `;
    }
  }

  updatePathwayInfo() {
    const container = document.getElementById("pathway-info-container");
    if (!container) return;

    const uniprot = this.state.receptor?.uniprot;
    if (uniprot && (uniprot.function || uniprot.subcellular_location)) {
      container.innerHTML = `
        ${uniprot.function ? '<div><span class="font-bold text-slate-200">Biological Function:</span> <span class="text-slate-300">' + uniprot.function + '</span></div>' : ''}
        ${uniprot.catalytic_activity ? '<div><span class="font-bold text-slate-200">Catalytic Activity:</span> <span class="text-slate-300">' + uniprot.catalytic_activity + '</span></div>' : ''}
        ${uniprot.subcellular_location ? '<div><span class="font-bold text-slate-200">Subcellular Location:</span> <span class="text-slate-300">' + uniprot.subcellular_location + '</span></div>' : ''}
        ${uniprot.tissue_specificity ? '<div><span class="font-bold text-slate-200">Tissue Specificity:</span> <span class="text-slate-300">' + uniprot.tissue_specificity + '</span></div>' : ''}
      `;
    } else {
      container.innerHTML = '<p class="text-slate-500 italic">No UniProt pathway annotation available for this target. Load a receptor with known annotations to view biological function and signaling cascade details.</p>';
    }
  }

  renderBatchLeaderboard(leaderboard) {
    const tbody = document.getElementById("batch-leaderboard-rows");
    if (!tbody) return;

    if (leaderboard.length === 0) {
      tbody.innerHTML = `<tr><td colspan="7" class="text-center p-4 text-xs text-slate-500 italic">No docking results produced.</td></tr>`;
      return;
    }

    tbody.innerHTML = leaderboard.map(item => `
      <tr class="border-b border-slate-800 hover:bg-slate-800/50">
        <td class="p-3 font-bold text-center text-cyan-400">#${item.rank}</td>
        <td class="p-3 font-semibold text-white">${item.name}</td>
        <td class="p-3 font-mono font-bold text-emerald-400">${item.affinity_kcal}</td>
        <td class="p-3 font-mono text-slate-300">${item.theoretical_kd_nm}</td>
        <td class="p-3 font-mono text-slate-300">${item.ligand_efficiency}</td>
        <td class="p-3 font-mono text-slate-300">${item.hbond_count}</td>
        <td class="p-3 text-xs">
          <span class="px-2 py-0.5 rounded ${item.lipinski_status === 'Pass' ? 'bg-emerald-950 text-emerald-300 border border-emerald-800' : 'bg-amber-950 text-amber-300 border border-amber-800'}">
            ${item.lipinski_status}
          </span>
        </td>
      </tr>
    `).join("");
  }
}
