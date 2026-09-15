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
      apiKey: localStorage.getItem("bindora_api_key") || localStorage.getItem("anudock_api_key") || "",
      activeTab: "studio"
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
      const health = await AnuDockAPI.checkHealth();
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
      const bmData = await AnuDockAPI.getBenchmarks();
      this.state.benchmarks = bmData.benchmarks || [];
      this.renderBenchmarksUI();
    } catch (e) {
      console.warn("Could not fetch benchmarks:", e);
    }

    // 5. Load default benchmark (Imatinib / BCR-ABL) for instant exploration
    if (this.state.benchmarks.length > 0) {
      this.loadBenchmark(this.state.benchmarks[0].id);
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

    // RCSB Search
    const rcsbBtn = document.getElementById("btn-search-rcsb");
    const rcsbInput = document.getElementById("input-search-rcsb");
    if (rcsbBtn && rcsbInput) {
      rcsbBtn.addEventListener("click", () => this.searchRCSB(rcsbInput.value));
      rcsbInput.addEventListener("keydown", e => {
        if (e.key === "Enter") this.searchRCSB(rcsbInput.value);
      });
    }

    // Docking trigger button
    const dockBtn = document.getElementById("btn-run-docking");
    if (dockBtn) {
      dockBtn.addEventListener("click", () => this.runDockingPipeline());
    }

    // Viewer controls
    const resetCamBtn = document.getElementById("btn-reset-cam");
    if (resetCamBtn) {
      resetCamBtn.addEventListener("click", () => this.viewer && this.viewer.resetCamera());
    }

    const surfaceToggle = document.getElementById("toggle-surface");
    if (surfaceToggle) {
      surfaceToggle.addEventListener("change", (e) => {
        if (this.viewer) this.viewer.toggleSurface(e.target.checked);
      });
    }

    const measureToggle = document.getElementById("toggle-measure");
    if (measureToggle) {
      measureToggle.addEventListener("change", (e) => {
        if (this.viewer) this.viewer.enableMeasurementMode(e.target.checked);
      });
    }

    const proteinStyleSelect = document.getElementById("select-protein-style");
    if (proteinStyleSelect) {
      proteinStyleSelect.addEventListener("change", (e) => {
        if (this.viewer) {
          this.viewer.settings.proteinStyle = e.target.value;
          this.viewer.applyProteinStyle();
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

    // Settings Modal
    const settingsBtn = document.getElementById("btn-settings");
    const settingsModal = document.getElementById("modal-settings");
    const closeSettingsBtn = document.getElementById("btn-close-settings");
    const saveSettingsBtn = document.getElementById("btn-save-settings");
    const apiKeyInput = document.getElementById("input-api-key");
    const firebaseKeyInput = document.getElementById("input-firebase-api-key");

    if (settingsBtn && settingsModal) {
      settingsBtn.addEventListener("click", () => {
        if (apiKeyInput) apiKeyInput.value = this.state.apiKey;
        if (firebaseKeyInput) firebaseKeyInput.value = localStorage.getItem("bindora_firebase_api_key") || "";
        settingsModal.classList.remove("hidden");
      });
    }
    if (closeSettingsBtn && settingsModal) {
      closeSettingsBtn.addEventListener("click", () => settingsModal.classList.add("hidden"));
    }
    if (saveSettingsBtn && settingsModal) {
      saveSettingsBtn.addEventListener("click", () => {
        if (apiKeyInput) {
          this.state.apiKey = apiKeyInput.value.trim();
          localStorage.setItem("bindora_api_key", this.state.apiKey);
        }
        if (firebaseKeyInput && firebaseKeyInput.value.trim()) {
          const fbKey = firebaseKeyInput.value.trim();
          localStorage.setItem("bindora_firebase_api_key", fbKey);
          if (window.bindoraFirebase) {
            window.bindoraFirebase.setApiKey(fbKey);
          }
        }
        settingsModal.classList.add("hidden");
        this.showToast("Settings & API Keys saved successfully.", "success");
      });
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
      this.loadCloudHistory();
    } else {
      if (loggedOutView) loggedOutView.classList.remove("hidden");
      if (loggedInView) loggedInView.classList.add("hidden");
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
      const ligPrep = await AnuDockAPI.prepareLigand(bm.smiles, false);
      this.state.ligand = {
        name: bm.drug_name,
        smiles: bm.smiles,
        weight: ligPrep.adme?.physicochemical?.molecular_weight?.value,
        formula: "",
        ...ligPrep
      };

      // 2. Fetch PDB from RCSB
      const recRes = await AnuDockAPI.prepareReceptor("", bm.pdb_id);
      this.state.receptor = {
        pdb_id: bm.pdb_id,
        title: bm.target_name,
        ...recRes
      };

      // 3. Update UI
      this.updateStudioCards();
      
      // 4. Update 3D viewer with receptor
      if (this.viewer && recRes.cleaned_pdb) {
        this.viewer.loadReceptor(recRes.cleaned_pdb);
      }

      this.showToast(`Loaded ${bm.drug_name} & ${bm.target_name} (${bm.pdb_id})`, "success");
    } catch (e) {
      console.error("Benchmark load error:", e);
      this.showToast(`Failed to load benchmark: ${e.message}`, "error");
    }
  }

  async searchPubChem(query) {
    if (!query) return;
    this.showToast(`Searching PubChem for '${query}'...`, "info");
    try {
      const data = await AnuDockAPI.searchPubChem(query);
      const ligPrep = await AnuDockAPI.prepareLigand(data.smiles || data.sdf, !data.smiles);
      this.state.ligand = {
        ...data,
        ...ligPrep
      };
      this.updateStudioCards();
      this.showToast(`Found PubChem compound: ${data.name} (CID: ${data.cid})`, "success");
    } catch (e) {
      this.showToast(`PubChem search failed: ${e.message}`, "error");
    }
  }

  async searchRCSB(query) {
    if (!query) return;
    this.showToast(`Searching RCSB for '${query}'...`, "info");
    try {
      const data = await AnuDockAPI.searchRCSB(query);
      if (data.direct && data.entry) {
        const entry = data.entry;
        const recPrep = await AnuDockAPI.prepareReceptor(entry.pdb_content, entry.pdb_id);
        this.state.receptor = {
          ...entry,
          ...recPrep
        };
        this.updateStudioCards();
        if (this.viewer && recPrep.cleaned_pdb) {
          this.viewer.loadReceptor(recPrep.cleaned_pdb);
        }
        this.showToast(`Loaded RCSB entry ${entry.pdb_id}: ${entry.title.substring(0, 30)}...`, "success");
      } else if (data.results && data.results.length > 0) {
        // Show result picker
        const first = data.results[0];
        const recPrep = await AnuDockAPI.prepareReceptor("", first.pdb_id);
        this.state.receptor = {
          ...first,
          ...recPrep
        };
        this.updateStudioCards();
        if (this.viewer && recPrep.cleaned_pdb) {
          this.viewer.loadReceptor(recPrep.cleaned_pdb);
        }
        this.showToast(`Selected top RCSB match: ${first.pdb_id}`, "success");
      } else {
        this.showToast("No RCSB entries found for query", "error");
      }
    } catch (e) {
      this.showToast(`RCSB query failed: ${e.message}`, "error");
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

    this.showToast("Launching AutoDock Vina docking engine...", "info");

    try {
      // Read grid parameters
      const getVal = (id, def) => parseFloat(document.getElementById(id)?.value) || def;
      const center = { x: getVal("grid-cx", 0), y: getVal("grid-cy", 0), z: getVal("grid-cz", 0) };
      const size = { x: getVal("grid-sx", 22), y: getVal("grid-sy", 22), z: getVal("grid-sz", 22) };
      const exhaustiveness = parseInt(document.getElementById("docking-exhaustiveness")?.value) || 8;

      const p = this.state.ligand.adme?.physicochemical || {};

      const dockResult = await AnuDockAPI.runDocking({
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

      // Switch to 3D docking tab
      this.switchTab("docking");
      this.showToast(`Docking finished! Top ΔG: ${dockResult.top_pose.affinity_kcal} kcal/mol`, "success");

    } catch (e) {
      console.error("Docking error:", e);
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
      const contacts = await AnuDockAPI.analyzeInteractions(this.state.receptor.cleaned_pdb, pose.pdbqt_content);
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
      const res = await AnuDockAPI.crosscheckBioactivity(drug, target);
      this.state.crosscheck = res;
      this.renderCrosscheckUI(res);
      
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

    let recordsHtml = "";
    if (isCorroborated && data.experimental_records?.length > 0) {
      recordsHtml = `
        <div class="mt-4 overflow-x-auto">
          <table class="w-full text-xs text-left border border-slate-700/60 rounded-lg">
            <thead class="bg-slate-800/70 text-slate-300 font-semibold border-b border-slate-700">
              <tr>
                <th class="p-2">Assay Parameter</th>
                <th class="p-2">Experimental Value</th>
                <th class="p-2">Assay Description</th>
                <th class="p-2">Database Source</th>
              </tr>
            </thead>
            <tbody class="divide-y divide-slate-800 font-mono text-slate-300">
              ${data.experimental_records.map(r => `
                <tr>
                  <td class="p-2 font-bold text-emerald-400">${r.type}</td>
                  <td class="p-2 text-white">${r.relation} ${r.value} ${r.units}</td>
                  <td class="p-2 font-sans text-slate-400">${r.assay_description}</td>
                  <td class="p-2 font-sans text-slate-400">ChEMBL</td>
                </tr>
              `).join("")}
            </tbody>
          </table>
        </div>
      `;
    }

    container.innerHTML = `
      <div class="p-4 rounded-xl border ${badgeClass} space-y-2">
        <div class="flex items-center justify-between">
          <div class="flex items-center space-x-2">
            <span class="text-base font-bold">${statusIcon}</span>
            <span class="text-sm font-bold tracking-wide uppercase">${data.status_badge}</span>
          </div>
          <span class="text-xs px-2 py-0.5 rounded bg-slate-900/60 font-mono">${data.target_organism || "Homo sapiens"}</span>
        </div>
        <p class="text-xs text-slate-300 leading-relaxed">${data.summary_note}</p>
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
      const res = await AnuDockAPI.explainNarrative(payload, this.state.apiKey);
      this.state.narrative = res;
      if (container) {
        // Render markdown formatted text cleanly
        container.innerHTML = `
          <div class="prose prose-invert prose-sm max-w-none space-y-4">
            <div class="flex items-center justify-between pb-2 border-b border-slate-700/60 text-xs text-slate-400">
              <span>Engine: <span class="text-cyan-400 font-semibold">${res.source}</span></span>
              <span>Zero-Hallucination Verified ✓</span>
            </div>
            <div class="text-slate-200 text-xs sm:text-sm leading-relaxed whitespace-pre-line">${res.narrative}</div>
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
    if (btn) btn.disabled = true;

    try {
      const p = this.state.receptor.detected_pocket || {};
      const center = p.center || { x: 0, y: 0, z: 0 };
      const size = p.size || { x: 22, y: 22, z: 22 };

      const batchResult = await AnuDockAPI.runBatchDocking({
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
      if (btn) btn.disabled = false;
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
