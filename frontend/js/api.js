// Client API module for Bindora backend communications

const API_BASE = ""; // Relative path to match same-origin or proxied API

class BindoraAPI {
  static async request(endpoint, options = {}) {
    const defaultHeaders = {
      "Content-Type": "application/json"
    };

    const config = {
      ...options,
      headers: {
        ...defaultHeaders,
        ...options.headers
      }
    };

    try {
      const response = await fetch(`${API_BASE}${endpoint}`, config);
      const data = await response.json();
      if (!response.ok) {
        throw new Error(data.error || `HTTP error! status: ${response.status}`);
      }
      return data;
    } catch (error) {
      console.error(`[API Error] ${endpoint}:`, error);
      throw error;
    }
  }

  static async checkHealth() {
    return this.request("/api/health");
  }

  static async getBenchmarks() {
    return this.request("/api/benchmarks");
  }

  static async searchPubChem(query) {
    return this.request(`/api/search/pubchem?query=${encodeURIComponent(query)}`);
  }

  static async searchPubChemByCID(cid) {
    return this.request(`/api/search/pubchem/cid/${encodeURIComponent(cid)}`);
  }

  static async searchRCSB(query) {
    return this.request(`/api/search/rcsb?query=${encodeURIComponent(query)}`);
  }

  static async prepareLigand(structure, isSdf = false) {
    return this.request("/api/structure/ligand", {
      method: "POST",
      body: JSON.stringify({ structure, is_sdf: isSdf })
    });
  }

  static async prepareReceptor(pdbContent = "", pdbId = "", targetChain = null) {
    return this.request("/api/structure/receptor", {
      method: "POST",
      body: JSON.stringify({
        pdb_content: pdbContent,
        pdb_id: pdbId,
        target_chain: targetChain
      })
    });
  }

  static async runDocking(params) {
    return this.request("/api/docking/run", {
      method: "POST",
      body: JSON.stringify(params)
    });
  }

  static async analyzeInteractions(receptorPdb, posePdbqt) {
    return this.request("/api/docking/interactions", {
      method: "POST",
      body: JSON.stringify({
        receptor_pdb: receptorPdb,
        pose_pdbqt: posePdbqt
      })
    });
  }

  static async crosscheckBioactivity(drug, target) {
    return this.request(`/api/pkpd/crosscheck?drug=${encodeURIComponent(drug)}&target=${encodeURIComponent(target)}`);
  }

  static async explainNarrative(reportData, provider = "auto") {
    return this.request("/api/narrative/explain", {
      method: "POST",
      body: JSON.stringify({
        ...reportData,
        provider: provider
      })
    });
  }

  static async runBatchDocking(batchParams) {
    return this.request("/api/batch/run", {
      method: "POST",
      body: JSON.stringify(batchParams)
    });
  }

  static async redockValidate(params) {
    return this.request("/api/docking/redock-validate", {
      method: "POST",
      body: JSON.stringify(params)
    });
  }

  static async getReproducibilityVersions() {
    return this.request("/api/reproducibility/versions");
  }
}

// Primary export + backward compatibility alias
window.BindoraAPI = BindoraAPI;
window.AnuDockAPI = BindoraAPI;
