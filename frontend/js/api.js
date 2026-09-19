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

  static async getValidationReport() {
    return this.request("/api/validation-report");
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

  static async analyzeInteractions(receptorPdb, posePdbqt, smiles = "") {
    return this.request("/api/docking/interactions", {
      method: "POST",
      body: JSON.stringify({
        receptor_pdb: receptorPdb,
        pose_pdbqt: posePdbqt,
        smiles: smiles
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

  static async startBatchDocking(batchParams) {
    return this.request("/api/batch/start", {
      method: "POST",
      body: JSON.stringify(batchParams)
    });
  }

  static async getBatchStatus(jobId) {
    return this.request(`/api/batch/status/${jobId}`);
  }

  static async cancelBatchDocking(jobId) {
    return this.request(`/api/batch/cancel/${jobId}`, {
      method: "POST"
    });
  }

  static async redockValidate(params) {
    return this.request("/api/docking/redock-validate", {
      method: "POST",
      body: JSON.stringify(params)
    });
  }

  static async analyzeInteractions(receptorPdb, posePdbqt, smiles = "") {
    return this.request("/api/docking/analyze-interactions", {
      method: "POST",
      body: JSON.stringify({
        receptor_pdb: receptorPdb,
        pose_pdbqt: posePdbqt,
        smiles: smiles
      })
    });
  }

  static async getInteractionDiagram(smiles, interactions) {
    return this.request("/api/docking/interaction-diagram", {
      method: "POST",
      body: JSON.stringify({ smiles, interactions })
    });
  }

  static async getReproducibilityVersions() {
    return this.request("/api/reproducibility/versions");
  }

  static async getSimilarCompounds(smiles, threshold = 85, maxRecords = 5) {
    return this.request(`/api/ligand/similar?smiles=${encodeURIComponent(smiles)}&threshold=${threshold}&max=${maxRecords}`);
  }

  static async getEnsembleStructures(accession, limit = 8) {
    return this.request(`/api/ensemble/structures?accession=${encodeURIComponent(accession)}&limit=${limit}`);
  }

  static async runEnsembleDocking(params) {
    return this.request("/api/ensemble/run", {
      method: "POST",
      body: JSON.stringify(params)
    });
  }

  static async getPharmacophoreActives(target = "", chemblId = "", pdbId = "", uniprotAcc = "", maxActives = 10) {
    let url = `/api/pharmacophore/actives?max_actives=${maxActives}`;
    if (target) url += `&target=${encodeURIComponent(target)}`;
    if (chemblId) url += `&chembl_id=${encodeURIComponent(chemblId)}`;
    if (pdbId) url += `&pdb_id=${encodeURIComponent(pdbId)}`;
    if (uniprotAcc) url += `&uniprot_acc=${encodeURIComponent(uniprotAcc)}`;
    return this.request(url);
  }

  static async screenPharmacophore(candidates, consensusProfile) {
    return this.request("/api/pharmacophore/screen", {
      method: "POST",
      body: JSON.stringify({
        candidates: candidates,
        consensus_profile: consensusProfile
      })
    });
  }
}

// Primary export
window.BindoraAPI = BindoraAPI;
