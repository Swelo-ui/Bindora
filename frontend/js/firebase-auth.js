// Firebase Authentication & Realtime Database Module for Bindora
// Project ID: bindora-1db62
// Database: https://bindora-1db62-default-rtdb.asia-southeast1.firebasedatabase.app

class BindoraFirebase {
  constructor() {
    this.auth = null;
    this.db = null;
    this.currentUser = null;
    this.isInitialized = false;

    this.defaultConfig = {
      apiKey: localStorage.getItem("bindora_firebase_api_key") || atob("QUl6YVN5QVh1dTVFWVhycFByYWtPOF9RLXJOaXdUbDZTOXZHZ3hZ"),
      authDomain: "bindora-1db62.firebaseapp.com",
      databaseURL: "https://bindora-1db62-default-rtdb.asia-southeast1.firebasedatabase.app",
      projectId: "bindora-1db62",
      storageBucket: "bindora-1db62.firebasestorage.app",
      messagingSenderId: "615327811618",
      appId: "1:615327811618:web:3d5e4eda9b24476f89c87b",
      measurementId: "G-ZN1340CMFE"
    };

    this.init();
  }

  init(customConfig = null) {
    if (!window.firebase) {
      console.warn("[Firebase] Firebase SDK not loaded.");
      return false;
    }

    const config = customConfig || this.defaultConfig;

    try {
      if (!firebase.apps.length) {
        if (config.apiKey) {
          firebase.initializeApp(config);
          console.log("[Firebase] Initialized with config for project:", config.projectId);
        } else {
          console.log("[Firebase] Awaiting Web API Key configuration.");
          return false;
        }
      }

      this.auth = firebase.auth();
      this.db = firebase.database();
      this.isInitialized = true;

      // Listen for auth state changes
      this.auth.onAuthStateChanged(user => {
        this.currentUser = user;
        console.log("[Firebase Auth] User state changed:", user ? user.email : "Logged Out");
        this.onUserChanged(user);
      });

      return true;
    } catch (err) {
      console.error("[Firebase Init Error]:", err);
      return false;
    }
  }

  setApiKey(key) {
    if (!key) return;
    this.defaultConfig.apiKey = key.trim();
    localStorage.setItem("bindora_firebase_api_key", this.defaultConfig.apiKey);
    if (!this.isInitialized) {
      this.init();
    }
  }

  // --- Authentication Methods ---

  async signInWithGoogle() {
    if (!this.ensureInitialized()) return;
    const provider = new firebase.auth.GoogleAuthProvider();
    provider.addScope("profile");
    provider.addScope("email");
    try {
      const result = await this.auth.signInWithPopup(provider);
      return result.user;
    } catch (error) {
      console.error("[Google Auth Error]:", error);
      throw error;
    }
  }

  async signInWithEmail(email, password) {
    if (!this.ensureInitialized()) return;
    try {
      const result = await this.auth.signInWithEmailAndPassword(email, password);
      return result.user;
    } catch (error) {
      console.error("[Email SignIn Error]:", error);
      throw error;
    }
  }

  async signUpWithEmail(email, password) {
    if (!this.ensureInitialized()) return;
    try {
      const result = await this.auth.createUserWithEmailAndPassword(email, password);
      return result.user;
    } catch (error) {
      console.error("[Email SignUp Error]:", error);
      throw error;
    }
  }

  async signOut() {
    if (!this.auth) return;
    await this.auth.signOut();
  }

  // --- Realtime Database Methods ---

  async saveDockingRun(runData) {
    if (!this.currentUser) {
      throw new Error("Please sign in to save your docking simulations to cloud.");
    }
    if (!this.db) {
      throw new Error("Firebase Realtime Database is not initialized.");
    }

    const userId = this.currentUser.uid;
    const timestamp = Date.now();
    const runRef = this.db.ref(`users/${userId}/docking_runs`).push();

    const payload = {
      id: runRef.key,
      savedAt: timestamp,
      userEmail: this.currentUser.email,
      drugName: runData.drugName || "Unknown Drug",
      targetName: runData.targetName || "Target Receptor",
      pdbId: runData.pdbId || "N/A",
      affinityKcal: runData.affinityKcal || 0.0,
      theoreticalKdNm: runData.theoreticalKdNm || "N/A",
      ligandEfficiency: runData.ligandEfficiency || 0.0,
      hbondCount: runData.hbondCount || 0,
      lipinskiStatus: runData.lipinskiStatus || "Pass",
      crosscheckBadge: runData.crosscheckBadge || "Computational Prediction Only",
      isCrossChecked: !!runData.isCrossChecked
    };

    await runRef.set(payload);
    console.log("[Firebase RTDB] Saved run to cloud:", runRef.key);
    return payload;
  }

  async fetchSavedRuns() {
    if (!this.currentUser || !this.db) return [];
    const userId = this.currentUser.uid;
    const snapshot = await this.db.ref(`users/${userId}/docking_runs`).once("value");
    const data = snapshot.val();
    if (!data) return [];

    return Object.values(data).sort((a, b) => (b.savedAt || 0) - (a.savedAt || 0));
  }

  async deleteSavedRun(runId) {
    if (!this.currentUser || !this.db || !runId) return;
    const userId = this.currentUser.uid;
    await this.db.ref(`users/${userId}/docking_runs/${runId}`).remove();
  }

  ensureInitialized() {
    if (!this.isInitialized) {
      // Try to initialize if apiKey is saved
      const savedKey = localStorage.getItem("bindora_firebase_api_key");
      if (savedKey) {
        this.setApiKey(savedKey);
      }
      if (!this.isInitialized) {
        // Prompt user to enter their Web API Key in settings
        const modal = document.getElementById("modal-firebase-auth");
        if (modal) modal.classList.remove("hidden");
        throw new Error("Please configure your Firebase Web API Key first in Settings.");
      }
    }
    return true;
  }

  onUserChanged(user) {
    const userBtn = document.getElementById("btn-user-auth");
    const userDisplay = document.getElementById("user-display-name");
    const userAvatar = document.getElementById("user-avatar-img");
    const saveCloudBtn = document.getElementById("btn-save-cloud");

    if (user) {
      if (userBtn) {
        userBtn.title = `Signed in as ${user.displayName || user.email}`;
      }
      if (userDisplay) {
        userDisplay.textContent = (user.displayName || user.email.split("@")[0]);
      }
      if (userAvatar) {
        userAvatar.src = user.photoURL || `https://api.dicebear.com/7.x/bottts/svg?seed=${user.uid}`;
        userAvatar.classList.remove("hidden");
      }
      const defaultIcon = document.getElementById("user-default-icon");
      if (defaultIcon) defaultIcon.classList.add("hidden");

      if (saveCloudBtn) {
        saveCloudBtn.classList.remove("hidden");
      }
    } else {
      if (userDisplay) userDisplay.textContent = "Sign In";
      if (userAvatar) userAvatar.classList.add("hidden");
      const defaultIcon = document.getElementById("user-default-icon");
      if (defaultIcon) defaultIcon.classList.remove("hidden");
      if (saveCloudBtn) saveCloudBtn.classList.add("hidden");
    }

    // Update cloud history panel if visible
    if (window.app && typeof window.app.loadCloudHistory === "function") {
      window.app.loadCloudHistory();
    }
  }
}

window.bindoraFirebase = new BindoraFirebase();
