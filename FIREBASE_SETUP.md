# Firebase Setup & Free Tier Deployment Guide for Bindora

This guide explains step-by-step how to attach Firebase's free tier, host Bindora securely, and protect your API keys (`sk-or-v1-...`) so they **never leak** to GitHub or public clients.

---

## 1. How Bindora Architecture Fits Firebase Free Tier

- **Frontend (Static Web App & 3D WebGL Viewer):**
  - Hosted directly on **Firebase Hosting** (100% Free on Spark Plan: 10 GB storage, 360 MB/day transfer, custom domain + SSL included).
- **Backend Compute (AutoDock Vina & RDKit):**
  - AutoDock Vina requires native C++ binary compute and Python dependencies.
  - Recommended free deployment: Run the Python backend locally, or host on a free Python/container tier (e.g. Render Free Tier, Railway, or Google Cloud Run with free monthly quota).
- **AI Explainer (DeepSeek via OpenRouter):**
  - Handled on backend with automatic local disk caching to prevent unnecessary API calls and conserve your balance.

---

## 2. Step-by-Step: Firebase Free Hosting Setup

### Step 1: Install Firebase CLI
If you don't already have the Firebase CLI installed:

```bash
npm install -g firebase-tools
```

### Step 2: Log into Firebase
```bash
firebase login
```
This opens your browser to authorize your Google account.

### Step 3: Create a Free Firebase Project
1. Go to [Firebase Console](https://console.firebase.google.com/).
2. Click **Add project** (or **Create a project**).
3. Project name enter karein: `bindora-app` (or any unique name).
4. Google Analytics enable/disable kar sakte hain (optional).
5. Default **Spark Plan (Free / $0/month)** select rehta hai. Click **Create Project**.

### Step 4: Link Project in Terminal
In the `Bindora` project directory:
```bash
firebase use --add
```
Select the project you just created in the list and set the alias to `default`.

### Step 5: Deploy Frontend to Firebase Hosting
```bash
firebase deploy --only hosting
```
Once completed, Firebase will output your live URL:
```
✔  Deploy complete!
Project Console: https://console.firebase.google.com/project/bindora-app/overview
Hosting URL: https://bindora-app.web.app
```

---

## 3. Protecting API Keys & Preventing Leaks

Your OpenRouter API Key (`sk-or-v1-********************************`):

### Why `.gitignore` is Critical
1. We have added `.env` to `.gitignore`. **NEVER** run `git add .env` or remove `.env` from `.gitignore`.
2. Anyone viewing your public GitHub repository (`https://github.com/Swelo-ui/Bindora`) will only see `.env.example` with dummy placeholders, ensuring your key remains confidential.
3. For cloud deployments, store the key in your hosting provider's **Environment Variables** dashboard (`OPENROUTER_API_KEY`), not in code.

### How Bindora Prevents Extra API Calls
Bindora features built-in **deterministic SHA-256 caching**:
- When you run a docking analysis for `Aspirin vs COX-2`, the AI narrative is generated once using `deepseek/deepseek-v4-flash-0731`.
- The result is stored in `data/cache/`.
- If you or anyone re-opens or re-clicks that result, Bindora immediately serves the cached explanation with `cached: true` without sending another request to OpenRouter!
- The API is invoked **only** when a new drug-target pair is docked.
- If the API key is missing or offline, the system automatically uses its built-in rule-based reasoning engine at **$0 recurring cost**.
