"""
Bindora Dock — Native Desktop Launcher (Windows)
Runs Flask backend via Waitress WSGI on a dynamic loopback port
and displays the UI seamlessly inside a native OS WebView2 window.
"""

import sys
import os
import time
import socket
import atexit
import logging
import threading
import urllib.request
import multiprocessing
from pathlib import Path

# Safe stdout/stderr redirection for windowed (noconsole) mode
appdata = os.environ.get("APPDATA")
log_dir = Path(appdata) / "Bindora" if appdata else Path.home() / ".bindora"
try:
    log_dir.mkdir(parents=True, exist_ok=True)
    launcher_log_path = log_dir / "desktop_launcher.log"
    log_file = open(launcher_log_path, "a", encoding="utf-8", buffering=1)
    if sys.stdout is None:
        sys.stdout = log_file
    if sys.stderr is None:
        sys.stderr = log_file
except Exception:
    launcher_log_path = None

# Configure logging
log_handlers = []
if launcher_log_path:
    try:
        log_handlers.append(logging.FileHandler(launcher_log_path, encoding="utf-8"))
    except Exception:
        pass
if sys.stdout and not getattr(sys, "frozen", False):
    log_handlers.append(logging.StreamHandler(sys.stdout))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=log_handlers or [logging.NullHandler()]
)
logger = logging.getLogger("BindoraLauncher")

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.critical("Uncaught launcher exception:", exc_info=(exc_type, exc_value, exc_traceback))

sys.excepthook = handle_exception

# Ensure root directory is on sys.path
BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

# Configure OpenMM plugins directory for frozen environment
if getattr(sys, "frozen", False):
    bundle_candidate = Path(getattr(sys, "_MEIPASS", Path(sys.executable).resolve().parent))
    if not (bundle_candidate / "OpenMM.libs").exists() and (bundle_candidate / "_internal" / "OpenMM.libs").exists():
        bundle_dir = bundle_candidate / "_internal"
    else:
        bundle_dir = bundle_candidate

    openmm_plugins = bundle_dir / "OpenMM.libs" / "lib" / "plugins"
    if openmm_plugins.exists():
        os.environ["OPENMM_PLUGIN_DIR"] = str(openmm_plugins)
    openmm_lib = bundle_dir / "OpenMM.libs" / "lib"
    if openmm_lib.exists() and hasattr(os, "add_dll_directory"):
        try:
            os.add_dll_directory(str(openmm_lib))
        except Exception:
            pass

import webview
from backend.services.docking import terminate_active_docking_processes

def free_port() -> int:
    """Find a guaranteed-free loopback port."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]

def start_backend(port: int):
    """Start the production-grade Waitress WSGI server."""
    try:
        from waitress import serve
        from backend.app import app
        logger.info(f"Starting Waitress WSGI server on 127.0.0.1:{port}")
        serve(app, host="127.0.0.1", port=port, threads=6, channel_timeout=180, _quiet=True)
    except Exception as e:
        logger.critical(f"Waitress server encountered a fatal error: {e}", exc_info=True)

def wait_for_backend(port: int, timeout: float = 25.0) -> bool:
    """Poll the backend health endpoint until it is accepting requests."""
    start = time.time()
    url = f"http://127.0.0.1:{port}/api/health"
    while time.time() - start < timeout:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "BindoraDesktop/1.0"})
            with urllib.request.urlopen(req, timeout=1.0) as resp:
                if resp.status == 200:
                    return True
        except Exception:
            time.sleep(0.15)
    return False

def cleanup():
    """Ensure in-flight Vina calculations are terminated on exit to prevent orphaned processes."""
    try:
        killed = terminate_active_docking_processes()
        if killed > 0:
            logger.info(f"Terminated {killed} in-flight docking process(es).")
    except Exception as e:
        logger.warning(f"Error during shutdown cleanup: {e}")

atexit.register(cleanup)

def main():
    multiprocessing.freeze_support()
    port = free_port()
    
    # 1. Start backend in daemon thread
    backend_thread = threading.Thread(target=start_backend, args=(port,), daemon=True)
    backend_thread.start()

    # 2. Wait until backend is fully initialized
    is_ready = wait_for_backend(port, timeout=25.0)
    if not is_ready:
        logger.error(f"Backend server failed to respond on port {port} within timeout period.")
    else:
        logger.info(f"Backend server is healthy and accepting connections on 127.0.0.1:{port}.")

    # 3. Create native desktop window
    app_url = f"http://127.0.0.1:{port}"
    if is_ready:
        window = webview.create_window(
            title="Bindora Dock — Molecular Docking & PK/PD Suite",
            url=app_url,
            width=1440,
            height=900,
            min_size=(1100, 700),
            resizable=True,
            text_select=True,
            background_color="#0b0f19"
        )
    else:
        err_html = f"""<!DOCTYPE html>
        <html style="background:#090a0f;color:#e2e8f0;font-family:sans-serif;display:flex;align-items:center;justify-content:center;height:100vh;margin:0;">
        <div style="text-align:center;max-width:550px;padding:32px;background:#13151f;border:1px solid #282b3d;border-radius:12px;">
            <h2 style="color:#f87171;margin-bottom:12px;">Engine Initialization Timeout</h2>
            <p style="color:#94a3b8;font-size:14px;line-height:1.6;">The internal Bindora backend service did not respond within 25 seconds on port {port}.</p>
            <p style="color:#64748b;font-size:12px;margin-top:16px;">Log file: <code style="color:#38bdf8;">{launcher_log_path}</code></p>
            <button onclick="window.location.reload()" style="margin-top:20px;padding:10px 20px;background:#2563eb;color:white;border:none;border-radius:8px;cursor:pointer;font-weight:600;">Retry Connection</button>
        </div>
        </html>"""
        window = webview.create_window(
            title="Bindora Dock — Connection Notice",
            html=err_html,
            width=800,
            height=500,
            resizable=False,
            background_color="#0b0f19"
        )

    if window:
        window.events.closed += cleanup

    # 4. Start GUI event loop
    logger.info("Launching native WebView2 desktop window...")
    webview.start()

if __name__ == "__main__":
    main()
