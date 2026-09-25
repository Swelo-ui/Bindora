"""
Bindora Dock — Native Desktop Launcher (Windows)
Runs Flask backend via Waitress WSGI on a dynamic loopback port
and displays the UI seamlessly inside a native OS WebView2 window.
"""

import sys
import time
import socket
import atexit
import logging
import threading
import urllib.request
from pathlib import Path

# Ensure root directory is on sys.path
BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

# Configure OpenMM plugins directory for frozen environment
if getattr(sys, "frozen", False):
    import os
    bundle_dir = Path(getattr(sys, "_MEIPASS", Path(sys.executable).resolve().parent))
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

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("BindoraLauncher")

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

def wait_for_backend(port: int, timeout: float = 15.0) -> bool:
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
    port = free_port()
    
    # 1. Start backend in daemon thread
    backend_thread = threading.Thread(target=start_backend, args=(port,), daemon=True)
    backend_thread.start()

    # 2. Wait until backend is fully initialized
    if not wait_for_backend(port, timeout=15.0):
        logger.error("Backend server failed to respond within timeout period.")
    else:
        logger.info("Backend server is healthy and accepting connections.")

    # 3. Create native desktop window
    app_url = f"http://127.0.0.1:{port}"
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

    if window:
        window.events.closed += cleanup

    # 4. Start GUI event loop
    logger.info("Launching native WebView2 desktop window...")
    webview.start()

if __name__ == "__main__":
    main()
