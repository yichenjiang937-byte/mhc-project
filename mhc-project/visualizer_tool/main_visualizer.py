"""Launcher for the Streamlit visualiser.

Starts a local Streamlit server on a loopback port and opens the default
browser once the server is up. Intended both for ``python main_visualizer.py``
and for the PyInstaller build produced by ``build_exe.bat``.
"""

from __future__ import annotations

import os
import sys
import threading
import webbrowser

import streamlit.web.bootstrap

DEFAULT_PORT = 8501
HOST = "127.0.0.1"


def get_script_dir() -> str:
    """Directory holding the bundled app, including inside a frozen build."""
    if getattr(sys, "frozen", False):
        return getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def main() -> int:
    port = int(os.environ.get("MHC_VISUALIZER_PORT", DEFAULT_PORT))

    os.environ.setdefault("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
    os.environ.setdefault("STREAMLIT_SERVER_HEADLESS", "true")
    os.environ["STREAMLIT_SERVER_PORT"] = str(port)
    os.environ["STREAMLIT_SERVER_ADDRESS"] = HOST

    script_path = os.path.join(get_script_dir(), "app_visualization.py")
    if not os.path.exists(script_path):
        print(f"app_visualization.py not found at {script_path}", file=sys.stderr)
        return 1

    threading.Timer(1.5, lambda: webbrowser.open(f"http://{HOST}:{port}")).start()
    streamlit.web.bootstrap.run(script_path, False, [], {})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
