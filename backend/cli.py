"""
PandaProtein CLI launcher.

Allows end users who install via:
    pip install git+https://github.com/pritampanda15/PandaProtein.git
to start the full web workbench and auto-open their browser with:
    pandaprotein
"""

from __future__ import annotations

import argparse
import socket
import threading
import time
import urllib.request
import webbrowser

import uvicorn

__version__ = "2.0.0"


def _find_free_port(host: str, preferred_port: int, max_tries: int = 20) -> int:
    """Return preferred_port if free, otherwise scan upward for the next open port."""
    for offset in range(max_tries):
        port = preferred_port + offset
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                s.bind((host, port))
                return port
            except OSError:
                continue
    return preferred_port


def _open_browser_when_ready(url: str, delay_timeout: float = 15.0) -> None:
    """Poll /api/health until Uvicorn is listening, then launch the user's browser."""
    health_url = f"{url.rstrip('/')}/api/health"
    deadline = time.time() + delay_timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(health_url, timeout=0.6) as resp:
                if resp.status == 200:
                    webbrowser.open(url)
                    return
        except Exception:
            time.sleep(0.2)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="pandaprotein",
        description="PandaProtein — Deep Protein Preparation & Crystallographic Workbench",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host address to bind (default: 127.0.0.1; use 0.0.0.0 for remote access)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to serve on (default: 8000; auto-increments if busy)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Start the server without automatically opening a browser tab",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"PandaProtein {__version__}",
    )
    args = parser.parse_args()

    from .engines.export_engine import check_external_tools
    from .main import app

    port = _find_free_port(args.host, args.port)
    display_host = "127.0.0.1" if args.host in ("0.0.0.0", "::") else args.host
    url = f"http://{display_host}:{port}"

    tools = check_external_tools()
    gmx_status = "detected" if tools.get("gromacs") else "not installed (optional)"
    ob_status = "detected" if tools.get("openbabel") else "not installed (optional)"

    print("")
    print(f"  PandaProtein v{__version__}")
    print("  Deep Protein Preparation · PDBFixer · PROPKA 3 · COOT · MOLE/HOLE")
    print(f"  URL       : {url}")
    print(f"  GROMACS   : {gmx_status}")
    print(f"  OpenBabel : {ob_status}")
    print("  Press Ctrl+C to stop.")
    print("")

    if not args.no_browser:
        threading.Thread(target=_open_browser_when_ready, args=(url,), daemon=True).start()

    uvicorn.run(app, host=args.host, port=port, log_level="warning")


if __name__ == "__main__":
    main()
