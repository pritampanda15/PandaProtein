"""
PandaProtein CLI launcher.

Allows end users who install via:
    pip install git+https://github.com/pritampanda15/PandaProtein.git
to start the full web workbench and auto-open their browser with:
    pandaprotein

Supports optional cluster/SSH & embed ergonomics:
    --token             Require random access token (?t=...) for localhost/SSH security
    --base-path /prefix Serve behind a reverse proxy under a URL path prefix
    --idle-exit MIN     Auto-stop after MIN minutes with no active browser tab
    --reuse             Open a new tab on an already-running local PandaProtein server
    --pdb 6GGD          Auto-fetch and open a PDB entry on launch
"""

from __future__ import annotations

import argparse
import json
import os
import secrets
import socket
import threading
import time
import urllib.parse
import urllib.request
import webbrowser
from pathlib import Path

import uvicorn

__version__ = "2.0.0"

_SERVER_INFO_PATH = Path.home() / ".pandaprotein" / "server.json"


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


def _check_running_instance() -> dict | None:
    """Return running server metadata if ~/.pandaprotein/server.json points to a live instance."""
    if not _SERVER_INFO_PATH.exists():
        return None
    try:
        info = json.loads(_SERVER_INFO_PATH.read_text(encoding="utf-8"))
        base = info.get("base_url", "").rstrip("/")
        if not base:
            return None
        with urllib.request.urlopen(f"{base}/api/health", timeout=0.6) as resp:
            if resp.status == 200:
                return info
    except Exception:
        return None
    return None


def _write_server_info(base_url: str, token: str, port: int) -> None:
    try:
        _SERVER_INFO_PATH.parent.mkdir(parents=True, exist_ok=True)
        _SERVER_INFO_PATH.write_text(
            json.dumps({"base_url": base_url, "token": token, "port": port, "pid": os.getpid()}),
            encoding="utf-8",
        )
        os.chmod(_SERVER_INFO_PATH, 0o600)
    except Exception:
        pass


def _open_browser_when_ready(health_url: str, open_url: str, delay_timeout: float = 15.0) -> None:
    """Poll /api/health until Uvicorn is listening, then launch the user's browser."""
    deadline = time.time() + delay_timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(health_url, timeout=0.6) as resp:
                if resp.status == 200:
                    webbrowser.open(open_url)
                    return
        except Exception:
            time.sleep(0.2)


def _start_idle_watchdog(idle_minutes: float) -> None:
    """Shut down the server cleanly if no HTTP request/heartbeat arrives for idle_minutes."""
    def _watch() -> None:
        from . import main as main_mod
        limit_sec = max(60.0, idle_minutes * 60.0)
        while True:
            time.sleep(15.0)
            if (time.time() - main_mod._LAST_ACTIVITY_TS) > limit_sec:
                print(f"\n  [PandaProtein] No active browser for {idle_minutes:g} min — exiting cleanly.")
                os._exit(0)

    threading.Thread(target=_watch, daemon=True).start()


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
        "--token",
        action="store_true",
        help="Generate a random access token (?t=...) for shared-machine / SSH security",
    )
    parser.add_argument(
        "--base-path",
        default=os.environ.get("PANDAPROTEIN_BASE_PATH", ""),
        help="Reverse-proxy URL prefix (e.g. /pandaprotein)",
    )
    parser.add_argument(
        "--idle-exit",
        type=float,
        default=0.0,
        help="Auto-stop server after N minutes of browser inactivity (0 = disabled)",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="Reuse an already-running local PandaProtein instance and open a new tab",
    )
    parser.add_argument(
        "--pdb",
        default="",
        help="Optional RCSB PDB ID to auto-fetch on browser open (e.g. --pdb 6GGD)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"PandaProtein {__version__}",
    )
    args = parser.parse_args()

    # Optional single-instance reuse
    if args.reuse:
        existing = _check_running_instance()
        if existing:
            params: dict[str, str] = {}
            if existing.get("token"):
                params["t"] = existing["token"]
            if args.pdb:
                params["pdb"] = args.pdb.strip().upper()
            q = f"?{urllib.parse.urlencode(params)}" if params else ""
            reuse_url = f"{existing['base_url']}/{q}"
            print(f"  Reusing running PandaProtein instance: {reuse_url}")
            if not args.no_browser:
                webbrowser.open(reuse_url)
            return

    token = secrets.token_urlsafe(18) if args.token else (os.environ.get("PANDAPROTEIN_TOKEN") or "")
    if token:
        os.environ["PANDAPROTEIN_TOKEN"] = token
    if args.base_path:
        os.environ["PANDAPROTEIN_BASE_PATH"] = "/" + args.base_path.strip("/")
    if args.idle_exit > 0:
        os.environ["PANDAPROTEIN_IDLE_EXIT"] = str(args.idle_exit)

    from .engines.export_engine import check_external_tools
    from .main import app

    port = _find_free_port(args.host, args.port)
    display_host = "127.0.0.1" if args.host in ("0.0.0.0", "::") else args.host
    prefix = os.environ.get("PANDAPROTEIN_BASE_PATH", "")
    base_url = f"http://{display_host}:{port}{prefix}"

    params = {}
    if token:
        params["t"] = token
    if args.pdb:
        params["pdb"] = args.pdb.strip().upper()
    query_str = f"?{urllib.parse.urlencode(params)}" if params else ""
    open_url = f"{base_url}/{query_str}"
    health_url = f"http://{display_host}:{port}/api/health"

    _write_server_info(base_url, token, port)

    tools = check_external_tools()
    gmx_status = "detected" if tools.get("gromacs") else "not installed (optional)"
    ob_status = "detected" if tools.get("openbabel") else "not installed (optional)"

    print("")
    print(f"  PandaProtein v{__version__}")
    print("  Deep Protein Preparation · PDBFixer · PROPKA 3 · COOT · ADapt Prep · MOLE/HOLE")
    print(f"  URL       : {open_url}")
    print(f"  GROMACS   : {gmx_status}")
    print(f"  OpenBabel : {ob_status}")
    if token:
        print("  Security  : Token-protected (?t=...)")
    if args.idle_exit > 0:
        print(f"  Idle Exit : {args.idle_exit:g} min")
        _start_idle_watchdog(args.idle_exit)
    print("  Press Ctrl+C to stop.")
    print("")

    if not args.no_browser:
        threading.Thread(target=_open_browser_when_ready, args=(health_url, open_url), daemon=True).start()

    uvicorn.run(app, host=args.host, port=port, log_level="warning")


if __name__ == "__main__":
    main()
