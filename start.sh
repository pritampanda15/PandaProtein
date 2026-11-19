#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON_BIN="$ROOT_DIR/.venv/bin/python"

if [ ! -x "$PYTHON_BIN" ]; then
  echo "Virtual environment not found at $ROOT_DIR/.venv"
  echo "Create it with: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
  exit 1
fi

if [ "${1:-}" = "--dev" ]; then
  echo "Starting PDBPrep Studio in DEV mode (FastAPI :8000 + Vite HMR :5173)..."
  "$PYTHON_BIN" -m uvicorn backend.main:app --reload --host 127.0.0.1 --port 8000 &
  BACK_PID=$!
  trap "kill $BACK_PID 2>/dev/null || true" EXIT
  npm --prefix "$ROOT_DIR/frontend" run dev
else
  echo "Starting PDBPrep Studio on http://127.0.0.1:8000 (auto-reload enabled) ..."
  exec "$PYTHON_BIN" -m uvicorn backend.main:app --reload --reload-dir "$ROOT_DIR/backend" --host 127.0.0.1 --port 8000
fi
