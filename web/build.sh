#!/usr/bin/env bash
# Build the browser-only version of EZplot (Shinylive / WebAssembly) into site/.
# Serve locally with:  python -m http.server -d site 8000
set -euo pipefail
cd "$(dirname "$0")/.."
stage=web/_app
rm -rf "$stage" site
mkdir -p "$stage"
rsync -a --exclude __pycache__ ezplot "$stage/"
cat > "$stage/app.py" <<'PY'
from ezplot.app import app  # noqa: F401
PY
# shinylive only scans app.py for imports, so list the package's dependencies;
# those outside Pyodide (adjustText) are fetched from PyPI by micropip
printf '%s\n' pandas numpy scipy matplotlib pillow adjustText > "$stage/requirements.txt"
uv run --frozen shinylive export "$stage" site
rm -rf "$stage"
echo "Built site/ ($(du -sh site | cut -f1))"
