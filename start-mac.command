#!/usr/bin/env bash
# EZplot launcher for macOS: double-click this file (first time: right-click > Open)
# First run installs uv (a small Python manager, in your home folder, no admin rights)
# and the app's packages (~2 minutes). Later runs start in a few seconds.
set -e
cd "$(dirname "$0")"
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
if ! command -v uv >/dev/null 2>&1; then
  echo "First run: installing uv (one time only)..."
  if command -v curl >/dev/null 2>&1; then
    curl -LsSf https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 sh
  else
    wget -qO- https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 sh
  fi
fi
echo "Starting EZplot (the first start downloads Python and packages, please wait)..."
uv run --frozen --no-dev ezplot "$@" || { echo; echo "EZplot stopped with an error (see above)."; read -r -p "Press Enter to close..."; }
