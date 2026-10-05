#!/usr/bin/env bash
# ============================================================================
# 16S Pipeline — start the app natively (after bash setup_native.sh)
#
# Usage:
#   bash run_native.sh              # http://localhost:<7000 + UID>
#   bash run_native.sh --reload     # restart automatically when app/ changes
#   PORT=9000 bash run_native.sh    # pick another port
#
# Uses ./data and ./microbiome.db in this checkout — separate from the Docker
# container's data volume. --reload kills running pipelines on every save.
# ============================================================================
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

if ! command -v conda >/dev/null 2>&1; then
    for base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
        [ -x "$base/bin/conda" ] && export PATH="$base/bin:$PATH" && break
    done
fi

# config.py locates envs (e.g. maaslin2_16S) under $CONDA_BASE/envs
export CONDA_BASE="$(conda info --base)"
export PORT="${PORT:-$((7000 + $(id -u)))}"

RELOAD=()
if [ "${1:-}" = "--reload" ]; then
    RELOAD=(--reload --reload-dir app)
fi

echo "16S Pipeline (native) → http://localhost:${PORT}"
exec conda run -n microbiome_16S --no-capture-output \
    uvicorn app.main:app --host 0.0.0.0 --port "$PORT" "${RELOAD[@]+"${RELOAD[@]}"}"
