#!/usr/bin/env bash
# ============================================================================
# 16S Pipeline — Docker entrypoint
# Initializes data directories and SILVA references on first run,
# then starts the uvicorn server.
# ============================================================================
set -e

# Activate conda
eval "$(conda shell.bash hook)"

# ── Ensure data directories exist (for fresh volume mounts) ──────────────
mkdir -p /app/data/uploads /app/data/datasets /app/data/combined \
    /app/data/exports /app/data/picrust2_runs /app/data/kegg_cache \
    /app/data/sra_cache /app/data/references

# ── Copy SILVA references from image into data volume if missing ─────────
for ref in silva_nr99_v138.1_train_set.fa.gz silva_species_assignment_v138.1.fa.gz ecoli_16S.fasta; do
    if [ ! -f "/app/data/references/${ref}" ] && [ -f "/opt/silva/${ref}" ]; then
        echo "Copying ${ref} to data volume..."
        cp "/opt/silva/${ref}" "/app/data/references/${ref}"
    fi
done

PORT="${PORT:-8016}"

echo "============================================="
echo "  16S Pipeline"
echo "  http://localhost:${PORT}"
echo "============================================="

exec conda run -n microbiome_16S --no-capture-output \
    uvicorn app.main:app --host 0.0.0.0 --port "${PORT}"
