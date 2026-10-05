#!/usr/bin/env bash
# ============================================================================
# 16S Pipeline — native (non-Docker) setup for development
#
# Builds the same 5 conda environments as the Dockerfile, downloads the SILVA
# references into data/references/, and verifies every tool and package the
# app calls. Safe to re-run: existing envs are kept and only checked.
#
# End users should use Docker (see README). This script is for editing the
# code and seeing changes without rebuilding the image.
#
# Usage:
#   bash setup_native.sh                  # create missing envs, verify all
#   bash setup_native.sh --check          # verify only, install nothing
#   bash setup_native.sh --recreate ENV   # delete and rebuild one env
#                                         # (repeatable, or --recreate all)
#
# Then start the app with:  bash run_native.sh
#
# Keep the package lists below in sync with the Dockerfile.
# ============================================================================
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

ALL_ENVS=(microbiome_16S dada2_16S analysis_16S maaslin2_16S picrust2_16S)
CHECK_ONLY=0
RECREATE=()

while [ $# -gt 0 ]; do
    case "$1" in
        --check) CHECK_ONLY=1 ;;
        --recreate)
            [ $# -ge 2 ] || { echo "--recreate needs an env name or 'all'"; exit 1; }
            if [ "$2" = "all" ]; then RECREATE=("${ALL_ENVS[@]}"); else RECREATE+=("$2"); fi
            shift ;;
        -h|--help) sed -n '2,21p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

info() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33mWARNING:\033[0m %s\n' "$*"; }
ok()   { printf '  \033[32mok\033[0m    %s\n' "$*"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$*"; FAILED=$((FAILED + 1)); }
FAILED=0

# ── Find conda / mamba ─────────────────────────────────────────────────────
if ! command -v conda >/dev/null 2>&1; then
    for base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
        [ -x "$base/bin/conda" ] && export PATH="$base/bin:$PATH" && break
    done
fi
command -v conda >/dev/null 2>&1 || {
    echo "conda not found. Install Miniforge first: https://github.com/conda-forge/miniforge"
    exit 1
}
if command -v mamba >/dev/null 2>&1; then MAMBA=mamba; else MAMBA=conda; fi
CONDA_BASE="$(conda info --base)"
info "Using $MAMBA ($CONDA_BASE)"

OS="$(uname -s)"
ARCH="$(uname -m)"
if [ "$OS" != "Linux" ]; then
    warn "Native setup is only tested on Linux. On $OS/$ARCH some bioconda packages may be unavailable."
fi

# Required env: stop the script if creation fails
require_env() { local rc=0; create_env "$@" || rc=$?; [ "$rc" -ne 2 ] || { echo "Failed to create $1"; exit 1; }; }

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }
should_recreate() { local e; for e in "${RECREATE[@]+"${RECREATE[@]}"}"; do [ "$e" = "$1" ] && return 0; done; return 1; }
rscript() { local env="$1"; shift; conda run -n "$env" Rscript -e "$@"; }

# Create an env only if it is missing (or was asked to be rebuilt).
# Returns 0 if the env was (re)created, 1 if it already existed, 2 on failure.
create_env() {
    local env="$1"; shift
    if env_exists "$env" && should_recreate "$env"; then
        info "Removing $env"
        conda env remove -n "$env" -y >/dev/null
    fi
    if env_exists "$env"; then
        info "$env exists — keeping it (use --recreate $env to rebuild)"
        return 1
    fi
    info "Creating $env"
    "$MAMBA" create -n "$env" "$@" -y || return 2
    return 0
}

# ── Install ────────────────────────────────────────────────────────────────
if [ "$CHECK_ONLY" -eq 0 ]; then

    # 1. microbiome_16S — Python web app + CLI tools.
    # Topped up even when it exists, so newly added packages get installed.
    require_env microbiome_16S -c conda-forge python=3.11
    info "Installing CLI tools into microbiome_16S"
    "$MAMBA" install -n microbiome_16S --override-channels -c conda-forge -c bioconda \
        fastqc cutadapt mafft fasttree bbmap sra-tools vsearch -y
    info "Installing Python packages into microbiome_16S"
    conda run -n microbiome_16S pip install --quiet \
        fastapi "uvicorn[standard]" dash dash-bootstrap-components plotly \
        sqlalchemy pandas numpy scipy scikit-bio python-multipart biom-format \
        matplotlib fpdf2 statsmodels dash-uploader matplotlib-venn scikit-learn

    # 2. dada2_16S — DADA2 + taxonomy
    require_env dada2_16S --override-channels -c conda-forge -c bioconda \
        bioconductor-dada2 r-optparse r-jsonlite

    # 3. analysis_16S — ALDEx2, DESeq2, ANCOM-BC2
    require_env analysis_16S --override-channels -c conda-forge -c bioconda \
        bioconductor-phyloseq bioconductor-ancombc bioconductor-deseq2 \
        bioconductor-aldex2 bioconductor-microbiome r-optparse r-jsonlite
    # conda-forge's "noarch" r-directlabels ships an x86_64 .so (see Dockerfile)
    if [ "$(uname -m)" != "x86_64" ]; then
        info "Rebuilding directlabels (CRAN) for $(uname -m)"
        rscript analysis_16S "install.packages('directlabels', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')"
    fi

    # 4. maaslin2_16S — MaAsLin2, LinDA, vegan
    rc=0
    create_env maaslin2_16S --override-channels -c conda-forge -c bioconda \
        bioconductor-maaslin2 r-optparse r-jsonlite \
        r-remotes r-ggrepel r-lme4 r-foreach r-rmutil r-fbasics || rc=$?
    [ "$rc" -ne 2 ] || { echo "Failed to create maaslin2_16S"; exit 1; }
    if [ "$rc" -eq 0 ]; then
        # modeest from CRAN: conda-forge has no linux-aarch64 r-stable (see Dockerfile)
        info "Installing vegan and modeest (CRAN) into maaslin2_16S"
        rscript maaslin2_16S "install.packages(c('vegan', 'modeest'), repos='https://cloud.r-project.org', INSTALL_opts='--no-lock', Ncpus=4)"
        info "Installing LinDA (GitHub) into maaslin2_16S"
        rscript maaslin2_16S "
            tryCatch(
                remotes::install_github('zhouhj1994/LinDA', upgrade='never', INSTALL_opts='--no-lock'),
                error = function(e) {
                    message('First attempt failed, retrying...'); Sys.sleep(10)
                    remotes::install_github('zhouhj1994/LinDA', upgrade='never', INSTALL_opts='--no-lock')
                })"
    fi

    # 5. picrust2_16S — optional; no bioconda package outside Linux x86_64
    if [ "$OS" = "Linux" ] && [ "$ARCH" = "x86_64" ]; then
        rc=0
        create_env picrust2_16S --override-channels -c conda-forge -c bioconda picrust2 || rc=$?
        [ "$rc" -ne 2 ] || warn "PICRUSt2 install failed — the app runs without it"
    else
        info "Skipping picrust2_16S on $OS/$ARCH (no bioconda package)"
    fi

    # ── Data directories + references ─────────────────────────────────────
    mkdir -p data/uploads data/datasets data/combined data/exports \
        data/picrust2_runs data/kegg_cache data/sra_cache data/references

    SILVA_URL="https://zenodo.org/record/4587955/files"
    for ref in silva_nr99_v138.1_train_set.fa.gz silva_species_assignment_v138.1.fa.gz \
               silva_nr99_v138.1_wSpecies_train_set.fa.gz; do
        dest="data/references/$ref"
        if [ -e "$dest" ]; then
            continue
        fi
        info "Downloading $ref"
        curl -fL --retry 3 -o "$dest.part" "$SILVA_URL/$ref"
        mv "$dest.part" "$dest"
    done
fi

# ── Verify ─────────────────────────────────────────────────────────────────
# Tools are checked inside the env's own bin/, because `conda run` also sees
# ~/bin and /usr/bin — a tool found there would be missing in Docker.
info "Verifying environments"

for env in microbiome_16S dada2_16S analysis_16S maaslin2_16S; do
    env_exists "$env" && ok "env $env" || bad "env $env is missing"
done

if env_exists microbiome_16S; then
    prefix="$(conda run -n microbiome_16S printenv CONDA_PREFIX)"
    for tool in fastqc cutadapt mafft FastTree bbduk.sh prefetch fasterq-dump vsearch uvicorn; do
        [ -x "$prefix/bin/$tool" ] && ok "$tool" || bad "$tool not in microbiome_16S"
    done
    missing="$(conda run -n microbiome_16S python -c "
import importlib
mods = ['fastapi', 'uvicorn', 'dash', 'dash_bootstrap_components', 'plotly',
        'sqlalchemy', 'pandas', 'numpy', 'scipy', 'skbio', 'multipart', 'biom',
        'matplotlib', 'fpdf', 'statsmodels', 'dash_uploader', 'matplotlib_venn',
        'sklearn']
missing = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception:
        missing.append(m)
print(' '.join(missing))" 2>/dev/null)" || missing="(python failed)"
    [ -z "$missing" ] && ok "Python packages" \
        || bad "Python packages missing: $missing (run without --check to install)"
fi

# One Rscript call per env (conda run is slow to start)
check_r() {
    local env="$1"; shift
    env_exists "$env" || return 0
    local pkgs missing
    pkgs="$(printf "'%s'," "$@")"
    missing="$(rscript "$env" "p <- c(${pkgs%,}); cat(p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)])" 2>/dev/null)" \
        || missing="(Rscript failed)"
    if [ -z "$missing" ]; then
        ok "R packages in $env"
    else
        bad "R packages missing in $env: $missing — try: bash setup_native.sh --recreate $env"
    fi
}
check_r dada2_16S    dada2 optparse jsonlite
check_r analysis_16S ANCOMBC microbiome DESeq2 ALDEx2 phyloseq optparse jsonlite
check_r maaslin2_16S Maaslin2 LinDA vegan optparse jsonlite

if env_exists picrust2_16S; then
    ok "env picrust2_16S (optional)"
else
    warn "picrust2_16S not installed — PICRUSt2, pathway and KEGG map pages will not work"
fi

for ref in silva_nr99_v138.1_train_set.fa.gz silva_species_assignment_v138.1.fa.gz \
           silva_nr99_v138.1_wSpecies_train_set.fa.gz ecoli_16S.fasta; do
    [ -e "data/references/$ref" ] && ok "$ref" || bad "data/references/$ref missing"
done

echo
if [ "$FAILED" -gt 0 ]; then
    echo "$FAILED check(s) failed."
    exit 1
fi
echo "All checks passed. Start the app with:  bash run_native.sh"
