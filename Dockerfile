# ============================================================================
# 16S Pipeline — Dockerfile
# Builds all 5 conda environments + SILVA references into a single image.
# Image size: ~10-15 GB (5 conda envs + reference databases)
#
# Build:   docker build -t 16s-pipeline .
# Run:     docker compose up -d
# Open:    http://localhost:8016
# ============================================================================

FROM condaforge/mambaforge:latest

LABEL maintainer="16S Pipeline"
LABEL description="End-to-end 16S rRNA microbiome analysis platform"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# ── System dependencies ─────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    wget \
    procps \
    lsof \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ── Conda environment 1: microbiome_16S (Python + CLI tools) ───────────────
RUN mamba create -n microbiome_16S -c conda-forge python=3.11 -y && \
    mamba clean -afy

# Install bioinformatics CLI tools
RUN mamba install -n microbiome_16S --override-channels -c conda-forge -c bioconda \
    fastqc cutadapt mafft fasttree bbmap sra-tools vsearch -y && \
    mamba clean -afy

# Install Python packages via pip
RUN conda run -n microbiome_16S pip install --no-cache-dir \
    fastapi \
    "uvicorn[standard]" \
    dash \
    dash-bootstrap-components \
    plotly \
    sqlalchemy \
    pandas \
    numpy \
    scipy \
    scikit-bio \
    python-multipart \
    biom-format \
    matplotlib \
    fpdf2 \
    statsmodels \
    dash-uploader \
    matplotlib-venn \
    scikit-learn

# ── Conda environment 2: dada2_16S (R + DADA2) ────────────────────────────
RUN mamba create -n dada2_16S --override-channels -c conda-forge -c bioconda \
    bioconductor-dada2 r-optparse r-jsonlite -y && \
    mamba clean -afy

# ── Conda environment 3: analysis_16S (R + DA tools) ──────────────────────
RUN mamba create -n analysis_16S --override-channels -c conda-forge -c bioconda \
    bioconductor-phyloseq bioconductor-ancombc bioconductor-deseq2 \
    bioconductor-aldex2 bioconductor-microbiome r-optparse r-jsonlite -y && \
    mamba clean -afy

# conda-forge's r-directlabels (needed by ALDEx2) is "noarch" but ships an
# x86_64 directlabels.so, so ALDEx2 can't load on arm64. Rebuild it from CRAN
# for the target architecture.
RUN conda run -n analysis_16S Rscript -e \
    "install.packages('directlabels', repos='https://cloud.r-project.org', INSTALL_opts='--no-lock')"

# ── Conda environment 4: maaslin2_16S (MaAsLin2 + vegan + LinDA) ─────────
# modeest (needed by LinDA) comes from CRAN below: conda-forge has no
# linux-aarch64 build of its dependency r-stable. rmutil and fBasics are its
# compiled dependencies that conda-forge does build for both architectures.
RUN mamba create -n maaslin2_16S --override-channels -c conda-forge -c bioconda \
    bioconductor-maaslin2 r-optparse r-jsonlite \
    r-remotes r-ggrepel r-lme4 r-foreach r-rmutil r-fbasics -y && \
    mamba clean -afy

# Install vegan and modeest from CRAN
RUN conda run -n maaslin2_16S Rscript -e \
    "install.packages(c('vegan', 'modeest'), repos='https://cloud.r-project.org', INSTALL_opts='--no-lock', Ncpus=4)"

# Install LinDA from a pinned GitHub archive (v0.2.0, its latest commit,
# 2023-12-16). A plain archive download, unlike remotes::install_github, isn't
# subject to GitHub's API rate limit, which failed an arm64 build. Its
# dependencies come from conda above; a failed install is only a warning in R,
# so check that LinDA actually loads, and retry in case of a network hiccup.
ARG LINDA_URL=https://github.com/zhouhj1994/LinDA/archive/af0f62fad83f25a0df272da9ccbec7db49149166.tar.gz
RUN for attempt in 1 2 3; do \
        conda run -n maaslin2_16S Rscript -e \
            "install.packages('${LINDA_URL}', repos=NULL, type='source', INSTALL_opts='--no-lock')" && \
        conda run -n maaslin2_16S Rscript -e "library(LinDA)" && break; \
        echo "LinDA install attempt $attempt failed"; sleep 20; \
    done; \
    conda run -n maaslin2_16S Rscript -e "library(LinDA)"

# Fail the build if any differential abundance package cannot be loaded
RUN conda run -n analysis_16S Rscript -e " \
        for (p in c('ANCOMBC', 'microbiome', 'DESeq2', 'ALDEx2', 'phyloseq')) \
            if (!requireNamespace(p, quietly=TRUE)) stop('missing R package: ', p)" && \
    conda run -n maaslin2_16S Rscript -e " \
        for (p in c('Maaslin2', 'LinDA', 'vegan', 'modeest')) \
            if (!requireNamespace(p, quietly=TRUE)) stop('missing R package: ', p)"

# ── Conda environment 5: picrust2_16S ───────────────────────────────────────
# Required on amd64: fail the build rather than publish an image without it.
# Its post-link script downloads from GitHub, which can fail transiently, so
# retry a few times first. arm64 has no bioconda package, so skip it there
# (the app handles missing PICRUSt2 gracefully).
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "amd64" ]; then \
        for attempt in 1 2 3; do \
            mamba create -n picrust2_16S --override-channels -c conda-forge -c bioconda \
                picrust2 -y && break; \
            echo "PICRUSt2 install attempt $attempt failed"; \
            rm -rf /opt/conda/envs/picrust2_16S; sleep 30; \
        done; \
        mamba clean -afy; \
        conda run -n picrust2_16S picrust2_pipeline.py --version \
        || { echo "ERROR: PICRUSt2 install failed on amd64"; exit 1; }; \
    else \
        echo "Skipping PICRUSt2 on $TARGETARCH (no bioconda package available)"; \
    fi

# ── Download SILVA 138.1 references to a staging location ─────────────────
# Stored in /opt/silva so they can be copied into the data volume at first run
RUN mkdir -p /opt/silva && \
    wget -q -O /opt/silva/silva_nr99_v138.1_train_set.fa.gz \
        "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz" && \
    wget -q -O /opt/silva/silva_species_assignment_v138.1.fa.gz \
        "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz" && \
    wget -q -O /opt/silva/silva_nr99_v138.1_wSpecies_train_set.fa.gz \
        "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz"

# ── Copy application code ─────────────────────────────────────────────────
COPY app/ /app/app/
COPY r_scripts/ /app/r_scripts/
COPY data/references/ecoli_16S.fasta /opt/silva/ecoli_16S.fasta
COPY docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# ── Environment variables ─────────────────────────────────────────────────
ENV PORT=8016
ENV CONDA_BASE=/opt/conda
# Store database inside the data volume so it persists
ENV DATABASE_PATH=/app/data/microbiome.db

EXPOSE 8016

ENTRYPOINT ["/app/docker-entrypoint.sh"]
