# MicrobiomeDash -- 16S rRNA Microbiome Analysis Dashboard

A web-based tool for processing, managing, and visualizing 16S rRNA amplicon sequencing data. Built with Plotly Dash + FastAPI + SQLite.

**Three integrated tools:**

| Tool | Purpose |
|------|---------|
| **Pipeline Engine** | FASTQ.gz -> DADA2 -> Taxonomy -> Phylogeny -> BIOM -> PICRUSt2 |
| **Data Manager** | Browse, download, combine, subsample datasets across studies |
| **Analysis Dashboard** | Alpha/beta diversity, differential abundance, pathway analysis, KEGG maps |

**Analysis capabilities:**

- **Alpha diversity** -- Shannon, Simpson, observed OTUs with Kruskal-Wallis / Mann-Whitney tests
- **Beta diversity** -- Bray-Curtis / Jaccard distance, PCoA, NMDS, PERMANOVA (pairwise + global)
- **Taxonomy** -- Stacked bar plots at any taxonomic level
- **Differential abundance** -- 5 tools: ALDEx2, DESeq2, ANCOM-BC2, LinDA, MaAsLin2; all-pairwise mode; volcano plots
- **Pathway analysis** -- PICRUSt2 output analysis with multi-tool DA, KO-to-KEGG aggregation, errorbar/heatmap/PCA plots (ggpicrust2-inspired)
- **KEGG Map** -- Targeted pathway inspection with DA-colored KEGG maps

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the App](#running-the-app)
- [Project Structure](#project-structure)
- [Windows WSL Setup](#windows-wsl-setup)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Prerequisites

- **Operating System**: Ubuntu/Debian Linux (native or WSL2 on Windows)
- **Conda**: Miniconda or Miniforge installed ([install guide](https://docs.conda.io/en/latest/miniconda.html))
- **RAM**: 8 GB minimum, 16 GB recommended (16 GB required for PICRUSt2)
- **Disk Space**: ~10 GB for software + reference databases
- **System libraries** (install if missing):

```bash
sudo apt install -y build-essential curl git wget libxml2-dev \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libtiff5-dev libjpeg-dev libpng-dev
```

---

## Installation

### Step 1: Clone the repository

```bash
git clone https://github.com/tatsu1207/microbiome-dashboard.git
cd microbiome-dashboard
```

### Step 2: Run the setup script

The setup script automatically:
- Creates the `microbiome` Conda environment with Python 3.11 + R 4.3
- Installs all Python packages (FastAPI, Dash, scikit-bio, biom-format, etc.)
- Installs R packages (DADA2, ALDEx2, DESeq2, ANCOM-BC2, MaAsLin2, LinDA, vegan)
- Installs bioinformatics CLI tools (FastQC, Cutadapt, MAFFT, FastTree, vsearch)
- Creates a separate `picrust2` Conda environment for PICRUSt2
- Downloads SILVA 138.1 reference databases (optional, prompted)
- Generates `app/config.py` with auto-detected paths

```bash
chmod +x setup_ubuntu.sh
./setup_ubuntu.sh
```

> Expected time: 15-30 minutes depending on internet speed and system.
> The R/Bioconductor package installation takes the longest.

### Step 3: Activate the environment

```bash
conda activate microbiome
```

---

## Running the App

```bash
conda activate microbiome
./run.sh
```

The app runs in the background. The port is auto-assigned based on your UID (7000 + UID). Open the URL shown in the terminal output.

To run manually in the foreground:

```bash
conda activate microbiome
uvicorn app.main:app --reload --reload-exclude data --host 0.0.0.0 --port 8050
```

---

## Project Structure

```
microbiome-dashboard/
├── app/
│   ├── main.py                  # FastAPI + Dash entry point
│   ├── config.py                # Auto-generated settings and paths
│   ├── api/                     # FastAPI REST endpoints
│   │   ├── pipeline.py          # Pipeline control API
│   │   └── upload.py            # File upload API
│   ├── pipeline/                # Pipeline Engine
│   │   ├── runner.py            # Pipeline orchestrator
│   │   ├── qc.py                # FastQC quality control
│   │   ├── trim.py              # Cutadapt adapter trimming
│   │   ├── dada2.py             # DADA2 denoising
│   │   ├── taxonomy.py          # Taxonomic assignment
│   │   ├── phylogeny.py         # Phylogenetic tree building
│   │   ├── biom_convert.py      # BIOM format conversion
│   │   └── picrust2.py          # PICRUSt2 functional prediction
│   ├── data_manager/            # Data Management
│   │   ├── biom_ops.py          # BIOM file operations
│   │   ├── mothur_convert.py    # Mothur format import
│   │   ├── rare_asv.py          # Rare ASV filtering
│   │   └── subsample.py         # Rarefaction subsampling
│   ├── analysis/                # Analysis Engine
│   │   ├── shared.py            # Shared BIOM/metadata helpers
│   │   ├── r_runner.py          # R subprocess wrapper
│   │   ├── alpha.py             # Alpha diversity (skbio + scipy)
│   │   ├── beta.py              # Beta diversity, PCoA, NMDS, PERMANOVA
│   │   ├── taxonomy.py          # Taxonomy aggregation
│   │   ├── diff_abundance.py    # Multi-tool DA dispatcher + volcano plots
│   │   ├── pathways.py          # PICRUSt2 pathway DA
│   │   ├── kegg_aggregation.py  # KO-to-KEGG aggregation + annotation
│   │   ├── kegg_map.py          # KEGG pathway map helpers
│   │   └── pathway_plots.py     # Errorbar, heatmap, PCA visualizations
│   ├── dashboard/               # Plotly Dash UI
│   │   ├── app.py               # Dash app initialization
│   │   ├── layout.py            # Sidebar nav + page routing
│   │   └── pages/               # One file per page
│   └── db/                      # SQLAlchemy models + database
│       ├── database.py          # Session management
│       └── models.py            # 11 ORM tables
├── r_scripts/                   # R analysis scripts
│   ├── run_dada2.R              # DADA2 pipeline
│   ├── run_taxonomy.R           # Taxonomy assignment
│   ├── run_nmds.R               # NMDS ordination (vegan)
│   ├── run_aldex2.R             # ALDEx2 DA
│   ├── run_deseq2.R             # DESeq2 DA
│   ├── run_ancombc.R            # ANCOM-BC2 DA
│   ├── run_linda.R              # LinDA DA
│   └── run_maaslin2.R           # MaAsLin2 DA
├── data/                        # Data storage (gitignored except placeholders)
│   ├── uploads/                 # User FASTQ uploads
│   ├── datasets/                # Processed pipeline outputs
│   ├── picrust2_runs/           # PICRUSt2 output directories
│   ├── kegg_cache/              # Cached KEGG API data (24h TTL)
│   ├── references/              # SILVA databases + E. coli reference
│   ├── combined/                # Combined/merged datasets
│   └── exports/                 # User exports
├── setup_ubuntu.sh              # One-command installation script
├── run.sh                       # Start the application
├── environment.yml              # Conda environment specification
├── requirements.txt             # Python dependencies (pip)
└── Makefile                     # Development commands
```

---

## Download Reference Databases

The setup script will offer to download SILVA automatically. If you skipped that step:

```bash
cd data/references

# SILVA 138.1 training set for DADA2 (~24 MB)
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz

# SILVA 138.1 species assignment (~77 MB)
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

cd ../..
```

The E. coli 16S reference (`ecoli_16S.fasta`) is included in the repository and used for primer orientation detection.

---

## Windows WSL Setup

If you're on Windows, you need WSL2 with Ubuntu first. If you already have it, skip to [Installation](#installation).

### Step 1: Enable WSL

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This installs WSL2 with Ubuntu by default. If WSL is already installed but you need Ubuntu:

```powershell
wsl --install -d Ubuntu-24.04
```

Restart your computer when prompted.

### Step 2: Initial Ubuntu Setup

After restart, Ubuntu will open automatically (or search for "Ubuntu" in the Start menu). It will ask you to create a username and password.

### Step 3: Install system dependencies and Conda

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential curl git wget libxml2-dev \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libtiff5-dev libjpeg-dev libpng-dev

# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

### Tips for WSL

- **Access WSL files from Windows**: Open File Explorer and go to `\\wsl$\Ubuntu\home\<your-username>`
- **Access Windows files from WSL**: Your C: drive is at `/mnt/c/`
- **VS Code integration**: Install the [Remote - WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) extension

---

## Troubleshooting

### "conda: command not found"

```bash
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

### R package installation fails (DADA2)

Make sure system libraries are installed (see [Prerequisites](#prerequisites)), then retry:

```bash
conda activate microbiome
Rscript -e 'BiocManager::install("dada2", force=TRUE)'
```

### "Permission denied" on setup_ubuntu.sh

```bash
chmod +x setup_ubuntu.sh
```

### Port already in use

```bash
# Find and kill the process using the port
lsof -ti:8050 | xargs kill -9

# Or use a different port
uvicorn app.main:app --reload --reload-exclude data --host 0.0.0.0 --port 8051
```

### PICRUSt2 installation fails

PICRUSt2 requires its own environment due to dependency conflicts:

```bash
conda create -n picrust2 -c bioconda -c conda-forge picrust2 -y
```

The setup script handles this automatically.

### WSL runs out of memory / PICRUSt2 OOM killed (exit code 137)

PICRUSt2's phylogenetic placement requires ~11 GB RAM. Create or edit `C:\Users\<YourUsername>\.wslconfig`:

```ini
[wsl2]
memory=16GB
swap=8GB
```

Then restart WSL from PowerShell: `wsl --shutdown`

> The pipeline (FastQC, Cutadapt, DADA2, taxonomy, phylogeny) works fine with 8 GB. Only PICRUSt2 requires 16 GB. If PICRUSt2 fails, the rest of the pipeline still completes -- you can re-run PICRUSt2 later from the Pipeline Status page.

---

## License

MIT License -- see [LICENSE](LICENSE) for details.
