# 16S Pipeline -- 16S rRNA Microbiome Analysis Platform

A web-based tool for processing, managing, and visualizing 16S rRNA amplicon sequencing data. Built with Plotly Dash + FastAPI + SQLite.

**Supported platforms**: Docker (Windows/macOS/Linux)

**Supported input**: Illumina paired-end or single-end amplicon FASTQ files targeting specific 16S variable regions (V1-V2, V3-V4, V4, V4-V5, V5-V6). Full-length 16S long reads (PacBio HiFi, Nanopore) are also supported -- auto-detected at upload, processed with DADA2 using platform-appropriate error models.

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
- **SRA Download** -- Fetch public datasets from NCBI SRA by accession
- **Analysis Report** -- Generate comprehensive PDF reports with all analysis results

---

## Table of Contents

- [Installation](#installation)
- [Project Structure](#project-structure)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Installation

Everything is packaged in a single Docker image -- no conda, R, or system libraries needed. See the [Tutorial](TUTORIAL.md) for a step-by-step guide with example data.

### Requirements

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows/macOS) or Docker Engine (Linux)
- **RAM**: 8 GB minimum, 16 GB recommended (PICRUSt2 needs 11 GB+). Docker Desktop defaults to only ~4 GB — increase it in Docker Desktop → **Settings** → **Resources** → **Memory**
- **Disk**: ~15 GB for the Docker image
- **Note**: PICRUSt2-related features are not available on Apple Silicon (ARM64)

### Windows

1. Download and install [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/). It will enable WSL2 automatically if needed. Restart your PC when prompted.

2. Open Docker Desktop and wait until it shows **"Docker Desktop is running"** (green icon in the system tray).

3. Open **PowerShell** and run:

```powershell
mkdir 16s-pipeline
cd 16s-pipeline
curl.exe -O https://raw.githubusercontent.com/tatsu1207/16S-Pipeline/main/docker-compose.yml
docker compose up -d
```

4. Open <a href="http://localhost:8016" target="_blank">http://localhost:8016</a> in your browser.

> The first run downloads the image (~15 GB), which may take 5-15 minutes depending on your internet speed. From the second time, you can either run `docker compose up -d` again or open Docker Desktop and start the container from the **Containers** tab -- it will start instantly.

### macOS

1. Download and install [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/) (supports both Apple Silicon and Intel).

2. Open **Terminal** and run:

```bash
mkdir 16s-pipeline && cd 16s-pipeline
curl -O https://raw.githubusercontent.com/tatsu1207/16S-Pipeline/main/docker-compose.yml
docker compose up -d
```

3. Open <a href="http://localhost:8016" target="_blank">http://localhost:8016</a> in your browser.

> The first run downloads the image, which may take 5-15 minutes. From the second time, you can either run `docker compose up -d` again or start the container from Docker Desktop.

> **Note**: PICRUSt2-related analysis (functional prediction, pathway analysis, KEGG maps) is not available on Apple Silicon (ARM64) due to lack of native support.

### Linux

1. Install [Docker Engine](https://docs.docker.com/engine/install/) if not already installed.

2. Run:

```bash
mkdir 16s-pipeline && cd 16s-pipeline
curl -O https://raw.githubusercontent.com/tatsu1207/16S-Pipeline/main/docker-compose.yml
docker compose up -d
```

3. Open <a href="http://localhost:8016" target="_blank">http://localhost:8016</a> in your browser.

> The first run downloads the image, which may take 5-15 minutes. Subsequent starts with `docker compose up -d` are instant.

### Managing the container

```bash
docker compose up -d         # Start
docker compose down          # Stop
docker compose logs -f       # View logs

# Use a different port
PORT=9000 docker compose up -d
```

### Where is my data stored?

All data (uploads, pipeline outputs, database) is stored in a Docker volume called `pipeline-data`. Your data persists across container restarts and updates.

```bash
# Back up your data
docker compose down
docker run --rm -v pipeline-data:/data -v $(pwd):/backup alpine tar czf /backup/16s-backup.tar.gz /data

# View volume info
docker volume inspect 16s-pipeline_pipeline-data
```

---

## Project Structure

```
16S-Pipeline/
├── app/
│   ├── main.py                  # FastAPI + Dash entry point
│   ├── config.py                # Auto-generated settings and paths
│   ├── api/                     # FastAPI REST endpoints
│   │   ├── pipeline.py          # Pipeline control API
│   │   └── upload.py            # File upload API
│   ├── pipeline/                # Pipeline Engine
│   │   ├── runner.py            # Pipeline orchestrator
│   │   ├── detect.py            # Auto-detect sequencing type + variable region
│   │   ├── quality.py           # Quality profiling + auto trunc_len detection
│   │   ├── qc.py                # FastQC quality control
│   │   ├── qc_pdf.py            # QC report PDF generation
│   │   ├── trim.py              # Cutadapt adapter trimming
│   │   ├── dada2.py             # DADA2 denoising (R wrapper)
│   │   ├── taxonomy.py          # Taxonomic assignment (R wrapper)
│   │   ├── phylogeny.py         # Phylogenetic tree building
│   │   ├── biom_convert.py      # BIOM format conversion
│   │   └── picrust2.py          # PICRUSt2 functional prediction
│   ├── data_manager/            # Data Management
│   │   ├── biom_ops.py          # BIOM region detection, extraction, combining
│   │   ├── mothur_convert.py    # Bidirectional BIOM/MOTHUR conversion
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
│   ├── report/                  # Report Generation
│   │   ├── report_generator.py  # PDF report builder (matplotlib + fpdf2)
│   │   └── methods_text.py      # Auto-generate Materials & Methods text
│   ├── sra/                     # SRA Integration
│   │   └── downloader.py        # NCBI SRA download via prefetch + fasterq-dump
│   ├── dashboard/               # Plotly Dash UI
│   │   ├── app.py               # Dash app initialization
│   │   ├── layout.py            # Sidebar nav + page routing
│   │   └── pages/               # One file per page
│   ├── utils/                   # Utility modules
│   │   ├── file_handler.py      # Register local FASTQ files into DB
│   │   └── metadata_parser.py   # Metadata CSV/TSV parser + validator
│   └── db/                      # SQLAlchemy models + database
│       ├── database.py          # Session management
│       └── models.py            # ORM tables
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
│   ├── sra_cache/               # Cached SRA downloads
│   ├── references/              # SILVA databases + E. coli reference
│   ├── combined/                # Combined/merged datasets
│   └── exports/                 # User exports
├── Dockerfile                   # Docker image build
├── docker-compose.yml           # One-command Docker deployment
└── docker-entrypoint.sh         # Docker container startup script
```

---

## Troubleshooting

### Container exits immediately

```bash
# Check the logs for error messages
docker compose logs

# Ensure Docker Desktop is running (Windows/macOS)
# Ensure you have enough RAM allocated to Docker
```

On Windows, Docker Desktop defaults to using half your system RAM. To increase it: Docker Desktop > Settings > Resources > Memory.

### Port 8016 already in use

```bash
# Use a different port (e.g., 9016)
PORT=9016 docker compose up -d
# Then open http://localhost:9016
```

### How to reset everything

```bash
docker compose down
docker volume rm 16s-pipeline_pipeline-data
docker compose up -d
```

This deletes all uploaded data, pipeline outputs, and the database.

### PICRUSt2 fails (exit code 137 / OOM)

PICRUSt2 requires ~11 GB RAM. Increase Docker's memory allocation: Docker Desktop > **Settings** > **Resources** > **Memory** (set to 16 GB).

> The main pipeline (FastQC, Cutadapt, DADA2, taxonomy, phylogeny) works fine with 8 GB. Only PICRUSt2 requires 16 GB. If PICRUSt2 fails, the rest of the pipeline still completes -- you can re-run PICRUSt2 later from the Pipeline Status page.

---

## Citation

If you use 16S Pipeline in your research, please cite:

> Unno T. (2026). 16S-Pipeline: A comprehensive web-based platform for end-to-end 16S rRNA amplicon sequencing analysis. *Journal of Microbiology*. Epub ahead of print. https://doi.org/10.71150/jm.2603014

```bibtex
@article{unno2026sixteenspipeline,
  title={16S-Pipeline: A comprehensive web-based platform for end-to-end 16S rRNA amplicon sequencing analysis},
  author={Unno, Tatsuya},
  journal={Journal of Microbiology},
  year={2026},
  doi={10.71150/jm.2603014},
  note={Epub ahead of print}
}
```

---

## License

MIT License -- see [LICENSE](LICENSE) for details.
