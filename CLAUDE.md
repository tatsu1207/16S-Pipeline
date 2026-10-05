# 16S Pipeline — Claude Reference Guide

Read this file before working on this project. It explains how the web tool works end-to-end.

## What This Tool Does

An end-to-end microbiome analysis platform: upload raw 16S rRNA FASTQ files (or fetch them from NCBI SRA) → denoise with DADA2 → assign taxonomy → build phylogenetic tree → run functional prediction → perform statistical analysis → generate a PDF report. All through a browser UI.

**Supported inputs**: Illumina paired-end/single-end (V1-V2, V3-V4, V4, V4-V5, V5-V6 regions), full-length 16S long reads (PacBio HiFi, Nanopore).

## Architecture

- **Frontend**: Plotly Dash (Python interactive dashboard)
- **Backend API**: FastAPI (REST endpoints)
- **Database**: SQLite via SQLAlchemy ORM (`microbiome.db` at project root natively; `/app/data/microbiome.db` in Docker via `DATABASE_PATH`)
- **R Integration**: Subprocess calls to conda environments (`conda run -n <env> Rscript ...`)
- **Pipeline execution**: Background Python threads with subprocess spawning
- **Server**: Uvicorn. Port = `$PORT` if set (Docker: 8016), otherwise 7000 + user UID (native)

### 5 Conda Environments

| Environment | Purpose |
|---|---|
| `microbiome_16S` | Web app (Python 3.11), FastQC, Cutadapt, MAFFT, FastTree, BBMap, vsearch, sra-tools |
| `dada2_16S` | R + DADA2 for denoising and taxonomy |
| `analysis_16S` | R + ALDEx2, DESeq2, ANCOM-BC2 (+ phyloseq, microbiome) |
| `maaslin2_16S` | R + MaAsLin2, LinDA, vegan (separate env due to R version conflicts) |
| `picrust2_16S` | PICRUSt2 functional prediction (optional; Linux x86_64 only) |

The `Dockerfile` is the source of truth for env contents. `setup_native.sh` mirrors it — keep the two package lists in sync. `config.py` falls back to `analysis_16S` for MaAsLin2/LinDA/NMDS if `maaslin2_16S` doesn't exist under `$CONDA_BASE/envs`.

## Project Structure

```
app/
├── main.py                    # FastAPI + Dash entry point, eagerly imports all pages
├── config.py                  # Paths, conda env names, CPU/thread limits, defaults, port
├── api/
│   ├── upload.py              # POST /api/upload (no read counting/region detection — UI doesn't use it)
│   └── pipeline.py            # POST/GET /api/pipeline/* — launch/cancel/status
├── pipeline/
│   ├── runner.py              # Main orchestrator, _run_pipeline()
│   ├── detect.py              # SE/PE, variable region, platform detection
│   ├── quality.py             # Auto-truncation parameter detection (Q20 sliding window)
│   ├── dada2.py               # DADA2 subprocess wrapper
│   ├── taxonomy.py            # SILVA taxonomy (+ long-read species: NB classifier, vsearch fallback)
│   ├── phylogeny.py           # MAFFT alignment + FastTree
│   ├── biom_convert.py        # ASV table → BIOM (HDF5) conversion
│   ├── qc.py / qc_pdf.py      # FastQC wrapper, QC report
│   ├── trim.py                # Cutadapt primer trimming
│   └── picrust2.py            # PICRUSt2 wrapper
├── sra/
│   ├── downloader.py          # prefetch + fasterq-dump by accession
│   ├── metadata_fetcher.py    # BioSample metadata from NCBI
│   └── submission.py          # SRA submission spreadsheet generation
├── report/
│   ├── report_generator.py    # Analysis report PDF
│   └── methods_text.py        # Auto-generated Materials & Methods paragraph
├── dashboard/
│   ├── app.py                 # Dash app initialization
│   ├── layout.py              # Sidebar navigation + URL routing (unknown paths → "Coming Soon")
│   └── pages/                 # Each page registers its own callbacks
│       ├── intro_page.py              # /
│       ├── file_manager.py            # /files — upload FASTQ, attach metadata, delete
│       ├── sra_download_page.py       # /sra-download
│       ├── sra_submit_page.py         # /sra-submit
│       ├── pipeline_status.py         # /pipeline — launch + monitor DADA2 pipeline
│       ├── biom_browser_page.py       # /biom-browser
│       ├── rare_asv_page.py           # /rare-asv — low-prevalence ASV removal
│       ├── subsampling_page.py        # /subsample — rarefaction + sample filtering
│       ├── combine_page.py            # /combine — merge BIOM files
│       ├── datasets_page.py           # /datasets — V-region extraction
│       ├── sample_tree_page.py        # /sample-tree — outlier detection
│       ├── alpha_page.py              # /alpha
│       ├── beta_page.py               # /beta — distances, PCoA, NMDS, PERMANOVA
│       ├── taxonomy_page.py           # /taxonomy — stacked bar plots
│       ├── diff_abundance_page.py     # /diff-abundance — 5-tool DA comparison
│       ├── picrust2_page.py           # /picrust2 — standalone PICRUSt2 runner
│       ├── pathways_page.py           # /pathways — pathway DA comparison
│       ├── kegg_map_page.py           # /kegg-map
│       ├── report_page.py             # /report
│       └── mothur_page.py             # not routed or imported (MOTHUR page was removed)
├── analysis/
│   ├── shared.py              # Helpers: BIOM↔DataFrame, metadata loading
│   ├── alpha.py, beta.py, taxonomy.py
│   ├── diff_abundance.py      # Dispatcher for 5 DA tools (calls R scripts)
│   ├── pathways.py, pathway_plots.py, kegg_aggregation.py, kegg_map.py
│   └── r_runner.py            # R subprocess wrapper (streams output, parses JSON result)
├── data_manager/
│   ├── biom_ops.py            # Region detection/extraction, dataset combining
│   ├── subsample.py           # Rarefaction to uniform depth
│   ├── rare_asv.py            # Filter ASVs by prevalence/abundance
│   └── mothur_convert.py      # unused since the MOTHUR page was removed
├── db/
│   ├── database.py            # Session management, init_db()
│   └── models.py              # 12 ORM tables
└── utils/
    ├── file_handler.py        # register_upload() — not called by the app
    └── metadata_parser.py     # CSV/TSV metadata parsing

r_scripts/                     # run_dada2, run_taxonomy, run_aldex2, run_ancombc,
                               # run_deseq2, run_linda, run_maaslin2, run_nmds (.R)

data/                          # All gitignored except .gitkeep files
├── uploads/                   # Raw FASTQ uploads
├── datasets/                  # Pipeline outputs (per dataset_id)
├── combined/                  # Merged datasets
├── picrust2_runs/             # PICRUSt2 outputs
├── references/                # SILVA 138.1 (+ wSpecies train set) + E. coli ref
├── sra_cache/                 # SRA downloads (can grow very large; safe to clear when idle)
├── kegg_cache/                # Cached KEGG API responses (24h TTL)
└── exports/                   # User-downloaded files

test_samples/                  # Tutorial FASTQs + metadata.tsv (tracked in git)
```

Removed features (don't re-add references): Random Forest, Correlation, Network, Association, Longitudinal analysis, MOTHUR conversion page.

## Pipeline Execution Flow

When a user launches a pipeline, these steps run sequentially in a background thread:

1. **FastQC** — Quality reports → `{dataset_dir}/qc/`
2. **Cutadapt** — Trim primers → `{dataset_dir}/trimmed/`
3. **Auto-truncation** — Analyze quality scores, find where 10bp sliding window mean < Q20. For PE: ensures `trunc_f + trunc_r >= insert_len + min_overlap`
4. **DADA2** — Denoise reads, remove chimeras → `asv_table.tsv`, `rep_seqs.fasta`, `track_reads.tsv`. Uses platform-specific error models for long reads.
5. **Taxonomy** — SILVA 138.1 assignment → `taxonomy.tsv`. Long reads: exact match, then NB classifier with the species-level train set, then vsearch fallback (≥ 99% identity).
6. **Phylogenetic tree** — MAFFT + FastTree → `tree.nwk`
7. **BIOM conversion** — Combine ASV table + taxonomy into HDF5 BIOM → `asv_table.biom`
8. **PICRUSt2** (optional) — Functional prediction in separate conda env

### Status Tracking

- Progress stored in `{dataset_dir}/status.json` with `current_step`, `progress_pct`, `steps_completed`, `pid`
- UI polls `GET /api/pipeline/status/{dataset_id}` via `dcc.Interval`
- Cancellation: sets `threading.Event` flag + kills subprocess via `os.killpg()`
- PID persists in `status.json` for recovery after server restart

## Data Flow Summary

```
Upload FASTQ (or SRA download) → detect SE/PE + region + platform → attach metadata (CSV/TSV)
→ launch pipeline → DADA2 denoising → taxonomy → tree → BIOM
→ optional: filter ASVs, rarefy, combine datasets, extract V-regions, outlier detection
→ analysis: alpha/beta diversity, taxonomy plots, differential abundance, pathways
→ export: PNG/SVG/PDF plots, CSV tables, PDF report with methods text
```

## Database (SQLite + SQLAlchemy)

Key tables in `app/db/models.py`:
- **projects** — Study grouping
- **uploads** — FASTQ upload batches (sequencing_type, variable_region, platform, study)
- **fastq_files** — Individual files (sample_name, read_direction, read_count, avg_read_length)
- **upload_metadata** — Key-value per (upload, sample)
- **datasets** — Pipeline outputs (status: pending/processing/complete/failed)
- **dataset_fastq_files** — M2M linking datasets to FASTQ files
- **samples** / **sample_metadata** — Post-pipeline sample data
- **qc_metrics**, **dataset_combinations**, **analysis_results**, **picrust2_runs**

The same sample name can exist in several uploads — identify rows by `(upload_id, sample_name)`, never by name alone (the File Manager's checkboxes use `"{upload_id}:{sample_name}"` keys).

Session pattern:
```python
from app.db.database import get_session
with get_session() as db:
    dataset = db.query(Dataset).filter(...).first()
```

## Key Code Patterns

### Dash Callbacks
```python
@dash_app.callback(
    Output("id-output", "children"),
    Input("id-button", "n_clicks"),
    [State("id-store", "data")],
)
def callback(n_clicks, store_data):
    if not n_clicks:
        return no_update
    return html.Div(...)
```

All page modules are eagerly imported in `main.py` so callbacks register before the server starts. Each page has a `get_layout()` function (`pipeline_status.py` uses a module-level `layout`). URL routing is a callback in `layout.py`; adding or removing a page means touching `main.py`, the sidebar and the route in `layout.py`, plus README/TUTORIAL.

### R Script Execution
```python
# In analysis/r_runner.py — spawns R in appropriate conda env
conda run -n {env_name} Rscript {script}.R --arg1 val1 --arg2 val2
```
Streams stdout/stderr to logger, parses JSON from last stdout line as result. Script → env mapping is `_R_SCRIPT_ENV_MAP` in `config.py`.

### DADA2 R Script (`r_scripts/run_dada2.R`)
- All output uses `log_msg()` (= `cat()` + `flush.console()`) to ensure output is flushed before potential crashes.
- Long-running steps (dereplication, denoising, merging, chimera removal) log elapsed time on completion.
- Auto-skips `filterAndTrim()` if `filtered/` directory already contains output files (also supported via explicit `--skip_filter` flag).
- The Python wrapper (`app/pipeline/dada2.py`) decodes signal names on crash (e.g., SIGSEGV, SIGKILL) and includes the last R output line in the error message.

### Pipeline Threading
- Global dicts: `_running_pipelines`, `_cancel_events`, `_active_procs`
- Each dataset_id gets its own cancel event
- Daemon threads auto-cleanup on shutdown

## Important Gotchas

1. **Truncation params (0, 0)** = auto-detect. The auto-detection in `quality.py` is critical for paired-end merge success.
2. **Long-read detection** (`detect.py`): median read length > 1000 bp. Median quality ≥ Q25 → PacBio HiFi, else Nanopore. Long reads use different DADA2 error models (PacBioErrfun vs loessErrfun).
3. **PICRUSt2 needs 11GB+ RAM** and has no ARM64 package. Can fail independently without blocking the main pipeline.
4. **Primer detection threshold**: >30% of reads must contain a known primer for it to be detected.
5. **BIOM format is HDF5** (binary). Requires biom-format library. Rep seqs are embedded as observation metadata.
6. **KEGG API has rate limits** — responses cached 24h in `data/kegg_cache/`.
7. **Taxonomy can be None** — several places must handle missing taxonomy gracefully (past bug source).
8. **Multi-region combining** has two modes: by-sequence (same region) or by-taxonomy (cross-region, uses E. coli alignment positions).
9. **Callback import order matters** — all page modules must be imported before Dash starts.
10. **File uploads** use dash-uploader (chunked) rather than standard dcc.Upload (which blocks the browser for large files).
11. **Pipeline crash recovery**: If the R subprocess dies (segfault, OOM, signal), the daemon thread also dies silently. The UI detects this via PID liveness check and marks the dataset as "failed". On re-run, DADA2 automatically skips filtering if `filtered/` files already exist.
12. **Host PATH can mask missing tools**: `conda run` also sees `~/bin` and `/usr/bin`, so a tool missing from an env may still work natively but fail in Docker (this hid the missing vsearch). `setup_native.sh --check` looks inside each env's own `bin/`.
13. **Native and Docker data are separate**: native uses `./data` + `./microbiome.db`; the container uses its `pipeline-data` volume.

## Config Essentials (`app/config.py`)

- `DATA_DIR` = `{PROJECT_DIR}/data` (not overridable); `DATABASE_PATH` env overrides the DB location
- `SILVA_TRAIN_SET` / `SILVA_SPECIES` / `SILVA_SPECIES_TRAIN_SET` = reference databases in `data/references/`
- `CONDA_BASE` env (default `~/miniforge3`) — used to detect `maaslin2_16S`
- `CPU_COUNT`, `MAX_THREADS` = all cores − 1 (no fixed cap). Pipeline page suggests 2 threads per sample up to `MAX_THREADS`.
- `R_DEFAULT_THREADS` = min(32, `MAX_THREADS`) — default (not cap) for R DA tools, which start one R worker per thread
- `picrust2_default_threads()` = RAM-aware (~2 GB per worker), capped at `MAX_THREADS`
- `DADA2_DEFAULTS` = `{trim_left_f: 0, trim_left_r: 0, trunc_len_f: 0, trunc_len_r: 0, min_overlap: 12, threads: MAX_THREADS}`
- `LONGREAD_DADA2_DEFAULTS` = `{min_len: 1000, max_len: 1600, max_ee: 10, band_size: 32}`
- `conda_cmd(args, env_name)` = wrapper that builds `["conda", "run", "-n", env_name, "--no-capture-output"] + args`

## Running the App

End users run Docker only (see README). For development, run natively — no image rebuild needed.

```bash
# Native development (this checkout)
bash setup_native.sh                 # build missing envs, download SILVA, verify everything
bash setup_native.sh --check         # verify only
bash setup_native.sh --recreate ENV  # rebuild one env (or: --recreate all)
bash run_native.sh --reload          # http://localhost:{7000+UID}; restarts on edits in app/
                                     # (--reload kills running pipelines on every save)

# Docker (what users run)
docker compose up -d                 # http://localhost:8016, image ghcr.io/tatsu1207/16s-pipeline:latest
docker compose pull && docker compose up -d   # update to the latest published image
```

### Releasing a new Docker image

`.github/workflows/docker-publish.yml` builds and pushes the image to GHCR **only when a GitHub release is published or the workflow is run manually** — pushing to `main` does not rebuild the image.
