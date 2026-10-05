"""
MicrobiomeDash — SILVA taxonomy assignment R script wrapper.
"""
import json
import logging
import subprocess
import tempfile
from collections import defaultdict
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from app.config import (
    CONDA_ENV_NAME, DADA2_ENV_NAME, LONGREAD_SPECIES_MIN_ID, R_SCRIPTS_DIR,
    SILVA_SPECIES, SILVA_SPECIES_TRAIN_SET, SILVA_TRAIN_SET, conda_cmd,
)

# Species epithets that do not name a species
_NON_SPECIES = {"sp.", "bacterium", "phage", "uncultured", "unidentified", "metagenome"}


def run_taxonomy(
    rep_seqs_path: Path,
    output_dir: Path,
    threads: int,
    logger: logging.Logger,
    proc_callback: Callable[[subprocess.Popen], None] | None = None,
    skip_species: bool = True,
    longread: bool = False,
) -> dict:
    """Run SILVA 138.1 taxonomy assignment via the R script.

    Short reads: genus-level Bayesian classification, optionally followed by
    exact-match species assignment (addSpecies).
    Long reads (``longread=True``): Bayesian classification with the SILVA
    species-level training set, overridden by a unique exact match where one
    exists, then a vsearch fallback (>= 99% identity) for ASVs still lacking
    a species.

    Returns:
        dict with: success (bool), taxonomy_path (str)
    """
    taxonomy_path = output_dir / "taxonomy.tsv"

    cmd_args = [
        "Rscript", str(R_SCRIPTS_DIR / "run_taxonomy.R"),
        "--rep_seqs", str(rep_seqs_path),
        "--output", str(taxonomy_path),
        "--silva_train", str(SILVA_TRAIN_SET),
        "--threads", str(threads),
    ]
    if skip_species:
        cmd_args.append("--skip_species")
    elif longread and SILVA_SPECIES_TRAIN_SET.exists():
        cmd_args.extend([
            "--species_train", str(SILVA_SPECIES_TRAIN_SET),
            "--silva_species", str(SILVA_SPECIES),
        ])
    else:
        if longread:
            logger.warning(
                f"{SILVA_SPECIES_TRAIN_SET.name} not found; "
                "falling back to exact-match species assignment"
            )
        cmd_args.extend(["--silva_species", str(SILVA_SPECIES)])

    cmd = conda_cmd(cmd_args, env_name=DADA2_ENV_NAME)

    logger.info("Running taxonomy assignment (SILVA 138.1)...")

    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        start_new_session=True,
    )
    if proc_callback:
        proc_callback(proc)

    last_line = ""
    for line in proc.stdout:
        stripped = line.rstrip()
        if stripped:
            logger.info(f"[Taxonomy] {stripped}")
            last_line = stripped

    proc.wait()

    if proc.returncode != 0:
        raise RuntimeError(f"Taxonomy R script failed (exit code {proc.returncode})")

    # Parse JSON status
    try:
        r_status = json.loads(last_line)
        if r_status.get("status") == "error":
            raise RuntimeError(f"Taxonomy error: {r_status.get('message', 'unknown')}")
    except json.JSONDecodeError:
        pass

    if not taxonomy_path.exists():
        raise FileNotFoundError(f"Taxonomy output not found: {taxonomy_path}")

    if longread and not skip_species:
        add_species_vsearch(rep_seqs_path, taxonomy_path, threads, logger, proc_callback)

    logger.info(f"Taxonomy assignment complete: {taxonomy_path}")

    return {"success": True, "taxonomy_path": str(taxonomy_path)}


def add_species_vsearch(
    rep_seqs_path: Path,
    taxonomy_path: Path,
    threads: int,
    logger: logging.Logger,
    proc_callback: Callable[[subprocess.Popen], None] | None = None,
    min_id: float = LONGREAD_SPECIES_MIN_ID,
) -> int:
    """Fill missing species by vsearch global alignment against SILVA.

    For each ASV with a genus but no species, the best hits (>= ``min_id``) in
    the SILVA species-assignment set are collected. A species is assigned only
    when every best hit agrees on a single species within the ASV's genus;
    ties between species are left unassigned. Updates ``taxonomy_path`` in
    place and returns the number of ASVs newly assigned.
    """
    tax = pd.read_csv(taxonomy_path, sep="\t", dtype=str, keep_default_na=False)
    tax = tax.replace({"NA": ""})
    todo = set(tax.loc[(tax["Genus"] != "") & (tax["Species"] == ""), "ASV_ID"])
    if not todo:
        return 0

    seqs = _read_fasta(rep_seqs_path)
    with tempfile.TemporaryDirectory() as tmp:
        query = Path(tmp) / "query.fasta"
        hits_path = Path(tmp) / "hits.tsv"
        query.write_text("".join(f">{a}\n{seqs[a]}\n" for a in todo if a in seqs))
        cmd = conda_cmd([
            "vsearch", "--usearch_global", str(query),
            "--db", str(SILVA_SPECIES),
            "--id", str(min_id),
            "--strand", "plus",
            "--maxaccepts", "100", "--maxrejects", "500",
            "--userout", str(hits_path),
            "--userfields", "query+target+id",
            "--notrunclabels",  # keep "Genus species" in target labels
            "--threads", str(threads),
        ], env_name=CONDA_ENV_NAME)
        logger.info(
            f"Species fallback: aligning {len(todo)} ASVs to SILVA "
            f"(vsearch, >= {min_id:.0%} identity)..."
        )
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            start_new_session=True,
        )
        if proc_callback:
            proc_callback(proc)
        output = proc.communicate()[0]
        if proc.returncode != 0:
            logger.warning(f"vsearch species fallback failed; skipping: {output[-500:]}")
            return 0
        hits = pd.read_csv(
            hits_path, sep="\t", names=["query", "target", "id"], dtype={"id": float},
        ) if hits_path.stat().st_size else pd.DataFrame(columns=["query", "target", "id"])

    genus_of = dict(zip(tax["ASV_ID"], tax["Genus"]))
    best = defaultdict(set)
    for asv, grp in hits.groupby("query"):
        top = grp[grp["id"] >= grp["id"].max() - 1e-9]
        for target in top["target"]:
            parts = target.split()
            if len(parts) < 3:
                continue
            genus, epithet = parts[1], parts[2]
            if epithet in _NON_SPECIES or genus not in _genus_names(genus_of[asv]):
                best[asv].add(None)
            else:
                best[asv].add(epithet)

    assigned = {a: s.pop() for a, s in best.items() if len(s) == 1 and None not in s}
    if assigned:
        mask = tax["ASV_ID"].isin(assigned)
        tax.loc[mask, "Species"] = tax.loc[mask, "ASV_ID"].map(assigned)
        tax.replace({"": "NA"}).to_csv(taxonomy_path, sep="\t", index=False)
    logger.info(
        f"Species fallback: {len(assigned)} of {len(todo)} ASVs assigned "
        f"({len(best) - len(assigned)} had hits but ambiguous or genus mismatch)"
    )
    return len(assigned)


def _genus_names(genus: str) -> set[str]:
    """SILVA genus labels like 'Escherichia-Shigella' cover each component."""
    return {genus, *genus.replace("/", "-").split("-")}


def _read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()
    return seqs
