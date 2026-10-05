"""
Methods Text Generator — Auto-generate a Materials & Methods paragraph
from pipeline parameters stored in the Dataset model.
"""
from app.db.database import SessionLocal
from app.db.models import Dataset, Sample


def generate_methods_text(dataset_id: int) -> str:
    """
    Generate a Materials & Methods paragraph for a completed dataset.
    Returns a formatted text string ready for publication.
    """
    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            return "Dataset not found."

        samples = (
            db.query(Sample)
            .filter(Sample.dataset_id == dataset_id)
            .all()
        )

        return _build_text(dataset, samples)
    finally:
        db.close()


def _build_text(dataset, samples: list) -> str:
    """Build the methods paragraph from dataset and sample data."""
    parts = []

    # Sequencing description
    region = dataset.variable_region or "16S rRNA"
    platform = _platform_name(dataset.platform)
    seq_type = dataset.sequencing_type or "paired-end"

    parts.append(
        f"Raw 16S rRNA gene amplicon sequences targeting the {region} region "
        f"were generated using {platform} {seq_type} sequencing."
    )

    # Primer trimming
    if dataset.custom_fwd_primer or dataset.custom_rev_primer:
        fwd = dataset.custom_fwd_primer or "N/A"
        rev = dataset.custom_rev_primer or "N/A"
        parts.append(
            f"Primer sequences (forward: {fwd}, reverse: {rev}) were removed "
            f"using Cutadapt (Martin, 2011)."
        )
    else:
        parts.append(
            "Region-specific primer sequences were removed using Cutadapt (Martin, 2011)."
        )

    # DADA2 processing
    dada2_text = "Quality filtering and amplicon sequence variant (ASV) inference were performed using DADA2 (Callahan et al., 2016)"

    if dataset.platform in ("pacbio", "nanopore"):
        error_model = "PacBioErrfun" if dataset.platform == "pacbio" else "loessErrfun"
        dada2_text += f" with the {error_model} error model for long-read data"
    elif seq_type == "paired-end":
        trunc_f = dataset.trunc_len_f
        trunc_r = dataset.trunc_len_r
        if trunc_f and trunc_r:
            dada2_text += (
                f". Forward reads were truncated at {trunc_f} bp "
                f"and reverse reads at {trunc_r} bp"
            )
            if dataset.min_overlap:
                dada2_text += f" with a minimum overlap of {dataset.min_overlap} bp for read merging"
    elif seq_type == "single-end":
        trunc_f = dataset.trunc_len_f
        if trunc_f:
            dada2_text += f". Reads were truncated at {trunc_f} bp"

    dada2_text += ". Chimeric sequences were removed using the consensus method."
    parts.append(dada2_text)

    # Taxonomy
    longread = dataset.platform in ("pacbio", "nanopore")
    if longread:
        parts.append(
            "Taxonomy was assigned to species level using the naive Bayesian classifier "
            "(Wang et al., 2007) against the SILVA NR99 v138.1 species-level training set "
            "(Quast et al., 2013). Where an ASV matched a single species exactly in the "
            "SILVA species-assignment database, that species took precedence. ASVs still "
            "lacking a species were aligned against the same database with VSEARCH "
            "(Rognes et al., 2016), and a species was assigned when all best hits at "
            "\u226599% identity agreed on a single species within the assigned genus."
        )
    else:
        parts.append(
            "Taxonomy was assigned to genus level using the naive Bayesian classifier "
            "(Wang et al., 2007) against the SILVA NR99 v138.1 reference database "
            "(Quast et al., 2013)."
        )

    # Phylogenetic tree
    parts.append(
        "Representative ASV sequences were aligned using MAFFT (Katoh and Standley, 2013), "
        "and a phylogenetic tree was constructed using FastTree 2 (Price et al., 2010) "
        "under the GTR+CAT model."
    )

    # Summary statistics
    if samples:
        nonchim = [s.read_count_nonchimeric for s in samples if s.read_count_nonchimeric]
        if nonchim:
            mean_reads = int(sum(nonchim) / len(nonchim))
            min_reads = min(nonchim)
            max_reads = max(nonchim)
            parts.append(
                f"A total of {dataset.asv_count or '?'} ASVs were identified "
                f"across {dataset.sample_count or len(samples)} samples, "
                f"with a mean of {mean_reads:,} non-chimeric reads per sample "
                f"(range: {min_reads:,}–{max_reads:,})."
            )
        else:
            parts.append(
                f"A total of {dataset.asv_count or '?'} ASVs were identified "
                f"across {dataset.sample_count or len(samples)} samples."
            )

    # PICRUSt2 (if run)
    if dataset.picrust_dir_path:
        parts.append(
            "Functional potential was predicted using PICRUSt2 (Douglas et al., 2020) "
            "based on phylogenetic placement of ASV sequences against a reference tree "
            "of sequenced genomes."
        )

    # Differential abundance. Which tools a user ran isn't recorded, so all
    # five are described one sentence each for the user to trim.
    da_text = (
        "\n\nDifferential abundance between groups was tested on raw counts with "
        "up to five complementary methods. ALDEx2 (Fernandes et al., 2014) was run "
        "with 128 Monte Carlo Dirichlet instances and a centred log-ratio "
        "transformation with scale uncertainty (\u03b3 = 0.5; Nixon et al., 2025); "
        "significance was assessed with the Wilcoxon rank-sum test, and effect sizes "
        "were reported. ANCOM-BC2 (Lin and Peddada, 2024) was run with default "
        "settings, which exclude taxa present in fewer than 10% of samples; taxa were "
        "considered differentially abundant only if they also passed the ANCOM-BC2 "
        "pseudo-count sensitivity analysis. DESeq2 (Love et al., 2014) was run with "
        "size factors estimated by the \u201cposcounts\u201d method and the Wald test. "
        "LinDA (Zhou et al., 2022) was run with default settings. MaAsLin2 (Mallick "
        "et al., 2021) was run with total-sum scaling, log transformation and a "
        "linear model, excluding features present in fewer than 10% of samples."
    )
    if dataset.picrust_dir_path:
        da_text += " The same methods were applied to predicted pathway abundances."
    da_text += (
        " For all methods, P values were adjusted with the Benjamini\u2013Hochberg "
        "procedure (Benjamini and Hochberg, 1995), and features with an adjusted "
        "P value (q) < 0.05 were considered significant."
    )
    parts.append(da_text)

    # Pipeline credit
    parts.append(
        "All analyses were performed using 16S-Pipeline "
        "(https://github.com/tatsu1207/16S-Pipeline), "
        "an open-source web-based platform for end-to-end 16S rRNA amplicon analysis."
    )

    # References
    refs = [
        "Benjamini Y, Hochberg Y. "
        "Controlling the false discovery rate: a practical and powerful approach to multiple testing. "
        "Journal of the Royal Statistical Society: Series B. 1995;57(1):289-300.",

        "Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. "
        "DADA2: High-resolution sample inference from Illumina amplicon data. "
        "Nature Methods. 2016;13(7):581-583.",

        "Douglas GM, Maffei VJ, Zaneveld JR, Yurgel SN, Brown JR, Taylor CM, Huttenhower C, Langille MGI. "
        "PICRUSt2 for prediction of metagenome functions. "
        "Nature Biotechnology. 2020;38(6):685-688.",

        "Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB. "
        "Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, "
        "16S rRNA gene sequencing and selective growth experiments by compositional data analysis. "
        "Microbiome. 2014;2:15.",

        "Katoh K, Standley DM. "
        "MAFFT multiple sequence alignment software version 7: improvements in performance and usability. "
        "Molecular Biology and Evolution. 2013;30(4):772-780.",

        "Lin H, Peddada SD. "
        "Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. "
        "Nature Methods. 2024;21(1):83-91.",

        "Love MI, Huber W, Anders S. "
        "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. "
        "Genome Biology. 2014;15:550.",

        "Mallick H, Rahnavard A, McIver LJ, Ma S, Zhang Y, Nguyen LH, Tickle TL, Weingart G, Ren B, "
        "Schwager EH, Chatterjee S, Thompson KN, Wilkinson JE, Subramanian A, Lu Y, Waldron L, "
        "Paulson JN, Franzosa EA, Bravo HC, Huttenhower C. "
        "Multivariable association discovery in population-scale meta-omics studies. "
        "PLoS Computational Biology. 2021;17(11):e1009442.",

        "Martin M. "
        "Cutadapt removes adapter sequences from high-throughput sequencing reads. "
        "EMBnet.journal. 2011;17(1):10-12.",

        "Nixon MP, Gloor GB, Silverman JD. "
        "Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization. "
        "Genome Biology. 2025;26:139.",

        "Price MN, Dehal PS, Arkin AP. "
        "FastTree 2 — approximately maximum-likelihood trees for large alignments. "
        "PLoS ONE. 2010;5(3):e9490.",

        "Quast C, Pruesse E, Yilmaz P, Gerber J, Schweer T, Yarza P, Peplies J, Glöckner FO. "
        "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. "
        "Nucleic Acids Research. 2013;41(D1):D590-D596.",

        *([
            "Rognes T, Flouri T, Nichols B, Quince C, Mahé F. "
            "VSEARCH: a versatile open source tool for metagenomics. "
            "PeerJ. 2016;4:e2584.",
        ] if longread else []),

        "Wang Q, Garrity GM, Tiedje JM, Cole JR. "
        "Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. "
        "Applied and Environmental Microbiology. 2007;73(16):5261-5267.",

        "Zhou H, He K, Chen J, Zhang X. "
        "LinDA: linear models for differential abundance analysis of microbiome compositional data. "
        "Genome Biology. 2022;23:95.",
    ]

    parts.append("\n\nReferences\n" + "\n".join(f"  {i+1}. {r}" for i, r in enumerate(refs)))

    return " ".join(parts)


def _platform_name(platform: str | None) -> str:
    """Convert internal platform name to publication-friendly name."""
    names = {
        "illumina": "Illumina",
        "pacbio": "PacBio HiFi",
        "nanopore": "Oxford Nanopore",
    }
    return names.get(platform or "illumina", "Illumina")
