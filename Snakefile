# ME/CFS gut 16S amplicon analysis
#
# Each rule is one script, each script is one argparse command with fixed inputs and
# outputs. The Snakefile only wires them together; no analysis logic lives here.
#
# TWO ENVIRONMENTS, and Snakemake itself lives in neither. Analysis needs QIIME 2, the
# download step needs sra-tools, and snakemake is usually installed somewhere else again.
# So every rule names its own environment through `conda run` instead of inheriting
# whichever environment you happened to launch from. That makes the pipeline independent
# of how it is invoked, and it is also why `qiime` is reachable from inside the rules:
# `conda run` sets CONDA_PREFIX, which the QIIME 2 CLI requires and a bare PATH entry
# does not provide.
#
# Environment names are in config/config.yaml under `environments`.
#
# Run from the project root, with snakemake on your PATH from any environment:
#     snakemake -n                 dry run, always first
#     snakemake --cores 4
#
# The scripts all skip work whose output already exists, so re-running is cheap and an
# interrupted run picks up where it stopped.

configfile: "config/config.yaml"

QIIME_DIR = "results/qiime"
TABLES = "results/tables"
FIGURES = "results/figures"

ENVS = config.get("environments", {})
CONDA = ENVS.get("conda_command", "conda")
QIIME_ENV = ENVS.get("qiime2", "qiime2-amplicon-2025.7")
SRA_ENV = ENVS.get("sra_tools", "bioenv")


def in_env(env, command):
    """Run a command inside a named conda environment."""
    return f"{CONDA} run -n {env} --no-capture-output {command}"


def qiime_script(name):
    return in_env(QIIME_ENV, f"python scripts/{name}")


rule all:
    input:
        f"{FIGURES}/read_quality_light.png",
        f"{FIGURES}/pcoa_weighted_unifrac_light.png",
        f"{FIGURES}/alpha_diversity_light.png",
        f"{FIGURES}/differential_abundance_light.png",
        f"{TABLES}/differential_abundance.tsv",


rule fetch_metadata:
    output:
        metadata="metadata/sample_metadata.tsv",
        manifest="data/download_manifest.tsv",
    log:
        "logs/01_fetch_metadata.log",
    shell:
        qiime_script("01_fetch_metadata.py") + " > {log} 2>&1"


rule download_reads:
    input:
        manifest=rules.fetch_metadata.output.manifest,
    output:
        # The real outputs are 126 FASTQ files whose names come from the manifest, which
        # Snakemake cannot know ahead of time. This record lists every one of them with
        # its spot count, so it is both the tracked output and the provenance file.
        record="data/download_record.tsv",
    log:
        "logs/02_download_reads.log",
    shell:
        in_env(SRA_ENV, "python scripts/02_download_reads.py") + " > {log} 2>&1"


rule read_diagnostics:
    input:
        record=rules.download_reads.output.record,
        metadata=rules.fetch_metadata.output.metadata,
    output:
        by_sample=f"{TABLES}/read_diagnostics_by_sample.tsv",
        by_position=f"{TABLES}/quality_by_position.tsv",
        light=f"{FIGURES}/read_quality_light.png",
        dark=f"{FIGURES}/read_quality_dark.png",
    log:
        "logs/03_read_diagnostics.log",
    shell:
        qiime_script("03_read_diagnostics.py") + " > {log} 2>&1"


rule denoise:
    input:
        diagnostics=rules.read_diagnostics.output.by_sample,
        metadata=rules.fetch_metadata.output.metadata,
    output:
        table=f"{QIIME_DIR}/table.qza",
        rep_seqs=f"{QIIME_DIR}/rep-seqs.qza",
        stats=f"{QIIME_DIR}/denoising-stats.qza",
        filtered_metadata="metadata/sample_metadata_filtered.tsv",
        stats_tsv=f"{TABLES}/denoising_stats.tsv",
    log:
        "logs/04_denoise.log",
    shell:
        qiime_script("04_denoise.py") + " > {log} 2>&1"


rule assign_taxonomy:
    input:
        rep_seqs=rules.denoise.output.rep_seqs,
    output:
        taxonomy=f"{QIIME_DIR}/taxonomy.qza",
        taxonomy_tsv=f"{TABLES}/taxonomy.tsv",
    log:
        "logs/05_assign_taxonomy.log",
    shell:
        qiime_script("05_assign_taxonomy.py") + " > {log} 2>&1"


rule filter_and_tree:
    input:
        table=rules.denoise.output.table,
        taxonomy=rules.assign_taxonomy.output.taxonomy,
        taxonomy_tsv=rules.assign_taxonomy.output.taxonomy_tsv,
    output:
        table=f"{QIIME_DIR}/table-filtered.qza",
        rep_seqs=f"{QIIME_DIR}/rep-seqs-filtered.qza",
        tree=f"{QIIME_DIR}/rooted-tree.qza",
        depth=f"{TABLES}/depth_after_filtering.tsv",
    log:
        "logs/06_filter_and_tree.log",
    shell:
        qiime_script("06_filter_and_tree.py") + " > {log} 2>&1"


rule diversity:
    input:
        table=rules.filter_and_tree.output.table,
        tree=rules.filter_and_tree.output.tree,
        depth=rules.filter_and_tree.output.depth,
        metadata=rules.denoise.output.filtered_metadata,
    output:
        # Named files rather than the output directory, because core-metrics-phylogenetic
        # refuses to write into a directory that already exists, and a directory() output
        # would have Snakemake remove and recreate it on every touch.
        weighted=f"{QIIME_DIR}/core-metrics/weighted_unifrac_distance_matrix.qza",
        unweighted=f"{QIIME_DIR}/core-metrics/unweighted_unifrac_distance_matrix.qza",
        bray=f"{QIIME_DIR}/core-metrics/bray_curtis_distance_matrix.qza",
        jaccard=f"{QIIME_DIR}/core-metrics/jaccard_distance_matrix.qza",
        pcoa=f"{QIIME_DIR}/core-metrics/weighted_unifrac_pcoa_results.qza",
        shannon=f"{QIIME_DIR}/core-metrics/shannon_vector.qza",
        rarefaction=f"{QIIME_DIR}/alpha-rarefaction.qzv",
        alpha=f"{TABLES}/alpha_diversity.tsv",
    log:
        "logs/07_diversity.log",
    shell:
        qiime_script("07_diversity.py") + " > {log} 2>&1"


rule statistics:
    input:
        weighted=rules.diversity.output.weighted,
        shannon=rules.diversity.output.shannon,
        table=rules.filter_and_tree.output.table,
        taxonomy=rules.assign_taxonomy.output.taxonomy,
        metadata=rules.denoise.output.filtered_metadata,
    output:
        differentials=f"{QIIME_DIR}/stats/ancombc-adjusted.qza",
        unadjusted=f"{QIIME_DIR}/stats/ancombc-unadjusted.qza",
        genus=f"{QIIME_DIR}/genus-table-filtered.qza",
        da_table=f"{TABLES}/differential_abundance.tsv",
        ancombc_metadata="metadata/sample_metadata_ancombc.tsv",
    log:
        "logs/08_statistics.log",
    shell:
        qiime_script("08_statistics.py") + " > {log} 2>&1"


rule figures:
    input:
        pcoa=rules.diversity.output.pcoa,
        alpha=rules.diversity.output.alpha,
        differentials=rules.statistics.output.differentials,
        genus=rules.statistics.output.genus,
    output:
        expand(
            f"{FIGURES}/{{name}}_{{mode}}.png",
            name=["pcoa_weighted_unifrac", "alpha_diversity", "differential_abundance"],
            mode=["light", "dark"],
        ),
    log:
        "logs/09_figures.log",
    shell:
        qiime_script("09_figures.py") + " > {log} 2>&1"
