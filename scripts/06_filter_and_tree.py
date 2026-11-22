#!/usr/bin/env python3
"""Curate the feature table, then build the phylogeny from what survives.

Order matters here. Taxonomy comes first (script 05), then filtering, then the tree, so
the phylogeny describes the features that actually get analysed rather than including
branches that are about to be thrown away.

What is removed:

  mitochondria and chloroplast   16S-like sequences from host and plant material. They
                                 amplify with bacterial primers and are not part of the
                                 bacterial community.

  no phylum-level assignment     at 134 bp, a sequence that cannot be placed even to
                                 phylum is far more likely an artefact than a novel
                                 lineage, and it would contribute to UniFrac through a
                                 branch nobody can interpret.

What is deliberately NOT removed: rare features. Prevalence filtering before diversity
would bias alpha diversity downward in exactly the samples that carry the most rare
taxa. That filter belongs at the differential abundance step and is applied there.

Then mafft, masking, FastTree and midpoint rooting, which UniFrac and Faith's PD need.

Finishes by reporting per-sample depth on the FILTERED table, because that is what sets
the rarefaction depth, not the pre-filter counts.
"""

import argparse
import collections
import csv
import subprocess
import sys
import zipfile
from pathlib import Path

import yaml


def run(command, description):
    print(f"\n>>> {description}")
    print("    " + " ".join(str(c) for c in command))
    result = subprocess.run([str(c) for c in command], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\nCommand failed with exit code {result.returncode}:", file=sys.stderr)
        print((result.stderr or result.stdout).strip()[-3000:], file=sys.stderr)
        sys.exit(1)
    return result.stdout


def count_sequences(qza_path):
    with zipfile.ZipFile(qza_path) as archive:
        matches = [n for n in archive.namelist() if n.endswith("dna-sequences.fasta")]
        if not matches:
            return None
        with archive.open(matches[0]) as handle:
            return sum(1 for line in handle if line.startswith(b">"))


def load_table(qza_path, workdir):
    """Export a FeatureTable artifact and load it with biom."""
    import biom
    out = workdir / qza_path.stem
    run(["qiime", "tools", "export", "--input-path", qza_path,
         "--output-path", out], f"export {qza_path.name}")
    return biom.load_table(out / "feature-table.biom")


def rank_present(lineage, prefix):
    """True only if the rank exists AND is populated, not an empty 'p__' placeholder."""
    for part in (p.strip() for p in lineage.split(";")):
        if part.startswith(prefix):
            return len(part) > len(prefix)
    return False


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--qiimedir", default="results/qiime", type=Path)
    parser.add_argument("--workdir", default="results/export", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    filt = config["filtering"]
    phylo = config["phylogeny"]

    table = args.qiimedir / "table.qza"
    rep_seqs = args.qiimedir / "rep-seqs.qza"
    taxonomy = args.qiimedir / "taxonomy.qza"
    for required in (table, rep_seqs, taxonomy):
        if not required.exists():
            sys.exit(f"{required} not found. Run the earlier scripts first.")

    args.workdir.mkdir(parents=True, exist_ok=True)
    tables_dir = Path("results/tables")
    tables_dir.mkdir(parents=True, exist_ok=True)

    exclude = ",".join(filt["exclude_taxa"])
    include = filt["require_rank"]

    # ---- what the filter will do, before doing it -----------------------------------
    taxonomy_tsv = tables_dir / "taxonomy.tsv"
    if taxonomy_tsv.exists():
        rows = [r for r in csv.DictReader(taxonomy_tsv.open(), delimiter="\t")
                if not list(r.values())[0].startswith("#")]
        taxon_key = next(k for k in rows[0] if k.lower() == "taxon")
        contaminant = sum(
            1 for r in rows
            if any(t in r[taxon_key].lower() for t in filt["exclude_taxa"]))
        no_phylum = sum(1 for r in rows if not rank_present(r[taxon_key], "p__"))
        print(f"Of {len(rows):,} ASVs:")
        print(f"  {contaminant:,} match {exclude}")
        print(f"  {no_phylum:,} have no populated phylum")

    # ---- filter ---------------------------------------------------------------------
    filtered_table = args.qiimedir / "table-filtered.qza"
    filtered_seqs = args.qiimedir / "rep-seqs-filtered.qza"

    if args.force or not filtered_table.exists():
        run(["qiime", "taxa", "filter-table",
             "--i-table", table,
             "--i-taxonomy", taxonomy,
             "--p-exclude", exclude,
             "--p-include", include,
             "--o-filtered-table", filtered_table],
            f"drop {exclude} and anything without a phylum, from the table")
    else:
        print(f"\n>>> filter-table: {filtered_table} exists, skipping")

    if args.force or not filtered_seqs.exists():
        run(["qiime", "taxa", "filter-seqs",
             "--i-sequences", rep_seqs,
             "--i-taxonomy", taxonomy,
             "--p-exclude", exclude,
             "--p-include", include,
             "--o-filtered-sequences", filtered_seqs],
            "apply the same filter to the representative sequences")
    else:
        print(f"\n>>> filter-seqs: {filtered_seqs} exists, skipping")

    before = count_sequences(rep_seqs)
    after = count_sequences(filtered_seqs)
    if before and after:
        print(f"\n  ASVs: {before:,} -> {after:,}  "
              f"({100 * (before - after) / before:.1f}% removed)")

    # ---- tree -----------------------------------------------------------------------
    rooted = args.qiimedir / "rooted-tree.qza"
    if args.force or not rooted.exists():
        print(f"\n    mafft on {after:,} sequences of 134 bp, then FastTree. A short")
        print(f"    alignment carries less phylogenetic signal than full-length 16S, so")
        print(f"    UniFrac here is coarser than usual. Method is still right.")
        run(["qiime", "phylogeny", "align-to-tree-mafft-fasttree",
             "--i-sequences", filtered_seqs,
             "--p-n-threads", phylo["threads"],
             "--o-alignment", args.qiimedir / "aligned-rep-seqs.qza",
             "--o-masked-alignment", args.qiimedir / "masked-aligned-rep-seqs.qza",
             "--o-tree", args.qiimedir / "unrooted-tree.qza",
             "--o-rooted-tree", rooted],
            "align, mask, build and root the tree")
    else:
        print(f"\n>>> phylogeny: {rooted} exists, skipping")

    # ---- depth on the filtered table, which sets the rarefaction depth --------------
    try:
        biom_table = load_table(filtered_table, args.workdir)
    except ImportError:
        print("\nbiom is not importable; skipping the depth summary.", file=sys.stderr)
        return

    sample_ids = list(biom_table.ids("sample"))
    depths = {s: int(biom_table.sum("sample")[i]) for i, s in enumerate(sample_ids)}

    metadata_path = Path("metadata/sample_metadata_filtered.tsv")
    groups = {}
    if metadata_path.exists():
        for row in csv.DictReader(metadata_path.open(), delimiter="\t"):
            if not row["sample-id"].startswith("#"):
                groups[row["sample-id"]] = row["group"]

    ordered = sorted(depths.items(), key=lambda kv: kv[1])
    values = [v for _, v in ordered]
    print(f"\nDepth after filtering, {len(values)} samples:")
    print(f"  min {values[0]:,}   median {values[len(values) // 2]:,}   "
          f"max {values[-1]:,}")
    print(f"\n  shallowest 5:")
    for sample, depth in ordered[:5]:
        print(f"    {sample:<16}{depth:>10,}  ({groups.get(sample, '?')})")

    print(f"\n  samples retained at each candidate rarefaction depth:")
    for candidate in (5000, 10000, 15000, 20000, 25000, 30000):
        kept = sum(1 for v in values if v >= candidate)
        by_group = collections.Counter(
            groups.get(s, "?") for s, v in depths.items() if v >= candidate)
        print(f"    {candidate:>7,}: {kept:>3}/{len(values)} samples  {dict(by_group)}")

    depth_table = tables_dir / "depth_after_filtering.tsv"
    with depth_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample-id", "group", "reads_after_filtering"])
        for sample, depth in ordered:
            writer.writerow([sample, groups.get(sample, ""), depth])
    print(f"\nWrote {depth_table}")

    print(f"\nPick the rarefaction depth from the table above: the aim is to keep every")
    print(f"sample if the shallowest is not an outlier, since dropping samples costs more")
    print(f"than the reads do. 07_diversity.py takes it as a parameter.")


if __name__ == "__main__":
    main()
