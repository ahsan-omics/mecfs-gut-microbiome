#!/usr/bin/env python3
"""Assign taxonomy to the ASVs against SILVA 138, using vsearch consensus classification.

Three steps, each cached so the script is safe to re-run:

  1. Download the SILVA 138 99% reference sequences and taxonomy.
  2. Extract the 341F/805R amplicon from those references and truncate it to the read
     length, so the reference describes the same stretch of gene as the queries.
  3. classify-consensus-vsearch.

The naive Bayes classifier is deliberately not used. It needs more RAM than this machine
has, while vsearch streams alignments and needs very little. Trimming the reference to
the amplicon takes it from about 700 Mbp to 67 Mbp, which is what makes this comfortable
rather than merely possible. The reasoning is in config.yaml under `taxonomy`.

Reports how deep the assignments go, which is the thing to look at: with 134 bp of V3,
genus-level resolution is the realistic ceiling and species-level calls should be rare.
"""

import argparse
import collections
import csv
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

import yaml

LEVELS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
USER_AGENT = "mecfs-gut-microbiome"


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def run(command, description):
    print(f"\n>>> {description}")
    print("    " + " ".join(str(c) for c in command))
    result = subprocess.run([str(c) for c in command], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\nCommand failed with exit code {result.returncode}:", file=sys.stderr)
        print((result.stderr or result.stdout).strip()[-3000:], file=sys.stderr)
        sys.exit(1)
    return result.stdout


def download(url, destination):
    if destination.exists() and destination.stat().st_size > 0:
        print(f"    {destination.name} already present ({human(destination.stat().st_size)})")
        return
    print(f"    fetching {destination.name} ...", end="", flush=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=900) as response, \
            open(temporary, "wb") as handle:
        while True:
            block = response.read(1 << 20)
            if not block:
                break
            handle.write(block)
    temporary.rename(destination)
    print(f" {human(destination.stat().st_size)}")


def extract_from_qza(qza_path, suffix, destination):
    with zipfile.ZipFile(qza_path) as archive:
        matches = [n for n in archive.namelist() if n.endswith(suffix)]
        if not matches:
            return None
        with archive.open(matches[0]) as src, open(destination, "wb") as dst:
            dst.write(src.read())
    return destination


def count_sequences(qza_path):
    with zipfile.ZipFile(qza_path) as archive:
        matches = [n for n in archive.namelist() if n.endswith("dna-sequences.fasta")]
        if not matches:
            return None
        with archive.open(matches[0]) as handle:
            return sum(1 for line in handle if line.startswith(b">"))


def depth_of(lineage):
    """How many ranks are actually populated, ignoring empty 'g__' style placeholders."""
    parts = [p.strip() for p in lineage.split(";")]
    depth = 0
    for part in parts:
        body = part.split("__", 1)[-1] if "__" in part else part
        if not body or body.lower() in {"unassigned", "unclassified", "uncultured"}:
            break
        depth += 1
    return depth


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--qiimedir", default="results/qiime", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    tax = config["taxonomy"]
    primers = config["primers"]
    trunc_len = int(config["denoising"]["trunc_len"])

    rep_seqs = args.qiimedir / "rep-seqs.qza"
    if not rep_seqs.exists():
        sys.exit(f"{rep_seqs} not found. Run 04_denoise.py first.")
    asv_count = count_sequences(rep_seqs)
    print(f"{asv_count:,} ASVs to classify, {trunc_len} bp each")

    reference_dir = Path(tax["reference_dir"])
    reference_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. reference ---------------------------------------------------------------
    print(f"\n>>> SILVA 138 99% reference into {reference_dir}/")
    seqs = reference_dir / "silva-138-99-seqs.qza"
    taxonomy_ref = reference_dir / "silva-138-99-tax.qza"
    download(tax["reference_seqs_url"], seqs)
    download(tax["reference_tax_url"], taxonomy_ref)

    # ---- 2. trim the reference to our amplicon ---------------------------------------
    search_seqs = seqs
    if tax.get("extract_region", True):
        extracted = reference_dir / f"silva-138-99-v3-{trunc_len}bp.qza"
        if args.force or not extracted.exists():
            print(f"\n    Trimming the reference to the {primers['forward']['name']}/"
                  f"{primers['reverse']['name']} amplicon and truncating to {trunc_len} bp.")
            print(f"    This is the slow step, and it only has to happen once.")
            run(["qiime", "feature-classifier", "extract-reads",
                 "--i-sequences", seqs,
                 "--p-f-primer", primers["forward"]["sequence"],
                 "--p-r-primer", primers["reverse"]["sequence"],
                 "--p-trunc-len", trunc_len,
                 "--p-min-length", tax["extract_min_length"],
                 "--p-max-length", tax["extract_max_length"],
                 "--p-n-jobs", tax["threads"],
                 "--o-reads", extracted],
                "extract the V3 region from SILVA")
        else:
            print(f"\n>>> extract-reads: {extracted} exists, skipping")
        search_seqs = extracted
        n = count_sequences(extracted)
        if n:
            print(f"    reference now {n:,} sequences x {trunc_len} bp "
                  f"(~{n * trunc_len / 1e6:.0f} Mbp)")

    # ---- 3. classify ------------------------------------------------------------------
    taxonomy_qza = args.qiimedir / "taxonomy.qza"
    # --o-search-results is a REQUIRED output on this action, not an optional extra.
    # It holds the top hits per query in BLAST6 format, which is also the thing to look
    # at if an assignment ever looks wrong.
    hits_qza = args.qiimedir / "vsearch-hits.qza"
    if args.force or not taxonomy_qza.exists():
        run(["qiime", "feature-classifier", "classify-consensus-vsearch",
             "--i-query", rep_seqs,
             "--i-reference-reads", search_seqs,
             "--i-reference-taxonomy", taxonomy_ref,
             "--p-perc-identity", tax["perc_identity"],
             "--p-query-cov", tax["query_cov"],
             "--p-maxaccepts", tax["maxaccepts"],
             "--p-min-consensus", tax["min_consensus"],
             "--p-threads", tax["threads"],
             "--o-classification", taxonomy_qza,
             "--o-search-results", hits_qza],
            "classify ASVs against SILVA")
    else:
        print(f"\n>>> classify: {taxonomy_qza} exists, skipping")

    # ---- 4. report -------------------------------------------------------------------
    tables = Path("results/tables")
    tables.mkdir(parents=True, exist_ok=True)
    taxonomy_tsv = tables / "taxonomy.tsv"
    if not extract_from_qza(taxonomy_qza, "taxonomy.tsv", taxonomy_tsv):
        sys.exit("could not read taxonomy.tsv out of the classification artifact")

    rows = [r for r in csv.DictReader(taxonomy_tsv.open(), delimiter="\t")
            if not list(r.values())[0].startswith("#")]
    key = list(rows[0].keys())
    taxon_key = next(k for k in key if k.lower() == "taxon")

    depths = collections.Counter()
    unassigned = 0
    phyla = collections.Counter()
    for row in rows:
        lineage = row[taxon_key]
        if lineage.strip().lower().startswith("unassigned"):
            unassigned += 1
            depths[0] += 1
            continue
        d = depth_of(lineage)
        depths[d] += 1
        if d >= 2:
            phyla[lineage.split(";")[1].strip()] += 1

    total = len(rows)
    print(f"\nClassified {total:,} ASVs")
    print(f"  unassigned: {unassigned:,} ({100 * unassigned / total:.1f}%)")
    print(f"\nHow deep the assignments reach (cumulative):")
    for level_index, level in enumerate(LEVELS, start=1):
        reached = sum(v for d, v in depths.items() if d >= level_index)
        print(f"  to {level:<8} {reached:>6,}  {100 * reached / total:5.1f}%")

    print(f"\nTop phyla by ASV count (a gut sanity check):")
    for name, count in phyla.most_common(8):
        print(f"  {name:<28} {count:>6,}  {100 * count / total:5.1f}%")

    print(f"\nWrote {taxonomy_tsv}")
    print(f"\nWith {trunc_len} bp of V3, genus is the realistic ceiling. A low species")
    print(f"percentage here is the sequencing design showing through, not a failure.")


if __name__ == "__main__":
    main()
