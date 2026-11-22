#!/usr/bin/env python3
"""Import into QIIME 2, trim primers, and denoise into ASVs.

Four steps, each skipped if its output already exists so the script is safe to re-run:

  1. Build a manifest of FORWARD reads only, excluding the failed libraries.
  2. qiime tools import       -> SampleData[SequencesWithQuality]
  3. qiime cutadapt trim-single -> removes 341F from the 5' end
  4. qiime dada2 denoise-single -> feature table, representative sequences, stats

Then it reads the denoising stats back and reports retention at every stage, because
that is the evidence for whether the binned quality scores have broken DADA2's error
model. See `denoising.expect_binned_quality_issue` in config.yaml.

Every qiime command is printed before it runs, so a failure can be reproduced by hand.
"""

import argparse
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


def read_metadata(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    # QIIME 2 metadata carries a #q2:types line that is not a sample.
    return [r for r in rows if not r["sample-id"].startswith("#")]


def extract_from_qza(qza_path, inner_name, destination):
    """Pull one file out of a .qza, which is just a zip with a uuid directory inside."""
    with zipfile.ZipFile(qza_path) as archive:
        matches = [n for n in archive.namelist() if n.endswith(inner_name)]
        if not matches:
            return None
        with archive.open(matches[0]) as src, open(destination, "wb") as dst:
            dst.write(src.read())
    return destination


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--rawdir", default="data/raw", type=Path)
    parser.add_argument("--metadata", default="metadata/sample_metadata.tsv", type=Path)
    parser.add_argument("--diagnostics",
                        default="results/tables/read_diagnostics_by_sample.tsv", type=Path)
    parser.add_argument("--outdir", default="results/qiime", type=Path)
    parser.add_argument("--force", action="store_true",
                        help="recompute steps whose outputs already exist")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    primers = config["primers"]
    cut = config["cutadapt"]
    den = config["denoising"]
    min_reads = int(config["sample_filtering"]["min_reads"])

    if not args.diagnostics.exists():
        sys.exit(f"{args.diagnostics} not found. Run 03_read_diagnostics.py first.")

    args.outdir.mkdir(parents=True, exist_ok=True)
    tables = Path("results/tables")
    tables.mkdir(parents=True, exist_ok=True)

    # ---- 1. manifest, with the failed libraries excluded ----------------------------
    depths = {r["sample"]: int(r["reads_forward"])
              for r in csv.DictReader(args.diagnostics.open(), delimiter="\t")}
    metadata = read_metadata(args.metadata)

    keep, dropped = [], []
    for row in metadata:
        sample = row["sample-id"]
        (keep if depths.get(sample, 0) >= min_reads else dropped).append(row)

    print(f"Samples: {len(metadata)} total, {len(keep)} kept, {len(dropped)} excluded "
          f"below {min_reads:,} reads")
    for row in dropped:
        print(f"  excluded {row['sample-id']}  {depths.get(row['sample-id'], 0):,} reads"
              f"  ({row['group']})")
    from collections import Counter
    print(f"  groups kept: {dict(Counter(r['group'] for r in keep))}")

    manifest_path = args.outdir / "manifest.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample-id", "absolute-filepath"])
        for row in keep:
            sample = row["sample-id"]
            path = (args.rawdir / f"{sample}_1.fastq.gz").resolve()
            if not path.exists():
                sys.exit(f"{path} not found.")
            writer.writerow([sample, path])
    print(f"\nWrote {manifest_path} (forward reads only)")

    # A metadata file matching the kept samples, for every downstream step.
    filtered_metadata = Path("metadata/sample_metadata_filtered.tsv")
    with filtered_metadata.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        columns = list(metadata[0].keys())
        writer.writerow(columns)
        types = next((r for r in csv.DictReader(args.metadata.open(), delimiter="\t")
                      if r["sample-id"].startswith("#")), None)
        if types:
            writer.writerow([types[c] for c in columns])
        for row in keep:
            writer.writerow([row[c] for c in columns])
    print(f"Wrote {filtered_metadata} ({len(keep)} samples)")

    # ---- 2. import ------------------------------------------------------------------
    demux = args.outdir / "demux.qza"
    if args.force or not demux.exists():
        run(["qiime", "tools", "import",
             "--type", "SampleData[SequencesWithQuality]",
             "--input-format", "SingleEndFastqManifestPhred33V2",
             "--input-path", manifest_path,
             "--output-path", demux],
            "import forward reads")
    else:
        print(f"\n>>> import: {demux} exists, skipping")

    # ---- 3. primer removal ----------------------------------------------------------
    trimmed = args.outdir / "trimmed.qza"
    if args.force or not trimmed.exists():
        command = ["qiime", "cutadapt", "trim-single",
                   "--i-demultiplexed-sequences", demux,
                   "--p-front", primers["forward"]["sequence"],
                   "--p-error-rate", cut["error_rate"],
                   "--p-cores", cut["cores"],
                   "--o-trimmed-sequences", trimmed]
        if cut.get("discard_untrimmed"):
            command.append("--p-discard-untrimmed")
        run(command, f"remove {primers['forward']['name']} from the 5' end")
    else:
        print(f"\n>>> cutadapt: {trimmed} exists, skipping")

    # ---- 4. denoise -----------------------------------------------------------------
    table = args.outdir / "table.qza"
    rep_seqs = args.outdir / "rep-seqs.qza"
    stats = args.outdir / "denoising-stats.qza"
    if args.force or not table.exists():
        print(f"\n    This is the slow step. 61 samples on {den['threads']} threads will")
        print(f"    take a while; DADA2 gives no progress output until it finishes.")
        run(["qiime", "dada2", "denoise-single",
             "--i-demultiplexed-seqs", trimmed,
             "--p-trunc-len", den["trunc_len"],
             "--p-trim-left", den["trim_left"],
             "--p-max-ee", den["max_ee"],
             "--p-trunc-q", den["trunc_q"],
             "--p-chimera-method", den["chimera_method"],
             "--p-pooling-method", den["pooling_method"],
             "--p-n-threads", den["threads"],
             "--o-table", table,
             "--o-representative-sequences", rep_seqs,
             "--o-denoising-stats", stats],
            f"denoise into ASVs (trunc-len {den['trunc_len']})")
    else:
        print(f"\n>>> dada2: {table} exists, skipping")

    # ---- 5. read the stats back -----------------------------------------------------
    stats_tsv = tables / "denoising_stats.tsv"
    extracted = extract_from_qza(stats, "stats.tsv", stats_tsv)
    if not extracted:
        print("\nCould not read denoising stats out of the .qza; inspect it with "
              "`qiime metadata tabulate`.", file=sys.stderr)
        return

    rows = [r for r in csv.DictReader(stats_tsv.open(), delimiter="\t")
            if not r[next(iter(r))].startswith("#")]
    if not rows:
        return

    def column(name):
        key = next((k for k in rows[0] if k.lower().replace("-", "-") == name), None)
        return [float(r[key]) for r in rows] if key else None

    header = list(rows[0].keys())
    numeric = [h for h in header[1:] if h.lower() not in ("sample-id",)]
    print(f"\n{'stage':<24}{'total':>14}{'% of input':>12}")
    values = {}
    for name in numeric:
        try:
            values[name] = sum(float(r[name]) for r in rows)
        except (ValueError, KeyError):
            continue
    baseline = values.get("input") or (list(values.values())[0] if values else 0)
    for name, total in values.items():
        if "percentage" in name.lower():
            continue
        share = 100 * total / baseline if baseline else 0
        print(f"  {name:<22}{total:>14,.0f}{share:>11.1f}%")

    print(f"\nPer-sample retention (input -> non-chimeric), worst 5:")
    try:
        ordered = sorted(
            rows, key=lambda r: float(r["non-chimeric"]) / max(1.0, float(r["input"])))
        for r in ordered[:5]:
            sample = r[header[0]]
            kept_pct = 100 * float(r["non-chimeric"]) / max(1.0, float(r["input"]))
            print(f"  {sample:<16}{float(r['input']):>10,.0f} -> "
                  f"{float(r['non-chimeric']):>10,.0f}  ({kept_pct:.1f}%)")
    except KeyError:
        pass

    print(f"\nWrote {stats_tsv}")
    print(f"\nArtifacts in {args.outdir}: table.qza, rep-seqs.qza, denoising-stats.qza")
    print(f"\nWhat to look at: retention through 'non-chimeric' should be well over half.")
    print(f"If it is poor, or the ASV count is implausibly high for 61 gut samples, that")
    print(f"is the binned quality scores breaking DADA2's error model, and we move the")
    print(f"denoising into R with a monotonic error function. See config.yaml.")


if __name__ == "__main__":
    main()
