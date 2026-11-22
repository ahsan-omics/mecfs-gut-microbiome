#!/usr/bin/env python3
"""Alpha and beta diversity, via QIIME 2 plugins.

Two plugin calls do the work:

  core-metrics-phylogenetic   rarefies to a common depth, then computes four alpha
                              metrics, four distance matrices and their PCoAs in one
                              pass.

  alpha-rarefaction           the rarefaction curve, which is the check on whether the
                              chosen depth sits past the point where richness stops
                              climbing. Without it, the depth is a guess.

Nothing here is analysis code written for this repository. These are the standard
implementations, which is the right choice on its own merits and also keeps anything from
the private BiomeDelta pipeline well clear of a public MIT repo.

Rarefying throws reads away on purpose. The reason is that richness and most distance
measures are sensitive to depth, so a sample sequenced twice as deeply would look more
diverse for no biological reason. The depth defaults to the shallowest sample so that
every sample survives, because losing a sample costs more than losing reads.
"""

import argparse
import csv
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from statistics import median

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


def read_alpha_vector(qza_path):
    """Pull the per-sample values out of an alpha diversity artifact."""
    with zipfile.ZipFile(qza_path) as archive:
        matches = [n for n in archive.namelist() if n.endswith("alpha-diversity.tsv")]
        if not matches:
            return {}
        text = archive.read(matches[0]).decode()
    values = {}
    for row in csv.reader(text.splitlines(), delimiter="\t"):
        if len(row) < 2 or row[0].startswith("#") or row[0].lower() in ("", "sample-id"):
            continue
        try:
            values[row[0]] = float(row[1])
        except ValueError:
            continue
    return values


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--qiimedir", default="results/qiime", type=Path)
    parser.add_argument("--metadata",
                        default="metadata/sample_metadata_filtered.tsv", type=Path)
    parser.add_argument("--depth", type=int, default=None,
                        help="rarefaction depth; overrides the config")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    div = config["diversity"]

    table = args.qiimedir / "table-filtered.qza"
    tree = args.qiimedir / "rooted-tree.qza"
    for required in (table, tree, args.metadata):
        if not required.exists():
            sys.exit(f"{required} not found. Run 06_filter_and_tree.py first.")

    # ---- depth ----------------------------------------------------------------------
    depth_table = Path("results/tables/depth_after_filtering.tsv")
    if not depth_table.exists():
        sys.exit(f"{depth_table} not found. Run 06_filter_and_tree.py first.")
    depths = {r["sample-id"]: int(r["reads_after_filtering"])
              for r in csv.DictReader(depth_table.open(), delimiter="\t")}
    groups = {r["sample-id"]: r["group"]
              for r in csv.DictReader(args.metadata.open(), delimiter="\t")
              if not r["sample-id"].startswith("#")}

    values = sorted(depths.values())
    depth = args.depth or div.get("sampling_depth") or values[0]
    kept = [s for s, v in depths.items() if v >= depth]
    lost = [s for s, v in depths.items() if v < depth]

    print(f"Rarefaction depth: {depth:,}")
    print(f"  sample depths: min {values[0]:,}  median {int(median(values)):,}  "
          f"max {values[-1]:,}")
    print(f"  samples retained: {len(kept)}/{len(depths)}")
    if lost:
        print(f"  DROPPED at this depth:")
        for sample in lost:
            print(f"    {sample}  {depths[sample]:,}  ({groups.get(sample, '?')})")
    from collections import Counter
    print(f"  groups retained: {dict(Counter(groups.get(s, '?') for s in kept))}")
    discarded = sum(v - depth for v in depths.values() if v >= depth)
    print(f"  reads discarded by rarefying: {discarded:,} of "
          f"{sum(depths.values()):,} ({100 * discarded / sum(depths.values()):.1f}%)")

    # ---- core metrics ---------------------------------------------------------------
    core = args.qiimedir / "core-metrics"

    # Completeness is judged by the files themselves, never by whether the directory is
    # there. An interrupted run leaves the directory populated but missing exactly the
    # outputs that were in flight, and a directory-existence check would then skip the
    # work, exit cleanly, and leave the pipeline broken several steps later.
    expected = [core / f"{m}_vector.qza" for m in div["alpha_metrics"]]
    expected += [core / f"{m}_distance_matrix.qza" for m in div["beta_metrics"]]
    expected += [core / f"{m}_pcoa_results.qza" for m in div["beta_metrics"]]
    missing = [p for p in expected if not p.exists()]

    if args.force or missing:
        if core.exists():
            if missing and not args.force:
                print(f"\n    {core} exists but {len(missing)} of {len(expected)} "
                      f"expected artifacts are missing, for example "
                      f"{missing[0].name}.")
                print(f"    Rebuilding it: core-metrics-phylogenetic cannot write into "
                      f"an existing directory,")
                print(f"    so a partial one has to be removed rather than topped up.")
            shutil.rmtree(core)
        run(["qiime", "diversity", "core-metrics-phylogenetic",
             "--i-phylogeny", tree,
             "--i-table", table,
             "--p-sampling-depth", depth,
             "--m-metadata-file", args.metadata,
             "--p-n-jobs-or-threads", div["threads"],
             "--output-dir", core],
            f"rarefy to {depth:,} and compute alpha, beta and PCoA")
    else:
        print(f"\n>>> core-metrics: all {len(expected)} artifacts present, skipping "
              f"(use --force to redo)")

    # ---- rarefaction curve ------------------------------------------------------------
    curve = args.qiimedir / "alpha-rarefaction.qzv"
    if args.force or not curve.exists():
        maximum = div.get("rarefaction_max_depth") or int(median(values))
        run(["qiime", "diversity", "alpha-rarefaction",
             "--i-table", table,
             "--i-phylogeny", tree,
             "--p-max-depth", maximum,
             "--p-steps", div["rarefaction_steps"],
             "--p-iterations", div["rarefaction_iterations"],
             "--m-metadata-file", args.metadata,
             "--o-visualization", curve],
            f"rarefaction curve to {maximum:,}")
    else:
        print(f"\n>>> alpha-rarefaction: {curve} exists, skipping")

    # ---- export the alpha vectors ----------------------------------------------------
    tables_dir = Path("results/tables")
    tables_dir.mkdir(parents=True, exist_ok=True)

    collected = {}
    for metric in div["alpha_metrics"]:
        vector = core / f"{metric}_vector.qza"
        if vector.exists():
            collected[metric] = read_alpha_vector(vector)
        else:
            print(f"  note: {vector.name} not found", file=sys.stderr)

    if collected:
        out = tables_dir / "alpha_diversity.tsv"
        metrics = list(collected)
        samples = sorted(next(iter(collected.values())))
        with out.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample-id", "group"] + metrics)
            for sample in samples:
                writer.writerow([sample, groups.get(sample, "")]
                                + [collected[m].get(sample, "") for m in metrics])
        print(f"\nWrote {out}")

        # Descriptive only. The tests are script 08's job; this is a look, not a result.
        print(f"\nAlpha diversity by group (medians, descriptive only, not a test):")
        labels = sorted(set(groups.values()))
        header = "  " + "metric".ljust(20) + "".join(l.ljust(16) for l in labels)
        print(header)
        for metric in metrics:
            cells = ""
            for label in labels:
                vals = sorted(v for s, v in collected[metric].items()
                              if groups.get(s) == label)
                cells += (f"{median(vals):.3f}" if vals else "-").ljust(16)
            print(f"  {metric.ljust(20)}{cells}")

    print(f"\nArtifacts in {core}:")
    print(f"  *_vector.qza            alpha diversity per sample")
    print(f"  *_distance_matrix.qza   beta diversity between samples")
    print(f"  *_pcoa_results.qza      ordinations")
    print(f"  *_emperor.qzv           interactive PCoA plots")
    print(f"\nOpen {curve} at https://view.qiime2.org to check the rarefaction curve.")
    print(f"What to look for: the richness curves should flatten before {depth:,}. If they")
    print(f"are still climbing there, the depth is undersampling the richer samples and")
    print(f"observed_features is measuring sequencing effort as much as biology.")


if __name__ == "__main__":
    main()
