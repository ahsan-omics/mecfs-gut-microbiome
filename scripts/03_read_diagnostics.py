#!/usr/bin/env python3
"""Read diagnostics: quality by position, primer content, and the overlap test.

Three questions, each of which has to be answered before DADA2 can be run sensibly.

  Where do we truncate?   The usual advice is to read it off the quality curve, but these
                          scores are binned to Q11/Q25/Q37 so there is no smooth curve.
                          The rule used instead is in config: keep the longest prefix in
                          which at least 80 percent of bases are Q25 or better. Two
                          alternatives are reported next to it so the choice is visible.

  Are the primers on?     R1 carries 341F and R2 carries 805R. Primer sequence is not
                          biological variation, and because these primers are degenerate,
                          leaving them on splits one organism into several ASVs.

  Can the mates merge?    V3-V4 is about 460 bp and the reads are 2 x 151, so there should
                          be no overlap. This tests it on real read pairs rather than
                          asserting it from arithmetic, because DADA2's merge step does not
                          error on non-overlapping reads. It silently discards almost
                          everything and reports success.

Outputs two tables and one figure, in light and dark versions for the README.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import yaml

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# Validated with the data-viz palette validator, both modes, all checks pass.
THEME = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e",
              "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7",
              "series": ["#2a78d6", "#eb6834"]},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7",
             "muted": "#898781", "grid": "#2c2c2a", "axis": "#383835",
             "series": ["#3987e5", "#d95926"]},
}


def iupac_regex(sequence):
    import re
    return re.compile("".join(IUPAC[b] for b in sequence.upper()))


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def scan_file(path, subsample, primer_pattern):
    """One pass: exact read count, plus per-position quality on a subsample."""
    max_length = 0
    qsum = np.zeros(512, dtype=np.int64)
    qn = np.zeros(512, dtype=np.int64)
    # Counts of bases at or above each of the observed bins, per position.
    at_least_25 = np.zeros(512, dtype=np.int64)
    at_least_37 = np.zeros(512, dtype=np.int64)
    primer_hits = 0
    examined = 0
    total = 0

    with gzip.open(path, "rt") as handle:
        for index, line in enumerate(handle):
            phase = index % 4
            if phase == 1:
                total += 1
                if examined < subsample and primer_pattern.match(line):
                    primer_hits += 1
            elif phase == 3 and examined < subsample:
                q = np.frombuffer(line.rstrip("\n").encode(), dtype=np.uint8) - 33
                n = len(q)
                max_length = max(max_length, n)
                qsum[:n] += q
                qn[:n] += 1
                at_least_25[:n] += (q >= 25)
                at_least_37[:n] += (q >= 37)
                examined += 1

    return {
        "reads": total, "examined": examined, "length": max_length,
        "qsum": qsum[:max_length], "qn": qn[:max_length],
        "ge25": at_least_25[:max_length], "ge37": at_least_37[:max_length],
        "primer_fraction": primer_hits / examined if examined else 0.0,
    }


def best_overlap(forward, reverse_rc, min_length, max_mismatch):
    """Longest suffix of R1 matching a prefix of revcomp(R2), within the mismatch budget."""
    limit = min(len(forward), len(reverse_rc))
    for length in range(limit, min_length - 1, -1):
        a = forward[-length:]
        b = reverse_rc[:length]
        mismatches = sum(1 for x, y in zip(a, b) if x != y)
        if mismatches <= max_mismatch * length:
            return length, mismatches
    return 0, 0


def overlap_test(pairs, settings):
    """Does any real overlap exist? With a 158 bp gap it should not."""
    found = 0
    lengths = []
    for forward, reverse in pairs:
        length, _ = best_overlap(
            forward, revcomp(reverse),
            settings["overlap_min_length"], settings["overlap_max_mismatch_fraction"])
        if length:
            found += 1
            lengths.append(length)
    return found, lengths


def read_pairs(forward_path, reverse_path, limit):
    out = []
    with gzip.open(forward_path, "rt") as f, gzip.open(reverse_path, "rt") as r:
        for index, (fl, rl) in enumerate(zip(f, r)):
            if index % 4 == 1:
                out.append((fl.strip(), rl.strip()))
                if len(out) >= limit:
                    break
    return out


def choose_truncation(fraction_ge, threshold):
    """Longest prefix where the fraction of bases at or above the bin stays above the bar."""
    below = np.where(fraction_ge < threshold)[0]
    return int(below[0]) if below.size else int(len(fraction_ge))


def make_figure(positions, mean_q, low_q, high_q, truncation, counts, groups,
                mode, out_path, config):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    t = THEME[mode]
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 7.2), gridspec_kw={"height_ratios": [1, 1], "hspace": 0.38})
    fig.patch.set_facecolor(t["surface"])

    for ax in (ax1, ax2):
        ax.set_facecolor(t["surface"])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(t["axis"])
            ax.spines[side].set_linewidth(1.0)
        ax.tick_params(colors=t["muted"], labelsize=9, length=3, width=1.0)
        ax.grid(True, color=t["grid"], linewidth=0.8, axis="y")
        ax.set_axisbelow(True)

    # --- panel A: quality by position. One series, so no legend; the title names it.
    ax1.fill_between(positions, low_q, high_q, color=t["series"][0], alpha=0.18,
                     linewidth=0)
    ax1.plot(positions, mean_q, color=t["series"][0], linewidth=2.0)
    ax1.axhline(25, color=t["muted"], linewidth=1.0, linestyle=(0, (4, 3)))
    ax1.text(positions[-1], 25, " Q25", color=t["muted"], fontsize=8,
             va="center", ha="left")
    # Only draw a truncation marker if there is one. Quality holding to the last base is
    # a result, and a line pinned to the right axis would misrepresent it as a cut point.
    if truncation < len(positions):
        ax1.axvline(truncation, color=t["series"][1], linewidth=2.0)
        ax1.annotate(f"truncate at {truncation}", xy=(truncation, ax1.get_ylim()[1]),
                     xytext=(-6, -10), textcoords="offset points",
                     color=t["series"][1], fontsize=9, ha="right", va="top")
    else:
        ax1.text(0.99, 0.06, "quality holds to the last base, no truncation needed",
                 transform=ax1.transAxes, color=t["series"][1], fontsize=9,
                 ha="right", va="bottom")

    ax1.set_title("Forward read quality by position", color=t["ink"], fontsize=11,
                  loc="left", pad=22)
    # Caption sits under the title, not under the axes, where it used to collide with
    # the second panel's title.
    ax1.text(0, 1.045, "line is the mean across all samples; band is the 10th to 90th "
                       "percentile of per-sample means",
             transform=ax1.transAxes, color=t["muted"], fontsize=8, va="bottom")
    ax1.set_xlabel("position in read (bp)", color=t["ink2"], fontsize=9)
    ax1.set_ylabel("mean Phred", color=t["ink2"], fontsize=9)
    ax1.set_xlim(1, positions[-1])

    # --- panel B: depth per sample, two series, so a legend is required.
    order = np.argsort(-counts)
    sorted_counts = counts[order]
    sorted_groups = [groups[i] for i in order]
    labels = sorted(set(groups))
    colours = {label: t["series"][i] for i, label in enumerate(labels)}
    ax2.bar(range(len(sorted_counts)), sorted_counts,
            color=[colours[g] for g in sorted_groups], width=0.78, linewidth=0)
    ax2.axhline(sorted_counts.min(), color=t["muted"], linewidth=1.0,
                linestyle=(0, (4, 3)))
    # In clear space above the descending bars, rather than on top of the short ones.
    ax2.text(0.30, 0.86,
             f"shallowest {sorted_counts.min():,}   median {int(np.median(sorted_counts)):,}"
             f"   deepest {sorted_counts.max():,}",
             transform=ax2.transAxes, color=t["muted"], fontsize=8,
             ha="left", va="center")
    ax2.set_title("Reads per sample", color=t["ink"], fontsize=11, loc="left", pad=10)
    ax2.set_xlabel("samples, deepest to shallowest", color=t["ink2"], fontsize=9)
    ax2.set_ylabel("reads", color=t["ink2"], fontsize=9)
    ax2.set_xlim(-1, len(sorted_counts))
    ax2.set_xticks([])
    handles = [Patch(facecolor=colours[l], label=l) for l in labels]
    legend = ax2.legend(handles=handles, frameon=False, fontsize=9,
                        loc="upper right", ncol=2)
    for text in legend.get_texts():
        text.set_color(t["ink2"])

    fig.savefig(out_path, dpi=200, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--rawdir", default="data/raw", type=Path)
    parser.add_argument("--metadata", default="metadata/sample_metadata.tsv", type=Path)
    parser.add_argument("--outdir", default="results", type=Path)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    settings = config["read_diagnostics"]
    primers = config["primers"]

    rows = [r for r in csv.DictReader(args.metadata.open(), delimiter="\t")
            if not r["sample-id"].startswith("#")]
    if not rows:
        sys.exit(f"no samples in {args.metadata}")
    groups = {r["sample-id"]: r["group"] for r in rows}

    forward_pattern = iupac_regex(primers["forward"]["sequence"])
    reverse_pattern = iupac_regex(primers["reverse"]["sequence"])

    print(f"Scanning {len(rows)} samples in {args.rawdir}")
    print(f"  per-position stats on the first {settings['subsample_reads']:,} reads per "
          f"file; read counts are exact\n")

    per_sample = []
    per_position_means = []      # kept from this pass so the files are read only once
    acc = None
    for index, sample in enumerate(sorted(groups), 1):
        forward_path = args.rawdir / f"{sample}_1.fastq.gz"
        reverse_path = args.rawdir / f"{sample}_2.fastq.gz"
        if not forward_path.exists():
            sys.exit(f"{forward_path} not found. Run 02_download_reads.py first.")

        f = scan_file(forward_path, settings["subsample_reads"], forward_pattern)
        r = scan_file(reverse_path, settings["subsample_reads"], reverse_pattern)

        if acc is None:
            acc = {k: np.zeros_like(f[k]) for k in ("qsum", "qn", "ge25", "ge37")}
        for k in acc:
            acc[k][:len(f[k])] += f[k]
        per_position_means.append(f["qsum"] / np.maximum(1, f["qn"]))

        per_sample.append({
            "sample": sample, "group": groups[sample],
            "reads_forward": f["reads"], "reads_reverse": r["reads"],
            "read_length": f["length"],
            "mean_phred": round(float(f["qsum"].sum() / max(1, f["qn"].sum())), 2),
            "primer_341F_fraction": round(f["primer_fraction"], 4),
            "primer_805R_fraction": round(r["primer_fraction"], 4),
        })
        if index % 10 == 0 or index == len(groups):
            print(f"  {index}/{len(groups)} samples scanned")

    counts = np.array([p["reads_forward"] for p in per_sample])
    group_list = [p["group"] for p in per_sample]

    # ---- failed libraries -----------------------------------------------------------
    # Reported loudly rather than left for someone to notice in a figure label. A sample
    # with a few hundred reads cannot be rarefied to any useful depth, and leaving it in
    # forces the rarefaction depth down to its level, which would discard most of the
    # reads in every other sample.
    minimum = int(config["sample_filtering"]["min_reads"])
    failed = sorted((p for p in per_sample if p["reads_forward"] < minimum),
                    key=lambda p: p["reads_forward"])
    kept = [p for p in per_sample if p["reads_forward"] >= minimum]
    print(f"\nSample depth:")
    print(f"  min {counts.min():,}   median {int(np.median(counts)):,}   "
          f"max {counts.max():,}")
    if failed:
        survivor_min = min(p["reads_forward"] for p in kept)
        print(f"\n  {len(failed)} sample(s) below the {minimum:,} read threshold:")
        for p in failed:
            print(f"    {p['sample']}  {p['reads_forward']:>8,} reads  ({p['group']})")
        print(f"  next shallowest kept: {survivor_min:,}, a "
              f"{survivor_min / max(1, failed[-1]['reads_forward']):.0f}-fold gap, so the "
              f"threshold is not\n  sensitive to where it is drawn")
        from collections import Counter as _C
        print(f"  group sizes after exclusion: "
              f"{dict(_C(p['group'] for p in kept))}")
        print(f"  usable rarefaction depth becomes {survivor_min:,}, keeping "
              f"{len(kept)} samples")

    # ---- quality summary and the truncation choice ----------------------------------
    mean_by_position = acc["qsum"] / np.maximum(1, acc["qn"])
    fraction_ge25 = acc["ge25"] / np.maximum(1, acc["qn"])
    fraction_ge37 = acc["ge37"] / np.maximum(1, acc["qn"])
    positions = np.arange(1, len(mean_by_position) + 1)

    threshold = float(settings["truncation_fraction_threshold"])
    chosen = choose_truncation(fraction_ge25, threshold)

    print(f"\nQuality by position, pooled over all samples:")
    for p in (1, 25, 50, 75, 100, 125, len(mean_by_position)):
        if p <= len(mean_by_position):
            i = p - 1
            print(f"  position {p:>3}: mean Q {mean_by_position[i]:5.1f}   "
                  f"Q25+ {100 * fraction_ge25[i]:5.1f}%   "
                  f"Q37 {100 * fraction_ge37[i]:5.1f}%")

    print(f"\nTruncation, by the rule in config "
          f"(longest prefix with >= {100 * threshold:.0f}% of bases at Q25+):")
    print(f"  chosen: {chosen} bp")
    for alternative in (0.90, 0.70):
        print(f"  for comparison, at >= {100 * alternative:.0f}%: "
              f"{choose_truncation(fraction_ge25, alternative)} bp")

    # ---- primers --------------------------------------------------------------------
    f_fraction = np.mean([p["primer_341F_fraction"] for p in per_sample])
    r_fraction = np.mean([p["primer_805R_fraction"] for p in per_sample])
    print(f"\nPrimers at the 5' end, mean across samples:")
    print(f"  {primers['forward']['name']} on R1: {100 * f_fraction:.1f}%")
    print(f"  {primers['reverse']['name']} on R2: {100 * r_fraction:.1f}%")
    if f_fraction > 0.5:
        print(f"  present, so they must be trimmed before denoising. Degenerate positions"
              f"\n  mean one organism would otherwise split across several ASVs.")

    # ---- the overlap test -----------------------------------------------------------
    print(f"\nOverlap test on {settings['overlap_test_samples']} samples x "
          f"{settings['overlap_test_pairs']} pairs")
    tested = 0
    found_total = 0
    all_lengths = []
    for sample in sorted(groups)[:int(settings["overlap_test_samples"])]:
        pairs = read_pairs(args.rawdir / f"{sample}_1.fastq.gz",
                           args.rawdir / f"{sample}_2.fastq.gz",
                           int(settings["overlap_test_pairs"]))
        found, lengths = overlap_test(pairs, settings)
        tested += len(pairs)
        found_total += found
        all_lengths.extend(lengths)
        print(f"  {sample}: {found}/{len(pairs)} pairs with a >= "
              f"{settings['overlap_min_length']} bp overlap")
    print(f"  total: {found_total}/{tested} "
          f"({100 * found_total / max(1, tested):.1f}%)")
    if found_total / max(1, tested) < 0.05:
        print(f"  essentially no overlap, as the amplicon length predicted. Merging is")
        print(f"  impossible; forward reads only. This is now shown, not assumed.")

    # ---- write ----------------------------------------------------------------------
    tables = args.outdir / "tables"
    figures = args.outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    sample_table = tables / "read_diagnostics_by_sample.tsv"
    with sample_table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_sample[0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(per_sample)

    position_table = tables / "quality_by_position.tsv"
    with position_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["position", "mean_phred", "fraction_ge_q25", "fraction_ge_q37",
                         "bases_counted"])
        for i, p in enumerate(positions):
            writer.writerow([int(p), round(float(mean_by_position[i]), 3),
                             round(float(fraction_ge25[i]), 5),
                             round(float(fraction_ge37[i]), 5), int(acc["qn"][i])])

    # Per-sample spread for the band in panel A, collected during the single scan above.
    stack = np.vstack([m[:len(mean_by_position)] for m in per_position_means])
    low = np.percentile(stack, 10, axis=0)
    high = np.percentile(stack, 90, axis=0)

    for mode in ("light", "dark"):
        out_path = figures / f"read_quality_{mode}.png"
        make_figure(positions, mean_by_position, low, high, chosen,
                    counts, group_list, mode, out_path, config)
        print(f"\nWrote {out_path}")

    print(f"Wrote {sample_table}")
    print(f"Wrote {position_table}")

    # The truncation figure above is measured on untrimmed reads. Primer removal comes
    # first in the pipeline, so what reaches DADA2 is shorter by the primer length and
    # the parameter has to be stated against that, not against the raw read.
    primer_length = len(primers["forward"]["sequence"])
    effective = chosen - primer_length
    print(f"\nParameters for denoising:")
    if chosen >= len(mean_by_position):
        print(f"  quality holds to the last base, so nothing is lost to truncation")
    else:
        print(f"  quality rule truncates the raw read at {chosen} of "
              f"{len(mean_by_position)} bp")
    print(f"  but cutadapt strips {primer_length} bp of "
          f"{primers['forward']['name']} first, so reads reach DADA2 at "
          f"{len(mean_by_position) - primer_length} bp")
    print(f"  --p-trunc-len {effective}      (not {chosen}; that would be measured "
          f"against the untrimmed read)")
    print(f"  --p-trim-left 0        (cutadapt has already removed the primer)")


if __name__ == "__main__":
    main()
