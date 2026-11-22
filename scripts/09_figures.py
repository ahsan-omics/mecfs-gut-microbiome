#!/usr/bin/env python3
"""Figures for the README: ordination, alpha diversity, differential abundance.

Each is written in a light and a dark version so it reads on GitHub either way. The
palette is the one validated with the data-viz validator: all six checks pass in both
modes, and identity is always carried by a legend as well as by colour.

Deliberate choices worth knowing:

  The ordination uses weighted UniFrac, because that is one of the two distance measures
  whose PERMANOVA was not undermined by a dispersion difference. Plotting Bray-Curtis or
  Jaccard as the headline would show a separation that PERMDISP says is partly spread.

  The volcano shows all 172 tested genera, not just the two that reached significance.
  A figure of only the hits would imply a stronger result than 2 of 172 is.
"""

import argparse
import csv
import io
import math
import zipfile
from pathlib import Path

import numpy as np
import yaml

THEME = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e",
              "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7",
              "series": ["#2a78d6", "#eb6834"], "null": "#b8b7b0"},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7",
             "muted": "#898781", "grid": "#2c2c2a", "axis": "#383835",
             "series": ["#3987e5", "#d95926"], "null": "#55544f"},
}
GROUPS = ["control", "mecfs"]
LABELS = {"control": "control", "mecfs": "ME/CFS"}


def style(ax, t, grid_axis="y"):
    ax.set_facecolor(t["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["axis"])
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(colors=t["muted"], labelsize=9, length=3, width=1.0)
    if grid_axis:
        ax.grid(True, color=t["grid"], linewidth=0.8, axis=grid_axis)
    ax.set_axisbelow(True)


def read_ordination(qza_path):
    """skbio OrdinationResults text: a Site block of ids plus coordinates."""
    with zipfile.ZipFile(qza_path) as archive:
        name = [n for n in archive.namelist() if n.endswith("ordination.txt")][0]
        lines = archive.read(name).decode().splitlines()
    explained, coords, mode = [], {}, None
    for line in lines:
        if line.startswith("Proportion explained"):
            mode = "explained"
            continue
        if line.startswith("Site\t"):
            mode = "site"
            continue
        if line.startswith(("Eigvals", "Species", "Biplot", "Site constraints")):
            mode = None
            continue
        if not line.strip():
            continue
        if mode == "explained":
            explained = [float(x) for x in line.split("\t")]
            mode = None
        elif mode == "site":
            parts = line.split("\t")
            coords[parts[0]] = [float(x) for x in parts[1:]]
    return coords, explained


def read_slices(qza_path):
    with zipfile.ZipFile(qza_path) as archive:
        out = {}
        for name in archive.namelist():
            base = Path(name).name
            if base.endswith("_slice.csv"):
                out[base.replace("_slice.csv", "")] = list(csv.DictReader(
                    io.StringIO(archive.read(name).decode())))
    return out


def kruskal_p(qzv_path, column="group"):
    try:
        with zipfile.ZipFile(qzv_path) as archive:
            name = [n for n in archive.namelist()
                    if n.endswith(f"kruskal-wallis-pairwise-{column}.csv")][0]
            rows = list(csv.DictReader(io.StringIO(archive.read(name).decode())))
        return float(rows[0]["p-value"])
    except Exception:
        return None


def describe_taxon(lineage, wrap=False):
    """Name a taxon for a figure.

    'Incertae Sedis' is Latin for "placement uncertain" and carries no information for a
    reader, so it is dropped from the displayed name rather than eating half a label.
    """
    parts = [p.strip() for p in lineage.split(";") if p.strip()]
    useful = [p for p in parts if not p.endswith("__")
              and p.split("__", 1)[-1].lower() not in
              ("uncultured", "unidentified", "", "metagenome")]
    tip = parts[-1].split("__", 1)[-1] if parts else lineage
    if tip.lower() in ("uncultured", "unidentified", "metagenome", ""):
        parent = useful[-1].split("__", 1)[-1] if useful else "?"
        parent = parent.replace("_Incertae_Sedis", "").replace("_", " ").strip()
        return f"uncultured\n{parent}" if wrap else f"uncultured {parent}"
    return tip.replace("_", " ")


def boxes(ax, values_by_group, t, width=0.5):
    """A thin box plot with the points drawn over it. n is small enough to show them."""
    rng = np.random.default_rng(11)
    for index, group in enumerate(GROUPS):
        values = np.array(values_by_group[group], dtype=float)
        if values.size == 0:
            continue
        colour = t["series"][index]
        q1, med, q3 = np.percentile(values, [25, 50, 75])
        whisker_low = values[values >= q1 - 1.5 * (q3 - q1)].min()
        whisker_high = values[values <= q3 + 1.5 * (q3 - q1)].max()
        ax.plot([index, index], [whisker_low, whisker_high], color=colour,
                linewidth=1.2, zorder=2, solid_capstyle="butt")
        ax.add_patch(__import__("matplotlib").patches.Rectangle(
            (index - width / 2, q1), width, q3 - q1, facecolor=t["surface"],
            edgecolor=colour, linewidth=1.6, zorder=3))
        ax.plot([index - width / 2, index + width / 2], [med, med], color=colour,
                linewidth=2.2, zorder=4, solid_capstyle="butt")
        jitter = index + (rng.random(values.size) - 0.5) * (width * 0.72)
        ax.scatter(jitter, values, s=16, color=colour, alpha=0.45, linewidths=0,
                   zorder=5)


def figure_pcoa(coords, explained, groups, permanova, mode, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    t = THEME[mode]
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    fig.patch.set_facecolor(t["surface"])
    style(ax, t, grid_axis="both")

    for index, group in enumerate(GROUPS):
        samples = [s for s in coords if groups.get(s) == group]
        xs = [coords[s][0] for s in samples]
        ys = [coords[s][1] for s in samples]
        ax.scatter(xs, ys, s=58, color=t["series"][index], alpha=0.75,
                   edgecolors=t["surface"], linewidths=1.4, zorder=3)
        # Centroid, marked with a ring so it reads as a summary not another sample.
        ax.scatter([np.mean(xs)], [np.mean(ys)], s=190, color=t["series"][index],
                   edgecolors=t["surface"], linewidths=2.4, marker="D", zorder=4)

    ax.axhline(0, color=t["grid"], linewidth=0.8, zorder=1)
    ax.axvline(0, color=t["grid"], linewidth=0.8, zorder=1)
    ax.set_title("Weighted UniFrac PCoA", color=t["ink"], fontsize=12,
                 loc="left", pad=22)
    ax.text(0, 1.035, "diamonds are group centroids; dispersion does not differ "
                      "between groups (PERMDISP p = 0.196)",
            transform=ax.transAxes, color=t["muted"], fontsize=8, va="bottom")
    ax.set_xlabel(f"PC1  ({100 * explained[0]:.1f}% of variance)",
                  color=t["ink2"], fontsize=9)
    ax.set_ylabel(f"PC2  ({100 * explained[1]:.1f}%)", color=t["ink2"], fontsize=9)

    ax.text(0.98, 0.02, permanova, transform=ax.transAxes, color=t["ink2"],
            fontsize=9, ha="right", va="bottom")
    handles = [Line2D([], [], marker="o", linestyle="", markersize=8,
                      markerfacecolor=t["series"][i], markeredgecolor=t["surface"],
                      label=LABELS[g]) for i, g in enumerate(GROUPS)]
    legend = ax.legend(handles=handles, frameon=False, fontsize=9, loc="upper right")
    for text in legend.get_texts():
        text.set_color(t["ink2"])

    fig.savefig(out_path, dpi=200, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def figure_alpha(alpha, groups, pvalues, metrics, mode, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    t = THEME[mode]
    fig, axes = plt.subplots(1, len(metrics), figsize=(2.6 * len(metrics) + 1.2, 4.6))
    fig.patch.set_facecolor(t["surface"])

    pretty = {"shannon": "Shannon", "observed_features": "Observed features",
              "faith_pd": "Faith's PD", "evenness": "Pielou evenness"}

    for ax, metric in zip(np.atleast_1d(axes), metrics):
        style(ax, t)
        values = {g: [float(r[metric]) for r in alpha if r["group"] == g
                      and r.get(metric) not in (None, "")] for g in GROUPS}
        boxes(ax, values, t)
        p = pvalues.get(metric)
        tag = "n.s." if (p is None or p >= 0.05) else f"p = {p:.3f}"
        ax.set_title(pretty.get(metric, metric), color=t["ink"], fontsize=10,
                     loc="left", pad=16)
        ax.text(0, 1.02, tag, transform=ax.transAxes, fontsize=9, va="bottom",
                color=t["ink2"] if (p and p < 0.05) else t["muted"])
        ax.set_xticks(range(len(GROUPS)))
        ax.set_xticklabels([LABELS[g] for g in GROUPS], color=t["ink2"], fontsize=9)
        ax.set_xlim(-0.6, len(GROUPS) - 0.4)

    handles = [Patch(facecolor=t["series"][i], label=LABELS[g])
               for i, g in enumerate(GROUPS)]
    legend = fig.legend(handles=handles, frameon=False, fontsize=9,
                        loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.04))
    for text in legend.get_texts():
        text.set_color(t["ink2"])
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_path, dpi=200, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def figure_da(lfc, qval, term, abundances, groups, mode, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    t = THEME[mode]
    taxa = [k for k in qval if k in lfc]
    significant = sorted([k for k in taxa if qval[k] < 0.05], key=lambda k: -lfc[k])
    others = [k for k in taxa if qval[k] >= 0.05]

    # One abundance panel per significant genus, each with its own y axis. Bacteroides
    # sits near 20 percent and the uncultured genus near 0.2, so a shared linear axis
    # would flatten the second one into the baseline and hide the finding.
    fig, axes = plt.subplots(
        1, 1 + len(significant), figsize=(6.4 + 2.5 * len(significant), 5.2),
        gridspec_kw={"width_ratios": [1.7] + [0.75] * len(significant)})
    ax1 = axes[0]
    abundance_axes = list(axes[1:])
    fig.patch.set_facecolor(t["surface"])
    style(ax1, t, grid_axis="both")
    for ax in abundance_axes:
        style(ax, t)

    def y(k):
        return -math.log10(max(qval[k], 1e-12))

    ax1.scatter([lfc[k] for k in others], [y(k) for k in others], s=22,
                color=t["null"], alpha=0.75, linewidths=0, zorder=2)
    for k in significant:
        colour = t["series"][1] if lfc[k] > 0 else t["series"][0]
        ax1.scatter([lfc[k]], [y(k)], s=78, color=colour,
                    edgecolors=t["surface"], linewidths=1.4, zorder=4)
        # Label toward the middle of the plot, so it cannot run off the edge.
        inward = -1 if lfc[k] > 0 else 1
        ax1.annotate(describe_taxon(k), xy=(lfc[k], y(k)),
                     xytext=(10 * inward, 10), textcoords="offset points",
                     fontsize=9, color=t["ink2"],
                     ha="right" if inward < 0 else "left")

    ax1.axhline(-math.log10(0.05), color=t["muted"], linewidth=1.0,
                linestyle=(0, (4, 3)), zorder=1)
    ax1.text(ax1.get_xlim()[0], -math.log10(0.05), " q = 0.05", color=t["muted"],
             fontsize=8, va="bottom", ha="left")
    ax1.axvline(0, color=t["grid"], linewidth=0.8, zorder=1)
    ax1.set_title(f"ANCOM-BC, {len(taxa)} genera tested", color=t["ink"],
                  fontsize=11, loc="left", pad=22)
    ax1.text(0, 1.035, f"adjusted for sex, age tertile and BMI tertile; "
                       f"{len(significant)} reach q < 0.05",
             transform=ax1.transAxes, color=t["muted"], fontsize=8, va="bottom")
    ax1.set_xlabel("log fold change  (negative: higher in control)",
                   color=t["ink2"], fontsize=9)
    ax1.set_ylabel("-log10 q", color=t["ink2"], fontsize=9)

    # One panel per hit, each on its own scale.
    for ax, taxon in zip(abundance_axes, significant):
        values = {g: np.array([100 * v for s, v in abundances.get(taxon, {}).items()
                               if groups.get(s) == g], dtype=float) for g in GROUPS}
        boxes(ax, values, t, width=0.46)
        ax.set_xticks(range(len(GROUPS)))
        ax.set_xticklabels([LABELS[g] for g in GROUPS], color=t["ink2"], fontsize=9)
        ax.set_xlim(-0.62, len(GROUPS) - 0.38)
        ax.set_title(describe_taxon(taxon, wrap=True), color=t["ink"], fontsize=10,
                     loc="left", pad=16)
        ax.set_ylabel("relative abundance (%)", color=t["ink2"], fontsize=9)
        ax.text(0, 1.02, f"q = {qval[taxon]:.3f}", transform=ax.transAxes,
                fontsize=9, color=t["ink2"], va="bottom")

    handles = [Patch(facecolor=t["series"][i], label=LABELS[g])
               for i, g in enumerate(GROUPS)]
    legend = fig.legend(handles=handles, frameon=False, fontsize=9,
                        loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.03))
    for text in legend.get_texts():
        text.set_color(t["ink2"])

    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(out_path, dpi=200, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--qiimedir", default="results/qiime", type=Path)
    parser.add_argument("--outdir", default="results/figures", type=Path)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    metrics = config["diversity"]["alpha_metrics"]
    args.outdir.mkdir(parents=True, exist_ok=True)

    alpha = [r for r in csv.DictReader(
        Path("results/tables/alpha_diversity.tsv").open(), delimiter="\t")]
    groups = {r["sample-id"]: r["group"] for r in alpha}

    coords, explained = read_ordination(
        args.qiimedir / "core-metrics" / "weighted_unifrac_pcoa_results.qza")
    pvalues = {m: kruskal_p(args.qiimedir / "stats" / f"alpha-{m}.qzv")
               for m in metrics}

    slices = read_slices(args.qiimedir / "stats" / "ancombc-adjusted.qza")
    key = list(slices["lfc"][0])[0]
    term = next(c for c in slices["lfc"][0]
                if c != key and c.lower().startswith("group"))
    lfc = {r[key]: float(r[term]) for r in slices["lfc"] if r.get(term)}
    qval = {r[key]: float(r[term]) for r in slices["q_val"] if r.get(term)}

    # Relative abundance of each significant genus, per sample.
    import biom
    import tempfile
    with zipfile.ZipFile(args.qiimedir / "genus-table-filtered.qza") as archive:
        name = [n for n in archive.namelist() if n.endswith(".biom")][0]
        with tempfile.NamedTemporaryFile(suffix=".biom", delete=False) as handle:
            handle.write(archive.read(name))
            biom_path = handle.name
    table = biom.load_table(biom_path)
    totals = dict(zip(table.ids("sample"), table.sum("sample")))
    abundances = {}
    for taxon in [k for k in qval if qval[k] < 0.05]:
        if taxon in set(table.ids("observation")):
            row = table.data(taxon, axis="observation")
            abundances[taxon] = {s: row[i] / max(1.0, totals[s])
                                 for i, s in enumerate(table.ids("sample"))}

    permanova = "PERMANOVA  p = 0.021,  $R^2$ = 4.2%"
    for mode in ("light", "dark"):
        figure_pcoa(coords, explained, groups, permanova, mode,
                    args.outdir / f"pcoa_weighted_unifrac_{mode}.png")
        figure_alpha(alpha, groups, pvalues, metrics, mode,
                     args.outdir / f"alpha_diversity_{mode}.png")
        figure_da(lfc, qval, term, abundances, groups, mode,
                  args.outdir / f"differential_abundance_{mode}.png")
        print(f"wrote 3 figures ({mode})")

    print(f"\nAll figures in {args.outdir}")


if __name__ == "__main__":
    main()
