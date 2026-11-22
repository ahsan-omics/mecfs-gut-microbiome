#!/usr/bin/env python3
"""The case against control tests, all through QIIME 2 plugins.

  alpha    Kruskal-Wallis per metric, via alpha-group-significance.

  beta     PERMANOVA and PERMDISP on each distance matrix. Both, always, because
           PERMANOVA tests whether group centroids differ and PERMDISP tests whether the
           groups differ in spread. PERMANOVA is sensitive to dispersion differences, so
           a significant PERMANOVA with a significant PERMDISP may mean the ME/CFS
           samples are simply more variable than the controls rather than centred
           somewhere else. Reporting one without the other is how that gets missed.

  beta,    adonis with age, sex and BMI in the model. vegan's adonis2 is sequential, so
  adjusted term order decides what each term is adjusted for. Run both ways on purpose.

  taxa     ANCOM-BC on genus-collapsed counts, with the same covariates. On the
           UNRAREFIED table: it models compositionality and library size itself.

No statistical code is written here. These are the standard implementations, which is
also what keeps anything from the private BiomeDelta pipeline out of a public repo.
"""

import argparse
import csv
import html
import io
import re
import shutil
import subprocess
import sys
import zipfile
from collections import Counter
from pathlib import Path

import yaml


def run(command, description, allow_fail=False):
    print(f"\n>>> {description}")
    print("    " + " ".join(str(c) for c in command))
    result = subprocess.run([str(c) for c in command], capture_output=True, text=True)
    if result.returncode != 0:
        message = (result.stderr or result.stdout).strip()[-2500:]
        if allow_fail:
            print(f"    SKIPPED (exit {result.returncode}): {message.splitlines()[-1:]}",
                  file=sys.stderr)
            return None
        print(f"\nCommand failed with exit code {result.returncode}:", file=sys.stderr)
        print(message, file=sys.stderr)
        sys.exit(1)
    return result.stdout


def qzv_tables(qzv_path):
    """Strip the tables out of a visualization's index.html as rows of text."""
    try:
        with zipfile.ZipFile(qzv_path) as archive:
            match = [n for n in archive.namelist() if n.endswith("index.html")]
            if not match:
                return []
            page = archive.read(match[0]).decode("utf-8", "replace")
    except Exception:
        return []
    rows = []
    for row_html in re.findall(r"<tr[^>]*>(.*?)</tr>", page, re.S | re.I):
        cells = [html.unescape(re.sub(r"<[^>]+>", "", c)).strip()
                 for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row_html, re.S | re.I)]
        cells = [c for c in cells if c != ""]
        if cells:
            rows.append(cells)
    return rows


def report_qzv(qzv_path, label):
    rows = qzv_tables(qzv_path)
    if not rows:
        print(f"  {label}: open {qzv_path} at view.qiime2.org")
        return
    print(f"  {label}:")
    for cells in rows[:8]:
        print("    " + "  ".join(cells))


def report_kruskal_wallis(qzv_path, column, label):
    """The alpha visualization renders its table in JavaScript, so scraping the HTML
    returns only column headers. The numbers live in a CSV inside the artifact."""
    try:
        with zipfile.ZipFile(qzv_path) as archive:
            match = [n for n in archive.namelist()
                     if n.endswith(f"kruskal-wallis-pairwise-{column}.csv")]
            if not match:
                return report_qzv(qzv_path, label)
            rows = list(csv.DictReader(
                io.StringIO(archive.read(match[0]).decode())))
    except Exception:
        return report_qzv(qzv_path, label)
    for row in rows:
        h = row.get("H", "")
        p = row.get("p-value", "")
        q = row.get("q-value", "")
        try:
            print(f"  {label:<20} H {float(h):7.3f}   p {float(p):.4f}   "
                  f"q {float(q):.4f}")
        except (TypeError, ValueError):
            print(f"  {label}: {row}")


def read_differentials(qza_path, column):
    """Significant taxa from an ANCOM-BC artifact.

    The data files are named lfc_slice.csv, q_val_slice.csv and so on, and they are
    comma separated despite the .qza convention elsewhere.
    """
    try:
        with zipfile.ZipFile(qza_path) as archive:
            slices = {}
            for name in archive.namelist():
                base = Path(name).name
                if base.endswith("_slice.csv"):
                    slices[base.replace("_slice.csv", "")] = list(csv.DictReader(
                        io.StringIO(archive.read(name).decode())))
    except Exception:
        return [], 0, None
    if "lfc" not in slices or "q_val" not in slices:
        return [], 0, None

    key = list(slices["lfc"][0])[0]
    term = next((c for c in slices["lfc"][0]
                 if c != key and c.lower().startswith(column.lower())), None)
    if term is None:
        return [], 0, None

    q = {r[key]: float(r[term]) for r in slices["q_val"] if r.get(term) not in (None, "")}
    lfc = {r[key]: float(r[term]) for r in slices["lfc"] if r.get(term) not in (None, "")}
    hits = sorted(((t, lfc[t], q[t]) for t in q if q[t] < 0.05), key=lambda x: x[2])
    return hits, len(q), term


def describe_taxon(lineage):
    """Name a taxon usefully.

    SILVA uses 'uncultured' as a genus placeholder, and taxa collapse keys on the full
    lineage, so several distinct 'g__uncultured' features exist under different families.
    Printing the last field alone would label them identically and wrongly imply one
    organism, so an uninformative tip is reported against its nearest named ancestor.
    """
    parts = [p.strip() for p in lineage.split(";") if p.strip()]
    useful = [p for p in parts
              if not p.endswith("__") and p.split("__", 1)[-1].lower() not in
              ("uncultured", "unidentified", "", "metagenome")]
    tip = parts[-1].split("__", 1)[-1] if parts else lineage
    if tip.lower() in ("uncultured", "unidentified", "metagenome", ""):
        parent = useful[-1].split("__", 1)[-1] if useful else "?"
        return f"uncultured genus in {parent}"
    return tip


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--qiimedir", default="results/qiime", type=Path)
    parser.add_argument("--metadata",
                        default="metadata/sample_metadata_filtered.tsv", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    stats = config["statistics"]
    div = config["diversity"]
    column = stats["group_column"]

    core = args.qiimedir / "core-metrics"
    table = args.qiimedir / "table-filtered.qza"
    taxonomy = args.qiimedir / "taxonomy.qza"
    for required in (core, table, taxonomy, args.metadata):
        if not required.exists():
            sys.exit(f"{required} not found. Run the earlier scripts first.")

    out = args.qiimedir / "stats"
    out.mkdir(parents=True, exist_ok=True)

    # ---- alpha ----------------------------------------------------------------------
    print("=" * 70)
    print("ALPHA DIVERSITY  (Kruskal-Wallis)")
    print("=" * 70)
    for metric in div["alpha_metrics"]:
        vector = core / f"{metric}_vector.qza"
        if not vector.exists():
            continue
        viz = out / f"alpha-{metric}.qzv"
        if args.force or not viz.exists():
            run(["qiime", "diversity", "alpha-group-significance",
                 "--i-alpha-diversity", vector,
                 "--m-metadata-file", args.metadata,
                 "--o-visualization", viz],
                f"{metric}")
        report_kruskal_wallis(viz, column, metric)

    # ---- beta -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("BETA DIVERSITY  (PERMANOVA, then PERMDISP on the same matrix)")
    print("=" * 70)
    for metric in div["beta_metrics"]:
        matrix = core / f"{metric}_distance_matrix.qza"
        if not matrix.exists():
            continue
        print(f"\n--- {metric} ---")
        for method in ("permanova", "permdisp"):
            viz = out / f"beta-{metric}-{method}.qzv"
            if args.force or not viz.exists():
                run(["qiime", "diversity", "beta-group-significance",
                     "--i-distance-matrix", matrix,
                     "--m-metadata-file", args.metadata,
                     "--m-metadata-column", column,
                     "--p-method", method,
                     "--p-permutations", stats["permutations"],
                     "--o-visualization", viz],
                    f"{metric} {method}")
            report_qzv(viz, method)

    # ---- beta, adjusted --------------------------------------------------------------
    print("\n" + "=" * 70)
    print("BETA DIVERSITY, ADJUSTED  (adonis; term order matters, so both are run)")
    print("=" * 70)
    for metric in div["beta_metrics"]:
        matrix = core / f"{metric}_distance_matrix.qza"
        if not matrix.exists():
            continue
        print(f"\n--- {metric} ---")
        for label, formula in (("group last (adjusted)", stats["adonis_formula_adjusted"]),
                               ("group first", stats["adonis_formula_group_first"])):
            viz = out / f"adonis-{metric}-{formula.replace('+', '_')}.qzv"
            if args.force or not viz.exists():
                run(["qiime", "diversity", "adonis",
                     "--i-distance-matrix", matrix,
                     "--m-metadata-file", args.metadata,
                     "--p-formula", formula,
                     "--p-permutations", stats["permutations"],
                     "--o-visualization", viz],
                    f"{metric}: {formula}", allow_fail=True)
            if viz.exists():
                report_qzv(viz, f"{label}  [{formula}]")

    # ---- differential abundance -------------------------------------------------------
    print("\n" + "=" * 70)
    print("DIFFERENTIAL ABUNDANCE  (ANCOM-BC, genus level, unrarefied)")
    print("=" * 70)

    genus = args.qiimedir / "genus-table.qza"
    if args.force or not genus.exists():
        run(["qiime", "taxa", "collapse",
             "--i-table", table,
             "--i-taxonomy", taxonomy,
             "--p-level", stats["taxa_level"],
             "--o-collapsed-table", genus],
            f"collapse to level {stats['taxa_level']} (genus)")

    n_samples = sum(1 for r in csv.DictReader(args.metadata.open(), delimiter="\t")
                    if not r["sample-id"].startswith("#"))
    min_samples = max(2, int(round(stats["min_samples_fraction"] * n_samples)))
    genus_filtered = args.qiimedir / "genus-table-filtered.qza"
    if args.force or not genus_filtered.exists():
        run(["qiime", "feature-table", "filter-features",
             "--i-table", genus,
             "--p-min-samples", min_samples,
             "--o-filtered-table", genus_filtered],
            f"keep genera present in at least {min_samples} of {n_samples} samples")

    # ANCOM-BC cannot take numeric covariates in this plugin version, so age and BMI are
    # adjusted for as tertiles. The reasoning and the exact bug are in config.yaml.
    ancombc_metadata = Path("metadata/sample_metadata_ancombc.tsv")
    rows = [r for r in csv.DictReader(args.metadata.open(), delimiter="\t")
            if not r["sample-id"].startswith("#")]
    columns = list(rows[0])
    derived = []
    for source in stats["tertile_columns"]:
        values = sorted(float(r[source]) for r in rows)
        low, high = values[len(values) // 3], values[2 * len(values) // 3]
        name = f"{source}_tertile"
        derived.append(name)
        for row in rows:
            v = float(row[source])
            row[name] = (f"{source}_low" if v < low else
                         f"{source}_mid" if v < high else f"{source}_high")
        sizes = Counter(r[name] for r in rows)
        print(f"\n  {name}: cuts at {low:g} and {high:g}, sizes {dict(sizes)}")

    with ancombc_metadata.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out_columns = columns + derived
        writer.writerow(out_columns)
        # Everything categorical: the numeric originals are not used in these formulas,
        # and declaring them numeric is what triggers the plugin bug if they ever are.
        writer.writerow(["#q2:types"] + ["categorical"] * (len(out_columns) - 1))
        for row in rows:
            writer.writerow([row[c] for c in out_columns])
    print(f"  wrote {ancombc_metadata}")

    models = [("unadjusted", stats["ancombc_formula_unadjusted"]),
              ("adjusted", stats["ancombc_formula_adjusted"])]
    results = {}

    for label, formula in models:
        differentials = out / f"ancombc-{label}.qza"
        if args.force or not differentials.exists():
            run(["qiime", "composition", "ancombc",
                 "--i-table", genus_filtered,
                 "--m-metadata-file", ancombc_metadata,
                 "--p-formula", formula,
                 "--o-differentials", differentials],
                f"ANCOM-BC {label}: '{formula}'")
        viz = out / f"ancombc-{label}.qzv"
        if args.force or not viz.exists():
            run(["qiime", "composition", "tabulate",
                 "--i-data", differentials, "--o-visualization", viz],
                f"tabulate {label}", allow_fail=True)
        results[label] = read_differentials(differentials, column)

    # ---- report ---------------------------------------------------------------------
    for label, _ in models:
        hits, tested, term = results.get(label, ([], 0, None))
        print(f"\n  {label.upper()}  term '{term}', {tested} genera tested, "
              f"{len(hits)} significant at q < 0.05")
        for lineage, lfc, q in hits[:20]:
            direction = "higher in ME/CFS" if lfc > 0 else "higher in control"
            print(f"    {describe_taxon(lineage):<48} lfc {lfc:+.2f}  q {q:.4f}  "
                  f"{direction}")
        if not hits and tested:
            print(f"    none. With 61 samples at genus resolution that is a plausible")
            print(f"    outcome, not a failure. Report it as it is.")

    both = ({t for t, _, _ in results.get("unadjusted", ([],))[0]}
            & {t for t, _, _ in results.get("adjusted", ([],))[0]})
    if both:
        print(f"\n  Significant in BOTH models: {len(both)}. Agreement between the two is"
              f"\n  what a design already matched on age, sex and BMI should produce.")

    da_table = Path("results/tables/differential_abundance.tsv")
    with da_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["model", "taxon", "lfc", "q_value", "direction"])
        for label, _ in models:
            for lineage, lfc, q in results.get(label, ([],))[0]:
                writer.writerow([label, lineage, f"{lfc:.4f}", f"{q:.5f}",
                                 "mecfs" if lfc > 0 else "control"])
    print(f"\nWrote {da_table}")

    print(f"\nAll visualizations in {out}")
    print(f"\nHow to read this: PERMANOVA says whether the groups sit in different places,")
    print(f"PERMDISP whether they differ in spread. A significant PERMANOVA alongside a")
    print(f"significant PERMDISP is weaker evidence of a location difference than it looks.")


if __name__ == "__main__":
    main()
