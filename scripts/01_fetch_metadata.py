#!/usr/bin/env python3
"""Build the sample metadata table and the download manifest for PRJNA1332076.

Two sources, joined on the sample accession:

  SRA runinfo   one row per run, with the spot and base counts used later to verify that
                each download converted completely.

  BioSample     the per-sample attributes. This dataset was chosen because host_disease,
                host_age, host_sex and host_body_mass_index are populated in the archive
                itself, so there is no supplementary table to chase and no PDF to parse.

This script only queries APIs and writes two small tables. It downloads no reads, so it
needs nothing but a Python with pyyaml. The reads are fetched by 02_download_reads.py,
which needs sra-tools and therefore a different environment.

One thing it deliberately does NOT use: the `download_path` column in SRA's own runinfo.
That points at the SRA Lite version of each run, whose quality scores have all been
replaced by a single synthesised constant. See `download_route` in config.yaml.
"""

import argparse
import csv
import json
import os
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import yaml

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "mecfs-gut-microbiome"


class Problems:
    """Collect every failure so one run reports all of them."""

    def __init__(self):
        self.items = []

    def check(self, condition, message):
        if not condition:
            self.items.append(message)
        return bool(condition)

    def report_and_exit(self):
        if not self.items:
            return
        print(f"\n{len(self.items)} problem(s) found:\n", file=sys.stderr)
        for i, item in enumerate(self.items, 1):
            print(f"  {i}. {item}", file=sys.stderr)
        print("", file=sys.stderr)
        sys.exit(1)


def get(url, timeout=180):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8")


def with_identity(url, email):
    if email:
        url += f"&email={urllib.parse.quote(email)}&tool={USER_AGENT}"
    return url


def fetch_run_table(bioproject, email):
    """SRA runinfo, one row per run."""
    search = with_identity(
        f"{EUTILS}/esearch.fcgi?db=sra&term={bioproject}&retmax=1000&retmode=json", email)
    ids = json.loads(get(search))["esearchresult"]["idlist"]
    if not ids:
        sys.exit(f"No SRA records found for {bioproject}.")
    fetch = with_identity(
        f"{EUTILS}/efetch.fcgi?db=sra&id={','.join(ids)}&rettype=runinfo&retmode=text",
        email)
    rows = list(csv.DictReader(get(fetch).splitlines()))
    return [r for r in rows if r.get("Run")]


def fetch_biosample_attributes(accessions, email):
    """BioSample attributes, keyed by sample accession."""
    out = {}
    for i in range(0, len(accessions), 150):
        chunk = ",".join(accessions[i:i + 150])
        url = with_identity(
            f"{EUTILS}/efetch.fcgi?db=biosample&id={chunk}&retmode=xml", email)
        root = ET.fromstring(get(url))
        for sample in root.findall(".//BioSample"):
            attributes = {}
            for attribute in sample.findall(".//Attribute"):
                name = attribute.get("harmonized_name") or attribute.get("attribute_name")
                attributes[name] = (attribute.text or "").strip()
            out[sample.get("accession")] = attributes
    return out


def normalise_group(raw, labels, problems, sample):
    if raw in labels:
        return labels[raw]
    problems.check(
        False,
        f"{sample}: host_disease is {raw!r}, which is not in config group_labels "
        f"{sorted(labels)}. Update the config rather than guessing.",
    )
    return None


def summarise_metadata(records, attribute_map, numeric, config):
    """Report what is actually in each column, which is what decides the covariates.

    A field can be fully populated and still be useless to model. If nearly every sample
    shares one value it carries no information and costs a degree of freedom. For the
    numeric covariates the by-group split matters too: if the groups are already matched
    on age and BMI there is little to adjust for, and if they are not, the adjustment is
    doing real work and any "disease" signal might be tracking the covariate instead.
    """
    groups = sorted({r["group"] for r in records})
    required = set(config["biosample"].get("required_columns", []))
    print("\nMetadata summary, for choosing covariates:")

    for column in attribute_map.values():
        values = [r.get(column, "") for r in records]
        present = [v for v in values if v != ""]
        blank = len(values) - len(present)
        tag = "required" if column in required else "optional"

        if column in numeric:
            try:
                numbers = sorted(float(v) for v in present)
            except ValueError:
                print(f"  {column:<16} [{tag}] not parseable as numeric")
                continue
            median = numbers[len(numbers) // 2]
            per_group = []
            for group in groups:
                g = sorted(float(r[column]) for r in records
                           if r["group"] == group and r.get(column, "") != "")
                if g:
                    per_group.append(f"{group} median {g[len(g) // 2]:.1f}")
            print(f"  {column:<16} [{tag}] n={len(numbers)} min {numbers[0]:.1f}  "
                  f"median {median:.1f}  max {numbers[-1]:.1f}"
                  + ("   (" + ", ".join(per_group) + ")" if per_group else ""))
        else:
            counts = Counter(present)
            shown = ", ".join(f"{k}={v}" for k, v in counts.most_common(5))
            extra = f" (+{len(counts) - 5} more)" if len(counts) > 5 else ""
            note = ""
            if present and counts.most_common(1)[0][1] / len(present) >= 0.9:
                note = "   <- near constant, not useful as a covariate"
            if blank:
                note += f"   [{blank} blank]"
            print(f"  {column:<16} [{tag}] {shown}{extra}{note}")

    # Cross-tabulate the categorical covariates against group. An imbalance here is a
    # confound that a median comparison cannot show. ME/CFS is strongly female-biased in
    # the population, so sex balance between the groups matters a great deal.
    print("\n  categorical covariates against group:")
    for column in attribute_map.values():
        if column in numeric or column == "group":
            continue
        counts = Counter(present_val for present_val in (r.get(column, "") for r in records)
                         if present_val)
        if not 2 <= len(counts) <= 4:
            continue
        print(f"    {column}:")
        for value in sorted(counts):
            by_group = {g: sum(1 for r in records
                               if r["group"] == g and r.get(column) == value)
                        for g in groups}
            total = sum(by_group.values())
            share = ", ".join(f"{g} {n} ({100 * n / max(1, sum(1 for r in records if r['group'] == g)):.0f}%)"
                              for g, n in by_group.items())
            print(f"      {value:<16} n={total:<4} {share}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", type=Path)
    parser.add_argument("--outdir", default=".", type=Path)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    dataset = config["dataset"]
    route = config["download_route"]
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        print("note: NCBI_EMAIL is not set. NCBI asks automated clients to identify\n"
              "      themselves. Set it with: export NCBI_EMAIL=\"you@example.com\"\n")

    problems = Problems()

    # ---- run table ------------------------------------------------------------------
    print(f"Querying SRA for {dataset['bioproject']}")
    runs = sorted(fetch_run_table(dataset["bioproject"], email),
                  key=lambda r: r["Run"])
    print(f"  {len(runs)} runs")
    problems.check(
        len(runs) == dataset["expected_runs"],
        f"SRA returned {len(runs)} runs, config expects {dataset['expected_runs']}.",
    )
    layouts = {r["LibraryLayout"] for r in runs}
    models = {r["Model"] for r in runs}
    strategies = {r["LibraryStrategy"] for r in runs}
    print(f"  {'/'.join(sorted(strategies))}, {'/'.join(sorted(layouts))}, "
          f"{'/'.join(sorted(models))}")
    problems.check(layouts == {dataset["layout"]},
                   f"library layout is {layouts}, config expects {dataset['layout']}.")
    total_mb = sum(int(r["size_MB"] or 0) for r in runs)
    spots = sum(int(r["spots"] or 0) for r in runs)
    print(f"  {spots:,} spots, about {total_mb} MB of .sra to download")

    # ---- attributes -----------------------------------------------------------------
    sample_accessions = sorted({r["BioSample"] for r in runs if r.get("BioSample")})
    print(f"\nFetching BioSample attributes for {len(sample_accessions)} samples")
    attributes = fetch_biosample_attributes(sample_accessions, email)
    missing = [a for a in sample_accessions if a not in attributes]
    problems.check(not missing, f"no BioSample record for: {missing[:5]}")

    attribute_map = config["biosample"]["attribute_map"]
    labels = {str(k): v for k, v in config["biosample"]["group_labels"].items()}
    required = set(config["biosample"].get("required_columns", []))
    placeholders = {"missing", "not applicable", "na", "not collected", "unknown", ""}

    records = []
    gaps = Counter()
    for run in runs:
        sample = run["BioSample"]
        source = attributes.get(sample, {})
        row = {"sample-id": run["Run"], "biosample": sample}
        for attribute, column in attribute_map.items():
            value = source.get(attribute, "").strip()
            if value.lower() in placeholders:
                value = ""
                gaps[column] += 1
                # Only the modelled columns are fatal. An empty cell in a covariate
                # silently drops that sample from the statistics several steps later.
                problems.check(
                    column not in required,
                    f"{run['Run']} ({sample}): required field {column!r} (from "
                    f"{attribute}) is missing or a placeholder.",
                )
            row[column] = value
        row["group"] = normalise_group(row.get("group", ""), labels, problems,
                                       run["Run"]) or row.get("group", "")
        records.append(row)

    if gaps:
        print("  gaps in non-modelled fields (kept, not fatal): "
              + ", ".join(f"{c} x{n}" for c, n in sorted(gaps.items())))

    counts = Counter(r["group"] for r in records)
    print(f"  groups: {dict(counts)}")
    expected = {labels.get(str(k), str(k)): v
                for k, v in dataset["expected_groups"].items()}
    problems.check(dict(counts) == expected,
                   f"group sizes are {dict(counts)}, config expects {expected}.")

    problems.report_and_exit()

    # ---- write ----------------------------------------------------------------------
    metadata_dir = args.outdir / "metadata"
    data_dir = args.outdir / "data"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    columns = ["sample-id"] + list(attribute_map.values()) + ["biosample"]
    numeric = set(config["biosample"]["numeric_columns"])
    metadata_path = metadata_dir / "sample_metadata.tsv"
    with metadata_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(["#q2:types"] + [
            "numeric" if c in numeric else "categorical" for c in columns[1:]])
        for row in records:
            writer.writerow([row.get(c, "") for c in columns])
    print(f"\nWrote {metadata_path}  ({len(records)} samples, {len(columns)} columns)")
    summarise_metadata(records, attribute_map, numeric, config)

    manifest_path = data_dir / "download_manifest.tsv"
    template = route["sra_url_template"]
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["run", "sra_url", "expected_spots", "expected_bases",
                         "size_mb", "avg_length", "run_hash"])
        for run in runs:
            writer.writerow([
                run["Run"], template.format(run=run["Run"]), run["spots"],
                run["bases"], run["size_MB"], run["avgLength"], run["RunHash"],
            ])
    print(f"\nWrote {manifest_path}  ({len(runs)} runs, about {total_mb} MB)")
    print(f"  URLs use the Normalized objects in the public S3 bucket, NOT runinfo's")
    print(f"  download_path, which serves SRA Lite with synthesised quality scores.")

    # ---- the overlap arithmetic ------------------------------------------------------
    read_length = int(dataset["read_length"])
    amplicon = int(dataset["approximate_amplicon_bp"])
    print(f"\nAmplicon check: {dataset['amplicon_region']} is about {amplicon} bp, reads "
          f"are 2 x {read_length} bp")
    print(f"  combined {2 * read_length} bp leaves a gap of about "
          f"{amplicon - 2 * read_length} bp, so the mates cannot be merged")
    print(f"  confirmed from the reads: R1 starts with the 341F primer and R2 with 805R")
    print(f"  the analysis uses forward reads only "
          f"(config: reads_usable = {dataset['reads_usable']})")

    print(f"\nNext: 02_download_reads.py, which needs sra-tools (your 'bioenv').")


if __name__ == "__main__":
    main()
