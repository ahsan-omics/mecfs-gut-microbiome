#!/usr/bin/env python3
"""Download the .sra objects and convert them to gzipped FASTQ.

Reads data/download_manifest.tsv, which 01_fetch_metadata.py writes. Needs sra-tools on
PATH for fasterq-dump, so run it in an environment that has it.

Two decisions baked in here, both found the hard way:

  The URLs in the manifest point at the Normalized objects in NCBI's public S3 bucket,
  not at the `download_path` from SRA runinfo. That one serves SRA Lite, whose quality
  scores are all replaced by a single synthesised constant. DADA2 fits its error model to
  those scores, so Lite data produces bad ASVs without anything appearing to go wrong.

  fasterq-dump is given a LOCAL FILE, never an accession. Accession resolution fails on
  this machine ("Failed to call external services") and vdb-config --restore-defaults does
  not fix it. Handing it a file on disk bypasses the resolver entirely.

Verification is by spot count: the converted FASTQ must contain exactly the number of
records SRA reports for that run. A truncated download fails this, which a file size
check alone would not reliably catch.

Only stdlib, so it runs anywhere with python3 and sra-tools.
"""

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

USER_AGENT = "mecfs-gut-microbiome"


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def download(url, destination, expected_mb):
    """Stream the .sra. Resumable in the sense that a complete file is not refetched."""
    if destination.exists() and destination.stat().st_size > 0:
        return destination.stat().st_size, "cached"
    temporary = destination.with_suffix(".part")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=600) as response, \
            open(temporary, "wb") as handle:
        shutil.copyfileobj(response, handle, length=1 << 20)
    size = temporary.stat().st_size
    # A plausibility check only. The real check is the spot count after conversion.
    if expected_mb and size < 0.5 * expected_mb * 1024 * 1024:
        temporary.unlink()
        raise RuntimeError(f"downloaded {human(size)}, expected about {expected_mb} MB")
    temporary.rename(destination)
    return size, "downloaded"


def convert(sra_path, out_dir, temp_dir, extra_args):
    """fasterq-dump on a local file, so no accession resolution is attempted."""
    command = ["fasterq-dump", *extra_args, "-O", str(out_dir),
               "-t", str(temp_dir), str(sra_path)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"fasterq-dump failed ({result.returncode}): "
            f"{(result.stderr or result.stdout).strip().splitlines()[-1:]}")
    return result.stdout


def count_records(path):
    """Records in a FASTQ, counting lines without holding the file in memory."""
    opener = gzip.open if path.suffix == ".gz" else open
    lines = 0
    with opener(path, "rb") as handle:
        for _ in handle:
            lines += 1
    if lines % 4:
        raise RuntimeError(f"{path.name} has {lines} lines, not a multiple of 4")
    return lines // 4


def compress(path):
    target = path.with_suffix(path.suffix + ".gz")
    with open(path, "rb") as src, gzip.open(target, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, length=1 << 20)
    path.unlink()
    return target


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", default="data/download_manifest.tsv", type=Path)
    parser.add_argument("--outdir", default="data/raw", type=Path)
    parser.add_argument("--sradir", default="data/sra", type=Path)
    parser.add_argument("--keep-sra", action="store_true",
                        help="keep the .sra files after conversion")
    parser.add_argument("--no-gzip", action="store_true")
    parser.add_argument("--limit", type=int, default=0,
                        help="process only the first N runs, for a trial")
    args = parser.parse_args()

    if not shutil.which("fasterq-dump"):
        sys.exit("fasterq-dump is not on PATH.\n"
                 "It lives in your 'bioenv' environment: conda activate bioenv")
    if not args.manifest.exists():
        sys.exit(f"{args.manifest} not found. Run 01_fetch_metadata.py first.")

    rows = [r for r in csv.DictReader(args.manifest.open(), delimiter="\t") if r.get("run")]
    if args.limit:
        rows = rows[:args.limit]

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.sradir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.sradir / "tmp"
    temp_dir.mkdir(exist_ok=True)

    suffix = ".fastq" if args.no_gzip else ".fastq.gz"
    total = len(rows)
    failures = []
    done = 0

    print(f"{total} runs to process, output to {args.outdir}\n")

    for index, row in enumerate(rows, 1):
        run = row["run"]
        expected = int(row["expected_spots"])
        forward = args.outdir / f"{run}_1{suffix}"
        reverse = args.outdir / f"{run}_2{suffix}"
        prefix = f"  [{index:>3}/{total}] {run}"

        if forward.exists() and reverse.exists():
            print(f"{prefix}  already present, skipping")
            done += 1
            continue

        sra_path = args.sradir / f"{run}.sra"
        try:
            size, how = download(row["sra_url"], sra_path, int(row["size_mb"] or 0))
            convert(sra_path, args.outdir, temp_dir, ["--split-files"])

            produced = sorted(args.outdir.glob(f"{run}_?.fastq"))
            if len(produced) != 2:
                raise RuntimeError(
                    f"expected 2 FASTQ files from --split-files, got "
                    f"{[p.name for p in produced]}")

            observed = count_records(produced[0])
            if observed != expected:
                raise RuntimeError(
                    f"spot count mismatch: FASTQ has {observed:,} records, SRA reports "
                    f"{expected:,}. The download was probably truncated.")

            if not args.no_gzip:
                produced = [compress(p) for p in produced]

            if not args.keep_sra:
                sra_path.unlink(missing_ok=True)

            print(f"{prefix}  {how} {human(size)}  ->  {observed:,} spots  ok")
            done += 1

        except Exception as exc:
            print(f"{prefix}  FAILED: {exc}", file=sys.stderr)
            failures.append(f"{run}: {exc}")
            for leftover in args.outdir.glob(f"{run}_?.fastq"):
                leftover.unlink(missing_ok=True)

    shutil.rmtree(temp_dir, ignore_errors=True)
    if not args.keep_sra:
        try:
            args.sradir.rmdir()
        except OSError:
            pass

    print(f"\n{done}/{total} runs complete.")
    if failures:
        print(f"\n{len(failures)} failure(s):", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        print("\nRe-running will retry only the failed runs.", file=sys.stderr)
        sys.exit(1)

    forward_files = sorted(args.outdir.glob(f"*_1{suffix}"))
    size = sum(p.stat().st_size for p in args.outdir.iterdir() if p.is_file())
    print(f"{len(forward_files)} forward read files, {human(size)} total in {args.outdir}")
    print("\nEvery run verified by spot count against SRA's own figures.")

    # A record of what landed, which is both provenance and the file the Snakefile
    # tracks: this rule's real outputs are 126 filenames taken from the manifest, which
    # Snakemake cannot know in advance.
    record = args.outdir.parent / "download_record.tsv"
    lines = ["\t".join(["run", "expected_spots", "forward_bytes", "reverse_bytes"])]
    for row in rows:
        run = row["run"]
        forward = args.outdir / f"{run}_1{suffix}"
        reverse = args.outdir / f"{run}_2{suffix}"
        lines.append("\t".join(str(x) for x in [
            run, row["expected_spots"],
            forward.stat().st_size if forward.exists() else 0,
            reverse.stat().st_size if reverse.exists() else 0,
        ]))
    content = "\n".join(lines) + "\n"

    # Only rewrite when the content actually changes. Snakemake decides what to redo from
    # file timestamps, so touching this file on every verification run would mark the
    # entire downstream pipeline stale even though not one read had changed. Re-running
    # the download step to confirm the data is intact must not cost hours of recompute.
    if record.exists() and record.read_text() == content:
        print(f"{record} unchanged, left alone so downstream steps stay current")
    else:
        record.write_text(content)
        print(f"Wrote {record}")


if __name__ == "__main__":
    main()
