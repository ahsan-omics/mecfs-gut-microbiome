#!/usr/bin/env python3
"""Draw the workflow diagram, light and dark, as PNGs for the README.

Needs a Chromium or Chrome binary to rasterise; it is not part of the analysis and is not
in environment.yml. Run it only when the diagram changes.

The diagram deliberately does not show nine numbered boxes in a row. A reader can get
that from the step table in the README. What it shows instead are the three things that
are hard to see from prose:

  two archive sources join      reads from SRA, phenotypes from BioSample
  the reverse reads are dropped with the measured reason on the arrow
  the table splits two ways     rarefied for diversity, unrarefied for ANCOM-BC

That last split is the one design decision in the pipeline most likely to be got wrong
elsewhere, so it is the element carrying the colour.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

THEME = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e",
              "muted": "#898781", "rule": "#c3c2b7", "panel": "#f2f1ec",
              "blue": "#2a78d6", "orange": "#eb6834"},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7",
             "muted": "#898781", "rule": "#383835", "panel": "#242422",
             "blue": "#3987e5", "orange": "#d95926"},
}

W, H = 880, 846

# x, y, width, height, title, subtitle, emphasis
NODES = [
    (50, 24, 280, 52, "NCBI SRA", "63 runs, paired 2 x 151 bp", "panel"),
    (550, 24, 280, 52, "NCBI BioSample", "disease, age, sex, BMI", "panel"),
    (255, 128, 370, 56, "01-02  fetch and verify", "spot counts checked against SRA", "node"),
    (255, 232, 370, 56, "03  read diagnostics", "quality, primers, overlap test", "node"),
    (255, 336, 370, 56, "04  cutadapt, then DADA2", "3,601 ASVs at 134 bp", "node"),
    (255, 440, 370, 56, "05  classify", "vsearch consensus", "node"),
    (668, 440, 180, 56, "SILVA 138", "trimmed to V3", "panel"),
    (255, 544, 370, 56, "06  filter and phylogeny", "2,908 ASVs, rooted tree", "node"),
    (60, 660, 320, 56, "07  diversity", "alpha, beta, PCoA", "node"),
    (500, 660, 320, 56, "08  statistics", "PERMANOVA, PERMDISP, ANCOM-BC", "node"),
    (255, 764, 370, 56, "09  figures and tables", "", "node"),
]


def svg(theme):
    t = THEME[theme]
    out = [
        f'<svg viewBox="0 0 {W} {H}" width="{W}" height="{H}" role="img" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'aria-label="Workflow: reads from SRA and phenotypes from BioSample join, the '
        f'reverse reads are dropped because the mates do not overlap, and the feature '
        f'table splits into a rarefied branch for diversity and an unrarefied branch for '
        f'differential abundance.">',
        '<defs>',
        f'<marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
        f'markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{t["ink2"]}"/></marker>',
        f'<marker id="ab" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
        f'markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{t["blue"]}"/></marker>',
        f'<marker id="ao" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
        f'markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{t["orange"]}"/></marker>',
        '</defs>',
        f'<rect width="{W}" height="{H}" fill="{t["surface"]}"/>',
    ]

    for x, y, w, h, title, subtitle, kind in NODES:
        fill = t["panel"] if kind == "panel" else t["surface"]
        dash = ' stroke-dasharray="5 4"' if kind == "panel" else ""
        out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{fill}" '
                   f'stroke="{t["rule"]}" stroke-width="1.5"{dash}/>')
        cx = x + w / 2
        ty = y + (33 if subtitle else 34)
        out.append(f'<text x="{cx}" y="{ty}" text-anchor="middle" font-size="14" '
                   f'font-weight="600" fill="{t["ink"]}">{title}</text>')
        if subtitle:
            out.append(f'<text x="{cx}" y="{y + 48}" text-anchor="middle" font-size="11.5" '
                       f'fill="{t["muted"]}">{subtitle}</text>')

    def line(points, colour, marker, width=1.8, dash=""):
        pts = " ".join(f"{a},{b}" for a, b in points)
        d = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(f'<polyline points="{pts}" fill="none" stroke="{colour}" '
                   f'stroke-width="{width}" marker-end="url(#{marker})"{d}/>')

    def label(x, y, text, colour, anchor="middle", weight="400", size=11.5):
        out.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
                   f'font-weight="{weight}" fill="{colour}">{text}</text>')

    ink2, muted, blue, orange = t["ink2"], t["muted"], t["blue"], t["orange"]

    # Two archive sources converge.
    # Labels sit centred on the horizontal run, short enough to clear the two vertical
    # drops on either side rather than crossing their own arrows.
    line([(190, 76), (190, 104), (340, 104), (340, 128)], ink2, "a")
    label(265, 97, "full quality scores", muted)
    line([(690, 76), (690, 104), (540, 104), (540, 128)], ink2, "a")
    label(615, 97, "per-sample phenotypes", muted)

    line([(440, 184), (440, 232)], ink2, "a")
    line([(440, 288), (440, 336)], ink2, "a")
    label(452, 316, "forward reads only: mates do not overlap", muted, anchor="start")

    line([(440, 392), (440, 440)], ink2, "a")
    line([(668, 468), (629, 468)], ink2, "a")
    line([(440, 496), (440, 544)], ink2, "a")

    # The split. This is the claim the figure exists to make, so it carries the colour.
    line([(440, 600), (440, 630), (220, 630), (220, 660)], blue, "ab", width=2.2)
    label(228, 623, "rarefied to 18,848", blue, anchor="start", weight="600")
    line([(440, 600), (440, 630), (660, 630), (660, 660)], orange, "ao", width=2.2)
    label(652, 623, "unrarefied counts", orange, anchor="end", weight="600")

    line([(220, 716), (220, 742), (436, 742), (436, 764)], blue, "ab")
    line([(660, 716), (660, 742), (444, 742), (444, 764)], orange, "ao")

    out.append('</svg>')
    return "\n".join(out)


def html(theme):
    t = THEME[theme]
    # overflow:hidden matters: without it the window is a few pixels short of the drawing,
    # Chromium shows scrollbars, and the screenshot captures them as grey bars.
    return (f'<!doctype html><meta charset="utf-8">'
            f'<style>html,body{{margin:0;padding:0;overflow:hidden;'
            f'background:{t["surface"]};'
            f'font-family:system-ui,-apple-system,"Segoe UI",sans-serif}}'
            f'svg{{display:block}}</style>'
            f'{svg(theme)}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="results/figures", type=Path)
    parser.add_argument("--scale", type=float, default=2.0)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    browser = next((shutil.which(b) for b in
                    ("chromium", "chromium-browser", "google-chrome", "chrome")
                    if shutil.which(b)), None)
    if not browser:
        sys.exit("No Chromium or Chrome binary found; cannot rasterise the diagram.")

    for theme in ("light", "dark"):
        # Written beside the output rather than in /tmp: a snap-packaged Chromium cannot
        # read /tmp, and silently produces a blank image instead of an error.
        page = args.outdir / f".workflow_{theme}.html"
        page.write_text(html(theme))
        png = args.outdir / f"workflow_{theme}.png"
        subprocess.run([browser, "--headless", "--disable-gpu", "--no-sandbox",
                        f"--screenshot={png.resolve()}",
                        f"--window-size={W},{H}",
                        f"--force-device-scale-factor={args.scale}",
                        f"--virtual-time-budget=2000",
                        page.resolve().as_uri()],
                       capture_output=True, text=True, timeout=180)
        page.unlink(missing_ok=True)
        if png.exists():
            print(f"  wrote {png}  ({png.stat().st_size / 1024:.0f} KB)")
        else:
            print(f"  FAILED to write {png}", file=sys.stderr)


if __name__ == "__main__":
    main()
