#!/usr/bin/env python3
"""Render the presentation's code screenshots to results/code_shots/*.png.

Each shot is a slide-ready PNG of real source lines, syntax-highlighted,
with the load-bearing lines marked and a footer naming exactly where the
code lives (file:lines) — the same provenance habit as the benchmark
figures, so a frame from the video can always be traced back to the repo.

Side-by-side shots put the "before" pane on the left and the "after" (or
the deliberately wrong control) on the right, because the diff between
neighbours IS the technique being shown (PLAN §4b).

Usage: .venv/bin/python scripts/code_shots.py [--out results/code_shots]
Needs pygments + pillow (scripts/requirements.txt, presentation only).
"""

import argparse
import io
import pathlib
import subprocess
import sys

from PIL import Image, ImageDraw, ImageFont
from pygments import highlight
from pygments.formatters import ImageFormatter
from pygments.lexers import CppLexer, BashLexer, TextLexer

ROOT = pathlib.Path(__file__).resolve().parent.parent

FONT_PATH = "/System/Library/Fonts/Menlo.ttc"
CODE_FONT_SIZE = 22
LABEL_FONT_SIZE = 20
CAPTION_FONT_SIZE = 24
BG = "#FFFFFF"
LABEL_COLOUR = "#6B7280"   # DARK_GREY, matches bench_plot.py annotations
CAPTION_COLOUR = "#374151"  # TEXT, matches bench_plot.py
HL_COLOUR = "#FFF3B8"       # soft yellow band behind the lines that matter
PAD = 24                    # outer margin and inter-pane gap, px


def lexer_for(path: str):
    if path.endswith((".cpp", ".h")):
        return CppLexer()
    if path.endswith(".sh"):
        return BashLexer()
    return TextLexer()


def slice_ranges(path: pathlib.Path, ranges):
    """Return (text, mapping) where text joins the requested 1-based line
    ranges with an elision marker, and mapping translates original line
    numbers to line numbers inside the joined text (for highlighting)."""
    lines = path.read_text().splitlines()
    out, mapping = [], {}
    for which, (lo, hi) in enumerate(ranges):
        if which:
            out.append("        // ...")  # elision between joined ranges
        for n in range(lo, hi + 1):
            mapping[n] = len(out) + 1
            out.append(lines[n - 1])
    return "\n".join(out), mapping


def render_pane(path: pathlib.Path, ranges, hl_original_lines, label):
    text, mapping = slice_ranges(path, ranges)
    hl = [mapping[n] for n in hl_original_lines if n in mapping]
    png = highlight(
        text,
        lexer_for(str(path)),
        ImageFormatter(
            font_name="Menlo",
            font_size=CODE_FONT_SIZE,
            style="xcode",
            line_numbers=False,
            hl_lines=hl,
            hl_color=HL_COLOUR,
            image_pad=14,
        ),
    )
    img = Image.open(io.BytesIO(png)).convert("RGB")

    # Pane label strip on top: which file, which lines. The pane widens if
    # the label is longer than the code, so provenance is never cut off.
    font = ImageFont.truetype(FONT_PATH, LABEL_FONT_SIZE)
    strip_h = LABEL_FONT_SIZE + 14
    where = ", ".join(f"{lo}-{hi}" for lo, hi in ranges)
    text = f"{label}   {path.relative_to(ROOT)}:{where}"
    text_w = int(ImageDraw.Draw(img).textlength(text, font=font)) + 8
    labelled = Image.new("RGB", (max(img.width, text_w),
                                 img.height + strip_h), BG)
    labelled.paste(img, (0, strip_h))
    ImageDraw.Draw(labelled).text((4, 4), text, fill=LABEL_COLOUR, font=font)
    return labelled


def render_text_pane(text: str, label: str):
    png = highlight(
        text, TextLexer(),
        ImageFormatter(font_name="Menlo", font_size=CODE_FONT_SIZE,
                       style="xcode", line_numbers=False, image_pad=14),
    )
    img = Image.open(io.BytesIO(png)).convert("RGB")
    font = ImageFont.truetype(FONT_PATH, LABEL_FONT_SIZE)
    strip_h = LABEL_FONT_SIZE + 14
    labelled = Image.new("RGB", (img.width, img.height + strip_h), BG)
    labelled.paste(img, (0, strip_h))
    ImageDraw.Draw(labelled).text((4, 4), label, fill=LABEL_COLOUR, font=font)
    return labelled


def wrap_caption(draw, caption, font, max_w):
    """Greedy word wrap so a caption longer than the panes folds onto
    extra lines instead of running off both edges of the canvas."""
    lines, current = [], ""
    for word in caption.split():
        trial = f"{current} {word}".strip()
        if not current or draw.textlength(trial, font=font) <= max_w:
            current = trial
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def compose(panes, caption, out_path: pathlib.Path):
    """Panes side by side, one caption centred underneath."""
    font = ImageFont.truetype(FONT_PATH, CAPTION_FONT_SIZE)
    line_h = CAPTION_FONT_SIZE + 8
    width = sum(p.width for p in panes) + PAD * (len(panes) + 1)
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    lines = wrap_caption(measure, caption, font, width - 2 * PAD)
    cap_h = line_h * len(lines) + 2 * PAD
    height = max(p.height for p in panes) + 2 * PAD + cap_h
    canvas = Image.new("RGB", (width, height), BG)
    x = PAD
    for p in panes:
        canvas.paste(p, (x, PAD))
        x += p.width + PAD
    draw = ImageDraw.Draw(canvas)
    for row, line in enumerate(lines):
        text_w = draw.textlength(line, font=font)
        draw.text(((width - text_w) / 2,
                   height - cap_h + PAD + row * line_h),
                  line, fill=CAPTION_COLOUR, font=font)
    canvas.save(out_path)
    print(f"wrote {out_path}  ({canvas.width}x{canvas.height})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(ROOT / "results" / "code_shots"))
    args = parser.parse_args()
    out_dir = pathlib.Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    kernels = ROOT / "src" / "solver_opt_kernels.cpp"
    baseline = ROOT / "src" / "solver_baseline.cpp"
    common = ROOT / "slurm" / "bench_common.sh"

    # 0. One shot per ladder rung, before on the left (the level below) and
    #    after on the right, highlighting exactly the lines that rung changes.
    #    Measured effects quote Table 2 of RESULTS.md (sweep E, huge pages).
    compose(
        [
            render_pane(baseline, [(80, 89)], [],
                        "baseline — the reference stencil"),
            render_pane(kernels, [(54, 67)], [],
                        "level 0 — the same arithmetic, line for line"),
        ],
        "No techniques applied. Measured: -0.1 % (Table 2), so the ladder "
        "harness itself costs nothing and changes no answers.",
        out_dir / "shot_L0_anchor.png",
    )

    compose(
        [
            render_pane(kernels, [(41, 61)],
                        [43, 45, 54, 55, 56, 57, 58, 59, 60, 61],
                        "level 0 — everything recomputed per cell"),
            render_pane(kernels, [(114, 121), (126, 129), (135, 139)],
                        [114, 115, 116, 117, 118, 120, 121,
                         126, 135, 136, 137, 138, 139],
                        "level 1 — hoisted per step and per row"),
        ],
        "Hoisting, a lookup table and CSE move work to the level where it "
        "varies. Names only, never reordering, so the answer stays "
        "bit-identical. Measured: +0.3 % (Table 2).",
        out_dir / "shot_L1_hoisting.png",
    )

    compose(
        [
            render_pane(kernels, [(227, 243)],
                        [235, 236, 237, 238, 239, 240, 241, 242, 243],
                        "level 2 — g.index arithmetic on every access"),
            render_pane(kernels, [(316, 337)],
                        [318, 319, 320, 329, 330, 331, 332, 333,
                         334, 335, 336, 337],
                        "level 3 — row bases, memory walked in layout order"),
        ],
        "Variance outside, stock inside, three row bases per row. Measured: "
        "+0.3 % (Table 2). The baseline already walked this way, and the "
        "ctl-order control prices the wrong way.",
        out_dir / "shot_L3_order.png",
    )

    compose(
        [
            render_pane(kernels, [(316, 328)], [316, 328],
                        "level 3 — bounds tested inside the loop condition"),
            render_pane(kernels, [(407, 416), (447, 456)],
                        [409, 410, 447, 448, 450],
                        "level 4 — plain locals, one loop per equation"),
        ],
        "Loop splitting and a branch-free interior. No cell tests which "
        "equation it follows. Measured: -0.2 % (Table 2), and the "
        "ctl-branch control times the fused version.",
        out_dir / "shot_L4_split.png",
    )

    compose(
        [
            render_pane(kernels, [(413, 433)],
                        [414, 415, 416, 425, 426, 427, 428, 429,
                         430, 431, 432, 433],
                        "level 4 — row + stock_i on every access"),
            render_pane(kernels, [(507, 530)],
                        [510, 511, 512, 513, 522, 523, 524, 525,
                         526, 527, 528, 529, 530],
                        "level 5 — raw row pointers, addressing does the work"),
        ],
        "Induction variable simplification. Ten index additions per cell "
        "become pointer indexing. Measured: +2.2 % (Table 2).",
        out_dir / "shot_L5_pointers.png",
    )

    compose(
        [
            render_pane(kernels, [(521, 541)], [521],
                        "level 5 — one cell per iteration"),
            render_pane(kernels, [(617, 627), (707, 710)],
                        [622, 623, 624, 625, 626, 709],
                        "level 6 — four independent cells per iteration"),
        ],
        "Unroll by four. Four dependency chains between one pair of "
        "loop-back branches. Measured: +0.5 % (Table 2). At -O2 the "
        "compiler already unrolls counted loops.",
        out_dir / "shot_L6_unroll.png",
    )

    # 1. The punchline shot: the one transformation the compiler is forbidden
    #    to make at -O2. Left, level 1's five divisions; right, level 2's
    #    per-step reciprocals and the multiply-only stencil.
    compose(
        [
            render_pane(kernels, [(154, 162)], [154, 155, 156, 157, 158],
                        "level 1 — five divisions per cell"),
            render_pane(kernels, [(206, 214), (244, 253)],
                        [208, 209, 210, 211, 212, 213,
                         244, 245, 246, 247, 248],
                        "level 2 — invert once per step, multiply per cell"),
        ],
        "Strength reduction: x/(2·ds) → x·inv_two_ds changes the last bits, "
        "so -O2 must not do it.  Measured: +148.7 % throughput (Table 2).",
        out_dir / "shot_L2_strength_reduction.png",
    )

    # 2. The loop-order negative control: identical arithmetic, loops swapped.
    compose(
        [
            render_pane(kernels, [(507, 521)], [507, 521],
                        "level 5 — variance outside, stock inside"),
            render_pane(kernels, [(806, 822)], [809, 810],
                        "ctl-order — swapped on purpose"),
        ],
        "Same arithmetic, same tables, loops swapped: 2.1× slower (Table 3). "
        "This is what the null level-3 rung was worth.",
        out_dir / "shot_ctl_order.png",
    )

    # 3. The iron rule: no benchmark without a just-passed test suite.
    compose(
        [render_pane(common, [(23, 32)], [28, 29, 30],
                     "build_and_validate — runs before every sweep")],
        "No number in RESULTS.md came from an unvalidated binary.",
        out_dir / "shot_iron_rule.png",
    )

    # 4. Provenance: the header every benchmark CSV carries.
    csvs = sorted((ROOT / "results").glob("bench_pages_*.csv"))
    if csvs:
        header = "\n".join(
            line for line in csvs[-1].read_text().splitlines()
            if line.startswith("#"))
        compose(
            [render_text_pane(header, f"results/{csvs[-1].name}, header")],
            "Stamped by the job itself: host, job id, compiler, flags, CPU, "
            "caches. The plot scripts refuse CSVs without it.",
            out_dir / "shot_csv_header.png",
        )
    else:
        print("skip shot_csv_header: no results/bench_pages_*.csv on disk",
              file=sys.stderr)

    # 5. The matching gate's own words: run the test binary, show its output.
    test_bin = ROOT / "tests" / "test_opt_matches"
    if test_bin.exists():
        out = subprocess.run([str(test_bin)], cwd=ROOT, capture_output=True,
                             text=True).stdout
        compose(
            [render_text_pane(out.strip(), "./test_opt_matches (local run)")],
            "Every rung and both controls, held to the baseline's answer.",
            out_dir / "shot_optmatches.png",
        )
    else:
        print("skip shot_optmatches: build test_opt_matches first (make)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
