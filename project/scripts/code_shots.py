#!/usr/bin/env python3
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""Render the presentation's code screenshots to results/m2/code_shots/*.png.

Each shot is a slide-ready PNG of real source lines, syntax highlighted,
with the lines that matter marked and a label naming the file and lines,
so a frame from the video can always be traced back to the repo.

Line ranges are found by searching for text in the source, not by line
number, so a shot follows its code when lines move.

The Milestone 1 shots in results/code_shots/ came from the code at git tag
milestone1 and are left as they are.

Usage: .venv/bin/python scripts/code_shots.py [--out results/m2/code_shots]
Needs pygments and pillow (scripts/requirements.txt, presentation only).
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
LABEL_COLOUR = "#6B7280"   # DARK_GREY, matches the figures
CAPTION_COLOUR = "#374151"  # TEXT, matches the figures
HL_COLOUR = "#FFF3B8"       # soft yellow band behind the lines that matter
PAD = 24                    # outer margin and gap between panes, px


def lexer_for(path: str):
    if path.endswith((".cpp", ".h")):
        return CppLexer()
    if path.endswith(".sh"):
        return BashLexer()
    return TextLexer()


def find_line(path: pathlib.Path, text: str, after: int = 0) -> int:
    """1-based number of the first line after `after` that contains text.
    Stops the script if the text is gone, so a stale shot is never drawn."""
    lines = path.read_text().splitlines()
    for number in range(after + 1, len(lines) + 1):
        if text in lines[number - 1]:
            return number
    sys.exit(f"code_shots: '{text}' not found in {path.relative_to(ROOT)}")


def span(path: pathlib.Path, first_text: str, last_text: str):
    """(first, last) line numbers from the line holding first_text to the
    next line holding last_text."""
    first = find_line(path, first_text)
    last = find_line(path, last_text, first - 1)
    return (first, last)


def lines_matching(path: pathlib.Path, ranges, texts):
    """Line numbers inside ranges whose line contains any of texts."""
    lines = path.read_text().splitlines()
    found = []
    for lo, hi in ranges:
        for number in range(lo, hi + 1):
            if any(text in lines[number - 1] for text in texts):
                found.append(number)
    return found


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

    # Pane label strip on top, which file and which lines. The pane widens if
    # the label is longer than the code, so the label is never cut off.
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


def code_pane(path, first_text, last_text, highlight_texts, label):
    """A pane from first_text to last_text with highlight_texts marked."""
    ranges = [span(path, first_text, last_text)]
    marked = lines_matching(path, ranges, highlight_texts)
    return render_pane(path, ranges, marked, label)


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
    parser.add_argument("--out", default=str(ROOT / "results" / "m2" / "code_shots"))
    args = parser.parse_args()
    out_dir = pathlib.Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    stencil = ROOT / "src" / "stencil.cpp"
    stencil_avx = ROOT / "src" / "stencil_avx.cpp"
    solver_simd = ROOT / "src" / "solver_simd.cpp"
    solver_omp = ROOT / "src" / "solver_omp.cpp"
    solver_mpi = ROOT / "src" / "solver_mpi.cpp"
    common = ROOT / "slurm" / "bench_common.sh"

    # The data structure. Two flat sheets, the index mapping and the loop nest that walks them.
    grid_h = ROOT / "include" / "grid.h"
    baseline = ROOT / "src" / "solver_baseline.cpp"
    compose(
        [code_pane(grid_h,
                   "Position of node (stock_i, var_j) in the flat buffer.",
                   "std::vector<double>& next() { return next_; }",
                   ["int index(", "current() const", "next() {"],
                   "the grid, two flat sheets"),
         code_pane(baseline,
                   "for (int step = 1; step <= num_timesteps; ++step) {",
                   "double south_west = current[grid.index(stock_i - 1, var_j - 1)];",
                   ["for (int", "double north = ", "double south = "],
                   "the loop nest")],
        "Each sheet is ns x nv doubles stored row after row. East and west are "
        "next door in memory, north and south are one row, ns doubles, away. "
        "At 2048 x 512 a row is 16 KiB, a sheet 8 MiB, and the two sheets 16 MiB.",
        out_dir / "shot_data_layout.png",
    )

    # Rung 1, serial. The scalar row update every later rung is checked against.
    compose(
        [code_pane(stencil,
                   "Central differences give each derivative",
                   "out_row[stock_i] = centre + dt * change;",
                   ["double change", "out_row[stock_i]"],
                   "serial, one cell per pass")],
        "One cell per pass. Five derivatives, five PDE terms, one step back "
        "in time. Measured 159 Mcell-updates/s, 2.8x the baseline.",
        out_dir / "shot_serial_row.png",
    )

    # Rung 2, simd. The scalar sum on the left, the same sum with FMA on the right.
    compose(
        [code_pane(stencil,
                   "The five terms of the Heston PDE",
                   "out_row[stock_i] = centre + dt * change;",
                   ["_term =", "double change", "out_row[stock_i]"],
                   "serial, one cell"),
         code_pane(stencil_avx,
                   "Add up the five PDE terms and the discounting.",
                   "_mm256_storeu_pd(out_row + stock_i, new_value);",
                   ["change = _mm256", "_mm256_fmadd_pd(dt", "_mm256_storeu_pd"],
                   "simd, four cells")],
        "The same sum on four neighbouring cells at once, each multiply and "
        "add fused into one instruction. Measured 2.4x over serial.",
        out_dir / "shot_simd_fma.png",
    )

    # Rung 3, omp. The serial row loop on the left, the same loop with one pragma on the right.
    compose(
        [code_pane(solver_simd,
                   "Update every interior row, four cells at a time.",
                   "}",
                   [],
                   "simd, one thread"),
         code_pane(solver_omp,
                   "Split the interior rows between the threads",
                   "}",
                   ["#pragma omp parallel for"],
                   "omp, every core")],
        "One line splits the rows between 24 threads. Each thread writes only "
        "its own rows, so no locks. Measured 15.1x over simd.",
        out_dir / "shot_omp_pragma.png",
    )

    # Rung 4, mpi. The edge row exchange, receive posted before send.
    compose(
        [code_pane(solver_mpi,
                   "The bottom process has nobody below it",
                   "MPI_Waitall(",
                   ["MPI_Irecv", "MPI_Isend", "MPI_Waitall"],
                   "mpi, swapping edge rows each step")],
        "Each process owns a block of rows and swaps one edge row with each "
        "neighbour every step. On one node 2x12 beat 1x24 by 1.05x.",
        out_dir / "shot_mpi_exchange.png",
    )

    # No benchmark without a just-passed test suite.
    compose(
        [code_pane(common,
                   "build_and_validate() {",
                   "make avx-check",
                   ["make test", "make avx-check"],
                   "runs before every sweep")],
        "Every timed binary has just passed the tests, or the job stops.",
        out_dir / "shot_iron_rule.png",
    )

    # The header every benchmark CSV carries.
    csvs = sorted((ROOT / "results").glob("bench_m2_ladder_*.csv"))
    if csvs:
        header = "\n".join(
            line for line in csvs[-1].read_text().splitlines()
            if line.startswith("#") and not line.startswith("# sweep,"))
        compose(
            [render_text_pane(header, f"results/{csvs[-1].name}, header")],
            "Written by the job itself. Host, job id, compiler, flags, "
            "binding, CPU and caches. The plot script refuses CSVs from "
            "anywhere but a compute node.",
            out_dir / "shot_csv_header.png",
        )
    else:
        print("skip shot_csv_header: no results/bench_m2_ladder_*.csv on disk",
              file=sys.stderr)

    # The matching test's own output from a compute node, where simd really is AVX2.
    # Only the call rows of the 24 thread run, the put rows say the same
    log = ROOT / "results" / "logs" / "slurm-651295.out"
    if log.exists():
        lines = log.read_text().splitlines()
        start = lines.index("== tests/test_ladder_matches")
        kept = []
        for line in lines[start + 1:]:
            if " call: " in line or line == "test_ladder_matches: PASS":
                kept.append(line)
            if line == "test_ladder_matches: PASS":
                break
        compose(
            [render_text_pane("\n".join(kept), "tests/test_ladder_matches on r730, 24 threads, job 651295")],
            "Every rung matches the baseline, and the parallel rungs match "
            "simd bit for bit.",
            out_dir / "shot_ladder_test.png",
        )
    else:
        print(f"skip shot_ladder_test: {log.relative_to(ROOT)} not on disk",
              file=sys.stderr)

if __name__ == "__main__":
    main()
