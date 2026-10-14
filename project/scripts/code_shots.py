#!/usr/bin/env python3
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""Render the presentation's code screenshots to results/m2/code_shots/*.png.

Each shot is a slide-ready PNG of real source lines, syntax highlighted,
with the lines that matter marked and a label naming the file and lines,
so a frame from the video can always be traced back to the repo.

The text on a slide is only as large as the slot divided by the longest
line, so each shot keeps to a few lines, removes the common indent, and
wraps a long line onto the next one. A wrapped preprocessor line ends in a
backslash, which is how C++ continues a directive, so every shot is still
valid code. Before and after panes are stacked, one above the other.

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
import sys
import textwrap

from PIL import Image, ImageDraw, ImageFont
from pygments import highlight
from pygments.formatters import ImageFormatter
from pygments.lexers import BashLexer, CLexer, CppLexer, GasLexer, TextLexer

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Menlo on the Mac, DejaVu Sans Mono on Linux and WSL
FONTS = [("/System/Library/Fonts/Menlo.ttc", "Menlo"),
         ("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", "DejaVu Sans Mono")]
FONT_PATH, FONT_NAME = next(((path, name) for path, name in FONTS if pathlib.Path(path).exists()), FONTS[-1])
CODE_FONT_SIZE = 26
LABEL_FONT_SIZE = 22
BG = "#FFFFFF"
LABEL_COLOUR = "#6B7280"   # DARK_GREY, matches the figures
HL_COLOUR = "#FFF3B8"       # soft yellow band behind the lines that matter
PAD = 18                    # outer margin and gap between panes, px
# Longest line on a full width code slide. Fewer characters means larger text on the slide
WRAP = 72


def lexer_for(path: str):
    if path.endswith((".cpp", ".h")):
        return CppLexer()
    if path.endswith(".c"):
        return CLexer()
    if path.endswith(".sh") or path.endswith("Makefile"):
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


def span(path: pathlib.Path, first_text: str, last_text: str, after: int = 0):
    """(first, last) line numbers from the line holding first_text to the
    next line holding last_text."""
    first = find_line(path, first_text, after)
    last = find_line(path, last_text, first - 1)
    return (first, last)


def wrap_line(line: str, width: int):
    """Breaks one long line at spaces. Continuation lines are indented four
    more spaces, a comment carries on as a comment, and a preprocessor line
    gets a trailing backslash."""
    if len(line) <= width:
        return [line]
    indent = line[:len(line) - len(line.lstrip())]
    directive = line.lstrip().startswith("#")
    # A wrapped comment carries on as a comment
    follow = "// " if line.lstrip().startswith("//") else "    "
    pieces = textwrap.wrap(line.strip(), width=width - len(indent) - (2 if directive else 0),
                           subsequent_indent=follow, break_long_words=False, break_on_hyphens=False)
    if directive:
        pieces = [piece + " \\" for piece in pieces[:-1]] + [pieces[-1]]
    return [indent + piece for piece in pieces]


def slice_ranges(path: pathlib.Path, ranges, marked, width, dropped=lambda text: False):
    """Joins the requested 1-based line ranges with an elision marker, removes the
    common indent and wraps long lines. Lines for which dropped is true are left out.
    Returns the text and the output lines to mark."""
    lines = path.read_text().splitlines()
    chosen = []
    for which, (lo, hi) in enumerate(ranges):
        if which:
            chosen.append((None, None))
        for number in range(lo, hi + 1):
            if not dropped(lines[number - 1]):
                chosen.append((number, lines[number - 1]))
    indents = [len(text) - len(text.lstrip()) for number, text in chosen if number and text.strip()]
    cut = min(indents) if indents else 0
    out, hl = [], []
    for number, text in chosen:
        if number is None:
            # The elision between joined ranges, at the indent of the line after it
            out.append("// ...")
            continue
        for piece in wrap_line(text[cut:], width):
            out.append(piece)
            if number in marked:
                hl.append(len(out))
    return "\n".join(out), hl


def label_strip(img: Image.Image, text: str):
    """Puts a grey label above a pane, widening the pane if the label is longer."""
    font = ImageFont.truetype(FONT_PATH, LABEL_FONT_SIZE)
    strip_h = LABEL_FONT_SIZE + 12
    text_w = int(ImageDraw.Draw(img).textlength(text, font=font)) + 10
    labelled = Image.new("RGB", (max(img.width, text_w), img.height + strip_h), BG)
    labelled.paste(img, (0, strip_h))
    ImageDraw.Draw(labelled).text((4, 3), text, fill=LABEL_COLOUR, font=font)
    return labelled


def render(text: str, lexer, hl):
    png = highlight(text, lexer, ImageFormatter(font_name=FONT_NAME, font_size=CODE_FONT_SIZE, style="xcode",
                                                line_numbers=False, hl_lines=hl, hl_color=HL_COLOUR, image_pad=10))
    return Image.open(io.BytesIO(png)).convert("RGB")


def code_pane(path, label, ranges, highlight_texts=(), width=WRAP, leave_out=(), bare=False, start_after=None):
    """A pane of one or more line ranges, each given as (first_text, last_text), with lines holding highlight_texts marked.

    leave_out drops every line holding one of its texts, and bare drops comment lines and blank
    lines, so a long loop fits on a slide. The label then has to say what was left out.
    start_after starts the search below the line holding that text, for a line that also
    appears earlier in the file."""
    def dropped(text):
        stripped = text.strip()
        if bare and (not stripped or stripped.startswith("//")):
            return True
        return any(t in text for t in leave_out)

    numbers = []
    after = find_line(path, start_after) if start_after else 0
    for first_text, last_text in ranges:
        lo, hi = span(path, first_text, last_text, after)
        numbers.append((lo, hi))
        after = hi
    lines = path.read_text().splitlines()
    marked = {n for lo, hi in numbers for n in range(lo, hi + 1) if any(t in lines[n - 1] for t in highlight_texts)}
    text, hl = slice_ranges(path, numbers, marked, width, dropped)
    where = ", ".join(f"{lo}-{hi}" if lo != hi else f"{lo}" for lo, hi in numbers)
    return label_strip(render(text, lexer_for(str(path)), hl), f"{label}   {path.relative_to(ROOT)}:{where}")


def text_pane(text: str, label: str, lexer=None):
    return label_strip(render(text, lexer or TextLexer(), []), label)


def compose(panes, out_path: pathlib.Path, stacked=True):
    """Panes one above the other (stacked) or side by side, on a white canvas."""
    if stacked:
        width = max(p.width for p in panes) + 2 * PAD
        height = sum(p.height for p in panes) + PAD * (len(panes) + 1)
    else:
        width = sum(p.width for p in panes) + PAD * (len(panes) + 1)
        height = max(p.height for p in panes) + 2 * PAD
    canvas = Image.new("RGB", (width, height), BG)
    x = y = PAD
    for p in panes:
        canvas.paste(p, (x, y))
        if stacked:
            y += p.height + PAD
        else:
            x += p.width + PAD
    canvas.save(out_path)
    print(f"wrote {out_path}  ({canvas.width}x{canvas.height})")


def disassembly_lines(path: pathlib.Path, function: str, count: int):
    """The first count load lines of one function in an objdump listing, without their addresses."""
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.endswith(f"<{function}>:"))
    kept = []
    for line in lines[start + 1:]:
        if line.endswith(">:"):
            break
        instruction = line.split(":", 1)[1].strip() if ":" in line else line.strip()
        if instruction.startswith(("vmovupd", "vinsertf128")):
            kept.append(" ".join(instruction.split()))
        if len(kept) == count:
            break
    return kept


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(ROOT / "results" / "m2" / "code_shots"))
    args = parser.parse_args()
    out_dir = pathlib.Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    grid_h = ROOT / "include" / "grid.h"
    grid_cpp = ROOT / "src" / "grid.cpp"
    baseline = ROOT / "src" / "solver_baseline.cpp"
    stencil = ROOT / "src" / "stencil.cpp"
    stencil_avx = ROOT / "src" / "stencil_avx.cpp"
    solver_simd = ROOT / "src" / "solver_simd.cpp"
    solver_omp = ROOT / "src" / "solver_omp.cpp"
    solver_mpi = ROOT / "src" / "solver_mpi.cpp"
    main_cpp = ROOT / "src" / "main.cpp"
    common = ROOT / "slurm" / "bench_common.sh"
    flush = ROOT / "slurm" / "flush_denormals.c"

    # The data structure. The index mapping and the two flat buffers.
    compose([code_pane(grid_h, "the grid", [
        ("// Position of node (stock_i, var_j) in the flat buffer.", "int index("),
        ("// The two buffers.", "std::vector<double> next_;")], ["int index(", "std::vector<double> current_;", "std::vector<double> next_;"])],
        out_dir / "shot_data_layout.png")

    # Rung 1, serial. The division in the baseline, the reciprocal worked out once, the multiplication in the serial row.
    compose([
        code_pane(baseline, "before, baseline", [("double dV_dS = (east - west) /", "double d2V_dS2 = (east")], ["/ ("]),
        code_pane(stencil, "once per solve", [("inv_two_ds = 1.0 /", "inv_ds_sq = 1.0 /")], ["= 1.0 /"]),
        code_pane(stencil, "after, serial", [("double dV_dS = (east - west) * inv_two_ds;", "double d2V_dS2 = (east")], ["* inv_"])],
        out_dir / "shot_serial_change.png")

    # Rung 2, simd. The serial loop and the AVX2 loop, line for line. Each scalar load, the update
    # and the store has a vector twin that works on four cells, and the counter moves by four.
    compose([
        code_pane(stencil, "before, serial, one cell per iteration, comments left out", [
            ("for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {", "double west = row_centre[stock_i - 1];"),
            ("out_row[stock_i] = centre + dt * change;", "}")], bare=True),
        code_pane(stencil_avx, "after, simd, four cells per iteration, comments left out", [
            ("while (stock_i + 3 <= last_interior) {", "__m256d west = _mm256_loadu_pd(row_centre + stock_i - 1);"),
            ("__m256d new_value = _mm256_fmadd_pd(dt, change, centre);", "    }")],
            ["_mm256_loadu_pd", "_mm256_fmadd_pd", "_mm256_storeu_pd", "stock_i + 4"], bare=True)],
        out_dir / "shot_simd_change.png")

    # The scalar sum and the FMA chain it became, for the appendix.
    compose([
        code_pane(stencil, "serial, one cell", [("// The five terms of the Heston PDE.", "out_row[stock_i] = centre + dt * change;")],
                  ["double change", "out_row[stock_i]"]),
        code_pane(stencil_avx, "simd, four cells", [("// Add up the five PDE terms and the discounting.", "_mm256_storeu_pd(out_row + stock_i, new_value);")],
                  ["change = _mm256", "_mm256_fmadd_pd(dt", "_mm256_storeu_pd"])],
        out_dir / "shot_simd_fma.png")

    # What the compiler made of one unaligned load, before and after the two flags
    e2 = ROOT / "results" / "e2_disasm_657421.txt"
    e2b = ROOT / "results" / "e2b_disasm_657498.txt"
    row = "update_row_avx(StencilConstants const&, double const*, double*, int, double const*, double const*)"
    if e2.exists() and e2b.exists():
        compose([text_pane("\n".join(disassembly_lines(e2, row, 4)), "g++ 8.5, -mavx2 -mfma, each load split into two 128 bit loads"),
                 text_pane("\n".join(disassembly_lines(e2b, row, 2)), "with -mno-avx256-split-unaligned-load, one load")],
                out_dir / "shot_split_loads.png")
    else:
        print("skip shot_split_loads: the E2 and E2b disassembly files are not on disk", file=sys.stderr)

    # Rung 3, omp. The simd row loop, and the same loop with one pragma.
    compose([
        code_pane(solver_simd, "before, simd", [("// Update every interior row, four cells at a time.", "}")]),
        code_pane(solver_omp, "after, omp", [("// Split the interior rows between the threads", "}")], ["#pragma omp parallel for"])],
        out_dir / "shot_omp_pragma.png")

    # Rung 4, mpi. The omp time loop with the main changes for MPI marked. Each process loops over
    # its own rows, updates only the boundaries it owns, exchanges its edge rows before the swap,
    # and the rows are gathered after the loop. The step timers are left out so the loop fits.
    compose([code_pane(solver_mpi, "mpi, main changes from omp marked, no timers or comments", [
        ("RowRange owned_rows = rows_for_rank(num_variance_nodes, rank, num_ranks);", "RowRange owned_rows = rows_for_rank("),
        ("for (int step = 1; step <= num_timesteps; ++step) {", "gather_rows(grid, num_stock_nodes")],
        ["RowRange owned_rows", "var_j = first_interior", "if (owns_zero_row)", "owned_rows.begin, owned_rows.end",
         "if (owns_top_row)", "exchange_edge_rows(", "gather_rows("],
        leave_out=("omp_get_wtime", "result."), bare=True, start_after="SolveResult MpiSolver::solve(")],
        out_dir / "shot_mpi_loop.png")

    # Rung 4, mpi. One neighbour's receive and send, then the wait.
    compose([code_pane(solver_mpi, "mpi, one neighbour of the edge row exchange", [
        ("// The bottom process has nobody below it.", "requests.push_back(send_below);"),
        ("// Wait for every send and receive to finish.", "MPI_Waitall(")], ["MPI_Irecv", "MPI_Isend", "MPI_Waitall"])],
        out_dir / "shot_mpi_exchange.png")

    # Appendix shots, for questions about the less usual code.
    compose([code_pane(flush, "loaded with LD_PRELOAD by the bench scripts", [("/* Runs before main.", "}")],
                       ["__attribute__", "_mm_setcsr"])], out_dir / "shot_flush_denormals.png")
    compose([code_pane(stencil, "copied into locals before the loop", [
        ("// Copy the constants into local variables first.", "int num_stock_nodes = constants.num_stock_nodes;")],
        ["// Otherwise"])], out_dir / "shot_locals.png")
    compose([
        code_pane(grid_cpp, "the swap", [("// Swaps which buffer is read and which is written.", "}")], ["std::swap"]),
        code_pane(solver_omp, "bound once, before the time loop", [("const std::vector<double>& current = grid.current();", "std::vector<double>& next = grid.next();")])],
        out_dir / "shot_swap.png")
    compose([code_pane(solver_mpi, "which rows a process owns", [("static RowRange rows_for_rank(", "    return rows;")], ["extra_rows"])],
            out_dir / "shot_rows_for_rank.png")
    compose([code_pane(solver_mpi, "E9, inside the parallel region", [
        ("// Both rows to send are finished once every thread reaches this line.", "post_edge_rows(next, num_stock_nodes, owned_rows, rank, num_ranks, requests);")],
        ["#pragma omp barrier", "#pragma omp master", "post_edge_rows("])], out_dir / "shot_overlap.png")
    compose([code_pane(main_cpp, "MPI started for the overlap rung", [("// Every solver runs inside MPI.", "MPI_Init_thread(")], ["MPI_Init_thread("])],
            out_dir / "shot_init_thread.png")
    compose([code_pane(common, "runs before every sweep", [("build_and_validate() {", "make avx-check")], ["make test", "make avx-check"])],
            out_dir / "shot_iron_rule.png")

    # The header every benchmark CSV carries.
    csvs = sorted((ROOT / "results").glob("bench_m2_ladder_*.csv"))
    if csvs:
        header = "\n".join(line for line in csvs[-1].read_text().splitlines()
                           if line.startswith("#") and not line.startswith("# sweep,"))
        compose([text_pane(header, f"results/{csvs[-1].name}, header written by the job")], out_dir / "shot_csv_header.png")
    else:
        print("skip shot_csv_header: no results/bench_m2_ladder_*.csv on disk", file=sys.stderr)


if __name__ == "__main__":
    main()
