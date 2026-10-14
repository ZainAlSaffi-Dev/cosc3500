# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""F10, who updates what. One picture per parallel rung on a small grid.

f10a_simd.png, four neighbouring cells of one row in one AVX2 register, and the
three overlapping windows each row is loaded as.
f10b_omp.png, OpenMP threads each owning a band of rows.
f10c_mpi.png, MPI ranks each owning a block of rows, with the edge rows swapped
between neighbours every step.
f10_decomposition.png, the three side by side, for documents.

Each slide picture is drawn at the size of the figure column on its slide.

    .venv/bin/python scripts/decomposition_diagram.py --out results/m2/f10_decomposition.png
"""

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import figstyle
from figstyle import BLUE, GREY, LIGHT_GREY, DARK_GREY, RED, ORANGE, TEXT

WORKER_COLOURS = [BLUE, ORANGE, "#10B981", "#8B5CF6"]
PALE_BLUE = "#BFD3FB"
# The figure column of a slide with text beside it, the same as HALF in m2_figures.py
HALF = (4.0, 2.6)
LABEL_FONT = 9
DPI = 250


def split(num_rows, workers):
    """First rows each worker owns, the first workers taking one extra, as the code does."""
    per = num_rows // workers
    extra = num_rows % workers
    bounds = []
    begin = 0
    for worker in range(workers):
        size = per + (1 if worker < extra else 0)
        bounds.append((begin, begin + size))
        begin += size
    return bounds


def draw_grid(ax, rows, columns, colour_of, y0=0.0):
    """One square per cell, coloured by colour_of(row, column), row 0 at the bottom."""
    for row in range(rows):
        for column in range(columns):
            ax.add_patch(Rectangle((column, y0 + row), 1, 1, facecolor=colour_of(row, column), edgecolor="white", lw=1.0))


def simd_picture(ax):
    """One row of the grid, four cells per register, and the three windows a row is loaded as."""
    columns = 11
    groups = [(1, 5), (5, 9)]

    def colour(row, column):
        if column == 0 or column == columns - 1:
            return GREY
        if row != 1:
            return LIGHT_GREY
        if groups[0][0] <= column < groups[0][1]:
            return BLUE
        if groups[1][0] <= column < groups[1][1]:
            return PALE_BLUE
        return RED

    y0 = 2.7
    draw_grid(ax, 3, columns, colour, y0=y0)
    ax.add_patch(Rectangle((1, y0 + 1), 4, 1, fill=False, edgecolor=TEXT, lw=2.2))
    for row, name in ((2, "row above"), (1, "this row"), (0, "row below")):
        ax.text(-0.3, y0 + row + 0.5, name, ha="right", va="center", fontsize=LABEL_FONT, color=TEXT)
    ax.text(3.0, y0 + 3.2, "one register, 4 cells", ha="center", va="bottom", fontsize=LABEL_FONT, color=TEXT, weight=600)
    ax.text(9.5, y0 + 3.2, "last cell,\nscalar", ha="center", va="bottom", fontsize=LABEL_FONT - 0.5, color=RED)
    # The three windows of four cells that this row is loaded as for the first group
    for offset, (start, name) in enumerate(((0, "west"), (1, "centre"), (2, "east"))):
        y = 1.85 - offset * 0.8
        ax.add_patch(Rectangle((start + 0.05, y), 3.9, 0.6, facecolor=PALE_BLUE, edgecolor=BLUE, lw=1.0))
        ax.text(start + 2.0, y + 0.3, name, ha="center", va="center", fontsize=LABEL_FONT - 0.5, color=TEXT)
    ax.text(6.6, 1.35, "each of the three rows\nis loaded as west, centre\nand east, so 9 overlapping\nloads per 4 cells",
            ha="left", va="center", fontsize=LABEL_FONT - 0.5, color=DARK_GREY)
    ax.set_xlim(-3.6, 14.6)
    ax.set_ylim(0.0, y0 + 4.4)
    ax.set_aspect("equal")
    ax.axis("off")


def omp_picture(ax):
    """The interior rows split into one band per thread, the boundaries in grey."""
    rows, columns, threads = 12, 16, 4
    bands = split(rows - 2, threads)

    def colour(row, column):
        if row == 0 or row == rows - 1 or column == 0 or column == columns - 1:
            return GREY
        for thread, (begin, end) in enumerate(bands):
            if begin <= row - 1 < end:
                return WORKER_COLOURS[thread]
        return LIGHT_GREY

    draw_grid(ax, rows, columns, colour)
    for thread, (begin, end) in enumerate(bands):
        ax.text(columns + 0.3, (begin + end) / 2 + 1, f"thread {thread}", va="center", fontsize=LABEL_FONT,
                color=WORKER_COLOURS[thread], weight=600)
    ax.text(columns + 0.3, 0.5, "boundaries,\nmain thread", va="center", fontsize=LABEL_FONT - 0.5, color=DARK_GREY)
    ax.text(columns + 0.3, rows - 0.5, "top row,\nmain thread", va="center", fontsize=LABEL_FONT - 0.5, color=DARK_GREY)
    ax.set_xlim(-0.3, columns + 5.2)
    ax.set_ylim(-0.3, rows + 0.3)
    ax.set_aspect("equal")
    ax.axis("off")


def mpi_picture(ax):
    """Every row split into one block per rank, with the two rows at each block edge swapped every step."""
    rows, columns, ranks = 12, 16, 3
    blocks = split(rows, ranks)

    def colour(row, column):
        for rank, (begin, end) in enumerate(blocks):
            if begin <= row < end:
                return WORKER_COLOURS[rank]
        return LIGHT_GREY

    draw_grid(ax, rows, columns, colour)
    for rank, (begin, end) in enumerate(blocks):
        ax.text(columns + 0.3, (begin + end) / 2, f"rank {rank}", va="center", fontsize=LABEL_FONT,
                color=WORKER_COLOURS[rank], weight=600)
        if rank > 0:
            ax.plot([0, columns], [begin, begin], color=TEXT, lw=2.4)
            # The last row of the block below and the first row of this block trade places
            ax.add_patch(FancyArrowPatch((-0.6, begin - 1.0), (-0.6, begin + 1.0), arrowstyle="<->",
                                         mutation_scale=14, color=RED, lw=2.0))
    ax.text(-1.1, rows / 2, "edge rows\nexchanged\nevery step", ha="right", va="center", fontsize=LABEL_FONT - 0.5, color=RED)
    ax.set_xlim(-5.2, columns + 3.6)
    ax.set_ylim(-0.3, rows + 0.3)
    ax.set_aspect("equal")
    ax.axis("off")


def save(fig, path):
    # Cropped to the drawing, so it fills its column on the slide
    fig.tight_layout(pad=0.3)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/m2/f10_decomposition.png")
    args = parser.parse_args()
    figstyle.apply_slide()
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    for name, draw in (("f10a_simd.png", simd_picture), ("f10b_omp.png", omp_picture), ("f10c_mpi.png", mpi_picture)):
        fig, ax = plt.subplots(figsize=HALF)
        draw(ax)
        save(fig, out.parent / name)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.2))
    for ax, draw in zip(axes, (simd_picture, omp_picture, mpi_picture)):
        draw(ax)
    save(fig, out)


if __name__ == "__main__":
    main()
