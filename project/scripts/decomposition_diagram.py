# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""F10, who updates what. One picture per parallel rung on a small grid.

Left, four neighbouring cells of one row in one AVX2 register.
Middle, OpenMP threads each owning a band of rows.
Right, MPI ranks each owning a slab, with the halo rows copied from the neighbours.

    arch -arm64 .venv/bin/python scripts/decomposition_diagram.py --out results/m2/f10_decomposition.png
"""

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import figstyle
from figstyle import BLUE, GREY, LIGHT_GREY, DARK_GREY, RED, ORANGE, TEXT

COLUMNS = 16
ROWS = 12
WORKER_COLOURS = [BLUE, ORANGE, "#10B981", "#8B5CF6"]


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


def draw_grid(ax, colour_of):
    for row in range(ROWS):
        for col in range(COLUMNS):
            ax.add_patch(Rectangle((col, row), 1, 1, facecolor=colour_of(row, col), edgecolor="white", lw=0.8))
    ax.set_xlim(-0.5, COLUMNS + 0.5)
    ax.set_ylim(-1.5, ROWS + 0.5)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/m2/f10_decomposition.png")
    args = parser.parse_args()
    figstyle.apply()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))

    # SIMD, one row, cells 1 to 4 in one register and the scalar remainder at the end
    ax = axes[0]
    simd_row = 5

    def simd_colour(row, col):
        if col == 0 or col == COLUMNS - 1:
            return GREY
        if row != simd_row:
            return LIGHT_GREY
        if col <= 12:
            return [BLUE, ORANGE, "#10B981"][(col - 1) // 4]
        return RED

    draw_grid(ax, simd_colour)
    for start in (1, 5, 9):
        ax.add_patch(Rectangle((start, simd_row), 4, 1, fill=False, edgecolor=TEXT, lw=2))
    ax.text(COLUMNS / 2, ROWS + 0.1, "SIMD, one AVX2 register = 4 neighbouring cells", ha="center", fontsize=11, weight=600)
    ax.text(COLUMNS / 2, -1.0, "boxed groups of 4 share one instruction, red cells are the scalar remainder,\ngrey columns are the stock edges", ha="center", fontsize=8.5, color=DARK_GREY)

    # OpenMP, interior rows 1 to ROWS - 2 split between 4 threads
    ax = axes[1]
    bands = split(ROWS - 2, 4)

    def omp_colour(row, col):
        if row == 0 or row == ROWS - 1 or col == 0 or col == COLUMNS - 1:
            return GREY
        for thread, (begin, end) in enumerate(bands):
            if begin <= row - 1 < end:
                return WORKER_COLOURS[thread]
        return LIGHT_GREY

    draw_grid(ax, omp_colour)
    for thread, (begin, end) in enumerate(bands):
        ax.text(COLUMNS + 0.2, (begin + end) / 2 + 1, f"thread {thread}", va="center", fontsize=8.5, color=WORKER_COLOURS[thread])
    ax.text(COLUMNS / 2, ROWS + 0.1, "OpenMP, each thread owns a band of rows", ha="center", fontsize=11, weight=600)
    ax.text(COLUMNS / 2, -1.0, "all threads read the old sheet and write only their own rows,\ngrey boundaries are done by the main thread after the barrier", ha="center", fontsize=8.5, color=DARK_GREY)

    # MPI, every row split between 3 ranks, halo rows outlined
    ax = axes[2]
    slabs = split(ROWS, 3)

    def mpi_colour(row, col):
        for rank, (begin, end) in enumerate(slabs):
            if begin <= row < end:
                return WORKER_COLOURS[rank]
        return LIGHT_GREY

    draw_grid(ax, mpi_colour)
    for rank, (begin, end) in enumerate(slabs):
        ax.text(COLUMNS + 0.2, (begin + end) / 2, f"rank {rank}", va="center", fontsize=8.5, color=WORKER_COLOURS[rank])
        if rank > 0:
            ax.plot([0, COLUMNS], [begin, begin], color=TEXT, lw=2.5)
            ax.annotate("", xy=(-0.3, begin - 0.5), xytext=(-0.3, begin + 0.5),
                        arrowprops=dict(arrowstyle="<->", color=RED, lw=1.5))
    ax.text(COLUMNS / 2, ROWS + 0.1, "MPI, each rank owns a slab of rows", ha="center", fontsize=11, weight=600)
    ax.text(COLUMNS / 2, -1.0, "each step the two rows either side of a black line are swapped,\nso every rank has the one row it reads from each neighbour", ha="center", fontsize=8.5, color=DARK_GREY)

    fig.tight_layout()
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
