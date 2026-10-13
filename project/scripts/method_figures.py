# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""F13 and F14, the two method figures for a viewer who has not seen Milestone 1.

F13, one time step. The current sheet at expiry holds the payoff, and the
nine cells around one cell give that cell's value in the next sheet, one dt
earlier. The next sheet fills row by row and its edges hold known values.
F14, the Black-Scholes check. With xi = 0 and rho = 0 the variance never
moves and Heston becomes Black-Scholes. The solver runs its full stencil on
results/collapse_work and its price today is drawn over the closed form.

    arch -arm64 .venv/bin/python scripts/method_figures.py --out results/m2
"""

import argparse
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import figstyle
from figstyle import BLUE, GREY, LIGHT_GREY, DARK_GREY, RED, TEXT, SPINE
from bs_collapse_anim import bs_call, read_variance_row

PALE_BLUE = "#BFD3FB"

# The contract and grid of results/collapse_work/collapse.cfg
STRIKE = 5250.0
MATURITY = 0.25
SPOT = 5200.0
RATE = 0.045
DIV_YIELD = 0.013
V0 = 0.04
NUM_STOCK = 526
NUM_VARIANCE = 101
S_MAX = 4.0 * STRIKE
V_MAX = 0.2


def cell(ax, col, row, colour, label="", label_colour=TEXT, size=1.0, fontsize=9):
    ax.add_patch(Rectangle((col * size, row * size), size, size, facecolor=colour, edgecolor="white", lw=1.0))
    if label:
        ax.text((col + 0.5) * size, (row + 0.5) * size, label, ha="center", va="center", fontsize=fontsize, color=label_colour)


def option_shades(stock_values, variance_values, time_to_expiry):
    """Colour for every cell of a small sheet, from the closed form at that time."""
    colours = matplotlib.colormaps["Blues"]
    largest = 0.5 * STRIKE
    shades = []
    for variance in variance_values:
        values = bs_call(stock_values, STRIKE, RATE, DIV_YIELD, math.sqrt(variance), time_to_expiry)
        shades.append([colours(0.12 + 0.8 * min(value / largest, 1.0)) for value in values])
    return shades


def grid_figure(out_dir):
    """One time step, the current sheet at expiry feeding one cell of the next sheet."""
    columns = 10
    rows = 7
    centre_col = 4
    centre_row = 3
    gap = 4.6
    next_left = columns + gap
    stock_values = np.linspace(0.6 * STRIKE, 1.5 * STRIKE, columns)
    variance_values = np.linspace(0.01, 0.12, rows)
    payoff_shades = option_shades(stock_values, variance_values, 0.0)

    fig, ax = plt.subplots(figsize=(11, 4.3))

    # Current sheet, t = T, every cell holds the payoff
    for row in range(rows):
        for col in range(columns):
            cell(ax, col, row, payoff_shades[row][col])
    ax.add_patch(Rectangle((centre_col - 1, centre_row - 1), 3, 3, fill=False, edgecolor=RED, lw=2.2))
    ax.text(2.0, rows - 0.5, "S < K, payoff 0", ha="center", va="center", fontsize=8.5, color=TEXT)
    ax.text(8.0, rows - 0.5, "payoff S - K", ha="center", va="center", fontsize=8.5, color="white")
    ax.text(columns / 2, rows + 1.05, "current, t = T", ha="center", fontsize=12, weight=600)
    ax.text(columns / 2, rows + 0.3, "known, every cell holds the payoff max(S - K, 0)", ha="center", fontsize=10, color=DARK_GREY)
    for offset, name in ((-1, "i-1"), (0, "i"), (1, "i+1")):
        ax.text(centre_col + offset + 0.5, -0.3, name, ha="center", va="top", fontsize=10, color=TEXT)
    for offset, name in ((-1, "j-1"), (0, "j"), (1, "j+1")):
        ax.text(-0.2, centre_row + offset + 0.5, name, ha="right", va="center", fontsize=10, color=TEXT)
    ax.add_patch(FancyArrowPatch((0, -1.2), (columns, -1.2), arrowstyle="-|>", mutation_scale=14, color=TEXT, lw=1.2))
    ax.text(columns / 2, -1.45, "stock price S, column i", ha="center", va="top", fontsize=10.5, color=TEXT)
    ax.add_patch(FancyArrowPatch((-1.4, 0), (-1.4, rows), arrowstyle="-|>", mutation_scale=14, color=TEXT, lw=1.2))
    ax.text(-1.65, rows / 2, "variance v, row j", ha="right", va="center", rotation=90, fontsize=10.5, color=TEXT)

    # Next sheet, t = T - dt, filled row by row along each row
    step_shades = option_shades(stock_values, variance_values, 0.02 * MATURITY)
    for row in range(rows):
        for col in range(columns):
            on_edge = row == 0 or row == rows - 1 or col == 0 or col == columns - 1
            written = row < centre_row or (row == centre_row and col < centre_col)
            if row == centre_row and col == centre_col:
                colour = BLUE
            elif on_edge:
                colour = GREY
            elif written:
                colour = step_shades[row][col]
            else:
                colour = LIGHT_GREY
            ax.add_patch(Rectangle((next_left + col, row), 1, 1, facecolor=colour, edgecolor="white", lw=1.0))
    ax.text(next_left + columns / 2, rows + 1.05, "next, t = T - dt", ha="center", fontsize=12, weight=600)
    ax.text(next_left + columns / 2, rows + 0.3, "filled row by row, along each row", ha="center", fontsize=10, color=DARK_GREY)
    ax.text(next_left + columns / 2, 0.5, "v = 0, drift terms only", ha="center", va="center", fontsize=8, color=TEXT)
    ax.text(next_left + columns / 2, rows - 0.5, "top v, copy of the row below", ha="center", va="center", fontsize=8, color=TEXT)
    ax.text(next_left + 0.5, rows / 2, "S = 0, value 0", ha="center", va="center", rotation=90, fontsize=8, color=TEXT)
    ax.text(next_left + columns - 0.5, rows / 2, "largest S, known", ha="center", va="center", rotation=-90, fontsize=8, color=TEXT)

    # The nine cells of the current sheet give the one new cell at the same i and j
    ax.add_patch(Rectangle((next_left + centre_col, centre_row), 1, 1, fill=False, edgecolor=RED, lw=2.2))
    arrow_y = rows - 1.6
    ax.add_patch(FancyArrowPatch((columns + 0.4, arrow_y), (next_left - 0.4, arrow_y), arrowstyle="-|>", mutation_scale=16, color=TEXT, lw=1.6))
    ax.text(columns + gap / 2, arrow_y + 0.3, "one step back, dt", ha="center", va="bottom", fontsize=10.5, color=TEXT)
    ax.text(columns + gap / 2, arrow_y - 0.3, "red block of nine\ngives the red cell", ha="center", va="top", fontsize=10, color=RED)

    # Key, which derivative each of the nine cells feeds
    key_left = columns + gap / 2 - 1.5
    key_bottom = 0.0
    key_labels = {(-1, -1): "$V_{Sv}$", (1, -1): "$V_{Sv}$", (-1, 1): "$V_{Sv}$", (1, 1): "$V_{Sv}$",
                  (-1, 0): "$V_S$\n$V_{SS}$", (1, 0): "$V_S$\n$V_{SS}$",
                  (0, -1): "$V_v$\n$V_{vv}$", (0, 1): "$V_v$\n$V_{vv}$", (0, 0): "all"}
    for (col_offset, row_offset), label in key_labels.items():
        ax.add_patch(Rectangle((key_left + 1 + col_offset, key_bottom + 1 + row_offset), 1, 1, facecolor="white", edgecolor=SPINE, lw=0.8))
        ax.text(key_left + 1.5 + col_offset, key_bottom + 1.5 + row_offset, label, ha="center", va="center", fontsize=7.5, color=TEXT, linespacing=1.1)
    ax.add_patch(Rectangle((key_left, key_bottom), 3, 3, fill=False, edgecolor=RED, lw=2.2))
    ax.text(columns + gap / 2, key_bottom - 0.3, "the derivatives\neach cell feeds", ha="center", va="top", fontsize=9, color=DARK_GREY)

    # The update for one cell, then the two swap and the step repeats back to today
    ax.text(next_left + columns / 2, -0.35, "$V^{new}_{i,j} = V_{i,j} + dt \\times$ (the PDE terms from the nine cells)", ha="center", va="top", fontsize=10.5, color=TEXT)
    ax.text(next_left + columns / 2, -1.15, "then current and next swap, and the same step\nrepeats back to t = 0, where the price is read", ha="center", va="top", fontsize=10, color=DARK_GREY)

    ax.set_xlim(-2.4, next_left + columns + 0.4)
    ax.set_ylim(-2.9, rows + 2.0)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    out = out_dir / "f13_grid.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", pad_inches=0.1)
    print(f"wrote {out}")


def collapse_figure(out_dir, dump_dir):
    """Today's PDE slice at v = v0 over the closed form, and the gap between them."""
    snapshots = sorted(dump_dir.glob("snap_*.csv"))
    if not snapshots:
        raise SystemExit(f"no snapshots in {dump_dir}, run make bs-collapse first")
    variance_spacing = V_MAX / (NUM_VARIANCE - 1)
    v0_row = int(round(V0 / variance_spacing))
    stock_axis = np.linspace(0.0, S_MAX, NUM_STOCK)
    # The last snapshot is the final step, which is today
    pde = read_variance_row(snapshots[-1], v0_row)
    closed_form = bs_call(stock_axis, STRIKE, RATE, DIV_YIELD, math.sqrt(V0), MATURITY)
    payoff = np.maximum(stock_axis - STRIKE, 0.0)

    window = (stock_axis >= 0.7 * STRIKE) & (stock_axis <= 1.4 * STRIKE)
    spot_col = int(round(SPOT / (S_MAX / (NUM_STOCK - 1))))
    pde_price = pde[spot_col]
    bs_price = closed_form[spot_col]
    rel_err = abs(pde_price - bs_price) / bs_price
    max_gap = np.max(np.abs(pde - closed_form)[window])
    print(f"at S = {SPOT:.0f}: PDE {pde_price:.4f}  BS {bs_price:.4f}  relative {rel_err:.2e}  max gap in window {max_gap:.4f}")

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True, gridspec_kw={"height_ratios": [2.6, 1.0]})
    s = stock_axis[window]
    top.plot(s, payoff[window], ":", color=DARK_GREY, lw=1.4, label="payoff at expiry")
    top.plot(s, closed_form[window], color=GREY, lw=7, alpha=0.55, solid_capstyle="round", label="Black-Scholes closed form")
    top.plot(s, pde[window], color=BLUE, lw=1.8, label="Heston solver with xi = 0, today")
    top.axvline(STRIKE, color=SPINE, lw=0.8)
    top.text(STRIKE + 40, top.get_ylim()[1] * 0.92, "strike 5250", fontsize=9, color=DARK_GREY)
    top.plot(SPOT, pde_price, "o", color=RED, ms=6)
    top.annotate(f"at S = 5200\nsolver {pde_price:.2f}, closed form {bs_price:.2f}", (SPOT, pde_price), xytext=(SPOT - 1500, pde_price + 380),
                 fontsize=9.5, color=TEXT, arrowprops=dict(arrowstyle="-", color=DARK_GREY, lw=0.8))
    top.set_ylabel("option value, dollars")
    top.legend(loc="upper left")
    top.set_title("With no volatility of volatility the solver lands on Black-Scholes")

    bottom.plot(s, (pde - closed_form)[window], color=BLUE, lw=1.4)
    bottom.axhline(0.0, color=GREY, lw=0.8)
    bottom.set_ylim(-0.5, 0.5)
    bottom.set_ylabel("solver minus\nclosed form, $")
    bottom.set_xlabel("stock price S")
    bottom.text(s[-1], 0.32, f"largest gap {max_gap:.3f} dollars", ha="right", fontsize=9.5, color=DARK_GREY)
    figstyle.finish_axes(top, bottom)
    fig.tight_layout()
    out = out_dir / "f14_bs_collapse.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/m2")
    parser.add_argument("--dump", default="results/collapse_work/dump")
    args = parser.parse_args()
    figstyle.apply()
    out_dir = pathlib.Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    grid_figure(out_dir)
    collapse_figure(out_dir, pathlib.Path(args.dump))


if __name__ == "__main__":
    main()
