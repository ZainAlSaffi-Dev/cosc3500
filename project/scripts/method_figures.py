# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""F13, F14 and F18, the method figures for a viewer who has not seen Milestone 1.

F13, one time step. The current sheet at expiry holds the payoff, and the
nine cells around one cell give that cell's value in the next sheet, one dt
earlier. The next sheet fills row by row and its edges hold known values.
F14, the Black-Scholes check. With xi = 0 and rho = 0 the variance never
moves and Heston becomes Black-Scholes. The solver runs its full stencil on
results/collapse_work and its price today is drawn over the closed form.
F18, the value map. The time value over share price and variance half way
through the solve of config/demo.cfg, from its snapshots in results/demo.

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


# Slide sizes in inches, the same as WIDE and HALF in m2_figures.py
WIDE = (6.6, 2.75)
HALF = (4.0, 2.6)
DEMO_SPOT = SPOT


def grid_figure(out_dir):
    """One time step, the current buffer at expiry feeding one cell of the next buffer."""
    columns = 8
    rows = 6
    centre_col = 3
    centre_row = 2
    gap = 5.0
    next_left = columns + gap
    stock_values = np.linspace(0.6 * STRIKE, 1.5 * STRIKE, columns)
    variance_values = np.linspace(0.01, 0.12, rows)
    payoff_shades = option_shades(stock_values, variance_values, 0.0)

    fig, ax = plt.subplots(figsize=WIDE)

    # Current buffer, t = T, every cell holds the payoff
    for row in range(rows):
        for col in range(columns):
            cell(ax, col, row, payoff_shades[row][col])
    ax.add_patch(Rectangle((centre_col - 1, centre_row - 1), 3, 3, fill=False, edgecolor=RED, lw=2.2))
    ax.text(columns / 2, rows + 1.0, "current, t = T", ha="center", fontsize=10.5, weight=600)
    ax.text(columns / 2, rows + 0.25, "every cell holds the payoff", ha="center", fontsize=9, color=DARK_GREY)
    for offset, name in ((-1, "i-1"), (0, "i"), (1, "i+1")):
        ax.text(centre_col + offset + 0.5, -0.25, name, ha="center", va="top", fontsize=8.5, color=TEXT)
    for offset, name in ((-1, "j-1"), (0, "j"), (1, "j+1")):
        ax.text(-0.15, centre_row + offset + 0.5, name, ha="right", va="center", fontsize=8.5, color=TEXT)
    ax.add_patch(FancyArrowPatch((0, -1.25), (columns, -1.25), arrowstyle="-|>", mutation_scale=12, color=TEXT, lw=1.1))
    ax.text(columns / 2, -1.5, "share price S", ha="center", va="top", fontsize=9, color=TEXT)
    ax.add_patch(FancyArrowPatch((-1.5, 0), (-1.5, rows), arrowstyle="-|>", mutation_scale=12, color=TEXT, lw=1.1))
    ax.text(-1.75, rows / 2, "variance v", ha="right", va="center", rotation=90, fontsize=9, color=TEXT)

    # Next buffer, t = T - dt, filled row by row
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
    ax.add_patch(Rectangle((next_left + centre_col, centre_row), 1, 1, fill=False, edgecolor=RED, lw=2.2))
    ax.text(next_left + columns / 2, rows + 1.0, "next, t = T - dt", ha="center", fontsize=10.5, weight=600)
    ax.text(next_left + columns / 2, rows + 0.25, "grey edges are known values", ha="center", fontsize=9, color=DARK_GREY)
    ax.text(next_left + columns / 2, -1.5, "then the buffers swap and the step repeats", ha="center", va="top", fontsize=9, color=DARK_GREY)

    # The nine cells of the current buffer give the one new cell at the same i and j
    arrow_y = rows - 0.9
    ax.add_patch(FancyArrowPatch((columns + 0.3, arrow_y), (next_left - 0.3, arrow_y), arrowstyle="-|>", mutation_scale=14, color=TEXT, lw=1.4))
    ax.text(columns + gap / 2, arrow_y + 0.25, "one step back", ha="center", va="bottom", fontsize=9, color=TEXT)

    # Key, which derivative each of the nine cells feeds
    key_left = columns + gap / 2 - 1.5
    key_bottom = 0.3
    key_labels = {(-1, -1): "$V_{Sv}$", (1, -1): "$V_{Sv}$", (-1, 1): "$V_{Sv}$", (1, 1): "$V_{Sv}$",
                  (-1, 0): "$V_S$", (1, 0): "$V_S$", (0, -1): "$V_v$", (0, 1): "$V_v$", (0, 0): "$V$"}
    for (col_offset, row_offset), label in key_labels.items():
        ax.add_patch(Rectangle((key_left + 1 + col_offset, key_bottom + 1 + row_offset), 1, 1, facecolor="white", edgecolor=SPINE, lw=0.8))
        ax.text(key_left + 1.5 + col_offset, key_bottom + 1.5 + row_offset, label, ha="center", va="center", fontsize=8, color=TEXT)
    ax.add_patch(Rectangle((key_left, key_bottom), 3, 3, fill=False, edgecolor=RED, lw=2.0))
    ax.text(columns + gap / 2, key_bottom - 0.2, "the nine-point stencil,\nand where each cell is used", ha="center", va="top", fontsize=8.5, color=RED)

    ax.set_xlim(-2.4, next_left + columns + 0.3)
    ax.set_ylim(-2.2, rows + 1.6)
    ax.set_aspect("equal")
    ax.axis("off")
    out = out_dir / "f13_grid.png"
    fig.savefig(out, dpi=250, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"wrote {out}")


def collapse_figure(out_dir, dump_dir):
    """Today's PDE slice at v = v0 over the closed form, and the gap between them."""
    snapshots = sorted(dump_dir.glob("snap_*.csv"))
    if not snapshots:
        print(f"skip f14_bs_collapse: no snapshots in {dump_dir}, run make collapse-anim first", file=sys.stderr)
        return
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

    fig, (top, bottom) = plt.subplots(2, 1, figsize=WIDE, sharex=True, gridspec_kw={"height_ratios": [2.4, 1.0]})
    s = stock_axis[window]
    top.plot(s, payoff[window], ":", color=DARK_GREY, lw=1.4, label="payoff at expiry")
    top.plot(s, closed_form[window], color=GREY, lw=6, alpha=0.55, solid_capstyle="round", label="Black-Scholes formula")
    top.plot(s, pde[window], color=BLUE, lw=1.6, label="solver with xi = 0, today")
    top.axvline(STRIKE, color=SPINE, lw=0.8)
    top.plot(SPOT, pde_price, "o", color=RED, ms=5)
    top.annotate(f"S = 5200, solver {pde_price:.2f},\nformula {bs_price:.2f}", (SPOT, pde_price), xytext=(SPOT + 650, pde_price - 60),
                 fontsize=8.5, color=TEXT, va="top", arrowprops=dict(arrowstyle="-", color=DARK_GREY, lw=0.8))
    top.set_ylabel("value ($)")
    top.legend(loc="upper left", fontsize=8.5)

    bottom.plot(s, (pde - closed_form)[window], color=BLUE, lw=1.3)
    bottom.axhline(0.0, color=GREY, lw=0.8)
    bottom.set_ylim(-0.5, 0.5)
    bottom.set_ylabel("gap ($)")
    bottom.set_xlabel("share price S")
    bottom.text(s[-1], 0.25, f"largest gap {max_gap:.3f}", ha="right", fontsize=8.5, color=DARK_GREY)
    figstyle.finish_axes(top, bottom)
    fig.tight_layout(pad=0.4)
    out = out_dir / "f14_bs_collapse.png"
    fig.savefig(out, dpi=250)
    plt.close(fig)
    print(f"wrote {out}")


def value_map_figure(out_dir, demo_dir):
    """F18, the time value over share price and variance half way through the solve of config/demo.cfg."""
    snapshots = sorted(demo_dir.glob("snap_*.csv"))
    if not snapshots:
        print(f"skip f18_value_map: no snapshots in {demo_dir}, run the demo solve first", file=sys.stderr)
        return
    # demo.cfg has the same contract and grid as the collapse run, with xi = 0.35 and rho = -0.7
    path = snapshots[max(len(snapshots) // 2 - 1, 0)]
    values = np.loadtxt(path, delimiter=",")
    stock = np.linspace(0.0, S_MAX, values.shape[1])
    variance = np.linspace(0.0, V_MAX, values.shape[0])
    # The time value, what the option is worth above what it would pay now. It shows how value
    # grows with variance near the strike, which the value itself hides under the payoff
    values = values - np.maximum(stock - STRIKE, 0.0)[np.newaxis, :]
    s_window = (stock >= 3600) & (stock <= 7400)
    v_window = variance <= 0.16
    fig, ax = plt.subplots(figsize=HALF)
    mesh = ax.pcolormesh(stock[s_window], variance[v_window], values[np.ix_(v_window, s_window)], cmap="viridis", shading="auto")
    colourbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    colourbar.set_label("time value, V minus payoff ($)", fontsize=9)
    colourbar.ax.tick_params(labelsize=8.5)
    ax.axvline(STRIKE, color="white", lw=0.9, ls="--")
    ax.text(STRIKE + 70, 0.15, "strike", color="white", fontsize=8.5, va="top")
    ax.plot(DEMO_SPOT, V0, "o", mfc="none", mec="white", mew=1.6, ms=8)
    ax.annotate("today's share price\nand variance", (DEMO_SPOT, V0), xytext=(DEMO_SPOT - 1500, V0 + 0.055), color="white",
                fontsize=8.5, arrowprops=dict(arrowstyle="-", color="white", lw=0.9))
    ax.set_xlabel("share price S")
    ax.set_ylabel("variance v")
    ax.grid(False)
    fig.tight_layout(pad=0.4)
    out = out_dir / "f18_value_map.png"
    fig.savefig(out, dpi=250)
    plt.close(fig)
    print(f"wrote {out} from {path.name}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/m2")
    parser.add_argument("--dump", default="results/collapse_work/dump")
    parser.add_argument("--demo", default="results/demo", help="snapshots of config/demo.cfg, for the value map")
    args = parser.parse_args()
    figstyle.apply_slide()
    out_dir = pathlib.Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    grid_figure(out_dir)
    collapse_figure(out_dir, pathlib.Path(args.dump))
    value_map_figure(out_dir, pathlib.Path(args.demo))


if __name__ == "__main__":
    main()
