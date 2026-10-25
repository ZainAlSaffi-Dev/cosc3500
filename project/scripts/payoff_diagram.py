#!/usr/bin/env python3
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

"""The payoff diagram for the presentation's first Introduction beat.

One curve, max(S - K, 0), the strike kink annotated, today's share price marked.
Styled to match the benchmark figures (same palette as bench_plot.py) so
the deck reads as one set of figures rather than two.

Usage: .venv/bin/python scripts/payoff_diagram.py [--out results/payoff.png]
"""

import argparse

import matplotlib
matplotlib.use("Agg")  # file output only, no display needed
from figstyle import BLUE, GREY, DARK_GREY, TEXT, GRID
import figstyle
figstyle.apply()
import matplotlib.pyplot as plt
import numpy as np

STRIKE = 5250.0
SPOT = 5200.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/payoff.png")
    args = parser.parse_args()

    stock = np.linspace(4600, 5900, 400)
    payoff = np.maximum(stock - STRIKE, 0.0)

    # Drawn at the size of the figure column on its slide, so the text keeps its size
    fig, ax = plt.subplots(figsize=(4.0, 2.6), dpi=250)
    ax.plot(stock, payoff, color=BLUE, linewidth=2.6, zorder=3,
            solid_capstyle="round")

    # The two prices that define the story, the strike where the kink sits,
    # and today's share price just below it.
    ax.axvline(STRIKE, color=GREY, linestyle="--", linewidth=1.2, zorder=2)
    ax.annotate(f"strike K = ${STRIKE:,.0f}\nbelow it the call\npays nothing",
                xy=(STRIKE, 0), xytext=(4620, 400),
                fontsize=9, color=DARK_GREY,
                arrowprops=dict(arrowstyle="->", color=DARK_GREY, lw=1.0))
    ax.plot([SPOT], [0], "o", color=TEXT, markersize=7, zorder=4)
    ax.annotate(f"today's share\nprice ${SPOT:,.0f}",
                xy=(SPOT, 0), xytext=(4620, 130),
                fontsize=9, color=TEXT,
                arrowprops=dict(arrowstyle="->", color=TEXT, lw=1.0))

    ax.set_xlabel("share price at expiry ($)", fontsize=9.5, color=TEXT)
    ax.set_ylabel("payoff ($)", fontsize=9.5, color=TEXT)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(colors=DARK_GREY, labelsize=8.5)
    ax.set_ylim(-40, 700)

    fig.tight_layout(pad=0.4)
    fig.savefig(args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
