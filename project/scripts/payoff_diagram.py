#!/usr/bin/env python3
"""The payoff diagram for the presentation's first Introduction beat.

One curve, max(S − K, 0), the strike kink annotated, today's spot marked.
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

    fig, ax = plt.subplots(figsize=(9, 5), dpi=160)
    ax.plot(stock, payoff, color=BLUE, linewidth=3, zorder=3,
            solid_capstyle="round")

    # The two prices that define the story: the strike where the kink sits,
    # and today's spot just below it.
    ax.axvline(STRIKE, color=GREY, linestyle="--", linewidth=1.4, zorder=2)
    ax.annotate(f"strike K = ${STRIKE:,.0f}\nthe kink: below it the option"
                "\nexpires worthless",
                xy=(STRIKE, 0), xytext=(STRIKE + 60, 260),
                fontsize=11, color=DARK_GREY,
                arrowprops=dict(arrowstyle="->", color=DARK_GREY, lw=1.2))
    ax.plot([SPOT], [0], "o", color=TEXT, markersize=9, zorder=4)
    ax.annotate(f"today's spot S = ${SPOT:,.0f}",
                xy=(SPOT, 0), xytext=(4680, 120),
                fontsize=11, color=TEXT,
                arrowprops=dict(arrowstyle="->", color=TEXT, lw=1.2))

    ax.set_title("European call at expiry: the right, not the obligation, "
                 "to buy at K", fontsize=13, color=TEXT, pad=12)
    ax.set_xlabel("share price at expiry ($)", fontsize=11, color=TEXT)
    ax.set_ylabel("option payoff ($)", fontsize=11, color=TEXT)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(colors=DARK_GREY)
    ax.set_ylim(-40, 700)

    fig.tight_layout()
    fig.savefig(args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
