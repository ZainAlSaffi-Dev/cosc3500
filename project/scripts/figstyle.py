# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

"""Shared look for every figure in the project.

Each plotting script imports this and calls apply() once before drawing, so
the ladder, the scaling sweep, the roofline, the convergence study and the
payoff diagram all read as one set of figures. The rules are the usual ones
for quiet charts. The data carries the ink, the scaffolding recedes. Grid
lines are faint, the top and right spines are gone, colour is reserved for
meaning (blue for the thing that worked, grey for the things that did not,
red for deliberately worse code) and text sits in dark grey rather than in
the series colour.
"""

import matplotlib

# Colour constants. These predate this module and every script already uses
# these names, so they live here now and the scripts import them.
BLUE = "#2563EB"        # the optimised solver, the rung that worked
GREY = "#9CA3AF"        # baseline, guides, reference lines
LIGHT_GREY = "#CDD3DC"  # null rungs: present, honest, visually quiet
DARK_GREY = "#6B7280"   # annotation text
RED = "#DC2626"         # negative controls, deliberately worse code
AMBER = "#D97706"       # regressions and warnings
ORANGE = "#F97316"      # the memory roof in the roofline
TEXT = "#374151"        # titles, labels, values
GRID = "#ECEEF2"        # fainter than the old #E5E7EB on purpose
SPINE = "#D1D5DB"


def apply():
    """Set the shared rcParams. Call once, before any figure is created."""
    matplotlib.rcParams.update({
        "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"],
        "font.size": 11,
        "text.color": TEXT,
        "axes.titlesize": 14,
        "axes.titleweight": 600,
        "axes.titlecolor": TEXT,
        "axes.labelsize": 11,
        "axes.labelcolor": TEXT,
        "axes.edgecolor": SPINE,
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.7,
        "xtick.color": DARK_GREY,
        "ytick.color": DARK_GREY,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.frameon": False,
        "legend.fontsize": 9.5,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def apply_slide():
    """The shared look at slide size, for figures drawn at the size they appear in the deck.

    A figure drawn at 12 inches and shrunk into a 5.5 inch slot loses half its
    point size, which made the first Milestone 2 deck unreadable. The slide
    figures are drawn at their slot size instead, so these sizes are close to
    the sizes on the slide.
    """
    apply()
    matplotlib.rcParams.update({
        "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"],
        "font.size": 10.5,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "lines.markersize": 5,
    })


def finish_axes(*axes):
    """Drop the top and right spines. The rcParams cannot express this."""
    for ax in axes:
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
