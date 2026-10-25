#!/usr/bin/env python3
# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

"""The animated diagrams and charts of the Quarto deck, presentation/deck_m2.qmd.

Each one is an inline SVG. The steps are reveal.js fragments, so the presenter
moves them on with the next key, and the charts draw themselves when their
slide opens. The deck holds each SVG between two marker comments, such as
<!-- svg:ladder --> and <!-- /svg:ladder -->, and this script rewrites what is
between them, so it can be run again after any change. The chart numbers come
from docs/RESULTS.md, Milestone 2. The full run diagram reads the solver
snapshots of the test grid in presentation/anim_m2/snapshots/dump, which git
ignores, so make them first as the comment at the top of
presentation/anim_m2/scenes.py says.

Usage: python3 scripts/deck_m2_svgs.py [--deck presentation/deck_m2.qmd]
"""

import argparse
import csv
import math
import os
import random
import re

# The deck palette, the same values as presentation/deck_m2.scss.
BLUE = "#2563EB"
LIGHT_BLUE = "#93C5FD"
PALE_BLUE = "#BFDBFE"
INK = "#1F2937"
SLATE = "#6B7280"
LINE = "#D1D5DB"
GREY = "#9CA3AF"
RED = "#DC2626"
AMBER = "#D97706"
UNKNOWN = "#CBD5E1"

# Payoff colours of the 8 columns in the grid diagram, cheap on the left.
PAYOFF_COLOURS = ["#DCE8F5", "#DCE8F5", "#DCE8F5", "#DCE8F5", "#B4D3EA", "#66A9D6", "#2B78BB", "#0B4A8C"]

# The ladder, job 659238, million cell updates per second on 2048 x 512.
LADDER = [("baseline", "", 57, "1x", ""), ("serial", "", 159, "2.8x", "x2.80"), ("simd", "", 456, "8.0x", "x2.87"),
          ("omp", "24 threads", 6367, "112x", "x14.0"), ("mpi", "2 x 12", 6984, "123x", "x1.10")]

# The roofline, jobs 651302, 659238 and 659241.
FLOPS_PER_BYTE = 31 / 24
SCALAR_PEAK = 5.12
AVX2_PEAK_CORE = 43.4
AVX2_PEAK_NODE = 971
TRIAD_CORE = 12.4
TRIAD_NODE = 67.9
SERIAL_GFLOPS = 4.9
SIMD_GFLOPS = 14.1
SIMD_GFLOPS_BIG = 12.1
OMP_GFLOPS = 195
OMP_GFLOPS_BIG = 93


def fragment(index=None, kind=""):
    # The class and step number of one reveal.js fragment.
    classes = ("fragment " + kind).strip()
    if index is None:
        return f'class="{classes}"'
    return f'class="{classes}" data-fragment-index="{index}"'


def rect(x, y, width, height, fill, stroke="#FFFFFF", stroke_width=3, extra=""):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{height:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" {extra}/>'


def text(x, y, words, size=22, anchor="middle", colour="", weight="400", extra=""):
    # colour is a class from deck_m2.scss, dim, red, blue or amber.
    colour_class = f' class="{colour}"' if colour else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}"{colour_class} {extra}>{words}</text>'


def line(x1, y1, x2, y2, colour, width=3, extra=""):
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{colour}" stroke-width="{width}" {extra}/>'


def arrow_head(name, colour):
    return (f'<defs><marker id="{name}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{colour}"/></marker></defs>')


def log_position(value, low, high, pixel_low, pixel_high):
    # Where value sits on a log axis that runs from low to high.
    fraction = (math.log10(value) - math.log10(low)) / (math.log10(high) - math.log10(low))
    return pixel_low + fraction * (pixel_high - pixel_low)


def svg_open(width, height, extra=""):
    return f'<svg class="anim" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" {extra}>'


# ---------------------------------------------------------------------------
# Diagrams with steps
# ---------------------------------------------------------------------------

def heston_paths(count, steps, seed):
    # Sample paths of the share price and the variance under the deck's Heston parameters.
    # Euler steps with the variance floored at zero, one step a trading day.
    generator = random.Random(seed)
    dt = 0.25 / steps
    paths = []
    for _ in range(count):
        price = 5200.0
        variance = 0.04
        prices = [price]
        variances = [variance]
        for _ in range(steps):
            shock_variance = generator.gauss(0, 1)
            shock_price = -0.7 * shock_variance + math.sqrt(1 - 0.7 * 0.7) * generator.gauss(0, 1)
            root = math.sqrt(max(variance, 0.0))
            price = price * math.exp((0.045 - 0.013 - 0.5 * variance) * dt + root * math.sqrt(dt) * shock_price)
            variance = max(variance + 1.5 * (0.04 - variance) * dt + 0.35 * root * math.sqrt(dt) * shock_variance, 0.0)
            prices.append(price)
            variances.append(variance)
        paths.append((prices, variances))
    return paths


def subscript(name, letters):
    # A derivative name such as V with SS as a subscript, for SVG text.
    return f'{name}<tspan baseline-shift="sub" font-size="70%">{letters}</tspan>'


def polyline(points, colour, width=3, extra=""):
    coordinates = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return f'<polyline points="{coordinates}" fill="none" stroke="{colour}" stroke-width="{width}" {extra}/>'


def closest_path(paths, target):
    # The sample path whose share price at expiry is nearest to target.
    return min(paths, key=lambda path: abs(path[0][-1] - target))


def option_paths():
    # One share price path that ends above the strike, one that ends below, then the payoff for every price.
    # Then eight more paths and the price today as the average discounted payoff.
    plot_left = 150
    plot_right = 800
    plot_top = 40
    plot_bottom = 400
    payoff_left = 850
    payoff_right = 1230
    low_price = 4400
    high_price = 6100
    price_y = lambda price: plot_bottom - (price - low_price) / (high_price - low_price) * (plot_bottom - plot_top)
    payoff_x = lambda price: payoff_left + max(price - 5250, 0) / 850 * (payoff_right - payoff_left)
    time_x = lambda k: plot_left + k / 63 * (plot_right - plot_left)
    strike_y = price_y(5250)
    out = [svg_open(1400, 560), arrow_head("option-arrow", INK), arrow_head("price-arrow", BLUE)]
    # The axes, the strike and today's price.
    out.append(line(plot_left, plot_bottom, plot_right + 10, plot_bottom, INK, 3, 'marker-end="url(#option-arrow)"'))
    out.append(line(plot_left, plot_bottom, plot_left, plot_top - 10, INK, 3, 'marker-end="url(#option-arrow)"'))
    out.append(text(plot_left - 16, plot_top + 10, "share price S", 20, "end"))
    out.append(text(plot_left, plot_bottom + 30, "today, S = 5200", 20))
    out.append(text(plot_right, plot_bottom + 30, "expiry, 3 months", 20))
    for price in (4800, 5600):
        out.append(text(plot_left - 12, price_y(price) + 7, str(price), 18, "end", "dim"))
    out.append(line(plot_left, strike_y, payoff_left, strike_y, RED, 2, 'stroke-dasharray="8 6"'))
    out.append(text(plot_left - 12, strike_y + 7, "K = 5250", 19, "end", "red"))
    out.append(f'<circle cx="{plot_left}" cy="{price_y(5200):.1f}" r="9" fill="{INK}"/>')
    out.append(line(payoff_left, plot_bottom, payoff_left, plot_top - 10, INK, 3, 'marker-end="url(#option-arrow)"'))
    out.append(line(payoff_left, plot_bottom, payoff_right + 10, plot_bottom, INK, 3, 'marker-end="url(#option-arrow)"'))
    out.append(text((payoff_left + payoff_right) / 2, plot_bottom + 30, "what the option pays", 20))

    def one_path(prices, colour, width):
        # The path, its end point and, above the strike, the bar of what it pays.
        end_y = price_y(prices[-1])
        drawn = polyline([(time_x(k), price_y(price)) for k, price in enumerate(prices)], colour, width)
        drawn += f'<circle cx="{plot_right}" cy="{end_y:.1f}" r="7" fill="{BLUE}"/>'
        if prices[-1] > 5250:
            drawn += line(plot_right, end_y, payoff_left, end_y, BLUE, 2, 'stroke-dasharray="4 5"')
            drawn += line(payoff_left, end_y, payoff_x(prices[-1]), end_y, BLUE, 7)
        else:
            drawn += line(plot_right, end_y, payoff_left, end_y, BLUE, 2, 'stroke-dasharray="4 5"')
            drawn += f'<circle cx="{payoff_left}" cy="{end_y:.1f}" r="7" fill="{BLUE}"/>'
        return drawn

    def caption(first, second, colour="blue", weight="400"):
        return text(plot_left, 470, first, 22, "start", colour, weight) + text(plot_left, 502, second, 22, "start", colour, weight)

    pool = heston_paths(4000, 63, 11)
    high = closest_path(pool, 5700)[0]
    low = closest_path(pool, 5000)[0]
    high_end = round(high[-1])
    low_end = round(low[-1])
    # Step 0, a path that ends above the strike.
    out.append(f'<g {fragment(0)}>' + one_path(high, BLUE, 3.5) + "</g>")
    out.append(f'<g {fragment(0, "current-visible")}>'
               + caption(f"ends at {high_end}, so buy at 5250 and sell at {high_end}",
                         f"the option pays {high_end} - 5250 = {high_end - 5250}") + "</g>")
    # Step 1, a path that ends below the strike.
    out.append(f'<g {fragment(1)}>' + one_path(low, BLUE, 3.5) + "</g>")
    out.append(f'<g {fragment(1, "current-visible")}>'
               + caption(f"ends at {low_end}, below the strike, so nobody buys at 5250", "the option pays 0") + "</g>")
    # Step 2, the payoff for every price at expiry.
    out.append(f'<g {fragment(2)}>'
               + polyline([(payoff_left, plot_bottom), (payoff_left, strike_y), (payoff_right, price_y(6100))], RED, 5)
               + text(payoff_left + 30, price_y(4900), "max(S - K, 0)", 24, "start", "red", "700") + "</g>")
    out.append(f'<g {fragment(2, "current-visible")}>'
               + caption("the payoff for every price at expiry", "S - K above the strike, 0 below it", "red") + "</g>")
    # Step 3, eight more paths.
    others = "".join(one_path(prices, LIGHT_BLUE, 2.5) for prices, _ in heston_paths(8, 63, 2))
    out.append(f'<g {fragment(3)}>' + others + "</g>")
    out.append(f'<g {fragment(3, "current-visible")}>'
               + caption("today nobody knows which path the share takes", "each path ends at its own payoff") + "</g>")
    # Step 4, the price today.
    out.append(f'<g {fragment(4)}>'
               + f'<path d="M {payoff_left - 16} {plot_top + 6} Q 480 {plot_top - 60} {plot_left + 10} {price_y(5200) - 10:.1f}" fill="none" stroke="{BLUE}" stroke-width="4" marker-end="url(#price-arrow)"/>'
               + caption("price today = the average payoff over all paths,", "discounted back to today, 196.17", "blue", "700") + "</g>")
    out.append("</svg>")
    return "\n".join(out)


def backwards_in_time(step):
    # Five small grids from today to expiry. On step `step` the expiry grid fills with the payoff
    # and the steps back from expiry towards today appear.
    cell = 30
    columns = 6
    rows = 4
    grid_width = cell * columns
    gap = 90
    left = 70
    top = 40
    grid_x = lambda k: left + k * (grid_width + gap)
    colours = ["#DCE8F5", "#DCE8F5", "#DCE8F5", "#66A9D6", "#2B78BB", "#0B4A8C"]
    out = [svg_open(1400, 230), arrow_head("back-arrow", BLUE)]
    for k in range(5):
        for r in range(rows):
            for c in range(columns):
                out.append(rect(grid_x(k) + c * cell, top + r * cell, cell, cell, UNKNOWN, "#FFFFFF", 2))
    out.append(text(grid_x(0) + grid_width / 2, top + rows * cell + 34, "today, t = 0", 21, "middle", "", "700"))
    out.append(text(grid_x(4) + grid_width / 2, top + rows * cell + 34, "expiry, t = T", 21, "middle", "", "700"))
    out.append(text(grid_x(2) + grid_width / 2, top + rows * cell + 34, "each grid holds V at every S and v", 19, "middle", "dim"))
    payoff = "".join(rect(grid_x(4) + c * cell, top + r * cell, cell, cell, colours[c], "#FFFFFF", 2) for r in range(rows) for c in range(columns))
    steps = [payoff, text(grid_x(4) + grid_width / 2, top - 12, "the payoff, known", 20, "middle", "red")]
    for k in range(4, 0, -1):
        start_x = grid_x(k) - 10
        end_x = grid_x(k - 1) + grid_width + 10
        steps.append(line(start_x, top + rows * cell / 2, end_x, top + rows * cell / 2, BLUE, 4, 'marker-end="url(#back-arrow)"'))
        steps.append(text((start_x + end_x) / 2, top + rows * cell / 2 - 14, "Δt", 20, "middle", "blue"))
    steps.append(text(grid_x(2) + grid_width / 2, top - 12, "solved one step back at a time", 20, "middle", "blue"))
    out.append(f'<g {fragment(step)}>' + "".join(steps) + "</g>")
    out.append("</svg>")
    return "\n".join(out)


# The colour of each PDE term and the derivative it uses, the same on every slide.
TERM_COLOURS = {"S": BLUE, "SS": AMBER, "v": "#059669", "vv": "#7C3AED", "Sv": "#DB2777", "r": "#0891B2"}


def finite_differences():
    # The cell at S = 5250, v = 0.04 of the test grid at expiry, with the value along its row and its 3 by 3 stencil.
    # Steps 0 to 3 work out the slope, the curvature, the two variance differences and the cross difference.
    curve_left = 70
    curve_right = 640
    axis_y = 330
    low_price = 5125
    high_price = 5375
    price_x = lambda price: curve_left + (price - low_price) / (high_price - low_price) * (curve_right - curve_left - 30)
    value_y = lambda value: axis_y - value * 2.0
    point = lambda price: (price_x(price), value_y(max(price - 5250, 0)))
    out = [svg_open(1400, 530), arrow_head("difference-arrow", INK)]
    out.append(text(curve_left, 30, "test grid at expiry, ΔS = 50, Δv = 0.02", 20, "start", "dim"))
    # The value along the row through the cell, the payoff max(S - 5250, 0).
    out.append(line(curve_left, axis_y, curve_right, axis_y, INK, 3, 'marker-end="url(#difference-arrow)"'))
    out.append(text(curve_right, axis_y + 30, "S", 22, "end"))
    out.append(text(curve_left, 70, "V along the row v = 0.04", 20, "start"))
    out.append(polyline([point(low_price), point(5250), point(high_price)], SLATE, 3))
    for price, name in ((5200, "west"), (5250, "centre"), (5300, "east")):
        x, y = point(price)
        out.append(line(x, axis_y, x, y, LINE, 2, 'stroke-dasharray="5 5"'))
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="9" fill="{INK}"/>')
        out.append(text(x, axis_y + 30, name, 20))
        out.append(text(x, axis_y + 54, str(price), 18, "middle", "dim"))
    # The stencil, with the share prices under the columns and the variances left of the rows.
    stencil_left = 860
    stencil_top = 60
    size = 104
    names = [["north west", "north", "north east"], ["west", "centre", "east"], ["south west", "south", "south east"]]
    values = ["0", "0", "50"]
    cell_x = lambda column: stencil_left + column * size
    cell_y = lambda row: stencil_top + row * size
    for row in range(3):
        for column in range(3):
            out.append(rect(cell_x(column), cell_y(row), size, size, "#FFFFFF", LINE, 2))
            out.append(text(cell_x(column) + size / 2, cell_y(row) + 24, names[row][column], 15, "middle", "dim"))
            out.append(text(cell_x(column) + size / 2, cell_y(row) + 70, values[column], 30, "middle", "", "700"))
    for column, price in enumerate(("5200", "5250", "5300")):
        out.append(text(cell_x(column) + size / 2, cell_y(3) + 28, price, 18, "middle", "dim"))
    for row, variance in enumerate(("0.06", "0.04", "0.02")):
        out.append(text(stencil_left - 12, cell_y(row) + size / 2 + 6, variance, 18, "end", "dim"))
    out.append(text(cell_x(1) + size / 2, cell_y(3) + 54, "every row holds the same payoff at expiry", 18, "middle", "dim"))

    def weights(cells, colour):
        # A coloured border and weight on each cell a difference reads. Each entry is row, column and weight.
        drawn = ""
        for row, column, weight in cells:
            drawn += rect(cell_x(column) + 3, cell_y(row) + 3, size - 6, size - 6, "none", colour, 5, 'rx="6"')
            drawn += text(cell_x(column) + size - 10, cell_y(row) + size - 10, weight, 20, "end", "", "700", f'style="fill:{colour}"')
        return drawn

    def formula(first, second, colour):
        return (text(curve_left, 476, first, 25, "start", "", "700", f'style="fill:{colour}"')
                + text(curve_left, 510, second, 21, "start", "dim"))

    west = point(5200)
    centre = point(5250)
    east = point(5300)
    # Step 0, the slope from west to east.
    out.append(f'<g {fragment(0)}>' + line(west[0] - 30, west[1] + 15, east[0] + 30, east[1] - 15, TERM_COLOURS["S"], 4) + "</g>")
    out.append(f'<g {fragment(0, "current-visible")}>' + weights([(1, 0, "-1"), (1, 2, "+1")], TERM_COLOURS["S"])
               + formula(f"slope {subscript('V', 'S')} ≈ (east - west) / 2ΔS = (50 - 0) / 100 = 0.5",
                         "the rise from west to east over two spacings", TERM_COLOURS["S"]) + "</g>")
    # Step 1, the change in slope across the centre.
    out.append(f'<g {fragment(1)}>' + line(west[0], west[1] - 8, centre[0], centre[1] - 8, TERM_COLOURS["SS"], 4)
               + line(centre[0], centre[1] - 8, east[0], east[1] - 8, TERM_COLOURS["SS"], 4)
               + text((west[0] + centre[0]) / 2, west[1] - 26, "slope 0", 18, "middle", "", "400", f'style="fill:{TERM_COLOURS["SS"]}"')
               + text((centre[0] + east[0]) / 2 - 40, (centre[1] + east[1]) / 2 - 20, "slope 1", 18, "middle", "", "400", f'style="fill:{TERM_COLOURS["SS"]}"')
               + "</g>")
    out.append(f'<g {fragment(1, "current-visible")}>' + weights([(1, 0, "+1"), (1, 1, "-2"), (1, 2, "+1")], TERM_COLOURS["SS"])
               + formula(f"curvature {subscript('V', 'SS')} ≈ (east - 2 centre + west) / ΔS² = 50 / 2500 = 0.02",
                         "the slope goes from 0 to 1 across one spacing of 50", TERM_COLOURS["SS"]) + "</g>")
    # Step 2, the two variance differences from north and south.
    out.append(f'<g {fragment(2, "current-visible")}>' + weights([(0, 1, "+1"), (1, 1, "-2"), (2, 1, "+1")], TERM_COLOURS["vv"])
               + formula(f"{subscript('V', 'vv')} ≈ (north - 2 centre + south) / Δv² = 0",
                         f"{subscript('V', 'v')} ≈ (north - south) / 2Δv = 0 the same way, from north and south alone", TERM_COLOURS["vv"]) + "</g>")
    # Step 3, the cross difference from the four corners.
    out.append(f'<g {fragment(3, "current-visible")}>' + weights([(0, 0, "-1"), (0, 2, "+1"), (2, 0, "+1"), (2, 2, "-1")], TERM_COLOURS["Sv"])
               + formula(f"{subscript('V', 'Sv')} ≈ (north east - north west - south east + south west) / 4ΔSΔv = 0",
                         "how the slope along S changes with v, so it needs the four corners", TERM_COLOURS["Sv"]) + "</g>")
    out.append("</svg>")
    return "\n".join(out)


def code_text(x, y, pieces, size=19):
    # One line of C++ in a monospace font. Each piece is words and a colour, or words alone in ink.
    spans = ""
    for piece in pieces:
        if isinstance(piece, tuple):
            spans += f'<tspan style="fill:{piece[1]}">{piece[0]}</tspan>'
        else:
            spans += piece
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="start" style="white-space:pre" font-family="Menlo, Consolas, monospace">{spans}</text>'


def pde_to_code():
    # The update of the cell at S = 5250, v = 0.04 in maths with real numbers, above the lines of src/stencil.cpp that run it.
    # The first group shows on open, steps 0 and 1 swap in the terms and the update, step 2 adds the stability limit.
    maths_left = 60
    code_left = 60
    code_top = 250
    line_gap = 30
    out = [svg_open(1400, 560)]

    def maths(y, pieces, size=23):
        # One line of maths. Each piece is words and a colour, or words alone in ink.
        spans = ""
        for piece in pieces:
            if isinstance(piece, tuple):
                spans += f'<tspan style="fill:{piece[1]}" font-weight="700">{piece[0]}</tspan>'
            else:
                spans += piece
        return f'<text x="{maths_left}" y="{y:.1f}" font-size="{size}" text-anchor="start">{spans}</text>'

    def code_card(file_line, lines):
        # A white card with the file name and line numbers on top and the code lines under it.
        height = 50 + line_gap * len(lines)
        drawn = rect(code_left - 20, code_top - 10, 1080, height, "#FFFFFF", LINE, 2, 'rx="8"')
        drawn += text(code_left, code_top + 16, file_line, 17, "start", "dim")
        for index, pieces in enumerate(lines):
            drawn += code_text(code_left, code_top + 50 + index * line_gap, pieces)
        return drawn

    S = TERM_COLOURS["S"]
    SS = TERM_COLOURS["SS"]
    v = TERM_COLOURS["v"]
    vv = TERM_COLOURS["vv"]
    Sv = TERM_COLOURS["Sv"]
    r = TERM_COLOURS["r"]
    sub = subscript
    differences = (maths(40, ["The five differences, each a multiplication by a reciprocal worked out once per solve"], 24)
                   + maths(84, [(sub("V", "S"), S), " = 0.5,   ", (sub("V", "SS"), SS), " = 0.02,   ",
                                (sub("V", "v"), v), " = 0,   ", (sub("V", "vv"), vv), " = 0,   ", (sub("V", "Sv"), Sv), " = 0"])
                   + maths(124, ["dividing by 2ΔS = 100 becomes multiplying by inv_two_ds = 1 / 100"], 21)
                   + code_card("src/stencil.cpp, lines 107 to 111, every cell of every step", [
                       ["double ", ("dV_dS", S), " = (east - west) * inv_two_ds;"],
                       ["double ", ("d2V_dS2", SS), " = (east - 2.0 * centre + west) * inv_ds_sq;"],
                       ["double ", ("dV_dv", v), " = (north - south) * inv_two_dv;"],
                       ["double ", ("d2V_dv2", vv), " = (north - 2.0 * centre + south) * inv_dv_sq;"],
                       ["double ", ("d2V_dS_dv", Sv), " = (north_east - north_west - south_east + south_west) * inv_four_ds_dv;"]]))
    terms = (maths(40, ["The terms, each a coefficient times one difference"], 24)
             + maths(84, [("½ v S² " + sub("V", "SS"), SS), " = 0.5 × 0.04 × 5250² × 0.02 = 11025"])
             + maths(124, [("(r - q) S " + sub("V", "S"), S), " = 0.032 × 5250 × 0.5 = 84"])
             + maths(164, [("ρξvS " + sub("V", "Sv"), Sv), " = 0,   ", ("½ ξ² v " + sub("V", "vv"), vv), " = 0,   ", ("κ(θ - v) " + sub("V", "v"), v), " = 0"])
             + code_card("src/stencil.cpp, lines 114 to 118", [
                 ["double ", ("stock_diffusion_term", SS), " = stock_diffusion * stock_price_sq[stock_i] * d2V_dS2;"],
                 ["double ", ("cross_term", Sv), " = cross_diffusion * stock_price[stock_i] * d2V_dS_dv;"],
                 ["double ", ("variance_diffusion_term", vv), " = variance_diffusion * d2V_dv2;"],
                 ["double ", ("stock_drift_term", S), " = stock_drift * stock_price[stock_i] * dV_dS;"],
                 ["double ", ("variance_drift_term", v), " = variance_drift * dV_dv;"]]))
    update = (maths(40, ["One step back, the old value plus Δt times the sum of the terms"], 24)
              + maths(84, ["change = ", ("11025", SS), " + ", ("84", S), " - ", ("0.045 × 0", r), " = 11109"])
              + maths(124, ["new V = 0 + Δt × change = 4.464 × 10⁻⁶ × 11109 = ", ("0.04959", BLUE)])
              + maths(164, ["the solver's own snapshot of this cell holds 0.04959375"], 21)
              + code_card("src/stencil.cpp, lines 121 and 122, line 121 wrapped at each plus sign", [
                  ["double change = ", ("stock_diffusion_term", SS), " + ", ("cross_term", Sv)],
                  ["               + ", ("variance_diffusion_term", vv), " + ", ("stock_drift_term", S)],
                  ["               + ", ("variance_drift_term", v), " - ", ("rate * centre", r), ";"],
                  ["out_row[stock_i] = centre + dt * change;"]]))
    stability = (text(code_left, 470, "Δt has a stability limit, 5.66 × 10⁻⁶ on this 421 x 51 grid, or the values grow without bound.", 21, "start", "amber", "700")
                 + text(code_left, 502, "The full 2048 x 512 grid has a limit of 3.66 × 10⁻⁷, so it needs 980,000 steps.", 21, "start", "amber", "700"))
    out.append(f'<g {fragment(0, "fade-out")}>' + differences + "</g>")
    out.append(f'<g {fragment(0, "current-visible")}>' + terms + "</g>")
    out.append(f'<g {fragment(1)}>' + update + "</g>")
    out.append(f'<g {fragment(2)}>' + stability + "</g>")
    out.append("</svg>")
    return "\n".join(out)


def full_run(snapshot_folder):
    # The value along the v = 0.04 row of the test grid, from the payoff at expiry back to today.
    # Steps 0 to 2 add the curves 1 and 9 days, 37 days and 91 days before expiry.
    plot_left = 150
    plot_right = 1000
    plot_top = 30
    plot_bottom = 520
    low_price = 4500
    high_price = 6000
    price_x = lambda price: plot_left + (price - low_price) / (high_price - low_price) * (plot_right - plot_left)
    value_y = lambda value: plot_bottom - value / 900 * (plot_bottom - plot_top)
    # The test grid has 421 stock nodes 50 apart, and variance row 2 is v = 0.04.
    prices = [50.0 * i for i in range(421)]
    shown = [i for i in range(421) if low_price <= prices[i] <= high_price]

    def curve(values, colour, width):
        return polyline([(price_x(prices[i]), value_y(values[i])) for i in shown], colour, width)

    def snapshot(step):
        path = os.path.join(snapshot_folder, f"snap_{step:06d}.csv")
        if not os.path.exists(path):
            raise SystemExit(f"{path} is missing, make the snapshots as the comment at the top of presentation/anim_m2/scenes.py says")
        with open(path) as snapshot_file:
            rows = list(csv.reader(snapshot_file))
        return [float(value) for value in rows[2]]

    out = [svg_open(1400, 580), arrow_head("run-arrow", INK)]
    out.append(line(plot_left, plot_bottom, plot_right + 10, plot_bottom, INK, 3, 'marker-end="url(#run-arrow)"'))
    out.append(line(plot_left, plot_bottom, plot_left, plot_top - 10, INK, 3, 'marker-end="url(#run-arrow)"'))
    out.append(text(plot_right, plot_bottom + 34, "share price S", 20, "end"))
    out.append(text(plot_left - 16, plot_top + 10, "V", 22, "end"))
    for price in (4750, 5250, 5750):
        out.append(text(price_x(price), plot_bottom + 34, str(price), 18, "middle", "dim"))
    for value in (250, 500, 750):
        out.append(text(plot_left - 12, value_y(value) + 7, str(value), 18, "end", "dim"))
    steps = [(0, [(560, "#BFDBFE"), (5600, "#93C5FD")], "1 and 9 days before expiry"),
             (1, [(22400, "#3B82F6")], "37 days before expiry"),
             (2, [(56000, "#1E3A8A")], "today, 91 days before expiry")]
    for index, curves, caption in steps:
        drawn = "".join(curve(snapshot(step), colour, 4) for step, colour in curves)
        out.append(f'<g {fragment(index)}>' + drawn + "</g>")
        out.append(f'<g {fragment(index, "current-visible")}>' + text(1040, 120, caption, 22, "start", "blue", "700") + "</g>")
    today = snapshot(56000)[104]
    out.append(f'<g {fragment(2)}>' + f'<circle cx="{price_x(5200):.1f}" cy="{value_y(today):.1f}" r="10" fill="{BLUE}"/>'
               + text(1040, 160, f"V at S = 5200 is {today:.2f}", 22, "start", "blue")
               + text(1040, 192, "196.17 on the full 2048 x 512 grid", 20, "start", "dim")
               + "</g>")
    # The payoff goes on top, so the later curves do not hide it.
    out.append(curve([max(price - 5250, 0) for price in prices], RED, 4))
    out.append(text(price_x(5520), value_y(120), "expiry, the payoff", 20, "start", "red"))
    out.append(text(1040, 60, "56,000 steps on the 421 x 51 test grid", 20, "start", "dim"))
    out.append("</svg>")
    return "\n".join(out)


def grid_step():
    # Two 8 by 6 grids. Column i, variance row j counted from the bottom.
    cell_width = 56
    cell_height = 52
    top = 96
    left_x = 120
    right_x = 860
    row_y = lambda j: top + (5 - j) * cell_height
    bottom = top + 6 * cell_height
    out = [svg_open(1400, 560), arrow_head("grid-arrow", INK), arrow_head("swap-arrow", BLUE)]
    out.append(text(left_x + 4 * cell_width, 40, "current, t = T", 28, weight="700"))
    out.append(text(left_x + 4 * cell_width, 72, "every cell holds the payoff", 20, "middle", "dim"))
    out.append(text(right_x + 4 * cell_width, 40, "next, t = T - dt", 28, weight="700"))
    out.append(text(right_x + 4 * cell_width, 72, "grey edges come from the boundary rules", 20, "middle", "dim"))
    for j in range(6):
        for i in range(8):
            out.append(rect(left_x + i * cell_width, row_y(j), cell_width, cell_height, PAYOFF_COLOURS[i]))
            on_edge = i in (0, 7) or j in (0, 5)
            out.append(rect(right_x + i * cell_width, row_y(j), cell_width, cell_height, GREY if on_edge else UNKNOWN))
    # The axes of the left grid.
    out.append(line(left_x, bottom + 34, left_x + 8 * cell_width, bottom + 34, INK, 3, 'marker-end="url(#grid-arrow)"'))
    out.append(text(left_x + 4 * cell_width, bottom + 66, "share price S, 2048 nodes in the real grid", 20))
    out.append(line(left_x - 30, bottom, left_x - 30, top, INK, 3, 'marker-end="url(#grid-arrow)"'))
    out.append(text(left_x - 46, top + 3 * cell_height, "variance v, 512 rows", 20, extra=f'transform="rotate(-90 {left_x - 46} {top + 3 * cell_height})"'))
    # The arrow between the grids and the stencil legend.
    out.append(line(left_x + 8 * cell_width + 30, 150, right_x - 30, 150, INK, 3, 'marker-end="url(#grid-arrow)"'))
    out.append(text((left_x + 8 * cell_width + right_x) / 2, 136, "one step back", 22))
    names = [["V_Sv", "V_v", "V_Sv"], ["V_S", "V", "V_S"], ["V_Sv", "V_v", "V_Sv"]]
    legend_x = 637
    legend_y = 222
    legend_cell = 42
    for r in range(3):
        for c in range(3):
            out.append(rect(legend_x + c * legend_cell, legend_y + r * legend_cell, legend_cell, legend_cell, "#FFFFFF", LINE, 1.5))
            out.append(text(legend_x + c * legend_cell + legend_cell / 2, legend_y + r * legend_cell + 27, names[r][c], 15))
    out.append(rect(legend_x, legend_y, 3 * legend_cell, 3 * legend_cell, "none", RED, 4))
    out.append(text(legend_x + 1.5 * legend_cell, legend_y + 3 * legend_cell + 30, "the nine-point stencil", 19, "middle", "red"))
    out.append(text(legend_x + 1.5 * legend_cell, legend_y + 3 * legend_cell + 54, "and the derivative", 19, "middle", "red"))
    out.append(text(legend_x + 1.5 * legend_cell, legend_y + 3 * legend_cell + 78, "each cell feeds", 19, "middle", "red"))
    # Steps 0 to 2. The stencil moves along row 1 and one new cell is written each time.
    for step, i in enumerate((1, 2, 3)):
        out.append(f'<g {fragment(step, "current-visible")}>'
                   + rect(left_x + (i - 1) * cell_width, row_y(2), 3 * cell_width, 3 * cell_height, "none", RED, 6)
                   + rect(right_x + i * cell_width, row_y(1), cell_width, cell_height, "none", RED, 6) + "</g>")
        out.append(f'<g {fragment(step)}>' + rect(right_x + i * cell_width, row_y(1), cell_width, cell_height, PAYOFF_COLOURS[i]) + "</g>")
    # Step 3, the rest of row 1.
    out.append(f'<g {fragment(3)}>' + "".join(rect(right_x + i * cell_width, row_y(1), cell_width, cell_height, PAYOFF_COLOURS[i]) for i in range(4, 7)) + "</g>")
    out.append(f'<g {fragment(3, "current-visible")}>' + text(right_x + 4 * cell_width, bottom + 40, "row 1 done, 2046 cells in the real grid", 21, "middle", "blue") + "</g>")
    # Step 4, every other interior row.
    out.append(f'<g {fragment(4)}>' + "".join(rect(right_x + i * cell_width, row_y(j), cell_width, cell_height, PAYOFF_COLOURS[i]) for j in range(2, 5) for i in range(1, 7)) + "</g>")
    out.append(f'<g {fragment(4, "current-visible")}>' + text(right_x + 4 * cell_width, bottom + 40, "rows 2 to 4 the same way, in any order", 21, "middle", "blue") + "</g>")
    # Step 5, the swap.
    out.append(f'<g {fragment(5)}>'
               + f'<path d="M {right_x + 4 * cell_width} {bottom + 20} Q 850 {bottom + 130} 600 {bottom + 62}" fill="none" stroke="{BLUE}" stroke-width="4" marker-end="url(#swap-arrow)"/>'
               + text(right_x + 4 * cell_width, bottom + 106, "std::swap, next becomes current", 22, "middle", "blue", "700")
               + text(right_x + 4 * cell_width, bottom + 134, "two pointers change, no doubles are copied", 20, "middle", "dim")
               + "</g>")
    out.append("</svg>")
    return "\n".join(out)


def data_layout():
    # A 6 by 3 grid that slides down into one line of memory.
    cell_width = 64
    cell_height = 46
    left = 130
    flat_y = 220
    row_top = {2: 10, 1: 60, 0: 110}
    row_fill = {0: "#DBEAFE", 1: "#BFDBFE", 2: LIGHT_BLUE}
    out = [svg_open(1300, 340), arrow_head("stride-arrow", BLUE)]
    for j in (2, 1, 0):
        out.append(text(left - 14, row_top[j] + 30, f"j = {j}", 20, "end", "dim"))
        for i in range(6):
            out.append(rect(left + i * cell_width, row_top[j], cell_width, cell_height, "none", LINE, 2, 'stroke-dasharray="6 5"'))
    out.append(text(left + 6 * cell_width + 40, 40, "3 variance rows of 6 stock nodes here,", 22, "start"))
    out.append(text(left + 6 * cell_width + 40, 70, "512 rows of 2048 nodes in the real grid.", 22, "start"))
    # Step 0, each row slides down into one line of memory.
    for j in (0, 1, 2):
        shift_x = j * 6 * cell_width
        shift_y = flat_y - row_top[j]
        cells = "".join(rect(left + i * cell_width, row_top[j], cell_width, cell_height, row_fill[j])
                        + text(left + i * cell_width + cell_width / 2, row_top[j] + 30, f"{i},{j}", 19) for i in range(6))
        out.append(f'<g {fragment(0, "go")} style="--tx:{shift_x}px; --ty:{shift_y}px">{cells}</g>')
    out.append(f'<g {fragment(0)}>' + text(left + 18 * cell_width, flat_y - 40, "one std::vector&lt;double&gt;, row 0, then row 1, then row 2", 20, "end", "dim") + "</g>")
    # Step 1, the index of each cell.
    out.append(f'<g {fragment(1)}>' + "".join(text(left + k * cell_width + cell_width / 2, flat_y - 10, str(k), 18, "middle", "blue") for k in range(18))
               + text(left + 6 * cell_width + 40, 120, "index(i, j) = j × 6 + i", 26, "start", "blue", "700")
               + text(left + 6 * cell_width + 40, 150, "j × 2048 + i in the real grid", 20, "start", "dim") + "</g>")
    # Step 2, east is the next double and north is one whole row on.
    start = left + 2 * cell_width + cell_width / 2
    end = left + 8 * cell_width + cell_width / 2
    out.append(f'<g {fragment(2)}>'
               + f'<path d="M {start} {flat_y + cell_height + 6} Q {(start + end) / 2} {flat_y + cell_height + 70} {end} {flat_y + cell_height + 6}" fill="none" stroke="{BLUE}" stroke-width="4" marker-end="url(#stride-arrow)"/>'
               + text(end + 30, flat_y + cell_height + 46, "north of (2,0) is 6 cells on, 2048 doubles or 16 KiB in the real grid", 20, "start", "blue")
               + "</g>")
    out.append("</svg>")
    return "\n".join(out)


def simd_row():
    # One row of 15 cells, updated 4 at a time with one cell left for the scalar row.
    cell_size = 60
    left = 20
    row_y = 200
    out = [svg_open(940, 340)]
    for i in range(15):
        on_edge = i in (0, 14)
        out.append(rect(left + i * cell_size, row_y, cell_size, cell_size, GREY if on_edge else "#E2E8F0"))
        out.append(text(left + i * cell_size + cell_size / 2, row_y + cell_size + 30, str(i), 20, "middle", "dim"))
    out.append(text(left + 7.5 * cell_size, row_y + cell_size + 66, "one row, 15 cells here, 2048 in the real grid", 22, "middle", "dim"))
    # Step 0, the west, centre and east loads for cells 1 to 4.
    loads = [("west", 0, 88), ("centre", 1, 118), ("east", 2, 148)]
    bars = "".join(rect(left + first * cell_size + 4, bar_y, 4 * cell_size - 8, 22, BLUE if name == "centre" else LIGHT_BLUE, "none", 0, 'rx="4"')
                   + text(left + first * cell_size + 4 * cell_size + 8, bar_y + 18, name, 20, "start", "blue") for name, first, bar_y in loads)
    out.append(f'<g {fragment(0, "current-visible")}>{bars}'
               + text(left + 3 * cell_size, 60, "three of the nine 4-wide loads", 22, "middle", "blue") + "</g>")
    # Steps 0 to 2, one group of four cells per pass.
    for step, first in ((0, 1), (1, 5), (2, 9)):
        box = rect(left + first * cell_size - 3, row_y - 3, 4 * cell_size + 6, cell_size + 6, "none", BLUE, 6, 'rx="6"')
        out.append(f'<g {fragment(step, "current-visible")}>{box}</g>')
        out.append(f'<g {fragment(step)}>' + "".join(rect(left + i * cell_size, row_y, cell_size, cell_size, LIGHT_BLUE) for i in range(first, first + 4)) + "</g>")
    out.append(f'<g {fragment(1, "current-visible")}>' + text(left + 7.5 * cell_size, 130, "stock_i = 5, cells 5 to 8 in the next pass", 24, "middle", "blue") + "</g>")
    out.append(f'<g {fragment(2, "current-visible")}>' + text(left + 7.5 * cell_size, 130, "stock_i = 9, the loop stops at stock_i = 13", 24, "middle", "blue") + "</g>")
    # Step 3, the one cell left over goes to the scalar row.
    out.append(f'<g {fragment(3)}>' + rect(left + 13 * cell_size, row_y, cell_size, cell_size, "#FCD34D") + "</g>")
    out.append(f'<g {fragment(3, "current-visible")}>' + rect(left + 13 * cell_size - 3, row_y - 3, cell_size + 6, cell_size + 6, "none", AMBER, 6, 'rx="6"')
               + text(left + 7.5 * cell_size, 110, "cell 13 goes to the scalar row", 24, "middle", "amber")
               + text(left + 7.5 * cell_size, 142, "2046 interior cells give 511 groups and 2 left over", 22, "middle", "amber") + "</g>")
    out.append("</svg>")
    return "\n".join(out)


def omp_bands():
    # Rows 1 to 510 split by schedule(static), 6 threads of 22 rows, then 18 of 21.
    grid_x = 170
    grid_width = 380
    grid_y = 20
    grid_height = 460
    row_height = grid_height / 512
    row_top = lambda r: grid_y + grid_height - (r + 1) * row_height
    out = [svg_open(960, 520)]
    out.append(rect(grid_x, grid_y, grid_width, grid_height, "#E5E7EB", LINE, 1))
    out.append(text(grid_x - 24, grid_y + grid_height / 2, "variance rows 0 to 511", 22, "middle", "dim",
                    extra=f'transform="rotate(-90 {grid_x - 24} {grid_y + grid_height / 2})"'))
    bands = []
    first_row = 1
    for thread in range(24):
        num_rows = 22 if thread < 6 else 21
        bands.append((thread, first_row, first_row + num_rows - 1))
        first_row = first_row + num_rows
    shapes = []
    for thread, first, last in bands:
        shapes.append(rect(grid_x, row_top(last), grid_width, (last - first + 1) * row_height, LIGHT_BLUE if thread % 2 == 0 else BLUE, "none", 0))
    band_label = lambda thread, first, last: text(grid_x + grid_width + 16, (row_top(last) + row_top(first) + row_height) / 2 + 7,
                                                  f"thread {thread}, rows {first} to {last}", 21, "start")
    shapes.append(band_label(*bands[0]))
    shapes.append(text(grid_x + grid_width + 16, grid_y + grid_height - 90, "⋮  threads 1 to 22", 21, "start", "dim"))
    shapes.append(band_label(*bands[23]))
    out.append(f'<g {fragment(0)}>' + "".join(shapes) + "</g>")
    note_box = lambda colour: rect(grid_x + grid_width + 10, 150, 390, 76, "#FFFFFF", colour, 3, 'rx="8"')
    note_x = grid_x + grid_width + 205
    out.append(f'<g {fragment(1, "current-visible")}>' + note_box(AMBER)
               + text(note_x, 182, "implicit barrier at the end of the loop", 21, "middle", "amber")
               + text(note_x, 210, "all 24 threads wait for the last one", 21, "middle", "amber") + "</g>")
    edges = (rect(grid_x, row_top(0), grid_width, 8, AMBER, "none", 0) + rect(grid_x, grid_y, grid_width, 8, AMBER, "none", 0)
             + rect(grid_x, grid_y, 7, grid_height, AMBER, "none", 0) + rect(grid_x + grid_width - 7, grid_y, 7, grid_height, AMBER, "none", 0))
    out.append(f'<g {fragment(2)}>{edges}</g>')
    out.append(f'<g {fragment(2, "current-visible")}>' + note_box(AMBER)
               + text(note_x, 182, "then the main thread alone updates", 21, "middle", "amber")
               + text(note_x, 210, "row 0, the stock edges, then row 511", 21, "middle", "amber") + "</g>")
    out.append(f'<g {fragment(3)}>' + note_box(BLUE)
               + text(note_x, 182, "swap the buffers", 21, "middle", "blue")
               + text(note_x, 210, "and the same split every step", 21, "middle", "blue") + "</g>")
    out.append("</svg>")
    return "\n".join(out)


def mpi_rows():
    # Two processes, each with a whole copy of the grid, drawn as 8 rows each.
    labels = ["511", "⋮", "257", "256", "255", "254", "⋮", "0"]
    bar_height = 40
    gap = 6
    first_y = 70
    width = 260
    left_of = {0: 110, 1: 590}
    bar_y = lambda k: first_y + k * (bar_height + gap)
    owned = {0: range(4, 8), 1: range(0, 4)}
    edge = {0: 4, 1: 3}
    copied = {0: 3, 1: 4}
    row_name = lambda k: f"row {labels[k]}" if labels[k] != "⋮" else "⋮"
    out = [svg_open(960, 500), arrow_head("row-arrow", BLUE)]
    out.append(text(left_of[0] + width / 2, 30, "process 0, owns rows 0 to 255", 22, weight="700"))
    out.append(text(left_of[1] + width / 2, 30, "process 1, owns rows 256 to 511", 22, weight="700"))
    for rank in (0, 1):
        for k in range(len(labels)):
            dashed = 'stroke-dasharray="7 5"' if k == copied[rank] else ""
            out.append(rect(left_of[rank], bar_y(k), width, bar_height, "#F3F4F6", GREY if k == copied[rank] else LINE, 2, f'rx="4" {dashed}'))
            out.append(text(left_of[rank] + width / 2, bar_y(k) + 27, row_name(k), 19, "middle", "dim"))
    # Step 0, the rows each process owns.
    shapes = []
    for rank in (0, 1):
        for k in owned[rank]:
            shapes.append(rect(left_of[rank], bar_y(k), width, bar_height, PALE_BLUE, "none", 0, 'rx="4"'))
            shapes.append(text(left_of[rank] + width / 2, bar_y(k) + 27, row_name(k), 19))
    out.append(f'<g {fragment(0)}>' + "".join(shapes) + "</g>")
    # Step 1, the edge rows.
    shapes = []
    for rank in (0, 1):
        k = edge[rank]
        shapes.append(rect(left_of[rank], bar_y(k), width, bar_height, BLUE, "none", 0, 'rx="4"'))
        shapes.append(text(left_of[rank] + width / 2, bar_y(k) + 27, f"row {labels[k]}, edge row", 19, extra='style="fill:#FFFFFF"'))
    out.append(f'<g {fragment(1)}>' + "".join(shapes) + "</g>")
    # Step 2, each edge row lands in the other process's copied row.
    row_255 = bar_y(edge[0]) + bar_height / 2
    row_256 = bar_y(edge[1]) + bar_height / 2
    middle = (left_of[0] + width + left_of[1]) / 2
    shapes = [line(left_of[0] + width + 6, row_255, left_of[1] - 8, row_255, BLUE, 4, 'marker-end="url(#row-arrow)"'),
              line(left_of[1] - 6, row_256, left_of[0] + width + 8, row_256, BLUE, 4, 'marker-end="url(#row-arrow)"'),
              rect(left_of[1], bar_y(copied[1]), width, bar_height, LIGHT_BLUE, "none", 0, 'rx="4"'),
              text(left_of[1] + width / 2, bar_y(copied[1]) + 27, "row 255, copied row", 19),
              rect(left_of[0], bar_y(copied[0]), width, bar_height, LIGHT_BLUE, "none", 0, 'rx="4"'),
              text(left_of[0] + width / 2, bar_y(copied[0]) + 27, "row 256, copied row", 19),
              text(middle, bar_y(7) + 4, "16 KiB", 22, "middle", "blue", "700"),
              text(middle, bar_y(7) + 30, "each way", 20, "middle", "blue")]
    out.append(f'<g {fragment(2)}>' + "".join(shapes) + "</g>")
    out.append(f'<g {fragment(3)}>' + text(480, bar_y(8) + 34, "MPI_Waitall, then each process swaps its buffers", 22, "middle", "blue", "700") + "</g>")
    out.append("</svg>")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Charts that draw themselves when the slide opens
# ---------------------------------------------------------------------------

def ladder_chart():
    # Log scale bars, one per rung, growing left to right.
    plot_left = 120
    plot_right = 880
    plot_top = 30
    plot_bottom = 450
    low = 10
    high = 20000
    y_of = lambda value: log_position(value, low, high, plot_bottom, plot_top)
    out = [svg_open(900, 540)]
    for tick, label in ((10, "10"), (100, "100"), (1000, "1,000"), (10000, "10,000")):
        out.append(line(plot_left, y_of(tick), plot_right, y_of(tick), LINE, 1.5))
        out.append(text(plot_left - 12, y_of(tick) + 7, label, 20, "end", "dim"))
    out.append(text(36, (plot_top + plot_bottom) / 2, "Mcell-updates/s, log scale", 20, "middle", "dim",
                    extra=f'transform="rotate(-90 36 {(plot_top + plot_bottom) / 2})"'))
    slot = (plot_right - plot_left) / len(LADDER)
    bar_width = slot * 0.62
    for k, (name, detail, rate, total, gain) in enumerate(LADDER):
        centre = plot_left + slot * (k + 0.5)
        delay = 0.15 + 0.3 * k
        colour = GREY if k == 0 else BLUE
        top = y_of(rate)
        out.append(rect(centre - bar_width / 2, top, bar_width, plot_bottom - top, colour, "none", 0, f'class="bar" style="--d:{delay:.2f}s"'))
        labels = text(centre, top - 34, f"{rate:,}", 22, weight="700") + text(centre, top - 10, total, 20, "middle", "dim")
        if gain:
            labels += text(centre, top + 30, gain, 21, weight="700", extra='style="fill:#FFFFFF"')
        out.append(f'<g class="pop" style="--d:{delay + 0.3:.2f}s">{labels}</g>')
        out.append(text(centre, plot_bottom + 30, name, 21))
        if detail:
            out.append(text(centre, plot_bottom + 56, detail, 19, "middle", "dim"))
    out.append(line(plot_left, plot_bottom, plot_right, plot_bottom, INK, 2))
    out.append("</svg>")
    return "\n".join(out)


def roofline(node):
    # Log-log roofline of one core, or of the whole node when node is True.
    plot_left = 110
    plot_right = 870
    plot_top = 30
    plot_bottom = 450
    x_low = 0.1
    x_high = 10
    y_low = 1
    y_high = 2000 if node else 100
    x_of = lambda value: log_position(value, x_low, x_high, plot_left, plot_right)
    y_of = lambda value: log_position(value, y_low, y_high, plot_bottom, plot_top)
    out = [svg_open(900, 520), arrow_head("rise-arrow", BLUE)]
    ticks_y = (1, 10, 100, 1000) if node else (1, 10, 100)
    for tick in ticks_y:
        out.append(line(plot_left, y_of(tick), plot_right, y_of(tick), LINE, 1.5))
        out.append(text(plot_left - 12, y_of(tick) + 7, f"{tick:,}", 20, "end", "dim"))
    for tick in (0.1, 1, 10):
        out.append(line(x_of(tick), plot_top, x_of(tick), plot_bottom, LINE, 1.5))
        out.append(text(x_of(tick), plot_bottom + 30, f"{tick:g}", 20, "middle", "dim"))
    out.append(text((plot_left + plot_right) / 2, plot_bottom + 62, "arithmetic intensity, flops per byte", 20, "middle", "dim"))
    out.append(text(34, (plot_top + plot_bottom) / 2, "GFLOP/s", 20, "middle", "dim", extra=f'transform="rotate(-90 34 {(plot_top + plot_bottom) / 2})"'))
    # The roofs, drawn one after another when the slide opens.
    bandwidth = TRIAD_NODE if node else TRIAD_CORE
    vector_peak = AVX2_PEAK_NODE if node else AVX2_PEAK_CORE
    ridge = min(vector_peak / bandwidth, x_high)
    memory_end = bandwidth * ridge
    roofs = [(x_low, bandwidth * x_low, ridge, memory_end, INK, f"memory, {bandwidth} GB/s", 0.1)]
    if ridge < x_high:
        roofs.append((ridge, vector_peak, x_high, vector_peak, INK, "", 0.35))
    else:
        roofs.append((x_low * 3, vector_peak, x_high, vector_peak, INK, "", 0.35))
    for x1, y1, x2, y2, colour, label, delay in roofs:
        out.append(f'<line x1="{x_of(x1):.1f}" y1="{y_of(y1):.1f}" x2="{x_of(x2):.1f}" y2="{y_of(y2):.1f}" stroke="{colour}" stroke-width="4" pathLength="1" class="draw" style="--d:{delay}s"/>')
    out.append(f'<g class="pop" style="--d:0.5s">' + text(plot_right - 6, y_of(vector_peak) - 12, f"AVX2 FMA peak, {vector_peak:g}", 20, "end") + "</g>")
    out.append(f'<g class="pop" style="--d:0.3s">' + text(x_of(x_low) + 12, y_of(bandwidth * x_low) - 30, f"memory roof, {bandwidth} GB/s", 20, "start",
                                                          extra=f'transform="rotate({-math.degrees(math.atan2(y_of(1) - y_of(10), x_of(10) - x_of(1))):.1f} {x_of(x_low) + 12:.1f} {y_of(bandwidth * x_low) - 30:.1f})"') + "</g>")
    if not node:
        scalar_start = SCALAR_PEAK / bandwidth
        out.append(f'<line x1="{x_of(scalar_start):.1f}" y1="{y_of(SCALAR_PEAK):.1f}" x2="{x_of(x_high):.1f}" y2="{y_of(SCALAR_PEAK):.1f}" stroke="{SLATE}" stroke-width="3" stroke-dasharray="10 7" class="pop" style="--d:0.6s"/>')
        out.append(f'<g class="pop" style="--d:0.6s">' + text(plot_right - 6, y_of(SCALAR_PEAK) + 26, f"scalar peak, {SCALAR_PEAK}", 20, "end", "dim") + "</g>")
    kernel_x = x_of(FLOPS_PER_BYTE)
    out.append(f'<line x1="{kernel_x:.1f}" y1="{plot_top}" x2="{kernel_x:.1f}" y2="{plot_bottom}" stroke="{GREY}" stroke-width="2.5" stroke-dasharray="6 6"/>')
    out.append(text(kernel_x + 8, plot_bottom - 12, "31 / 24 = 1.29", 19, "start", "dim"))
    dot = lambda y, fill, stroke=BLUE: f'<circle cx="{kernel_x:.1f}" cy="{y:.1f}" r="11" fill="{fill}" stroke="{stroke}" stroke-width="4"/>'
    if not node:
        # Step 0, serial on the scalar peak. Step 1, simd rises above it. Step 2, the big grid just below simd.
        serial_y = y_of(SERIAL_GFLOPS)
        simd_y = y_of(SIMD_GFLOPS)
        big_y = y_of(SIMD_GFLOPS_BIG)
        out.append(f'<g {fragment(0)}>' + dot(serial_y, "#FFFFFF", SLATE) + text(kernel_x - 20, serial_y + 7, f"serial {SERIAL_GFLOPS}", 21, "end") + "</g>")
        out.append(f'<g {fragment(1)}><g {fragment(1, "go")} style="--tx:0px; --ty:{simd_y - serial_y:.1f}px">' + dot(serial_y, BLUE) + "</g>"
                   + text(kernel_x - 20, simd_y - 8, f"simd {SIMD_GFLOPS}", 21, "end", "blue", "700")
                   + text(kernel_x - 20, (serial_y + simd_y) / 2 + 7, "2.9x", 21, "end", "blue") + "</g>")
        out.append(f'<g {fragment(2)}><circle cx="{kernel_x:.1f}" cy="{big_y:.1f}" r="7" fill="#FFFFFF" stroke="{BLUE}" stroke-width="3"/>'
                   + text(kernel_x + 20, big_y + 12, f"simd, 4096 x 1024, {SIMD_GFLOPS_BIG}", 20, "start", "blue") + "</g>")
    else:
        # Step 0, one core of simd for reference. Step 1, 24 threads rise to 195. Step 2, the big grid sits on the memory roof.
        simd_y = y_of(SIMD_GFLOPS)
        omp_y = y_of(OMP_GFLOPS)
        big_y = y_of(OMP_GFLOPS_BIG)
        out.append(f'<g {fragment(0)}>' + dot(simd_y, "#FFFFFF", SLATE) + text(kernel_x - 20, simd_y + 7, f"simd, one core, {SIMD_GFLOPS}", 21, "end") + "</g>")
        out.append(f'<g {fragment(1)}><g {fragment(1, "go")} style="--tx:0px; --ty:{omp_y - simd_y:.1f}px">' + dot(simd_y, BLUE) + "</g>"
                   + text(kernel_x - 20, omp_y + 7, f"omp, 2048 x 512, {OMP_GFLOPS}", 21, "end", "blue", "700") + "</g>")
        out.append(f'<g {fragment(2)}>' + dot(big_y, "#FFFFFF") + text(kernel_x + 22, big_y + 30, f"omp, 4096 x 1024, {OMP_GFLOPS_BIG}, on the roof", 20, "start", "blue") + "</g>")
    out.append(line(plot_left, plot_bottom, plot_right, plot_bottom, INK, 2))
    out.append("</svg>")
    return "\n".join(out)


def impact(rows):
    # Bars for the rung before and after a code change. Each row is name, rate, label, grows.
    # A row that grows is a fragment, so it fills in after the code has changed.
    label_width = 220
    bar_left = label_width + 10
    bar_space = 1000 - bar_left - 190
    biggest = max(rate for _, rate, _, _ in rows)
    row_height = 44
    out = [svg_open(1000, row_height * len(rows) + 8)]
    for k, (name, rate, label, grows) in enumerate(rows):
        y = 4 + k * row_height
        width = bar_space * rate / biggest
        out.append(text(label_width, y + 28, name, 21, "end"))
        bar = rect(bar_left, y + 4, width, row_height - 12, BLUE if grows else GREY, "none", 0, 'rx="4"')
        caption = text(bar_left + width + 12, y + 28, label, 21, "start", "blue" if grows else "dim", "700" if grows else "400")
        if grows:
            out.append(f'<g {fragment(None)}><g class="grow-x">{bar}</g>{caption}</g>')
        else:
            out.append(bar + caption)
    out.append("</svg>")
    return "\n".join(out)


# The snapshots of the 421 by 51 test grid, every 560th step.
SNAPSHOT_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "presentation", "anim_m2", "snapshots", "dump")


BUILDERS = {
    "option": option_paths,
    "backwards": lambda: backwards_in_time(3),
    "differences": finite_differences,
    "pde-code": pde_to_code,
    "grid": grid_step,
    "full-run": lambda: full_run(SNAPSHOT_FOLDER),
    "layout": data_layout,
    "simd": simd_row,
    "omp": omp_bands,
    "mpi": mpi_rows,
    "ladder": ladder_chart,
    "roof-core": lambda: roofline(False),
    "roof-node": lambda: roofline(True),
    "impact-serial": lambda: impact([("baseline", 57, "57", False), ("serial", 159, "159, 2.8x", True)]),
    "impact-simd": lambda: impact([("serial", 159, "159", False), ("simd", 456, "456, 2.9x", True)]),
    "impact-omp": lambda: impact([("simd", 456, "456", False), ("omp, 24 threads", 6367, "6367, 14x", True)]),
    "impact-mpi": lambda: impact([("omp, 24 threads", 6367, "6367", False), ("mpi, one node", 6984, "6984, 1.10x", True)]),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--deck", default="presentation/deck_m2.qmd")
    args = parser.parse_args()
    with open(args.deck) as deck_file:
        deck = deck_file.read()
    for name, build in BUILDERS.items():
        # Everything between the two marker comments is replaced.
        pattern = re.compile(rf"(<!-- svg:{name} -->\n).*?(\n<!-- /svg:{name} -->)", re.S)
        if not pattern.search(deck):
            raise SystemExit(f"no <!-- svg:{name} --> marker in {args.deck}")
        deck = pattern.sub(lambda match: match.group(1) + build() + match.group(2), deck)
    with open(args.deck, "w") as deck_file:
        deck_file.write(deck)
    print(f"rewrote {len(BUILDERS)} diagrams in {args.deck}")


if __name__ == "__main__":
    main()
