# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Four manim animations for the recorded Milestone 2 video, 1600 by 900 at 30 frames a second.
# Render commands, with MEDIA as any scratch folder and LaTeX on the PATH:
#   export PATH=/Library/TeX/texbin:$PATH
#   uv run --python 3.12 --with manim manim -qh -r 1600,900 --fps 30 --save_sections --media_dir MEDIA scenes.py ModelToPde
#   uv run --python 3.12 --with manim manim -qh -r 1600,900 --fps 30 --save_sections --media_dir MEDIA scenes.py PdeToCode
#   uv run --python 3.12 --with manim manim -qh -r 1600,900 --fps 30 --save_sections --media_dir MEDIA scenes.py Differences
#   uv run --python 3.12 --with manim manim -qh -r 1600,900 --fps 30 --save_sections --media_dir MEDIA scenes.py FullRun
# FullRun reads the solver snapshots of the test grid from snapshots/ next to this file, or from the
# folder in HESTON_SNAPSHOTS. snapshots/early holds steps 1 to 60 and snapshots/dump every 560th step,
# one CSV of 51 variance rows by 421 stock columns each. From the project folder they are made by:
#   ./heston --config presentation/anim_m2/full_run.cfg --solver serial --dump-every 560 --dump-dir presentation/anim_m2/snapshots/dump
#   ./heston --config presentation/anim_m2/full_run.cfg --solver serial --nt 60 --maturity 0.00026785714285714287 --dump-every 1 --dump-dir presentation/anim_m2/snapshots/early
# The second run keeps the same dt, 0.25 / 56000, so its 60 steps are the first 60 steps of the full run.
# Then copy each video here and write its stop times, for example:
#   cp MEDIA/videos/scenes/900p30/ModelToPde.mp4 model_pde.mp4
#   uv run --python 3.12 --with manim python scenes.py MEDIA/videos/scenes/900p30/sections/ModelToPde.json model_pde.json
# The other three are PdeToCode to pde_code, Differences to differences and FullRun to full_run.

import json
import os
import sys

import numpy as np
from manim import *

config.background_color = "#F0F3F9"

# Colours from the deck palette.
INK = "#1F2937"
GREY_TEXT = "#6B7280"
BLUE_ACCENT = "#2563EB"
RED_ACCENT = "#DC2626"
AMBER_ACCENT = "#D97706"
GREEN_ACCENT = "#059669"
CARD_FILL = "#FFFFFF"
CARD_BORDER = "#D1D5DB"
HIGHLIGHT_FILL = "#F6C9C9"

# One colour for each of the six PDE terms, the same in every video.
STOCK_DRIFT_COLOUR = BLUE_ACCENT
STOCK_DIFFUSION_COLOUR = AMBER_ACCENT
VARIANCE_DRIFT_COLOUR = GREEN_ACCENT
VARIANCE_DIFFUSION_COLOUR = "#7C3AED"
CROSS_COLOUR = "#DB2777"
DISCOUNT_COLOUR = "#0891B2"

# LaTeX text in Helvetica to match the deck. The maths stays in Computer Modern.
SANS_TEMPLATE = TexTemplate()
SANS_TEMPLATE.add_to_preamble(r"\usepackage[scaled]{helvet}")
SANS_TEMPLATE.add_to_preamble(r"\renewcommand{\familydefault}{\sfdefault}")

MONO_FONT = "Menlo"
CAPTION_SIZE = 32

# Still time at the end of every section, where the deck pauses the video.
SECTION_END_PAUSE = 0.4

# The webcam covers the right quarter of the width and the bottom quarter of the height.
# Nothing may be drawn there.
WEBCAM_LEFT = config.frame_width / 4
WEBCAM_TOP = -config.frame_height / 4

# The six PDE terms in the order they appear in the PDE.
# Each entry is the name, the LaTeX and the colour.
PDE_TERMS = [
    ("stock_diffusion", r"\tfrac{1}{2} v S^2 V_{SS}", STOCK_DIFFUSION_COLOUR),
    ("cross", r"\rho \xi v S V_{Sv}", CROSS_COLOUR),
    ("variance_diffusion", r"\tfrac{1}{2} \xi^2 v V_{vv}", VARIANCE_DIFFUSION_COLOUR),
    ("stock_drift", r"(r - q) S V_S", STOCK_DRIFT_COLOUR),
    ("variance_drift", r"\kappa(\theta - v) V_v", VARIANCE_DRIFT_COLOUR),
    ("discount", r"r V", DISCOUNT_COLOUR),
]

# The five central differences, in the order of the PDE terms they feed.
# Each entry is the derivative, its colour, the numerator, the denominator, the cell weights and the reciprocal name.
# A cell is (column offset, row offset) from the centre cell, with rows counted upwards.
CENTRAL_DIFFERENCES = [
    ("V_{SS}", STOCK_DIFFUSION_COLOUR, r"V_{i+1,j} - 2V_{i,j} + V_{i-1,j}", r"\Delta S^2", {(1, 0): "+1", (0, 0): "-2", (-1, 0): "+1"}, "inv_ds_sq"),
    ("V_{Sv}", CROSS_COLOUR, r"V_{i+1,j+1} - V_{i-1,j+1} - V_{i+1,j-1} + V_{i-1,j-1}", r"4\Delta S \Delta v", {(1, 1): "+1", (-1, 1): "-1", (1, -1): "-1", (-1, -1): "+1"}, "inv_four_ds_dv"),
    ("V_{vv}", VARIANCE_DIFFUSION_COLOUR, r"V_{i,j+1} - 2V_{i,j} + V_{i,j-1}", r"\Delta v^2", {(0, 1): "+1", (0, 0): "-2", (0, -1): "+1"}, "inv_dv_sq"),
    ("V_S", STOCK_DRIFT_COLOUR, r"V_{i+1,j} - V_{i-1,j}", r"2\Delta S", {(1, 0): "+1", (-1, 0): "-1"}, "inv_two_ds"),
    ("V_v", VARIANCE_DRIFT_COLOUR, r"V_{i,j+1} - V_{i,j-1}", r"2\Delta v", {(0, 1): "+1", (0, -1): "-1"}, "inv_two_dv"),
]

# The reciprocals in the order of src/stencil.cpp lines 28 to 32.
# Each entry is the name and its denominator.
RECIPROCALS = [
    ("inv_two_ds", r"2\Delta S"),
    ("inv_ds_sq", r"\Delta S^2"),
    ("inv_two_dv", r"2\Delta v"),
    ("inv_dv_sq", r"\Delta v^2"),
    ("inv_four_ds_dv", r"4\Delta S \Delta v"),
]

# Lines 121 and 122 of src/stencil.cpp, with line 121 wrapped at each plus sign.
# The tokens are exactly the ones in the file.
CODE_LINES = [
    "double change = stock_diffusion_term",
    "              + cross_term",
    "              + variance_diffusion_term",
    "              + stock_drift_term",
    "              + variance_drift_term",
    "              - rate * centre;",
    "out_row[stock_i] = centre + dt * change;",
]

# The colour of each term name in the code, matching the maths.
CODE_COLOURS = {
    "stock_diffusion_term": STOCK_DIFFUSION_COLOUR,
    "cross_term": CROSS_COLOUR,
    "variance_diffusion_term": VARIANCE_DIFFUSION_COLOUR,
    "stock_drift_term": STOCK_DRIFT_COLOUR,
    "variance_drift_term": VARIANCE_DRIFT_COLOUR,
    "rate * centre": DISCOUNT_COLOUR,
}

# Names and values of the nine cells in the Differences video, top row first.
# Every row holds the payoff max(S - K, 0) at S = 5200, 5250 and 5300 with K = 5250.
CELL_NAMES = [
    ["north_west", "north", "north_east"],
    ["west", "centre", "east"],
    ["south_west", "south", "south_east"],
]
CELL_VALUES = [
    ["0", "0", "50"],
    ["0", "0", "50"],
    ["0", "0", "50"],
]

# The five differences with real numbers.
# Each entry is the derivative, its colour, the numerator pieces, the denominator and the result.
# A numerator piece named after a cell gets its number from that cell.
NUMBER_DIFFERENCES = [
    ("V_S", STOCK_DRIFT_COLOUR, [("east", "50"), ("minus", "-"), ("west", "0")], r"2 \times 50", "0.5"),
    ("V_{SS}", STOCK_DIFFUSION_COLOUR, [("east", "50"), ("minus", "-"), ("two", r"2 \times"), ("centre", "0"), ("plus", "+"), ("west", "0")], r"50^2", "0.02"),
    ("V_v", VARIANCE_DRIFT_COLOUR, [("north", "0"), ("minus", "-"), ("south", "0")], r"2 \times 0.02", "0"),
    ("V_{vv}", VARIANCE_DIFFUSION_COLOUR, [("north", "0"), ("minus", "-"), ("two", r"2 \times"), ("centre", "0"), ("plus", "+"), ("south", "0")], r"0.02^2", "0"),
    ("V_{Sv}", CROSS_COLOUR, [("north_east", "50"), ("minus", "-"), ("north_west", "0"), ("minus_two", "-"), ("south_east", "50"), ("plus", "+"), ("south_west", "0")], r"4 \times 50 \times 0.02", "0"),
]

# Solver snapshots for the FullRun video.
SNAPSHOT_FOLDER = os.environ.get("HESTON_SNAPSHOTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "snapshots"))

# The part of the test grid the heatmap shows, S from 4000 to 6600 and v from 0 to 0.2.
FIRST_COLUMN = 80
LAST_COLUMN = 132
TOP_ROW = 10
STOCK_SPACING = 50.0
STRIKE = 5250.0

# Colour scale of the heatmap, white at no time value and dark blue at 450.
# The colour follows the square root of the time value, because time value grows like the square root of tau.
MAX_TIME_VALUE = 450.0
COLOUR_STOP_POSITIONS = np.array([0.0, 0.08, 0.4, 1.0])
COLOUR_STOPS = np.array([[255, 255, 255], [191, 219, 254], [59, 130, 246], [30, 58, 138]], dtype=float)


# A bold title in the top left corner of the frame.
# Every section opens by setting one.
def make_title(text):
    title = Tex(r"\textbf{" + text + "}", tex_template=SANS_TEMPLATE, font_size=42, color=INK)
    title.to_corner(UL, buff=0.45)
    return title


# A grey caption in the bottom left corner of the frame.
# A double backslash in the text starts a second line.
def make_caption(text, font_size=CAPTION_SIZE):
    caption = Tex(text, tex_template=SANS_TEMPLATE, font_size=font_size, color=GREY_TEXT, tex_environment="flushleft")
    caption.to_corner(DL, buff=0.45)
    return caption


# One line of maths built from named pieces.
# Each piece is (name, LaTeX) or (name, LaTeX, colour) and can move or morph on its own.
def make_equation(pieces, font_size):
    tex_strings = []
    for piece in pieces:
        tex_strings.append(piece[1])
    line = MathTex(*tex_strings, font_size=font_size, color=INK)
    parts = {}
    for index in range(len(pieces)):
        piece = pieces[index]
        if len(piece) == 3:
            line[index].set_color(piece[2])
        parts[piece[0]] = line[index]
    return line, parts


# Animations that turn one set of named pieces into another.
# Pieces with the same name move and cross fade into each other, the others fade out or fade in.
def morph_pieces(old_parts, new_parts):
    animations = []
    for name, new_piece in new_parts.items():
        if name in old_parts:
            animations.append(FadeTransform(old_parts[name], new_piece))
        else:
            animations.append(FadeIn(new_piece))
    for name, old_piece in old_parts.items():
        if name not in new_parts:
            animations.append(FadeOut(old_piece))
    return animations


# The six PDE terms as named pieces on two lines, each term after its sign.
# The prefix goes in front of line one and the suffix after line two.
def pde_pieces(prefix, suffix, sign_before_first_term):
    line_one = list(prefix)
    line_two = []
    for index in range(len(PDE_TERMS)):
        name, tex, colour = PDE_TERMS[index]
        if index < 3:
            line = line_one
        else:
            line = line_two
        if name == "discount":
            sign = "-"
        else:
            sign = "+"
        # An empty group in front keeps the spacing of a binary plus at the start of line two.
        if index == 3:
            sign = "{}" + sign
        if index > 0 or sign_before_first_term:
            line.append(("sign_" + name, sign))
        line.append((name, tex, colour))
    line_two.extend(suffix)
    return line_one, line_two


# Two centred lines of maths, one above the other.
# Returns the group and one dictionary with the named pieces of both lines.
def make_two_lines(line_one_pieces, line_two_pieces, font_size, centre_y, line_gap):
    line_one, parts_one = make_equation(line_one_pieces, font_size)
    line_two, parts_two = make_equation(line_two_pieces, font_size)
    line_one.move_to(UP * (centre_y + line_gap / 2))
    line_two.move_to(UP * (centre_y - line_gap / 2))
    parts = {}
    parts.update(parts_one)
    parts.update(parts_two)
    return VGroup(line_one, line_two), parts


# A formula with a left side, a fraction and an optional right side.
# The fraction bar sits level with the relation sign, as LaTeX would set it.
def make_fraction_formula(left_pieces, numerator_pieces, denominator_pieces, right_pieces, font_size):
    left, left_parts = make_equation(left_pieces, font_size)
    numerator, numerator_parts = make_equation(numerator_pieces, font_size)
    denominator, denominator_parts = make_equation(denominator_pieces, font_size)

    # Bar a little wider than the wider of the two lines.
    bar_width = max(numerator.width, denominator.width) + 0.16
    bar = Line(LEFT * bar_width / 2, RIGHT * bar_width / 2, stroke_width=2.5, color=INK)
    numerator.next_to(bar, UP, buff=0.1)
    denominator.next_to(bar, DOWN, buff=0.1)
    fraction = VGroup(numerator, bar, denominator)

    # Fraction to the right of the left side, bar level with the relation sign.
    fraction.next_to(left, RIGHT, buff=0.2)
    fraction.shift(UP * (left[-1].get_y() - bar.get_y()))
    group = VGroup(left, fraction)

    # Optional right side, its first sign level with the bar.
    right = None
    right_parts = {}
    if right_pieces:
        right, right_parts = make_equation(right_pieces, font_size)
        right.next_to(fraction, RIGHT, buff=0.2)
        right.shift(UP * (bar.get_y() - right[0].get_y()))
        group.add(right)

    return {
        "group": group,
        "left": left_parts,
        "numerator": numerator,
        "numerator_parts": numerator_parts,
        "bar": bar,
        "denominator": denominator,
        "right": right,
        "right_parts": right_parts,
    }


# Moves a formula so its fraction bar is at height y and its left edge at x.
# Formulas in a column then line up on their bars.
def place_formula(formula, left_x, bar_y):
    group = formula["group"]
    group.shift(RIGHT * (left_x - group.get_left()[0]))
    group.shift(UP * (bar_y - formula["bar"].get_y()))


# A 3 by 3 block of cells, top row first.
# Index row * 3 + column picks a cell.
def make_patch(cell_width, cell_height):
    cells = VGroup()
    for row in range(3):
        for column in range(3):
            cell = Rectangle(width=cell_width, height=cell_height, stroke_color=INK, stroke_width=3, fill_color=CARD_FILL, fill_opacity=1)
            cells.add(cell)
    cells.arrange_in_grid(rows=3, cols=3, buff=0)
    return cells


# The cell at a column and row offset from the centre cell.
# Rows are counted upwards, so (1, 1) is the top right cell.
def cell_at(cells, column_offset, row_offset):
    row = 1 - row_offset
    column = 1 + column_offset
    return cells[row * 3 + column]


# A white card with a grey border behind a group of mobjects.
# It matches the code blocks of the deck.
def make_card(content, padding):
    card = RoundedRectangle(corner_radius=0.15, width=content.width + 2 * padding, height=content.height + 2 * padding, fill_color=CARD_FILL, fill_opacity=1, stroke_color=CARD_BORDER, stroke_width=2)
    card.move_to(content)
    # Drawn behind everything, also behind mobjects added before it.
    card.set_z_index(-1)
    return card


# Turns a grid of time values into an RGBA image, highest variance row at the top.
# The colour is the square root of the time value over 450, read off the colour stops.
def time_value_colours(time_value):
    scaled = np.sqrt(np.clip(time_value, 0.0, MAX_TIME_VALUE) / MAX_TIME_VALUE)
    channels = []
    for channel in range(3):
        channels.append(np.interp(scaled, COLOUR_STOP_POSITIONS, COLOUR_STOPS[:, channel]))
    channels.append(np.full_like(scaled, 255.0))
    image = np.dstack(channels).astype(np.uint8)
    return image[::-1]


# Shared setup for the four videos.
# Holds the title and caption and checks the webcam corner at the end of every section.
class DeckScene(Scene):
    def setup(self):
        self.title = None
        self.caption = None

    # Swaps the title for a new one.
    # Returns the animations so they can run with others.
    def swap_title(self, text):
        animations = []
        if self.title is not None:
            animations.append(FadeOut(self.title))
        self.title = make_title(text)
        animations.append(FadeIn(self.title))
        return animations

    # Swaps the caption for a new one.
    # Returns the animations so they can run with others.
    def swap_caption(self, text, font_size=CAPTION_SIZE):
        animations = []
        if self.caption is not None:
            animations.append(FadeOut(self.caption))
        self.caption = make_caption(text, font_size)
        animations.append(FadeIn(self.caption))
        return animations

    # Holds the last frame of a section for the deck to pause on.
    # First checks that nothing on screen reaches into the webcam corner.
    def end_section(self):
        for mobject in self.mobjects:
            for part in mobject.family_members_with_points():
                if part.get_right()[0] > WEBCAM_LEFT and part.get_bottom()[1] < WEBCAM_TOP:
                    raise ValueError("Something is drawn in the webcam corner at " + str(part.get_corner(DR)))
        self.wait(SECTION_END_PAUSE)


# From the Heston model to the pricing PDE.
# Four sections, the model, the expectation, Ito's lemma and the PDE.
class ModelToPde(DeckScene):
    def construct(self):
        self.model_section()
        self.expectation_section()
        self.ito_section()
        self.pde_section()

    # Section 1 shows the two random processes with the parameter values beside them.
    # The three equations line up on their equals signs.
    def model_section(self):
        self.next_section("The Heston model")

        # The share price, the variance and the correlation of the two random drivers.
        stock_equation, self.stock_parts = make_equation([
            ("lhs", "dS"), ("eq", "="),
            ("drift", r"(r - q) S", STOCK_DRIFT_COLOUR), ("dt", r"\, dt"), ("plus", "+"),
            ("vol", r"\sqrt{v}\, S", STOCK_DIFFUSION_COLOUR), ("dW", r"\, dW_S"),
        ], 44)
        variance_equation, self.variance_parts = make_equation([
            ("lhs", "dv"), ("eq", "="),
            ("drift", r"\kappa(\theta - v)", VARIANCE_DRIFT_COLOUR), ("dt", r"\, dt"), ("plus", "+"),
            ("vol", r"\xi \sqrt{v}", VARIANCE_DIFFUSION_COLOUR), ("dW", r"\, dW_v"),
        ], 44)
        correlation_equation, self.correlation_parts = make_equation([
            ("lhs", r"dW_S\, dW_v"), ("eq", "="), ("rho", r"\rho", CROSS_COLOUR), ("dt", r"\, dt"),
        ], 44)
        self.model = VGroup(stock_equation, variance_equation, correlation_equation)

        # Equals signs on one vertical line, the parameter values to the right of each equation.
        self.values = VGroup(
            MathTex(r"r = 0.045, \quad q = 0.013", font_size=36, color=INK),
            MathTex(r"\kappa = 1.5, \quad \theta = 0.04, \quad \xi = 0.35", font_size=36, color=INK),
            MathTex(r"\rho = -0.7", font_size=36, color=INK),
        )
        rows = [(stock_equation, self.stock_parts, 1.9), (variance_equation, self.variance_parts, 0.8), (correlation_equation, self.correlation_parts, -0.3)]
        for index in range(3):
            equation, parts, y = rows[index]
            equation.shift(UP * (y - parts["eq"].get_y()))
            equation.shift(RIGHT * (-4.3 - parts["eq"].get_x()))
            self.values[index].set_y(y)
            self.values[index].align_to(RIGHT * 1.2, LEFT)

        # Spot and variance today and the contract.
        self.contract = MathTex(r"S_0 = 5200, \quad v_0 = 0.04, \quad K = 5250, \quad T = 0.25 \text{ years}", font_size=40, color=INK)
        self.contract.move_to(DOWN * 1.45)
        self.contract.align_to(LEFT * 6.4, LEFT)

        self.play(*self.swap_title("The Heston model"), LaggedStart(*[Write(equation) for equation in self.model], lag_ratio=0.3), run_time=1.6)
        self.play(FadeIn(self.values, lag_ratio=0.3), run_time=0.8)
        self.play(FadeIn(self.contract), *self.swap_caption("The share price and its variance both move randomly."), run_time=0.6)
        self.wait(0.8)
        self.end_section()

    # Section 2 writes the price as an expected discounted payoff.
    # The model shrinks to the top of the frame and stays there.
    def expectation_section(self):
        self.next_section("The price as an expectation")
        stock_equation, variance_equation, correlation_equation = self.model
        self.play(
            FadeOut(self.values), FadeOut(self.contract),
            stock_equation.animate.scale(0.78).move_to(LEFT * 3.4 + UP * 2.45),
            variance_equation.animate.scale(0.78).move_to(RIGHT * 3.4 + UP * 2.45),
            correlation_equation.animate.scale(0.78).move_to(UP * 1.65),
            *self.swap_title("The price as an expectation"),
            run_time=0.8,
        )

        # The option value as an expected discounted payoff.
        self.expectation, self.expectation_parts = make_equation([
            ("value", r"V(S, v, t) ="), ("open", r"\mathbb{E}\big[\,"),
            ("discount", r"e^{-r(T-t)}", DISCOUNT_COLOUR), ("payoff", r"\max(S_T - K,\, 0)"),
            ("given", r"\;\big|\; S_t = S,\ v_t = v \,\big]"),
        ], 44)
        self.expectation.move_to(DOWN * 0.1)
        self.play(Write(self.expectation), *self.swap_caption(r"The price is the expected discounted payoff.\\The Feynman-Kac theorem turns this expectation into a PDE."), run_time=1.0)
        self.wait(1.0)
        self.end_section()

    # Section 3 expands dV with Ito's lemma.
    # The three products of random terms become dt terms.
    def ito_section(self):
        self.next_section("Ito's lemma")
        self.play(FadeOut(self.expectation), *self.swap_title(r"It\^{o}'s lemma"), run_time=0.5)

        # dV to second order on two lines.
        line_one, self.ito_one = make_equation([
            ("dV", "dV"), ("eq", "="), ("time", r"V_t\, dt"), ("plus_one", "+"),
            ("stock", r"V_S\, dS"), ("plus_two", "+"), ("variance", r"V_v\, dv"),
        ], 44)
        line_two, ito_two = make_equation([
            ("sign_stock_diffusion", "{}+"), ("stock_diffusion", r"\tfrac{1}{2} V_{SS}\, (dS)^2"),
            ("sign_cross", "+"), ("cross", r"V_{Sv}\, dS\, dv"),
            ("sign_variance_diffusion", "+"), ("variance_diffusion", r"\tfrac{1}{2} V_{vv}\, (dv)^2"),
        ], 44)
        line_one.move_to(UP * 0.55)
        line_two.move_to(DOWN * 0.45)
        self.play(Write(line_one), Write(line_two), *self.swap_caption(r"It\^{o}'s lemma expands $dV$ to second order."), run_time=1.0)

        # The three products of random terms, each in its term colour.
        self.rules = VGroup(
            MathTex(r"(dS)^2 = v S^2\, dt", font_size=38, color=STOCK_DIFFUSION_COLOUR),
            MathTex(r"(dv)^2 = \xi^2 v\, dt", font_size=38, color=VARIANCE_DIFFUSION_COLOUR),
            MathTex(r"dS\, dv = \rho \xi v S\, dt", font_size=38, color=CROSS_COLOUR),
        )
        self.rules.arrange(RIGHT, buff=0.9)
        self.rules.move_to(DOWN * 1.55)
        self.play(FadeIn(self.rules, lag_ratio=0.3), run_time=0.8)
        self.wait(0.5)

        # Each product becomes its dt term.
        new_line_two, self.ito_two = make_equation([
            ("sign_stock_diffusion", "{}+"), ("stock_diffusion", r"\tfrac{1}{2} v S^2 V_{SS}\, dt", STOCK_DIFFUSION_COLOUR),
            ("sign_cross", "+"), ("cross", r"\rho \xi v S V_{Sv}\, dt", CROSS_COLOUR),
            ("sign_variance_diffusion", "+"), ("variance_diffusion", r"\tfrac{1}{2} \xi^2 v V_{vv}\, dt", VARIANCE_DIFFUSION_COLOUR),
        ], 44)
        new_line_two.move_to(DOWN * 0.45)
        self.play(*morph_pieces(ito_two, self.ito_two), *self.swap_caption("Each product of random terms becomes a $dt$ term."), run_time=0.9)
        self.ito_lines = VGroup(line_one, new_line_two)
        self.wait(0.6)
        self.end_section()

    # Section 4 collects the dt terms, sets them equal to rV and moves rV across.
    # The result is the Heston PDE in the six term colours.
    def pde_section(self):
        self.next_section("The Heston PDE")

        # The dt terms of dV, equal to rV.
        drift_line, drift_parts = make_equation([
            ("time_derivative", "V_t"),
            ("sign_stock_drift", "+"), ("stock_drift", r"(r - q) S V_S", STOCK_DRIFT_COLOUR),
            ("sign_variance_drift", "+"), ("variance_drift", r"\kappa(\theta - v) V_v", VARIANCE_DRIFT_COLOUR),
            ("sign_stock_diffusion", "+"), ("stock_diffusion", r"\tfrac{1}{2} v S^2 V_{SS}", STOCK_DIFFUSION_COLOUR),
            ("sign_cross", "+"), ("cross", r"\rho \xi v S V_{Sv}", CROSS_COLOUR),
            ("sign_variance_diffusion", "+"), ("variance_diffusion", r"\tfrac{1}{2} \xi^2 v V_{vv}", VARIANCE_DIFFUSION_COLOUR),
            ("eq", "="), ("discount", "r V", DISCOUNT_COLOUR),
        ], 40)
        drift_line.move_to(DOWN * 1.55)

        # Each dt term flies down from dV, the two drift terms with their drift from the model.
        sources = {
            "time_derivative": [self.ito_one["time"]],
            "stock_drift": [self.stock_parts["drift"], self.ito_one["stock"]],
            "variance_drift": [self.variance_parts["drift"], self.ito_one["variance"]],
            "stock_diffusion": [self.ito_two["stock_diffusion"]],
            "cross": [self.ito_two["cross"]],
            "variance_diffusion": [self.ito_two["variance_diffusion"]],
        }
        animations = [FadeOut(self.rules)]
        for name, part in drift_parts.items():
            if name in sources:
                animations.append(TransformFromCopy(VGroup(*sources[name]), part))
            else:
                animations.append(FadeIn(part))
        self.play(*animations, *self.swap_title("The Heston PDE"), *self.swap_caption(r"Under the pricing measure $e^{-rt}V$ has no drift,\\so the $dt$ terms of $dV$ add up to $rV$."), run_time=1.0)
        self.wait(1.0)

        # rV moves across the equals sign, giving the PDE on two lines.
        line_one, line_two = pde_pieces([("time_derivative", "V_t")], [("eq", "="), ("zero", "0")], True)
        pde, pde_parts = make_two_lines(line_one, line_two, 50, 0.3, 1.0)
        self.play(
            FadeOut(self.model), FadeOut(self.ito_lines),
            *morph_pieces(drift_parts, pde_parts),
            *self.swap_caption("Moving $rV$ to the left gives the Heston PDE."),
            run_time=1.0,
        )
        self.wait(0.8)
        self.end_section()


# From the PDE to the two C++ lines that step one cell.
# Four sections, the time step, the differences, the reciprocals and the code.
class PdeToCode(DeckScene):
    def construct(self):
        self.time_step_section()
        self.differences_section()
        self.reciprocals_section()
        self.code_section()

    # Section 1 writes one explicit Euler step and its stability limit.
    # Two grids show the limit with real numbers.
    def time_step_section(self):
        self.next_section("Explicit Euler steps")

        # Time to expiry, the PDE in that time and one Euler step.
        tau_line = MathTex(r"\tau = T - t, \qquad V_\tau = L V", font_size=46, color=INK)
        tau_line.move_to(UP * 2.45)
        tau_line.align_to(LEFT * 6.4, LEFT)
        euler = MathTex(r"V^{n+1} = V^n + \Delta t \cdot L V^n", font_size=46, color=INK)
        euler.move_to(UP * 1.4)
        euler.align_to(LEFT * 6.4, LEFT)
        euler_box = SurroundingRectangle(euler, color=BLUE_ACCENT, buff=0.15, stroke_width=4)

        # The largest stable step, set by the corner cell with the biggest diffusion.
        stability = MathTex(r"\Delta t_{\text{stable}} = \frac{1}{v_{\max} S_{\max}^2 / \Delta S^2 + \xi^2 v_{\max} / \Delta v^2 + r}", font_size=42, color=INK)
        stability.move_to(UP * 0.0)
        stability.align_to(LEFT * 6.4, LEFT)

        # The limit on the test grid and on the full grid.
        grid_rows = VGroup()
        row_texts = [
            (r"\text{test grid } 421 \times 51", r"\Delta t = 4.464 \times 10^{-6}, \quad \text{limit } 5.66 \times 10^{-6}"),
            (r"\text{full grid } 2048 \times 512", r"\text{limit } 3.66 \times 10^{-7}, \quad 980{,}000 \text{ steps}"),
        ]
        for index in range(2):
            label = MathTex(row_texts[index][0], font_size=38, color=GREY_TEXT)
            value = MathTex(row_texts[index][1], font_size=38, color=INK)
            label.move_to(UP * (-1.05 - 0.65 * index))
            label.align_to(LEFT * 6.4, LEFT)
            value.set_y(label.get_y())
            value.align_to(LEFT * 1.8, LEFT)
            grid_rows.add(VGroup(label, value))

        self.play(*self.swap_title("Explicit Euler steps"), Write(tau_line), run_time=0.7)
        self.play(Write(euler), Create(euler_box), run_time=0.8)
        self.play(Write(stability), *self.swap_caption("Each step is explicit, so no equation is solved per step."), run_time=0.9)
        self.play(FadeIn(grid_rows, lag_ratio=0.4), *self.swap_caption(r"Each step is explicit, so no equation is solved per step.\\The work is one stencil pass over the grid."), run_time=0.8)
        self.wait(0.6)
        self.end_section()
        self.time_step_mobjects = VGroup(tau_line, euler, euler_box, stability, grid_rows)

    # Section 2 shows each derivative as a central difference on a 3 by 3 stencil.
    # The cells a difference reads light up with their weights.
    def differences_section(self):
        self.next_section("Central differences on the grid")

        # The stencil with column labels below and row labels on the left.
        self.patch = make_patch(1.1, 1.1)
        self.patch.move_to(LEFT * 4.0 + UP * 0.25)
        column_labels = VGroup()
        for column, text in enumerate(["i-1", "i", "i+1"]):
            label = MathTex(text, font_size=36, color=INK)
            label.next_to(self.patch[6 + column], DOWN, buff=0.18)
            column_labels.add(label)
        row_labels = VGroup()
        for row, text in enumerate(["j+1", "j", "j-1"]):
            label = MathTex(text, font_size=36, color=INK)
            label.next_to(self.patch[row * 3], LEFT, buff=0.18)
            row_labels.add(label)
        stock_axis = Tex(r"share price $S$, step $\Delta S$", tex_template=SANS_TEMPLATE, font_size=30, color=GREY_TEXT)
        stock_axis.next_to(column_labels, DOWN, buff=0.15)
        variance_axis = Tex(r"variance $v$, step $\Delta v$", tex_template=SANS_TEMPLATE, font_size=30, color=GREY_TEXT)
        variance_axis.next_to(self.patch, UP, buff=0.2)
        variance_axis.align_to(row_labels, LEFT)
        self.patch_labels = VGroup(column_labels, row_labels, stock_axis, variance_axis)

        self.play(
            FadeOut(self.time_step_mobjects),
            *self.swap_title("Central differences on the grid"),
            *self.swap_caption(r"$V_{i,j}$ is the value at share price node $i$ and variance node $j$.", 30),
            FadeIn(self.patch), FadeIn(self.patch_labels),
            run_time=0.6,
        )

        # One formula per derivative, bars lined up in a column to the right of the stencil.
        self.formulas = []
        weights = VGroup()
        lit_cells = []
        bar_heights = [2.45, 1.3, 0.15, -1.0, -2.15]
        for index in range(len(CENTRAL_DIFFERENCES)):
            derivative, colour, numerator, denominator, cell_weights, reciprocal = CENTRAL_DIFFERENCES[index]
            formula = make_fraction_formula([("derivative", derivative, colour), ("approx", r"\approx")], [("numerator", numerator)], [("denominator", denominator)], [], 38)
            place_formula(formula, -1.9, bar_heights[index])
            self.formulas.append(formula)

            # Cells this difference reads, each with its weight.
            new_weights = VGroup()
            new_cells = []
            for (column_offset, row_offset), weight in cell_weights.items():
                cell = cell_at(self.patch, column_offset, row_offset)
                new_cells.append(cell)
                weight_label = MathTex(weight, font_size=40, color=INK)
                weight_label.move_to(cell)
                new_weights.add(weight_label)

            animations = [FadeOut(weights), FadeIn(new_weights), Write(formula["group"])]
            for cell in lit_cells:
                animations.append(cell.animate.set_fill(CARD_FILL))
            for cell in new_cells:
                animations.append(cell.animate.set_fill(HIGHLIGHT_FILL))
            self.play(*animations, run_time=0.6)
            self.wait(0.2)
            weights = new_weights
            lit_cells = new_cells

        # All nine cells together make the stencil.
        animations = [FadeOut(weights)]
        for cell in self.patch:
            animations.append(cell.animate.set_fill(HIGHLIGHT_FILL).set_stroke(RED_ACCENT))
        self.play(*animations, *self.swap_caption("The cross derivative needs the four corners, so the stencil has nine cells.", 30), run_time=0.6)
        self.wait(0.3)
        self.end_section()

    # Section 3 turns each division into a multiplication by a reciprocal.
    # The reciprocals are worked out once per solve.
    def reciprocals_section(self):
        self.next_section("Divisions become multiplications")

        # Each fraction unfolds into one line with a slash.
        division_lines = []
        animations = []
        for index in range(len(CENTRAL_DIFFERENCES)):
            derivative, colour, numerator, denominator, cell_weights, reciprocal = CENTRAL_DIFFERENCES[index]
            line, parts = make_equation([("derivative", derivative, colour), ("approx", r"\approx"), ("numerator", "(" + numerator + ")"), ("slash", "/"), ("denominator", "(" + denominator + ")")], 38)
            line.move_to(UP * (2.6 - 0.65 * index))
            line.align_to(LEFT * 6.4, LEFT)
            division_lines.append(parts)
            formula = self.formulas[index]
            animations.append(ReplacementTransform(formula["left"]["derivative"], parts["derivative"]))
            animations.append(ReplacementTransform(formula["left"]["approx"], parts["approx"]))
            animations.append(TransformMatchingShapes(formula["numerator"], parts["numerator"]))
            animations.append(ReplacementTransform(formula["bar"], parts["slash"]))
            animations.append(TransformMatchingShapes(formula["denominator"], parts["denominator"]))
        self.play(
            FadeOut(self.patch), FadeOut(self.patch_labels),
            *self.swap_title("Divisions become multiplications"),
            *self.swap_caption("The baseline divides five times at every cell in every step."),
            *animations,
            run_time=0.9,
        )
        self.wait(0.4)

        # The five reciprocals in a card, in the order of src/stencil.cpp lines 28 to 32.
        definitions = VGroup()
        definition_denominators = {}
        definition_rest = VGroup()
        for name, denominator in RECIPROCALS:
            name_text = Text(name, font=MONO_FONT, font_size=26, color=RED_ACCENT)
            value, value_parts = make_equation([("eq", r"= 1 /"), ("denominator", "(" + denominator + ")")], 34)
            value.next_to(name_text, RIGHT, buff=0.15)
            definitions.add(VGroup(name_text, value))
            definition_denominators[name] = value_parts["denominator"]
            definition_rest.add(name_text, value_parts["eq"])
        first_row = VGroup(*definitions[:3]).arrange(RIGHT, buff=0.55)
        second_row = VGroup(*definitions[3:]).arrange(RIGHT, buff=0.55)
        rows = VGroup(first_row, second_row).arrange(DOWN, buff=0.3)
        rows.move_to(DOWN * 1.15)
        card = make_card(rows, 0.2)

        animations = [FadeIn(card), FadeIn(definition_rest)]
        for index in range(len(CENTRAL_DIFFERENCES)):
            reciprocal = CENTRAL_DIFFERENCES[index][5]
            animations.append(TransformFromCopy(division_lines[index]["denominator"], definition_denominators[reciprocal]))
        self.play(*animations, *self.swap_caption(r"Each reciprocal is worked out once per solve,\\at \texttt{src/stencil.cpp} lines 28 to 32."), run_time=0.9)
        self.wait(0.4)

        # The slash and denominator of each line become a dot and the reciprocal name.
        animations = []
        for index in range(len(CENTRAL_DIFFERENCES)):
            reciprocal = CENTRAL_DIFFERENCES[index][5]
            parts = division_lines[index]
            dot = MathTex(r"\cdot", font_size=38, color=INK)
            dot.next_to(parts["numerator"], RIGHT, buff=0.15)
            dot.set_y(parts["approx"].get_y())
            name_text = Text(reciprocal, font=MONO_FONT, font_size=28, color=RED_ACCENT)
            name_text.next_to(dot, RIGHT, buff=0.15)
            animations.append(FadeTransform(parts["slash"], dot))
            animations.append(FadeTransform(parts["denominator"], name_text))
        self.play(*animations, *self.swap_caption(r"The time loop multiplies by the reciprocals,\\at \texttt{src/stencil.cpp} lines 107 to 111."), run_time=0.9)
        self.wait(0.5)
        self.end_section()

    # Section 4 morphs the update into the two C++ lines that run it.
    # Each maths line turns into the code line beside it.
    def code_section(self):
        self.next_section("The update in C++")

        # The code as one Text, so the leading spaces keep the plus signs lined up.
        code = Text("\n".join(CODE_LINES), font=MONO_FONT, font_size=28, color=INK, t2c=CODE_COLOURS, line_spacing=1.2)
        code.move_to(RIGHT * 0.45 + UP * 0.42)
        code_line_groups = []
        start = 0
        for line in CODE_LINES:
            glyph_count = len(line.replace(" ", ""))
            code_line_groups.append(VGroup(*code[start:start + glyph_count]))
            start = start + glyph_count

        # Line numbers 121 and 122 in a grey gutter on the left.
        first_number = Text("121", font=MONO_FONT, font_size=28, color=GREY_TEXT)
        first_number.next_to(code_line_groups[0], LEFT, buff=0.45)
        second_number = Text("122", font=MONO_FONT, font_size=28, color=GREY_TEXT)
        second_number.next_to(code_line_groups[6], LEFT, buff=0.45)
        second_number.align_to(first_number, LEFT)
        card = make_card(VGroup(code, first_number, second_number), 0.25)

        # The update in maths, one line per term, each line starting where its code line starts.
        maths_pieces = [
            [("change", r"\text{change}"), ("eq", "="), ("term", r"\tfrac{1}{2} v S^2 V_{SS}", STOCK_DIFFUSION_COLOUR)],
            [("sign", "{}+"), ("term", r"\rho \xi v S V_{Sv}", CROSS_COLOUR)],
            [("sign", "{}+"), ("term", r"\tfrac{1}{2} \xi^2 v V_{vv}", VARIANCE_DIFFUSION_COLOUR)],
            [("sign", "{}+"), ("term", r"(r - q) S V_S", STOCK_DRIFT_COLOUR)],
            [("sign", "{}+"), ("term", r"\kappa(\theta - v) V_v", VARIANCE_DRIFT_COLOUR)],
            [("sign", "{}-"), ("term", r"r V_{i,j}", DISCOUNT_COLOUR)],
            [("new", r"V^{n+1}_{i,j}"), ("eq", "="), ("old", r"V^n_{i,j}"), ("plus", "+"), ("dt", r"\Delta t \cdot \text{change}")],
        ]
        maths_lines = VGroup()
        for index in range(len(maths_pieces)):
            line, parts = make_equation(maths_pieces[index], 38)
            line.set_y(code_line_groups[index].get_y())
            line.align_to(code_line_groups[index], LEFT)
            maths_lines.add(line)

        others = []
        for mobject in self.mobjects:
            if mobject is not self.title and mobject is not self.caption:
                others.append(mobject)
        self.play(
            *[FadeOut(mobject) for mobject in others],
            *self.swap_title("The update in C++"),
            *self.swap_caption("The update for one cell, with one line for each term."),
            FadeIn(card),
            run_time=0.6,
        )
        self.play(LaggedStart(*[Write(line) for line in maths_lines], lag_ratio=0.1), run_time=1.0)
        self.wait(0.5)

        # Each maths line morphs into its line of code.
        morphs = []
        for index in range(len(CODE_LINES)):
            morphs.append(FadeTransform(maths_lines[index], code_line_groups[index]))
        self.play(*morphs, run_time=1.0)
        self.play(
            FadeIn(first_number), FadeIn(second_number),
            *self.swap_caption(r"\texttt{src/stencil.cpp} lines 121 and 122,\\with line 121 wrapped at each plus sign."),
            run_time=0.5,
        )
        self.wait(0.5)
        self.end_section()


# One cell of the test grid worked through its first step with real numbers.
# Four sections, from the 3 by 3 patch to the new value of the cell and its neighbours.
class Differences(DeckScene):
    def construct(self):
        self.patch_section()
        self.differences_section()
        self.terms_section()
        self.update_section()

    # Section 1 shows the 3 by 3 patch around S = 5250, v = 0.04 at expiry.
    # The patch sits top right and the grid spacing on the left.
    def patch_section(self):
        self.next_section("One cell, first step back")

        # Nine cells, each with its name at the top and its value in the middle.
        self.patch = make_patch(1.8, 1.2)
        self.patch.move_to(RIGHT * 3.75 + UP * 0.75)
        names = VGroup()
        self.cell_values = {}
        for row in range(3):
            for column in range(3):
                cell = self.patch[row * 3 + column]
                name = Text(CELL_NAMES[row][column], font=MONO_FONT, font_size=19, color=GREY_TEXT)
                name.move_to(cell.get_top() + DOWN * 0.22)
                names.add(name)
                value = MathTex(CELL_VALUES[row][column], font_size=46, color=INK)
                value.move_to(cell.get_center() + DOWN * 0.12)
                self.cell_values[CELL_NAMES[row][column]] = value
        values = VGroup(*self.cell_values.values())
        centre_outline = Rectangle(width=1.8, height=1.2, stroke_color=RED_ACCENT, stroke_width=7)
        centre_outline.move_to(self.patch[4])

        # Share prices under the columns and variances left of the rows.
        labels = VGroup()
        for column, text in enumerate(["5200", "5250", "5300"]):
            label = MathTex(text, font_size=34, color=INK)
            label.next_to(self.patch[6 + column], DOWN, buff=0.15)
            labels.add(label)
        stock_letter = MathTex("S", font_size=38, color=GREY_TEXT)
        stock_letter.next_to(labels[2], RIGHT, buff=0.3)
        row_labels = VGroup()
        for row, text in enumerate(["0.06", "0.04", "0.02"]):
            label = MathTex(text, font_size=34, color=INK)
            label.next_to(self.patch[row * 3], LEFT, buff=0.15)
            row_labels.add(label)
        variance_letter = MathTex("v", font_size=38, color=GREY_TEXT)
        variance_letter.next_to(row_labels[0], UP, buff=0.3)
        labels.add(stock_letter, row_labels, variance_letter)

        # Test grid spacing and the strike on the left.
        self.parameters = VGroup(
            MathTex(r"\Delta S = 50, \quad \Delta v = 0.02", font_size=38, color=INK),
            MathTex(r"\Delta t = 0.25 / 56000 \approx 4.464 \times 10^{-6}", font_size=38, color=INK),
            MathTex(r"K = 5250", font_size=38, color=INK),
        )
        self.parameters.arrange(DOWN, aligned_edge=LEFT, buff=0.4)
        self.parameters.move_to(UP * 1.1)
        self.parameters.align_to(LEFT * 6.4, LEFT)

        self.play(*self.swap_title("One cell, first step back"), Create(self.patch), FadeIn(labels), run_time=0.6)
        self.play(FadeIn(names), LaggedStart(*[FadeIn(value, scale=1.3) for value in values], lag_ratio=0.05), run_time=0.8)
        self.play(
            Create(centre_outline), FadeIn(self.parameters, lag_ratio=0.2),
            *self.swap_caption(r"The cell at $S = 5250$, $v = 0.04$ and its eight neighbours.\\At expiry every cell holds the payoff $\max(S - K,\, 0)$."),
            run_time=0.7,
        )
        self.wait(0.6)
        self.end_section()
        self.patch_mobjects = VGroup(self.patch, names, values, centre_outline, labels)

    # Section 2 works out the five differences.
    # Each number flies from its cell into the formula.
    def differences_section(self):
        self.next_section("The five differences")
        self.play(FadeOut(self.parameters), *self.swap_title("The five differences"), *self.swap_caption("Each number comes from its cell in the patch."), run_time=0.4)

        cell_lookup = {}
        for row in range(3):
            for column in range(3):
                cell_lookup[CELL_NAMES[row][column]] = self.patch[row * 3 + column]

        self.formulas = []
        lit_cells = []
        bar_heights = [2.3, 1.2, 0.1, -1.0, -2.1]
        for index in range(len(NUMBER_DIFFERENCES)):
            derivative, colour, numerator_pieces, denominator, result = NUMBER_DIFFERENCES[index]
            formula = make_fraction_formula([("derivative", derivative, colour), ("eq", "=")], numerator_pieces, [("denominator", denominator)], [("eq", "="), ("result", result)], 40)
            place_formula(formula, -6.4, bar_heights[index])
            self.formulas.append(formula)

            # Numbers named after a cell fly in from that cell, the rest fade in.
            animations = []
            new_cells = []
            for name, part in formula["numerator_parts"].items():
                if name in self.cell_values:
                    animations.append(TransformFromCopy(self.cell_values[name], part))
                    new_cells.append(cell_lookup[name])
                else:
                    animations.append(FadeIn(part))
            for cell in lit_cells:
                animations.append(cell.animate.set_fill(CARD_FILL))
            for cell in new_cells:
                animations.append(cell.animate.set_fill(HIGHLIGHT_FILL))
            animations.append(FadeIn(formula["group"][0]))
            animations.append(FadeIn(formula["bar"]))
            animations.append(FadeIn(formula["denominator"]))
            self.play(*animations, run_time=0.5)
            self.play(Write(formula["right"]), run_time=0.25)
            lit_cells = new_cells

        # The three variance differences are zero.
        boxes = VGroup()
        for index in [2, 3, 4]:
            boxes.add(SurroundingRectangle(self.formulas[index]["right"], color=AMBER_ACCENT, buff=0.1, stroke_width=4))
        animations = []
        for cell in lit_cells:
            animations.append(cell.animate.set_fill(CARD_FILL))
        self.play(Create(boxes), *animations, *self.swap_caption("Every variance difference is 0 because every row holds the same payoff.", 30), run_time=0.5)
        self.wait(0.5)
        self.end_section()
        self.zero_boxes = boxes

    # Section 3 multiplies each difference by its coefficient.
    # Only the two share price terms are nonzero.
    def terms_section(self):
        self.next_section("The six terms")

        # The five results collapse into one summary line at the top.
        summary, summary_parts = make_equation([
            ("V_S", "V_S", STOCK_DRIFT_COLOUR), ("eq_s", "="), ("value_s", "0.5"), ("comma_s", r",\quad"),
            ("V_SS", "V_{SS}", STOCK_DIFFUSION_COLOUR), ("eq_ss", "="), ("value_ss", "0.02"), ("comma_ss", r",\quad"),
            ("V_v", "V_v", VARIANCE_DRIFT_COLOUR), ("eq_v", "="), ("V_vv", "V_{vv}", VARIANCE_DIFFUSION_COLOUR), ("eq_vv", "="),
            ("V_Sv", "V_{Sv}", CROSS_COLOUR), ("eq_sv", "="), ("zero", "0"),
        ], 44)
        summary.move_to(UP * 2.4)
        summary.align_to(LEFT * 6.4, LEFT)
        summary_keys = [("V_S", "eq_s", "value_s"), ("V_SS", "eq_ss", "value_ss"), ("V_v", "eq_v", None), ("V_vv", "eq_vv", None), ("V_Sv", "eq_sv", "zero")]
        animations = [FadeOut(self.patch_mobjects), FadeOut(self.zero_boxes)]
        for index in range(len(self.formulas)):
            formula = self.formulas[index]
            derivative_key, eq_key, value_key = summary_keys[index]
            animations.append(ReplacementTransform(formula["left"]["derivative"], summary_parts[derivative_key]))
            animations.append(ReplacementTransform(formula["right_parts"]["eq"], summary_parts[eq_key]))
            if value_key is None:
                animations.append(FadeOut(formula["right_parts"]["result"]))
            else:
                animations.append(ReplacementTransform(formula["right_parts"]["result"], summary_parts[value_key]))
            animations.append(FadeOut(formula["left"]["eq"]))
            animations.append(FadeOut(formula["numerator"]))
            animations.append(FadeOut(formula["bar"]))
            animations.append(FadeOut(formula["denominator"]))
        animations.append(FadeIn(summary_parts["comma_s"]))
        animations.append(FadeIn(summary_parts["comma_ss"]))
        self.play(*animations, *self.swap_title("The six terms"), *self.swap_caption("Each term multiplies one difference or the cell value by its coefficient.", 30), run_time=0.8)
        self.summary = summary

        # One line per term, left aligned under the summary.
        stock_diffusion_line, stock_diffusion_parts = make_equation([
            ("term", r"\tfrac{1}{2} v S^2 V_{SS}", STOCK_DIFFUSION_COLOUR), ("eq", "="),
            ("product", r"0.5 \times 0.04 \times 5250^2 \times"), ("difference", "0.02"),
            ("eq_result", "="), ("result", "11025"),
        ], 44)
        stock_drift_line, stock_drift_parts = make_equation([
            ("term", r"(r - q) S V_S", STOCK_DRIFT_COLOUR), ("eq", "="),
            ("product", r"0.032 \times 5250 \times"), ("difference", "0.5"),
            ("eq_result", "="), ("result", "84"),
        ], 44)
        zero_line, zero_parts = make_equation([
            ("cross", r"\rho \xi v S V_{Sv}", CROSS_COLOUR), ("zero_cross", r"= 0, \quad"),
            ("variance_diffusion", r"\tfrac{1}{2} \xi^2 v V_{vv}", VARIANCE_DIFFUSION_COLOUR), ("zero_vd", r"= 0, \quad"),
            ("variance_drift", r"\kappa(\theta - v) V_v", VARIANCE_DRIFT_COLOUR), ("zero_v", "= 0"),
        ], 44)
        discount_line, discount_parts = make_equation([
            ("term", r"-r V", DISCOUNT_COLOUR), ("eq", "="), ("product", r"-0.045 \times 0"), ("eq_result", "="), ("result", "0"),
        ], 44)
        term_lines = VGroup(stock_diffusion_line, stock_drift_line, zero_line, discount_line)
        for index in range(len(term_lines)):
            term_lines[index].move_to(UP * (1.3 - 1.05 * index))
            term_lines[index].align_to(LEFT * 6.4, LEFT)

        # The two share price terms, each taking its difference from the summary.
        for line_parts, value_key in [(stock_diffusion_parts, "value_ss"), (stock_drift_parts, "value_s")]:
            self.play(FadeIn(line_parts["term"]), FadeIn(line_parts["eq"]), FadeIn(line_parts["product"]), TransformFromCopy(summary_parts[value_key], line_parts["difference"]), run_time=0.6)
            self.play(Write(VGroup(line_parts["eq_result"], line_parts["result"])), run_time=0.4)
        self.play(Write(zero_line), Write(discount_line), *self.swap_caption("Only the two share price terms are nonzero."), run_time=0.7)
        self.wait(0.5)
        self.end_section()
        self.stock_diffusion_parts = stock_diffusion_parts
        self.stock_drift_parts = stock_drift_parts
        self.term_lines = term_lines

    # Section 4 adds the terms, takes one step of size dt and compares with the solver.
    # The cells above and below take different steps, because the curvature term grows with v.
    def update_section(self):
        self.next_section("The new value")

        # The two nonzero terms move into the sum.
        change_line, change_parts = make_equation([
            ("change", r"\text{change}"), ("eq", "="), ("first", "11025", STOCK_DIFFUSION_COLOUR), ("plus", "+"),
            ("second", "84", STOCK_DRIFT_COLOUR), ("eq_sum", "="), ("sum", "11109"),
        ], 46)
        change_line.move_to(UP * 2.4)
        change_line.align_to(LEFT * 6.4, LEFT)
        animations = [FadeOut(self.summary)]
        for line_parts, key in [(self.stock_diffusion_parts, "first"), (self.stock_drift_parts, "second")]:
            animations.append(ReplacementTransform(line_parts["result"], change_parts[key]))
            for name, part in line_parts.items():
                if name != "result":
                    animations.append(FadeOut(part))
        animations.append(FadeOut(self.term_lines[2]))
        animations.append(FadeOut(self.term_lines[3]))
        animations.append(FadeIn(change_parts["change"]))
        animations.append(FadeIn(change_parts["eq"]))
        animations.append(FadeIn(change_parts["plus"]))
        self.play(*animations, *self.swap_title("The new value"), run_time=0.7)
        self.play(Write(VGroup(change_parts["eq_sum"], change_parts["sum"])), run_time=0.4)

        # One explicit Euler step with the numbers in place.
        update_line = MathTex(r"V^{n+1}_{i,j} = 0 + 4.464 \times 10^{-6} \times 11109 \approx", "0.04959", font_size=46, color=INK)
        update_line[1].set_color(RED_ACCENT)
        update_line.move_to(UP * 1.4)
        update_line.align_to(LEFT * 6.4, LEFT)
        snapshot = Tex(r"The solver's own snapshot of this cell holds $0.04959375$.", tex_template=SANS_TEMPLATE, font_size=34, color=GREEN_ACCENT)
        snapshot.move_to(UP * 0.5)
        snapshot.align_to(LEFT * 6.25, LEFT)
        snapshot_box = SurroundingRectangle(snapshot, color=GREEN_ACCENT, buff=0.15, stroke_width=3, corner_radius=0.1)
        self.play(Write(update_line), *self.swap_caption("A cell worth nothing at expiry is worth about 5 cents one step earlier."), run_time=0.7)
        self.play(FadeIn(snapshot), Create(snapshot_box), run_time=0.5)
        self.wait(0.4)

        # The column of three cells after the same step.
        table_texts = [
            [r"v", r"\tfrac{1}{2} v S^2 V_{SS}", r"\text{new value}"],
            [r"0.06", r"16537.5", r"0.07420"],
            [r"0.04", r"11025", r"0.04959"],
            [r"0.02", r"5512.5", r"0.02498"],
        ]
        table = VGroup()
        for row in range(4):
            for column in range(3):
                entry = MathTex(table_texts[row][column], font_size=38, color=INK)
                if row == 0:
                    entry.set_color(GREY_TEXT)
                if row == 2:
                    entry.set_color(RED_ACCENT)
                if row > 0 and column == 1:
                    entry.set_color(STOCK_DIFFUSION_COLOUR)
                table.add(entry)
        table.arrange_in_grid(rows=4, cols=3, buff=(0.9, 0.22), col_alignments="lrr")
        table.move_to(DOWN * 1.05)
        table.align_to(LEFT * 6.25, LEFT)
        self.play(
            FadeIn(table, lag_ratio=0.05),
            *self.swap_caption(r"In the same step the cell above becomes 0.07420 and the one below 0.02498,\\because the $\tfrac{1}{2} v S^2 V_{SS}$ term grows with $v$.", 28),
            run_time=0.7,
        )
        self.wait(0.6)
        self.end_section()


# The whole solve on the test grid, drawn from the solver's own snapshots.
# Three sections, the first 60 steps, the full run and the price read off today.
class FullRun(DeckScene):
    def construct(self):
        self.load_snapshots()
        self.first_steps_section()
        self.whole_run_section()
        self.price_section()

    # Reads the snapshots, keeps the part of the grid the heatmap shows and builds one image per snapshot.
    # Snapshot k is step k + 1 for the first 60, then step 560 times (k - 59).
    def load_snapshots(self):
        paths = []
        for step in range(1, 61):
            paths.append(os.path.join(SNAPSHOT_FOLDER, "early", "snap_%06d.csv" % step))
        for dump in range(1, 101):
            paths.append(os.path.join(SNAPSHOT_FOLDER, "dump", "snap_%06d.csv" % (560 * dump)))

        stock_prices = np.arange(FIRST_COLUMN, LAST_COLUMN + 1) * STOCK_SPACING
        payoff = np.maximum(stock_prices - STRIKE, 0.0)
        self.stock_prices = stock_prices
        self.payoff = payoff
        self.value_rows = []
        self.colour_frames = []
        for path in paths:
            grid = np.loadtxt(path, delimiter=",")
            crop = grid[0:TOP_ROW + 1, FIRST_COLUMN:LAST_COLUMN + 1]
            self.value_rows.append(crop[2])
            self.colour_frames.append(time_value_colours(crop - payoff))
        self.final_grid = grid

    # The step number of snapshot k.
    # The first 60 snapshots are every step, the rest every 560 steps.
    def step_of(self, index):
        if index < 60:
            return index + 1
        return 560 * (index - 59)

    # The snapshot index the tracker points at.
    # The tracker moves smoothly, so it is rounded to a whole snapshot.
    def current_index(self):
        return int(round(self.index_tracker.get_value()))

    # Centre of one heatmap cell, counted from the left column and the bottom row of the crop.
    # Used for the outlines of the worked cell and of today's node.
    def heatmap_cell_centre(self, column, row):
        cell_width = self.heatmap.width / (LAST_COLUMN - FIRST_COLUMN + 1)
        cell_height = self.heatmap.height / (TOP_ROW + 1)
        x = self.heatmap.get_left()[0] + (column + 0.5) * cell_width
        y = self.heatmap.get_bottom()[1] + (row + 0.5) * cell_height
        return np.array([x, y, 0.0])

    # Section 1 runs the first 60 steps on the heatmap.
    # The worked cell at S = 5250, v = 0.04 counts up from 0.0496.
    def first_steps_section(self):
        self.next_section("The solver running")
        self.index_tracker = ValueTracker(0)

        # The heatmap of time value, one cell per grid node, scaled up without smoothing.
        self.heatmap = ImageMobject(self.colour_frames[0])
        self.heatmap.set_resampling_algorithm(RESAMPLING_ALGORITHMS["nearest"])
        self.heatmap.stretch_to_fit_width(6.3)
        self.heatmap.stretch_to_fit_height(2.25)
        self.heatmap.move_to(RIGHT * (-5.75 + 3.15) + UP * (0.35 + 1.125))

        # Each frame shows the snapshot the tracker points at.
        def show_snapshot(image):
            image.pixel_array = self.colour_frames[self.current_index()]
        self.heatmap.add_updater(show_snapshot)
        border = SurroundingRectangle(self.heatmap, color=INK, buff=0, stroke_width=2)

        # Ticks for S under the heatmap and for v on its left.
        ticks = VGroup()
        for text, column in [("4000", 0), ("5250", 25), ("6600", 52)]:
            tick = MathTex(text, font_size=26, color=INK)
            tick.move_to(self.heatmap_cell_centre(column, 0) * RIGHT + UP * 0.1)
            ticks.add(tick)
        for text, row in [("0", 0), ("0.1", 5), ("0.2", 10)]:
            tick = MathTex(text, font_size=26, color=INK)
            tick.move_to(self.heatmap_cell_centre(0, row) * UP)
            tick.align_to(RIGHT * -5.9, RIGHT)
            ticks.add(tick)
        stock_letter = MathTex("S", font_size=32, color=GREY_TEXT)
        stock_letter.move_to(RIGHT * 1.0 + UP * 0.1)
        variance_letter = MathTex("v", font_size=32, color=GREY_TEXT)
        variance_letter.move_to(RIGHT * -6.15 + UP * 2.85)
        heading = Tex(r"time value $V - \max(S - K,\,0)$", tex_template=SANS_TEMPLATE, font_size=28, color=GREY_TEXT)
        heading.move_to(UP * 2.85)
        heading.align_to(self.heatmap, LEFT)

        # Colour bar to the right of the heatmap, square root scale from 0 to 450.
        bar_values = np.linspace(MAX_TIME_VALUE, 0.0, 200).reshape(200, 1)
        colour_bar = ImageMobject(time_value_colours(bar_values)[::-1])
        colour_bar.stretch_to_fit_width(0.18)
        colour_bar.stretch_to_fit_height(2.25)
        colour_bar.move_to(RIGHT * 0.9 + UP * self.heatmap.get_y())
        colour_bar_border = SurroundingRectangle(colour_bar, color=INK, buff=0, stroke_width=1.5)
        colour_ticks = VGroup()
        for value in [0, 50, 200, 450]:
            tick = MathTex(str(value), font_size=22, color=INK)
            y = colour_bar.get_bottom()[1] + np.sqrt(value / MAX_TIME_VALUE) * colour_bar.height
            tick.move_to(UP * y)
            tick.align_to(RIGHT * (colour_bar.get_right()[0] + 0.1), LEFT)
            colour_ticks.add(tick)
        self.heatmap_frame = Group(border, ticks, stock_letter, variance_letter, heading, colour_bar, colour_bar_border, colour_ticks)

        # The worked cell outlined in red.
        cell_width = self.heatmap.width / (LAST_COLUMN - FIRST_COLUMN + 1)
        cell_height = self.heatmap.height / (TOP_ROW + 1)
        self.worked_outline = Rectangle(width=cell_width, height=cell_height, stroke_color=RED_ACCENT, stroke_width=4)
        self.worked_outline.move_to(self.heatmap_cell_centre(105 - FIRST_COLUMN, 2))

        # Read-outs under the heatmap, each number kept just right of its label.
        step_label = MathTex(r"\text{step } n =", font_size=36, color=INK)
        tau_label = MathTex(r"\tau = n \times 4.464 \times 10^{-6}\ \text{years} =", font_size=36, color=INK)
        cell_label = MathTex(r"V(5250,\, 0.04) =", font_size=36, color=RED_ACCENT)
        labels = [step_label, tau_label, cell_label]
        for index in range(3):
            labels[index].move_to(UP * (-0.55 - 0.62 * index))
            labels[index].align_to(LEFT * 6.4, LEFT)
        step_number = Integer(1, font_size=36, color=INK)
        tau_number = DecimalNumber(4.464e-6, num_decimal_places=6, font_size=36, color=INK)
        cell_number = DecimalNumber(self.value_rows[0][105 - FIRST_COLUMN], num_decimal_places=4, font_size=36, color=RED_ACCENT)

        # Height of a digit centre above the centre of an equals sign, so numbers sit on the label baseline.
        reference = MathTex("=", "0", font_size=36)
        digit_offset = reference[1].get_y() - reference[0].get_y()

        # Each number follows the tracker and sits on the baseline of its label.
        def follow(number, label, value_of_index):
            def update(mobject):
                mobject.set_value(value_of_index(self.current_index()))
                mobject.next_to(label, RIGHT, buff=0.15)
                mobject.set_y(label[0][-1].get_y() + digit_offset)
            number.add_updater(update)
            update(number)
        follow(step_number, step_label, lambda index: self.step_of(index))
        follow(tau_number, tau_label, lambda index: self.step_of(index) * 0.25 / 56000)
        follow(cell_number, cell_label, lambda index: self.value_rows[index][105 - FIRST_COLUMN])
        self.readouts = VGroup(step_label, tau_label, cell_label, step_number, tau_number, cell_number)

        self.play(
            *self.swap_title("The solver running"),
            FadeIn(self.heatmap), FadeIn(self.heatmap_frame), Create(self.worked_outline), FadeIn(self.readouts),
            *self.swap_caption(r"The first 60 steps. Time value appears at the strike first,\\and the worked cell holds 0.0496 after step 1."),
            run_time=0.6,
        )
        self.wait(0.4)
        self.play(self.index_tracker.animate.set_value(59), rate_func=linear, run_time=3.5)
        self.wait(0.3)
        self.end_section()

    # Section 2 runs the whole solve to step 56,000 beside a plot of V at v = 0.04.
    # The caption states the rule on each of the four edges of the grid.
    def whole_run_section(self):
        self.next_section("The whole run to today")

        # Axes for V against S at v = 0.04, with hand placed labels.
        self.axes = Axes(x_range=[4000, 6600, 2600], y_range=[0, 1400, 1400], x_length=3.9, y_length=2.25, tips=False, axis_config={"color": INK, "stroke_width": 2, "include_ticks": False})
        self.axes.move_to(RIGHT * (2.55 + 1.95) + UP * (0.35 + 1.125))
        axis_labels = VGroup()
        for value in [4000, 5250, 6600]:
            label = MathTex(str(value), font_size=24, color=INK)
            label.move_to(self.axes.c2p(value, 0) + DOWN * 0.25)
            axis_labels.add(label)
        for value in [0, 700, 1400]:
            label = MathTex(str(value), font_size=24, color=INK)
            label.next_to(self.axes.c2p(4000, value), LEFT, buff=0.1)
            axis_labels.add(label)
        plot_heading = Tex(r"$V$ against $S$ at $v = 0.04$", tex_template=SANS_TEMPLATE, font_size=28, color=GREY_TEXT)
        plot_heading.move_to(UP * 2.85)
        plot_heading.align_to(self.axes, LEFT)

        # The payoff as a dashed grey kink.
        payoff_line = VMobject(stroke_color=GREY_TEXT, stroke_width=3)
        payoff_line.set_points_as_corners([self.axes.c2p(4000, 0), self.axes.c2p(STRIKE, 0), self.axes.c2p(6600, 6600 - STRIKE)])
        payoff_dashes = DashedVMobject(payoff_line, num_dashes=40)

        # The value curve and the worked cell on it, redrawn from the current snapshot.
        def value_curve():
            row = self.value_rows[self.current_index()]
            points = []
            for index in range(len(self.stock_prices)):
                points.append(self.axes.c2p(self.stock_prices[index], max(row[index], 0.0)))
            curve = VMobject(stroke_color=BLUE_ACCENT, stroke_width=4)
            curve.set_points_as_corners(points)
            return curve
        def worked_dot():
            value = self.value_rows[self.current_index()][105 - FIRST_COLUMN]
            return Dot(self.axes.c2p(STRIKE, value), radius=0.06, color=RED_ACCENT)
        self.curve = always_redraw(value_curve)
        self.curve_dot = always_redraw(worked_dot)
        self.plot = VGroup(self.axes, axis_labels, plot_heading, payoff_dashes)

        self.play(
            *self.swap_title("The whole run to today"),
            FadeIn(self.plot), FadeIn(self.curve), FadeIn(self.curve_dot),
            *self.swap_caption(r"At the edges the $v = 0$ row keeps only the drift terms,\\a call at $S = 0$ is worth 0, the largest $S$ is worth $S e^{-q\tau} - K e^{-r\tau}$,\\and the top variance row copies the row below.", 28),
            run_time=0.6,
        )
        self.play(self.index_tracker.animate.set_value(159), rate_func=linear, run_time=4.5)
        self.wait(0.4)
        self.end_section()

    # Section 3 reads today's price at S0 = 5200, v0 = 0.04 with the Greeks from central differences.
    # The two cross-checks follow in the caption.
    def price_section(self):
        self.next_section("The price today")

        # Price and Greeks at node i = 104, j = 2 of the final snapshot.
        grid = self.final_grid
        price = grid[2, 104]
        delta = (grid[2, 105] - grid[2, 103]) / (2.0 * 50.0)
        gamma = (grid[2, 105] - 2.0 * price + grid[2, 103]) / (50.0 * 50.0)
        vega = (grid[3, 104] - grid[1, 104]) / (2.0 * 0.02)

        # Today's node outlined in green on the heatmap and marked on the plot.
        node_outline = self.worked_outline.copy().set_stroke(GREEN_ACCENT)
        node_outline.move_to(self.heatmap_cell_centre(104 - FIRST_COLUMN, 2))
        node_dot = Dot(self.axes.c2p(5200, price), radius=0.08, color=GREEN_ACCENT)
        node_label = MathTex("%.4f" % price, font_size=28, color=GREEN_ACCENT)
        node_label.next_to(node_dot, LEFT, buff=0.12).shift(UP * 0.2)

        # The read-outs become the price and the Greeks.
        results = VGroup(
            MathTex(r"V(5200,\, 0.04) = %.4f" % price, font_size=38, color=GREEN_ACCENT),
            MathTex(r"\Delta = V_S = %.4f, \quad \Gamma = V_{SS} = %.6f" % (delta, gamma), font_size=36, color=INK),
            MathTex(r"\text{vega} = V_v = %.1f \text{ per unit of variance}" % vega, font_size=36, color=INK),
        )
        for index in range(3):
            results[index].move_to(UP * (-0.55 - 0.62 * index))
            results[index].align_to(LEFT * 6.4, LEFT)

        self.readouts.clear_updaters()
        for mobject in self.readouts:
            mobject.clear_updaters()
        self.play(
            *self.swap_title("The price today"),
            ReplacementTransform(self.worked_outline, node_outline),
            FadeOut(self.curve_dot), FadeIn(node_dot), FadeIn(node_label),
            FadeOut(self.readouts), FadeIn(results, lag_ratio=0.2),
            *self.swap_caption(r"The price today is read at the grid node $S_0 = 5200$, $v_0 = 0.04$,\\with delta, gamma and vega from central differences."),
            run_time=0.9,
        )
        self.wait(1.2)
        self.play(*self.swap_caption(r"The 2048 by 512 grid gives 196.1684,\\and Fourier inversion gives 196.1692."), run_time=0.5)
        self.wait(1.0)
        self.end_section()


# Writes the stop times of one rendered video to a JSON file.
# Each stop is the end time in seconds of one section, read from manim's section index.
def write_stops(section_index_path, output_path):
    with open(section_index_path) as file:
        sections = json.load(file)
    stops = []
    total_frames = 0
    for section in sections:
        total_frames = total_frames + int(section["nb_frames"])
        stops.append(round(total_frames / 30, 3))
    with open(output_path, "w") as file:
        json.dump({"stops": stops}, file)
    print(stops)


if __name__ == "__main__":
    write_stops(sys.argv[1], sys.argv[2])
