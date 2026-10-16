# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.

# Two manim animations for the Milestone 2 deck, rendered at 1600 by 900 and 30 frames a second.
# Render commands, with MEDIA as any scratch folder and LaTeX on the PATH:
#   export PATH=/Library/TeX/texbin:$PATH
#   uv run --python 3.12 --with manim manim -qh -r 1600,900 --fps 30 --save_sections --media_dir MEDIA scenes.py Derivation
#   uv run --python 3.12 --with manim manim -qh -r 1600,900 --fps 30 --save_sections --media_dir MEDIA scenes.py Differences
# Then copy the two videos here and write the stop times of each one:
#   cp MEDIA/videos/scenes/900p30/Derivation.mp4 derivation.mp4
#   cp MEDIA/videos/scenes/900p30/Differences.mp4 differences.mp4
#   uv run --python 3.12 --with manim python scenes.py MEDIA/videos/scenes/900p30/sections/Derivation.json derivation.json
#   uv run --python 3.12 --with manim python scenes.py MEDIA/videos/scenes/900p30/sections/Differences.json differences.json

import json
import sys

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

# One colour for each of the six PDE terms, the same in both videos.
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

# Still time at the end of every section, where the deck pauses the video.
SECTION_END_PAUSE = 0.6

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

# The five central differences of the first video.
# Each entry is the derivative, its colour, the numerator, the denominator, the cell weights and the reciprocal name.
# A cell is (column offset, row offset) from the centre cell, with rows counted upwards.
CENTRAL_DIFFERENCES = [
    ("V_S", STOCK_DRIFT_COLOUR, r"V_{i+1,j} - V_{i-1,j}", r"2\Delta S", {(1, 0): "+1", (-1, 0): "-1"}, "inv_two_ds"),
    ("V_{SS}", STOCK_DIFFUSION_COLOUR, r"V_{i+1,j} - 2V_{i,j} + V_{i-1,j}", r"\Delta S^2", {(1, 0): "+1", (0, 0): "-2", (-1, 0): "+1"}, "inv_ds_sq"),
    ("V_v", VARIANCE_DRIFT_COLOUR, r"V_{i,j+1} - V_{i,j-1}", r"2\Delta v", {(0, 1): "+1", (0, -1): "-1"}, "inv_two_dv"),
    ("V_{vv}", VARIANCE_DIFFUSION_COLOUR, r"V_{i,j+1} - 2V_{i,j} + V_{i,j-1}", r"\Delta v^2", {(0, 1): "+1", (0, 0): "-2", (0, -1): "+1"}, "inv_dv_sq"),
    ("V_{Sv}", CROSS_COLOUR, r"V_{i+1,j+1} - V_{i-1,j+1} - V_{i+1,j-1} + V_{i-1,j-1}", r"4\Delta S \Delta v", {(1, 1): "+1", (-1, 1): "-1", (1, -1): "-1", (-1, -1): "+1"}, "inv_four_ds_dv"),
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

# Names and values of the nine cells in the second video, top row first.
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

# The five differences of the second video with real numbers.
# Each entry is the derivative, its colour, the numerator pieces, the denominator and the result.
# A numerator piece named after a cell gets its number from that cell.
NUMBER_DIFFERENCES = [
    ("V_S", STOCK_DRIFT_COLOUR, [("east", "50"), ("minus", "-"), ("west", "0")], r"2 \times 50", "0.5"),
    ("V_{SS}", STOCK_DIFFUSION_COLOUR, [("east", "50"), ("minus", "-"), ("two", r"2 \times"), ("centre", "0"), ("plus", "+"), ("west", "0")], r"50^2", "0.02"),
    ("V_v", VARIANCE_DRIFT_COLOUR, [("north", "0"), ("minus", "-"), ("south", "0")], r"2 \times 0.02", "0"),
    ("V_{vv}", VARIANCE_DIFFUSION_COLOUR, [("north", "0"), ("minus", "-"), ("two", r"2 \times"), ("centre", "0"), ("plus", "+"), ("south", "0")], r"0.02^2", "0"),
    ("V_{Sv}", CROSS_COLOUR, [("north_east", "50"), ("minus", "-"), ("north_west", "0"), ("minus_two", "-"), ("south_east", "50"), ("plus", "+"), ("south_west", "0")], r"4 \times 50 \times 0.02", "0"),
]


# A bold title in the top left corner of the frame.
# Every section opens by setting one.
def make_title(text):
    title = Tex(r"\textbf{" + text + "}", tex_template=SANS_TEMPLATE, font_size=42, color=INK)
    title.to_corner(UL, buff=0.45)
    return title


# A grey caption in the bottom left corner of the frame.
# A double backslash in the text starts a second line.
def make_caption(text, font_size=34):
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


# From the Heston model to the two C++ lines that step one cell.
# Six sections, one for each step of the derivation.
class Derivation(Scene):
    def construct(self):
        self.model_section()
        self.pde_section()
        self.time_step_section()
        self.differences_section()
        self.reciprocals_section()
        self.code_section()

    # Swaps the title and the caption for new ones.
    # Returns the two animations so they can run with others.
    def swap_title_and_caption(self, title_text, caption_text):
        new_title = make_title(title_text)
        new_caption = make_caption(caption_text)
        animations = [FadeOut(self.title), FadeIn(new_title), FadeOut(self.caption), FadeIn(new_caption)]
        self.title = new_title
        self.caption = new_caption
        return animations

    # Swaps only the caption for a new one.
    # Returns the two animations so they can run with others.
    def swap_caption(self, caption_text):
        new_caption = make_caption(caption_text)
        animations = [FadeOut(self.caption), FadeIn(new_caption)]
        self.caption = new_caption
        return animations

    # Section 1 shows the two random processes of the Heston model.
    # The three equations line up on their equals signs.
    def model_section(self):
        self.next_section("The Heston model")
        self.title = make_title("The Heston model")
        self.caption = make_caption("The share price and its variance both move randomly.")

        # The share price, the variance and the correlation of the two random drivers.
        stock_equation, self.stock_parts = make_equation([
            ("lhs", "dS"), ("eq", "="),
            ("drift", r"(r - q) S", STOCK_DRIFT_COLOUR), ("dt", r"\, dt"), ("plus", "+"),
            ("vol", r"\sqrt{v}\, S", STOCK_DIFFUSION_COLOUR), ("dW", r"\, dW_S"),
        ], 52)
        variance_equation, self.variance_parts = make_equation([
            ("lhs", "dv"), ("eq", "="),
            ("drift", r"\kappa(\theta - v)", VARIANCE_DRIFT_COLOUR), ("dt", r"\, dt"), ("plus", "+"),
            ("vol", r"\xi \sqrt{v}", VARIANCE_DIFFUSION_COLOUR), ("dW", r"\, dW_v"),
        ], 52)
        correlation_equation, self.correlation_parts = make_equation([
            ("lhs", r"dW_S\, dW_v"), ("eq", "="), ("rho", r"\rho", CROSS_COLOUR), ("dt", r"\, dt"),
        ], 52)
        self.model = VGroup(stock_equation, variance_equation, correlation_equation)

        # Equals signs on one vertical line, a grey label to the right of each equation.
        labels = VGroup()
        rows = [(stock_equation, self.stock_parts, 1.5, "share price"), (variance_equation, self.variance_parts, 0.3, "variance"), (correlation_equation, self.correlation_parts, -0.9, "correlation")]
        for equation, parts, y, label_text in rows:
            equation.shift(UP * (y - parts["eq"].get_y()))
            equation.shift(RIGHT * (-2.6 - parts["eq"].get_x()))
            label = Tex(label_text, tex_template=SANS_TEMPLATE, font_size=34, color=GREY_TEXT)
            label.move_to(UP * y)
            label.align_to(RIGHT * 3.7, LEFT)
            labels.add(label)
        self.model_labels = labels

        self.play(FadeIn(self.title), run_time=0.6)
        for index in range(3):
            self.play(Write(self.model[index]), FadeIn(labels[index]), run_time=1.3)
        self.play(FadeIn(self.caption), run_time=0.6)
        self.wait(SECTION_END_PAUSE)

    # Section 2 builds the PDE one term at a time.
    # Each term comes from the model piece boxed in the same colour.
    def pde_section(self):
        self.next_section("From the model to a PDE")

        # Model shrinks to the top of the frame.
        stock_equation, variance_equation, correlation_equation = self.model
        self.play(
            stock_equation.animate.scale(0.75).move_to(LEFT * 3.4 + UP * 2.45),
            variance_equation.animate.scale(0.75).move_to(RIGHT * 3.4 + UP * 2.45),
            correlation_equation.animate.scale(0.75).move_to(UP * 1.6),
            FadeOut(self.model_labels),
            *self.swap_title_and_caption("From the model to a PDE", "The Feynman-Kac theorem turns this expected discounted payoff into a PDE."),
            run_time=1.2,
        )

        # The option value as an expected discounted payoff.
        self.expectation, self.expectation_parts = make_equation([
            ("value", r"V(S, v, t) ="), ("open", r"\mathbb{E}\big[\,"),
            ("discount", r"e^{-r(T-t)}", DISCOUNT_COLOUR), ("payoff", r"\max(S_T - K,\, 0)"),
            ("given", r"\;\big|\; S_t = S,\ v_t = v \,\big]"),
        ], 42)
        self.expectation.move_to(UP * 0.5)
        self.play(Write(self.expectation), run_time=1.5)

        # The PDE with only the time derivative and the zero showing at first.
        line_one, line_two = pde_pieces([("time_derivative", "V_t")], [("eq", "="), ("zero", "0")], True)
        self.pde, self.pde_parts = make_two_lines(line_one, line_two, 50, -1.2, 1.0)
        self.play(FadeIn(self.pde_parts["time_derivative"]), FadeIn(self.pde_parts["eq"]), FadeIn(self.pde_parts["zero"]), run_time=0.8)

        # Each step names the term, the model pieces it comes from and the caption.
        term_colours = {name: colour for name, tex, colour in PDE_TERMS}
        build_steps = [
            ("stock_drift", [self.stock_parts["drift"]], r"The drift of $S$ gives $(r - q)\,S\,V_S$."),
            ("stock_diffusion", [self.stock_parts["vol"]], r"The variance $v S^2$ of $S$ gives $\tfrac{1}{2} v S^2 V_{SS}$."),
            ("variance_drift", [self.variance_parts["drift"]], r"The drift of $v$ gives $\kappa(\theta - v)\,V_v$."),
            ("variance_diffusion", [self.variance_parts["vol"]], r"The variance $\xi^2 v$ of $v$ gives $\tfrac{1}{2} \xi^2 v\,V_{vv}$."),
            ("cross", [self.correlation_parts["rho"], self.stock_parts["vol"], self.variance_parts["vol"]], r"The correlation $\rho$ links the two and gives $\rho\,\xi\,v S\,V_{Sv}$."),
            ("discount", [self.expectation_parts["discount"]], r"Discounting at the rate $r$ gives $-rV$."),
        ]
        boxes = VGroup()
        for name, sources, caption_text in build_steps:
            new_boxes = VGroup()
            for source in sources:
                new_boxes.add(SurroundingRectangle(source, color=term_colours[name], buff=0.08, stroke_width=3))
            self.play(
                FadeOut(boxes),
                Create(new_boxes),
                TransformFromCopy(VGroup(*sources), self.pde_parts[name]),
                FadeIn(self.pde_parts["sign_" + name]),
                *self.swap_caption(caption_text),
                run_time=1.2,
            )
            boxes = new_boxes
            self.wait(0.6)

        self.play(FadeOut(boxes), *self.swap_caption("Each coloured term comes from one piece of the model."), run_time=0.8)
        self.wait(SECTION_END_PAUSE)

    # Section 3 turns the PDE into one explicit Euler step backwards from expiry.
    # The six terms become the operator L V.
    def time_step_section(self):
        self.next_section("Step back from expiry")

        # Model and expectation leave, then the PDE becomes V_tau = terms.
        self.play(
            FadeOut(self.model), FadeOut(self.expectation),
            *self.swap_title_and_caption("Step back from expiry", r"The solver starts at expiry and steps back in the time to expiry $\tau = T - t$."),
            run_time=0.7,
        )
        line_one, line_two = pde_pieces([("time_derivative", r"V_\tau"), ("eq", "=")], [], False)
        tau_pde, tau_parts = make_two_lines(line_one, line_two, 50, 0.6, 1.0)
        self.play(*morph_pieces(self.pde_parts, tau_parts), run_time=1.4)
        self.wait(0.8)

        # V_tau = L V on top, with L V spelled out under it.
        top_line, top_parts = make_equation([("time_derivative", r"V_\tau"), ("eq", "="), ("operator", "L V")], 50)
        top_line.move_to(UP * 2.2)
        line_one, line_two = pde_pieces([("operator_name", "L V"), ("operator_eq", "=")], [], False)
        operator_lines, operator_parts = make_two_lines(line_one, line_two, 50, 0.45, 1.0)
        new_parts = {}
        new_parts.update(top_parts)
        new_parts.update(operator_parts)
        self.play(*morph_pieces(tau_parts, new_parts), *self.swap_caption("The operator $L$ collects the six terms."), run_time=1.5)
        self.wait(0.8)

        # The time derivative becomes a forward difference, then one explicit Euler step.
        euler_difference = MathTex(r"\frac{V^{n+1} - V^n}{\Delta t}", "=", r"L V^n", font_size=50, color=INK)
        euler_difference.move_to(DOWN * 1.75)
        self.play(TransformMatchingShapes(top_line.copy(), euler_difference), *self.swap_caption(r"One explicit Euler step moves every cell back by $\Delta t$."), run_time=1.5)
        self.wait(0.5)
        euler_update = MathTex(r"V^{n+1}", "=", r"V^n", "+", r"\Delta t \cdot L V^n", font_size=54, color=INK)
        euler_update.move_to(DOWN * 1.75)
        self.play(TransformMatchingShapes(euler_difference, euler_update), run_time=1.5)
        euler_box = SurroundingRectangle(euler_update, color=BLUE_ACCENT, buff=0.2, stroke_width=4)
        self.play(Create(euler_box), run_time=0.6)
        self.wait(SECTION_END_PAUSE)
        self.time_step_mobjects = VGroup(top_line, operator_lines, euler_update, euler_box)

    # Section 4 shows each derivative as a central difference on a 3 by 3 stencil.
    # The cells a difference reads light up with their weights.
    def differences_section(self):
        self.next_section("Central differences on the grid")

        # The stencil with column labels below and row labels on the left.
        self.patch = make_patch(1.1, 1.1)
        self.patch.move_to(LEFT * 4.0 + UP * 0.05)
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
        stock_axis.next_to(column_labels, DOWN, buff=0.18)
        variance_axis = Tex(r"variance $v$, step $\Delta v$", tex_template=SANS_TEMPLATE, font_size=30, color=GREY_TEXT)
        variance_axis.next_to(self.patch, UP, buff=0.2)
        variance_axis.align_to(row_labels, LEFT)
        self.patch_labels = VGroup(column_labels, row_labels, stock_axis, variance_axis)

        self.play(
            FadeOut(self.time_step_mobjects),
            *self.swap_title_and_caption("Central differences on the grid", r"$V_{i,j}$ is the value at share price node $i$ and variance node $j$."),
            FadeIn(self.patch), FadeIn(self.patch_labels),
            run_time=1.2,
        )

        # One formula per derivative, bars lined up in a column to the right of the stencil.
        self.formulas = []
        weights = VGroup()
        lit_cells = []
        bar_heights = [2.35, 1.2, 0.05, -1.1, -2.25]
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
            self.play(*animations, run_time=1.4)
            self.wait(0.5)
            weights = new_weights
            lit_cells = new_cells

        # All nine cells together make the stencil.
        animations = [FadeOut(weights)]
        for cell in self.patch:
            animations.append(cell.animate.set_fill(HIGHLIGHT_FILL).set_stroke(RED_ACCENT))
        self.play(*animations, *self.swap_caption("The cross derivative needs the four corners, so the stencil has nine cells."), run_time=1.0)
        self.wait(SECTION_END_PAUSE)

    # Section 5 turns each division into a multiplication by a reciprocal.
    # The reciprocals are worked out once per solve.
    def reciprocals_section(self):
        self.next_section("Divisions become multiplications")

        # Each fraction unfolds into one line with a slash.
        division_lines = []
        animations = []
        for index in range(len(CENTRAL_DIFFERENCES)):
            derivative, colour, numerator, denominator, cell_weights, reciprocal = CENTRAL_DIFFERENCES[index]
            line, parts = make_equation([("derivative", derivative, colour), ("approx", r"\approx"), ("numerator", "(" + numerator + ")"), ("slash", "/"), ("denominator", "(" + denominator + ")")], 38)
            line.move_to(UP * (2.55 - 0.75 * index))
            line.align_to(LEFT * 6.6, LEFT)
            division_lines.append(parts)
            formula = self.formulas[index]
            animations.append(ReplacementTransform(formula["left"]["derivative"], parts["derivative"]))
            animations.append(ReplacementTransform(formula["left"]["approx"], parts["approx"]))
            animations.append(TransformMatchingShapes(formula["numerator"], parts["numerator"]))
            animations.append(ReplacementTransform(formula["bar"], parts["slash"]))
            animations.append(TransformMatchingShapes(formula["denominator"], parts["denominator"]))
        self.play(
            FadeOut(self.patch), FadeOut(self.patch_labels),
            *self.swap_title_and_caption("Divisions become multiplications", "The baseline divides five times at every cell in every step."),
            *animations,
            run_time=1.6,
        )
        self.wait(0.8)

        # The five reciprocals in a card, three on the first row and two on the second.
        definitions = VGroup()
        definition_parts = []
        for index in range(len(CENTRAL_DIFFERENCES)):
            derivative, colour, numerator, denominator, cell_weights, reciprocal = CENTRAL_DIFFERENCES[index]
            name = Text(reciprocal, font=MONO_FONT, font_size=26, color=RED_ACCENT)
            value, value_parts = make_equation([("eq", r"= 1 /"), ("denominator", "(" + denominator + ")")], 36)
            value.next_to(name, RIGHT, buff=0.15)
            definitions.add(VGroup(name, value))
            definition_parts.append((name, value_parts))
        first_row = VGroup(*definitions[:3]).arrange(RIGHT, buff=0.6)
        second_row = VGroup(*definitions[3:]).arrange(RIGHT, buff=0.6)
        rows = VGroup(first_row, second_row).arrange(DOWN, buff=0.35)
        rows.move_to(DOWN * 1.75)
        card = make_card(rows, 0.25)

        animations = [FadeIn(card)]
        for index in range(len(CENTRAL_DIFFERENCES)):
            name, value_parts = definition_parts[index]
            animations.append(TransformFromCopy(division_lines[index]["denominator"], value_parts["denominator"]))
            animations.append(FadeIn(name))
            animations.append(FadeIn(value_parts["eq"]))
        self.play(*animations, *self.swap_caption(r"Each reciprocal is worked out once per solve, \texttt{src/stencil.cpp} lines 28 to 32."), run_time=1.6)
        self.wait(0.8)

        # The slash and denominator of each line become a dot and the reciprocal name.
        animations = []
        for index in range(len(CENTRAL_DIFFERENCES)):
            derivative, colour, numerator, denominator, cell_weights, reciprocal = CENTRAL_DIFFERENCES[index]
            parts = division_lines[index]
            dot = MathTex(r"\cdot", font_size=38, color=INK)
            dot.next_to(parts["numerator"], RIGHT, buff=0.15)
            dot.set_y(parts["approx"].get_y())
            name = Text(reciprocal, font=MONO_FONT, font_size=28, color=RED_ACCENT)
            name.next_to(dot, RIGHT, buff=0.15)
            animations.append(FadeTransform(parts["slash"], dot))
            animations.append(FadeTransform(parts["denominator"], name))
        self.play(*animations, *self.swap_caption(r"The time loop multiplies by the reciprocals, \texttt{src/stencil.cpp} lines 107 to 111."), run_time=1.6)
        self.wait(SECTION_END_PAUSE)

    # Section 6 morphs the update into the two C++ lines that run it.
    # Each maths line turns into the code line beside it.
    def code_section(self):
        self.next_section("The update in C++")

        # The code as one Text, so the leading spaces keep the plus signs lined up.
        code = Text("\n".join(CODE_LINES), font=MONO_FONT, font_size=30, color=INK, t2c=CODE_COLOURS, line_spacing=1.2)
        code.move_to(RIGHT * 0.55 + DOWN * 0.05)
        code_line_groups = []
        start = 0
        for line in CODE_LINES:
            glyph_count = len(line.replace(" ", ""))
            code_line_groups.append(VGroup(*code[start:start + glyph_count]))
            start = start + glyph_count

        # Line numbers 121 and 122 in a grey gutter on the left.
        first_number = Text("121", font=MONO_FONT, font_size=30, color=GREY_TEXT)
        first_number.next_to(code_line_groups[0], LEFT, buff=0.5)
        second_number = Text("122", font=MONO_FONT, font_size=30, color=GREY_TEXT)
        second_number.next_to(code_line_groups[6], LEFT, buff=0.5)
        second_number.align_to(first_number, LEFT)
        card = make_card(VGroup(code, first_number, second_number), 0.35)

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
            line, parts = make_equation(maths_pieces[index], 44)
            line.set_y(code_line_groups[index].get_y())
            line.align_to(code_line_groups[index], LEFT)
            maths_lines.add(line)

        self.play(
            *[FadeOut(mobject) for mobject in self.mobjects if mobject is not self.title and mobject is not self.caption],
            *self.swap_title_and_caption("The update in C++", "The update for one cell, with one line for each term."),
            FadeIn(card),
            run_time=1.0,
        )
        self.play(LaggedStart(*[Write(line) for line in maths_lines], lag_ratio=0.25), run_time=2.4)
        self.wait(0.8)

        # Each maths line morphs into its line of code.
        morphs = []
        for index in range(len(CODE_LINES)):
            morphs.append(FadeTransform(maths_lines[index], code_line_groups[index]))
        self.play(*morphs, run_time=2.2)
        self.play(
            FadeIn(first_number), FadeIn(second_number),
            *self.swap_caption(r"\texttt{src/stencil.cpp}, lines 121 and 122. Line 121 is wrapped at each plus sign."),
            run_time=0.8,
        )
        self.wait(SECTION_END_PAUSE)


# One cell of the test grid worked through its first step with real numbers.
# Four sections, from the 3 by 3 patch to the new value of the cell.
class Differences(Scene):
    def construct(self):
        self.patch_section()
        self.differences_section()
        self.terms_section()
        self.update_section()

    # Swaps the title and the caption for new ones.
    # Returns the animations so they can run with others.
    def swap_title_and_caption(self, title_text, caption_text, font_size=34):
        new_title = make_title(title_text)
        new_caption = make_caption(caption_text, font_size)
        animations = [FadeOut(self.title), FadeIn(new_title), FadeOut(self.caption), FadeIn(new_caption)]
        self.title = new_title
        self.caption = new_caption
        return animations

    # Swaps only the caption for a new one.
    # Returns the two animations so they can run with others.
    def swap_caption(self, caption_text, font_size=34):
        new_caption = make_caption(caption_text, font_size)
        animations = [FadeOut(self.caption), FadeIn(new_caption)]
        self.caption = new_caption
        return animations

    # Section 1 shows the 3 by 3 patch around S = 5250, v = 0.04 at expiry.
    # The test grid parameters sit on the right.
    def patch_section(self):
        self.next_section("One cell, first step back from expiry")
        self.title = make_title("One cell, first step back from expiry")
        self.caption = make_caption(r"The cell at $S = 5250$ and $v = 0.04$ with its eight neighbours.\\At expiry every cell holds the payoff $\max(S - K,\, 0)$.", 32)

        # Nine cells, each with its name at the top and its value in the middle.
        self.patch = make_patch(1.9, 1.3)
        self.patch.move_to(LEFT * 2.95 + UP * 0.05)
        names = VGroup()
        self.cell_values = {}
        for row in range(3):
            for column in range(3):
                cell = self.patch[row * 3 + column]
                name = Text(CELL_NAMES[row][column], font=MONO_FONT, font_size=20, color=GREY_TEXT)
                name.move_to(cell.get_top() + DOWN * 0.24)
                names.add(name)
                value = MathTex(CELL_VALUES[row][column], font_size=50, color=INK)
                value.move_to(cell.get_center() + DOWN * 0.13)
                self.cell_values[CELL_NAMES[row][column]] = value
        values = VGroup(*self.cell_values.values())
        centre_outline = Rectangle(width=1.9, height=1.3, stroke_color=RED_ACCENT, stroke_width=7)
        centre_outline.move_to(self.patch[4])

        # Share prices under the columns and variances left of the rows.
        labels = VGroup()
        for column, text in enumerate(["5200", "5250", "5300"]):
            label = MathTex(text, font_size=36, color=INK)
            label.next_to(self.patch[6 + column], DOWN, buff=0.18)
            labels.add(label)
        stock_letter = MathTex("S", font_size=40, color=GREY_TEXT)
        stock_letter.next_to(labels[2], RIGHT, buff=0.35)
        row_labels = VGroup()
        for row, text in enumerate(["0.06", "0.04", "0.02"]):
            label = MathTex(text, font_size=36, color=INK)
            label.next_to(self.patch[row * 3], LEFT, buff=0.18)
            row_labels.add(label)
        variance_letter = MathTex("v", font_size=40, color=GREY_TEXT)
        variance_letter.next_to(row_labels[0], UP, buff=0.35)
        labels.add(stock_letter, row_labels, variance_letter)

        # Test grid parameters, one group per line.
        parameter_lines = VGroup(
            MathTex(r"\Delta S = 50, \quad \Delta v = 0.02", font_size=38, color=INK),
            MathTex(r"\Delta t = 0.25 / 56000 \approx 4.464 \times 10^{-6}", font_size=38, color=INK),
            MathTex(r"r = 0.045, \quad q = 0.013", font_size=38, color=INK),
            MathTex(r"\kappa = 1.5, \quad \theta = 0.04", font_size=38, color=INK),
            MathTex(r"\xi = 0.35, \quad \rho = -0.7", font_size=38, color=INK),
            MathTex(r"K = 5250", font_size=38, color=INK),
        )
        parameter_lines.arrange(DOWN, aligned_edge=LEFT, buff=0.38)
        parameter_lines.move_to(RIGHT * 3.55 + UP * 0.25)
        self.parameters = parameter_lines

        self.play(FadeIn(self.title), Create(self.patch), FadeIn(labels), run_time=1.2)
        self.play(FadeIn(names), LaggedStart(*[FadeIn(value, scale=1.3) for value in values], lag_ratio=0.08), run_time=1.6)
        self.play(Create(centre_outline), FadeIn(self.caption), run_time=0.8)
        self.wait(1.0)
        self.play(FadeIn(self.parameters, lag_ratio=0.15), run_time=1.2)
        self.wait(0.6)
        self.wait(SECTION_END_PAUSE)
        self.patch_mobjects = VGroup(self.patch, names, values, centre_outline, labels)

    # Section 2 works out the five differences.
    # Each number flies from its cell into the formula.
    def differences_section(self):
        self.next_section("The five differences")
        self.play(FadeOut(self.parameters), *self.swap_title_and_caption("The five differences", "Each difference takes its numbers from the cells of the patch."), run_time=0.8)

        cell_lookup = {}
        for row in range(3):
            for column in range(3):
                cell_lookup[CELL_NAMES[row][column]] = self.patch[row * 3 + column]

        self.formulas = []
        lit_cells = []
        bar_heights = [2.4, 1.25, 0.1, -1.05, -2.2]
        for index in range(len(NUMBER_DIFFERENCES)):
            derivative, colour, numerator_pieces, denominator, result = NUMBER_DIFFERENCES[index]
            formula = make_fraction_formula([("derivative", derivative, colour), ("eq", "=")], numerator_pieces, [("denominator", denominator)], [("eq", "="), ("result", result)], 42)
            place_formula(formula, 0.65, bar_heights[index])
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
            self.play(*animations, run_time=1.8)
            self.play(Write(formula["right"]), run_time=0.7)
            self.wait(0.6)
            lit_cells = new_cells

        # The three variance differences are zero.
        boxes = VGroup()
        for index in [2, 3, 4]:
            boxes.add(SurroundingRectangle(self.formulas[index]["right"], color=AMBER_ACCENT, buff=0.12, stroke_width=4))
        animations = []
        for cell in lit_cells:
            animations.append(cell.animate.set_fill(CARD_FILL))
        self.play(Create(boxes), *animations, *self.swap_caption("Every variance difference is 0 because the payoff is the same in every row.", 32), run_time=1.0)
        self.wait(SECTION_END_PAUSE)
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
        summary.move_to(UP * 2.45)
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
        self.play(*animations, *self.swap_title_and_caption("The six terms", "Each term multiplies one difference or the cell value by its coefficient."), run_time=1.6)
        self.summary = summary

        # One line per term, left aligned under the summary.
        stock_diffusion_line, stock_diffusion_parts = make_equation([
            ("term", r"\tfrac{1}{2} v S^2 V_{SS}", STOCK_DIFFUSION_COLOUR), ("eq", "="),
            ("product", r"0.5 \times 0.04 \times 5250^2 \times"), ("difference", "0.02"),
            ("eq_result", "="), ("result", "11025"),
        ], 46)
        stock_drift_line, stock_drift_parts = make_equation([
            ("term", r"(r - q) S V_S", STOCK_DRIFT_COLOUR), ("eq", "="),
            ("product", r"0.032 \times 5250 \times"), ("difference", "0.5"),
            ("eq_result", "="), ("result", "84"),
        ], 46)
        zero_line, zero_parts = make_equation([
            ("cross", r"\rho \xi v S V_{Sv}", CROSS_COLOUR), ("zero_cross", r"= 0, \quad"),
            ("variance_diffusion", r"\tfrac{1}{2} \xi^2 v V_{vv}", VARIANCE_DIFFUSION_COLOUR), ("zero_vd", r"= 0, \quad"),
            ("variance_drift", r"\kappa(\theta - v) V_v", VARIANCE_DRIFT_COLOUR), ("zero_v", "= 0"),
        ], 46)
        discount_line, discount_parts = make_equation([
            ("term", r"-r V", DISCOUNT_COLOUR), ("eq", "="), ("product", r"-0.045 \times 0"), ("eq_result", "="), ("result", "0"),
        ], 46)
        term_lines = VGroup(stock_diffusion_line, stock_drift_line, zero_line, discount_line)
        for index in range(len(term_lines)):
            term_lines[index].move_to(UP * (1.3 - 1.1 * index))
            term_lines[index].align_to(LEFT * 6.4, LEFT)

        # The two share price terms, each taking its difference from the summary.
        for line_parts, value_key in [(stock_diffusion_parts, "value_ss"), (stock_drift_parts, "value_s")]:
            self.play(FadeIn(line_parts["term"]), FadeIn(line_parts["eq"]), FadeIn(line_parts["product"]), TransformFromCopy(summary_parts[value_key], line_parts["difference"]), run_time=1.3)
            self.play(Write(VGroup(line_parts["eq_result"], line_parts["result"])), run_time=0.8)
            self.wait(0.6)
        self.play(Write(zero_line), run_time=1.2)
        self.play(Write(discount_line), *self.swap_caption("Only the two share price terms are nonzero."), run_time=1.0)
        self.wait(SECTION_END_PAUSE)
        self.stock_diffusion_parts = stock_diffusion_parts
        self.stock_drift_parts = stock_drift_parts
        self.term_lines = term_lines

    # Section 4 adds the terms and takes one step of size dt.
    # The result matches the solver's own snapshot of the cell.
    def update_section(self):
        self.next_section("The new value")

        # The two nonzero terms move into the sum.
        change_line, change_parts = make_equation([
            ("change", r"\text{change}"), ("eq", "="), ("first", "11025", STOCK_DIFFUSION_COLOUR), ("plus", "+"),
            ("second", "84", STOCK_DRIFT_COLOUR), ("eq_sum", "="), ("sum", "11109"),
        ], 50)
        change_line.move_to(UP * 1.7)
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
        self.play(*animations, *self.swap_title_and_caption("The new value", r"One explicit Euler step with $\Delta t = 4.464 \times 10^{-6}$ years."), run_time=1.5)
        self.play(Write(VGroup(change_parts["eq_sum"], change_parts["sum"])), run_time=0.8)
        self.wait(0.4)

        # The Euler step first in symbols, then with the numbers in place.
        update_line, update_parts = make_equation([
            ("new", r"V^{n+1}_{i,j}"), ("eq", "="), ("old", r"V^n_{i,j}"), ("plus", "+"),
            ("dt", r"\Delta t"), ("times", r"\cdot"), ("change", r"\text{change}"),
        ], 50)
        update_line.move_to(UP * 0.25)
        self.play(Write(update_line), run_time=1.2)
        self.wait(0.8)
        number_line, number_parts = make_equation([
            ("new", r"V^{n+1}_{i,j}"), ("eq", "="), ("old", "0"), ("plus", "+"),
            ("dt", r"4.464 \times 10^{-6}"), ("times", r"\times"), ("change", "11109"),
            ("approx", r"\approx"), ("result", "0.04959", RED_ACCENT),
        ], 50)
        number_line.move_to(UP * 0.25)
        number_without_result = {}
        for name, part in number_parts.items():
            if name != "approx" and name != "result":
                number_without_result[name] = part
        self.play(*morph_pieces(update_parts, number_without_result), run_time=1.4)
        self.play(Write(VGroup(number_parts["approx"], number_parts["result"])), run_time=0.8)
        self.wait(0.6)

        # The solver's own snapshot of the same cell.
        snapshot = Tex(r"The solver's own snapshot of this cell holds $0.04959375$.", tex_template=SANS_TEMPLATE, font_size=40, color=GREEN_ACCENT)
        snapshot.move_to(DOWN * 1.25)
        snapshot_box = SurroundingRectangle(snapshot, color=GREEN_ACCENT, buff=0.2, stroke_width=3, corner_radius=0.1)
        self.play(
            FadeIn(snapshot), Create(snapshot_box),
            *self.swap_caption(r"A cell worth nothing at expiry is worth about 5 cents one step earlier.\\Almost all of it comes from the curvature at the strike.", 32),
            run_time=1.0,
        )
        self.wait(SECTION_END_PAUSE)


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
