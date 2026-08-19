#!/usr/bin/env python3
"""Memory-behaviour animations (PLAN P7): why loop order and layout dominate.

Two styles, both driven by the same flat layout the solver actually uses
(idx = var_j*ns + stock_i, stock contiguous, see include/grid.h):

  --style cache    The array as memory. Two panels side by side sweep the
                   9-point stencil across it in the two possible orders:
                   row-major (what levels 0-6 do) and the ctl-order control
                   (loops swapped). A modelled fully-associative LRU cache of
                   64-byte lines runs underneath with a live hit/miss counter,
                   so the difference shows up as a number on screen instead of
                   a claim. Pass --slowdown with the measured ctl-order
                   slowdown from the benchmark so the caption quotes a real
                   result.

  --style buffers  The double-buffer ping-pong: read `current`, write `next`,
                   then std::swap exchanges two pointers in O(1) without
                   copying a single element. After that the per-level data
                   structures appear as the ladder climbs: row constants, the
                   S and S^2 tables, row base offsets, row pointers, and the
                   4-wide unroll.

Styling follows scripts/weather_map.py: dark background, amber #ffb000 for
the live stencil, cyan #7fd4ff for cache-resident memory, H.264 via imageio.

Usage:
    python3 scripts/memory_anim.py --style cache --out results/cache.mp4 \\
        --slowdown 31.4
    python3 scripts/memory_anim.py --style buffers --out results/buffers.mp4
"""

import argparse
import sys

import numpy as np

# Palette, shared with weather_map.py so the video looks like one artefact.
AMBER = "#ffb000"      # the live stencil / the cell being written
CYAN = "#7fd4ff"       # memory currently resident in cache
DIM = "#2a2a35"        # cold memory
MID = "#4b5563"        # grid lines, secondary text
WHITE = "#f5f5f5"
RED = "#ff5c5c"        # misses, the wrong-order panel
GREEN = "#5cd97b"      # hits

DOUBLES_PER_LINE = 8   # 64-byte cache line / 8 bytes per double


class ModelCache:
    """Fully-associative LRU cache of fixed-size lines.

    Deliberately tiny (a few tens of lines) so the contrast between the two
    traversals is visible in a 30-second clip. The real machine's L1 is 512
    lines; the MECHANISM being shown is identical, only the scale differs.
    """

    def __init__(self, num_lines):
        self.num_lines = num_lines
        self.lines = []      # most-recently-used last
        self.hits = 0
        self.misses = 0      # == lines fetched from memory
        self.cells = 0       # cell updates completed, for bytes-per-cell

    def access(self, flat_index):
        line = flat_index // DOUBLES_PER_LINE
        if line in self.lines:
            self.hits += 1
            self.lines.remove(line)     # refresh recency
            self.lines.append(line)
            return True
        self.misses += 1
        self.lines.append(line)
        if len(self.lines) > self.num_lines:
            self.lines.pop(0)           # evict least-recently-used
        return False

    @property
    def hit_rate(self):
        total = self.hits + self.misses
        return self.hits / total if total else 0.0

    @property
    def bytes_per_cell(self):
        """The number that actually drives the roofline: every miss drags a
        whole 64-byte line across the bus, however few of its 8 doubles get
        used before it is evicted."""
        if not self.cells:
            return 0.0
        return self.misses * DOUBLES_PER_LINE * 8 / self.cells

    def resident_doubles(self, total_doubles):
        """Flat indices currently held, for drawing."""
        held = np.zeros(total_doubles, dtype=bool)
        for line in self.lines:
            start = line * DOUBLES_PER_LINE
            held[start:min(start + DOUBLES_PER_LINE, total_doubles)] = True
        return held


def stencil_indices(stock_i, var_j, num_stock):
    """The nine cells the interior update reads, as flat indices. These are
    the ten accesses solver_opt_kernels.cpp makes: nine reads and one write,
    at the same flat position as the centre."""
    out = []
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            out.append((var_j + dj) * num_stock + (stock_i + di))
    return out


def visit_order(num_stock, num_var, row_major):
    """Cell visit sequence. row_major=True is the solver's order (outer
    variance, inner stock); False is the ctl-order control, which walks the
    variance direction on the inside and therefore jumps num_stock doubles
    per step."""
    if row_major:
        for var_j in range(1, num_var - 1):
            for stock_i in range(1, num_stock - 1):
                yield stock_i, var_j
    else:
        for stock_i in range(1, num_stock - 1):
            for var_j in range(1, num_var - 1):
                yield stock_i, var_j


def draw_panel(ax, num_stock, num_var, resident, live, centre, title,
               cache, colour):
    """One memory panel: the flat array drawn the way it is laid out, so
    reading left to right and then wrapping, which is row-major order."""
    import matplotlib.patches as patches

    total = num_stock * num_var
    image = np.zeros((num_var, num_stock, 3))
    cold = np.array([0.16, 0.16, 0.21])
    image[:, :] = cold
    resident_2d = resident[:total].reshape(num_var, num_stock)
    image[resident_2d] = np.array([0.30, 0.62, 0.80])   # cyan-ish: in cache
    for flat in live:
        if 0 <= flat < total:
            image[flat // num_stock, flat % num_stock] = np.array(
                [1.0, 0.69, 0.0])                        # amber: live stencil
    if centre is not None:
        image[centre // num_stock, centre % num_stock] = np.array(
            [1.0, 1.0, 1.0])                             # white: being written

    ax.imshow(image, interpolation="nearest", aspect="auto", origin="lower")
    # A cache line holds 8 doubles, and that is the unit memory moves in.
    for x in range(0, num_stock + 1, DOUBLES_PER_LINE):
        ax.axvline(x - 0.5, color=MID, linewidth=0.4, alpha=0.55)
    # pad clears the two counter lines drawn above the axes at y=1.045/1.115.
    ax.set_title(title, color=colour, fontsize=11.5, pad=46)
    ax.set_xlabel("stock_i  →  contiguous in memory (64 B = 8 doubles per tick)",
                  color=MID, fontsize=7.5)
    ax.set_ylabel("var_j  →  each row is num_stock doubles further on",
                  color=MID, fontsize=7.5)
    ax.tick_params(colors=MID, labelsize=6.5)
    for spine in ax.spines.values():
        spine.set_color(MID)

    rate = cache.hit_rate
    # The ideal for this stencil is 16 bytes per cell, meaning one double of
    # `cur` read and one double of `next` written with everything else served
    # from cache. Green means the traversal is achieving roughly that, and red
    # means it is dragging whole cache lines across the bus to use one double
    # out of each of them.
    good = cache.bytes_per_cell < 16.0
    ax.text(0.02, 1.115,
            f"bytes moved per cell update  {cache.bytes_per_cell:6.1f} B",
            transform=ax.transAxes, color=GREEN if good else RED,
            fontsize=10.5, fontfamily="monospace", fontweight="bold")
    ax.text(0.02, 1.045,
            f"lines fetched {cache.misses:,}   hits {cache.hits:,}   "
            f"hit rate {100 * rate:5.1f}%",
            transform=ax.transAxes, color=MID, fontsize=8.5,
            fontfamily="monospace")
    # A bar showing bytes/cell against a 64-byte reference (one whole line
    # per cell), so the two panels can be compared at a glance.
    fraction = min(1.0, cache.bytes_per_cell / 128.0)
    ax.add_patch(patches.Rectangle((0.02, -0.175), 0.96, 0.035,
                                   transform=ax.transAxes, facecolor=DIM,
                                   clip_on=False))
    ax.add_patch(patches.Rectangle((0.02, -0.175), 0.96 * fraction, 0.035,
                                   transform=ax.transAxes,
                                   facecolor=GREEN if good else RED,
                                   clip_on=False))


def render_cache(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import imageio.v3 as iio

    plt.style.use("dark_background")

    num_stock, num_var = args.ns, args.nv
    total = num_stock * num_var
    orders = [
        (True, "row-major  (solver levels 0-6)", CYAN),
        (False, "ctl-order control  (loops swapped)", RED),
    ]
    caches = [ModelCache(args.cache_lines) for _ in orders]
    walkers = [visit_order(num_stock, num_var, row_major)
               for row_major, _, _ in orders]

    # Both panels advance the same number of CELLS per frame, so the clip
    # compares equal work rather than equal time.
    cells_total = (num_stock - 2) * (num_var - 2)
    per_frame = max(1, cells_total // args.frames)

    frames = []
    states = [(None, [])] * len(orders)
    for _ in range(args.frames):
        for panel, walker in enumerate(walkers):
            for _ in range(per_frame):
                try:
                    stock_i, var_j = next(walker)
                except StopIteration:
                    break
                touched = stencil_indices(stock_i, var_j, num_stock)
                for flat in touched:
                    caches[panel].access(flat)
                caches[panel].cells += 1
                states[panel] = (var_j * num_stock + stock_i, touched)

        fig, axes = plt.subplots(1, 2, figsize=(13, 6.0), dpi=120)
        for panel, (_, title, colour) in enumerate(orders):
            centre, live = states[panel]
            draw_panel(axes[panel], num_stock, num_var,
                       caches[panel].resident_doubles(total), live, centre,
                       title, caches[panel], colour)

        caption = (f"modelled {args.cache_lines}-line LRU cache, 64 B lines   |   "
                   f"grid {num_stock}x{num_var}, flat idx = var_j*{num_stock} "
                   f"+ stock_i   |   the model simulates capacity only")
        if args.slowdown:
            traffic = (caches[1].bytes_per_cell /
                       max(1e-9, caches[0].bytes_per_cell))
            # Be upfront about the gap between traffic and time. On the
            # rangpur node the two buffers come to 16 MiB against a 30 MB L3,
            # so the extra traffic is absorbed by cache and costs much less
            # time than the byte count on its own suggests. "4.7x traffic so
            # 4.7x slower" would be the easy thing to claim here, and wrong.
            caption += (
                f"\nThe swapped order moves {traffic:.1f}x the memory, but "
                f"measured on rangpur it costs only {args.slowdown:.2f}x the "
                f"time. At this footprint the extra traffic is served by a "
                f"30 MB L3, not DRAM, so the kernel is not bandwidth-bound "
                f"and the bytes are cheap. Layout still decides who pays.")
        fig.suptitle("Same 9-point stencil, same arithmetic, two traversal "
                     "orders", color=WHITE, fontsize=14, y=0.985)
        fig.text(0.5, 0.018, caption, color=MID, fontsize=8.5, ha="center")
        fig.tight_layout(rect=(0, 0.10, 1, 0.945))
        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
        plt.close(fig)

    # Hold the final state so the closing counters are readable.
    frames.extend([frames[-1]] * args.fps * 2)
    write_movie(frames, args)
    print(f"row-major: {caches[0].bytes_per_cell:.1f} B/cell, "
          f"hit rate {100 * caches[0].hit_rate:.1f}%")
    print(f"ctl-order: {caches[1].bytes_per_cell:.1f} B/cell, "
          f"hit rate {100 * caches[1].hit_rate:.1f}%")
    print(f"traffic ratio: {caches[1].bytes_per_cell / max(1e-9, caches[0].bytes_per_cell):.1f}x "
          f"more memory moved for identical arithmetic")


# ---------------------------------------------------------------------------
# --style buffers
# ---------------------------------------------------------------------------

# What each ladder rung adds, and the data structure that appears with it.
LADDER_SCENES = [
    ("level 0", "no techniques, every weight recomputed per cell", []),
    ("level 1", "hoisting + lookup tables + CSE",
     ["row constants: half_v, rho_xi_v, half_xi2_v, mean_rev  (once per ROW)",
      "S table: stock[i] = i*ds  (once per STEP, reused by every row)"]),
    ("level 2", "strength reduction",
     ["reciprocals: inv_two_ds, inv_ds_sq, inv_two_dv, inv_dv_sq, "
      "inv_four_ds_dv",
      "S^2 table: stock_sq[i] = S*S"]),
    ("level 3", "traversal order made explicit",
     ["row base offsets: row, row_above, row_below  (once per ROW)"]),
    ("level 4", "loop splitting", ["hoisted bounds: last_stock, last_var",
                                   "v=0 row peeled into its own loop"]),
    ("level 5", "induction variable",
     ["row POINTERS: row_mid, row_above, row_below, out_row",
      "ten var_j*ns + stock_i evaluations per cell  ->  zero"]),
    ("level 6", "unrolling x4",
     ["four independent cells per loop iteration",
      "one compare-and-branch amortised over four cells"]),
]


def draw_buffers_scene(ax, phase, swap_progress):
    """The two sheets and the two pointers, with the swap animated."""
    import matplotlib.patches as patches

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    left_x, right_x = 1.2, 6.0
    for x, name, colour in ((left_x, "buffer A", CYAN), (right_x, "buffer B",
                                                         AMBER)):
        ax.add_patch(patches.FancyBboxPatch(
            (x, 1.4), 2.8, 2.6, boxstyle="round,pad=0.08",
            facecolor=DIM, edgecolor=colour, linewidth=1.6))
        ax.text(x + 1.4, 1.05, name, color=colour, ha="center", fontsize=9)
        # A hint of the grid inside each buffer.
        for row in range(6):
            for col in range(7):
                ax.add_patch(patches.Rectangle(
                    (x + 0.25 + col * 0.33, 1.7 + row * 0.36), 0.28, 0.30,
                    facecolor=colour, alpha=0.12 + 0.03 * row, edgecolor="none"))

    # The two pointer labels slide across each other during the swap.
    cur_x = left_x + 1.4 + swap_progress * (right_x - left_x)
    next_x = right_x + 1.4 - swap_progress * (right_x - left_x)
    ax.text(cur_x, 4.55, "current()", color=WHITE, ha="center", fontsize=11,
            fontweight="bold")
    ax.annotate("", xy=(cur_x, 4.1), xytext=(cur_x, 4.45),
                arrowprops={"arrowstyle": "-|>", "color": WHITE, "lw": 1.6})
    ax.text(next_x, 5.25, "next()", color=AMBER, ha="center", fontsize=11,
            fontweight="bold")
    ax.annotate("", xy=(next_x, 4.1), xytext=(next_x, 5.15),
                arrowprops={"arrowstyle": "-|>", "color": AMBER, "lw": 1.6})

    if phase == "read":
        ax.text(5.0, 0.45, "read the finished sheet, write the new one",
                color=MID, ha="center", fontsize=10)
    else:
        ax.text(5.0, 0.45,
                "std::swap exchanges two pointers: O(1), not one double copied",
                color=GREEN, ha="center", fontsize=10)


def render_buffers(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import imageio.v3 as iio

    plt.style.use("dark_background")
    frames = []

    # Scene 1: the ping-pong, three full cycles.
    for _ in range(3):
        for step in range(args.fps):
            progress = step / args.fps
            phase = "read" if progress < 0.6 else "swap"
            swap = 0.0 if phase == "read" else (progress - 0.6) / 0.4
            fig, ax = plt.subplots(figsize=(11, 5.4), dpi=120)
            draw_buffers_scene(ax, phase, swap)
            fig.suptitle("Two buffers, swapped and never copied",
                         color=WHITE, fontsize=14, y=0.95)
            fig.tight_layout()
            fig.canvas.draw()
            frames.append(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
            plt.close(fig)

    # Scene 2: the ladder's data structures accumulating.
    accumulated = []
    for name, technique, additions in LADDER_SCENES:
        accumulated.extend(additions)
        for _ in range(int(args.fps * 1.6)):
            fig, ax = plt.subplots(figsize=(11, 5.4), dpi=120)
            ax.set_xlim(0, 10)
            ax.set_ylim(0, 6)
            ax.axis("off")
            ax.text(0.3, 5.5, name, color=AMBER, fontsize=17,
                    fontweight="bold")
            ax.text(2.3, 5.55, technique, color=WHITE, fontsize=12)
            # The list grows every rung, so the line spacing shrinks to fit
            # rather than running off the bottom of the frame.
            spacing = min(0.42, 4.05 / max(1, len(accumulated)))
            size = 10.5 if spacing > 0.36 else 9.0
            for i, line in enumerate(accumulated):
                fresh = line in additions
                ax.text(0.5, 4.75 - i * spacing, "•  " + line,
                        color=AMBER if fresh else MID,
                        fontsize=size if fresh else size - 1,
                        fontfamily="monospace")
            ax.text(0.3, 0.25,
                    "each rung adds to the one below it: the ladder is "
                    "cumulative, which is what makes the ablation readable",
                    color=MID, fontsize=9)
            fig.suptitle("What the optimised kernel keeps in memory, rung by "
                         "rung", color=WHITE, fontsize=14, y=0.97)
            fig.tight_layout()
            fig.canvas.draw()
            frames.append(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
            plt.close(fig)

    frames.extend([frames[-1]] * args.fps * 2)
    write_movie(frames, args)


def write_movie(frames, args):
    import imageio.v3 as iio

    # imageio's libx264 needs even dimensions.
    height, width = frames[0].shape[:2]
    if height % 2 or width % 2:
        frames = [f[: height - height % 2, : width - width % 2] for f in frames]
    iio.imwrite(args.out, frames, fps=args.fps, codec="libx264")
    print(f"wrote {args.out} ({len(frames)} frames @ {args.fps} fps)")
    if args.gif:
        gif_path = args.out.rsplit(".", 1)[0] + ".gif"
        iio.imwrite(gif_path, frames[::2], duration=2000 / args.fps, loop=0)
        print(f"wrote {gif_path}")


def measured_slowdown(paths):
    """ctl-order's slowdown against its PAIRED reference, read from the bench
    CSVs. Only a comparison made inside one sweep means anything, because the
    sweeps differ in time-loop length and in how memory was allocated. Sweep E
    is preferred (full ladder plus both controls, huge pages, which is what a
    real run gets), then D (the dedicated controls job), then A."""
    import statistics

    buckets = {}
    for path in paths:
        try:
            handle = open(path)
        except OSError:
            continue
        with handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split(",")
                if len(fields) != 12:
                    continue
                buckets.setdefault((fields[0], fields[9]), []).append(
                    float(fields[11]))

    for sweep in ("E", "D", "A"):
        control = buckets.get((sweep, "opt-ctl-order"))
        # Compare against the fastest ladder rung measured in the same sweep.
        rungs = [statistics.median(values)
                 for (s, solver), values in buckets.items()
                 if s == sweep and solver.startswith("opt-L")]
        if control and rungs:
            return max(rungs) / statistics.median(control)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--style", choices=("cache", "buffers"),
                        default="cache")
    parser.add_argument("--out", default="results/memory_cache.mp4")
    parser.add_argument("--ns", type=int, default=128,
                        help="stock nodes in the toy grid (cache style)")
    parser.add_argument("--nv", type=int, default=48,
                        help="variance nodes in the toy grid (cache style)")
    parser.add_argument("--cache-lines", type=int, default=48,
                        help="lines in the modelled LRU cache. The default is\n                             exactly 3 rows at the default --ns, i.e.\n                             just enough to hold what the stencil needs")
    parser.add_argument("--frames", type=int, default=140)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--slowdown", type=float, default=None,
                        help="measured ctl-order slowdown, for the caption")
    parser.add_argument("--bench-csv", nargs="*", default=None,
                        help="rangpur bench CSVs; the ctl-order slowdown is "
                             "read straight out of them, so the caption cannot "
                             "drift from the measurement it quotes")
    parser.add_argument("--gif", action="store_true")
    args = parser.parse_args()

    if args.bench_csv and args.slowdown is None:
        args.slowdown = measured_slowdown(args.bench_csv)
        if args.slowdown:
            print(f"ctl-order slowdown read from the benchmark CSVs: "
                  f"{args.slowdown:.1f}x")
        else:
            print("no ctl-order rows in the given CSVs, so the caption will "
                  "leave out the measured slowdown", file=sys.stderr)

    if args.style == "cache":
        if args.ns % DOUBLES_PER_LINE:
            print(f"--ns must be a multiple of {DOUBLES_PER_LINE} so rows "
                  f"start on a cache-line boundary", file=sys.stderr)
            return 1
        render_cache(args)
    else:
        render_buffers(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
