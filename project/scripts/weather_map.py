#!/usr/bin/env python3
"""Weather map (PLAN P6): results/<run>/snap_*.csv -> animated value surface.

Each snapshot CSV: nv rows (variance), ns columns (stock price) — layout
fixed by src/io.cpp dump_snapshot. The animation runs in SOLVER order:
expiry payoff -> today.

Styles and colour modes (pick by what the data's job is):
  --style flat (default)  heatmap over the (S, v) plane.
  --style surface         3-D sheet V(S, v) with a slow camera drift — the
                          "spreadsheet" itself, kink smoothing out.
  sequential (default)    dark background + inferno (black->red->yellow):
                          low values recede, heat glows. Magnitude story.
  --diverging             white background + red/blue around 0 for the
                          instability clip, where polarity IS the story.
                          NaN/inf cells render black.
  --time-value            subtract the expiry payoff, plotting V - payoff:
                          the time value is what actually diffuses — raw V
                          wastes the colour range on the linear deep-ITM
                          ramp. (Call payoff; strike via --strike.)

Usage:
    python3 scripts/weather_map.py results/smooth --out results/weather_map.mp4
    python3 scripts/weather_map.py results/smooth --time-value --out tv.mp4
    python3 scripts/weather_map.py results/smooth --style surface --out s.mp4
    python3 scripts/weather_map.py results/blowup --diverging --max-frames 90
"""

import argparse
import glob
import os
import re

import numpy as np

S_MAX_MULT = 4.0  # matches every config: S spans [0, 4*strike]


def load_snapshots(snapshot_dir: str):
    """Sorted (step, 2-D array) pairs, ascending step = solver time order."""
    paths = sorted(glob.glob(os.path.join(snapshot_dir, "snap_*.csv")))
    if not paths:
        raise SystemExit(f"no snap_*.csv files in {snapshot_dir}")
    snaps = []
    for p in paths:
        step = int(re.search(r"snap_(\d+)", p).group(1))
        snaps.append((step, np.loadtxt(p, delimiter=",")))
    return snaps


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot_dir")
    parser.add_argument("--out", default="results/weather_map.mp4")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--style", choices=("flat", "surface"), default="flat")
    parser.add_argument("--diverging", action="store_true",
                        help="two-hue scale around 0 for the instability clip")
    parser.add_argument("--time-value", action="store_true",
                        help="plot V - expiry payoff (call) instead of V")
    parser.add_argument("--strike", type=float, default=5250.0,
                        help="strike for --time-value payoff and the S axis")
    parser.add_argument("--maturity", type=float, default=0.25,
                        help="option maturity in years — converts solver "
                             "steps into days-to-expiry for the titles")
    parser.add_argument("--spot", type=float, default=5200.0,
                        help="today's stock price, marked on flat frames")
    parser.add_argument("--v0", type=float, default=0.04,
                        help="today's variance, marked on flat frames")
    parser.add_argument("--caption", default=None,
                        help="one/two-line explainer under the axes "
                             "(default: auto text per style)")
    parser.add_argument("--gif", action="store_true",
                        help="also write a .gif next to the .mp4")
    parser.add_argument("--max-frames", type=int, default=0,
                        help="render only the first N snapshots (0 = all); "
                             "the blow-up clip's tail is all-NaN frames")
    parser.add_argument("--title", default="Heston option value V(S, v)")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import imageio.v3 as iio

    # Sequential mode is film material: dark stage so the heat glows.
    # Diverging keeps a white stage — a diverging scale's midpoint is light
    # by design, and calm regions should read quiet, not glowing.
    if not args.diverging:
        plt.style.use("dark_background")

    snaps = load_snapshots(args.snapshot_dir)
    # Total step count must come from the FULL dump list — computing it
    # after --max-frames truncation would mislabel every remaining frame's
    # step fraction and days-to-expiry.
    total_steps = max(step for step, _ in snaps)
    if args.max_frames > 0:
        snaps = snaps[: args.max_frames]

    s_max = S_MAX_MULT * args.strike
    num_stock = snaps[0][1].shape[1]
    stock_axis = np.linspace(0.0, s_max, num_stock)

    value_label = "option value ($)"
    if args.time_value:
        # Time value = V - intrinsic. The deep-ITM ramp cancels out and the
        # colour range is spent entirely on what the PDE actually moves.
        payoff = np.maximum(stock_axis - args.strike, 0.0)
        snaps = [(step, vals - payoff[np.newaxis, :]) for step, vals in snaps]
        value_label = "time value V − payoff ($)"

    # One colour scale across ALL frames, or the animation "breathes" and
    # instability frames clip. Non-finite values (blow-up) are masked.
    finite_vals = np.concatenate(
        [a[np.isfinite(a)].ravel() for _, a in snaps])
    if finite_vals.size == 0:
        raise SystemExit("every cell is NaN/inf — nothing to render")
    if args.diverging:
        # Polarity story: symmetric limits around zero, neutral midpoint.
        # 99.5th percentile, not max: one saturated cell must not wash out
        # the palette for every earlier frame.
        limit = np.percentile(np.abs(finite_vals), 99.5) or 1.0
        vmin, vmax, cmap_name, bad_colour = -limit, limit, "RdBu_r", "black"
    else:
        # Magnitude story: perceptually-uniform heat ramp. 98th percentile
        # cap keeps the late-time high-v glow from compressing the early
        # near-strike structure; hotter cells simply saturate to yellow.
        vmin = min(0.0, float(finite_vals.min()))
        vmax = float(np.percentile(finite_vals, 98.0)) or 1.0
        cmap_name, bad_colour = "inferno", "magenta"
    cmap = matplotlib.colormaps[cmap_name].copy()
    cmap.set_bad(bad_colour)  # broken cells must be unmistakable

    variance_axis = np.linspace(0.0, 1.0, snaps[0][1].shape[0])
    stock_mesh, variance_mesh = np.meshgrid(stock_axis, variance_axis)

    # Markers don't think in solver steps — translate to calendar time.
    # Step n holds the sheet n*dt before expiry; the last step is "today".
    # (Assumes the final dump lands on the final step, true whenever
    # dump_every divides nt.)
    days_total = args.maturity * 365.0
    days_per_step = days_total / total_steps

    frames = []
    if args.style == "surface":
        fig = plt.figure(figsize=(8, 4.5), dpi=120)
        ax = fig.add_subplot(projection="3d")
    else:
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)

    # Static explainer for a reader who has never seen Heston: what the
    # axes mean and which direction the solver moves. Drawn once at figure
    # level, so ax.clear() per frame leaves it alone.
    caption = args.caption
    if caption is None and args.style == "surface":
        caption = ("The sheet of what-if option values over stock price and "
                   "variance. The kinked expiry payoff is known exactly;\n"
                   "the solver smooths it backwards through time until it "
                   "reaches today's sheet.")
    elif caption is None:
        caption = ("Each pixel is one what-if scenario: stock price S "
                   "(x-axis) under variance v (y-axis); colour = the "
                   "option's value there.\nThe expiry payoff is known "
                   "exactly; the solver fills the sheet backwards in time "
                   "until it reaches today.")
    fig.subplots_adjust(bottom=0.22 if args.style == "flat" else 0.14)
    fig.text(0.5, 0.02, caption, ha="center", va="bottom", fontsize=7.5,
             color="0.65" if not args.diverging else "0.35")
    for index, (step, grid_vals) in enumerate(snaps):
        ax.clear()
        masked = np.ma.masked_invalid(grid_vals)
        broken = int(np.sum(~np.isfinite(grid_vals)))
        note = f"  [{broken} broken cells]" if broken else ""
        if args.style == "surface":
            ax.plot_surface(stock_mesh, variance_mesh,
                            np.where(np.isfinite(grid_vals), grid_vals, 0.0),
                            cmap=cmap, vmin=vmin, vmax=vmax,
                            rcount=64, ccount=128, linewidth=0,
                            antialiased=False)
            ax.set_zlim(vmin, float(finite_vals.max()))
            ax.set_zlabel(value_label, fontsize=8)
            # Slow camera drift: ~25 degrees of azimuth over the whole clip.
            ax.view_init(elev=25,
                         azim=-60 + 25.0 * index / max(1, len(snaps) - 1))
        else:
            image = ax.imshow(masked, origin="lower", aspect="auto",
                              cmap=cmap, vmin=vmin, vmax=vmax,
                              extent=(0.0, s_max, 0.0, 1.0))
            if index == 0:
                fig.colorbar(image, ax=ax, label=value_label)
            # Landmarks for a reader who doesn't know the contract: the
            # strike (where the payoff kinks) and today's market point —
            # the one cell whose value becomes the quoted price.
            accent = "0.75" if not args.diverging else "0.35"
            ax.axvline(args.strike, color=accent, linestyle="--",
                       linewidth=1.0)
            ax.text(args.strike + 0.012 * s_max, 0.95,
                    f"strike K = {args.strike:g}", color=accent, fontsize=8)
            ax.plot(args.spot, args.v0, "o", markersize=8,
                    markerfacecolor="none", markeredgecolor=accent,
                    markeredgewidth=1.5)
            ax.annotate(f"today's market (S = {args.spot:g}, v = {args.v0:g})"
                        "\n= the price we quote",
                        xy=(args.spot, args.v0),
                        xytext=(args.spot + 0.10 * s_max, args.v0 + 0.16),
                        fontsize=8, color=accent,
                        arrowprops=dict(arrowstyle="-", color=accent,
                                        linewidth=0.8))
        ax.set_xlabel("stock price S", fontsize=9)
        ax.set_ylabel("variance v", fontsize=9)
        # Headline in calendar time; the solver step drops to a second,
        # smaller line for anyone cross-referencing the code.
        days_left = step * days_per_step
        ax.set_title(f"{args.title} — {days_left:.0f} of {days_total:.0f} "
                     f"days before expiry\n"
                     f"(solver step {step}/{total_steps}){note}",
                     fontsize=9)
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        frames.append(frame)
    plt.close(fig)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    # H.264 mp4 (course submission format); imageio-ffmpeg ships ffmpeg.
    iio.imwrite(args.out, frames, fps=args.fps, codec="libx264")
    print(f"wrote {args.out} ({len(frames)} frames @ {args.fps} fps)")
    if args.gif:
        gif_path = os.path.splitext(args.out)[0] + ".gif"
        iio.imwrite(gif_path, frames, duration=1000 / args.fps, loop=0)
        print(f"wrote {gif_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
