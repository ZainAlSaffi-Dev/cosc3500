#!/usr/bin/env python3
"""Weather map (PLAN P6): results/<run>/snap_*.csv -> animated value surface.

Each snapshot CSV: nv rows (variance), ns columns (stock price) — layout
fixed by src/io.cpp dump_snapshot. Frames colour option value over the
(S, v) plane; the animation runs in SOLVER order: expiry payoff -> today.

Two colour modes (pick by what the data's job is, not by taste):
  sequential (default) — one hue light->dark for magnitude; smooth runs.
  diverging           — two hues around a neutral zero for the instability
                        clip, where the checkerboard swings +/- and polarity
                        IS the story. NaN/inf cells render black.

Usage:
    python3 scripts/weather_map.py results/run1 --out results/weather_map.mp4
    python3 scripts/weather_map.py results/blowup --diverging --fps 12
"""

import argparse
import glob
import os
import re
import sys

import numpy as np


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
    parser.add_argument("--diverging", action="store_true",
                        help="two-hue scale around 0 for the instability clip")
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

    snaps = load_snapshots(args.snapshot_dir)
    if args.max_frames > 0:
        snaps = snaps[: args.max_frames]

    # One colour scale across ALL frames, or the animation "breathes" and
    # instability frames clip. Non-finite values (blow-up) are masked black.
    finite_vals = np.concatenate(
        [a[np.isfinite(a)].ravel() for _, a in snaps])
    if finite_vals.size == 0:
        raise SystemExit("every cell is NaN/inf — nothing to render")
    if args.diverging:
        # Polarity story: symmetric limits around zero, neutral midpoint.
        # 99.5th percentile, not max: one saturated cell must not wash out
        # the palette for every earlier frame.
        limit = np.percentile(np.abs(finite_vals), 99.5) or 1.0
        vmin, vmax, cmap_name = -limit, limit, "RdBu_r"
    else:
        # Magnitude story: one hue, light -> dark.
        vmin, vmax, cmap_name = 0.0, float(finite_vals.max()), "Blues"
    cmap = matplotlib.colormaps[cmap_name].copy()
    cmap.set_bad("black")  # non-finite cells: unmistakably "broken"

    total_steps = max(step for step, _ in snaps)
    frames = []
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
    for step, grid_vals in snaps:
        ax.clear()
        masked = np.ma.masked_invalid(grid_vals)
        # origin="lower": v=0 at the bottom, matching how we drew the grid
        # on paper. extent maps pixel indices to (S, v) coordinates.
        image = ax.imshow(masked, origin="lower", aspect="auto",
                          cmap=cmap, vmin=vmin, vmax=vmax,
                          extent=(0.0, 21000.0, 0.0, 1.0))
        ax.set_xlabel("stock price S")
        ax.set_ylabel("variance v")
        broken = int(np.sum(~np.isfinite(grid_vals)))
        note = f"   [{broken} cells non-finite]" if broken else ""
        ax.set_title(f"{args.title} — step {step}/{total_steps}{note}",
                     fontsize=10)
        if len(frames) == 0:
            fig.colorbar(image, ax=ax, label="option value ($)")
        fig.canvas.draw()
        # Grab the rendered canvas as an RGB image array for the encoder.
        frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        frames.append(frame)
    plt.close(fig)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    # H.264 mp4 (course submission format). macro_block_size pads odd frame
    # sizes; imageio-ffmpeg ships its own ffmpeg binary.
    iio.imwrite(args.out, frames, fps=args.fps, codec="libx264")
    print(f"wrote {args.out} ({len(frames)} frames @ {args.fps} fps)")
    if args.gif:
        gif_path = os.path.splitext(args.out)[0] + ".gif"
        iio.imwrite(gif_path, frames, duration=1000 / args.fps, loop=0)
        print(f"wrote {gif_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
