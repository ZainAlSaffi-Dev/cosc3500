#!/usr/bin/env python3
"""Weather map (PLAN P6): results/<run>/snap_*.csv -> animated value surface.

Each snapshot CSV: nv rows (variance), ns columns (stock price) — layout
fixed by src/io.cpp dump_snapshot. Frames colour option value over the
(S, v) plane; animation runs backwards in time, expiry payoff -> today.

Usage:
    python3 scripts/weather_map.py results/run1 --out results/weather_map.mp4

TODO(P6):
  - glob + sort snap_*.csv by step number (descending = solver order)
  - fixed colour scale across frames (min/max over all snapshots, else the
    animation "breathes" and instability frames clip)
  - matplotlib imshow/pcolormesh per frame, S on x, v on y, colourbar
  - annotate step number + calendar time remaining
  - ffmpeg H.264 output (course submission format), ~20 fps
  - also emit a .gif for the video editor / README
"""

import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot_dir")
    parser.add_argument("--out", default="results/weather_map.mp4")
    parser.add_argument("--fps", type=int, default=20)
    args = parser.parse_args()
    print(f"TODO(P6): render {args.snapshot_dir} -> {args.out}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
