#!/usr/bin/env python3
"""Benchmark figures (PLAN P7): results/bench_*.csv -> before/after plots.

Input: rangpur CSVs from slurm/bench_serial.sh (PLAN §3 lines + '#' header
comments carrying host/flags/date). Never plot login-node or laptop numbers
as results.

Outputs:
  - bar chart: cell-updates/sec, baseline vs opt, per grid size (median of
    reps, min/max whiskers)
  - speedup-vs-grid-size line plot
  - results/bench_summary.md table for the video script

Usage:
    python3 scripts/bench_plot.py results/bench_1234.csv --out results/

TODO(P7): parse, aggregate median/min/max per (grid, solver), render both
figures, write the summary table.
"""

import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("--out", default="results/")
    args = parser.parse_args()
    print(f"TODO(P7): plot {args.csv} -> {args.out}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
