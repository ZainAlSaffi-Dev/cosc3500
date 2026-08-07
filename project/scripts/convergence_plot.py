#!/usr/bin/env python3
"""Grid-refinement convergence study (PLAN P4) — runs the sweep AND plots it.

Doubles resolution stepwise via CLI overrides (--ns/--nv/--nt), collects the
PLAN §3 CSV line per run, plots |price - finest| vs resolution on log-log
axes. Expected: first order in dt. The figure is a Milestone 1 deliverable.

Usage:
    python3 scripts/convergence_plot.py --binary ./heston \
        --config config/reference.cfg --out results/convergence.png

TODO(P4):
  - ladder e.g. (256,64,500) -> (512,128,1000) -> ... -> (4096,1024,8000)
    keeping dt/ds^2-style ratio stable (explicit scheme constraint)
  - subprocess each run, parse the CSV line (price = field 0)
  - reference = finest grid; plot abs diff vs nt, annotate slope
  - dump raw numbers to results/convergence.csv alongside the png
"""

import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="./heston")
    parser.add_argument("--config", default="config/reference.cfg")
    parser.add_argument("--out", default="results/convergence.png")
    args = parser.parse_args()
    print(f"TODO(P4): sweep {args.binary} -> {args.out}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
