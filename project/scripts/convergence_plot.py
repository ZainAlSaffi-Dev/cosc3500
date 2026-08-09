#!/usr/bin/env python3
"""Grid-refinement convergence study (PLAN P4) — runs the sweep AND plots it.

Halves both spacings per rung (ds, dv) and rescales nt to hold the explicit
stability margin fixed (dt scales with ds^2, so dt quarters per rung). The
error |price_k - price_finest| should then shrink ~4x per rung: second order
in the spacings == FIRST ORDER IN DT, the rate the video quotes.

Every rung is node-ALIGNED: ns/nv are chosen so spot 5200 and v0 0.04 land
exactly on grid nodes (extract_result reads the nearest cell — a half-cell
snap would bury the convergence trend under sawtooth noise).

Usage:
    python3 scripts/convergence_plot.py --binary ./heston \
        --config config/reference.cfg --out results/convergence.png
"""

import argparse
import csv
import math
import pathlib
import subprocess
import sys

# The refinement ladder. Spacings: ds = 21000/(ns-1), dv = 1/(nv-1).
# ns-1 must divide both 5200 and 21000 into integers (multiples of 105 work:
# spacing 200, 100, 50, 25); nv-1 must make 0.04 a node (multiples of 25).
# nt keeps dt ~0.79x the explicit stability bound at each rung.
LADDER = [
    # (ns,  nv,  nt)      ds     dv     dt
    (106,  26,   3500),   # 200  0.04   7.14e-5
    (211,  51,  14000),   # 100  0.02   1.79e-5
    (421, 101,  56000),   #  50  0.01   4.46e-6
    (841, 201, 224000),   #  25  0.005  1.12e-6
]
MATURITY_YEARS = 0.25  # matches every config; only used to report dt


def run_solver(binary: str, config: str, ns: int, nv: int, nt: int) -> dict:
    """One solve via the PLAN §3 CLI; returns the parsed CSV result line."""
    cmd = [binary, "--config", config,
           "--ns", str(ns), "--nv", str(nv), "--nt", str(nt)]
    print(f"  running ns={ns} nv={nv} nt={nt} ...", file=sys.stderr, flush=True)
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    fields = out.stdout.strip().splitlines()[-1].split(",")
    return {
        "price": float(fields[0]),
        "ns": ns, "nv": nv, "nt": nt,
        "dt": MATURITY_YEARS / nt,
        "seconds": float(fields[8]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="./heston")
    parser.add_argument("--config", default="config/reference.cfg")
    parser.add_argument("--out", default="results/convergence.png")
    parser.add_argument("--max-rungs", type=int, default=len(LADDER),
                        help="use fewer rungs for a quick check")
    args = parser.parse_args()

    rungs = LADDER[: args.max_rungs]
    if len(rungs) < 3:
        print("need >= 3 rungs to show a rate", file=sys.stderr)
        return 1

    results = [run_solver(args.binary, args.config, *r) for r in rungs]
    finest = results[-1]

    # Error proxy: distance from the finest rung's answer. The finest rung
    # itself has no error bar, so it anchors the comparison and is not plotted.
    for r in results[:-1]:
        r["abs_err"] = abs(r["price"] - finest["price"])

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Raw numbers next to the figure — the video quotes these directly.
    csv_path = out_path.with_suffix(".csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    # Observed order in dt between successive rungs: dt quarters per rung,
    # so err ratio ~4 <=> order 1 in dt (log base 4 of the ratio).
    print("\nrung results:")
    for r in results:
        err = f"{r['abs_err']:.6f}" if "abs_err" in r else "(anchor)"
        print(f"  ns={r['ns']:4d} nv={r['nv']:3d} nt={r['nt']:6d} "
              f"dt={r['dt']:.3e}  price={r['price']:.6f}  err={err}")
    orders = []
    for a, b in zip(results[:-2], results[1:-1]):
        ratio = a["abs_err"] / b["abs_err"]
        order = math.log(ratio) / math.log(a["dt"] / b["dt"])
        orders.append(order)
        print(f"  err ratio {ratio:.2f} -> observed order in dt: {order:.2f}")

    import matplotlib
    matplotlib.use("Agg")  # file output only, no display needed
    import matplotlib.pyplot as plt

    plotted = results[:-1]
    dts = [r["dt"] for r in plotted]
    errs = [r["abs_err"] for r in plotted]

    fig, ax = plt.subplots(figsize=(7, 5), dpi=160)
    # Single series: the title names it, so no legend box is needed.
    ax.loglog(dts, errs, "o-", color="#2563EB", linewidth=2, markersize=8,
              zorder=3)
    # Neutral dashed guide with exact slope 1, anchored at the finest point.
    guide_x = [dts[-1], dts[0]]
    guide_y = [errs[-1], errs[-1] * (dts[0] / dts[-1])]
    ax.loglog(guide_x, guide_y, "--", color="#9CA3AF", linewidth=1.5, zorder=2)
    ax.annotate("slope 1 (first order in dt)", xy=(guide_x[1], guide_y[1]),
                xytext=(8, -14), textcoords="offset points",
                color="#6B7280", fontsize=9)
    for r in plotted:  # direct labels: which grid produced each point
        ax.annotate(f"{r['ns']}×{r['nv']}", xy=(r["dt"], r["abs_err"]),
                    xytext=(8, 6), textcoords="offset points",
                    color="#374151", fontsize=9)
    ax.set_xlabel("timestep dt (years)")
    ax.set_ylabel(f"|price − finest| (finest = {finest['ns']}×{finest['nv']}, "
                  f"nt={finest['nt']})")
    ax.set_title("Grid-refinement convergence — Heston explicit FD, "
                 "reference call")
    ax.grid(True, which="both", color="#E5E7EB", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path)
    print(f"\nwrote {out_path} and {csv_path}")

    mean_order = sum(orders) / len(orders)
    if not 0.7 <= mean_order <= 1.5:
        print(f"WARNING: observed order {mean_order:.2f} far from 1 — "
              f"investigate before using the figure", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
