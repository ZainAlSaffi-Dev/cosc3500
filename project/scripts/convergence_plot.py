#!/usr/bin/env python3
"""Grid-refinement convergence study (PLAN P4). Runs the sweep and plots it.

Each rung halves both spacings (ds, dv) and rescales nt to keep the explicit
stability margin fixed. dt goes like ds^2, so dt quarters per rung.

The error is measured against an EXTERNAL anchor -- the semi-analytic Heston
price by Fourier inversion, imported from monte_carlo_check.py -- rather than
against this study's own finest grid. Self-anchoring gives the finest rung no
error bar, spends a rung of compute on a reference instead of a data point,
and silently subtracts the anchor's own error from every other rung.

Every rung must be node-aligned, in THREE senses, all checked at run time by
check_alignment() rather than asserted in a comment:
  - spot 5200 lands on a stock node, so the readout is not interpolated
  - the STRIKE lands on a stock node, so the payoff kink is discretised the
    same way on every rung instead of leaving a sawtooth on the trend
  - v0 0.04 lands on a variance node, at node 2 or deeper

Usage:
    python3 scripts/convergence_plot.py --binary ./heston \
        --config config/reference.cfg --out results/convergence.png
"""

import argparse
import csv
import datetime
import math
import pathlib
import socket
import subprocess
import sys

# The Fourier anchor lives in the Monte Carlo checker. Importing it rather than
# re-deriving it means there is exactly one implementation of Heston's
# characteristic function in the repository, so the two studies cannot drift
# apart and quietly disagree about what the "true" price is.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from monte_carlo_check import (  # noqa: E402
    load_config as mc_load_config,
    semi_analytic_call,
)

# The refinement ladder. Every rung must be NODE-ALIGNED, meaning the quote
# point (spot 5200, v0 0.04) lands exactly on a grid node, because otherwise
# extract_result's bilinear blend contributes an interpolation error of its
# own and buries the discretisation trend this study is trying to measure.
#
# ds = s_max_mult*strike/(ns-1) = 21000/(ns-1), so ns-1 must divide 5200 and
# 21000 into whole numbers, which multiples of 105 do.
#
# dv = v_max/(nv-1), so nv-1 must make 0.04 land on a node. THIS DEPENDS ON
# v_max AND IT CHANGED. The original ladder used nv-1 as a multiple of 25,
# which is correct only for v_max = 1.0. P7 moved config/reference.cfg to
# v_max = 0.64 for stability, and 0.04/(0.64/25) = 1.5625 — so every rung of
# the old ladder silently went off-node. With v_max = 0.64 the requirement is
# that nv-1 be a multiple of 16, and the ladder below starts at 32 so that v0
# sits at node 2 or deeper on every rung, clear of the special v=0 row that
# the vega stencil would otherwise reach into.
#
# check_alignment() below verifies this against the config at run time rather
# than trusting this comment, because trusting the comment is precisely how
# the ladder went stale in the first place.
# THE STRIKE MUST LAND ON A NODE TOO. The payoff has a kink at the strike, and
# a kink sitting mid-cell is discretised differently on each rung. That shows up
# as a sawtooth riding on top of the convergence trend, and it is why an earlier
# version of this ladder produced an error sequence that went DOWN and then back
# UP. ds must therefore divide gcd(5200, 5250) = 50, and with s_max = 21000 the
# coarsest grid that manages it is ns = 421 (ds = 50). Hence the ladder starts
# there rather than at something cheaper.
LADDER = [
    # (ns,   nv,   nt)      ds     dv      dt        dt/dt_stable  ~cost
    (421,   33,   56000),   # 50   0.02    4.46e-6   0.50          1 s
    (841,   65,  224000),   # 25   0.01    1.12e-6   0.50          20 s
    (1681, 129,  896000),   # 12.5 0.005   2.79e-7   0.50          5 min
    (3361, 257, 3584000),   # 6.25 0.0025 6.98e-8    0.50          1.3 h
]
MATURITY_YEARS = 0.25  # matches every config; only used to report dt


def read_config(path: str) -> dict:
    """The 'key = value' lines of a .cfg, as floats where they parse."""
    values = {}
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if "=" not in line:
                continue
            key, _, value = line.partition("=")
            try:
                values[key.strip()] = float(value.strip())
            except ValueError:
                values[key.strip()] = value.strip()
    return values


def check_alignment(config_path: str, rungs) -> list:
    """Complain about any rung whose quote point misses a grid node.

    Reads the geometry out of the config instead of assuming it, so that
    changing s_max_mult or v_max in the .cfg cannot silently invalidate the
    ladder the way the P7 v_max change did.
    """
    cfg = read_config(config_path)
    strike = cfg["option.strike"]
    s_max = cfg["grid.s_max_mult"] * strike
    v_max = cfg["grid.v_max"]
    spot, v0 = cfg["market.spot"], cfg["heston.v0"]
    complaints = []
    for ns, nv, _ in rungs:
        stock_spacing = s_max / (ns - 1)
        variance_spacing = v_max / (nv - 1)
        stock_node = spot / stock_spacing
        strike_node = strike / stock_spacing
        variance_node = v0 / variance_spacing
        if abs(stock_node - round(stock_node)) > 1e-9:
            complaints.append(
                f"ns={ns}: spot {spot} sits at node {stock_node:.4f}, not on a "
                f"node (ds={stock_spacing}). ns-1 must divide {spot:g} and "
                f"{s_max:g} evenly.")
        if abs(strike_node - round(strike_node)) > 1e-9:
            complaints.append(
                f"ns={ns}: strike {strike} sits at node {strike_node:.4f}, not "
                f"on a node (ds={stock_spacing}). The payoff kink then lands "
                f"mid-cell and its discretisation error rides on top of the "
                f"convergence trend as a sawtooth.")
        if abs(variance_node - round(variance_node)) > 1e-9:
            complaints.append(
                f"nv={nv}: v0 {v0} sits at node {variance_node:.4f}, not on a "
                f"node (dv={variance_spacing}). With v_max={v_max:g}, nv-1 must "
                f"be a multiple of {round(v0 and v_max / v0) or '?'}.")
        elif round(variance_node) < 2:
            complaints.append(
                f"nv={nv}: v0 lands on node {round(variance_node)}, so the vega "
                f"stencil reaches the special v=0 transport row. Use a finer "
                f"variance axis.")
    return complaints


def run_solver(binary: str, config: str, ns: int, nv: int, nt: int) -> dict:
    """One solve through the PLAN §3 CLI. Returns the parsed CSV result."""
    cmd = [binary, "--config", config,
           "--ns", str(ns), "--nv", str(nv), "--nt", str(nt)]
    print(f"  running ns={ns} nv={nv} nt={nt} ...", file=sys.stderr, flush=True)
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    # The solver warns on stderr when the quote point is off-node. capture_output
    # swallows that stream, so it is re-raised here rather than discarded — this
    # warning was being emitted and thrown away for the whole of P7.
    if "off-node" in out.stderr:
        raise SystemExit(
            f"solver reports an off-node quote at ns={ns} nv={nv}:\n"
            f"{out.stderr.strip()}")
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
    parser.add_argument("--replot", action="store_true",
                        help="redraw the figure from the existing CSV without "
                             "re-solving. The finest rung costs minutes, so a "
                             "change to a label should not cost a re-run.")
    args = parser.parse_args()

    rungs = LADDER[: args.max_rungs]
    if len(rungs) < 3:
        print("need >= 3 rungs to show a rate", file=sys.stderr)
        return 1

    # Refuse before spending two minutes of CPU on a ladder that cannot
    # measure what it claims to measure.
    complaints = check_alignment(args.config, rungs)
    if complaints:
        print(f"FAIL: refinement ladder is not node-aligned for "
              f"{args.config}", file=sys.stderr)
        for complaint in complaints:
            print(f"  - {complaint}", file=sys.stderr)
        return 1

    # The anchor is the semi-analytic Heston price obtained by Fourier
    # inversion, which shares no code and no discretisation with the solver.
    #
    # The earlier version of this study used the FINEST RUNG as the anchor,
    # which has two problems. The finest rung has no error bar of its own, so
    # it cannot be plotted and one rung of expensive compute is spent producing
    # a reference rather than a data point. Worse, once the coarse rungs get
    # close, the differences being measured are between two numbers that are
    # both wrong by a comparable amount, and the apparent order stops meaning
    # anything. An external anchor fixes both: every rung gets a real error.
    reference_price = semi_analytic_call(mc_load_config(args.config))
    print(f"anchor: semi-analytic (Fourier) call = {reference_price:.7f}",
          file=sys.stderr)

    if args.replot:
        csv_in = pathlib.Path(args.out).with_suffix(".csv")
        results = []
        with open(csv_in) as f:
            for row in csv.DictReader(
                    line for line in f if not line.startswith("#")):
                results.append({
                    "price": float(row["price"]), "ns": int(row["ns"]),
                    "nv": int(row["nv"]), "nt": int(row["nt"]),
                    "dt": float(row["dt"]), "seconds": float(row["seconds"]),
                })
        print(f"replotting {len(results)} rungs from {csv_in}",
              file=sys.stderr)
    else:
        results = [run_solver(args.binary, args.config, *r) for r in rungs]
    for r in results:
        r["abs_err"] = abs(r["price"] - reference_price)

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Raw numbers saved next to the figure, since the video quotes them.
    csv_path = out_path.with_suffix(".csv")
    with open(csv_path, "w", newline="") as f:
        # Stamp where this ran, in the same '# key=value' style the slurm jobs
        # use. This is an ACCURACY study, not a throughput one, so it is
        # legitimately a local run and its `seconds` column is not a benchmark.
        # Saying so in the file matters because convergence.csv lands in the
        # same results/ directory as the rangpur CSVs, and an unlabelled
        # timing column sitting next to stamped ones invites being quoted as
        # if it were one.
        f.write(f"# study=convergence (accuracy, not throughput)\n")
        f.write(f"# host={socket.gethostname()}\n")
        f.write(f"# date={datetime.datetime.now().astimezone().isoformat()}\n")
        f.write("# NOTE: run locally. The `seconds` column is a local timing "
                "and is NOT a benchmark result; every quotable timing comes "
                "from slurm/bench_*.sh on a rangpur compute node.\n")
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    # Observed order in dt between neighbouring rungs. dt quarters per rung,
    # so an error ratio near 4 means order 1 in dt (log base 4 of the ratio).
    print("\nrung results:")
    for r in results:
        err = f"{r['abs_err']:.6f}" if "abs_err" in r else "(anchor)"
        print(f"  ns={r['ns']:4d} nv={r['nv']:3d} nt={r['nt']:6d} "
              f"dt={r['dt']:.3e}  price={r['price']:.6f}  err={err}")
    orders = []
    # Every rung now has a real error against the external anchor, so every
    # consecutive pair yields an order. The old version had to drop the finest
    # rung because it WAS the anchor.
    for a, b in zip(results[:-1], results[1:]):
        ratio = a["abs_err"] / b["abs_err"]
        order = math.log(ratio) / math.log(a["dt"] / b["dt"])
        orders.append(order)
        print(f"  err ratio {ratio:.2f} -> observed order in dt: {order:.2f}")

    import matplotlib
    matplotlib.use("Agg")  # file output only, no display needed
    import figstyle
    figstyle.apply()
    import matplotlib.pyplot as plt

    plotted = results
    dts = [r["dt"] for r in plotted]
    errs = [r["abs_err"] for r in plotted]

    fig, ax = plt.subplots(figsize=(8, 5.2), dpi=160)
    # There is only one series and the title names it, so no legend is needed.
    ax.loglog(dts, errs, "o-", color=figstyle.BLUE, linewidth=2.2,
              markersize=9, zorder=3,
              markeredgecolor="white", markeredgewidth=1.2)
    # Neutral dashed guide with exact slope 1, anchored at the finest point.
    guide_x = [dts[-1], dts[0]]
    guide_y = [errs[-1], errs[-1] * (dts[0] / dts[-1])]
    ax.loglog(guide_x, guide_y, "--", color=figstyle.GREY, linewidth=1.5,
              zorder=2)
    # A REFERENCE line, not a fit. Labelling it "first order in dt" invited the
    # reader to assume the data follows it; the data actually steepens then
    # flattens onto a floor, so the annotation states the measured orders.
    # The label sits at the guide's midpoint, below the line, where it cannot
    # clip the right edge the way an endpoint label did.
    mid_x = (guide_x[0] * guide_x[1]) ** 0.5
    mid_y = errs[-1] * (mid_x / dts[-1])
    ax.annotate("a slope of exactly one, for reference", xy=(mid_x, mid_y),
                xytext=(0, -26), textcoords="offset points", ha="center",
                color=figstyle.DARK_GREY, fontsize=9)
    ax.annotate(
        "the observed order in dt is "
        + " and then ".join(f"{o:.2f}" for o in orders)
        + ".\nThe error flattens onto a floor rather than following one "
        "power law.",
        # Bottom-right: the data runs bottom-left to top-right, so this corner
        # is the only one that stays clear of the points and their labels.
        xy=(0.98, 0.03), xycoords="axes fraction", color=figstyle.DARK_GREY,
        fontsize=9, va="bottom", ha="right")
    for r in plotted:  # direct labels: which grid produced each point
        ax.annotate(f"{r['ns']}x{r['nv']}", xy=(r["dt"], r["abs_err"]),
                    xytext=(-10, 8), textcoords="offset points", ha="right",
                    color=figstyle.TEXT, fontsize=9.5, fontweight=600)
    ax.set_xlabel("timestep dt (years)")
    ax.set_ylabel("error in the price, dollars\n(against the Fourier price "
                  f"{reference_price:.4f})")
    ax.set_title("The answer under grid refinement")
    ax.grid(True, which="minor", color=figstyle.GRID, linewidth=0.4)
    figstyle.finish_axes(ax)
    fig.tight_layout()
    fig.savefig(out_path)
    print(f"\nwrote {out_path} and {csv_path}")

    # Report EVERY pairwise order, not just the mean. A mean can hide a pair
    # that barely converged: 2.2 and 0.4 average to a respectable 1.3 while
    # describing a curve that is flattening onto an error floor, which is a
    # materially different claim from "first order throughout".
    mean_order = sum(orders) / len(orders)
    print(f"\npairwise orders in dt: "
          f"{', '.join(f'{o:.2f}' for o in orders)}  (mean {mean_order:.2f})")
    weak = [o for o in orders if o < 0.5]
    if weak:
        print(f"NOTE: {len(weak)} rung pair(s) show an order below 0.5, which "
              f"means the error has stopped falling with the grid and is "
              f"sitting on a floor. Say so rather than quoting the mean.",
              file=sys.stderr)
    if not 0.3 <= mean_order <= 2.5:
        print(f"WARNING: mean observed order {mean_order:.2f} is implausible, "
              f"check this before using the figure", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
