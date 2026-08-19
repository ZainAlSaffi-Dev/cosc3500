#!/usr/bin/env python3
"""Roofline for the Heston kernel (PLAN P7): is it compute- or memory-bound?

The optimisation story rests on the claim that the kernel is memory-bound, so
layout and loop order are what matter. This script checks that claim rather
than just asserting it:

  1. Count flops per cell update, term by term, from the actual level-6
     kernel in src/solver_opt_kernels.cpp. The table gets printed so it can
     be checked against the source line by line.
  2. Count the bytes that have to cross the memory bus per cell update.
  3. Divide the two to get arithmetic intensity in flops per byte.
  4. Take the measured cell-updates/sec from the rangpur CSVs, multiply by
     the flop count for achieved GFLOP/s, and put that point under the roof.

Two ceilings get drawn:
  - a spec-sheet DRAM bandwidth roof, if --peak-bandwidth is given;
  - a streaming ceiling worked out from our own ctl-order control. That
    control walks the grid with an ns-element stride, so it uses about one
    double out of every 64-byte line it fetches, and its measured throughput
    implies a bandwidth the machine actually sustained. It is a lower bound
    on the real ceiling, but at least it is measured rather than quoted.

Usage:
    python3 scripts/roofline.py results/bench_*.csv --out results/roofline.png
    python3 scripts/roofline.py results/bench_*.csv --peak-bandwidth 68
"""

import argparse
import collections
import pathlib
import statistics
import sys

# The palette and rcParams shared by every figure in the project.
from figstyle import BLUE, GREY, DARK_GREY, RED, AMBER, TEXT
import figstyle

# The sweep-B points are still level 6, just at other grid sizes, so they
# wear a lighter shade of the level-6 blue rather than a new hue.
LIGHT_BLUE = "#94B5F6"

# Flops per interior cell update, read straight off step_level6(). Each entry
# is (expression as it appears in the kernel, flop count, note).
# Common subexpressions the compiler will certainly eliminate are counted only
# once, since `2.0 * V` appears in both V_SS and V_vv but is one multiply.
FLOP_TABLE = [
    ("2.0 * V", 1, "shared by V_SS and V_vv (CSE)"),
    ("(east - west) * inv_two_ds", 2, "V_S"),
    ("(east - 2V + west) * inv_ds_sq", 3, "V_SS"),
    ("(north - south) * inv_two_dv", 2, "V_v"),
    ("(north - 2V + south) * inv_dv_sq", 3, "V_vv"),
    ("(ne - nw - se + sw) * inv_four_ds_dv", 4, "V_Sv"),
    ("half_v * stock_sq[i] * V_SS", 2, "S-diffusion"),
    ("rho_xi_v * stock[i] * V_Sv", 2, "cross term"),
    ("half_xi2_v * V_vv", 1, "v-diffusion"),
    ("carry * stock[i] * V_S", 2, "S-drift"),
    ("mean_rev * V_v", 1, "v-drift"),
    ("rate * V", 1, "discounting"),
    ("A + B + C + D + E - F", 5, "combine the six PDE terms"),
    ("dt * (...)", 1, "explicit Euler step"),
    ("V + ...", 1, "add to the old value"),
]

# Bytes that have to cross the bus per cell update, assuming perfect
# row-major reuse. The 9-point stencil reads three rows, but walking along a
# row each new cell only pulls in one new double from each of them, so over a
# whole timestep it comes to reading every cell of cur once.
BYTE_TABLE = [
    ("read cur", 8, "each cell read once from memory; the other 8 of its 9 "
                    "stencil uses hit cache"),
    ("write next", 8, "each cell written once"),
    ("write-allocate", 8, "the store misses, so its line gets read before "
                          "being overwritten (no non-temporal stores here)"),
    ("stock[] / stock_sq[] tables", 0, "2*ns doubles, re-read per row but "
                                       "they live in L1, so about 0"),
]

# Bytes ctl-order moves per cell. Its inner loop strides ns doubles, so the
# three stencil rows land in three different cache lines and hardly any of
# the other 7 doubles in each line get used before it is evicted.
CTL_ORDER_BYTES_PER_CELL = 4 * 64  # 3 lines read + 1 line for the store


class Row:
    def __init__(self, fields):
        self.sweep = fields[0]
        self.build = fields[1]
        self.ns = int(fields[6])
        self.nv = int(fields[7])
        self.solver = fields[9]
        self.cups = float(fields[11])

    @property
    def grid(self):
        return f"{self.ns}x{self.nv}"


def read_csvs(paths):
    rows, provenance = [], []
    for path in paths:
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    text = line.lstrip("# ")
                    if not text.startswith("sweep,"):
                        provenance.append(text)
                    continue
                fields = line.split(",")
                if len(fields) == 12:
                    try:
                        rows.append(Row(fields))
                    except ValueError:
                        pass
    return rows, provenance


def median_cups(rows, solver, sweeps):
    values = [r.cups for r in rows if r.solver == solver and r.sweep in sweeps]
    return statistics.median(values) if values else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("--out", default="results/roofline.png")
    parser.add_argument("--peak-bandwidth", type=float, default=None,
                        help="spec-sheet DRAM bandwidth in GB/s (optional, "
                             "and quoted rather than measured, so label it "
                             "that way)")
    parser.add_argument("--peak-flops", type=float, default=None,
                        help="scalar double peak in GFLOP/s (optional)")
    args = parser.parse_args()

    rows, provenance = read_csvs(args.csv)
    if not rows:
        print("no benchmark rows found", file=sys.stderr)
        return 1

    # Same gate as bench_plot.py, for the same reason. This figure places a
    # measured throughput against a measured ceiling, so a laptop CSV wandering
    # into the results/bench_*.csv glob would put a point on the plot that the
    # ceiling does not apply to at all.
    hosts = sorted({text.partition("=")[2].strip().strip("'")
                    for text in provenance if text.startswith("host=")})
    if not hosts:
        print("FAIL: no '# host=' header in any CSV — cannot prove these "
              "numbers came from a compute node", file=sys.stderr)
        return 1
    strays = [h for h in hosts if not h.endswith(".compute.eait.uq.edu.au")]
    if strays:
        print(f"FAIL: {strays} are not rangpur compute nodes. Login-node and "
              f"laptop numbers are not results.", file=sys.stderr)
        return 1

    flops_per_cell = sum(count for _, count, _ in FLOP_TABLE)
    bytes_per_cell = sum(count for _, count, _ in BYTE_TABLE)
    intensity = flops_per_cell / bytes_per_cell

    print("flops per interior cell update (from step_level6):")
    for expression, count, note in FLOP_TABLE:
        print(f"  {count:2d}  {expression:<40s} {note}")
    print(f"  {'':2s}  {'':40s} ---")
    print(f"  {flops_per_cell:2d}  TOTAL\n")
    print("bytes across the bus per cell update:")
    for name, count, note in BYTE_TABLE:
        print(f"  {count:2d}  {name:<30s} {note}")
    print(f"  {bytes_per_cell:2d}  TOTAL")
    print(f"\narithmetic intensity = {flops_per_cell}/{bytes_per_cell} "
          f"= {intensity:.2f} flop/byte\n")

    # Measured points. Sweep A is the ladder at the reference footprint, and
    # sweep B repeats baseline and level 6 across four grid sizes. Every point
    # sits at the same arithmetic intensity, since it is the same kernel and
    # the same stencil, so any vertical spread between them can't be an
    # arithmetic difference. That is the argument the figure is making.
    points = []
    # Sweep E is the huge-page case, which is what a single production solve
    # actually gets, so it is the headline point when it's there.
    has_pages = any(r.sweep == "E" for r in rows)
    for solver, label, colour in (("baseline", "baseline", GREY),
                                  ("opt-L6", "opt L6", BLUE)):
        if has_pages:
            cups = median_cups(rows, solver, {"E"})
            if cups:
                points.append((f"{label}, huge pages", cups,
                               cups * flops_per_cell / 1e9, intensity, colour))
        cups = median_cups(rows, solver, {"A"})
        if cups:
            points.append((f"{label}, 4 KiB pages", cups,
                           cups * flops_per_cell / 1e9, intensity, colour))

    # Sweep B follows level 6 as the working set grows past each cache level.
    scaling = {}
    for row in rows:
        if row.sweep == "B" and row.solver == "opt-L6":
            scaling.setdefault(row.grid, []).append(row.cups)
    for grid in sorted(scaling,
                       key=lambda g: int(g.split("x")[0]) * int(g.split("x")[1])):
        ns, nv = (int(x) for x in grid.split("x"))
        mib = 2 * 8 * ns * nv / 2**20
        cups = statistics.median(scaling[grid])
        points.append((f"L6 at {grid}, {mib:.0f} MiB", cups,
                       cups * flops_per_cell / 1e9, intensity, LIGHT_BLUE))

    # Empirical streaming ceiling from the ctl-order control.
    ctl_cups = median_cups(rows, "opt-ctl-order", {"D"}) or \
        median_cups(rows, "opt-ctl-order", {"A"})
    empirical_bandwidth = None
    if ctl_cups:
        empirical_bandwidth = ctl_cups * CTL_ORDER_BYTES_PER_CELL / 1e9
        print(f"ctl-order sustained {ctl_cups / 1e6:.1f} Mcell-updates/s at "
              f"~{CTL_ORDER_BYTES_PER_CELL} B/cell")
        print(f"  => {empirical_bandwidth:.1f} GB/s of load/store traffic "
              f"sustained at this footprint.")
        print("  note: at 2048x512 the two buffers come to 16 MiB and this "
              "node's L3 is 30 MB,\n        so that traffic is served by "
              "cache, not DRAM. It is a lower bound on the\n"
              "        memory hierarchy's throughput, not a measurement of "
              "main memory bandwidth.")

    if points:
        best = max(points, key=lambda p: p[2])
        required = best[1] * bytes_per_cell / 1e9
        print(f"\n{best[0]}: {best[1] / 1e6:.1f} Mcell-updates/s "
              f"= {best[2]:.2f} GFLOP/s, demanding {required:.1f} GB/s")
        if empirical_bandwidth:
            print(f"  that is {100 * required / empirical_bandwidth:.0f}% of "
                  f"the empirically demonstrated bandwidth")
        if args.peak_bandwidth:
            print(f"  and {100 * required / args.peak_bandwidth:.0f}% of the "
                  f"quoted {args.peak_bandwidth:.0f} GB/s peak")

    # ---- the figure ----
    import matplotlib
    matplotlib.use("Agg")
    figstyle.apply()
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(9, 5.8), dpi=160)
    # The interesting window is a decade either side of the kernel's own
    # arithmetic intensity. The old four-decade axis was mostly empty canvas.
    x = np.logspace(-1.6, 1.35, 300)

    roofs = []
    if empirical_bandwidth:
        roofs.append((empirical_bandwidth, AMBER, "-",
                      f"memory roof demonstrated by ctl-order "
                      f"({empirical_bandwidth:.0f} GB/s)"))
    if args.peak_bandwidth:
        roofs.append((args.peak_bandwidth, GREY, "--",
                      f"quoted DRAM peak ({args.peak_bandwidth:.0f} GB/s)"))

    for bandwidth, colour, style, label in roofs:
        line = bandwidth * x
        if args.peak_flops:
            line = np.minimum(line, args.peak_flops)
        ax.loglog(x, line, style, color=colour, linewidth=2, label=label,
                  zorder=2)
    if args.peak_flops:
        ax.axhline(args.peak_flops, color=RED, linestyle=":", linewidth=1.6,
                   label=f"measured scalar peak "
                         f"({args.peak_flops:.1f} GFLOP/s)",
                   zorder=2)

    # Every point has the same arithmetic intensity, so they all stack up on
    # one vertical line. Rather than scattering labels around the dots and
    # hoping they miss each other, the labels form one tidy column to the
    # right, each connected to its dot by a thin leader line, with a minimum
    # vertical gap enforced in log space.
    ordered = sorted(points, key=lambda p: p[2], reverse=True)
    label_x = intensity * 2.4
    min_ratio = 1.24
    label_ys = []
    for _, _, gflops, _, _ in ordered:
        y = gflops
        if label_ys and label_ys[-1] / y < min_ratio:
            y = label_ys[-1] / min_ratio
        label_ys.append(y)
    for (label, _, gflops, ai, colour), label_y in zip(ordered, label_ys):
        ax.loglog([ai], [gflops], "o", color=colour, markersize=10, zorder=4,
                  markeredgecolor="white", markeredgewidth=1.2)
        ax.annotate(f"{gflops:.2f}   {label}",
                    xy=(ai, gflops), xytext=(label_x, label_y),
                    textcoords="data", ha="left", va="center",
                    color=TEXT, fontsize=9,
                    arrowprops={"arrowstyle": "-", "color": "#C7CDD6",
                                "lw": 0.8,
                                "shrinkA": 4, "shrinkB": 6})

    ax.axvline(intensity, color=DARK_GREY, linewidth=0.9, linestyle=":",
               zorder=1)
    # The intensity label sits at the bottom of its guide line, well away
    # from the label column and the legend.
    ax.annotate(f"this kernel: {intensity:.2f} flop per byte\n"
                f"({flops_per_cell} flops, {bytes_per_cell} bytes per cell "
                f"update)",
                xy=(intensity, 0), xycoords=ax.get_xaxis_transform(),
                xytext=(-8, 10), textcoords="offset points",
                color=DARK_GREY, fontsize=9, ha="right", va="bottom")

    ax.set_xlabel("arithmetic intensity (flop / byte)")
    ax.set_ylabel("performance (GFLOP/s)")
    # The title reports what was measured, not what I expected. The point of
    # the figure is to test the "memory-bound" claim, so a title that assumed
    # the answer would defeat the whole exercise.
    verdict = "roof not supplied"
    if points and args.peak_flops:
        best = max(points, key=lambda p: p[2])
        share = best[2] / args.peak_flops
        verdict = (f"{100 * share:.0f}% of scalar compute peak"
                   if share > 0.4 else
                   f"only {100 * share:.0f}% of scalar compute peak")
    ax.set_title(f"Roofline for the Heston stencil: {verdict}")
    ax.grid(True, which="minor", color=figstyle.GRID, linewidth=0.4)
    figstyle.finish_axes(ax)
    if roofs or args.peak_flops:
        ax.legend(fontsize=8.5, loc="lower right")

    # The five CSVs all repeat the same machine details, so show each once.
    seen, bits = set(), []
    for text in provenance:
        if text.startswith(("host=", "cpu=", "compiler=")):
            value = text.split("=", 1)[1].strip("'")
            if value not in seen:
                seen.add(value)
                bits.append(value)
    caption = "  |  ".join(bits)
    fig.text(0.01, 0.005, caption, fontsize=6.5, color=DARK_GREY)
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    print(f"\nwrote {out_path}")

    if not points:
        print("FAIL: no sweep-A rows to place on the roofline", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
