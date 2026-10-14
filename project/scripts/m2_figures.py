# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""Milestone 2 benchmark figures for the deck, one PNG per figure.

Reads the bench_m2_*.csv files from rangpur and draws each figure at the size
it has on its slide, so the text keeps its point size. Refuses any CSV whose
host is not a rangpur compute node, as in Milestone 1.

    .venv/bin/python scripts/m2_figures.py results/bench_m2_*.csv \
        --test-log results/logs/slurm-652647.out --peak results/peak_651302.txt --out results/m2

The ladder, scaling, size, roofline, step breakdown, placement, hybrid
layouts, correctness, -O levels and profile figures come from the main
sweeps. The compiler (E1), split load (E2 and E2b) and overlap (E9) figures
are drawn when their CSVs are passed too.
"""

import argparse
import collections
import math
import pathlib
import re
import statistics
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import figstyle
from figstyle import BLUE, GREY, LIGHT_GREY, DARK_GREY, RED, AMBER, ORANGE, TEXT

COMPUTE_NODE_SUFFIX = ".compute.eait.uq.edu.au"
COLUMNS = ["sweep", "label", "price", "delta", "gamma", "vega", "ns", "nv", "nt",
           "solver", "seconds", "cups", "threads", "ranks", "exchange_seconds",
           "boundary_seconds"]
# Cache sizes on the r730 nodes, from lscpu in the CSV headers
L2_MIB = 0.25
L3_MIB = 30.0

# Figure sizes in inches, matched to the slots in presentation/deck_m2.tex.
# WIDE fills the width of a slide, HALF fills the figure column of a slide with text beside it
WIDE = (6.6, 2.6)
HALF = (4.0, 2.6)
# Point sizes for numbers on bars, notes inside a chart and the provenance line
VALUE_FONT = 9
NOTE_FONT = 8.5
CAPTION_FONT = 6
DPI = 250


def read_csvs(paths):
    """Every row as a dict, plus the provenance lines of each file."""
    rows = []
    provenance = {}
    for path in paths:
        lines = pathlib.Path(path).read_text().splitlines()
        provenance[path] = [line.lstrip("# ") for line in lines
                            if line.startswith("#") and not line.startswith("# sweep,")]
        for line in lines:
            if not line or line.startswith("#"):
                continue
            fields = line.split(",")
            if len(fields) != len(COLUMNS):
                print(f"WARNING: skipped malformed line in {path}: {line[:60]}", file=sys.stderr)
                continue
            row = dict(zip(COLUMNS, fields))
            for key in ("ns", "nv", "nt", "threads", "ranks"):
                row[key] = int(row[key])
            for key in ("seconds", "cups", "exchange_seconds", "boundary_seconds"):
                row[key] = float(row[key])
            row["file"] = path
            rows.append(row)
    return rows, provenance


def check_provenance(provenance):
    """Returns the files whose host is not a rangpur compute node."""
    bad = []
    for path, lines in provenance.items():
        hosts = [line.split("=", 1)[1] for line in lines if line.startswith("host=")]
        if not hosts or not hosts[0].endswith(COMPUTE_NODE_SUFFIX):
            bad.append(path)
    return bad


def fact(provenance, path, key):
    for line in provenance[path]:
        if line.startswith(key + "="):
            return line.split("=", 1)[1]
    return "?"


def caption(provenance, paths):
    """One footer naming the node, CPU, compiler, flags, page regime and jobs, kept short."""
    first = paths[0]
    jobs = ", ".join(sorted({fact(provenance, p, "jobid") for p in paths}))
    nodes = ", ".join(sorted({fact(provenance, p, "nodes") for p in paths}))
    cpu = fact(provenance, first, "cpu")
    match = re.search(r"Xeon\(R\) CPU (\S+ \S+)", cpu)
    if match:
        cpu = "Xeon " + match.group(1)
    compiler = fact(provenance, first, "compiler")
    match = re.search(r"(\d+\.\d+\.\d+)", compiler)
    if match:
        compiler = "g++ " + match.group(1)
    level = re.search(r"-O\d", fact(provenance, first, "flags"))
    avx = fact(provenance, first, "avx_flags")
    avx_note = "-mavx2 -mfma"
    if "no-avx256-split" in avx:
        avx_note += ", loads not split"
    return (f"rangpur {nodes}, {cpu}, {compiler} {level.group(0) if level else ''}, AVX2 row {avx_note}, "
            f"huge pages, median of {fact(provenance, first, 'reps')} reps, jobs {jobs}")


def medians(rows, sweep):
    """Median, min and max rate per label, and one example row per label."""
    groups = collections.defaultdict(list)
    for row in rows:
        if row["sweep"] == sweep:
            groups[row["label"]].append(row)
    out = {}
    for label, group in groups.items():
        rates = [r["cups"] for r in group]
        out[label] = {
            "median": statistics.median(rates),
            "min": min(rates),
            "max": max(rates),
            "row": group[0],
            "rows": group,
            "step_ms": statistics.median([1000 * r["seconds"] / r["nt"] for r in group]),
            "boundary_ms": statistics.median([1000 * r["boundary_seconds"] / r["nt"] for r in group]),
            "exchange_ms": statistics.median([1000 * r["exchange_seconds"] / r["nt"] for r in group]),
        }
    return out


def save(fig, out_path, footer, top=1.0):
    """Lays the figure out above its footer and writes the PNG. top leaves room for a legend above the axes."""
    width = fig.get_figwidth()
    # A character of DejaVu Sans is about 0.6 of the point size wide
    lines = textwrap.wrap(footer, width=int(width * 72 / (CAPTION_FONT * 0.6))) if footer else []
    bottom = (len(lines) * CAPTION_FONT * 1.35 + 4) / 72 / fig.get_figheight()
    fig.tight_layout(rect=(0, bottom, 1, top), pad=0.4)
    if lines:
        fig.text(0.005, 0.012, "\n".join(lines), fontsize=CAPTION_FONT, color=DARK_GREY, va="bottom")
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    print(f"wrote {out_path}")


def ratio_text(value):
    """1x, 1.35x, 2.8x, 101x, with fewer digits as the number grows."""
    if abs(value - 1.0) < 1e-9:
        return "1x"
    if value >= 10:
        return f"{value:.0f}x"
    return f"{value:.2f}x" if value < 2 else f"{value:.1f}x"


def gain_text(value):
    """x2.81, x15.1, the gain of one rung over the rung below."""
    return f"x{value:.1f}" if value >= 10 else f"x{value:.2f}"


def share_inside(ax, x, height, share, width_note=None):
    """Writes a share in white inside the top of a bar, when the bar is tall enough to hold it."""
    if share >= 0.5:
        ax.text(x, height * 0.93, f"{share:.0f}%", ha="center", va="top", fontsize=VALUE_FONT - 0.5,
                color="white", weight=600)


def amdahl_fit(threads, speedups):
    """Least squares serial fraction f in 1/S = f + (1 - f)/p, which is linear in f."""
    top = 0.0
    bottom = 0.0
    for p, s in zip(threads, speedups):
        x = 1.0 - 1.0 / p
        y = 1.0 / s - 1.0 / p
        top += x * y
        bottom += x * x
    return top / bottom


def figure_ladder(ladder, provenance, paths, out_path):
    """F1, the five rungs on one node, each with its gain over the rung below."""
    order = [("baseline", "baseline"), ("serial", "serial"), ("simd", "simd"),
             ("omp24", "omp\n24 threads"), ("mpi2x12", "mpi\n2 x 12")]
    order = [(key, name) for key, name in order if key in ladder]
    fig, ax = plt.subplots(figsize=HALF)
    base = ladder["baseline"]["median"]
    floor = 20
    previous = None
    for index, (key, name) in enumerate(order):
        value = ladder[key]["median"] / 1e6
        ax.bar(index, value, color=GREY if key == "baseline" else BLUE, width=0.62)
        ax.text(index, value * 1.15, f"{value:,.0f}\n{ratio_text(ladder[key]['median'] / base)}", ha="center",
                va="bottom", fontsize=VALUE_FONT, color=TEXT, linespacing=1.1)
        if previous is not None:
            # The gain over the rung below, written inside the bar halfway up on the log axis
            ax.text(index, (floor * value) ** 0.5, gain_text(ladder[key]["median"] / previous), ha="center",
                    va="center", fontsize=VALUE_FONT - 0.5, color="white", weight=600)
        previous = ladder[key]["median"]
    ax.set_yscale("log")
    ax.set_ylim(floor, max(v["median"] for v in ladder.values()) / 1e6 * 12)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([name for _, name in order])
    ax.set_ylabel("Mcell-updates/s, log scale")
    ax.grid(axis="x", visible=False)
    figstyle.finish_axes(ax)
    save(fig, out_path, caption(provenance, paths))


def figure_tried(e3, e4, e5, e7, e8, e11, provenance, paths, out_path):
    """F1b, every variant that was tried against the version it would have replaced."""
    choices = []
    # E5 runs in the same job as E4, so like every other bar it is set against that job's omp rung
    # (RESULTS M2.6). Its own static row, the same solver with OMP_SCHEDULE=static, is the fallback
    reference = e4.get("omp_2048") or e5.get("static")
    if reference:
        for key, name in [("dynamic_1", "dynamic, chunk 1"), ("dynamic_8", "dynamic, chunk 8"), ("guided", "guided")]:
            if key in e5:
                choices.append((name, e5[key]["median"] / reference["median"], "E5"))
    if "persistent_2048" in e4 and "omp_2048" in e4:
        choices.append(("one region per solve", e4["persistent_2048"]["median"] / e4["omp_2048"]["median"], "E4"))
    if "edges_2048" in e8 and "omp_2048" in e8:
        choices.append(("boundaries in parallel", e8["edges_2048"]["median"] / e8["omp_2048"]["median"], "E8"))
    memory = []
    if "default_2048" in e11 and "pinned_2048" in e11:
        memory.append(("no huge pages", e11["default_2048"]["median"] / e11["pinned_2048"]["median"], "E11"))
    for size, grid in (("2048x512", "2048"), ("4096x1024", "4096")):
        if f"interleave_24_{size}" in e7 and f"default_24_{size}" in e7:
            memory.append((f"interleaved, {grid}", e7[f"interleave_24_{size}"]["median"] / e7[f"default_24_{size}"]["median"], "E7"))
    for size in ("2048", "4096"):
        if f"float_24_{size}" in e3 and f"omp_24_{size}" in e3:
            memory.append((f"float, {size}", e3[f"float_24_{size}"]["median"] / e3[f"omp_24_{size}"]["median"], "E3"))
    fig, axes = plt.subplots(1, 2, figsize=WIDE)
    for ax, entries, title in ((axes[0], choices, "OpenMP choices"), (axes[1], memory, "Memory, precision")):
        widest = max([ratio for _, ratio, _ in entries] + [1.3])
        for index, (name, ratio, experiment) in enumerate(entries):
            # Faster but rejected for accuracy is amber, slower is red, within 2 percent is grey
            if experiment == "E3":
                colour = AMBER
            elif ratio < 0.98:
                colour = RED
            elif ratio > 1.02:
                colour = BLUE
            else:
                colour = LIGHT_GREY
            ax.barh(index, ratio - 1.0, left=1.0, color=colour, height=0.62)
            if ratio > 2.0:
                ax.text(ratio - 0.06, index, f"{ratio:.2f}x", ha="right", va="center", fontsize=VALUE_FONT, color="white", weight=600)
            else:
                ax.text(max(ratio, 1.0) + 0.03 * widest, index, f"{ratio:.2f}x", va="center", fontsize=VALUE_FONT)
        ax.axvline(1.0, color=TEXT, lw=1)
        ax.set_yticks(range(len(entries)))
        ax.set_yticklabels([f"{name} ({experiment})" for name, _, experiment in entries], fontsize=9)
        ax.invert_yaxis()
        ax.set_xlim(0.4, widest + 0.35 * widest)
        ax.set_xlabel("rate against omp")
        ax.set_title(title, loc="left")
        ax.grid(axis="y", visible=False)
    figstyle.finish_axes(*axes)
    save(fig, out_path, caption(provenance, paths) + ", 24 threads")


def figure_strong(strong, provenance, paths, out_path):
    """F2, speedup and efficiency against threads with the Amdahl fit."""
    fig, (ax, eff) = plt.subplots(1, 2, figsize=WIDE)
    fits = {}
    for solver, colour, name in [("omp", BLUE, "omp, AVX2 row"), ("omp_scalar", GREY, "omp_scalar, scalar row")]:
        points = sorted((v["row"]["threads"], v["median"]) for k, v in strong.items() if v["row"]["solver"] == solver)
        if not points:
            continue
        one = points[0][1]
        threads = [p for p, _ in points]
        speedups = [rate / one for _, rate in points]
        f = amdahl_fit(threads[1:], speedups[1:])
        fits[solver] = f
        ax.plot(threads, speedups, "o-", color=colour, label=f"{solver}, f = {f:.3f}")
        model = [1.0 / (f + (1.0 - f) / p) for p in threads]
        ax.plot(threads, model, "--", color=colour, lw=1, alpha=0.7)
        ax.text(threads[-1] + 0.6, speedups[-1], f"{speedups[-1]:.1f}x", va="center", fontsize=VALUE_FONT, color=colour)
        eff.plot(threads, [s / p for s, p in zip(speedups, threads)], "o-", color=colour, label=name)
    top = max(max(v["row"]["threads"] for v in strong.values()), 1)
    ax.plot([1, top], [1, top], color=LIGHT_GREY, lw=1, label="ideal")
    ax.set_xlim(0, top + 4)
    ax.set_xlabel("OpenMP threads")
    ax.set_ylabel("speedup over 1 thread")
    ax.set_title("Speedup", loc="left")
    # Every point sits below the ideal line, so the upper left is empty
    ax.legend(loc="upper left", fontsize=8.5)
    eff.set_xlabel("OpenMP threads")
    eff.set_ylabel("speedup / threads")
    eff.set_ylim(0, 1.15)
    eff.set_title("Parallel efficiency", loc="left")
    eff.legend(loc="lower left")
    figstyle.finish_axes(ax, eff)
    save(fig, out_path, caption(provenance, paths))
    return fits


def figure_weak(weak, provenance, paths, out_path):
    """F3, time per step with fixed work per thread, and scaled speedup against Gustafson."""
    points = sorted((v["row"]["threads"], v) for v in weak.values())
    threads = [p for p, _ in points]
    step_ms = [v["step_ms"] for _, v in points]
    fig, (ax, sp) = plt.subplots(1, 2, figsize=WIDE)
    ax.plot(threads, step_ms, "o-", color=BLUE)
    ax.set_ylim(0, max(step_ms) * 1.3)
    ax.set_xlabel("OpenMP threads, 21 interior rows each")
    ax.set_ylabel("time per step (ms)")
    ax.set_title("Time per step", loc="left")
    # Scaled speedup is p times the one thread time over the p thread time
    scaled = [p * step_ms[0] / t for p, t in zip(threads, step_ms)]
    sp.plot(threads, scaled, "o-", color=BLUE, label="measured")
    sp.plot(threads, threads, color=LIGHT_GREY, lw=1, label="ideal")
    # Gustafson's f is the serial share of the parallel run, measured on the most threads
    widest = points[-1][1]
    serial_fraction = widest["boundary_ms"] / widest["step_ms"]
    gustafson = [p - serial_fraction * (p - 1) for p in threads]
    sp.plot(threads, gustafson, "--", color=GREY, label=f"Gustafson, f = {serial_fraction:.3f}")
    sp.text(threads[-1] + 0.6, scaled[-1], f"{scaled[-1]:.1f}", va="center", fontsize=VALUE_FONT, color=BLUE)
    sp.set_xlim(0, threads[-1] + 4)
    sp.set_xlabel("OpenMP threads")
    sp.set_ylabel("scaled speedup")
    sp.set_title("Scaled speedup", loc="left")
    sp.legend(loc="upper left")
    figstyle.finish_axes(ax, sp)
    save(fig, out_path, caption(provenance, paths))
    return serial_fraction, scaled[-1]


def figure_size(size, provenance, paths, out_path):
    """F4, throughput against working set with the cache bands, and omp over serial at each size."""
    fig, ax = plt.subplots(figsize=WIDE)
    lines = {}
    # The one grid whose two buffers sit exactly 2^24 bytes apart. Every rung runs slow there
    # for an address reason (RESULTS M2.8), so those points are drawn hollow and left off the lines
    odd_grid = "2896x724"
    for prefix, colour, name in [("serial", GREY, "serial"), ("simd", BLUE, "simd"), ("omp24", ORANGE, "omp, 24 threads")]:
        points = []
        for label, v in size.items():
            if label.startswith(prefix + "_"):
                row = v["row"]
                mib = 2 * row["ns"] * row["nv"] * 8 / 2**20
                points.append((mib, v["median"] / 1e6, f"{row['ns']}x{row['nv']}"))
        points.sort()
        lines[prefix] = points
        normal = [p for p in points if p[2] != odd_grid]
        ax.plot([p[0] for p in normal], [p[1] for p in normal], "o-", color=colour, label=name)
        for mib, rate, grid in points:
            if grid == odd_grid:
                ax.plot(mib, rate, "o", mfc="white", mec=colour, mew=1.5)
    if lines.get("serial") and lines.get("omp24"):
        serial = {grid: rate for _, rate, grid in lines["serial"]}
        for mib, rate, grid in lines["omp24"]:
            if grid != odd_grid and grid in serial:
                ax.text(mib, rate * 1.25, f"{rate / serial[grid]:.0f}x", ha="center", fontsize=VALUE_FONT, color=ORANGE)
    slow = [(rate, mib) for prefix in ("serial", "simd") for mib, rate, grid in lines.get(prefix, []) if grid == odd_grid]
    if slow:
        rate, mib = min(slow)
        ax.annotate("hollow, buffers exactly\n2^24 bytes apart", (mib, rate), xytext=(mib * 2.4, rate * 1.1),
                    fontsize=NOTE_FONT, color=DARK_GREY, va="center",
                    arrowprops=dict(arrowstyle="->", color=DARK_GREY, lw=0.8))
    bands = [(0.3, L2_MIB, ""), (L2_MIB, L3_MIB, "fits one L3"), (L3_MIB, 2 * L3_MIB, "fits both L3s"), (2 * L3_MIB, 400, "main memory")]
    shades = ["#F3F4F6", "#E5E7EB", "#F3F4F6", "#E5E7EB"]
    for (low, high, name), shade in zip(bands, shades):
        ax.axvspan(low, high, color=shade, zorder=0)
        if name:
            ax.text((max(low, 0.5) * high) ** 0.5, 0.97, name, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=NOTE_FONT, color=DARK_GREY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 300)
    ax.set_ylim(30, 30000)
    ax.set_xlabel("working set, both buffers (MiB)")
    ax.set_ylabel("Mcell-updates/s")
    ax.legend(loc="lower left", ncol=3)
    ax.grid(axis="x", visible=False)
    figstyle.finish_axes(ax)
    save(fig, out_path, caption(provenance, paths))


def figure_breakdown(strong, mpi, provenance, paths, out_path):
    """F6, each step split into the parallel update, the serial boundaries and the edge row exchange.

    Drawn as shares of the step, because at 24 threads the serial part is too thin to see on a time axis.
    """
    bars = []
    for threads in (1, 12, 24):
        key = f"omp_{threads}"
        if key in strong:
            v = strong[key]
            bars.append((f"{threads} thread{'s' if threads > 1 else ''}", v["step_ms"] - v["boundary_ms"], v["boundary_ms"], 0.0))
    for key, name in [("1node_1x24_2048x512", "1 node\n1 x 24"), ("1node_2x12_2048x512", "1 node\n2 x 12"),
                      ("2node_1x24_2048x512", "2 nodes\n1 x 24"), ("2node_2x12_2048x512", "2 nodes\n2 x 12")]:
        if key in mpi:
            v = mpi[key]
            bars.append(("mpi " + name, v["step_ms"] - v["boundary_ms"] - v["exchange_ms"], v["boundary_ms"], v["exchange_ms"]))
    fig, axes = plt.subplots(1, 2, figsize=WIDE, gridspec_kw={"width_ratios": [3, 4]}, sharey=True)
    panels = [[b for b in bars if not b[0].startswith("mpi")], [b for b in bars if b[0].startswith("mpi")]]
    for ax, panel, title in zip(axes, panels, ("OpenMP, one process", "MPI layouts")):
        labels = []
        for index, (name, update, boundary, exchange) in enumerate(panel):
            total = update + boundary + exchange
            parts = [100 * update / total, 100 * boundary / total, 100 * exchange / total]
            ax.bar(index, parts[0], color=BLUE, width=0.62)
            ax.bar(index, parts[1], bottom=parts[0], color=AMBER, width=0.62)
            ax.bar(index, parts[2], bottom=parts[0] + parts[1], color=RED, width=0.62)
            ax.text(index, 102, f"{parts[1] + parts[2]:.0f}%", ha="center", va="bottom", fontsize=VALUE_FONT, weight=600)
            labels.append(f"{name.replace('mpi ', '')}\n{total:.3f} ms")
        ax.set_xticks(range(len(panel)))
        ax.set_xticklabels(labels, fontsize=8.5)
        ax.set_ylim(0, 116)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_title(title, loc="left")
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("share of each step (%)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (BLUE, AMBER, RED)]
    fig.legend(handles, ["interior rows, parallel", "boundaries, serial", "edge row exchange"],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0), fontsize=9)
    figstyle.finish_axes(*axes)
    save(fig, out_path, caption(provenance, paths) + ", 2048 x 512, the number above each bar is the share that is not parallel",
         top=0.9)


def figure_placement(e6, e7, mpi, provenance, paths, out_path):
    """F7, which socket the threads run on and which socket's memory holds the pages."""
    fig, (left, right) = plt.subplots(1, 2, figsize=WIDE, gridspec_kw={"width_ratios": [1, 1.25]})
    # Left, 12 threads. spread puts them all on socket 0, close alternates sockets on this node
    groups = [("2048x512", "2048 x 512\n16 MiB"), ("4096x1024", "4096 x 1024\n64 MiB")]
    for index, (size, name) in enumerate(groups):
        one = e6.get(f"spread_12_{size}")
        both = e6.get(f"close_12_{size}")
        if not one or not both:
            continue
        for offset, value, colour, label in ((-0.2, one["median"], GREY, "one socket"), (0.2, both["median"], BLUE, "both sockets")):
            left.bar(index + offset, value / 1e6, width=0.38, color=colour, label=label if index == 0 else None)
            left.text(index + offset, value / 1e6 + 80, f"{value / 1e6:,.0f}", ha="center", va="bottom", fontsize=VALUE_FONT)
    left.set_xticks(range(len(groups)))
    left.set_xticklabels([name for _, name in groups])
    left.set_ylabel("Mcell-updates/s")
    # The headroom keeps the legend clear of the value labels
    left.set_ylim(0, 6600)
    left.set_title("12 threads, one socket or both", loc="left")
    left.legend(loc="upper right")
    # Right, 24 threads on the big grid. Each pair comes from one job, so it compares like with like
    pairs = [(e7.get("default_24_4096x1024"), e7.get("interleave_24_4096x1024"),
              "pages on\nsocket 0", "pages\ninterleaved", "numactl"),
             (mpi.get("1node_1x24_4096x1024"), mpi.get("1node_2x12_4096x1024"),
              "1 process\nof 24", "1 process\nper socket", "MPI")]
    ticks = []
    names = []
    for index, (before, after, before_name, after_name, how) in enumerate(pairs):
        if not before or not after:
            continue
        base = index * 2.4
        for offset, v, colour, name in ((0, before, GREY, before_name), (1, after, BLUE, after_name)):
            right.bar(base + offset, v["median"] / 1e6, width=0.8, color=colour)
            right.text(base + offset, v["median"] / 1e6 + 80, f"{v['median'] / 1e6:,.0f}", ha="center", va="bottom", fontsize=VALUE_FONT)
            ticks.append(base + offset)
            names.append(name)
        right.text(base + 0.5, after["median"] / 1e6 + 900, f"{how}, {after['median'] / before['median']:.2f}x",
                   ha="center", fontsize=VALUE_FONT, weight=600, color=BLUE)
    right.set_xticks(ticks)
    right.set_xticklabels(names, fontsize=8.5)
    right.set_ylim(0, 8200)
    right.set_title("24 threads, 4096 x 1024, pages", loc="left")
    for ax in (left, right):
        ax.grid(axis="x", visible=False)
    figstyle.finish_axes(left, right)
    save(fig, out_path, caption(provenance, paths))


def figure_layouts(mpi, provenance, paths, out_path):
    """F8, hybrid layouts on one and two nodes, with the edge row exchange share of each step."""
    layouts = [("1node_1x24", "1 node\n1 x 24"), ("1node_2x12", "1 node\n2 x 12"),
               ("2node_1x24", "2 nodes\n1 x 24"), ("2node_2x12", "2 nodes\n2 x 12")]
    fig, axes = plt.subplots(1, 2, figsize=WIDE)
    for ax, (size, grid) in zip(axes, (("2048x512", "2048 x 512"), ("4096x1024", "4096 x 1024"))):
        for index, (key, name) in enumerate(layouts):
            v = mpi.get(f"{key}_{size}")
            if not v:
                continue
            value = v["median"] / 1e6
            ax.bar(index, value, color=BLUE if key.startswith("1node") else ORANGE, width=0.62,
                   label=("one node" if key.startswith("1node") else "two nodes") if index in (0, 2) else None)
            ax.text(index, value + 150, f"{value:,.0f}", ha="center", va="bottom", fontsize=VALUE_FONT - 0.5)
            share_inside(ax, index, value, 100 * v["exchange_ms"] / v["step_ms"])
        ax.set_xticks(range(len(layouts)))
        ax.set_xticklabels([name for _, name in layouts])
        ax.set_ylabel("Mcell-updates/s")
        ax.set_ylim(0, 13500)
        ax.set_title(grid, loc="left")
        ax.grid(axis="x", visible=False)
    axes[0].legend(loc="upper left", title="white figure, exchange share", title_fontsize=8.5)
    figstyle.finish_axes(*axes)
    save(fig, out_path, caption(provenance, paths))


def figure_correctness(test_log, out_path):
    """F9, the worst relative difference from the baseline per rung, from a rangpur test log."""
    text = pathlib.Path(test_log).read_text()
    worst = collections.defaultdict(float)
    for name, diff in re.findall(r"^(\w+)\s+(?:call|put)\s*: worst relative diff ([0-9.e+-]+)", text, re.M):
        worst[name] = max(worst[name], float(diff))
    names = [n for n in ("serial", "simd", "omp", "mpi", "mpi_overlap", "omp_float") if n in worst]
    fig, ax = plt.subplots(figsize=HALF)
    floor = 1e-17
    for index, name in enumerate(names):
        value = max(worst[name], floor)
        # The float experiment is held to its own limit, so it is drawn in the rejected colour
        colour = RED if name == "omp_float" else BLUE
        ax.bar(index, value, color=colour, width=0.62)
        ax.text(index, value * 2.5, f"{worst[name]:.1e}", ha="center", fontsize=VALUE_FONT - 0.5)
    ax.axhline(1e-11, color=RED, lw=1, label="tolerance, 1e-11")
    if "omp_float" in worst:
        ax.axhline(1e-4, color=RED, lw=1, ls="--", label="float tolerance, 1e-4")
    ax.axhline(2.2e-16, color=GREY, lw=1, ls=":", label="one rounding, 2.2e-16")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_yscale("log")
    ax.set_ylim(floor, 1e0 if "omp_float" in worst else 1e-9)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=9)
    ax.set_ylabel("worst relative difference")
    ax.grid(axis="x", visible=False)
    figstyle.finish_axes(ax)
    save(fig, out_path, f"price, delta, gamma and vega, call and put, 128 x 32 grid, 5200 steps, from {pathlib.Path(test_log).name}")


def figure_optlevels(opt, provenance, paths, out_path):
    """F11, the baseline and serial rungs built at -O0 to -O3, one thread."""
    levels = ["O0", "O1", "O2", "O3"]
    fig, ax = plt.subplots(figsize=HALF)
    width = 0.38
    for offset, (solver, colour) in zip((-width / 2, width / 2), (("baseline", GREY), ("serial", BLUE))):
        for index, level in enumerate(levels):
            value = opt[f"{solver}_{level}"]["median"] / 1e6
            ax.bar(index + offset, value, width=width, color=colour, label=solver if index == 0 else None)
            ax.text(index + offset, value + 4, f"{value:.0f}", ha="center", va="bottom", fontsize=VALUE_FONT - 0.5)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"-{level}" for level in levels])
    ax.set_ylabel("Mcell-updates/s, 1 thread")
    ax.set_ylim(0, max(v["median"] for v in opt.values()) / 1e6 * 1.25)
    ax.legend(loc="upper left")
    ax.grid(axis="x", visible=False)
    figstyle.finish_axes(ax)
    save(fig, out_path, caption(provenance, paths) + ", each pair built at its own level, every level gives the -O2 answer")


def read_flat_profile(path):
    """The rows of a gprof flat profile as (percent of time, function name)."""
    rows = []
    for line in pathlib.Path(path).read_text().splitlines():
        match = re.match(r"^\s*([\d.]+)\s+[\d.]+\s+[\d.]+\s+(?:\d+\s+[\d.]+\s+[\d.]+\s+)?(\S.*)$", line)
        if match:
            rows.append((float(match.group(1)), match.group(2).strip()))
    return rows


def short_name(name):
    """A gprof function name without its argument list or template arguments, so it fits on an axis."""
    name = name.split("(")[0]
    # Remove the innermost <...> until none are left, std::vector<double, ...> becomes std::vector
    while re.search(r"<[^<>]*>", name):
        name = re.sub(r"<[^<>]*>", "", name)
    return name


def figure_profile(profile_paths, out_path):
    """F12, the gprof flat profile of the baseline, serial and simd rungs on one thread.

    One bar per rung, split into the hottest function and everything else.
    """
    fig, ax = plt.subplots(figsize=(WIDE[0], 1.9))
    solvers = []
    for index, path in enumerate(profile_paths):
        solvers.append(pathlib.Path(path).name.split("_")[1])
        share, name = read_flat_profile(path)[0]
        ax.barh(index, share, color=BLUE, height=0.58)
        ax.barh(index, 100 - share, left=share, color=AMBER, height=0.58)
        ax.text(2, index, f"{short_name(name)}, {share:.1f}%", va="center", fontsize=VALUE_FONT, color="white", weight=600)
        ax.text(101.5, index, f"the rest {100 - share:.1f}%", va="center", fontsize=VALUE_FONT, color=TEXT)
    ax.set_yticks(range(len(solvers)))
    ax.set_yticklabels(solvers)
    ax.invert_yaxis()
    ax.set_xlim(0, 122)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("share of run time (%)")
    ax.grid(axis="y", visible=False)
    figstyle.finish_axes(ax)
    jobs = sorted({pathlib.Path(path).stem.split("_")[-1] for path in profile_paths})
    save(fig, out_path, f"rangpur r730-2, g++ -O2 -pg, sampled every 10 ms, 2048 x 512, one thread, job {', '.join(jobs)}")


# Flops and bytes per interior cell update, the same count as scripts/roofline.py.
# FMA counts as two flops, so the AVX2 row does the same 31 flops as the scalar row
FLOPS_PER_CELL = 31
BYTES_PER_CELL = 24


def read_peaks(path):
    """The measured ceilings from bench_peak.sh, keyed by kind and thread count."""
    peaks = {}
    for line in pathlib.Path(path).read_text().splitlines():
        match = re.match(r"(scalar_peak_gflops|avx2_fma_peak_gflops|triad_bandwidth_gbs) ([0-9.]+)(?:\s+threads (\d+))?", line)
        if match:
            kind, value, threads = match.groups()
            key = (kind, int(threads) if threads else 1)
            # The script runs some cases twice, keep the best
            peaks[key] = max(peaks.get(key, 0.0), float(value))
    return peaks


def draw_roofline(ax, flat_roofs, bandwidth_roofs, points, y_range, labels_right=False):
    """One roofline. Flat roofs are compute ceilings, sloped ones are memory bandwidth.

    labels_right puts the point labels right of the points and the memory label under its line,
    for the single core plot where the left side is crowded.
    """
    intensity = FLOPS_PER_CELL / BYTES_PER_CELL
    xs = [0.1, 20]
    backing = dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.9)
    for value, name, colour, style, side in flat_roofs:
        ax.plot(xs, [value] * 2, color=colour, lw=1.3, ls=style)
        x_text = xs[0] * 1.12 if side == "left" else xs[1] * 0.9
        ax.text(x_text, value * 1.12, f"{name}, {value:.0f} GFLOP/s", ha=side, va="bottom", fontsize=NOTE_FONT, color=colour)
    for bandwidth, name in bandwidth_roofs:
        ax.plot(xs, [bandwidth * x for x in xs], color=ORANGE, lw=1.3)
        # The label lies along the sloped line, low on the left. On log axes a line through the origin
        # points along (x, y) from any point on or near it, so that direction gives the angle in data units
        # and matplotlib turns it into the angle on screen
        x_label = max(xs[0] * 1.15, y_range[0] * 1.1 / bandwidth)
        # Under the line, the text hangs below its baseline, so it never crosses the line
        y_label = bandwidth * x_label * (0.8 if labels_right else 1.25)
        ax.text(x_label, y_label, f"{name}, {bandwidth:.0f} GB/s", fontsize=NOTE_FONT, color=ORANGE,
                va="top" if labels_right else "baseline",
                rotation=math.degrees(math.atan2(y_label, x_label)), rotation_mode="anchor", transform_rotates_text=True)
    ax.axvline(intensity, color=LIGHT_GREY, lw=1, zorder=0)
    ax.text(intensity * 1.08, y_range[0] * 1.15, f"this kernel, {intensity:.2f} flop/byte", va="bottom", fontsize=NOTE_FONT,
            color=DARK_GREY)
    for gflops, name, colour in points:
        ax.plot(intensity, gflops, "o", color=colour, ms=6.5, zorder=3)
        ax.text(intensity * (1.18 if labels_right else 0.84), gflops, f"{name}  {gflops:.1f}", ha="left" if labels_right else "right",
                va="center", fontsize=VALUE_FONT, color=TEXT, bbox=backing, zorder=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xs)
    ax.set_ylim(*y_range)
    ax.set_xlabel("flop per byte of memory traffic")
    ax.set_ylabel("GFLOP/s")


def figure_roofline_core(ladder, peaks, provenance, paths, out_path):
    """F5a, one core. The serial rungs against the scalar and AVX2 ceilings of one core."""
    fig, ax = plt.subplots(figsize=HALF)
    flat = []
    if ("scalar_peak_gflops", 1) in peaks:
        flat.append((peaks[("scalar_peak_gflops", 1)], "scalar peak", GREY, "--", "left"))
    if ("avx2_fma_peak_gflops", 1) in peaks:
        flat.append((peaks[("avx2_fma_peak_gflops", 1)], "AVX2 FMA peak", BLUE, "-", "left"))
    sloped = []
    if ("triad_bandwidth_gbs", 1) in peaks:
        sloped.append((peaks[("triad_bandwidth_gbs", 1)], "memory"))
    points = [(ladder[key]["median"] * FLOPS_PER_CELL / 1e9, key, GREY if key == "baseline" else BLUE)
              for key in ("baseline", "serial", "simd") if key in ladder]
    draw_roofline(ax, flat, sloped, points, (0.5, 120), labels_right=True)
    figstyle.finish_axes(ax)
    save(fig, out_path, caption(provenance, paths))


def figure_roofline_node(ladder, size, peaks, provenance, paths, out_path):
    """F5b, the whole node. The 24 thread rung on two grids against the node's ceilings."""
    fig, ax = plt.subplots(figsize=HALF)
    flat = []
    if ("avx2_fma_peak_gflops", 24) in peaks:
        flat.append((peaks[("avx2_fma_peak_gflops", 24)], "AVX2 FMA peak", BLUE, "-", "left"))
    sloped = []
    if ("triad_bandwidth_gbs", 24) in peaks:
        sloped.append((peaks[("triad_bandwidth_gbs", 24)], "memory"))
    points = []
    for label, name in (("omp24_2048x512", "omp, 2048"), ("omp24_4096x1024", "omp, 4096")):
        if label in size:
            points.append((size[label]["median"] * FLOPS_PER_CELL / 1e9, name, BLUE if "2048" in label else ORANGE))
    draw_roofline(ax, flat, sloped, points, (10, 3000))
    figstyle.finish_axes(ax)
    save(fig, out_path, caption(provenance, paths))


def figure_compiler(e1, provenance, paths, out_path):
    """F15, experiment E1. What g++ makes of each rung under stronger flags, one thread."""
    # The four builds as the flags are typed, which are the rows RESULTS M2.6 quotes. The builds
    # that also keep the compiler's own loads whole are in the full E1 table in the appendix
    builds = [("O2", "-O2\nshipped"), ("O3", "-O3"), ("O3avx2", "-O3\n-mavx2 -mfma"),
              ("O3avx2fast", "-O3 -mavx2 -mfma\n-ffast-math")]
    builds = [(key, name) for key, name in builds if f"baseline_{key}" in e1]
    rungs = [("baseline", GREY), ("serial", "#93B4F5"), ("simd", BLUE)]
    fig, ax = plt.subplots(figsize=WIDE)
    width = 0.26
    for offset, (rung, colour) in zip((-width, 0, width), rungs):
        for index, (key, _) in enumerate(builds):
            v = e1.get(f"{rung}_{key}")
            if not v:
                continue
            value = v["median"] / 1e6
            ax.bar(index + offset, value, width=width * 0.95, color=colour, label=rung if index == 0 else None)
            ax.text(index + offset, value + 8, f"{value:.0f}", ha="center", va="bottom", fontsize=VALUE_FONT - 0.5)
    ax.set_xticks(range(len(builds)))
    ax.set_xticklabels([name for _, name in builds])
    ax.set_ylabel("Mcell-updates/s, 1 thread")
    ax.set_ylim(0, max(v["median"] for v in e1.values()) / 1e6 * 1.2)
    ax.legend(loc="upper left", ncol=3)
    ax.grid(axis="x", visible=False)
    figstyle.finish_axes(ax)
    save(fig, out_path, caption(provenance, paths) + ", whole program rebuilt with each set of flags")


def figure_split(e2, e2b, provenance, paths, out_path):
    """F16, experiments E2 and E2b. The same source with g++ splitting each unaligned load and without."""
    fig, axes = plt.subplots(1, 2, figsize=HALF)
    for ax, key, title in ((axes[0], "simd_1_2048", "simd, 1 thread"), (axes[1], "omp_24_2048", "omp, 24 threads")):
        before = e2.get(key)
        after = e2b.get(key)
        if not before or not after:
            continue
        for index, (v, colour) in enumerate(((before, GREY), (after, BLUE))):
            ax.bar(index, v["median"] / 1e6, color=colour, width=0.62)
            ax.text(index, v["median"] / 1e6 * 1.02, f"{v['median'] / 1e6:,.0f}", ha="center", va="bottom", fontsize=VALUE_FONT)
        ax.text(0.5, after["median"] / 1e6 * 1.2, f"{after['median'] / before['median']:.2f}x", ha="center",
                fontsize=VALUE_FONT + 0.5, weight=600, color=BLUE)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["split\nloads", "unsplit\nloads"])
        ax.set_ylim(0, after["median"] / 1e6 * 1.4)
        ax.set_title(title, loc="left")
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("Mcell-updates/s")
    figstyle.finish_axes(*axes)
    save(fig, out_path, caption(provenance, paths) + ", 2048 x 512")


def figure_overlap(e9, provenance, paths, out_path):
    """F17, experiment E9. The mpi rung against mpi_overlap in every hybrid layout."""
    layouts = [("1node_1x24", "1 node\n1 x 24"), ("1node_2x12", "1 node\n2 x 12"),
               ("2node_1x24", "2 nodes\n1 x 24"), ("2node_2x12", "2 nodes\n2 x 12")]
    fig, axes = plt.subplots(1, 2, figsize=WIDE)
    for ax, (size, grid) in zip(axes, (("2048x512", "2048 x 512"), ("4096x1024", "4096 x 1024"))):
        for index, (key, name) in enumerate(layouts):
            for offset, solver, colour in ((-0.22, "mpi", GREY), (0.22, "mpi_overlap", BLUE)):
                v = e9.get(f"{key}_{size}_{solver}")
                if not v:
                    continue
                value = v["median"] / 1e6
                ax.bar(index + offset, value, width=0.38, color=colour, label=solver if index == 0 else None)
                # The bars are too narrow to hold the share, so it sits just above each one
                share = 100 * v["exchange_ms"] / v["step_ms"]
                if share >= 0.5:
                    ax.text(index + offset, value + 120, f"{share:.0f}%", ha="center", va="bottom", fontsize=7.5, color=DARK_GREY)
        ax.set_xticks(range(len(layouts)))
        ax.set_xticklabels([name for _, name in layouts])
        ax.set_ylabel("Mcell-updates/s")
        ax.set_ylim(0, 12500)
        ax.set_title(grid, loc="left")
        ax.grid(axis="x", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (GREY, BLUE)]
    fig.legend(handles, ["mpi", "mpi_overlap"], loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0), fontsize=9)
    figstyle.finish_axes(*axes)
    save(fig, out_path, caption(provenance, paths) + ", the number above each bar is the exchange share", top=0.9)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("--test-log", help="slurm output holding a test_ladder_matches run")
    parser.add_argument("--peak", help="peak_<jobid>.txt from slurm/bench_peak.sh")
    parser.add_argument("--profile", nargs="*", default=[], help="profile_<solver>_<jobid>.txt from slurm/profile_m2.sh")
    parser.add_argument("--out", default="results/m2")
    args = parser.parse_args()

    rows, provenance = read_csvs(args.csv)
    bad = check_provenance(provenance)
    if bad:
        print("FAIL: these CSVs did not come from a rangpur compute node: " + ", ".join(bad), file=sys.stderr)
        return 1

    figstyle.apply_slide()
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    files = collections.defaultdict(list)
    for row in rows:
        if row["file"] not in files[row["sweep"]]:
            files[row["sweep"]].append(row["file"])

    sweeps = {name: medians(rows, name) for name in
              ("L", "S", "W", "Z", "M", "OPT", "E1", "E2", "E2b", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E11")}
    ladder, strong, weak, size, mpi = sweeps["L"], sweeps["S"], sweeps["W"], sweeps["Z"], sweeps["M"]

    summary = []
    if ladder:
        figure_ladder(ladder, provenance, files["L"], out / "f1_ladder.png")
        for key, v in ladder.items():
            summary.append(f"ladder {key}: median {v['median'] / 1e6:.1f} Mcell-updates/s (min {v['min'] / 1e6:.1f}, max {v['max'] / 1e6:.1f})")
    tried_files = files["E3"] + files["E4"] + files["E11"] + files["E6"]
    if tried_files:
        figure_tried(sweeps["E3"], sweeps["E4"], sweeps["E5"], sweeps["E7"], sweeps["E8"], sweeps["E11"],
                     provenance, tried_files, out / "f1b_tried.png")
    if strong:
        fits = figure_strong(strong, provenance, files["S"], out / "f2_strong.png")
        for solver, f in fits.items():
            summary.append(f"strong {solver}: Amdahl serial fraction {f:.4f}")
        for key, v in sorted(strong.items(), key=lambda item: (item[1]["row"]["solver"], item[1]["row"]["threads"])):
            summary.append(f"strong {key}: {v['median'] / 1e6:.1f} Mcell-updates/s, step {v['step_ms']:.4f} ms, boundary {v['boundary_ms']:.4f} ms")
    if weak:
        serial_fraction, scaled = figure_weak(weak, provenance, files["W"], out / "f3_weak.png")
        summary.append(f"weak: Gustafson serial share {serial_fraction:.4f}, scaled speedup on the most threads {scaled:.2f}")
        for key, v in sorted(weak.items(), key=lambda item: item[1]["row"]["threads"]):
            summary.append(f"weak {key}: step {v['step_ms']:.4f} ms, {v['median'] / 1e6:.1f} Mcell-updates/s")
    if size:
        figure_size(size, provenance, files["Z"], out / "f4_size.png")
        for key, v in sorted(size.items()):
            summary.append(f"size {key}: {v['median'] / 1e6:.1f} Mcell-updates/s")
    if strong or mpi:
        figure_breakdown(strong, mpi, provenance, files["S"] + files["M"], out / "f6_breakdown.png")
    if sweeps["E6"]:
        figure_placement(sweeps["E6"], sweeps["E7"], mpi, provenance, files["E6"] + files["M"], out / "f7_placement.png")
    for sweep in ("E1", "E2", "E2b", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E11"):
        for key, v in sorted(sweeps[sweep].items()):
            line = f"experiment {sweep} {key}: {v['median'] / 1e6:.1f} Mcell-updates/s (min {v['min'] / 1e6:.1f}, max {v['max'] / 1e6:.1f})"
            if v["exchange_ms"] > 0:
                line += f", exchange share {100 * v['exchange_ms'] / v['step_ms']:.1f}%"
            summary.append(line)
    if mpi:
        figure_layouts(mpi, provenance, files["M"], out / "f8_layouts.png")
        for key, v in sorted(mpi.items()):
            summary.append(f"mpi {key}: {v['median'] / 1e6:.1f} Mcell-updates/s, step {v['step_ms']:.4f} ms, exchange {v['exchange_ms']:.4f} ms, boundary {v['boundary_ms']:.4f} ms")
    if args.peak and ladder:
        peaks = read_peaks(args.peak)
        figure_roofline_core(ladder, peaks, provenance, files["L"], out / "f5a_roofline_core.png")
        if size:
            figure_roofline_node(ladder, size, peaks, provenance, files["Z"], out / "f5b_roofline_node.png")
        for (kind, threads), value in sorted(peaks.items()):
            summary.append(f"peak {kind} threads {threads}: {value:.2f}")
        for key, v in ladder.items():
            summary.append(f"roofline ladder {key}: {v['median'] * FLOPS_PER_CELL / 1e9:.1f} GFLOP/s")
        for key, v in size.items():
            if key.startswith("omp24_"):
                summary.append(f"roofline size {key}: {v['median'] * FLOPS_PER_CELL / 1e9:.1f} GFLOP/s")
    if args.test_log:
        figure_correctness(args.test_log, out / "f9_correctness.png")
    if sweeps["OPT"]:
        figure_optlevels(sweeps["OPT"], provenance, files["OPT"], out / "f11_optlevels.png")
        for key, v in sorted(sweeps["OPT"].items()):
            summary.append(f"optlevel {key}: {v['median'] / 1e6:.1f} Mcell-updates/s")
    if args.profile:
        figure_profile(args.profile, out / "f12_profile.png")
        for path in args.profile:
            for share, name in read_flat_profile(path)[:5]:
                summary.append(f"profile {pathlib.Path(path).stem}: {share:.2f}% {name}")
    if sweeps["E1"]:
        figure_compiler(sweeps["E1"], provenance, files["E1"], out / "f15_compiler.png")
    if sweeps["E2"] and sweeps["E2b"]:
        figure_split(sweeps["E2"], sweeps["E2b"], provenance, files["E2"] + files["E2b"], out / "f16_split.png")
    if sweeps["E9"]:
        figure_overlap(sweeps["E9"], provenance, files["E9"], out / "f17_overlap.png")

    (out / "summary.txt").write_text("\n".join(summary) + "\n")
    print(f"wrote {out / 'summary.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
