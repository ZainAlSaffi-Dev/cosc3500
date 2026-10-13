# AI declaration
# I used Claude, an AI assistant made by Anthropic, while working on this file.
# It helped me debug the code, draft pseudocode and implement it.
"""Milestone 2 benchmark figures, F1 to F9, F11 and F12.

Reads the bench_m2_*.csv files from rangpur and draws one PNG per figure.
Refuses any CSV whose host is not a rangpur compute node, as in Milestone 1.

    arch -arm64 .venv/bin/python scripts/m2_figures.py results/bench_m2_*.csv \
        --test-log results/logs/slurm-651147.out --out results/m2
"""

import argparse
import collections
import pathlib
import re
import statistics
import sys

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
    """One footer line naming the node, CPU, compiler, flags, page regime and jobs."""
    first = paths[0]
    jobs = ", ".join(sorted({fact(provenance, p, "jobid") for p in paths}))
    nodes = ", ".join(sorted({fact(provenance, p, "nodes") for p in paths}))
    return (f"rangpur {nodes}, {fact(provenance, first, 'cpu')}, "
            f"{fact(provenance, first, 'compiler').split(' 2021')[0]}, "
            f"{fact(provenance, first, 'flags')} plus {fact(provenance, first, 'avx_flags')} on the AVX2 row, "
            f"huge pages pinned, median of {fact(provenance, first, 'reps')}, jobs {jobs}")


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
            "halo_ms": statistics.median([1000 * r["exchange_seconds"] / r["nt"] for r in group]),
        }
    return out


def save(fig, provenance, paths, out_path):
    fig.text(0.01, 0.005, caption(provenance, paths), fontsize=6.5, color=DARK_GREY)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"wrote {out_path}")


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


def figure_ladder(ladder, exp_a, exp_b, provenance, paths, out_path):
    """F1, the ladder bars with cumulative and per rung speedup, and the rejected variants."""
    order = [("baseline", "baseline"), ("serial", "serial"), ("simd", "simd\n(AVX2)"),
             ("omp24", "omp\n(24 threads)"), ("mpi2x12", "mpi\n(2 ranks x 12)")]
    fig, (ax, side) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [3, 2]})
    base = ladder["baseline"]["median"]
    previous = None
    for index, (key, name) in enumerate(order):
        value = ladder[key]["median"] / 1e6
        colour = GREY if key == "baseline" else BLUE
        # A rung within 3 percent of the one below is a null, drawn quiet as in Milestone 1
        if previous is not None and abs(ladder[key]["median"] / previous - 1.0) < 0.03:
            colour = LIGHT_GREY
        ax.bar(index, value, color=colour, width=0.65)
        ax.errorbar(index, value, yerr=[[value - ladder[key]["min"] / 1e6], [ladder[key]["max"] / 1e6 - value]],
                    color=DARK_GREY, capsize=3, lw=0.8)
        label = f"{value:,.0f}\n{ladder[key]['median'] / base:.1f}x total"
        ax.text(index, value * 1.08, label, ha="center", va="bottom", fontsize=9, color=TEXT)
        if previous is not None:
            ax.text(index, value * 0.5, f"x{ladder[key]['median'] / previous:.2f}", ha="center",
                    va="center", fontsize=9, color="white", weight=600)
        previous = ladder[key]["median"]
    ax.set_yscale("log")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([name for _, name in order], fontsize=9.5)
    ax.set_ylabel("Mcell-updates/s (log scale)")
    ax.set_ylim(top=ax.get_ylim()[1] * 3)
    ax.set_title("Each rung adds one technique, 2048 x 512 grid")

    # The rejected variants against the omp rung they would have replaced
    variants = []
    if exp_a:
        reference = exp_a.get("omp_2048")
        if reference:
            for key, name in [("persistent_2048", "one region per solve (E4)"),
                              ("edges_2048", "parallel boundaries (E8)"),
                              ("dynamic_1", "schedule dynamic,1 (E5)"),
                              ("dynamic_8", "schedule dynamic,8 (E5)"),
                              ("guided", "schedule guided (E5)")]:
                if key in exp_a:
                    variants.append((name, exp_a[key]["median"] / reference["median"]))
            if "default_2048_ratio" in exp_a:
                variants.append(("default page regime (E11)", exp_a["default_2048_ratio"]))
            for size in ("2048", "4096"):
                if f"float_{size}_ratio" in exp_a:
                    variants.append((f"float sheets, {size}, less accurate (E3)", exp_a[f"float_{size}_ratio"]))
    if exp_b:
        for size in ("2048x512", "4096x1024"):
            reference = exp_b.get(f"close_24_{size}")
            if reference:
                for bind in ("false", "spread"):
                    key = f"{bind}_24_{size}"
                    if key in exp_b:
                        variants.append((f"bind {bind}, {size} (E6)", exp_b[key]["median"] / reference["median"]))
            reference = exp_b.get(f"default_24_{size}")
            key = f"interleave_24_{size}"
            if reference and key in exp_b:
                variants.append((f"numactl interleave, {size} (E7)", exp_b[key]["median"] / reference["median"]))
    for index, (name, ratio) in enumerate(variants):
        colour = RED if ratio < 0.98 else (BLUE if ratio > 1.02 else LIGHT_GREY)
        side.barh(index, ratio - 1.0, left=1.0, color=colour, height=0.6)
        side.text(max(ratio, 1.0) + 0.005, index, f"{ratio:.3f}x", va="center", fontsize=8.5)
    side.axvline(1.0, color=GREY, lw=1)
    side.set_yticks(range(len(variants)))
    side.set_yticklabels([name for name, _ in variants], fontsize=8.5)
    side.invert_yaxis()
    side.set_xlabel("throughput relative to the omp rung, 24 threads")
    side.set_title("Tried, measured, compared")
    figstyle.finish_axes(ax, side)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


def figure_strong(strong, provenance, paths, out_path):
    """F2, speedup and efficiency against threads with the Amdahl fit."""
    fig, (ax, eff) = plt.subplots(1, 2, figsize=(12, 5))
    fits = {}
    for solver, colour, name in [("omp", BLUE, "AVX2 row"), ("omp_scalar", GREY, "scalar row")]:
        points = sorted((v["row"]["threads"], v["median"]) for k, v in strong.items() if v["row"]["solver"] == solver)
        if not points:
            continue
        one = points[0][1]
        threads = [p for p, _ in points]
        speedups = [rate / one for _, rate in points]
        f = amdahl_fit(threads[1:], speedups[1:])
        fits[solver] = f
        ax.plot(threads, speedups, "o-", color=colour, label=f"{name}, Amdahl fit f = {f:.3f}")
        model = [1.0 / (f + (1.0 - f) / p) for p in threads]
        ax.plot(threads, model, "--", color=colour, lw=1, alpha=0.7)
        eff.plot(threads, [s / p for s, p in zip(speedups, threads)], "o-", color=colour, label=name)
    top = max(max(v["row"]["threads"] for v in strong.values()), 1)
    ax.plot([1, top], [1, top], color=LIGHT_GREY, lw=1, label="ideal")
    for axis in (ax, eff):
        axis.set_xlabel("OpenMP threads, OMP_PROC_BIND=close")
    # Measured by slurm/diag_placement.sh, close fills CPUs 0, 1, 2 and the node numbers sockets even and odd
    ax.text(24, 2.0, "on this node close alternates sockets,\nso every even count is split evenly", ha="right", color=DARK_GREY, fontsize=8.5)
    ax.set_ylabel("speedup over the same rung on 1 thread")
    ax.set_title("Strong scaling, 2048 x 512")
    ax.legend(loc="upper left")
    eff.set_ylabel("parallel efficiency (speedup / threads)")
    eff.set_ylim(0, 1.1)
    eff.set_title("Efficiency")
    eff.legend()
    figstyle.finish_axes(ax, eff)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)
    return fits


def figure_weak(weak, provenance, paths, out_path):
    """F3, time per step with fixed work per thread, and scaled speedup against Gustafson."""
    points = sorted((v["row"]["threads"], v) for v in weak.values())
    threads = [p for p, _ in points]
    step_ms = [v["step_ms"] for _, v in points]
    fig, (ax, sp) = plt.subplots(1, 2, figsize=(12, 5))
    ax.plot(threads, step_ms, "o-", color=BLUE)
    ax.set_ylim(0, max(step_ms) * 1.3)
    ax.set_xlabel("OpenMP threads (21 interior rows each)")
    ax.set_ylabel("time per step (ms)")
    ax.set_title("Weak scaling, the same work per thread")
    # Scaled speedup is p times the one thread time over the p thread time
    scaled = [p * step_ms[0] / t for p, t in zip(threads, step_ms)]
    sp.plot(threads, scaled, "o-", color=BLUE, label="measured scaled speedup")
    sp.plot(threads, threads, color=LIGHT_GREY, lw=1, label="ideal")
    # Gustafson's f is the serial share of the parallel run, measured on the most threads
    widest = points[-1][1]
    serial_fraction = widest["boundary_ms"] / widest["step_ms"]
    gustafson = [p - serial_fraction * (p - 1) for p in threads]
    sp.plot(threads, gustafson, "--", color=GREY, label=f"Gustafson, measured serial share f = {serial_fraction:.4f}")
    sp.set_xlabel("OpenMP threads")
    sp.set_ylabel("scaled speedup")
    sp.set_title("Gustafson")
    sp.legend(loc="upper left")
    figstyle.finish_axes(ax, sp)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


def figure_size(size, provenance, paths, out_path):
    """F4, throughput against working set with the cache bands, and omp over serial per size."""
    fig, (ax, sp) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [3, 2]})
    lines = {}
    for prefix, colour, name in [("serial", GREY, "serial"), ("simd", BLUE, "simd"), ("omp24", ORANGE, "omp, 24 threads")]:
        points = []
        for label, v in size.items():
            if label.startswith(prefix + "_"):
                row = v["row"]
                mib = 2 * row["ns"] * row["nv"] * 8 / 2**20
                points.append((mib, v["median"] / 1e6, f"{row['ns']}x{row['nv']}"))
        points.sort()
        lines[prefix] = points
        ax.plot([p[0] for p in points], [p[1] for p in points], "o-", color=colour, label=name)
    # The one grid whose two buffers sit exactly 2^24 bytes apart, measured slow by slurm/diag_size_m2.sh.
    # Serial and simd both run slow there, and the note points at the lower of the two.
    slow = [(rate, mib) for prefix in ("serial", "simd") for mib, rate, grid in lines[prefix] if grid == "2896x724"]
    if slow:
        rate, mib = min(slow)
        ax.annotate("serial and simd slow here,\nbuffers exactly 16 MiB apart,\nsee RESULTS.md", (mib, rate), xytext=(mib * 0.08, rate * 1.25),
                    fontsize=8, color=DARK_GREY, arrowprops=dict(arrowstyle="->", color=DARK_GREY, lw=0.8))
    bands = [(0, L2_MIB, "fits in one L2"), (L2_MIB, L3_MIB, "fits in one L3"),
             (L3_MIB, 2 * L3_MIB, "fits in both L3s"), (2 * L3_MIB, 1e4, "main memory")]
    shades = ["#F3F4F6", "#E5E7EB", "#F3F4F6", "#E5E7EB"]
    for (low, high, name), shade in zip(bands, shades):
        ax.axvspan(max(low, 0.1), high, color=shade, zorder=0)
        ax.text((max(low, 0.4) * min(high, 400)) ** 0.5, 1.02, name, transform=ax.get_xaxis_transform(),
                ha="center", fontsize=8, color=DARK_GREY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.3, 400)
    ax.set_xlabel("working set, both sheets (MiB)")
    ax.set_ylabel("Mcell-updates/s")
    ax.set_title("Throughput against grid size", pad=18)
    ax.legend(loc="upper left")
    if lines.get("serial") and lines.get("omp24"):
        names = [p[2] for p in lines["serial"]]
        ratio = [o[1] / s[1] for s, o in zip(lines["serial"], lines["omp24"])]
        # At 2896 by 724 the serial rate is the slow point, so its ratio is not a real speedup.
        colours = [GREY if name == "2896x724" else BLUE for name in names]
        sp.bar(range(len(names)), ratio, color=colours)
        for index, value in enumerate(ratio):
            if names[index] == "2896x724":
                sp.text(index, value * 0.5, f"{value:.0f}x\nserial\nslow", ha="center", va="center", fontsize=9, color="white")
            else:
                sp.text(index, value * 1.02, f"{value:.0f}x", ha="center", fontsize=9)
        sp.set_xticks(range(len(names)))
        sp.set_xticklabels(names, rotation=30, fontsize=8.5)
        sp.set_ylabel("omp on 24 threads over serial")
        sp.set_title("Combined speedup per size")
    figstyle.finish_axes(ax, sp)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


def figure_breakdown(strong, mpi, provenance, paths, out_path):
    """F6, time per step split into the parallel update, the serial boundaries and the edge row exchange."""
    bars = []
    for threads in (1, 12, 24):
        key = f"omp_{threads}"
        if key in strong:
            v = strong[key]
            bars.append((f"omp\n{threads} thr", v["step_ms"] - v["boundary_ms"], v["boundary_ms"], 0.0))
    for key, name in [("1node_1x24_2048x512", "mpi 1 node\n1 x 24"), ("1node_2x12_2048x512", "mpi 1 node\n2 x 12"),
                      ("2node_1x24_2048x512", "mpi 2 nodes\n1 x 24"), ("2node_2x12_2048x512", "mpi 2 nodes\n2 x 12")]:
        if key in mpi:
            v = mpi[key]
            bars.append((name, v["step_ms"] - v["boundary_ms"] - v["halo_ms"], v["boundary_ms"], v["halo_ms"]))
    # One panel for the omp thread counts and one for the mpi layouts, each on a linear axis
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [3, 4]})
    panels = [[b for b in bars if b[0].startswith("omp")], [b for b in bars if b[0].startswith("mpi")]]
    for ax, panel, title in zip(axes, panels, ("OpenMP, one process", "MPI layouts, 2 ranks per node or 1")):
        for index, (name, update, boundary, halo) in enumerate(panel):
            ax.bar(index, update, color=BLUE, label="interior rows, parallel" if index == 0 else None)
            ax.bar(index, boundary, bottom=update, color=AMBER, label="boundaries, serial" if index == 0 else None)
            if halo > 0:
                ax.bar(index, halo, bottom=update + boundary, color=RED, label="edge row exchange" if index == 0 else None)
            total = update + boundary + halo
            share = 100 * (boundary + halo) / total
            ax.text(index, total * 1.02, f"{total:.3f} ms\n{share:.1f}% not parallel", ha="center", fontsize=8.5)
        ax.set_xticks(range(len(panel)))
        ax.set_xticklabels([b[0] for b in panel], fontsize=9)
        ax.set_ylim(0, max(b[1] + b[2] + b[3] for b in panel) * 1.25)
        ax.set_ylabel("time per step (ms)")
        ax.set_title(title)
        ax.legend(loc="upper right", fontsize=8.5)
    fig.suptitle("Where the time goes, 2048 x 512", fontsize=14, weight=600, color=TEXT)
    figstyle.finish_axes(*axes)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


def figure_placement(exp_b, provenance, paths, out_path):
    """F7, binding and page placement at 6, 12 and 24 threads on both grids."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    for ax, size in zip(axes, ("2048x512", "4096x1024")):
        groups = [6, 12, 24]
        binds = [("false", LIGHT_GREY), ("close", BLUE), ("spread", ORANGE)]
        # What each binding means on this node, measured by slurm/diag_placement.sh
        meaning = {"false": "not bound", "close": "alternates sockets", "spread": "socket 0 only below 24"}
        for offset, (bind, colour) in enumerate(binds):
            values = []
            for threads in groups:
                v = exp_b.get(f"{bind}_{threads}_{size}")
                values.append(v["median"] / 1e9 if v else 0.0)
            ax.bar([g + (offset - 1) * 1.2 for g in groups], values, width=1.2, color=colour, label=f"OMP_PROC_BIND={bind}, {meaning[bind]}")
        v = exp_b.get(f"interleave_24_{size}")
        if v:
            ax.bar(24 + 2.4, v["median"] / 1e9, width=1.2, color=RED, label="close + numactl --interleave=all")
        ax.set_xticks(groups)
        ax.set_xlabel("OpenMP threads")
        ax.set_ylabel("Gcell-updates/s")
        ax.set_title(f"Thread and page placement, {size}")
        ax.legend(fontsize=8)
    figstyle.finish_axes(*axes)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


def figure_layouts(mpi, provenance, paths, out_path):
    """F8, hybrid layouts on one and two nodes, with the edge row exchange share of each step."""
    layouts = [("1node_1x24", "1 node\n1 rank x 24"), ("1node_2x12", "1 node\n2 ranks x 12"),
               ("2node_1x24", "2 nodes\n1 rank x 24 each"), ("2node_2x12", "2 nodes\n2 ranks x 12 each")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, size in zip(axes, ("2048x512", "4096x1024")):
        for index, (key, name) in enumerate(layouts):
            v = mpi.get(f"{key}_{size}")
            if not v:
                continue
            value = v["median"] / 1e9
            colour = BLUE if key.startswith("1node") else ORANGE
            ax.bar(index, value, color=colour)
            halo_share = 100 * v["halo_ms"] / v["step_ms"]
            ax.text(index, value * 1.02, f"{value:.2f}\nexchange {halo_share:.1f}%", ha="center", fontsize=8.5)
        ax.set_xticks(range(len(layouts)))
        ax.set_xticklabels([name for _, name in layouts], fontsize=8.5)
        ax.set_ylabel("Gcell-updates/s")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
        ax.set_title(f"Hybrid MPI and OpenMP, {size}")
    figstyle.finish_axes(*axes)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


def figure_correctness(test_log, out_path):
    """F9, the worst relative difference from the baseline per rung, from a rangpur test log."""
    text = pathlib.Path(test_log).read_text()
    worst = collections.defaultdict(float)
    for name, diff in re.findall(r"^(\w+)\s+(?:call|put)\s*: worst relative diff ([0-9.e+-]+)", text, re.M):
        worst[name] = max(worst[name], float(diff))
    names = [n for n in ("serial", "simd", "omp", "mpi", "omp_float") if n in worst]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    floor = 1e-17
    for index, name in enumerate(names):
        value = max(worst[name], floor)
        # The float experiment is held to its own limit, so it is drawn in the rejected colour
        colour = RED if name == "omp_float" else BLUE
        ax.bar(index, value, color=colour)
        ax.text(index, value * 1.5, f"{worst[name]:.2e}", ha="center", fontsize=9)
    ax.axhline(1e-11, color=RED, lw=1, label="test tolerance, 1e-11")
    if "omp_float" in worst:
        ax.axhline(1e-4, color=RED, lw=1, ls="--", label="float tolerance, 1e-4")
    ax.axhline(2.2e-16, color=GREY, lw=1, ls=":", label="double precision rounding, 2.2e-16")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_yscale("log")
    if "omp_float" in worst:
        ax.set_ylim(floor, 1e-1)
    else:
        ax.set_ylim(floor, 1e-9)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_ylabel("worst relative difference from baseline")
    ax.set_title("Every rung gives the same answer")
    figstyle.finish_axes(ax)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.text(0.01, 0.005, f"price, delta, gamma and vega, call and put, 128 x 32 grid, 5200 steps, from {pathlib.Path(test_log).name}",
             fontsize=6.5, color=DARK_GREY)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"wrote {out_path}")


def figure_optlevels(opt, provenance, paths, out_path):
    """F11, the baseline and serial rungs built at -O0 to -O3, one thread."""
    levels = ["O0", "O1", "O2", "O3"]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    width = 0.38
    for offset, (solver, colour) in zip((-width / 2, width / 2), (("baseline", GREY), ("serial", BLUE))):
        slowest = opt[f"{solver}_O0"]["median"]
        for index, level in enumerate(levels):
            value = opt[f"{solver}_{level}"]["median"] / 1e6
            ax.bar(index + offset, value, width=width, color=colour, label=solver if index == 0 else None)
            ax.text(index + offset, value + 3, f"{value:.0f}\n{opt[f'{solver}_{level}']['median'] / slowest:.1f}x",
                    ha="center", va="bottom", fontsize=8.5, color=TEXT)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"-{level}" for level in levels])
    ax.set_ylabel("Mcell-updates/s, one thread")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.2)
    ax.legend(loc="upper left")
    # The footer shows the default flags, so say that each pair replaced -O2 with its own level
    ax.text(0.99, 0.97, "each pair built with its own -O level in place of -O2,\nevery level gives exactly the -O2 answer",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=DARK_GREY)
    ax.set_title("What the compiler gives before any change of ours, 2048 x 512")
    figstyle.finish_axes(ax)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


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
    """F12, the gprof flat profile of the baseline, serial and simd rungs on one thread."""
    fig, axes = plt.subplots(1, len(profile_paths), figsize=(5 * len(profile_paths), 3.8))
    for ax, path in zip(axes, profile_paths):
        solver = pathlib.Path(path).name.split("_")[1]
        rows = read_flat_profile(path)[:5]
        names = [short_name(name) for _, name in rows]
        shares = [share for share, _ in rows]
        ax.barh(range(len(rows)), shares, color=BLUE, height=0.6)
        for index, share in enumerate(shares):
            ax.text(share + 1, index, f"{share:.1f}%", va="center", fontsize=8.5)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(names, fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlim(0, 115)
        ax.set_xlabel("share of run time (%)")
        ax.set_title(solver)
    figstyle.finish_axes(*axes)
    fig.suptitle("gprof flat profile, one thread, 2048 x 512", fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    jobs = sorted({pathlib.Path(path).stem.split("_")[-1] for path in profile_paths})
    fig.text(0.01, 0.005, f"rangpur r730-2, g++ -O2 -pg, sampled every 10 ms, job {', '.join(jobs)}", fontsize=6.5, color=DARK_GREY)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"wrote {out_path}")


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


def figure_roofline(ladder, size, peaks, provenance, paths, out_path):
    """F5, the measured compute and bandwidth roofs with every rung and the omp size sweep."""
    intensity = FLOPS_PER_CELL / BYTES_PER_CELL
    fig, ax = plt.subplots(figsize=(10, 6))
    xs = [0.05, 50]
    roofs = [(("scalar_peak_gflops", 1), "scalar peak, 1 core", GREY),
             (("avx2_fma_peak_gflops", 1), "AVX2 FMA peak, 1 core", BLUE),
             (("avx2_fma_peak_gflops", 24), "AVX2 FMA peak, 24 cores", BLUE)]
    for key, name, colour in roofs:
        if key in peaks:
            ax.plot(xs, [peaks[key]] * 2, color=colour, lw=1.2, ls="-" if key[1] == 24 else "--")
            ax.text(xs[1], peaks[key] * 1.05, f"{name} {peaks[key]:.0f} GFLOP/s", ha="right", fontsize=8.5, color=colour)
    for threads, style in ((1, "--"), (24, "-")):
        key = ("triad_bandwidth_gbs", threads)
        if key in peaks:
            bandwidth = peaks[key]
            ax.plot(xs, [bandwidth * x for x in xs], color=ORANGE, lw=1.2, ls=style)
            ax.text(0.15, bandwidth * 0.15 * 1.4, f"DRAM triad, {threads} thread{'s' if threads > 1 else ''}, {bandwidth:.0f} GB/s",
                    fontsize=8.5, color=ORANGE, rotation=0)
    ax.axvline(intensity, color=LIGHT_GREY, lw=1)
    ax.text(intensity * 1.05, 0.12, f"this kernel, {intensity:.2f} flop/byte", fontsize=8.5, color=DARK_GREY)
    names = {"baseline": "baseline", "serial": "serial", "simd": "simd", "omp24": "omp 24", "mpi2x12": "mpi 2x12"}
    for key, name in names.items():
        if key in ladder:
            gflops = ladder[key]["median"] * FLOPS_PER_CELL / 1e9
            ax.plot(intensity, gflops, "o", color=GREY if key == "baseline" else BLUE, ms=7)
            # omp and mpi land almost on top of each other, so one label goes above and one below
            nudge = {"omp24": 0.8, "mpi2x12": 1.25}.get(key, 1.0)
            ax.text(intensity * 0.8, gflops * nudge, f"{name} {gflops:.1f}", ha="right", va="center", fontsize=8.5)
    for label, v in size.items():
        if label.startswith("omp24_"):
            gflops = v["median"] * FLOPS_PER_CELL / 1e9
            ax.plot(intensity * 1.25, gflops, "s", color=ORANGE, ms=5)
            ax.text(intensity * 1.4, gflops, label.replace("omp24_", "omp 24, "), va="center", fontsize=7.5, color=DARK_GREY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xs)
    ax.set_ylim(0.1, 3000)
    ax.set_xlabel("arithmetic intensity (flop per byte of DRAM traffic)")
    ax.set_ylabel("GFLOP/s")
    ax.set_title("Roofline, measured ceilings on r730-2")
    figstyle.finish_axes(ax)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, provenance, paths, out_path)


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

    figstyle.apply()
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    files = collections.defaultdict(list)
    for row in rows:
        if row["file"] not in files[row["sweep"]]:
            files[row["sweep"]].append(row["file"])

    ladder = medians(rows, "L")
    exp_a = medians(rows, "E4") | medians(rows, "E5") | medians(rows, "E8") | medians(rows, "E11")
    # E11 ran in its own job, so its omp reference is its own pinned run
    if "pinned_2048" in exp_a and "default_2048" in exp_a:
        exp_a["default_2048_ratio"] = exp_a["default_2048"]["median"] / exp_a["pinned_2048"]["median"]
    # E3 ran in its own job, so each float run is compared with that job's omp run on the same grid
    exp_float = medians(rows, "E3")
    for size in ("2048", "4096"):
        if f"float_24_{size}" in exp_float and f"omp_24_{size}" in exp_float:
            exp_a[f"float_{size}_ratio"] = exp_float[f"float_24_{size}"]["median"] / exp_float[f"omp_24_{size}"]["median"]
    exp_b = medians(rows, "E6") | medians(rows, "E7")
    strong = medians(rows, "S")
    weak = medians(rows, "W")
    size = medians(rows, "Z")
    mpi = medians(rows, "M")

    summary = []
    if ladder:
        paths = files["L"] + files["E4"] + files["E11"] + files["E3"] + files["E6"]
        figure_ladder(ladder, exp_a, exp_b, provenance, [p for p in paths if p], out / "f1_ladder.png")
        for key, v in ladder.items():
            summary.append(f"ladder {key}: median {v['median'] / 1e6:.1f} Mcell-updates/s (min {v['min'] / 1e6:.1f}, max {v['max'] / 1e6:.1f})")
    if strong:
        fits = figure_strong(strong, provenance, files["S"], out / "f2_strong.png")
        for solver, f in fits.items():
            summary.append(f"strong {solver}: Amdahl serial fraction {f:.4f}")
        for key, v in sorted(strong.items(), key=lambda item: (item[1]["row"]["solver"], item[1]["row"]["threads"])):
            summary.append(f"strong {key}: {v['median'] / 1e6:.1f} Mcell-updates/s, step {v['step_ms']:.4f} ms, boundary {v['boundary_ms']:.4f} ms")
    if weak:
        figure_weak(weak, provenance, files["W"], out / "f3_weak.png")
        for key, v in sorted(weak.items(), key=lambda item: item[1]["row"]["threads"]):
            summary.append(f"weak {key}: step {v['step_ms']:.4f} ms, {v['median'] / 1e6:.1f} Mcell-updates/s")
    if size:
        figure_size(size, provenance, files["Z"], out / "f4_size.png")
        for key, v in sorted(size.items()):
            summary.append(f"size {key}: {v['median'] / 1e6:.1f} Mcell-updates/s")
    if strong or mpi:
        figure_breakdown(strong, mpi, provenance, files["S"] + files["M"], out / "f6_breakdown.png")
    if exp_b:
        figure_placement(exp_b, provenance, files["E6"], out / "f7_placement.png")
    for key, v in sorted((exp_a | exp_b | exp_float).items()):
        if not isinstance(v, dict):
            continue
        summary.append(f"experiment {key}: {v['median'] / 1e6:.1f} Mcell-updates/s (min {v['min'] / 1e6:.1f}, max {v['max'] / 1e6:.1f})")
    if mpi:
        figure_layouts(mpi, provenance, files["M"], out / "f8_layouts.png")
        for key, v in sorted(mpi.items()):
            summary.append(f"mpi {key}: {v['median'] / 1e6:.1f} Mcell-updates/s, step {v['step_ms']:.4f} ms, exchange {v['halo_ms']:.4f} ms, boundary {v['boundary_ms']:.4f} ms")
    if args.peak and ladder:
        peaks = read_peaks(args.peak)
        figure_roofline(ladder, size, peaks, provenance, files["L"] + files["Z"], out / "f5_roofline.png")
        for (kind, threads), value in sorted(peaks.items()):
            summary.append(f"peak {kind} threads {threads}: {value:.2f}")
    if args.test_log:
        figure_correctness(args.test_log, out / "f9_correctness.png")
    opt = medians(rows, "OPT")
    if opt:
        figure_optlevels(opt, provenance, files["OPT"], out / "f11_optlevels.png")
        for key, v in sorted(opt.items()):
            summary.append(f"optlevel {key}: {v['median'] / 1e6:.1f} Mcell-updates/s")
    if args.profile:
        figure_profile(args.profile, out / "f12_profile.png")
        for path in args.profile:
            for share, name in read_flat_profile(path)[:5]:
                summary.append(f"profile {pathlib.Path(path).stem}: {share:.2f}% {name}")

    (out / "summary.txt").write_text("\n".join(summary) + "\n")
    print(f"wrote {out / 'summary.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
