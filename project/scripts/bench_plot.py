#!/usr/bin/env python3
"""Benchmark figures (PLAN P7 / §4b): rangpur bench CSVs -> ladder + scaling.

Input: the CSVs written by slurm/bench_{ladder,scaling,controls}.sh. Each has
'#' header comments carrying host/compiler/flags/date/cache sizes, then rows:

    sweep,build,price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec

Sweeps:  A = ablation ladder    B = scaling      C = no-vectorise control
         D = negative controls

Never plot login-node or laptop numbers as results. The provenance header is
reproduced on every figure so a reader can check where the numbers came from.

Outputs (into --out):
  bench_ladder.png    cumulative ladder bars, median with min/max whiskers,
                      the two negative controls plotted separately and
                      labelled as controls, and a footnote saying why the
                      null rungs came out null
  bench_scaling.png   speedup vs grid size, plus achieved throughput
  bench_summary.md    the table the video script reads from

Usage:
    python3 scripts/bench_plot.py results/bench_*.csv --out results/
"""

import argparse
import collections
import pathlib
import statistics
import sys

# The palette and rcParams shared by every figure in the project.
from figstyle import (BLUE, GREY, LIGHT_GREY, DARK_GREY, RED, AMBER,
                      TEXT, GRID)
import figstyle

# Ladder rung -> the technique it adds, in the lecturer's vocabulary (L03).
LADDER_TECHNIQUE = {
    "baseline": "reference solver",
    "opt-L0": "none (anchor)",
    "opt-L1": "hoisting + lookup + CSE",
    "opt-L2": "strength reduction",
    "opt-L3": "traversal order",
    "opt-L4": "loop splitting",
    "opt-L5": "induction variable",
    "opt-L6": "unrolling x4",
    "opt-ctl-order": "CONTROL: loops swapped",
    "opt-ctl-branch": "CONTROL: v=0 row fused",
}

# Shorter names for the ladder's x-axis, so horizontal two-line tick labels
# do not run into their neighbours. The full names stay in the table.
LADDER_TICK = {
    "baseline": "reference",
    "opt-L0": "anchor",
    "opt-L1": "hoisting, CSE",
    "opt-L2": "strength\nreduction",
    "opt-L3": "traversal\norder",
    "opt-L4": "loop\nsplitting",
    "opt-L5": "induction\nvariables",
    "opt-L6": "unroll\nby four",
}
LADDER_ORDER = ["baseline", "opt-L0", "opt-L1", "opt-L2", "opt-L3",
                "opt-L4", "opt-L5", "opt-L6"]
CONTROLS = ["opt-ctl-order", "opt-ctl-branch"]

# A rung within this fraction of the previous one gets reported as a null
# result instead of a gain. 3% is about the run-to-run spread we see across
# five reps on a shared node.
NULL_BAND = 0.03


class Row:
    """One --bench rep. Plain record, one per CSV line."""

    def __init__(self, fields):
        self.sweep = fields[0]
        self.build = fields[1]
        self.price = float(fields[2])
        self.ns = int(fields[6])
        self.nv = int(fields[7])
        self.nt = int(fields[8])
        self.solver = fields[9]
        self.seconds = float(fields[10])
        self.cups = float(fields[11])

    @property
    def grid(self):
        return f"{self.ns}x{self.nv}"


def read_csvs(paths):
    """Returns (rows, provenance-lines). Malformed lines get reported rather
    than quietly dropped, so a truncated CSV from a killed job is visible."""
    rows, provenance, skipped = [], [], 0
    for path in paths:
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    text = line.lstrip("# ")
                    if not text.startswith("sweep,"):  # column header, not info
                        provenance.append(text)
                    continue
                fields = line.split(",")
                if len(fields) != 12:
                    skipped += 1
                    continue
                try:
                    rows.append(Row(fields))
                except ValueError:
                    skipped += 1
    if skipped:
        print(f"WARNING: skipped {skipped} malformed CSV line(s)",
              file=sys.stderr)
    return rows, provenance


def aggregate(rows, key):
    """Group rows by `key(row)` -> dict of median/min/max cell-updates/sec."""
    buckets = collections.defaultdict(list)
    for row in rows:
        buckets[key(row)].append(row.cups)
    return {
        name: {
            "median": statistics.median(values),
            "min": min(values),
            "max": max(values),
            "reps": len(values),
        }
        for name, values in buckets.items()
    }


def provenance_caption(provenance):
    """One short line for the figure corner: host, cpu, compiler, date."""
    facts = {}
    for text in provenance:
        if "=" in text:
            key, _, value = text.partition("=")
            facts.setdefault(key.strip(), value.strip())
    bits = []
    for key in ("host", "cpu", "compiler", "flags_default", "date"):
        if key in facts:
            bits.append(facts[key].strip("'"))
    return "  |  ".join(bits)


# A compute node on rangpur is named <node>.compute.eait.uq.edu.au. The LOGIN
# node is rangpur.eait.uq.edu.au and a laptop is anything else, so requiring
# this substring rejects both. Kept as a constant so the check below and the
# error message can never disagree about what counts as valid.
COMPUTE_NODE_SUFFIX = ".compute.eait.uq.edu.au"


def provenance_facts(provenance):
    """The '# key=value' header lines, collected into a dict of key -> values."""
    facts = collections.defaultdict(list)
    for text in provenance:
        if "=" in text:
            key, _, value = text.partition("=")
            facts[key.strip()].append(value.strip().strip("'"))
    return facts


def check_provenance(provenance):
    """The real gate: every CSV must name a rangpur COMPUTE node.

    The module docstring has always promised not to plot login-node or laptop
    numbers, but the only check used to be "did any rows parse", which a file
    full of laptop timings passes just as happily. Since the figures are
    rendered from a results/bench_*.csv glob and `make fetch` merges rather
    than replaces, a stray local CSV dropped in that directory would have been
    silently folded into the published tables. This makes the promise real.

    Returns a list of complaint strings; empty means the data is clean.
    """
    facts = provenance_facts(provenance)
    hosts = facts.get("host", [])
    complaints = []
    if not hosts:
        complaints.append(
            "no '# host=' header in any CSV — cannot prove these numbers came "
            "from a compute node. Regenerate with slurm/bench_*.sh, which "
            "stamps the header via bench_common.sh:write_header.")
    for host in sorted(set(hosts)):
        if not host.endswith(COMPUTE_NODE_SUFFIX):
            complaints.append(
                f"host={host!r} is not a rangpur compute node (expected a name "
                f"ending {COMPUTE_NODE_SUFFIX}). Login-node and laptop numbers "
                f"are not results.")
    if not facts.get("jobid"):
        complaints.append("no '# jobid=' header — no sbatch job to trace back to.")
    return complaints


def plot_ladder(ladder, controls, control_reference, provenance, out_path,
                regime=""):
    """7-bar cumulative ladder + the two controls, plotted separately."""
    import matplotlib
    matplotlib.use("Agg")  # file output only, no display needed
    figstyle.apply()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    names = [n for n in LADDER_ORDER if n in ladder]
    if not names:
        return None
    values = [ladder[n]["median"] / 1e6 for n in names]
    lows = [values[i] - ladder[n]["min"] / 1e6 for i, n in enumerate(names)]
    highs = [ladder[n]["max"] / 1e6 - values[i] for i, n in enumerate(names)]

    control_names = [n for n in CONTROLS if n in controls]
    control_values = [controls[n]["median"] / 1e6 for n in control_names]

    fig, (ax, ax_ctl) = plt.subplots(
        1, 2, figsize=(12.8, 5.5), dpi=160,
        gridspec_kw={"width_ratios": [len(names), max(2, len(control_names))]})

    # --- left panel, the ladder itself ---
    # Colour follows meaning and nothing else. The rungs that moved the
    # needle are blue, the null rungs are a quiet grey so they read as
    # honestly flat rather than as alarms, and the baseline is a darker grey
    # anchor. Each bar carries its value and its speedup over the baseline;
    # the per-step percentage only appears where a step actually happened,
    # because writing "+0% step" on six bars in a row is noise.
    colours = [GREY]
    steps = [None]
    for i in range(1, len(names)):
        gain = (values[i] - values[i - 1]) / values[i - 1]
        steps.append(gain)
        colours.append(BLUE if abs(gain) > NULL_BAND else LIGHT_GREY)
    bars = ax.bar(range(len(names)), values, color=colours, zorder=3,
                  width=0.66,
                  yerr=[lows, highs], capsize=2,
                  error_kw={"ecolor": DARK_GREY, "linewidth": 0.9})
    for i, bar in enumerate(bars):
        centre = bar.get_x() + bar.get_width() / 2
        speedup = values[i] / values[0]
        ax.annotate(f"{values[i]:.0f}", xy=(centre, values[i]),
                    xytext=(0, 16), textcoords="offset points",
                    ha="center", fontsize=10.5, fontweight=600, color=TEXT)
        ax.annotate(f"{speedup:.2f}x", xy=(centre, values[i]),
                    xytext=(0, 5), textcoords="offset points",
                    ha="center", fontsize=8.5, color=DARK_GREY)
        if steps[i] is not None and abs(steps[i]) > NULL_BAND:
            ax.annotate(f"{100 * steps[i]:+.0f}%", ha="center", va="center",
                        xy=(centre, values[i] / 2),
                        fontsize=11, fontweight=700, color="white")
    legend_handles = [
        Patch(facecolor=BLUE, label="a real gain (more than 3%)"),
        Patch(facecolor=LIGHT_GREY, label="null, within run-to-run spread"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", ncol=2,
              handlelength=1.2, columnspacing=1.2, fontsize=9)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(
        [f"{n.replace('opt-', '')}\n{LADDER_TICK.get(n, '')}"
         for n in names], fontsize=8.5)
    ax.set_ylabel("throughput (million cell-updates / sec)")
    ax.set_title("The optimisation ladder, one technique per rung"
                 + (f"\n{regime}" if regime else ""))
    ax.set_ylim(0, max(values) * 1.24)

    # --- right panel, the negative controls ---
    # Separate axes, red, and labelled CONTROL, because these are deliberately
    # worse code and shouldn't be mistaken for ladder rungs.
    if control_names and control_reference:
        # The reference has to come from the controls' own sweep. The controls
        # job runs a shorter time loop than the ladder job, and on this node a
        # shorter run spends more of itself at the turbo clock, so comparing a
        # sweep-D control against a sweep-A rung would be measuring the run
        # length instead of the technique. That is why the controls job
        # re-runs levels 5 and 6 alongside them.
        reference = control_reference
        ctl_bars = ax_ctl.bar(range(len(control_names)), control_values,
                              color=RED, zorder=3, width=0.6)
        ax_ctl.axhline(reference, color=BLUE, linestyle="--", linewidth=1.3,
                       zorder=2)
        ax_ctl.annotate("level 5, same job",
                        xy=(len(control_names) - 0.55, reference),
                        xytext=(0, 5), textcoords="offset points",
                        ha="right", fontsize=9, color=BLUE)
        for i, bar in enumerate(ctl_bars):
            centre = bar.get_x() + bar.get_width() / 2
            slow = reference / control_values[i] if control_values[i] else 0
            ax_ctl.annotate(f"{control_values[i]:.0f}",
                            xy=(centre, control_values[i]),
                            xytext=(0, 16), textcoords="offset points",
                            ha="center", fontsize=10.5, fontweight=600,
                            color=TEXT)
            ax_ctl.annotate(f"{slow:.1f}x slower",
                            xy=(centre, control_values[i]),
                            xytext=(0, 5), textcoords="offset points",
                            ha="center", fontsize=8.5, color=DARK_GREY)
        ax_ctl.set_xticks(range(len(control_names)))
        ax_ctl.set_xticklabels(
            [n.replace("opt-", "") for n in control_names], fontsize=9)
        ax_ctl.set_ylim(0, max(reference, max(control_values)) * 1.28)
    ax_ctl.set_title("Negative controls\nwhat levels 3 and 4\nare worth")

    figstyle.finish_axes(ax, ax_ctl)
    fig.text(0.01, 0.005, provenance_caption(provenance), fontsize=6.5,
             color=DARK_GREY)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_scaling(scaling, provenance, out_path):
    """Speedup and achieved throughput against grid size (sweep B)."""
    import matplotlib
    matplotlib.use("Agg")
    figstyle.apply()
    import matplotlib.pyplot as plt

    grids = sorted({g for (g, _) in scaling},
                   key=lambda g: int(g.split("x")[0]) * int(g.split("x")[1]))
    pairs = [(g, scaling.get((g, "baseline")), scaling.get((g, "opt-L6")))
             for g in grids]
    pairs = [p for p in pairs if p[1] and p[2]]
    if not pairs:
        return None

    labels = [p[0] for p in pairs]
    base = [p[1]["median"] / 1e6 for p in pairs]
    opt = [p[2]["median"] / 1e6 for p in pairs]
    speedup = [o / b for o, b in zip(opt, base)]
    # The working set is the two buffers together, measured in MiB. The point
    # where it passes the last cache level is where being memory-bound starts
    # to show up in the numbers.
    footprint = [2 * 8 * int(g.split("x")[0]) * int(g.split("x")[1]) / 2**20
                 for g in labels]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.5), dpi=160)

    ax.plot(range(len(labels)), base, "o-", color=GREY, linewidth=2.2,
            markersize=8, zorder=3)
    ax.plot(range(len(labels)), opt, "o-", color=BLUE, linewidth=2.2,
            markersize=8, zorder=3)
    # Direct labels at the line ends instead of a legend box. Two series is
    # few enough that naming them where they finish reads faster.
    ax.annotate("level 6", xy=(len(labels) - 1, opt[-1]),
                xytext=(8, 0), textcoords="offset points", va="center",
                fontsize=10, fontweight=600, color=BLUE)
    ax.annotate("baseline", xy=(len(labels) - 1, base[-1]),
                xytext=(8, 0), textcoords="offset points", va="center",
                fontsize=10, fontweight=600, color=GREY)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([f"{g}\n{f:.0f} MiB" for g, f in zip(labels, footprint)],
                       fontsize=9)
    ax.set_xlim(-0.4, len(labels) - 0.25)
    ax.set_ylabel("throughput (million cell-updates / sec)")
    ax.set_title("Throughput against grid size")

    ax2.plot(range(len(labels)), speedup, "o-", color=BLUE, linewidth=2.2,
             markersize=8, zorder=3)
    ax2.axhline(1.0, color=GREY, linestyle="--", linewidth=1.2, zorder=2)
    for i, value in enumerate(speedup):
        ax2.annotate(f"{value:.2f}x", xy=(i, value), xytext=(0, 9),
                     textcoords="offset points", ha="center", fontsize=10,
                     fontweight=600, color=TEXT)
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("level 6 over baseline")
    ax2.set_title("Speedup against grid size")
    ax2.set_ylim(0, max(speedup) * 1.35)

    # Show the page-regime split as two quiet background bands rather than
    # as labels on the points, which collided with the tick labels. The
    # re-run (job 556834) pins MALLOC_MMAP_THRESHOLD_ so every rep of a row
    # is served the same way, but a buffer needs a 2 MiB-aligned 2 MiB
    # stretch before transparent huge pages can back it, and in practice
    # only the 8 MiB and 32 MiB buffers got them. The two small grids
    # measure at sweep A's 4 KiB throughput and the two big ones at the
    # pinned huge-page throughput of sweeps E and F, so the split drawn
    # here records the measurement rather than the mmap arithmetic. Without
    # it the shape of the curve reads as a cache story, and part of it is
    # a page story.
    split = next((i for i, g in enumerate(labels)
                  if 8 * int(g.split("x")[0]) * int(g.split("x")[1])
                  >= 8 * 2**20), len(labels)) - 0.5
    left, right = ax.get_xlim()
    ax.axvspan(split, right, color="#EFF3FD", zorder=0)
    ax.text(0.5 * (left + split), 0.03, "4 KiB pages",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=9, fontstyle="italic", color=DARK_GREY)
    ax.text(0.5 * (split + right), 0.03, "2 MiB huge pages",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=9, fontstyle="italic", color=DARK_GREY)
    ax.set_xlim(left, right)

    figstyle.finish_axes(ax, ax2)
    fig.text(0.01, 0.005, provenance_caption(provenance), fontsize=6.5,
             color=DARK_GREY)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def write_summary(ladders, controls, control_reference, control_sweep,
                  scaling, novec, provenance, out_path):
    """The markdown table the video script reads from."""
    lines = ["# Serial optimisation results (PLAN P7)", ""]
    lines.append("Provenance (from the CSV headers, deduplicated):")
    lines.append("")
    seen = set()
    for text in provenance:
        if text in seen:
            continue
        seen.add(text)
        lines.append(f"- `{text}`")
    lines.append("")

    REGIME = {
        "A": ("Sweep A: ablation ladder, default allocator (4 KiB pages)",
              "This is what median-of-5 measures with the allocator left "
              "alone: rep 1 gets freshly mmap'd, huge-page-backed memory and "
              "reps 2-5 are served from the reused glibc heap on 4 KiB pages."),
        "E": ("Sweep E: ablation ladder, huge pages",
              "MALLOC_MMAP_THRESHOLD_ pinned below the buffer size so every "
              "rep gets mmap'd, transparent-huge-page-backed memory. **A "
              "single production solve is a fresh process, so this is the "
              "case a user actually gets.**"),
    }

    ladder = ladders.get("E") or ladders.get("A") or {}
    for sweep in ("E", "A"):
        entries = ladders.get(sweep)
        if not entries:
            continue
        names = [n for n in LADDER_ORDER if n in entries]
        if not names:
            continue
        title, blurb = REGIME[sweep]
        anchor = entries[names[0]]["median"]
        lines += [f"## {title}", "", blurb, "",
                  "| level | technique added | Mcell-updates/s "
                  "(median) | min-max | vs baseline | step gain |",
                  "|---|---|---|---|---|---|"]
        previous = None
        for name in names:
            entry = entries[name]
            step = ("-" if previous is None
                    else f"{100 * (entry['median'] - previous) / previous:+.1f}%")
            lines.append(
                f"| {name} | {LADDER_TECHNIQUE.get(name, '')} "
                f"| {entry['median'] / 1e6:.1f} "
                f"| {entry['min'] / 1e6:.1f}-{entry['max'] / 1e6:.1f} "
                f"| {entry['median'] / anchor:.2f}x | {step} |")
            previous = entry["median"]
        lines.append("")
        nulls, previous = [], None
        for name in names[1:]:
            if previous is not None:
                gain = (entries[name]["median"] - previous) / previous
                if gain <= NULL_BAND:
                    nulls.append(name)
            previous = entries[name]["median"]
        if nulls:
            lines += ["**Null rungs: " + ", ".join(nulls) + ".** Levels 3 and 4"
                      " are null because of how the baseline was written:"
                      " BaselineSolver already looped outer-variance/"
                      "inner-stock and already peeled the v=0 row into its own"
                      " loop, so there was nothing left for those rungs to win"
                      " back. What they are worth is measured by the negative"
                      " controls below instead. A null anywhere else means the"
                      " compiler had already done it at -O2, which is the main"
                      " finding here: the only rung that pays is level 2, and"
                      " it is also the only one on the list the compiler is"
                      " not allowed to do for you, because turning"
                      " `x/(2*ds)` into `x*(1/(2*ds))` changes the answer."
                      " The task sheet gives credit for reported null results,"
                      " and L03 says they show methodology.", ""]

    # The page-size effect deserves its own row-by-row comparison.
    if ladders.get("A") and ladders.get("E"):
        lines += ["## The page-size effect", "",
                  "Same binary, same node, same flags, same nt. The only"
                  " difference is whether the two 8 MiB grid buffers are"
                  " backed by 2 MiB transparent huge pages or by 4 KiB pages."
                  " 16 MiB on 4 KiB pages needs 4096 TLB entries; on 2 MiB"
                  " pages it needs 8.", "",
                  "| solver | 4 KiB pages | huge pages | gain from pages |",
                  "|---|---|---|---|"]
        for name in LADDER_ORDER:
            a, e = ladders["A"].get(name), ladders["E"].get(name)
            if not (a and e):
                continue
            lines.append(f"| {name} | {a['median'] / 1e6:.1f} "
                         f"| {e['median'] / 1e6:.1f} "
                         f"| {e['median'] / a['median']:.2f}x |")
        lines.append("")

    if controls:
        top = control_reference * 1e6 if control_reference else None
        lines.append(f"*(controls measured in sweep {control_sweep}.)*")
        lines.append("")
        lines += [f"## Sweep {control_sweep}: negative controls", "",
                  "Deliberately worse code that does the same arithmetic as"
                  " the top of the ladder (test_opt_matches holds them to the"
                  " same tolerance). These are the measurements that give"
                  " levels 3 and 4 their meaning. Both are compared against"
                  " the best ladder rung measured in the same job at the same"
                  " nt, since a cross-job comparison would be measuring run"
                  " length rather than technique.", "",
                  "| control | what it ablates | Mcell-updates/s | "
                  "slowdown vs paired reference |", "|---|---|---|---|"]
        for name in CONTROLS:
            if name not in controls:
                continue
            entry = controls[name]
            slow = (f"{top / entry['median']:.1f}x" if top else "-")
            lines.append(f"| {name} | {LADDER_TECHNIQUE.get(name, '')} "
                         f"| {entry['median'] / 1e6:.1f} | {slow} |")
        lines.append("")

    if scaling:
        grids = sorted({g for (g, _) in scaling},
                       key=lambda g: int(g.split("x")[0]) * int(g.split("x")[1]))
        lines += ["## Sweep B: scaling (pinned re-run, job 556834)", "",
                  "MALLOC_MMAP_THRESHOLD_ is pinned, so every rep of every"
                  " row is served the same way. The rows still span two page"
                  " regimes, as physics rather than accident: a buffer"
                  " smaller than one 2 MiB huge page cannot be THP-backed,"
                  " so the two small grids measure on 4 KiB pages (their"
                  " throughput matches sweep A) while the 8 MiB and 32 MiB"
                  " buffers get huge pages (matching sweeps E and F)."
                  " Within-regime comparisons are the clean ones."
                  " RESULTS.md §3 carries the full reading.", "",
                  "| grid | working set | baseline (Mcell-updates/s) |"
                  " level 6 (Mcell-updates/s) | speedup |",
                  "|---|---|---|---|---|"]
        for grid in grids:
            b = scaling.get((grid, "baseline"))
            o = scaling.get((grid, "opt-L6"))
            if not (b and o):
                continue
            ns, nv = (int(x) for x in grid.split("x"))
            mib = 2 * 8 * ns * nv / 2**20
            lines.append(f"| {grid} | {mib:.0f} MiB | {b['median'] / 1e6:.1f} "
                         f"| {o['median'] / 1e6:.1f} "
                         f"| {o['median'] / b['median']:.2f}x |")
        lines.append("")

    if novec:
        lines += ["## Sweep C: auto-vectorisation control "
                  "(-O2 -fno-tree-vectorize)", "",
                  "Does the ladder's speedup survive with the compiler's"
                  " vectoriser switched off? If it does, the gains come from"
                  " the techniques and not from hidden auto-SIMD. Compared"
                  " against **sweep A**, because both ran with the default"
                  " allocator. Pairing it against the huge-page sweep would"
                  " credit the vectoriser with the page-size effect.", "",
                  "| solver | -O2 | -O2 -fno-tree-vectorize | vectoriser's"
                  " share |", "|---|---|---|---|"]
        for name in ("baseline", "opt-L6"):
            # Sweep C rebuilt and re-ran with the default allocator, so the
            # only sweep it can be paired with is A. Pairing it with sweep E
            # would hand the 1.66x page-size effect to the vectoriser.
            with_vec = ladders.get("A", {}).get(name)
            without = novec.get(name)
            if not (with_vec and without):
                continue
            share = with_vec["median"] / without["median"]
            lines.append(f"| {name} | {with_vec['median'] / 1e6:.1f} "
                         f"| {without['median'] / 1e6:.1f} | {share:.2f}x |")
        lines.append("")

    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("--out", default="results/")
    args = parser.parse_args()

    rows, provenance = read_csvs(args.csv)
    if not rows:
        print("no benchmark rows found. did the sbatch jobs finish?",
              file=sys.stderr)
        return 1

    # Refuse BEFORE anything is written, so a bad run cannot leave a
    # half-updated figure or summary behind for someone to quote.
    complaints = check_provenance(provenance)
    if complaints:
        print("FAIL: provenance check", file=sys.stderr)
        for complaint in complaints:
            print(f"  - {complaint}", file=sys.stderr)
        return 1

    out_dir = pathlib.Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Two ladders, measured with two different memory setups. Sweep A ran with
    # the default allocator, where glibc serves reps 2-5 out of the reused
    # heap on 4 KiB pages. Sweep E pinned MALLOC_MMAP_THRESHOLD_ so every rep
    # gets freshly mmap'd, transparent-huge-page-backed memory. A single
    # production solve is a fresh process, so sweep E is what a user sees.
    ladders = {
        "A": aggregate([r for r in rows if r.sweep == "A"], lambda r: r.solver),
        "E": aggregate([r for r in rows if r.sweep == "E"], lambda r: r.solver),
    }
    LADDER_LABEL = {
        "A": ("default allocator, 4 KiB pages", "bench_ladder.png"),
        "E": ("huge pages (MALLOC_MMAP_THRESHOLD_ pinned)",
              "bench_ladder_hugepage.png"),
    }

    # The controls and their paired reference have to come from one sweep,
    # measured in one job at one nt. Prefer E (full ladder and both controls
    # together), then D (the dedicated controls job), then A.
    control_sweep = next((sw for sw in ("E", "D", "A")
                          if any(r.sweep == sw and r.solver in CONTROLS
                                 for r in rows)), None)
    controls, control_reference = {}, None
    if control_sweep:
        controls = aggregate([r for r in rows if r.sweep == control_sweep
                              and r.solver in CONTROLS], lambda r: r.solver)
        paired = aggregate([r for r in rows if r.sweep == control_sweep
                            and r.solver.startswith("opt-L")],
                           lambda r: r.solver)
        control_reference = max((e["median"] / 1e6 for e in paired.values()),
                               default=None)

    scaling = aggregate([r for r in rows if r.sweep == "B"],
                        lambda r: (r.grid, r.solver))
    novec = aggregate([r for r in rows if r.sweep == "C"], lambda r: r.solver)

    print(f"parsed {len(rows)} rows: "
          f"ladderA={len(ladders['A'])} ladderE={len(ladders['E'])} "
          f"controls={len(controls)} (sweep {control_sweep}) "
          f"scaling={len(scaling)} novec={len(novec)}")

    written = []
    for sweep, ladder in ladders.items():
        if not ladder:
            continue
        regime, filename = LADDER_LABEL[sweep]
        figure = plot_ladder(ladder, controls, control_reference, provenance,
                             out_dir / filename, regime)
        if figure:
            written.append(figure)
    figure = plot_scaling(scaling, provenance, out_dir / "bench_scaling.png")
    if figure:
        written.append(figure)
    written.append(write_summary(ladders, controls, control_reference,
                                 control_sweep, scaling, novec, provenance,
                                 out_dir / "bench_summary.md"))

    for path in written:
        print(f"wrote {path}")

    # Exit non-zero as a gate. A figure nobody can trace back to a compute
    # node isn't a result.
    if not any(ladders.values()):
        print("FAIL: no ladder rows in any sweep, nothing to plot",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
