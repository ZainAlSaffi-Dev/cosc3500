#!/usr/bin/env python3
"""Weather map (PLAN P6): results/<run>/snap_*.csv -> animated value surface.

Each snapshot CSV has nv rows (variance) and ns columns (stock price). That
layout is fixed by dump_snapshot in src/io.cpp. The animation runs in solver
order, from the expiry payoff back to today.

Styles and colour modes (pick whichever suits what the data has to show):
  --style flat (default)  heatmap over the (S, v) plane.
  --style surface         3-D sheet V(S, v) with a slow camera drift. The
                          "spreadsheet" itself, with the kink smoothing out.
  sequential (default)    dark background + inferno (black->red->yellow):
                          low values recede and heat glows. Use when
                          magnitude is the point.
  --diverging             white background + red/blue around 0 for the
                          instability clip, where the sign is the point.
                          NaN/inf cells render black.
  --time-value            subtract the expiry payoff and plot V - payoff.
                          The time value is the part that actually diffuses,
                          and raw V wastes the colour range on the linear
                          deep-ITM ramp. (Call payoff; strike via --strike.)

Making the motion visible (the reference parameters hide it, see below):
  --zoom                  crop both axes to the box that actually holds the
                          action. With s_max = 4K and v_max = 1.0 but a
                          three-month option at 20% vol, the live region is
                          a few percent of the frame; everything else is a
                          flat dead zone that never changes.
  --s-range LO HI         explicit stock-axis crop (overrides --zoom).
  --v-range LO HI         explicit variance-axis crop (overrides --zoom).
  --gamma G               colour compression, e.g. 0.5. The heat grows like
                          sqrt(tau), so a linear scale leaves the first half
                          of the clip nearly black on a fixed scale.
  --warp sqrt             pick frames so that sqrt(tau) advances evenly
                          rather than tau. Diffusion widths grow like
                          sqrt(tau), so this keeps the apparent motion
                          steady instead of front-loaded. Needs a dense
                          dump to make any difference.
  --frames N              how many frames to render after warping.
  --envelope              overlay the +/-1 sigma diffusion cone
                          S = K*exp(+/- sqrt(v*tau)), i.e. the region whose
                          outcome is still undecided. It flares open frame
                          by frame, which is the clearest cue of motion in
                          the whole clip.
  --price-trace           inset line plot of the quoted price at (spot, v0)
                          accumulating as the solver walks back to today,
                          plus a live numeric readout. The numbers move even
                          when the pixels look still.

Pass --config so the axes can't lie: strike, maturity, spot, v0, the model
parameters and the grid extents then all come from the same .cfg the solver
read. Without it the flags fall back to reference.cfg's values, and a dump
made with a different v_max gets mislabelled without any warning.

Usage:
    python3 scripts/weather_map.py results/smooth --out results/weather_map.mp4
    python3 scripts/weather_map.py results/smooth --time-value --out tv.mp4
    python3 scripts/weather_map.py results/smooth --style surface --out s.mp4
    python3 scripts/weather_map.py results/blowup --diverging --max-frames 90
    python3 scripts/weather_map.py results/demo --config config/demo.cfg \
        --time-value --zoom --gamma 0.5 --cap 99.9 --warp sqrt --frames 80 \
        --fps 12 --envelope --price-trace --out results/demo.mp4
"""

import argparse
import glob
import os
import re

import numpy as np

S_MAX_MULT = 4.0  # matches every config: S spans [0, 4*strike]


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


def warp_frame_indices(count: int, wanted: int, warp: str):
    """Which snapshot indices to render, in order.

    'linear' is even spacing in solver time. 'sqrt' spaces the frames evenly
    in sqrt(tau) instead: a diffusion front spreads like sqrt(tau), so equal
    steps in sqrt(tau) are equal steps of visible motion. Duplicates are
    dropped, so a sparse dump silently yields fewer frames than asked for.
    """
    if wanted <= 0 or wanted >= count:
        wanted = count
    fractions = np.linspace(0.0, 1.0, wanted)
    if warp == "sqrt":
        fractions = fractions ** 2  # tau ~ k^2  =>  sqrt(tau) ~ k
    picks = np.unique(np.round(fractions * (count - 1)).astype(int))
    return list(picks)


def read_config(path):
    """Parse the same "key = value" .cfg the C++ reads (see src/params.cpp).

    Every axis label, marker and zoom window in this script depends on the
    grid the dump was actually produced with. Typing those in by hand is how
    you end up with a frame claiming v runs to 1.0 when the solver stopped at
    0.2, so the default is to read them out of the config instead.
    """
    values = {}
    with open(path) as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, value = (part.strip() for part in line.split("=", 1))
            values[key] = value
    return values


def apply_config_defaults(args, parser):
    """Fill unset CLI flags from --config. Explicit flags always win."""
    cfg = read_config(args.config)
    # (cfg key, argparse dest), only the fields that shape the picture.
    mapping = [("option.strike", "strike"),
               ("option.maturity_years", "maturity"),
               ("market.spot", "spot"),
               ("heston.v0", "v0"),
               ("heston.kappa", "kappa"),
               ("heston.theta", "theta"),
               ("heston.xi", "xi"),
               ("grid.s_max_mult", "s_max_mult"),
               ("grid.v_max", "v_max")]
    for cfg_key, dest in mapping:
        if cfg_key not in cfg:
            continue
        # parser.get_default tells us whether the user actually typed the
        # flag; if they did, their value stands.
        if getattr(args, dest) == parser.get_default(dest):
            setattr(args, dest, float(cfg[cfg_key]))
    print(f"axes from {args.config}: strike={args.strike:g} "
          f"T={args.maturity:g} s_max={args.s_max_mult * args.strike:g} "
          f"v_max={args.v_max:g}")


def cir_variance_sd(v0, kappa, theta, xi, years):
    """Standard deviation of the Heston variance v_T, in closed form.

    The variance process is a CIR (Cox-Ingersoll-Ross) process, and its
    distribution at a future time is known exactly:
        Var(v_T) = v0*xi^2/kappa * (e^-kT - e^-2kT)
                 + theta*xi^2/(2*kappa) * (1 - e^-kT)^2
    We only need it to answer "which variance levels can the market
    plausibly reach?", which is what decides how much of the v axis is worth
    filming.
    """
    if kappa <= 0.0 or years <= 0.0:
        return 0.0
    decay = np.exp(-kappa * years)
    var = (v0 * xi * xi / kappa * (decay - decay * decay) +
           theta * xi * xi / (2.0 * kappa) * (1.0 - decay) ** 2)
    return float(np.sqrt(max(var, 0.0)))


def reachable_window(args, stock_axis, variance_axis):
    """The (S, v) box the market can plausibly reach, i.e. the shot worth
    filming.

    The grid deliberately runs way past anything realistic (v up to 1.0, so
    100% vol, and S up to 4x the strike) so the artificial boundaries sit far
    away from the answer. That is right for the numerics and awful for film:
    nearly the whole frame shows scenarios the model never visits, and a
    fixed colour scale computed over that dead zone flattens everything that
    does move.

    So crop to +/-3 standard deviations of where the two processes actually
    go over the option's life: variance from the CIR formula above, stock
    price from the lognormal spread exp(+/-3*sqrt(v_hi*T)) around the strike.
    Nothing is recomputed or faked, this only picks a viewport.
    """
    sd_v = cir_variance_sd(args.v0, args.kappa, args.theta, args.xi,
                           args.maturity)
    v_hi = args.v0 + 3.0 * sd_v
    # When xi is zero the variance is frozen, and this keeps a band around v0
    # in shot rather than collapsing the crop to a single row.
    v_hi = max(v_hi, 3.0 * args.v0, variance_axis[1])
    v_hi = min(v_hi, float(variance_axis[-1]))

    sigma_log = np.sqrt(max(v_hi, 0.0) * args.maturity)
    s_lo = args.strike * np.exp(-3.0 * sigma_log)
    s_hi = args.strike * np.exp(+3.0 * sigma_log)
    # Today's market always has to be in shot, since it's the cell we quote.
    s_lo = min(s_lo, args.spot * 0.95)
    s_hi = max(s_hi, args.spot * 1.05)
    return (max(float(stock_axis[0]), s_lo), min(float(stock_axis[-1]), s_hi),
            float(variance_axis[0]), v_hi)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot_dir")
    parser.add_argument("--config", default=None,
                        help="the .cfg the dump was produced with; supplies "
                             "strike/maturity/spot/v0/kappa/theta/xi/"
                             "s_max_mult/v_max so the axes cannot be "
                             "mislabelled. Explicit flags still win.")
    parser.add_argument("--out", default="results/weather_map.mp4")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--style", choices=("flat", "surface"), default="flat")
    parser.add_argument("--diverging", action="store_true",
                        help="two-hue scale around 0 for the instability clip")
    parser.add_argument("--time-value", action="store_true",
                        help="plot V - expiry payoff (call) instead of V")
    parser.add_argument("--strike", type=float, default=5250.0,
                        help="strike for --time-value payoff and the S axis")
    parser.add_argument("--maturity", type=float, default=0.25,
                        help="option maturity in years, used to turn solver "
                             "steps into days-to-expiry for the titles")
    parser.add_argument("--spot", type=float, default=5200.0,
                        help="today's stock price, marked on flat frames")
    parser.add_argument("--v0", type=float, default=0.04,
                        help="today's variance, marked on flat frames")
    # The three model parameters below never touch the pixels directly. They
    # tell --zoom how far the variance process can plausibly wander, which
    # decides how much of the v axis is worth showing.
    parser.add_argument("--kappa", type=float, default=1.5,
                        help="heston.kappa the dump was produced with (--zoom)")
    parser.add_argument("--theta", type=float, default=0.04,
                        help="heston.theta the dump was produced with (--zoom)")
    parser.add_argument("--xi", type=float, default=0.35,
                        help="heston.xi the dump was produced with (--zoom)")
    parser.add_argument("--s-max-mult", type=float, default=S_MAX_MULT,
                        help="grid.s_max_mult the dump was produced with; "
                             "the S axis runs [0, s_max_mult*strike]")
    parser.add_argument("--v-max", type=float, default=1.0,
                        help="grid.v_max the dump was produced with")
    parser.add_argument("--zoom", action="store_true",
                        help="auto-crop both axes to the live region")
    parser.add_argument("--s-range", type=float, nargs=2, default=None,
                        metavar=("LO", "HI"), help="explicit S-axis crop")
    parser.add_argument("--v-range", type=float, nargs=2, default=None,
                        metavar=("LO", "HI"), help="explicit v-axis crop")
    parser.add_argument("--gamma", type=float, default=1.0,
                        help="colour power-law: <1 lifts small values so the "
                             "early, faint frames are readable")
    parser.add_argument("--warp", choices=("linear", "sqrt"), default="linear",
                        help="frame spacing: even in tau, or even in sqrt(tau)")
    parser.add_argument("--frames", type=int, default=0,
                        help="frames to render after warping (0 = all)")
    parser.add_argument("--envelope", action="store_true",
                        help="overlay the +/-1 sigma cone S = K*exp(+/-sqrt(v*tau))")
    parser.add_argument("--price-trace", action="store_true",
                        help="inset: quoted price at (spot, v0) vs time")
    parser.add_argument("--caption", default=None,
                        help="one/two-line explainer under the axes "
                             "(default: auto text per style)")
    parser.add_argument("--gif", action="store_true",
                        help="also write a .gif next to the .mp4")
    parser.add_argument("--max-frames", type=int, default=0,
                        help="render only the first N snapshots (0 = all); "
                             "the blow-up clip's tail is all-NaN frames")
    parser.add_argument("--cap", type=float, default=98.0,
                        help="percentile that sets the top of the colour "
                             "scale. 98 keeps the late high-v glow from "
                             "compressing early structure; in --time-value "
                             "mode 99.9 is better (98 clips the final frame).")
    parser.add_argument("--title", default="Heston option value V(S, v)")
    args = parser.parse_args()
    if args.config:
        apply_config_defaults(args, parser)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import imageio.v3 as iio

    # Sequential mode is the one that looks like film, so it gets a dark stage
    # for the heat to glow against. Diverging mode keeps a white stage,
    # because its midpoint is light by design and calm regions ought to look
    # quiet rather than lit up.
    if not args.diverging:
        plt.style.use("dark_background")

    snaps = load_snapshots(args.snapshot_dir)
    # The total step count has to come from the full dump list. Working it
    # out after --max-frames has truncated things would mislabel every
    # remaining frame's step fraction and days-to-expiry.
    total_steps = max(step for step, _ in snaps)
    if args.max_frames > 0:
        snaps = snaps[: args.max_frames]
    # Frames get selected before the colour scale is fixed, so the scale
    # describes exactly the frames that end up being rendered.
    picks = warp_frame_indices(len(snaps), args.frames, args.warp)
    if len(picks) < len(snaps):
        print(f"frame selection: {len(picks)} of {len(snaps)} snapshots "
              f"({args.warp} spacing)")
    snaps = [snaps[i] for i in picks]

    s_max = args.s_max_mult * args.strike
    num_stock = snaps[0][1].shape[1]
    stock_axis = np.linspace(0.0, s_max, num_stock)
    variance_axis = np.linspace(0.0, args.v_max, snaps[0][1].shape[0])

    # The price we quote is one cell of the raw sheet; capture it before any
    # time-value subtraction or cropping can move the indices around.
    spot_col = int(np.argmin(np.abs(stock_axis - args.spot)))
    v0_row = int(np.argmin(np.abs(variance_axis - args.v0)))
    price_trace = [float(vals[v0_row, spot_col]) for _, vals in snaps]

    value_label = "option value ($)"
    if args.time_value:
        # Time value = V - intrinsic. The deep-ITM ramp cancels out and the
        # colour range is spent entirely on what the PDE actually moves.
        payoff = np.maximum(stock_axis - args.strike, 0.0)
        snaps = [(step, vals - payoff[np.newaxis, :]) for step, vals in snaps]
        value_label = "time value V − payoff ($)"

    # Crop to the live region. This is not cosmetic, because a three-month
    # option at 20% vol leaves roughly 97% of the full frame permanently
    # unchanged, and a fixed colour scale computed over that dead zone
    # flattens everything that does move.
    if args.s_range or args.v_range or args.zoom:
        if args.zoom and not (args.s_range and args.v_range):
            auto = reachable_window(args, stock_axis, variance_axis)
        else:
            auto = (stock_axis[0], stock_axis[-1],
                    variance_axis[0], variance_axis[-1])
        s_lo, s_hi = args.s_range if args.s_range else auto[0:2]
        v_lo, v_hi = args.v_range if args.v_range else auto[2:4]
        s_keep = np.flatnonzero((stock_axis >= s_lo) & (stock_axis <= s_hi))
        v_keep = np.flatnonzero((variance_axis >= v_lo) & (variance_axis <= v_hi))
        if s_keep.size < 2 or v_keep.size < 2:
            raise SystemExit("crop window keeps fewer than 2 nodes on an axis")
        stock_axis = stock_axis[s_keep]
        variance_axis = variance_axis[v_keep]
        snaps = [(step, vals[np.ix_(v_keep, s_keep)]) for step, vals in snaps]
        print(f"cropped to S in [{stock_axis[0]:.0f}, {stock_axis[-1]:.0f}], "
              f"v in [{variance_axis[0]:.4f}, {variance_axis[-1]:.4f}] "
              f"({len(v_keep)}x{len(s_keep)} cells)")
        # A crop down to a handful of nodes shows the grid rather than the
        # physics. The fix is a config whose v_max isn't mostly unreachable,
        # not smoother interpolation here.
        if len(v_keep) < 24 or len(s_keep) < 24:
            print("  warning: few nodes survive the crop, the frame will "
                  "look blocky. Re-run the solver with a v_max/s_max_mult "
                  "closer to the reachable region (see config/demo.cfg).")

    # One colour scale across every frame, otherwise the animation "breathes"
    # and the instability frames clip. Non-finite values (blow-up) are masked.
    finite_vals = np.concatenate(
        [a[np.isfinite(a)].ravel() for _, a in snaps])
    if finite_vals.size == 0:
        raise SystemExit("every cell is NaN/inf, nothing to render")
    if args.diverging:
        # When the sign is what matters the limits are symmetric around zero
        # with a neutral midpoint. Using the 99.5th percentile rather than the
        # maximum stops one saturated cell washing out every earlier frame.
        limit = np.percentile(np.abs(finite_vals), 99.5) or 1.0
        vmin, vmax, cmap_name, bad_colour = -limit, limit, "RdBu_r", "black"
    else:
        # When magnitude is what matters a perceptually uniform heat ramp is
        # used instead. Capping at the 98th percentile stops the late glow at
        # high variance from compressing the early structure near the strike,
        # and any hotter cells simply saturate to yellow.
        vmin = min(0.0, float(finite_vals.min()))
        vmax = float(np.percentile(finite_vals, args.cap)) or 1.0
        cmap_name, bad_colour = "inferno", "magenta"
    cmap = matplotlib.colormaps[cmap_name].copy()
    cmap.set_bad(bad_colour)  # broken cells should be impossible to miss

    # A power-law colour scale buys back the early frames. Heat grows like
    # sqrt(tau), so on a linear scale the first half of the clip sits in the
    # bottom few percent of the ramp and just looks black.
    norm = None
    if args.gamma != 1.0:
        if args.diverging:
            raise SystemExit("--gamma is for the sequential scale, not "
                             "--diverging (a power law has no sign)")
        # PowerNorm needs a non-negative floor.
        norm = matplotlib.colors.PowerNorm(args.gamma, vmin=max(0.0, vmin),
                                           vmax=vmax)
    # imshow/plot_surface take either (vmin, vmax) or a norm, never both.
    colour_kwargs = ({"norm": norm} if norm is not None
                     else {"vmin": vmin, "vmax": vmax})

    stock_mesh, variance_mesh = np.meshgrid(stock_axis, variance_axis)

    # Nobody thinks in solver steps, so convert to calendar time. Step n
    # holds the sheet n*dt before expiry, and the last step is "today".
    # (Assumes the final dump lands on the final step, true whenever
    # dump_every divides nt.)
    days_total = args.maturity * 365.0
    days_per_step = days_total / total_steps

    frames = []
    if args.style == "surface":
        fig = plt.figure(figsize=(8, 4.5), dpi=120)
        ax = fig.add_subplot(projection="3d")
    else:
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)

    # Inset price ticker, so there is at least one obviously moving element.
    # The heatmap can look static frame to frame while the quoted price
    # climbs from 0 to its value today, and this draws that climb.
    # It sits low and to the right, clear of the colour bar and of the hot
    # near-strike column the heat map lights up.
    trace_ax = fig.add_axes([0.545, 0.33, 0.165, 0.20]) if args.price_trace \
        else None

    # A static explainer for someone who has never seen Heston, saying what
    # the axes mean and which way the solver moves. It is drawn once at figure
    # level so that the per-frame ax.clear() leaves it alone.
    caption = args.caption
    if caption is None and args.style == "surface":
        caption = ("The sheet of what-if option values over stock price and "
                   "variance. The kinked expiry payoff is known exactly;\n"
                   "the solver smooths it backwards through time until it "
                   "reaches today's sheet.")
    elif caption is None:
        caption = ("Each pixel is one what-if scenario: stock price S "
                   "(x-axis) under variance v (y-axis); colour = the "
                   "option's value there.\nThe expiry payoff is known "
                   "exactly; the solver fills the sheet backwards in time "
                   "until it reaches today.")
    fig.subplots_adjust(bottom=0.22 if args.style == "flat" else 0.14)
    fig.text(0.5, 0.02, caption, ha="center", va="bottom", fontsize=7.5,
             color="0.65" if not args.diverging else "0.35")
    for index, (step, grid_vals) in enumerate(snaps):
        ax.clear()
        masked = np.ma.masked_invalid(grid_vals)
        broken = int(np.sum(~np.isfinite(grid_vals)))
        note = f"  [{broken} broken cells]" if broken else ""
        # How far back from expiry this frame sits. Drives the title and the
        # width of the diffusion cone.
        tau_years = args.maturity * step / total_steps
        if args.style == "surface":
            ax.plot_surface(stock_mesh, variance_mesh,
                            np.where(np.isfinite(grid_vals), grid_vals, 0.0),
                            cmap=cmap, rcount=64, ccount=128, linewidth=0,
                            antialiased=False, **colour_kwargs)
            # A little headroom, or the tallest frame's peak gets sliced off
            # flat by the top of the axes box and reads as a plateau.
            ax.set_zlim(vmin, 1.06 * float(finite_vals.max()))
            ax.set_zlabel(value_label, fontsize=8)
            # The camera drifts about 25 degrees over the whole clip.
            ax.view_init(elev=25,
                         azim=-60 + 25.0 * index / max(1, len(snaps) - 1))
        else:
            image = ax.imshow(masked, origin="lower", aspect="auto",
                              cmap=cmap, extent=(stock_axis[0], stock_axis[-1],
                                                 variance_axis[0],
                                                 variance_axis[-1]),
                              **colour_kwargs)
            if index == 0:
                fig.colorbar(image, ax=ax, label=value_label, pad=0.02)
            # Two landmarks for a reader who does not know the contract. One
            # is the strike, where the payoff kinks, and the other is today's
            # market point, the single cell that becomes the quoted price.
            accent = "0.75" if not args.diverging else "0.35"
            # Labels sit on top of the hottest pixels in the late frames, so
            # every one gets a semi-opaque plate to stay readable.
            plate = dict(facecolor="black" if not args.diverging else "white",
                         alpha=0.55, edgecolor="none", pad=1.6)
            ax.axvline(args.strike, color=accent, linestyle="--",
                       linewidth=1.0)
            ax.text(args.strike + 0.012 * (stock_axis[-1] - stock_axis[0]),
                    variance_axis[0] + 0.93 * (variance_axis[-1] -
                                               variance_axis[0]),
                    f"strike K = {args.strike:g}", color=accent, fontsize=8,
                    bbox=plate)
            if args.envelope and tau_years > 0.0:
                # The +/-1 sigma cone. At variance v, a lognormal stock
                # wanders sqrt(v*tau) in log-price over the remaining life,
                # so the scenarios whose moneyness is still genuinely in
                # doubt are the ones inside K*exp(+/-sqrt(v*tau)). The cone
                # flares open as the solver walks back, which is the clearest
                # single cue that time is passing.
                sigma_log = np.sqrt(np.maximum(variance_axis, 0.0) * tau_years)
                for sign in (+1.0, -1.0):
                    ax.plot(args.strike * np.exp(sign * sigma_log),
                            variance_axis, color="#7fd4ff", linewidth=1.1,
                            alpha=0.85)
                ax.text(0.985, 0.055,
                        "cyan: ±1σ cone, where the outcome is still in doubt",
                        transform=ax.transAxes, ha="right", fontsize=7,
                        color="#7fd4ff", bbox=plate)
            ax.plot(args.spot, args.v0, "o", markersize=8,
                    markerfacecolor="none", markeredgecolor=accent,
                    markeredgewidth=1.5)
            ax.annotate(f"today's market (S = {args.spot:g}, v = {args.v0:g})"
                        "\n= the price we quote",
                        xy=(args.spot, args.v0),
                        xytext=(0.045, 0.80), textcoords="axes fraction",
                        fontsize=8, color=accent, bbox=plate,
                        arrowprops=dict(arrowstyle="-", color=accent,
                                        linewidth=0.8))
        if trace_ax is not None:
            # Redraw the trace so far, with a dot marking where it has got to.
            trace_ax.clear()
            trace_ax.set_facecolor("black" if not args.diverging else "white")
            trace_ax.patch.set_alpha(0.72)
            trace_ax.tick_params(labelsize=6, length=2)
            trace_ax.set_title("price at today's market ($)", fontsize=6.5,
                               pad=3)
            trace_ax.set_xlabel("days before expiry", fontsize=6, labelpad=1)
            days_axis = [s * days_per_step for s, _ in snaps]
            # x runs high to low, since the clip starts at expiry and walks
            # back to today.
            trace_ax.set_xlim(days_total, 0.0)
            span = (max(price_trace) - min(price_trace)) or 1.0
            trace_ax.set_ylim(min(price_trace) - 0.05 * span,
                              max(price_trace) + 0.25 * span)
            trace_ax.plot(days_axis[: index + 1], price_trace[: index + 1],
                          color="#ffb000", linewidth=1.4)
            trace_ax.plot(days_axis[index], price_trace[index], "o",
                          color="#ffb000", markersize=4)
            trace_ax.text(0.05, 0.82, f"${price_trace[index]:,.2f}",
                          transform=trace_ax.transAxes, fontsize=9,
                          color="#ffb000", fontweight="bold")
        ax.set_xlabel("stock price S", fontsize=9)
        ax.set_ylabel("variance v", fontsize=9)
        # Headline in calendar time, with the solver step on a second smaller
        # line for anyone cross-referencing the code.
        days_left = step * days_per_step
        ax.set_title(f"{args.title}: {days_left:.0f} of {days_total:.0f} "
                     f"days before expiry\n"
                     f"(solver step {step}/{total_steps}){note}",
                     fontsize=9)
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
        frames.append(frame)
    plt.close(fig)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    # H.264 mp4 (course submission format); imageio-ffmpeg ships ffmpeg.
    iio.imwrite(args.out, frames, fps=args.fps, codec="libx264")
    print(f"wrote {args.out} ({len(frames)} frames @ {args.fps} fps)")
    if args.gif:
        gif_path = os.path.splitext(args.out)[0] + ".gif"
        iio.imwrite(gif_path, frames, duration=1000 / args.fps, loop=0)
        print(f"wrote {gif_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
