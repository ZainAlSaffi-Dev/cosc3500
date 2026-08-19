#!/usr/bin/env python3
"""Model collapse (reading A): watch the Heston PDE turn into Black-Scholes.

The Heston PDE has three extra terms compared with Black-Scholes, and every
one of them carries a factor of xi (the vol-of-vol):

    0.5*xi^2*v*V_vv        variance diffuses
    rho*xi*v*S*V_Sv        stock and variance move together
    kappa*(theta - v)*V_v  variance is pulled back to its long-run level

Set xi = 0 and rho = 0 and the first two vanish outright. Start the variance
exactly at its long-run level, v0 = theta, and the third one vanishes too,
because its coefficient kappa*(theta - v) is exactly zero on that row. What
is left on the v = v0 row is

    V_t + 0.5*v0*S^2*V_SS + (r - q)*S*V_S - r*V = 0

which is the Black-Scholes equation with a constant volatility sigma =
sqrt(v0). So the whole v = v0 row of the solver's sheet has to reproduce the
closed-form Black-Scholes curve at every moment of solver time, not just at
the end. That is what this script films.

The animation runs in solver order (expiry -> today), same as
scripts/weather_map.py:
  top panel    the solver's slice V(S, v0, tau) drawn over the closed-form
               Black-Scholes curve for the same tau, plus the expiry payoff
               for reference
  bottom panel the difference in dollars, on a fixed scale so the viewer can
               see it is flat noise and not a drifting bias

Nothing here is fitted or tuned. The only inputs are the solver's own CSV
snapshots (src/io.cpp: row = variance level, column = stock node) and the
textbook Black-Scholes formula.

Grid notes (why the defaults are what they are):
  * v0 = theta = 0.04 has to land exactly on a variance node, otherwise the
    row we read isn't the row where the collapse is exact. v_max = 0.2 with
    nv = 101 gives dv = 0.002 and v0 = node 20. The script checks this.
  * spot = 5200 has to land exactly on a stock node or the quoted-price
    readout is off: s_max = 4*5250 = 21000 with ns = 526 gives ds = 40 and
    spot = node 130. The script checks this too.
  * dt has to stay under the explicit-scheme bound
        dt_stable ~ 1 / (v_max*(ns-1)^2 + xi^2*(nv-1)^2/v_max + r).
    With xi = 0 the middle term disappears, so the bound relaxes to
    1/(0.2*525^2 + 0.045) = 1.814e-5. nt = 24000 gives dt = 1.042e-5, a
    ratio of 0.57. The same nt is also safe at xi = 0.35 (bound 1.633e-5,
    ratio 0.64), which is why xi_sweep_smile.py can reuse it unchanged.

Usage:
    arch -arm64 .venv/bin/python bs_collapse_anim.py --out collapse.mp4
    arch -arm64 .venv/bin/python bs_collapse_anim.py --frames 120 --gif
    arch -arm64 .venv/bin/python bs_collapse_anim.py --skip-solve   # reuse dumps
"""

import argparse
import glob
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

# Where the compiled solver and its .venv live, so the defaults just work.
PROJECT = "/Users/zer0/Documents/cosc3500/project"


# ---------------------------------------------------------------- maths ----

def norm_cdf(x):
    """Standard normal CDF, vectorised. Same identity as src/black_scholes.cpp:
    N(x) = erfc(-x/sqrt(2))/2, which is exact rather than an approximation."""
    from math import sqrt
    # math.erfc is scalar-only, so use the numpy-friendly rescaling of erf.
    # np has no erf either without scipy, so fall back to a vectorised erfc
    # built from the same math.erfc, applied elementwise.
    return 0.5 * np.vectorize(math.erfc)(-np.asarray(x, dtype=float) / sqrt(2.0))


def bs_call(spot, strike, rate, div_yield, sigma, tau):
    """Closed-form Black-Scholes call, vectorised over `spot`.

    tau is time remaining to expiry (the solver's tau), so tau = 0 is expiry
    and the formula degenerates to the payoff. Handled explicitly because
    d1/d2 divide by sigma*sqrt(tau).
    """
    spot = np.asarray(spot, dtype=float)
    if tau <= 0.0 or sigma <= 0.0:
        return np.maximum(spot - strike, 0.0)
    stock_leg = spot * math.exp(-div_yield * tau)
    strike_leg = strike * math.exp(-rate * tau)
    vol_sqrt_tau = sigma * math.sqrt(tau)
    # S = 0 makes log(S/K) = -inf. numpy would warn; sidestep it and set the
    # value by hand, since a call on a worthless stock is worthless.
    safe_spot = np.where(spot > 0.0, spot, 1.0)
    d1 = (np.log(safe_spot / strike) +
          (rate - div_yield + 0.5 * sigma * sigma) * tau) / vol_sqrt_tau
    d2 = d1 - vol_sqrt_tau
    price = stock_leg * norm_cdf(d1) - strike_leg * norm_cdf(d2)
    return np.where(spot > 0.0, price, 0.0)


# ------------------------------------------------------- solver plumbing ----

CFG_TEMPLATE = """\
# Written by bs_collapse_anim.py. This sits in the corner of the Heston
# parameter space where the model collapses to Black-Scholes, which needs the
# vol of vol at zero, rho at zero and v0 equal to theta.
# There are no CLI flags for the Heston parameters, so a .cfg is the only
# way to set them; this file is regenerated on every run.
option.strike         = {strike!r}
option.maturity_years = {maturity!r}
option.type           = call

market.spot      = {spot!r}
market.rate      = {rate!r}
market.div_yield = {div_yield!r}

heston.v0    = {v0!r}
heston.kappa = {kappa!r}
heston.theta = {v0!r}
heston.xi    = 0.0
heston.rho   = 0.0

grid.ns         = {ns:d}
grid.nv         = {nv:d}
grid.nt         = {nt:d}
grid.s_max_mult = {s_max_mult!r}
grid.v_max      = {v_max!r}
"""


def divisor_near(total, target):
    """Largest-but-closest divisor of `total` to `target`.

    The final snapshot must land exactly on step nt, otherwise the last frame
    is not "today" and every frame's tau label is off. Snapshots only happen
    when n % dump_every == 0, so dump_every has to divide nt.
    """
    target = max(1, int(round(target)))
    best = 1
    for candidate in range(1, total + 1):
        if total % candidate == 0 and abs(candidate - target) < abs(best - target):
            best = candidate
    return best


def run_solver(binary, cfg_path, dump_every, dump_dir):
    """Run ./heston and return (price, stderr_line). stdout is the one CSV line."""
    cmd = [binary, "--config", cfg_path,
           "--dump-every", str(dump_every), "--dump-dir", dump_dir]
    done = subprocess.run(cmd, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"solver failed:\n{done.stderr}")
    return done.stdout.strip(), done.stderr.strip()


def read_variance_row(path, row_index):
    """Pull a single row out of a snapshot CSV.

    A snapshot is nv rows (variance) by ns columns (stock), a layout fixed
    by dump_snapshot() in src/io.cpp. The collapse only ever needs the
    v = v0 row, so reading the single line instead of np.loadtxt on the whole
    sheet is ~100x less work.
    """
    with open(path) as handle:
        for index, line in enumerate(handle):
            if index == row_index:
                return np.array(line.split(","), dtype=float)
    raise SystemExit(f"{path}: no row {row_index}")


# -------------------------------------------------------------- the film ----

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", default=os.path.join(PROJECT, "heston"))
    parser.add_argument("--work", default="collapse_work",
                        help="scratch dir for the generated .cfg and snapshots")
    parser.add_argument("--out", default="bs_collapse.mp4")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--gif", action="store_true",
                        help="also write a .gif next to the .mp4")
    parser.add_argument("--frames", type=int, default=100,
                        help="target number of snapshots (rounded to a divisor of nt)")
    parser.add_argument("--skip-solve", action="store_true",
                        help="reuse the snapshots already in --work")
    # Contract and market, taken from config/reference.cfg.
    parser.add_argument("--strike", type=float, default=5250.0)
    parser.add_argument("--maturity", type=float, default=0.25)
    parser.add_argument("--spot", type=float, default=5200.0)
    parser.add_argument("--rate", type=float, default=0.045)
    parser.add_argument("--div-yield", type=float, default=0.013)
    # v0 doubles as theta here, because the collapse needs the two to be equal
    # so that the mean-reversion term switches itself off on the row we read.
    parser.add_argument("--v0", type=float, default=0.04)
    parser.add_argument("--kappa", type=float, default=1.5)
    # Grid.
    parser.add_argument("--ns", type=int, default=526)
    parser.add_argument("--nv", type=int, default=101)
    parser.add_argument("--nt", type=int, default=24000)
    parser.add_argument("--s-max-mult", type=float, default=4.0)
    parser.add_argument("--v-max", type=float, default=0.2)
    # Where the comparison is scored. Outside this band a call is either
    # worthless or pure intrinsic value, and a relative error there says
    # nothing about the numerics.
    parser.add_argument("--window", type=float, nargs=2, default=(0.7, 1.4),
                        metavar=("LO", "HI"),
                        help="scoring/viewing band as multiples of the strike")
    parser.add_argument("--rel-floor", type=float, default=1.0,
                        help="only score RELATIVE error where the BS value "
                             "exceeds this many dollars")
    parser.add_argument("--refine-check", action="store_true",
                        help="also solve on a twice-finer stock grid to show "
                             "the residual gap is discretisation error")
    args = parser.parse_args()

    sigma = math.sqrt(args.v0)
    stock_spacing = args.s_max_mult * args.strike / (args.ns - 1)
    variance_spacing = args.v_max / (args.nv - 1)

    # --- node-alignment checks -------------------------------------------
    # Reading the collapse off a row that is not exactly v = theta would show
    # a real (not numerical) mismatch, so refuse rather than mislead.
    v0_row = int(round(args.v0 / variance_spacing))
    v0_offset = abs(v0_row * variance_spacing - args.v0)
    if v0_offset > 1e-12 * max(1.0, args.v0):
        raise SystemExit(
            f"v0 = {args.v0} is not on a variance node (dv = {variance_spacing:g}, "
            f"nearest node {v0_row} sits at {v0_row * variance_spacing:g}). "
            "Choose nv/v_max so v0/dv is a whole number.")
    spot_col = int(round(args.spot / stock_spacing))
    spot_offset = abs(spot_col * stock_spacing - args.spot)
    if spot_offset > 1e-9 * args.spot:
        print(f"warning: spot {args.spot} is {spot_offset:.3g} off node "
              f"{spot_col} (ds = {stock_spacing:g}), so the quoted-price "
              "readout carries an O(ds) offset error")

    # --- stability budget -------------------------------------------------
    # xi = 0 kills the variance-diffusion term, so the bound is looser here
    # than for the same grid at xi = 0.35. Printed so the ratio is on record.
    dt = args.maturity / args.nt
    dt_stable = 1.0 / (args.v_max * (args.ns - 1) ** 2 + args.rate)
    print(f"grid {args.ns}x{args.nv}x{args.nt}   ds = {stock_spacing:g}   "
          f"dv = {variance_spacing:g}")
    print(f"dt = {dt:.4e}   dt_stable(xi=0) ~ {dt_stable:.4e}   "
          f"ratio = {dt / dt_stable:.3f}")
    if dt >= dt_stable:
        raise SystemExit("dt is over the stability bound, the run would blow up")

    os.makedirs(args.work, exist_ok=True)
    cfg_path = os.path.join(args.work, "collapse.cfg")
    dump_dir = os.path.join(args.work, "dump")
    dump_every = divisor_near(args.nt, args.nt / max(1, args.frames))

    if not args.skip_solve:
        with open(cfg_path, "w") as handle:
            handle.write(CFG_TEMPLATE.format(
                strike=args.strike, maturity=args.maturity, spot=args.spot,
                rate=args.rate, div_yield=args.div_yield, v0=args.v0,
                kappa=args.kappa, ns=args.ns, nv=args.nv, nt=args.nt,
                s_max_mult=args.s_max_mult, v_max=args.v_max))
        shutil.rmtree(dump_dir, ignore_errors=True)
        print(f"solving (dump every {dump_every} of {args.nt} steps) ...")
        csv_line, status = run_solver(args.binary, cfg_path, dump_every, dump_dir)
        print("  " + status)
        print("  " + csv_line)
        quoted_price = float(csv_line.split(",")[0])
    else:
        quoted_price = None

    paths = sorted(glob.glob(os.path.join(dump_dir, "snap_*.csv")))
    if not paths:
        raise SystemExit(f"no snapshots in {dump_dir}")

    stock_axis = np.linspace(0.0, args.s_max_mult * args.strike, args.ns)
    payoff = np.maximum(stock_axis - args.strike, 0.0)

    # Frame 0 is the expiry payoff itself, the sheet the solver starts from
    # (Grid::init_payoff). The solver never dumps step 0, so it gets rebuilt
    # here from the same formula. It makes the clip open on the kink, which
    # is what the PDE then spends the whole animation smoothing out.
    steps = [0]
    slices = [payoff.copy()]
    for path in paths:
        steps.append(int(re.search(r"snap_(\d+)", path).group(1)))
        slices.append(read_variance_row(path, v0_row))
    print(f"loaded {len(paths)} snapshots (+ the synthetic tau = 0 payoff frame)")

    # --- scoring ----------------------------------------------------------
    win_lo, win_hi = args.window[0] * args.strike, args.window[1] * args.strike
    in_window = (stock_axis >= win_lo) & (stock_axis <= win_hi)
    print(f"scoring window: S in [{win_lo:.0f}, {win_hi:.0f}] "
          f"({int(in_window.sum())} nodes)")

    # A finite-difference grid cannot see the payoff kink until the diffusion
    # has smeared it over several cells. The smear width in stock terms is
    # sigma*sqrt(tau)*K, so the first frames, where that width is under a few
    # ds, carry a real and expected numerical error that has nothing to do
    # with the Heston-to-Black-Scholes question. The frames get split on that
    # boundary so the headline number isn't dominated by it.
    cells_needed = 4.0
    tau_resolved = (cells_needed * stock_spacing /
                    (args.strike * sigma)) ** 2
    print(f"payoff kink is resolved (smear >= {cells_needed:g} cells) once "
          f"tau >= {tau_resolved:.5f} y = {tau_resolved * 365:.1f} days")

    reference = []          # closed-form BS curve per frame
    taus = []
    max_abs_err = []        # $ , over the window
    max_rel_err = []        # |diff| / local BS value, where BS >= --rel-floor
    for step, pde in zip(steps, slices):
        tau = args.maturity * step / args.nt
        taus.append(tau)
        exact = bs_call(stock_axis, args.strike, args.rate, args.div_yield,
                        sigma, tau)
        reference.append(exact)
        diff = np.abs(pde - exact)
        max_abs_err.append(float(diff[in_window].max()))
        scorable = in_window & (exact >= args.rel_floor)
        max_rel_err.append(float((diff[scorable] / exact[scorable]).max())
                           if scorable.any() else 0.0)
    taus = np.array(taus)
    max_abs_err = np.array(max_abs_err)
    max_rel_err = np.array(max_rel_err)
    resolved = taus >= tau_resolved

    final_pde, final_bs = slices[-1], reference[-1]
    exact_at_spot = float(bs_call(np.array([args.spot]), args.strike, args.rate,
                                  args.div_yield, sigma, args.maturity)[0])
    pde_at_spot = float(final_pde[spot_col])
    # This is the scale-free number a practitioner would quote, being the
    # worst dollar disagreement anywhere in the window divided by the
    # at-the-money price.
    rel_to_atm = max_abs_err / exact_at_spot

    print("\n=== measured agreement, PDE slice vs closed-form Black-Scholes ===")
    print(f"sigma used in the closed form: sqrt(v0) = {sigma:.6f}")
    print(f"{'':44}{'ALL frames':>14}{'resolved only':>16}")
    print(f"{'max |PDE - BS| in window ($)':44}"
          f"{max_abs_err.max():14.6f}{max_abs_err[resolved].max():16.6f}")
    print(f"{'  ... as a fraction of the ATM price':44}"
          f"{rel_to_atm.max():14.3e}{rel_to_atm[resolved].max():16.3e}")
    print(f"{f'max |PDE - BS| / BS  (where BS >= ${args.rel_floor:g})':44}"
          f"{max_rel_err.max():14.3e}{max_rel_err[resolved].max():16.3e}")
    print(f"\nfinal frame (tau = T, i.e. today):")
    print(f"  max |PDE - BS| in window   : "
          f"${np.abs(final_pde - final_bs)[in_window].max():.6f}  "
          f"({rel_to_atm[-1]:.3e} of the ATM price)")
    print(f"  max |PDE - BS| / BS        : {max_rel_err[-1]:.3e}")
    print(f"  at today's market (S = {args.spot:g}): PDE {pde_at_spot:.7f}  "
          f"BS {exact_at_spot:.7f}  rel "
          f"{abs(pde_at_spot - exact_at_spot) / exact_at_spot:.3e}")
    if quoted_price is not None:
        print(f"  (solver's own quoted-price line: {quoted_price:.7f})")

    if args.refine_check:
        # Halving ds, and shrinking dt to stay stable, shows that what is left
        # is discretisation error rather than a genuine gap between models.
        refine_ns = 2 * args.ns - 1          # keeps every old node, adds midpoints
        refine_dt_stable = 1.0 / (args.v_max * (refine_ns - 1) ** 2 + args.rate)
        refine_nt = int(math.ceil(args.maturity / (0.6 * refine_dt_stable) / 1000)) * 1000
        print(f"\nrefinement check: ns {args.ns} -> {refine_ns}, "
              f"nt {args.nt} -> {refine_nt} (ds halved)")
        refine_cfg = os.path.join(args.work, "collapse_refine.cfg")
        with open(refine_cfg, "w") as handle:
            handle.write(CFG_TEMPLATE.format(
                strike=args.strike, maturity=args.maturity, spot=args.spot,
                rate=args.rate, div_yield=args.div_yield, v0=args.v0,
                kappa=args.kappa, ns=refine_ns, nv=args.nv, nt=refine_nt,
                s_max_mult=args.s_max_mult, v_max=args.v_max))
        refine_dir = os.path.join(args.work, "dump_refine")
        shutil.rmtree(refine_dir, ignore_errors=True)
        line, _ = run_solver(args.binary, refine_cfg, refine_nt, refine_dir)
        refine_row = read_variance_row(
            os.path.join(refine_dir, f"snap_{refine_nt:06d}.csv"), v0_row)
        refine_axis = np.linspace(0.0, args.s_max_mult * args.strike, refine_ns)
        refine_exact = bs_call(refine_axis, args.strike, args.rate,
                               args.div_yield, sigma, args.maturity)
        refine_win = (refine_axis >= win_lo) & (refine_axis <= win_hi)
        refine_err = float(np.abs(refine_row - refine_exact)[refine_win].max())
        coarse_err = float(np.abs(final_pde - final_bs)[in_window].max())
        print(f"  coarse max |PDE - BS| : ${coarse_err:.6f}")
        print(f"  refined max |PDE - BS|: ${refine_err:.6f}")
        print(f"  ratio                 : {coarse_err / refine_err:.2f}x "
              f"(4x would be clean second order in ds)")
        print(f"  refined quoted price  : {line.split(',')[0]}  vs BS "
              f"{exact_at_spot:.7f}")

    # --- render -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import imageio.v3 as iio

    plt.style.use("dark_background")   # same stage as scripts/weather_map.py

    view = np.flatnonzero(in_window)
    stock_view = stock_axis[view]
    value_top = 1.05 * float(reference[-1][view].max())
    # One fixed error scale is used for the whole clip, because rescaling per
    # frame would hide the fact that the error never grows.
    err_limit = 1.25 * float(max_abs_err.max()) or 1e-6

    fig, (ax_value, ax_err) = plt.subplots(
        2, 1, figsize=(8, 5.5), dpi=120, sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.35], "hspace": 0.12})
    fig.subplots_adjust(left=0.10, right=0.975, top=0.86, bottom=0.155)
    fig.text(0.5, 0.015,
             "xi = 0 and rho = 0 freeze the variance at v0 = theta, and the Heston PDE becomes Black-Scholes with sigma = sqrt(v0) = "
             f"{100 * sigma:.0f}%.\nThe solver is told none of this: it runs the full 2-D Heston stencil and lands on the closed form anyway, at every instant of solver time.",
             ha="center", va="bottom", fontsize=7.5, color="0.62")

    days_total = args.maturity * 365.0
    frames = []
    for index, (step, pde) in enumerate(zip(steps, slices)):
        exact = reference[index]
        tau = args.maturity * step / args.nt
        ax_value.clear()
        ax_err.clear()

        ax_value.plot(stock_view, payoff[view], color="0.42", linewidth=0.9,
                      linestyle=":", label="expiry payoff max(S - K, 0)")
        # Closed form drawn thick and pale underneath the solver's thin
        # bright line, so if they agree the orange sits inside the blue band.
        ax_value.plot(stock_view, exact[view], color="#7fd4ff", linewidth=4.5,
                      alpha=0.55, solid_capstyle="round",
                      label="closed-form Black-Scholes")
        ax_value.plot(stock_view, pde[view], color="#ffb000", linewidth=1.3,
                      label="Heston PDE slice V(S, v0, tau)")
        ax_value.axvline(args.strike, color="0.5", linestyle="--", linewidth=0.8)
        ax_value.set_ylim(-0.02 * value_top, value_top)
        ax_value.set_ylabel("option value ($)", fontsize=9)
        ax_value.legend(loc="upper left", fontsize=7.5, framealpha=0.25)
        ax_value.text(args.strike + 0.008 * (stock_view[-1] - stock_view[0]),
                      0.05 * value_top, f"strike K = {args.strike:g}",
                      color="0.6", fontsize=7.5)

        difference = pde[view] - exact[view]
        ax_err.axhline(0.0, color="0.45", linewidth=0.8)
        ax_err.fill_between(stock_view, 0.0, difference, color="#ff5c5c",
                            alpha=0.75, linewidth=0)
        ax_err.set_ylim(-err_limit, err_limit)
        ax_err.set_ylabel("PDE − BS ($)", fontsize=9)
        ax_err.set_xlabel("stock price S", fontsize=9)
        ax_err.text(0.985, 0.86,
                    f"max |error| in window: ${max_abs_err[index]:.4f}"
                    f"   ({1e4 * max_abs_err[index] / exact_at_spot:.2f} "
                    "basis points of the at-the-money price)",
                    transform=ax_err.transAxes, ha="right", va="top",
                    fontsize=7.5, color="#ff9a9a")
        if 0.0 < tau < tau_resolved:
            # The opening frames need saying out loud, because that spike is
            # the grid failing to resolve a corner rather than the two models
            # disagreeing with each other.
            ax_err.text(0.015, 0.86,
                        "kink still narrower than a few grid cells",
                        transform=ax_err.transAxes, ha="left", va="top",
                        fontsize=7.5, color="0.6")

        days_left = days_total * (1.0 - step / args.nt)
        ax_value.set_title(
            f"Heston with the volatility-of-volatility switched off:  "
            f"{tau * 365.0:5.1f} of {days_total:.0f} days before expiry\n"
            f"(solver step {step}/{args.nt};  {days_left:.0f} days of solving still to go)",
            fontsize=9.5)
        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
    plt.close(fig)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    iio.imwrite(args.out, frames, fps=args.fps, codec="libx264")
    print(f"\nwrote {args.out} ({len(frames)} frames @ {args.fps} fps)")
    if args.gif:
        gif_path = os.path.splitext(args.out)[0] + ".gif"
        iio.imwrite(gif_path, frames, duration=1000 / args.fps, loop=0)
        print(f"wrote {gif_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
