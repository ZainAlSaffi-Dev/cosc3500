#!/usr/bin/env python3
"""xi sweep (reading B): the volatility smile flattening into Black-Scholes.

Black-Scholes says every strike on the same underlying trades at one single
volatility. The market says otherwise: quote every strike, back out the
volatility Black-Scholes would need to reproduce that quote, and you get a
curved, tilted line, the smile or skew. Heston produces that curve because
volatility is itself random (xi = vol-of-vol) and correlated with the stock
(rho).

So run the Heston solver at xi = 0.35, turn its prices into Black-Scholes
implied volatilities, and the skew shows up. Turn xi down to 0 and the curve
has to collapse onto the flat line sigma = sqrt(v0) = 20%, because with
frozen variance the model is Black-Scholes. That collapse is the animation,
and it shows that Heston contains Black-Scholes rather than competing
with it.

HOW A WHOLE SMILE COMES OUT OF ONE SOLVE
    The solver's final sheet gives V as a function of spot for one fixed
    strike, but a smile needs V across strikes at one fixed spot. Two ways
    round that:

      (i)  homogeneity. A European call is homogeneous of degree 1 in
           (S, K): scale both and the price scales,  V(lam*S, lam*K) =
           lam*V(S, K). So one solve gives every strike:
                V(spot, K') = (spot/S_i) * V(S_i, K)   with  K' = spot*K/S_i
      (ii) re-run the solver once per strike.

    This script uses (i), and here it holds for a stronger reason than the
    usual continuum argument. In the continuum the Heston PDE's S terms are
    all of the form S*d/dS or S^2*d^2/dS^2, which don't change under
    S -> lam*S, and the payoff max(S-K,0) is homogeneous, so homogeneity is
    exact. In this solver it is exact at the discrete level too: the grid
    spans [0, s_max_mult*K] with ns nodes, so ds = s_max_mult*K/(ns-1) scales
    with K, and the stencil terms reduce to
        0.5*v*S^2*V_SS -> 0.5*v*i^2*(E - 2V + W)     (ds cancels)
        (r-q)*S*V_S    -> (r-q)*i*(E - W)/2          (ds cancels)
    i.e. the discrete operator depends only on the node index, never on ds.
    Rescaling the strike rescales the grid onto itself.

    --verify N checks this against reading (ii) by re-running the solver with
    option.strike = K' for N strikes and comparing prices. Because node index
    m maps to strike spot*K/(m*ds), the re-run's spot lands exactly on node m
    as well, so the comparison carries no interpolation error at all.

WHERE THE INVERSION IS TRUSTWORTHY
    Outside a band around the money a call is either worthless or pure
    intrinsic value; vega collapses, and dividing a fraction of a cent of
    grid error by a near-zero vega produces garbage implied vols. The script
    inverts only over --window (default 0.7K to 1.4K) and masks any strike
    whose price sits at or below intrinsic rather than crashing. The xi = 0
    end of the sweep has a known answer of exactly 20%, so it also reports
    the measured deviation from 20% across that window, which says how far
    the window can reasonably be pushed.

STABILITY
    dt = T/nt must stay under
        dt_stable ~ 1/(v_max*(ns-1)^2 + xi^2*(nv-1)^2/v_max + r).
    xi only appears in the middle term, so xi = 0 relaxes the bound and one
    nt sized for the largest xi is safe for the whole sweep. Defaults:
    ns=526, nv=101, v_max=0.2 give bound 1.633e-5 at xi = 0.35 and 1.814e-5
    at xi = 0; nt = 24000 gives dt = 1.042e-5, ratio 0.64 at the worst end.

Usage:
    arch -arm64 .venv/bin/python xi_sweep_smile.py --out xi_smile.mp4
    arch -arm64 .venv/bin/python xi_sweep_smile.py --steps 30 --verify 5
    arch -arm64 .venv/bin/python xi_sweep_smile.py --skip-solve   # reuse dumps
"""

import argparse
import math
import os
import shutil
import subprocess
import sys

import numpy as np

PROJECT = "/Users/zer0/Documents/cosc3500/project"


# ------------------------------------------------- Black-Scholes + invert ----

def norm_cdf(x):
    """Standard normal CDF via erfc, same identity as src/black_scholes.cpp."""
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def bs_call(spot, strike, rate, div_yield, sigma, tau):
    """Closed-form Black-Scholes call price (scalar; the root finder calls it
    thousands of times, so it stays plain Python floats)."""
    if tau <= 0.0 or sigma <= 0.0:
        return max(spot - strike, 0.0)
    vol_sqrt_tau = sigma * math.sqrt(tau)
    d1 = (math.log(spot / strike) +
          (rate - div_yield + 0.5 * sigma * sigma) * tau) / vol_sqrt_tau
    d2 = d1 - vol_sqrt_tau
    return (spot * math.exp(-div_yield * tau) * norm_cdf(d1) -
            strike * math.exp(-rate * tau) * norm_cdf(d2))


def brent(f, lo, hi, tol=1e-12, max_iter=200):
    """Brent's method: bisection's guaranteed bracketing plus the speed of
    inverse quadratic interpolation. Returns NaN if [lo, hi] doesn't bracket
    a sign change, and the caller treats that as "no implied vol here"
    instead of an exception.

    (scipy.optimize.brentq would do this in one line, but the project's
    virtualenv has numpy/matplotlib/imageio only, so it is written out.)
    """
    eps = np.finfo(float).eps
    a, b = lo, hi
    fa, fb = f(a), f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        return float("nan")
    c, fc = a, fa
    d = e = b - a
    for _ in range(max_iter):
        if fb * fc > 0.0:          # c must stay on the opposite side of the root
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):      # keep b as the best estimate so far
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol1 = 2.0 * eps * abs(b) + 0.5 * tol
        half_gap = 0.5 * (c - b)
        if abs(half_gap) <= tol1 or fb == 0.0:
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            # Interpolate, using the secant method when only two points are
            # distinct and inverse quadratic when three are.
            s = fb / fa
            if a == c:
                p, q = 2.0 * half_gap * s, 1.0 - s
            else:
                q, r = fa / fc, fb / fc
                p = s * (2.0 * half_gap * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                q = -q
            p = abs(p)
            if 2.0 * p < min(3.0 * half_gap * q - abs(tol1 * q), abs(e * q)):
                e, d = d, p / q            # accept the interpolated step
            else:
                d = e = half_gap           # reject it, bisect instead
        else:
            d = e = half_gap
        a, fa = b, fb
        b += d if abs(d) > tol1 else math.copysign(tol1, half_gap)
        fb = f(b)
    return b


def implied_vol(price, spot, strike, rate, div_yield, tau,
                lo=1e-6, hi=5.0):
    """Volatility that makes Black-Scholes reproduce `price`, or NaN.

    NaN (not an exception) is returned whenever no such volatility exists:
      * price at or below the no-arbitrage floor
        max(S*e^-qT - K*e^-rT, 0). sigma -> 0 already gives that value, so
        there is nothing left for volatility to explain;
      * price at or above the ceiling S*e^-qT, which even infinite volatility
        cannot exceed.
    Deep in/out of the money that floor is essentially the whole price, which
    is exactly why the strike window has to be limited.
    """
    floor = max(spot * math.exp(-div_yield * tau) -
                strike * math.exp(-rate * tau), 0.0)
    ceiling = spot * math.exp(-div_yield * tau)
    # This leaves a small relative cushion, because prices this close to the
    # floor carry no usable volatility information on any finite grid.
    if not (floor + 1e-10 * max(1.0, spot) < price < ceiling):
        return float("nan")
    return brent(lambda sigma: bs_call(spot, strike, rate, div_yield, sigma,
                                       tau) - price, lo, hi)


# ------------------------------------------------------- solver plumbing ----

CFG_TEMPLATE = """\
# Written by xi_sweep_smile.py, one rung of the vol-of-vol ladder.
# The Heston parameters have no CLI flags, so a .cfg is the only way to set
# xi; this file is regenerated for every rung of the sweep.
option.strike         = {strike!r}
option.maturity_years = {maturity!r}
option.type           = call

market.spot      = {spot!r}
market.rate      = {rate!r}
market.div_yield = {div_yield!r}

heston.v0    = {v0!r}
heston.kappa = {kappa!r}
heston.theta = {theta!r}
heston.xi    = {xi!r}
heston.rho   = {rho!r}

grid.ns         = {ns:d}
grid.nv         = {nv:d}
grid.nt         = {nt:d}
grid.s_max_mult = {s_max_mult!r}
grid.v_max      = {v_max!r}
"""


def run_solver(binary, cfg_path, dump_every=0, dump_dir=None):
    """Run ./heston once; return its single stdout CSV line."""
    cmd = [binary, "--config", cfg_path]
    if dump_every > 0:
        cmd += ["--dump-every", str(dump_every), "--dump-dir", dump_dir]
    done = subprocess.run(cmd, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"solver failed:\n{done.stderr}")
    return done.stdout.strip()


def read_variance_row(path, row_index):
    """One row (one variance level) out of a snapshot CSV.

    Snapshot layout is fixed by dump_snapshot() in src/io.cpp: nv rows
    (variance), ns columns (stock node). Only the v = v0 row is ever needed,
    so reading that single line beats np.loadtxt on the whole sheet.
    """
    with open(path) as handle:
        for index, line in enumerate(handle):
            if index == row_index:
                return np.array(line.split(","), dtype=float)
    raise SystemExit(f"{path}: no row {row_index}")


# ------------------------------------------------------------- the sweep ----

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", default=os.path.join(PROJECT, "heston"))
    parser.add_argument("--work", default="xi_work",
                        help="scratch dir for the generated .cfg files and dumps")
    parser.add_argument("--out", default="xi_smile.mp4")
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--hold", type=int, default=12,
                        help="extra frames held on the first and last rung so "
                             "a viewer can read them")
    parser.add_argument("--gif", action="store_true")
    parser.add_argument("--steps", type=int, default=30,
                        help="how many xi values from --xi-max down to 0")
    parser.add_argument("--xi-max", type=float, default=0.35)
    parser.add_argument("--skip-solve", action="store_true",
                        help="reuse the snapshots already in --work")
    parser.add_argument("--verify", type=int, default=0,
                        help="re-run the solver at N rescaled strikes to check "
                             "the homogeneity shortcut (0 = skip)")
    parser.add_argument("--refine-check", action="store_true",
                        help="re-solve the largest-xi rung with ds halved, and "
                             "again with the variance axis doubled, to measure "
                             "how far the recovered smile is grid-converged")
    # Contract and market, from config/reference.cfg.
    parser.add_argument("--strike", type=float, default=5250.0)
    parser.add_argument("--maturity", type=float, default=0.25)
    parser.add_argument("--spot", type=float, default=5200.0)
    parser.add_argument("--rate", type=float, default=0.045)
    parser.add_argument("--div-yield", type=float, default=0.013)
    # Model. v0 = theta so that xi = 0 gives a variance frozen at v0 and the
    # collapse target is the clean flat line sigma = sqrt(v0).
    parser.add_argument("--v0", type=float, default=0.04)
    parser.add_argument("--kappa", type=float, default=1.5)
    parser.add_argument("--theta", type=float, default=0.04)
    parser.add_argument("--rho", type=float, default=-0.70)
    # Grid.
    parser.add_argument("--ns", type=int, default=526)
    parser.add_argument("--nv", type=int, default=101)
    parser.add_argument("--nt", type=int, default=24000)
    parser.add_argument("--s-max-mult", type=float, default=4.0)
    parser.add_argument("--v-max", type=float, default=0.2)
    parser.add_argument("--window", type=float, nargs=2, default=(0.7, 1.4),
                        metavar=("LO", "HI"),
                        help="strike window to invert over, as multiples of "
                             "the base strike")
    args = parser.parse_args()

    flat_sigma = math.sqrt(args.v0)          # the xi = 0 answer, known exactly
    stock_spacing = args.s_max_mult * args.strike / (args.ns - 1)
    variance_spacing = args.v_max / (args.nv - 1)

    # v0 must sit exactly on a variance node, or the row we read is not the
    # row the price lives on.
    v0_row = int(round(args.v0 / variance_spacing))
    if abs(v0_row * variance_spacing - args.v0) > 1e-12 * max(1.0, args.v0):
        raise SystemExit(f"v0 = {args.v0} is not on a variance node "
                         f"(dv = {variance_spacing:g}), adjust nv/v_max")

    # One nt for the whole sweep, sized for the largest xi. xi only appears
    # in the middle term, so every smaller xi is strictly safer.
    dt = args.maturity / args.nt
    def dt_stable(xi):
        return 1.0 / (args.v_max * (args.ns - 1) ** 2 +
                      xi * xi * (args.nv - 1) ** 2 / args.v_max + args.rate)
    print(f"grid {args.ns}x{args.nv}x{args.nt}   ds = {stock_spacing:g}   "
          f"dv = {variance_spacing:g}   v0 on row {v0_row}")
    print(f"dt = {dt:.4e}")
    print(f"  dt_stable(xi={args.xi_max:g}) = {dt_stable(args.xi_max):.4e}  "
          f"ratio {dt / dt_stable(args.xi_max):.3f}   <- the binding case")
    print(f"  dt_stable(xi=0)    = {dt_stable(0.0):.4e}  "
          f"ratio {dt / dt_stable(0.0):.3f}")
    if dt >= dt_stable(args.xi_max):
        raise SystemExit("dt is over the stability bound at the largest xi")

    os.makedirs(args.work, exist_ok=True)
    xi_values = np.linspace(args.xi_max, 0.0, args.steps)

    # --- solve one rung per xi -------------------------------------------
    rows = []
    for rung, xi in enumerate(xi_values):
        dump_dir = os.path.join(args.work, f"xi_{rung:02d}")
        snap = os.path.join(dump_dir, f"snap_{args.nt:06d}.csv")
        if not args.skip_solve:
            cfg_path = os.path.join(args.work, f"xi_{rung:02d}.cfg")
            with open(cfg_path, "w") as handle:
                handle.write(CFG_TEMPLATE.format(
                    strike=args.strike, maturity=args.maturity, spot=args.spot,
                    rate=args.rate, div_yield=args.div_yield, v0=args.v0,
                    kappa=args.kappa, theta=args.theta, xi=float(xi),
                    rho=args.rho, ns=args.ns, nv=args.nv, nt=args.nt,
                    s_max_mult=args.s_max_mult, v_max=args.v_max))
            shutil.rmtree(dump_dir, ignore_errors=True)
            # Dumping every nt steps writes exactly one sheet, the final one.
            line = run_solver(args.binary, cfg_path, args.nt, dump_dir)
            print(f"  xi = {xi:.4f}  price at (spot, v0) = {line.split(',')[0]}")
        rows.append(read_variance_row(snap, v0_row))

    # --- prices across strikes, by homogeneity ---------------------------
    # Node i of the base grid carries V(S_i, K), and homogeneity turns that
    # into the price of a different strike at today's spot, since
    #     V(spot, K') = (spot/S_i) * V(S_i, K),   K' = spot*K/S_i
    # Only interior nodes are used, because node 0 and the last node hold
    # boundary values that were imposed rather than solved for.
    stock_axis = np.linspace(0.0, args.s_max_mult * args.strike, args.ns)
    node_index = np.arange(1, args.ns - 1)
    equiv_strike = args.spot * args.strike / stock_axis[node_index]
    scale = args.spot / stock_axis[node_index]
    keep = ((equiv_strike >= args.window[0] * args.strike) &
            (equiv_strike <= args.window[1] * args.strike))
    node_index, equiv_strike, scale = (node_index[keep], equiv_strike[keep],
                                       scale[keep])
    order = np.argsort(equiv_strike)         # plot low strike -> high strike
    node_index, equiv_strike, scale = (node_index[order], equiv_strike[order],
                                       scale[order])
    print(f"\nstrike window: K' in [{equiv_strike[0]:.0f}, "
          f"{equiv_strike[-1]:.0f}] = [{equiv_strike[0] / args.strike:.3f}K, "
          f"{equiv_strike[-1] / args.strike:.3f}K]  "
          f"({equiv_strike.size} strikes, one per stock node)")

    smiles = []
    for row in rows:
        prices = scale * row[node_index]
        smiles.append(np.array([
            implied_vol(float(p), args.spot, float(k), args.rate,
                        args.div_yield, args.maturity)
            for p, k in zip(prices, equiv_strike)]))
    smiles = np.array(smiles)

    # --- honesty check on the window -------------------------------------
    # The xi = 0 rung has a known exact answer, because implied vol is flat at
    # sqrt(v0) there. Whatever deviation shows up is the pipeline's own
    # numerical error, which calibrates how wide the window can reasonably be.
    flat = smiles[-1]
    flat_dev = np.abs(flat - flat_sigma)
    print(f"\n=== window calibration from the xi = 0 rung "
          f"(exact answer: {100 * flat_sigma:.4f}% at every strike) ===")
    print(f"  recovered implied vol range : [{100 * np.nanmin(flat):.4f}%, "
          f"{100 * np.nanmax(flat):.4f}%]")
    print(f"  max deviation from {100 * flat_sigma:.0f}%      : "
          f"{100 * np.nanmax(flat_dev):.4f} vol points "
          f"(at K' = {equiv_strike[int(np.nanargmax(flat_dev))]:.0f} "
          f"= {equiv_strike[int(np.nanargmax(flat_dev))] / args.strike:.2f}K)")
    for tolerance in (0.001, 0.002, 0.005):
        good = np.flatnonzero(flat_dev <= tolerance)
        if good.size:
            print(f"  within {100 * tolerance:.1f} vol points over K' in "
                  f"[{equiv_strike[good[0]] / args.strike:.3f}K, "
                  f"{equiv_strike[good[-1]] / args.strike:.3f}K]")
        else:
            print(f"  within {100 * tolerance:.1f} vol points: nowhere")
    masked = int(np.sum(~np.isfinite(smiles)))
    print(f"  strikes masked (no implied vol exists): {masked} "
          f"of {smiles.size}")

    # --- homogeneity verification against per-strike re-runs -------------
    if args.verify > 0:
        print(f"\n=== homogeneity check: {args.verify} strikes re-solved from "
              "scratch at xi = %.4f ===" % xi_values[0])
        picks = np.unique(np.linspace(0, node_index.size - 1,
                                      args.verify).astype(int))
        worst = 0.0
        for pick in picks:
            m = int(node_index[pick])
            strike_prime = float(equiv_strike[pick])
            predicted = float(scale[pick] * rows[0][m])
            cfg_path = os.path.join(args.work, f"verify_{m:03d}.cfg")
            with open(cfg_path, "w") as handle:
                handle.write(CFG_TEMPLATE.format(
                    strike=strike_prime, maturity=args.maturity,
                    spot=args.spot, rate=args.rate, div_yield=args.div_yield,
                    v0=args.v0, kappa=args.kappa, theta=args.theta,
                    xi=float(xi_values[0]), rho=args.rho, ns=args.ns,
                    nv=args.nv, nt=args.nt, s_max_mult=args.s_max_mult,
                    v_max=args.v_max))
            direct = float(run_solver(args.binary, cfg_path).split(",")[0])
            rel = abs(predicted - direct) / max(abs(direct), 1e-300)
            worst = max(worst, rel)
            print(f"  K' = {strike_prime:9.3f} ({strike_prime / args.strike:.3f}K)"
                  f"  homogeneity {predicted:14.8f}   re-run {direct:14.8f}"
                  f"   rel {rel:.2e}")
        print(f"  worst relative disagreement: {worst:.2e}  "
              "(the solver prints only 10 significant figures, so ~1e-10 is "
              "exact agreement)")

    # --- is the SMILE itself grid-converged? -----------------------------
    # The xi = 0 calibration above only proves the pipeline reproduces a flat
    # 20%. The wings at large xi are a harder test, because prices out there
    # are fractions of a cent, so a fixed dollar of grid error buys a lot of
    # implied vol. Two independent refinements answer that, one halving ds and
    # the other doubling the variance axis to see whether truncating v_max at
    # 0.2 is biting.
    if args.refine_check:
        print(f"\n=== smile convergence at xi = {xi_values[0]:.4f} ===")
        variants = []
        fine_ns = 2 * args.ns - 1
        fine_nt = int(math.ceil(args.maturity /
                                (0.6 / (args.v_max * (fine_ns - 1) ** 2 +
                                        args.xi_max ** 2 * (args.nv - 1) ** 2 /
                                        args.v_max + args.rate)) / 1000)) * 1000
        variants.append(("ds/2", fine_ns, args.nv, fine_nt, args.v_max))
        tall_nv = 2 * args.nv - 1                 # same dv, twice the v range
        tall_v_max = 2.0 * args.v_max
        tall_nt = int(math.ceil(args.maturity /
                                (0.6 / (tall_v_max * (args.ns - 1) ** 2 +
                                        args.xi_max ** 2 * (tall_nv - 1) ** 2 /
                                        tall_v_max + args.rate)) / 1000)) * 1000
        variants.append(("v_max*2", args.ns, tall_nv, tall_nt, tall_v_max))

        base_smile = smiles[0]
        for tag, r_ns, r_nv, r_nt, r_v_max in variants:
            cfg_path = os.path.join(args.work, f"refine_{tag.replace('/', '')}.cfg")
            with open(cfg_path, "w") as handle:
                handle.write(CFG_TEMPLATE.format(
                    strike=args.strike, maturity=args.maturity, spot=args.spot,
                    rate=args.rate, div_yield=args.div_yield, v0=args.v0,
                    kappa=args.kappa, theta=args.theta, xi=float(xi_values[0]),
                    rho=args.rho, ns=r_ns, nv=r_nv, nt=r_nt,
                    s_max_mult=args.s_max_mult, v_max=r_v_max))
            dump = os.path.join(args.work, f"refine_{tag.replace('/', '')}")
            shutil.rmtree(dump, ignore_errors=True)
            line = run_solver(args.binary, cfg_path, r_nt, dump)
            r_row = read_variance_row(os.path.join(dump, f"snap_{r_nt:06d}.csv"),
                                      int(round(args.v0 / (r_v_max / (r_nv - 1)))))
            r_axis = np.linspace(0.0, args.s_max_mult * args.strike, r_ns)
            # The refined grid has different nodes, so the equivalent stock
            # prices have to be interpolated onto it (linear is ample at 2x ds).
            equiv_stock = args.spot * args.strike / equiv_strike
            r_prices = (args.spot / equiv_stock) * np.interp(equiv_stock,
                                                             r_axis, r_row)
            r_smile = np.array([
                implied_vol(float(p), args.spot, float(k), args.rate,
                            args.div_yield, args.maturity)
                for p, k in zip(r_prices, equiv_strike)])
            drift = np.abs(r_smile - base_smile)
            inner = ((equiv_strike >= 0.85 * args.strike) &
                     (equiv_strike <= 1.25 * args.strike))
            print(f"  {tag:8s} ns={r_ns} nv={r_nv} nt={r_nt} v_max={r_v_max}  "
                  f"price={line.split(',')[0]}")
            print(f"    smile drift vs base: worst {100 * np.nanmax(drift):.3f} "
                  f"vol points over the full window, "
                  f"{100 * np.nanmax(drift[inner]):.3f} over 0.85K-1.25K")
            worst_at = equiv_strike[int(np.nanargmax(drift))]
            print(f"    worst at K' = {worst_at:.0f} "
                  f"= {worst_at / args.strike:.2f}K")

    # --- headline smile numbers ------------------------------------------
    print("\n=== smile per rung (implied vol over the window) ===")
    print(f"{'xi':>8}{'min IV':>10}{'max IV':>10}{'spread':>10}"
          f"{'IV @ K=spot':>13}{'skew 0.9K->1.1K':>18}")
    atm_col = int(np.argmin(np.abs(equiv_strike - args.spot)))
    lo_col = int(np.argmin(np.abs(equiv_strike - 0.9 * args.strike)))
    hi_col = int(np.argmin(np.abs(equiv_strike - 1.1 * args.strike)))
    spreads = np.nanmax(smiles, axis=1) - np.nanmin(smiles, axis=1)
    for rung, xi in enumerate(xi_values):
        smile = smiles[rung]
        if rung % max(1, args.steps // 10) and rung not in (0, args.steps - 1):
            continue
        print(f"{xi:8.4f}{100 * np.nanmin(smile):9.3f}%"
              f"{100 * np.nanmax(smile):9.3f}%{100 * spreads[rung]:9.3f}%"
              f"{100 * smile[atm_col]:12.3f}%"
              f"{100 * (smile[hi_col] - smile[lo_col]):17.3f}%")

    # --- render -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import imageio.v3 as iio

    plt.style.use("dark_background")

    iv_lo = float(np.nanmin(smiles))
    iv_hi = float(np.nanmax(smiles))
    pad = 0.12 * (iv_hi - iv_lo) or 0.01

    fig, (ax_smile, ax_spread) = plt.subplots(
        1, 2, figsize=(9.6, 4.8), dpi=100,
        gridspec_kw={"width_ratios": [2.15, 1.0], "wspace": 0.28})
    fig.subplots_adjust(left=0.075, right=0.975, top=0.855, bottom=0.20)
    fig.text(0.5, 0.015,
             "Every point is a Black-Scholes implied volatility backed out of a Heston PDE price. "
             "Black-Scholes can only draw a flat line;\n"
             "Heston's randomness in the volatility itself is what bends it. Turn that randomness off and Heston becomes Black-Scholes.",
             ha="center", va="bottom", fontsize=7.5, color="0.62")

    # xi runs high -> low, so the trace panel is drawn right-to-left.
    frames = []
    order_of_play = ([0] * args.hold + list(range(args.steps)) +
                     [args.steps - 1] * args.hold)
    for position, rung in enumerate(order_of_play):
        xi = float(xi_values[rung])
        ax_smile.clear()
        ax_spread.clear()

        # Keeping ghosts of the rungs already shown leaves a visible fan
        # behind the collapse instead of a single line jumping about.
        for older in range(rung):
            ax_smile.plot(equiv_strike, 100 * smiles[older], color="#ffb000",
                          linewidth=0.6, alpha=0.10)
        ax_smile.plot(equiv_strike, 100 * smiles[rung], color="#ffb000",
                      linewidth=2.2, label="Heston implied volatility",
                      zorder=3)
        # Drawn on top of the Heston line (higher zorder) so that on the last
        # rung the dashes show through it. That way the viewer sees the two
        # curves coincide instead of one just hiding the other.
        ax_smile.axhline(100 * flat_sigma, color="#7fd4ff", linestyle="--",
                         linewidth=1.3, zorder=4,
                         label=f"Black-Scholes: one flat vol = "
                               f"{100 * flat_sigma:.0f}%")
        ax_smile.axvline(args.spot, color="0.45", linestyle=":", linewidth=0.9)
        ax_smile.text(args.spot, 100 * (iv_hi + 0.55 * pad),
                      f" spot {args.spot:g}", color="0.55", fontsize=7.5,
                      ha="left", va="top")
        ax_smile.set_xlim(equiv_strike[0], equiv_strike[-1])
        ax_smile.set_ylim(100 * (iv_lo - pad), 100 * (iv_hi + pad))
        ax_smile.set_xlabel("strike K", fontsize=9)
        ax_smile.set_ylabel("Black-Scholes implied volatility (%)", fontsize=9)
        ax_smile.legend(loc="upper right", fontsize=8, framealpha=0.25)
        ax_smile.set_title("the smile", fontsize=9.5, color="0.8")

        ax_spread.plot(xi_values[: rung + 1], 100 * spreads[: rung + 1],
                       color="#ffb000", linewidth=1.8)
        ax_spread.plot(xi, 100 * spreads[rung], "o", color="#ffb000",
                       markersize=6)
        ax_spread.set_xlim(args.xi_max * 1.03, -0.01)   # sweep runs right->left
        ax_spread.set_ylim(-0.05 * 100 * spreads.max(), 110 * spreads.max())
        ax_spread.set_xlabel("xi  (volatility of volatility)", fontsize=9)
        ax_spread.set_ylabel("smile spread: max IV − min IV (%)", fontsize=9)
        ax_spread.set_title("how curved the smile is", fontsize=9.5,
                            color="0.8")
        ax_spread.text(0.96, 0.92, f"{100 * spreads[rung]:.2f}%",
                       transform=ax_spread.transAxes, ha="right", va="top",
                       fontsize=13, color="#ffb000", fontweight="bold")

        tail = "   <- Black-Scholes" if xi == 0.0 else ""
        fig.suptitle(f"Heston contains Black-Scholes:  vol-of-vol xi = "
                     f"{xi:.4f}{tail}", fontsize=11.5, y=0.965)
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
