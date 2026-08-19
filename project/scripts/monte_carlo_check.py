#!/usr/bin/env python3
"""Monte-Carlo cross-check of the Heston PDE price by simulating the SDEs.

The C++ solver takes the deterministic route: it walks a partial differential
equation backwards from expiry to today. This script takes the stochastic
route instead: draw random paths for the stock and its variance, work out
what the option pays on each path, and average. Same model, two completely
different numerical machines. If both land on the same number, the PDE is
handling the model's randomness correctly.

The model (Heston 1993), under the risk-neutral measure:

    dS = (r - q) S dt + sqrt(v) S dW1        stock price
    dv = kappa (theta - v) dt + xi sqrt(v) dW2   its instantaneous variance
    corr(dW1, dW2) = rho

In plain English: the stock drifts at the risk-free rate minus the dividend
yield (that is what "risk-neutral" gets us, nobody's risk appetite appears),
and its wobbliness sqrt(v) is itself random. v is pulled back towards a long
run level theta at speed kappa, gets kicked around with strength xi, and its
kicks are correlated with the stock's kicks through rho. Negative rho is the
leverage effect: markets fall and volatility spikes at the same time.

Discretisation (both choices matter, see the comments at each site):
  * log-Euler for S: exact in the drift, and S = exp(...) can never go
    negative the way a naive Euler step can.
  * full-truncation Euler for v: clamps v at zero inside the coefficients.
    The reference parameters break the Feller condition (2 kappa theta
    = 0.12 < xi^2 = 0.1225), so v really does reach zero, and a scheme that
    is careless there gives NaNs or a biased price.

Variance reduction: antithetic variates (every random path is simulated
alongside its mirror image). Greeks: delta by bump-and-revalue with common
random numbers.

There is a third opinion available too: Heston is one of the few
stochastic-volatility models with a semi-analytic price (Fourier inversion of
its characteristic function). This script computes that as well, in numpy
with no extra dependencies, so if MC and the PDE disagree there is something
to settle it. Three routes, one number.

Usage:
    python3 monte_carlo_check.py --config config/smoke.cfg
    python3 monte_carlo_check.py --config config/smoke.cfg \
        --paths 200000 --steps 400 --pde-price 205.3242902
    python3 monte_carlo_check.py --config config/reference.cfg \
        --paths 1000000 --steps 800 --seed 7 --step-ladder 100,200,400,800
    python3 monte_carlo_check.py --config config/smoke.cfg --self-test
"""

import argparse
import math
import sys

import numpy as np

# 97.5th percentile of the standard normal -> two-sided 95% confidence band.
Z_95 = 1.959963984540054

# numpy renamed trapz -> trapezoid in 2.0; accept either so the script runs
# on whatever numpy the marking environment happens to have.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


# --------------------------------------------------------------------------
# Config files use the same "key = value" format the C++ reads, which is
# handled in src/params.cpp.
# --------------------------------------------------------------------------

def load_config(path: str) -> dict:
    """Parse the .cfg the C++ solver uses, so both routes price the same deal.

    Format is one 'key = value' per line, '#' starts a comment. Keys are the
    dotted names from src/params.cpp, e.g. option.strike, market.rate,
    heston.v0. Values are all numbers except option.type (call/put).
    """
    values = {}
    with open(path, "r") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if "=" not in line:
                raise SystemExit(f"{path}:{line_no}: expected 'key = value'")
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()

    def number(key: str) -> float:
        if key not in values:
            raise SystemExit(f"{path}: missing required key '{key}'")
        try:
            return float(values[key])
        except ValueError:
            raise SystemExit(f"{path}: '{key}' is not a number: {values[key]}")

    option_type = values.get("option.type", "call")
    if option_type not in ("call", "put"):
        raise SystemExit(f"{path}: option.type must be call or put")

    # The grid keys are deliberately ignored, because they describe the PDE's
    # mesh and Monte Carlo has its own knobs in --paths and --steps.
    return {
        "strike": number("option.strike"),
        "maturity": number("option.maturity_years"),
        "is_call": option_type == "call",
        "spot": number("market.spot"),
        "rate": number("market.rate"),
        "div_yield": number("market.div_yield"),
        "v0": number("heston.v0"),
        "kappa": number("heston.kappa"),
        "theta": number("heston.theta"),
        "xi": number("heston.xi"),
        "rho": number("heston.rho"),
    }


# --------------------------------------------------------------------------
# The simulation itself.
# --------------------------------------------------------------------------

def simulate_log_growth(cfg: dict, num_pairs: int, num_steps: int, seed: int):
    """Simulate the Heston SDEs; return log(S_T / S_0) for every path.

    Returns (growth, diagnostics) where growth has shape (2, num_pairs):
    row 0 is the ordinary path, row 1 is its antithetic mirror.

    Returning the growth rather than S_T itself: under log-Euler the spot
    only ever enters as an additive constant log(S_0), so
    S_T = S_0 * exp(growth). Keeping the growth lets the delta bump below
    re-price at a different spot without re-drawing a single random number,
    which is common random numbers in their purest form.
    """
    dt = cfg["maturity"] / num_steps
    sqrt_dt = math.sqrt(dt)
    kappa, theta, xi, rho = cfg["kappa"], cfg["theta"], cfg["xi"], cfg["rho"]
    drift = cfg["rate"] - cfg["div_yield"]

    # Cholesky factor of the 2x2 correlation matrix [[1, rho], [rho, 1]].
    # Its lower-triangular Cholesky factor L is [[1, 0], [rho, sqrt(1-rho^2)]],
    # so feeding L two independent standard normals (z1, z_indep) produces
    # (z1, rho*z1 + rho_bar*z_indep) with exactly correlation rho and unit
    # variances. This one line is where the "correlated random noise" of the
    # model actually enters the simulation.
    rho_bar = math.sqrt(1.0 - rho * rho)

    rng = np.random.default_rng(seed)

    # Both members of each antithetic pair are held at once, so these arrays
    # have two rows.
    # Row 0 uses the draws as-is, row 1 uses their negatives (see `sign`).
    log_s = np.zeros((2, num_pairs))          # log(S_t / S_0), starts at 0
    variance = np.full((2, num_pairs), cfg["v0"])
    sign = np.array([[1.0], [-1.0]])          # broadcasts down the rows

    zero_hits = 0                             # times v was clamped at 0
    first_step_corr = float("nan")

    for step in range(num_steps):
        z1_base = rng.standard_normal(num_pairs)
        z_indep_base = rng.standard_normal(num_pairs)

        # Antithetic variates mean running each draw and its negative. Both
        # are perfectly good standard normals, so the estimator stays
        # unbiased, but averaging the pair cancels most of the linear part of
        # the payoff's response to the noise, which roughly halves the error
        # bar for free.
        z1 = sign * z1_base
        z2 = rho * z1 + rho_bar * (sign * z_indep_base)

        if step == 0:
            first_step_corr = float(np.corrcoef(z1[0], z2[0])[0, 1])

        # Full truncation means every coefficient sees max(v, 0), while the
        # state v itself is allowed to go slightly negative and recover.
        # Feller is broken for the reference parameters, so v really does
        # hit zero. sqrt() of a raw negative v would be NaN, and simply
        # reflecting v to |v| would inject variance the model never had.
        v_plus = np.maximum(variance, 0.0)
        zero_hits += int(np.count_nonzero(variance < 0.0))
        vol = np.sqrt(v_plus)

        # log-Euler for the stock. Ito's lemma turns dS = (r-q) S dt +
        # sqrt(v) S dW1 into d(log S) = (r - q - v/2) dt + sqrt(v) dW1. The
        # The -v/2 is the Ito correction, which is there because the log of an
        # average is not the average of the log. Stepping the log rather than
        # the price keeps S strictly positive however large a draw arrives.
        log_s += (drift - 0.5 * v_plus) * dt + vol * sqrt_dt * z1

        # Euler-Maruyama for the variance, coefficients evaluated at v_plus.
        variance = (variance
                    + kappa * (theta - v_plus) * dt
                    + xi * vol * sqrt_dt * z2)

    total_states = 2 * num_pairs * num_steps
    diagnostics = {
        "zero_hit_fraction": zero_hits / total_states,
        "first_step_corr": first_step_corr,
        "final_variance_mean": float(np.mean(np.maximum(variance, 0.0))),
        "dt": dt,
    }
    return log_s, diagnostics


def price_from_growth(growth, cfg: dict, spot: float):
    """Discounted expected payoff, its standard error, and the plain-MC error.

    growth[0] and growth[1] are an antithetic pair, so one piece of
    independent evidence is a pair average, not an individual path. The
    confidence interval has to be built from the num_pairs pair averages;
    treating the 2*num_pairs paths as independent would understate the error.
    """
    terminal_spot = spot * np.exp(growth)
    if cfg["is_call"]:
        payoff = np.maximum(terminal_spot - cfg["strike"], 0.0)
    else:
        payoff = np.maximum(cfg["strike"] - terminal_spot, 0.0)

    discount = math.exp(-cfg["rate"] * cfg["maturity"])
    payoff *= discount

    pair_mean = 0.5 * (payoff[0] + payoff[1])
    num_pairs = pair_mean.size

    price = float(np.mean(pair_mean))
    std_err = float(np.std(pair_mean, ddof=1) / math.sqrt(num_pairs))

    # What the error bar would have been without antithetics, for the same
    # number of simulated paths. Shows the variance reduction is real.
    plain_err = float(np.std(payoff, ddof=1) / math.sqrt(payoff.size))
    return price, std_err, plain_err


def delta_by_bump(growth, cfg: dict, spot: float, bump_fraction: float):
    """Central-difference delta, re-priced on the same random paths.

    dPrice/dSpot ~= (P(S+h) - P(S-h)) / (2h). The whole trick is that both
    re-prices reuse `growth`, i.e. identical random numbers (common random
    numbers). Fresh draws on each side would leave two prices each carrying
    a few tenths of noise, and their difference divided by a small h would
    be pure garbage. With shared draws the noise cancels almost exactly.
    """
    bump = bump_fraction * spot
    price_up, err_up, _ = price_from_growth(growth, cfg, spot + bump)
    price_down, err_down, _ = price_from_growth(growth, cfg, spot - bump)
    delta = (price_up - price_down) / (2.0 * bump)
    # This error bar is crude on purpose. The two legs are strongly correlated
    # with each other, so it over-estimates rather than claiming a tight
    # statistical bound.
    delta_err = (err_up + err_down) / (2.0 * bump)
    return delta, delta_err


# --------------------------------------------------------------------------
# The third opinion, which is the semi-analytic Heston price in plain numpy.
# --------------------------------------------------------------------------

def semi_analytic_call(cfg: dict, num_nodes: int = 200000, u_max: float = 400.0):
    """Heston's closed-form-ish call price by Fourier inversion.

    Heston's trick: you cannot write down the distribution of S_T, but you CAN
    write down its characteristic function E[exp(i*u*log S_T)] in closed form.
    The price then falls out of two probability-like integrals,

        Call = S0 e^{-qT} P1  -  K e^{-rT} P2,
        Pj   = 1/2 + (1/pi) Int_0^inf Re[ e^{-i u log K} f_j(u) / (i u) ] du,

    which is just Black-Scholes' N(d1)/N(d2) structure with the two normal
    probabilities replaced by numerically integrated ones.

    Written in the Albrecher et al. "little Heston trap" form: the textbook
    version takes a complex logarithm that crosses its branch cut for longer
    maturities and silently returns nonsense. Flipping the sign convention
    (using g = num/den rather than den/num) keeps the log on its principal
    branch. This is a genuine trap, not a stylistic preference.

    Accuracy is set by the trapezoid grid, not by randomness: this is the
    referee, so it must be tighter than either competitor. Validate with
    --self-test, which drives xi -> 0 and checks the answer collapses onto
    Black-Scholes (no volatility-of-volatility means constant volatility).
    """
    spot, strike, maturity = cfg["spot"], cfg["strike"], cfg["maturity"]
    rate, div_yield = cfg["rate"], cfg["div_yield"]
    v0, kappa, theta = cfg["v0"], cfg["kappa"], cfg["theta"]
    xi, rho = cfg["xi"], cfg["rho"]

    log_spot = math.log(spot)
    mean_reversion_term = kappa * theta
    # This starts just off zero because the integrand has a removable
    # singularity at the origin.
    u = np.linspace(1e-10, u_max, num_nodes)
    imag = 1j

    def probability(j: int) -> float:
        # The two integrals differ only in these two constants.
        if j == 1:
            u_j, b_j = 0.5, kappa - rho * xi
        else:
            u_j, b_j = -0.5, kappa

        d = np.sqrt((rho * xi * imag * u - b_j) ** 2
                    - xi ** 2 * (2.0 * u_j * imag * u - u ** 2))
        numer = b_j - rho * xi * imag * u - d
        denom = b_j - rho * xi * imag * u + d
        g = numer / denom                       # "little trap" orientation
        decay = np.exp(-d * maturity)

        c_term = ((rate - div_yield) * imag * u * maturity
                  + (mean_reversion_term / xi ** 2)
                  * (numer * maturity
                     - 2.0 * np.log((1.0 - g * decay) / (1.0 - g))))
        d_term = (numer / xi ** 2) * (1.0 - decay) / (1.0 - g * decay)

        char_fn = np.exp(c_term + d_term * v0 + imag * u * log_spot)
        integrand = np.real(
            np.exp(-imag * u * math.log(strike)) * char_fn / (imag * u))
        return 0.5 + float(_trapezoid(integrand, u)) / math.pi

    call = (spot * math.exp(-div_yield * maturity) * probability(1)
            - strike * math.exp(-rate * maturity) * probability(2))
    if cfg["is_call"]:
        return call
    # Put-call parity gives P = C - S e^{-qT} + K e^{-rT}. This is the same
    # identity the C++ test suite checks, so a put config is covered too.
    return (call - spot * math.exp(-div_yield * maturity)
            + strike * math.exp(-rate * maturity))


def black_scholes_call(spot, strike, maturity, rate, div_yield, vol):
    """Closed-form BS call, used only by --self-test."""
    sqrt_t = math.sqrt(maturity)
    d1 = ((math.log(spot / strike)
           + (rate - div_yield + 0.5 * vol * vol) * maturity) / (vol * sqrt_t))
    d2 = d1 - vol * sqrt_t
    normal_cdf = lambda z: 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
    return (spot * math.exp(-div_yield * maturity) * normal_cdf(d1)
            - strike * math.exp(-rate * maturity) * normal_cdf(d2))


def run_self_test(cfg: dict):
    """Prove the semi-analytic referee before trusting it to judge anything.

    With xi -> 0 the variance stops moving, so Heston degenerates into
    Black-Scholes at constant volatility sqrt(v0). If the Fourier code does
    not reproduce the Black-Scholes number to many decimals, it is wrong and
    nothing else in this script's verdict means anything.
    """
    print("Self-test: xi -> 0 must collapse Heston onto Black-Scholes")
    flat = dict(cfg)
    flat["theta"] = cfg["v0"]        # no drift in v either
    flat["rho"] = 0.0
    reference = black_scholes_call(cfg["spot"], cfg["strike"], cfg["maturity"],
                                   cfg["rate"], cfg["div_yield"],
                                   math.sqrt(cfg["v0"]))
    print(f"  Black-Scholes, sigma = sqrt(v0) = {math.sqrt(cfg['v0']):.4f}"
          f" : {reference:.6f}")
    limit_error = 0.0
    for xi in (1e-1, 1e-2, 1e-3, 1e-4):
        flat["xi"] = xi
        got = semi_analytic_call(flat)
        limit_error = got - reference     # the last (smallest xi) is the limit
        print(f"  semi-analytic, xi = {xi:<7g} : {got:.6f}   "
              f"error {limit_error:+.6f}")
    verdict = "PASS" if abs(limit_error) < 1e-3 else "FAIL"
    print(f"  -> {verdict} (limit error {limit_error:+.2e})")
    print()


# --------------------------------------------------------------------------
# Reporting.
# --------------------------------------------------------------------------

def run_once(cfg: dict, num_paths: int, num_steps: int, seed: int,
             bump_fraction: float, want_delta: bool):
    num_pairs = max(1, num_paths // 2)
    growth, diagnostics = simulate_log_growth(cfg, num_pairs, num_steps, seed)
    price, std_err, plain_err = price_from_growth(growth, cfg, cfg["spot"])
    result = {
        "steps": num_steps,
        "paths": 2 * num_pairs,
        "price": price,
        "std_err": std_err,
        "plain_err": plain_err,
        "half_width": Z_95 * std_err,
        "diagnostics": diagnostics,
    }
    if want_delta:
        delta, delta_err = delta_by_bump(growth, cfg, cfg["spot"], bump_fraction)
        result["delta"] = delta
        result["delta_err"] = delta_err
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Monte-Carlo Heston price to cross-check the PDE solver.")
    parser.add_argument("--config", required=True,
                        help="path to the .cfg the C++ solver reads")
    parser.add_argument("--paths", type=int, default=200000,
                        help="total simulated paths (halved into antithetic "
                             "pairs; default 200000)")
    parser.add_argument("--steps", type=int, default=400,
                        help="time steps per path (default 400)")
    parser.add_argument("--seed", type=int, default=12345,
                        help="numpy default_rng seed; fixes the answer exactly")
    parser.add_argument("--pde-price", type=float, default=None,
                        help="price from ./heston, to test against the CI")
    parser.add_argument("--pde-delta", type=float, default=None,
                        help="delta from ./heston, for the same comparison")
    parser.add_argument("--bump", type=float, default=0.01,
                        help="delta bump as a fraction of spot (default 0.01)")
    parser.add_argument("--no-delta", action="store_true",
                        help="skip the bump-and-revalue delta")
    parser.add_argument("--step-ladder", default=None,
                        help="comma-separated step counts, e.g. 100,200,400 -- "
                             "runs each so the discretisation bias in v is "
                             "visible as a trend, not guessed at")
    parser.add_argument("--no-semi-analytic", action="store_true",
                        help="skip the Fourier-inversion reference price")
    parser.add_argument("--self-test", action="store_true",
                        help="validate the semi-analytic price against "
                             "Black-Scholes in the xi -> 0 limit, then exit")
    parser.add_argument("--gate", type=float, default=None, metavar="REL_TOL",
                        help="exit non-zero if |PDE - semi-analytic|/analytic "
                             "exceeds this. The only check in the repo that "
                             "constrains the xi and rho terms; see the note at "
                             "the bottom of main().")
    args = parser.parse_args()

    cfg = load_config(args.config)

    if args.self_test:
        run_self_test(cfg)
        return

    feller = 2.0 * cfg["kappa"] * cfg["theta"]
    xi_sq = cfg["xi"] ** 2

    print("=" * 74)
    print("Heston Monte-Carlo cross-check")
    print("=" * 74)
    print(f"config          : {args.config}")
    print(f"option          : {'call' if cfg['is_call'] else 'put'}  "
          f"K={cfg['strike']:g}  T={cfg['maturity']:g}y")
    print(f"market          : S0={cfg['spot']:g}  r={cfg['rate']:g}  "
          f"q={cfg['div_yield']:g}")
    print(f"heston          : v0={cfg['v0']:g}  kappa={cfg['kappa']:g}  "
          f"theta={cfg['theta']:g}  xi={cfg['xi']:g}  rho={cfg['rho']:g}")
    feller_verdict = "VIOLATED" if feller < xi_sq else "satisfied"
    print(f"Feller 2*k*theta={feller:.4f} vs xi^2={xi_sq:.4f} -> {feller_verdict}"
          "  (v can reach 0; full truncation required)")
    print(f"scheme          : log-Euler for S, full-truncation Euler for v, "
          f"antithetic pairs")
    print(f"seed            : {args.seed}")
    print()

    ladder = ([int(s) for s in args.step_ladder.split(",")]
              if args.step_ladder else [args.steps])

    header = f"{'steps':>7} {'dt':>10} {'MC price':>12} {'+/- 95% CI':>12} " \
             f"{'CI low':>11} {'CI high':>11}"
    print(header)
    print("-" * len(header))

    results = []
    for steps in ladder:
        res = run_once(cfg, args.paths, steps, args.seed, args.bump,
                       not args.no_delta)
        results.append(res)
        low = res["price"] - res["half_width"]
        high = res["price"] + res["half_width"]
        print(f"{steps:>7d} {res['diagnostics']['dt']:>10.2e} "
              f"{res['price']:>12.4f} {res['half_width']:>12.4f} "
              f"{low:>11.4f} {high:>11.4f}")
    print()

    final = results[-1]
    diag = final["diagnostics"]
    print("Noise diagnostics (finest run)")
    print(f"  sample corr(z1, z2) on step 1 : {diag['first_step_corr']:+.4f}  "
          f"(target rho = {cfg['rho']:+.4f})")
    print(f"  states with v < 0 (truncated) : "
          f"{100.0 * diag['zero_hit_fraction']:.3f}%  "
          "-- the Feller violation, live")
    print(f"  mean v at maturity            : "
          f"{diag['final_variance_mean']:.5f}  (theta = {cfg['theta']:g})")
    print(f"  paths                         : {final['paths']:,}")
    print(f"  std error with antithetics    : {final['std_err']:.4f}")
    print(f"  std error without (same paths): {final['plain_err']:.4f}  "
          f"({final['plain_err'] / final['std_err']:.2f}x wider)")
    print()

    reference = None if args.no_semi_analytic else semi_analytic_call(cfg)

    print("Comparison")
    print("-" * 74)
    low = final["price"] - final["half_width"]
    high = final["price"] + final["half_width"]
    print(f"  MC price        : {final['price']:.4f} "
          f"+/- {final['half_width']:.4f}   [{low:.4f}, {high:.4f}]")

    if reference is not None:
        inside = low <= reference <= high
        gap = reference - final["price"]
        sigmas = gap / final["std_err"] if final["std_err"] > 0 else 0.0
        print(f"  semi-analytic   : {reference:.4f}   "
              f"(Fourier inversion, no randomness, no mesh)")
        print(f"    vs MC         : {gap:+.4f}  ({sigmas:+.2f} MC sigma)  "
              f"inside 95% CI: {'YES' if inside else 'NO'}")

    if args.pde_price is not None:
        inside = low <= args.pde_price <= high
        gap = args.pde_price - final["price"]
        sigmas = gap / final["std_err"] if final["std_err"] > 0 else 0.0
        print(f"  PDE price       : {args.pde_price:.4f}")
        print(f"    vs MC         : {gap:+.4f}  ({sigmas:+.2f} MC sigma)  "
              f"inside 95% CI: {'YES' if inside else 'NO'}")
        if reference is not None:
            pde_gap = args.pde_price - reference
            print(f"    vs analytic   : {pde_gap:+.4f}  "
                  f"({100.0 * pde_gap / reference:+.2f}% of price)")

    if "delta" in final:
        print(f"  MC delta        : {final['delta']:.6f} "
              f"+/- {final['delta_err']:.6f}   "
              f"(bump {100 * args.bump:g}% of spot, common random numbers)")
        if args.pde_delta is not None:
            print(f"  PDE delta       : {args.pde_delta:.6f}   "
                  f"difference {args.pde_delta - final['delta']:+.6f}")
    print("-" * 74)

    if len(results) > 1:
        print("Reading the step ladder: a price that still trends as dt shrinks")
        print("is Euler bias in v; a price that has gone flat has converged, and")
        print("any gap left is the OTHER method's error, not this one's.")
    if reference is not None and args.pde_price is not None:
        print("Three-way reading: MC and semi-analytic are independent of each")
        print("other (random vs deterministic). If those two agree and the PDE")
        print("does not, the PDE mesh is the thing that needs refining.")

    # ---- the gate -----------------------------------------------------------
    #
    # WHY THIS EXISTS. The three C++ tests do not constrain the two terms that
    # make this a Heston solver rather than a Black-Scholes solver:
    #
    #   test_bs_collapse sets xi = 0, which zeroes BOTH the vol-of-vol
    #     diffusion (0.5*xi^2*v*V_vv) and the cross term (rho*xi*v*S*V_Sv).
    #     It therefore cannot see either of them.
    #
    #   test_parity checks C - P against S*exp(-qT) - K*exp(-rT). That identity
    #     holds in EVERY model, so it is blind to the model by construction.
    #     Measured, not assumed: holding the grid fixed and sweeping
    #     xi = 0.01 -> 0.90 and rho = 0.00 -> +0.90 moves the call price by
    #     more than $30 while the parity residual stays at +5.377e-06 to four
    #     significant figures. It is a linearity test, not a Heston test.
    #
    #   test_opt_matches only ever compares the solver against itself, so a
    #     systematic error in those terms cancels on both sides.
    #
    # Comparing the PDE price against the Fourier-inversion price at the FULL
    # parameter set is the only check in the repository that touches them, and
    # until now it only printed. --gate makes it fail.
    if args.gate is not None:
        if reference is None:
            print("\nGATE: cannot run with --no-semi-analytic", file=sys.stderr)
            return 1
        if args.pde_price is None:
            print("\nGATE: --gate needs --pde-price", file=sys.stderr)
            return 1
        rel = abs(args.pde_price - reference) / reference
        ok = rel <= args.gate
        print()
        print("=" * 74)
        print(f"GATE  PDE vs semi-analytic at the full parameter set "
              f"(xi={cfg['xi']:g}, rho={cfg['rho']:+g})")
        print(f"      PDE {args.pde_price:.7f}   analytic {reference:.7f}")
        print(f"      relative gap {rel:.3e}   tolerance {args.gate:.3e}   "
              f"{'PASS' if ok else 'FAIL'}")
        print("=" * 74)
        if not ok:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
