// Validation gate 2 (PLAN P3): put-call parity
//   C - P = S*exp(-q*T) - K*exp(-r*T)
// holds in ANY model — no Heston formula needed. Runs the FULL parameter set
// (xi != 0), exercising the cross-term and both variance boundaries that the
// BS-collapse test (xi = 0) cannot reach.
//
// Why the tolerance can be dollar-tight: the parity portfolio's payoff S - K
// is linear in S, central differences are exact on linear functions, and
// every boundary formula treats it exactly — the only survivor is the
// compounding gap (1 - r*dt)^nt vs e^(-rT), of order 1e-5 dollars here
// (STUDY_GUIDE §10 P3 notes). A genuine bug shows up dollars-sized.

#include <cmath>
#include <cstdio>

#include "params.h"
#include "solver.h"

// Node-aligned stable test grid (same shape in test_bs_collapse.cpp):
// spot must land exactly on a stock node because the readout takes the
// nearest cell — otherwise C - P compares against the wrong S.
//   s_max = 4*5250 = 21000, ns = 421 -> spacing 50: spot 5200 = node 104.
//   v_max = 1, nv = 51 -> spacing 0.02: v0 = 0.04 = node 2.
// nt = 56000 keeps dt ~0.79x the explicit stability bound (~5.66e-6).
static Config aligned_test_config() {
    Config cfg;  // params.h defaults are the reference-run market/model
    cfg.grid.num_stock_nodes = 421;
    cfg.grid.num_variance_nodes = 51;
    cfg.grid.num_timesteps = 56000;
    return cfg;
}

int main() {
    Config cfg = aligned_test_config();

    // Same grid, same parameters, only the payoff differs between the legs.
    BaselineSolver solver;
    cfg.option.is_call = true;
    const SolveResult call_result = solver.solve(cfg);
    cfg.option.is_call = false;
    const SolveResult put_result = solver.solve(cfg);

    // Guard: the comparison is meaningless if the runs were unstable.
    const double dt = cfg.option.maturity_years / cfg.grid.num_timesteps;
    if (dt > call_result.dt_stable_estimate) {
        std::printf("test_parity: UNSTABLE dt %.3e > bound %.3e\n", dt,
                    call_result.dt_stable_estimate);
        return 1;
    }

    const double maturity = cfg.option.maturity_years;
    const double parity_rhs =
        cfg.market.spot * std::exp(-cfg.market.div_yield * maturity) -
        cfg.option.strike * std::exp(-cfg.market.rate * maturity);
    const double parity_lhs = call_result.price - put_result.price;
    const double abs_gap = std::fabs(parity_lhs - parity_rhs);

    // Measured on this grid (2026-08-09): gap 5.4e-6 dollars, as the
    // compounding analysis above predicts. Tolerance 20x that — headroom for
    // compiler FP differences on rangpur, still ~5 orders below bug-sized.
    const double abs_tolerance = 1e-4;

    std::printf("parity: call %10.4f  put %10.4f  C-P %10.6f  "
                "S*e^-qT - K*e^-rT %10.6f  gap %.3e  (tol %.1e)\n",
                call_result.price, put_result.price, parity_lhs, parity_rhs,
                abs_gap, abs_tolerance);
    const bool ok = abs_gap < abs_tolerance;
    std::printf("test_parity: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
