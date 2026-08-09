// Validation gate 1 (PLAN P3): with xi = 0 the volatility of variance is
// zero, so variance follows its deterministic mean-reversion path; starting
// it AT the long-run mean (v0 = theta) pins it there forever — Heston
// degenerates to Black-Scholes with sigma = sqrt(theta). The PDE solver
// never sees the closed form, so agreement here validates payoff, stencil,
// boundaries and readout end to end.

#include <cmath>
#include <cstdio>

#include "black_scholes.h"
#include "params.h"
#include "solver.h"

// Node-aligned stable test grid (same shape in test_parity.cpp):
// extract_result reads the cell NEAREST (spot, v0) with no interpolation, so
// spot and v0 must land exactly on nodes — a half-cell snap would swamp the
// discretisation error this test measures (STUDY_GUIDE §10 P3 notes).
//   s_max = 4*5250 = 21000, ns = 421 -> spacing 50: spot 5200 = node 104,
//   strike 5250 = node 105.  v_max = 1, nv = 51 -> spacing 0.02: v0 = node 2.
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
    // The collapse: kill vol-of-vol, start variance at its resting point.
    cfg.heston.xi = 0.0;
    cfg.heston.v0 = cfg.heston.theta;

    // Measured on this grid (2026-08-09): rel err 1.161e-3 (call),
    // 1.116e-3 (put) — pure discretisation error at spacing 50.
    // Tolerance ~2x measured: a regression trips it, grid noise does not.
    const double rel_tolerance = 2.5e-3;

    const double sigma = std::sqrt(cfg.heston.theta);
    // Solver lives on the stack — it IS the object here, freed at scope exit;
    // no pointer needed because there is no polymorphic dispatch to do.
    BaselineSolver solver;
    bool all_ok = true;
    for (int leg = 0; leg < 2; ++leg) {
        cfg.option.is_call = (leg == 0);
        const SolveResult result = solver.solve(cfg);

        // Guard: the comparison is meaningless if the run was unstable.
        const double dt = cfg.option.maturity_years / cfg.grid.num_timesteps;
        if (dt > result.dt_stable_estimate) {
            std::printf("test_bs_collapse: UNSTABLE dt %.3e > bound %.3e\n",
                        dt, result.dt_stable_estimate);
            return 1;
        }

        const double closed_form = bs_price(
            cfg.option.is_call, cfg.market.spot, cfg.option.strike,
            cfg.market.rate, cfg.market.div_yield, sigma,
            cfg.option.maturity_years);
        const double rel_err =
            std::fabs(result.price - closed_form) / closed_form;
        std::printf("bs_collapse %-4s: solver %10.4f  closed-form %10.4f  "
                    "rel_err %.3e  (tol %.1e)\n",
                    cfg.option.is_call ? "call" : "put", result.price,
                    closed_form, rel_err, rel_tolerance);
        if (rel_err > rel_tolerance) all_ok = false;
    }
    std::printf("test_bs_collapse: %s\n", all_ok ? "PASS" : "FAIL");
    return all_ok ? 0 : 1;
}
