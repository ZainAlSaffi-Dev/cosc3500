// The second validation gate checks put-call parity, which says that
//   C - P = S*exp(-q*T) - K*exp(-r*T)
// This holds in any model at all, so no Heston formula is needed to test it.
// It runs with the full parameter set, which exercises the cross term and
// both variance boundaries that the Black-Scholes collapse test never
// reaches, because that one has to switch the vol of vol off.
//
// The tolerance can be this tight because the parity portfolio pays S - K,
// which is linear in S. Central differences are exact on linear functions and
// every boundary formula handles this one exactly, so the only error left is
// the small gap between compounding step by step and the exact exponential.
// A real bug would show up as dollars rather than fractions of a cent.

#include <cmath>
#include <cstdio>

#include "params.h"
#include "solver.h"

// A node-aligned and stable test grid, the same shape test_bs_collapse.cpp
// uses. The spot has to land exactly on a stock node, because otherwise C - P
// would be compared against a slightly different S. With s_max at 21000 and
// 421 nodes the spacing is 50, so a spot of 5200 is node 104, and with v_max
// at 1 and 51 nodes the spacing is 0.02, so v0 of 0.04 is node 2. The
// timestep count keeps dt at roughly 0.79 of the stability bound.
static Config aligned_test_config() {
    Config cfg;  // params.h defaults are the reference-run market/model
    cfg.grid.num_stock_nodes = 421;
    cfg.grid.num_variance_nodes = 51;
    cfg.grid.num_timesteps = 56000;
    return cfg;
}

int main() {
    Config cfg = aligned_test_config();

    // The two legs share a grid and parameters and differ only in payoff.
    BaselineSolver solver;
    cfg.option.is_call = true;
    const SolveResult call_result = solver.solve(cfg);
    cfg.option.is_call = false;
    const SolveResult put_result = solver.solve(cfg);

    // The comparison means nothing if the runs were unstable, so check.
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

    // Measured on this grid the gap is 5.4e-6 dollars, which is what the
    // compounding argument above predicts. The tolerance is twenty times that
    // so there is room for floating-point differences from rangpur's
    // compiler, and it still sits far below anything a real bug would cause.
    const double abs_tolerance = 1e-4;

    std::printf("parity: call %10.4f  put %10.4f  C-P %10.6f  "
                "S*e^-qT - K*e^-rT %10.6f  gap %.3e  (tol %.1e)\n",
                call_result.price, put_result.price, parity_lhs, parity_rhs,
                abs_gap, abs_tolerance);
    const bool ok = abs_gap < abs_tolerance;
    std::printf("test_parity: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
