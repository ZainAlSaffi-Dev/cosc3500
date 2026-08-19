// The first validation gate. Setting the vol of vol to zero leaves the
// variance following its deterministic mean-reversion path, and starting it
// at the long-run mean pins it there for good. Heston then collapses to
// Black-Scholes with a volatility of sqrt(theta). The PDE solver never sees
// the closed form, so agreement between the two checks the payoff, the
// stencil, the boundaries and the readout all at once.

#include <cmath>
#include <cstdio>

#include "black_scholes.h"
#include "params.h"
#include "solver.h"

// A node-aligned and stable test grid, the same shape test_parity.cpp uses.
// Both the spot and v0 have to land exactly on nodes, because being half a
// cell off would swamp the discretisation error this test is trying to
// measure. With s_max at 21000 and 421 nodes the spacing is 50, so a spot of
// 5200 is node 104 and the strike of 5250 is node 105, and with v_max at 1
// and 51 nodes the spacing is 0.02, so v0 is node 2.
static Config aligned_test_config() {
    Config cfg;  // params.h defaults are the reference-run market/model
    cfg.grid.num_stock_nodes = 421;
    cfg.grid.num_variance_nodes = 51;
    cfg.grid.num_timesteps = 56000;
    return cfg;
}

int main() {
    Config cfg = aligned_test_config();
    // These two lines are the collapse itself.
    cfg.heston.xi = 0.0;
    cfg.heston.v0 = cfg.heston.theta;

    // Measured on this grid the relative error is 1.161e-3 for the call and
    // 1.116e-3 for the put, which is discretisation error at a spacing of 50
    // and nothing more. The tolerance is about twice that, so a regression
    // trips it but ordinary grid noise does not.
    const double rel_tolerance = 2.5e-3;

    const double sigma = std::sqrt(cfg.heston.theta);
    // The solver lives on the stack. The variable is the object itself, not
    // a handle to one, and it's freed at scope exit. No pointer needed here
    // because nothing is being dispatched polymorphically.
    BaselineSolver solver;
    bool all_ok = true;
    for (int leg = 0; leg < 2; ++leg) {
        cfg.option.is_call = (leg == 0);
        const SolveResult result = solver.solve(cfg);

        // The comparison means nothing if the run was unstable, so check.
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
