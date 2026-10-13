// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <cmath>
#include <cstdio>

#include "params.h"
#include "solver.h"

// A grid where today's stock price and v0 land exactly on nodes.
// Same shape as test_bs_collapse, with a stable timestep.
static Config aligned_test_config() {
    Config cfg;
    cfg.grid.num_stock_nodes = 421;
    cfg.grid.num_variance_nodes = 51;
    cfg.grid.num_timesteps = 56000;
    return cfg;
}

// Checks put-call parity, C - P = S e^(-qT) - K e^(-rT), on the full model.
// Parity holds in any model, so a gap means a bug.
int main() {
    Config cfg = aligned_test_config();

    // Same grid and parameters for both. Only the payoff differs.
    BaselineSolver solver;
    cfg.option.is_call = true;
    SolveResult call_result = solver.solve(cfg);
    cfg.option.is_call = false;
    SolveResult put_result = solver.solve(cfg);

    // Fail on an unstable run.
    double dt = cfg.option.maturity_years / cfg.grid.num_timesteps;
    if (dt > call_result.dt_stable_estimate) {
        std::printf("test_parity: UNSTABLE dt %.3e > bound %.3e\n", dt, call_result.dt_stable_estimate);
        return 1;
    }

    // The two sides of the parity equation.
    double maturity = cfg.option.maturity_years;
    double parity_rhs = cfg.market.spot * std::exp(-cfg.market.div_yield * maturity) - cfg.option.strike * std::exp(-cfg.market.rate * maturity);
    double parity_lhs = call_result.price - put_result.price;
    double abs_gap = std::fabs(parity_lhs - parity_rhs);

    // Measured gap is 5.4e-6 dollars. The tolerance is about twenty times that.
    double abs_tolerance = 1e-4;

    std::printf("parity: call %10.4f  put %10.4f  C-P %10.6f  S*e^-qT - K*e^-rT %10.6f  gap %.3e  (tol %.1e)\n",
                call_result.price, put_result.price, parity_lhs, parity_rhs, abs_gap, abs_tolerance);

    if (abs_gap < abs_tolerance) {
        std::printf("test_parity: PASS\n");
        return 0;
    } else {
        std::printf("test_parity: FAIL\n");
        return 1;
    }
}
