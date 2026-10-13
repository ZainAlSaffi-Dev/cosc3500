// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <cmath>
#include <cstdio>

#include "black_scholes.h"
#include "params.h"
#include "solver.h"

// A grid where today's stock price, the strike and v0 all land exactly on nodes.
// Stock spacing 50, variance spacing 0.02, stable timestep.
static Config aligned_test_config() {
    Config cfg;
    cfg.grid.num_stock_nodes = 421;
    cfg.grid.num_variance_nodes = 51;
    cfg.grid.num_timesteps = 56000;
    return cfg;
}

// Switches off the random variance and checks the solver matches Black-Scholes.
// With xi = 0 and v0 = theta the variance never moves.
int main() {
    Config cfg = aligned_test_config();

    // Setting xi to 0 and v0 to theta turns Heston into Black-Scholes.
    cfg.heston.xi = 0.0;
    cfg.heston.v0 = cfg.heston.theta;

    // Measured error is about 1.1e-3. The tolerance is about twice that.
    double rel_tolerance = 2.5e-3;
    double sigma = std::sqrt(cfg.heston.theta);

    // The solver lives on the stack and is freed at the end of main.
    BaselineSolver solver;
    bool all_ok = true;

    for (int leg = 0; leg < 2; ++leg) {
        // First leg is the call, second leg is the put.
        cfg.option.is_call = (leg == 0);
        SolveResult result = solver.solve(cfg);

        // Fail on an unstable run.
        double dt = cfg.option.maturity_years / cfg.grid.num_timesteps;
        if (dt > result.dt_stable_estimate) {
            std::printf("test_bs_collapse: UNSTABLE dt %.3e > bound %.3e\n", dt, result.dt_stable_estimate);
            return 1;
        }

        double closed_form = bs_price(cfg.option.is_call, cfg.market.spot, cfg.option.strike, cfg.market.rate, cfg.market.div_yield, sigma, cfg.option.maturity_years);
        double rel_err = std::fabs(result.price - closed_form) / closed_form;

        const char* leg_name = "put";
        if (cfg.option.is_call) {
            leg_name = "call";
        }
        std::printf("bs_collapse %-4s: solver %10.4f  closed-form %10.4f  rel_err %.3e  (tol %.1e)\n",
                    leg_name, result.price, closed_form, rel_err, rel_tolerance);
        if (rel_err > rel_tolerance) {
            all_ok = false;
        }
    }

    if (all_ok) {
        std::printf("test_bs_collapse: PASS\n");
        return 0;
    } else {
        std::printf("test_bs_collapse: FAIL\n");
        return 1;
    }
}
