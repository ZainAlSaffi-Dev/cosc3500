// OptSolver produces the "after" half of the optimisation comparison from
// PLAN §4b. It follows BaselineSolver::solve almost line for line, and the
// only real difference is the update kernel, which gets chosen once before
// the time loop starts. Choosing it per cell with something like
// "if (level >= k)" would put that branch inside the very timings the ladder
// is supposed to be measuring.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "io.h"
#include "solver.h"
#include "solver_opt_kernels.h"

SolveResult OptSolver::solve(const Config& cfg) {
    // A level with no kernel written yet throws here, before any work starts,
    // because §1b keeps exceptions at the edges and out of the numerics.
    const StepKernel step = kernel_for_level(cfg.opt_level);
    if (step == nullptr)
        throw std::runtime_error(
            "opt-level " + std::to_string(cfg.opt_level) +
            ": kernel not written yet (P7 ladder in progress)");

    SolveResult result;

    Grid g(cfg.grid, cfg.option);
    g.init_payoff(cfg.option);

    const double strike = cfg.option.strike;
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const int num_timesteps = cfg.grid.num_timesteps;

    const double dt = cfg.option.maturity_years / num_timesteps;

    // The same stability estimate and Feller report the baseline prints,
    // tagged with the ladder level so the benchmark logs record which kernel
    // produced them.
    const double s_max = g.stock_price(g.num_stock_nodes() - 1);
    const double v_max = g.variance(g.num_variance_nodes() - 1);
    result.dt_stable_estimate = 1.0 /
        (v_max * s_max * s_max / (g.stock_spacing() * g.stock_spacing()) +
         xi * xi * v_max / (g.variance_spacing() * g.variance_spacing()) +
         rate);
    const bool feller_ok = 2.0 * kappa * theta >= xi * xi;
    std::fprintf(stderr,
                 "[opt-L%d] feller %s (2*kappa*theta=%.4f, xi^2=%.4f)  "
                 "dt=%.3e  dt_stable~%.3e\n",
                 cfg.opt_level, feller_ok ? "OK" : "VIOLATED",
                 2.0 * kappa * theta, xi * xi, dt, result.dt_stable_estimate);

    // These references survive swap_buffers() for the reason given in
    // solver_baseline.cpp, so cur is always the sheet being read.
    const std::vector<double>& cur = g.current();
    std::vector<double>& next = g.next();
    for (int n = 1; n <= num_timesteps; ++n) {
        const std::chrono::steady_clock::time_point step_start =
            std::chrono::steady_clock::now();

        // This call is the part being measured. It updates the interior and
        // the v=0 row using whichever techniques this level applies.
        step(g, cur, next, cfg, dt);

        // The boundaries are handled identically for every level and for the
        // baseline, and solver_baseline.cpp explains what each one means.
        const double tau = n * dt;
        const double s_zero_value =
            cfg.option.is_call ? 0.0 : strike * std::exp(-rate * tau);
        const double s_max_value =
            cfg.option.is_call
                ? s_max * std::exp(-div_yield * tau) -
                      strike * std::exp(-rate * tau)
                : 0.0;
        for (int var_j = 0; var_j < g.num_variance_nodes(); ++var_j) {
            next[g.index(0, var_j)] = s_zero_value;
            next[g.index(g.num_stock_nodes() - 1, var_j)] = s_max_value;
        }
        const int top_row = g.num_variance_nodes() - 1;
        for (int stock_i = 0; stock_i < g.num_stock_nodes(); ++stock_i) {
            next[g.index(stock_i, top_row)] =
                next[g.index(stock_i, top_row - 1)];
        }

        g.swap_buffers();
        result.seconds += std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - step_start)
                              .count();

        // Snapshots are written after the timer stops, as in the baseline.
        if (cfg.dump_every > 0 && n % cfg.dump_every == 0) {
            dump_snapshot(cfg.dump_dir, n, g);
        }
    }

    extract_result(g, cfg, result);
    result.cell_updates_per_sec = static_cast<double>(g.num_stock_nodes()) *
                                  g.num_variance_nodes() * num_timesteps /
                                  result.seconds;
    return result;
}
