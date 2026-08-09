// OptSolver — the measured half of the optimisation story (PLAN §4b).
// The structure mirrors BaselineSolver::solve line for line; the ONE
// difference is the update kernel, chosen ONCE before the time loop from
// cfg.opt_level. Never a per-cell "if (level >= k)" — that would poison
// the very timings the ladder exists to measure.

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
    // Dispatch once, at the edge: an unwritten ladder level throws here,
    // before any work (§1b — exceptions only at the edges, never in
    // numerics).
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

    // Same stability estimate and Feller report as baseline (PLAN §1/§1c),
    // tagged with the ladder level so bench logs are self-describing.
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

    // cur/next bindings survive swap_buffers() — swap exchanges the
    // vectors' contents, so cur is always the read sheet.
    const std::vector<double>& cur = g.current();
    std::vector<double>& next = g.next();
    for (int n = 1; n <= num_timesteps; ++n) {
        const std::chrono::steady_clock::time_point step_start =
            std::chrono::steady_clock::now();

        // The ladder under test: interior + v=0 row, per-level technique.
        step(g, cur, next, cfg, dt);

        // Shared boundary plumbing, identical for every level (and to
        // baseline) — see solver_baseline.cpp for the full narration.
        const double tau = n * dt;
        // S=0 boundary (Dirichlet): call worthless, put pays discounted K.
        const double s_zero_value =
            cfg.option.is_call ? 0.0 : strike * std::exp(-rate * tau);
        // S=s_max boundary (Dirichlet): deep-ITM forward value / zero.
        const double s_max_value =
            cfg.option.is_call
                ? s_max * std::exp(-div_yield * tau) -
                      strike * std::exp(-rate * tau)
                : 0.0;
        for (int var_j = 0; var_j < g.num_variance_nodes(); ++var_j) {
            next[g.index(0, var_j)] = s_zero_value;
            next[g.index(g.num_stock_nodes() - 1, var_j)] = s_max_value;
        }
        // v=v_max boundary (Neumann dV/dv = 0): copy the completed row.
        const int top_row = g.num_variance_nodes() - 1;
        for (int stock_i = 0; stock_i < g.num_stock_nodes(); ++stock_i) {
            next[g.index(stock_i, top_row)] =
                next[g.index(stock_i, top_row - 1)];
        }

        g.swap_buffers();
        result.seconds += std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - step_start)
                              .count();

        // Snapshots outside the timed region, as in baseline.
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
