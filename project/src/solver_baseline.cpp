// BaselineSolver is the simple version, written to be correct rather than
// fast. It recomputes every stencil weight per cell and hoists nothing. This
// is the "before" number for the optimisation section and the answer every
// other solver has to match, so nothing in here should be optimised.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "io.h"
#include "solver.h"

SolveResult BaselineSolver::solve(const Config& cfg) {
    SolveResult result;

    Grid g(cfg.grid, cfg.option);
    g.init_payoff(cfg.option);

    // Copying the parameters into short locals lets the arithmetic below look
    // like the PDE as it is written on paper.
    const double strike = cfg.option.strike;
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const int num_timesteps = cfg.grid.num_timesteps;
    const double rho = cfg.heston.rho;

    const double dt = cfg.option.maturity_years / num_timesteps;

    // An explicit scheme is only stable below a certain timestep. The worst
    // cell is the far corner where the stock price and the variance are both
    // at their largest, because that is where the diffusion coefficients are
    // biggest, so the estimate comes from there. PLAN §1 works this through.
    const double s_max = g.stock_price(g.num_stock_nodes() - 1);
    const double v_max = g.variance(g.num_variance_nodes() - 1);
    result.dt_stable_estimate = 1.0 /
        (v_max * s_max * s_max / (g.stock_spacing() * g.stock_spacing()) +
         xi * xi * v_max / (g.variance_spacing() * g.variance_spacing()) +
         rate);

    // The Feller condition says that when 2*kappa*theta is at least xi
    // squared, the variance can never reach zero. It is reported here for the
    // write-up and does not change anything the solver does.
    const bool feller_ok = 2.0 * kappa * theta >= xi * xi;
    std::fprintf(stderr,
                 "[baseline] feller %s (2*kappa*theta=%.4f, xi^2=%.4f)  "
                 "dt=%.3e  dt_stable~%.3e\n",
                 feller_ok ? "OK" : "VIOLATED", 2.0 * kappa * theta, xi * xi,
                 dt, result.dt_stable_estimate);

    // These two references stay valid across swap_buffers(), because the swap
    // exchanges the contents of the vectors rather than the references
    // themselves. So cur is always the sheet being read and next is always the
    // sheet being written.
    const std::vector<double>& cur = g.current();
    std::vector<double>& next = g.next();

    // n counts steps taken backwards from expiry towards today.
    for (int n = 1; n <= num_timesteps; ++n) {
        const std::chrono::steady_clock::time_point step_start =
            std::chrono::steady_clock::now();

        for (int var_j = 1; var_j < g.num_variance_nodes() - 1; ++var_j) {
            for (int stock_i = 1; stock_i < g.num_stock_nodes() - 1;
                 ++stock_i) {
                    double S = g.stock_price(stock_i);
                    double V = cur[g.index(stock_i, var_j)];
                    double v = g.variance(var_j);
                    double east = cur[g.index(stock_i + 1, var_j)];
                    double west = cur[g.index(stock_i - 1, var_j)];
                    double north = cur[g.index(stock_i, var_j + 1)];
                    double south = cur[g.index(stock_i, var_j - 1)];
                    double ne = cur[g.index(stock_i + 1, var_j + 1)];
                    double nw = cur[g.index(stock_i - 1, var_j + 1)];
                    double se = cur[g.index(stock_i + 1, var_j - 1)];
                    double sw = cur[g.index(stock_i - 1, var_j - 1)];
                    // finite differences declarations (what we use to calc the pde)
                    double V_S = (east - west) / (2.0 * g.stock_spacing());
                    double V_SS = (east - 2.0 * V + west) / (g.stock_spacing() * g.stock_spacing());
                    double V_v = (north - south) / (2.0 * g.variance_spacing());
                    double V_vv = (north - 2.0 * V + south) / (g.variance_spacing() * g.variance_spacing());
                    double V_Sv = (ne - nw - se + sw) / (4.0 * g.stock_spacing() * g.variance_spacing());
                    // heston pde update, its the model the simulation is based on
                    double pde = V + dt * (0.5 * v * S * S * V_SS + rho * xi * v * S * V_Sv + 0.5 * xi * xi * v * V_vv
                        + (rate - div_yield) * S * V_S + kappa * (theta - v) * V_v - rate * V);
                    next[g.index(stock_i, var_j)] = pde;
            }
        }

        // Every diffusion term carries a factor of v, so along the v=0 row
        // they all vanish and a transport equation is left behind. Solving
        // that is better than imposing a value by hand. The derivative in v
        // has to be a forward one-sided difference rather than a central one,
        // because the drift kappa*theta is positive and therefore points into
        // the grid, and a central difference would need a row at negative
        // variance that does not exist. PLAN §1c covers this.
        for (int stock_i = 1; stock_i < g.num_stock_nodes() - 1; ++stock_i) {
            double S = g.stock_price(stock_i);
            double V = cur[g.index(stock_i, 0)];
            double V_S = (cur[g.index(stock_i + 1, 0)] - cur[g.index(stock_i - 1, 0)]) / (2.0 * g.stock_spacing());
            double V_v = (cur[g.index(stock_i, 1)] - cur[g.index(stock_i, 0)]) / g.variance_spacing();
            double pde = V + dt * ((rate - div_yield) * S * V_S + kappa * theta * V_v - rate * V);
            next[g.index(stock_i, 0)] = pde;
        }

        // How far back from expiry this step has reached, which is what the
        // boundary values below are discounted over.
        const double tau = n * dt;

        // A stock that reaches zero stays there, so at S=0 the call is
        // worthless and the put is certain to pay the strike. At the top of
        // the stock axis the option is so deep in the money that the call
        // behaves like a forward and the put is worth nothing. Both edges are
        // therefore known values rather than something to solve for.
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

        // Vega has flattened out by the top of the variance axis, so that row
        // is given zero slope in v by copying the row underneath it. This runs
        // last because the row underneath has to be finished first.
        const int top_row = g.num_variance_nodes() - 1;
        for (int stock_i = 0; stock_i < g.num_stock_nodes(); ++stock_i) {
            next[g.index(stock_i, top_row)] =
                next[g.index(stock_i, top_row - 1)];
        }

        g.swap_buffers();
        result.seconds += std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - step_start)
                              .count();

        // Snapshots are written after the timer stops, so that file I/O never
        // lands inside the measured time.
        if (cfg.dump_every > 0 && n % cfg.dump_every == 0) {
            dump_snapshot(cfg.dump_dir, n, g);
        }
    }

    extract_result(g, cfg, result);
    // Throughput rather than raw seconds, so that runs on different grid sizes
    // can be compared against each other.
    result.cell_updates_per_sec = static_cast<double>(g.num_stock_nodes()) *
                                  g.num_variance_nodes() * num_timesteps /
                                  result.seconds;
    return result;
}
