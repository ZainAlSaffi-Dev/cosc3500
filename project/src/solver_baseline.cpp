// BaselineSolver — correct-first, deliberately straightforward:
// stencil weights recomputed per cell, no hoisting, no cleverness.
// This is the "before" of the optimisation story and the reference
// every other solver must match. Resist optimising this file.

#include <algorithm>  // std::clamp
#include <chrono>     // steady_clock timing
#include <cmath>      // std::exp
#include <cstdio>     // fprintf (stderr status line; stdout stays CSV-clean)
#include <stdexcept>
#include <vector>


#include "io.h"
#include "solver.h"

SolveResult BaselineSolver::solve(const Config& cfg) {
    SolveResult result;

    Grid g(cfg.grid, cfg.option);
    g.init_payoff(cfg.option);

    // Short local names so the numerics below read like the PDE.
    const double strike = cfg.option.strike;
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const int num_timesteps = cfg.grid.num_timesteps;
    const double rho = cfg.heston.rho;

    const double dt = cfg.option.maturity_years / num_timesteps;

    // Explicit-scheme speed limit (PLAN §1): the binding cell is the far
    // (s_max, v_max) corner, where the diffusion coefficients are largest.
    const double s_max = g.stock_price(g.num_stock_nodes() - 1);
    const double v_max = g.variance(g.num_variance_nodes() - 1);
    result.dt_stable_estimate = 1.0 /
        (v_max * s_max * s_max / (g.stock_spacing() * g.stock_spacing()) +
         xi * xi * v_max / (g.variance_spacing() * g.variance_spacing()) +
         rate);

    // Feller check (PLAN §1c): 2*kappa*theta >= xi^2 decides whether the
    // variance process can touch zero. Report-only; the scheme stays valid.
    const bool feller_ok = 2.0 * kappa * theta >= xi * xi;
    std::fprintf(stderr,
                 "[baseline] feller %s (2*kappa*theta=%.4f, xi^2=%.4f)  "
                 "dt=%.3e  dt_stable~%.3e\n",
                 feller_ok ? "OK" : "VIOLATED", 2.0 * kappa * theta, xi * xi,
                 dt, result.dt_stable_estimate);

    // Bind the two sheets once. References stay valid across swap_buffers()
    // because swap exchanges the vectors' contents, not these bindings —
    // cur is always the read sheet, next always the write sheet.
    const std::vector<double>& cur = g.current();
    std::vector<double>& next = g.next();
    // Time loop: n counts steps taken backwards from expiry toward today.
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

        // Every diffusion term carries a factor v and vanishes on this row;
        // solve the degenerate transport equation instead of imposing values:
        // V_t + (rate - div_yield)*S*V_S + kappa*theta*V_v - rate*V = 0
        // stock_i = 1 to num_stock_nodes-2 (the row's two endpoints belong to the S boundaries below)
        // V_S: central difference along row var_j = 0 of cur
        // V_v = FORWARD one-sided, (row 1 - row 0)/variance_spacing —
        // upwind matches the inward drift kappa*theta > 0; central would
        // need a meaningless v<0 ghost row (PLAN §1c forbids it)
        // next = V + dt*((rate - div_yield)*S*V_S + kappa*theta*V_v - rate*V)
        for (int stock_i = 1; stock_i < g.num_stock_nodes() - 1; ++stock_i) {
            double S = g.stock_price(stock_i);
            double V = cur[g.index(stock_i, 0)];
            double V_S = (cur[g.index(stock_i + 1, 0)] - cur[g.index(stock_i - 1, 0)]) / (2.0 * g.stock_spacing());
            double V_v = (cur[g.index(stock_i, 1)] - cur[g.index(stock_i, 0)]) / g.variance_spacing();
            double pde = V + dt * ((rate - div_yield) * S * V_S + kappa * theta * V_v - rate * V);
            next[g.index(stock_i, 0)] = pde;
        }

        // Discounting horizon: how far back from expiry we are after this
        // step.
        const double tau = n * dt;

        // S=0 boundary (Dirichlet): the stock is absorbed at zero — a call
        // is worthless, a put pays the discounted strike for certain.
        const double s_zero_value =
            cfg.option.is_call ? 0.0 : strike * std::exp(-rate * tau);
        // S=s_max boundary (Dirichlet): so deep in the money the option is
        // effectively a forward: call = S*e^(-q*tau) - K*e^(-r*tau); put = 0.
        const double s_max_value =
            cfg.option.is_call
                ? s_max * std::exp(-div_yield * tau) -
                      strike * std::exp(-rate * tau)
                : 0.0;
        for (int var_j = 0; var_j < g.num_variance_nodes(); ++var_j) {
            next[g.index(0, var_j)] = s_zero_value;
            next[g.index(g.num_stock_nodes() - 1, var_j)] = s_max_value;
        }

        // v=v_max boundary (Neumann dV/dv = 0): vega has saturated up here;
        // copy the fully-updated row below so the top row has zero v-slope.
        // Runs last so that row is complete (interior + S boundaries).
        const int top_row = g.num_variance_nodes() - 1;
        for (int stock_i = 0; stock_i < g.num_stock_nodes(); ++stock_i) {
            next[g.index(stock_i, top_row)] =
                next[g.index(stock_i, top_row - 1)];
        }

        g.swap_buffers();
        result.seconds += std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - step_start)
                              .count();

        // Snapshots outside the timed region: I/O must not pollute timings.
        if (cfg.dump_every > 0 && n % cfg.dump_every == 0) {
            dump_snapshot(cfg.dump_dir, n, g);
        }
    }

    extract_result(g, cfg, result);
    // Throughput metric — comparable across grid sizes (PLAN §1).
    result.cell_updates_per_sec = static_cast<double>(g.num_stock_nodes()) *
                                  g.num_variance_nodes() * num_timesteps /
                                  result.seconds;
    return result;
}

void Solver::extract_result(const Grid& g, const Config& cfg,
                            SolveResult& out) const {
    const std::vector<double>& V = g.current();
    // Price today = the finished sheet's value at the cell nearest (spot, v0).
    int spot_i = g.nearest_stock_index(cfg.market.spot);
    int v0_j = g.nearest_variance_index(cfg.heston.v0);
    // Greeks need a neighbour on each side; clamp one node in from the edges.
    spot_i = std::clamp(spot_i, 1, g.num_stock_nodes() - 2);
    v0_j = std::clamp(v0_j, 1, g.num_variance_nodes() - 2);
    out.price = V[g.index(spot_i, v0_j)];

    const double east = V[g.index(spot_i + 1, v0_j)];
    const double west = V[g.index(spot_i - 1, v0_j)];
    const double north = V[g.index(spot_i, v0_j + 1)];
    const double south = V[g.index(spot_i, v0_j - 1)];
    // delta and gamma: central differences in the stock direction — free
    // off the finished sheet.
    out.delta = (east - west) / (2.0 * g.stock_spacing());
    out.gamma = (east - 2.0 * out.price + west) /
                (g.stock_spacing() * g.stock_spacing());
    // vega here is per unit VARIANCE; the market's per-volatility vega is
    // this times 2*sqrt(v0) (chain rule — noted in the write-up).
    out.vega = (north - south) / (2.0 * g.variance_spacing());
}

std::unique_ptr<Solver> make_solver(const std::string& name) {
    // std::make_unique wraps `new` so ownership is never loose — the
    // returned unique_ptr deletes the solver when it goes out of scope.
    if (name == "baseline") return std::make_unique<BaselineSolver>();
    if (name == "opt") return std::make_unique<OptSolver>();
    throw std::runtime_error("unknown solver: " + name);
}
