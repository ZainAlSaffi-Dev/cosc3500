// The §4b ablation ladder, one kernel per level. Duplication between
// kernels is ACCEPTED by design — the ladder is the ablation instrument,
// and diffing two adjacent kernels must show exactly one technique.

#include "solver_opt_kernels.h"

// Level 0 — no techniques. Line-for-line the BaselineSolver arithmetic
// (interior stencil + v=0 degenerate row): the sanity anchor proving the
// ladder harness itself adds no cost and changes no answer.
static void step_level0(const Grid& g, const std::vector<double>& cur,
                        std::vector<double>& next, const Config& cfg,
                        double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    for (int var_j = 1; var_j < g.num_variance_nodes() - 1; ++var_j) {
        for (int stock_i = 1; stock_i < g.num_stock_nodes() - 1; ++stock_i) {
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
            double V_S = (east - west) / (2.0 * g.stock_spacing());
            double V_SS = (east - 2.0 * V + west) /
                          (g.stock_spacing() * g.stock_spacing());
            double V_v = (north - south) / (2.0 * g.variance_spacing());
            double V_vv = (north - 2.0 * V + south) /
                          (g.variance_spacing() * g.variance_spacing());
            double V_Sv = (ne - nw - se + sw) /
                          (4.0 * g.stock_spacing() * g.variance_spacing());
            double pde = V + dt * (0.5 * v * S * S * V_SS +
                                   rho * xi * v * S * V_Sv +
                                   0.5 * xi * xi * v * V_vv +
                                   (rate - div_yield) * S * V_S +
                                   kappa * (theta - v) * V_v - rate * V);
            next[g.index(stock_i, var_j)] = pde;
        }
    }

    // v=0 boundary: degenerate transport row — central in S, forward
    // (upwind) one-sided in v, matching the inward drift kappa*theta > 0.
    for (int stock_i = 1; stock_i < g.num_stock_nodes() - 1; ++stock_i) {
        double S = g.stock_price(stock_i);
        double V = cur[g.index(stock_i, 0)];
        double V_S = (cur[g.index(stock_i + 1, 0)] -
                      cur[g.index(stock_i - 1, 0)]) /
                     (2.0 * g.stock_spacing());
        double V_v = (cur[g.index(stock_i, 1)] - cur[g.index(stock_i, 0)]) /
                     g.variance_spacing();
        double pde = V + dt * ((rate - div_yield) * S * V_S +
                               kappa * theta * V_v - rate * V);
        next[g.index(stock_i, 0)] = pde;
    }
}

// Level 1 — per-row stencil-weight lookup table (PLAN §4 technique 1;
// L03 names: hoisting + lookup tables + common subexpression elimination).
// *** AUTHOR-WRITTEN KERNEL (working agreement §1d) — skeleton hints: ***
//   - Start from a copy of step_level0.
//   - Everything v-dependent is constant along a row: hoist v, 0.5*xi*xi*v,
//     rho*xi*v, kappa*(theta - v) out of the inner loop, computed once per
//     var_j.
//   - Everything S-dependent repeats identically for every row: S and S*S
//     per stock_i can live in std::vector lookup tables built once at the
//     top of the kernel (O(ns) work per step vs O(ns*nv) saved).
//   - The spacings never change: name 1/(2*ds), 1/(ds*ds), ... once.
//   - Keep the ARITHMETIC identical in value — test_opt_matches must stay
//     green at 1e-12 against baseline.
//
// static void step_level1(const Grid& g, const std::vector<double>& cur,
//                         std::vector<double>& next, const Config& cfg,
//                         double dt) { ... }

// Levels 2..6 (strength reduction, traversal order, loop splitting,
// induction variable, unrolling) land here after level 1, one commit each.

StepKernel kernel_for_level(int level) {
    switch (level) {
        case 0: return &step_level0;
        default: return nullptr;  // unwritten rungs — P7 in progress
    }
}
