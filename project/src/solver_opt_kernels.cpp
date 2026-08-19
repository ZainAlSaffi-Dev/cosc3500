// The optimisation ladder from PLAN §4b, written as one kernel per level.
//
// These kernels repeat most of each other on purpose. Each one is a complete
// copy of the one below it with a single technique added, so the diff between
// two neighbours shows that technique and nothing else, and any difference in
// their timings can only have come from it. Sharing the setup in a helper
// would be shorter but would make a rung able to change what its neighbour
// measures, which is exactly what an ablation has to rule out. That is also
// why this file is far longer than the ~200 lines §1b asks for.
//
// Levels 0 to 6 are the ladder itself. Levels 7 and 8 are negative controls,
// meaning deliberately worse code that does identical arithmetic. They exist
// to measure what levels 3 and 4 would have been worth had the baseline not
// already been written correctly. They are labelled opt-ctl-order and
// opt-ctl-branch rather than opt-L7 and opt-L8 so that nobody reads a CSV row
// as the ladder carrying on upwards.
//
// Levels 0 and 1 return bit-identical answers to BaselineSolver, because
// level 1 only names subexpressions and never reorders or rewrites them.
// Level 2 introduces reciprocals and so genuinely changes the arithmetic.
// From there up, levels 2 to 8 are bit-identical to each other and agree with
// the baseline to a measured tolerance, which tests/test_opt_matches checks.

#include <vector>

#include "solver_opt_kernels.h"

// Level 0 applies no techniques at all and repeats BaselineSolver's
// arithmetic line for line. It anchors the ladder by showing that the harness
// around these kernels costs nothing and changes no answers.
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

    // The v=0 row is pure transport, so it uses a central difference in S and
    // a forward one-sided difference in v to match the inward drift.
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

// Level 1 adds loop-invariant code motion, a lookup table and common
// subexpression elimination, which are techniques 7, 9 and 13 from L03. Level
// 0 recomputes quantities in every one of a million cells that really only
// change once per row, or never, so this rung works each one out at the level
// where it actually varies.
//
// Hoisting here may only give a subexpression a name, never change the order
// it was multiplied in, because floating-point multiplication is not
// associative and this rung has to stay bit-identical. Divisions stay as
// divisions for the same reason, and replacing them is level 2's job.
static void step_level1(const Grid& g, const std::vector<double>& cur,
                        std::vector<double>& next, const Config& cfg,
                        double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    // These spacing products are the same for every cell of every row of
    // every timestep, and level 0 recomputed all five of them per cell.
    const double two_ds = 2.0 * stock_spacing;
    const double ds_sq = stock_spacing * stock_spacing;
    const double two_dv = 2.0 * variance_spacing;
    const double dv_sq = variance_spacing * variance_spacing;
    const double four_ds_dv = 4.0 * stock_spacing * variance_spacing;
    // The cost of carry and the v=0 drift are constants of the market.
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    // S is the same on every row, so the column gets built once per step.
    // That trades ns multiplications for the ns*nv the loop would otherwise
    // do. The vector frees itself when the kernel returns.
    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        stock[stock_i] = stock_i * stock_spacing;
    }

    for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
        // Every coefficient here depends only on v, so each one is worked out
        // once per row rather than once per cell. On the reference grid that
        // is around two thousand fewer evaluations of each.
        const double v = var_j * variance_spacing;
        const double half_v = 0.5 * v;
        const double rho_xi_v = (rho * xi) * v;
        const double half_xi2_v = ((0.5 * xi) * xi) * v;
        const double mean_rev = kappa * (theta - v);

        for (int stock_i = 1; stock_i < num_stock_nodes - 1; ++stock_i) {
            // Each neighbour is loaded once into a name and then reused
            // across two or three of the derivatives below.
            const double S = stock[stock_i];
            const double V = cur[g.index(stock_i, var_j)];
            const double east = cur[g.index(stock_i + 1, var_j)];
            const double west = cur[g.index(stock_i - 1, var_j)];
            const double north = cur[g.index(stock_i, var_j + 1)];
            const double south = cur[g.index(stock_i, var_j - 1)];
            const double ne = cur[g.index(stock_i + 1, var_j + 1)];
            const double nw = cur[g.index(stock_i - 1, var_j + 1)];
            const double se = cur[g.index(stock_i + 1, var_j - 1)];
            const double sw = cur[g.index(stock_i - 1, var_j - 1)];
            const double V_S = (east - west) / two_ds;
            const double V_SS = (east - 2.0 * V + west) / ds_sq;
            const double V_v = (north - south) / two_dv;
            const double V_vv = (north - 2.0 * V + south) / dv_sq;
            const double V_Sv = (ne - nw - se + sw) / four_ds_dv;
            next[g.index(stock_i, var_j)] =
                V + dt * ((half_v * S) * S * V_SS + (rho_xi_v * S) * V_Sv +
                          half_xi2_v * V_vv + (carry * S) * V_S +
                          mean_rev * V_v - rate * V);
        }
    }

    // The v=0 transport row again, with the same quantities hoisted out.
    for (int stock_i = 1; stock_i < num_stock_nodes - 1; ++stock_i) {
        const double S = stock[stock_i];
        const double V = cur[g.index(stock_i, 0)];
        const double V_S =
            (cur[g.index(stock_i + 1, 0)] - cur[g.index(stock_i - 1, 0)]) /
            two_ds;
        const double V_v =
            (cur[g.index(stock_i, 1)] - cur[g.index(stock_i, 0)]) /
            variance_spacing;
        next[g.index(stock_i, 0)] =
            V + dt * ((carry * S) * V_S + kappa_theta * V_v - rate * V);
    }
}

// Level 2 adds strength reduction, technique 15 from L03, which here means
// getting the divisions out of the hot loop. A double division costs roughly
// 5 to 14 cycles and does not pipeline, while a multiply costs about one, so
// each of the five constants is inverted once per step and multiplied by
// instead.
//
// This is the rung where the last bits start to move, because x/(2*ds) and
// x*(1/(2*ds)) are different doubles once the reciprocal has itself been
// rounded. That is why test_opt_matches stops demanding bit-identical answers
// from here upwards.
static void step_level2(const Grid& g, const std::vector<double>& cur,
                        std::vector<double>& next, const Config& cfg,
                        double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    // Five divisions per cell become five reciprocals per step, and the loop
    // itself is left with plain multiplies.
    const double inv_two_ds = 1.0 / (2.0 * stock_spacing);
    const double inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    const double inv_two_dv = 1.0 / (2.0 * variance_spacing);
    const double inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    const double inv_four_ds_dv =
        1.0 / (4.0 * stock_spacing * variance_spacing);
    const double inv_dv = 1.0 / variance_spacing;
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    // Both S and S squared are tabulated per column, once per step.
    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    std::vector<double> stock_sq(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        const double S = stock_i * stock_spacing;
        stock[stock_i] = S;
        stock_sq[stock_i] = S * S;
    }

    for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
        const double v = var_j * variance_spacing;
        const double half_v = 0.5 * v;
        const double rho_xi_v = (rho * xi) * v;
        const double half_xi2_v = ((0.5 * xi) * xi) * v;
        const double mean_rev = kappa * (theta - v);

        for (int stock_i = 1; stock_i < num_stock_nodes - 1; ++stock_i) {
            const double V = cur[g.index(stock_i, var_j)];
            const double east = cur[g.index(stock_i + 1, var_j)];
            const double west = cur[g.index(stock_i - 1, var_j)];
            const double north = cur[g.index(stock_i, var_j + 1)];
            const double south = cur[g.index(stock_i, var_j - 1)];
            const double ne = cur[g.index(stock_i + 1, var_j + 1)];
            const double nw = cur[g.index(stock_i - 1, var_j + 1)];
            const double se = cur[g.index(stock_i + 1, var_j - 1)];
            const double sw = cur[g.index(stock_i - 1, var_j - 1)];
            const double V_S = (east - west) * inv_two_ds;
            const double V_SS = (east - 2.0 * V + west) * inv_ds_sq;
            const double V_v = (north - south) * inv_two_dv;
            const double V_vv = (north - 2.0 * V + south) * inv_dv_sq;
            const double V_Sv = (ne - nw - se + sw) * inv_four_ds_dv;
            next[g.index(stock_i, var_j)] =
                V + dt * (half_v * stock_sq[stock_i] * V_SS +
                          rho_xi_v * stock[stock_i] * V_Sv +
                          half_xi2_v * V_vv + carry * stock[stock_i] * V_S +
                          mean_rev * V_v - rate * V);
        }
    }

    // The v=0 transport row, with its divisions removed as well.
    for (int stock_i = 1; stock_i < num_stock_nodes - 1; ++stock_i) {
        const double V = cur[g.index(stock_i, 0)];
        const double V_S =
            (cur[g.index(stock_i + 1, 0)] - cur[g.index(stock_i - 1, 0)]) *
            inv_two_ds;
        const double V_v =
            (cur[g.index(stock_i, 1)] - cur[g.index(stock_i, 0)]) * inv_dv;
        next[g.index(stock_i, 0)] =
            V + dt * (carry * stock[stock_i] * V_S + kappa_theta * V_v -
                      rate * V);
    }
}

// Level 3 makes the cache-friendly traversal explicit, technique 22 from L03,
// by walking memory in the order it is actually laid out and addressing the
// three stencil rows from a base worked out once per row.
//
// This rung is expected to measure as nothing, because the baseline was
// already written with the right loop order and so has no lost performance to
// win back. What getting the order wrong actually costs is measured by the
// ctl-order control further down, which is also what the cache animation
// illustrates.
static void step_level3(const Grid& g, const std::vector<double>& cur,
                        std::vector<double>& next, const Config& cfg,
                        double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    const double inv_two_ds = 1.0 / (2.0 * stock_spacing);
    const double inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    const double inv_two_dv = 1.0 / (2.0 * variance_spacing);
    const double inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    const double inv_four_ds_dv =
        1.0 / (4.0 * stock_spacing * variance_spacing);
    const double inv_dv = 1.0 / variance_spacing;
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    std::vector<double> stock_sq(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        const double S = stock_i * stock_spacing;
        stock[stock_i] = S;
        stock_sq[stock_i] = S * S;
    }

    // With variance on the outside and stock on the inside, consecutive
    // iterations touch consecutive doubles, so every 64-byte cache line that
    // arrives has all eight of its doubles used before it is evicted.
    for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
        // The three rows the stencil reads, each addressed once per row.
        const int row = var_j * num_stock_nodes;
        const int row_above = row + num_stock_nodes;
        const int row_below = row - num_stock_nodes;

        const double v = var_j * variance_spacing;
        const double half_v = 0.5 * v;
        const double rho_xi_v = (rho * xi) * v;
        const double half_xi2_v = ((0.5 * xi) * xi) * v;
        const double mean_rev = kappa * (theta - v);

        for (int stock_i = 1; stock_i < num_stock_nodes - 1; ++stock_i) {
            const double V = cur[row + stock_i];
            const double east = cur[row + stock_i + 1];
            const double west = cur[row + stock_i - 1];
            const double north = cur[row_above + stock_i];
            const double south = cur[row_below + stock_i];
            const double ne = cur[row_above + stock_i + 1];
            const double nw = cur[row_above + stock_i - 1];
            const double se = cur[row_below + stock_i + 1];
            const double sw = cur[row_below + stock_i - 1];
            const double V_S = (east - west) * inv_two_ds;
            const double V_SS = (east - 2.0 * V + west) * inv_ds_sq;
            const double V_v = (north - south) * inv_two_dv;
            const double V_vv = (north - 2.0 * V + south) * inv_dv_sq;
            const double V_Sv = (ne - nw - se + sw) * inv_four_ds_dv;
            next[row + stock_i] =
                V + dt * (half_v * stock_sq[stock_i] * V_SS +
                          rho_xi_v * stock[stock_i] * V_Sv +
                          half_xi2_v * V_vv + carry * stock[stock_i] * V_S +
                          mean_rev * V_v - rate * V);
        }
    }

    // The v=0 row, addressed the same way as the interior above.
    const int row_zero = 0;
    const int row_one = num_stock_nodes;
    for (int stock_i = 1; stock_i < num_stock_nodes - 1; ++stock_i) {
        const double V = cur[row_zero + stock_i];
        const double V_S = (cur[row_zero + stock_i + 1] -
                            cur[row_zero + stock_i - 1]) * inv_two_ds;
        const double V_v =
            (cur[row_one + stock_i] - cur[row_zero + stock_i]) * inv_dv;
        next[row_zero + stock_i] =
            V + dt * (carry * stock[stock_i] * V_S + kappa_theta * V_v -
                      rate * V);
    }
}

// Level 4 adds loop splitting and a branch-free interior, techniques 18 and
// 20 from L03. Which equation a cell follows depends on its row rather than
// on the cell, so the v=0 row and the interior get separate loops with no
// conditional between them, and the loop bounds become plain locals.
//
// This is another expected null, because the baseline already peels that row
// out. The fused version that tests every cell is measured by the ctl-branch
// control below.
static void step_level4(const Grid& g, const std::vector<double>& cur,
                        std::vector<double>& next, const Config& cfg,
                        double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    const double inv_two_ds = 1.0 / (2.0 * stock_spacing);
    const double inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    const double inv_two_dv = 1.0 / (2.0 * variance_spacing);
    const double inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    const double inv_four_ds_dv =
        1.0 / (4.0 * stock_spacing * variance_spacing);
    const double inv_dv = 1.0 / variance_spacing;
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    std::vector<double> stock_sq(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        const double S = stock_i * stock_spacing;
        stock[stock_i] = S;
        stock_sq[stock_i] = S * S;
    }

    // The bounds are worked out once, so each loop condition is a plain
    // integer comparison and neither body contains a conditional.
    const int last_stock = num_stock_nodes - 1;
    const int last_var = num_variance_nodes - 1;

    // Every cell in the interior follows the full PDE with no special cases.
    for (int var_j = 1; var_j < last_var; ++var_j) {
        const int row = var_j * num_stock_nodes;
        const int row_above = row + num_stock_nodes;
        const int row_below = row - num_stock_nodes;

        const double v = var_j * variance_spacing;
        const double half_v = 0.5 * v;
        const double rho_xi_v = (rho * xi) * v;
        const double half_xi2_v = ((0.5 * xi) * xi) * v;
        const double mean_rev = kappa * (theta - v);

        for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
            const double V = cur[row + stock_i];
            const double east = cur[row + stock_i + 1];
            const double west = cur[row + stock_i - 1];
            const double north = cur[row_above + stock_i];
            const double south = cur[row_below + stock_i];
            const double ne = cur[row_above + stock_i + 1];
            const double nw = cur[row_above + stock_i - 1];
            const double se = cur[row_below + stock_i + 1];
            const double sw = cur[row_below + stock_i - 1];
            const double V_S = (east - west) * inv_two_ds;
            const double V_SS = (east - 2.0 * V + west) * inv_ds_sq;
            const double V_v = (north - south) * inv_two_dv;
            const double V_vv = (north - 2.0 * V + south) * inv_dv_sq;
            const double V_Sv = (ne - nw - se + sw) * inv_four_ds_dv;
            next[row + stock_i] =
                V + dt * (half_v * stock_sq[stock_i] * V_SS +
                          rho_xi_v * stock[stock_i] * V_Sv +
                          half_xi2_v * V_vv + carry * stock[stock_i] * V_S +
                          mean_rev * V_v - rate * V);
        }
    }

    // The transport row is peeled out into its own loop rather than being
    // picked out by a branch on every cell.
    const int row_one = num_stock_nodes;
    for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
        const double V = cur[stock_i];
        const double V_S = (cur[stock_i + 1] - cur[stock_i - 1]) * inv_two_ds;
        const double V_v = (cur[row_one + stock_i] - cur[stock_i]) * inv_dv;
        next[stock_i] = V + dt * (carry * stock[stock_i] * V_S +
                                  kappa_theta * V_v - rate * V);
    }
}

// Level 5 adds induction variable simplification, technique 14 from L03.
// Level 4 still computes `row + stock_i` for each of the ten accesses a cell
// makes, so this rung holds a raw pointer to the start of each row instead
// and lets the hardware's addressing mode do the work.
//
// Calling `.data()` hands out a raw pointer into the vector's buffer. PLAN
// §1b allows that only here, as locals inside the hot kernel that are never
// stored and never own anything, since the vectors still own every byte.
static void step_level5(const Grid& g, const std::vector<double>& cur,
                        std::vector<double>& next, const Config& cfg,
                        double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    const double inv_two_ds = 1.0 / (2.0 * stock_spacing);
    const double inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    const double inv_two_dv = 1.0 / (2.0 * variance_spacing);
    const double inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    const double inv_four_ds_dv =
        1.0 / (4.0 * stock_spacing * variance_spacing);
    const double inv_dv = 1.0 / variance_spacing;
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    std::vector<double> stock_sq(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        const double S = stock_i * stock_spacing;
        stock[stock_i] = S;
        stock_sq[stock_i] = S * S;
    }

    const int last_stock = num_stock_nodes - 1;
    const int last_var = num_variance_nodes - 1;
    // Base pointers into the two sheets. These are locals, and ownership of
    // the memory never moves anywhere.
    const double* const cur_base = cur.data();
    double* const next_base = next.data();

    for (int var_j = 1; var_j < last_var; ++var_j) {
        // The row bases move on by exactly one row per iteration, which turns
        // the per-cell index arithmetic into plain pointer indexing.
        const double* const row_mid = cur_base + var_j * num_stock_nodes;
        const double* const row_above = row_mid + num_stock_nodes;
        const double* const row_below = row_mid - num_stock_nodes;
        double* const out_row = next_base + var_j * num_stock_nodes;

        const double v = var_j * variance_spacing;
        const double half_v = 0.5 * v;
        const double rho_xi_v = (rho * xi) * v;
        const double half_xi2_v = ((0.5 * xi) * xi) * v;
        const double mean_rev = kappa * (theta - v);

        for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
            const double V = row_mid[stock_i];
            const double east = row_mid[stock_i + 1];
            const double west = row_mid[stock_i - 1];
            const double north = row_above[stock_i];
            const double south = row_below[stock_i];
            const double ne = row_above[stock_i + 1];
            const double nw = row_above[stock_i - 1];
            const double se = row_below[stock_i + 1];
            const double sw = row_below[stock_i - 1];
            const double V_S = (east - west) * inv_two_ds;
            const double V_SS = (east - 2.0 * V + west) * inv_ds_sq;
            const double V_v = (north - south) * inv_two_dv;
            const double V_vv = (north - 2.0 * V + south) * inv_dv_sq;
            const double V_Sv = (ne - nw - se + sw) * inv_four_ds_dv;
            out_row[stock_i] =
                V + dt * (half_v * stock_sq[stock_i] * V_SS +
                          rho_xi_v * stock[stock_i] * V_Sv +
                          half_xi2_v * V_vv + carry * stock[stock_i] * V_S +
                          mean_rev * V_v - rate * V);
        }
    }

    // The v=0 row, with rows 0 and 1 held as pointers.
    const double* const row_zero = cur_base;
    const double* const row_one = cur_base + num_stock_nodes;
    double* const out_zero = next_base;
    for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
        const double V = row_zero[stock_i];
        const double V_S =
            (row_zero[stock_i + 1] - row_zero[stock_i - 1]) * inv_two_ds;
        const double V_v = (row_one[stock_i] - row_zero[stock_i]) * inv_dv;
        out_zero[stock_i] = V + dt * (carry * stock[stock_i] * V_S +
                                      kappa_theta * V_v - rate * V);
    }
}

// Level 6 unrolls the inner loop by four, technique 17 from L03, so the loop
// overhead is spread over four cells and the processor gets four independent
// dependency chains to interleave rather than one. Whether that helps is
// worth measuring rather than assuming, since at -O2 the compiler already
// unrolls simple counted loops.
//
// This is the only rung whose per-cell arithmetic is unchanged from the one
// below, so levels 5 and 6 have to agree bit for bit, which doubles as a
// check that the unrolled body was copied out correctly.
static void step_level6(const Grid& g, const std::vector<double>& cur,
                        std::vector<double>& next, const Config& cfg,
                        double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    const double inv_two_ds = 1.0 / (2.0 * stock_spacing);
    const double inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    const double inv_two_dv = 1.0 / (2.0 * variance_spacing);
    const double inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    const double inv_four_ds_dv =
        1.0 / (4.0 * stock_spacing * variance_spacing);
    const double inv_dv = 1.0 / variance_spacing;
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    std::vector<double> stock_sq(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        const double S = stock_i * stock_spacing;
        stock[stock_i] = S;
        stock_sq[stock_i] = S * S;
    }

    const int last_stock = num_stock_nodes - 1;
    const int last_var = num_variance_nodes - 1;
    const double* const cur_base = cur.data();
    double* const next_base = next.data();

    for (int var_j = 1; var_j < last_var; ++var_j) {
        const double* const row_mid = cur_base + var_j * num_stock_nodes;
        const double* const row_above = row_mid + num_stock_nodes;
        const double* const row_below = row_mid - num_stock_nodes;
        double* const out_row = next_base + var_j * num_stock_nodes;

        const double v = var_j * variance_spacing;
        const double half_v = 0.5 * v;
        const double rho_xi_v = (rho * xi) * v;
        const double half_xi2_v = ((0.5 * xi) * xi) * v;
        const double mean_rev = kappa * (theta - v);

        int stock_i = 1;
        // Four independent cells per loop iteration. Each block is level 5's
        // body copied out with its index and local names numbered, and the
        // only real difference is that four of them now sit between one pair
        // of loop-back branches instead of one.
        for (; stock_i + 3 < last_stock; stock_i += 4) {
            const int i0 = stock_i;
            const int i1 = stock_i + 1;
            const int i2 = stock_i + 2;
            const int i3 = stock_i + 3;

            const double V0 = row_mid[i0];
            const double east0 = row_mid[i0 + 1];
            const double west0 = row_mid[i0 - 1];
            const double north0 = row_above[i0];
            const double south0 = row_below[i0];
            const double ne0 = row_above[i0 + 1];
            const double nw0 = row_above[i0 - 1];
            const double se0 = row_below[i0 + 1];
            const double sw0 = row_below[i0 - 1];
            const double V_S0 = (east0 - west0) * inv_two_ds;
            const double V_SS0 = (east0 - 2.0 * V0 + west0) * inv_ds_sq;
            const double V_v0 = (north0 - south0) * inv_two_dv;
            const double V_vv0 = (north0 - 2.0 * V0 + south0) * inv_dv_sq;
            const double V_Sv0 = (ne0 - nw0 - se0 + sw0) * inv_four_ds_dv;
            out_row[i0] = V0 + dt * (half_v * stock_sq[i0] * V_SS0 +
                                     rho_xi_v * stock[i0] * V_Sv0 +
                                     half_xi2_v * V_vv0 +
                                     carry * stock[i0] * V_S0 +
                                     mean_rev * V_v0 - rate * V0);

            const double V1 = row_mid[i1];
            const double east1 = row_mid[i1 + 1];
            const double west1 = row_mid[i1 - 1];
            const double north1 = row_above[i1];
            const double south1 = row_below[i1];
            const double ne1 = row_above[i1 + 1];
            const double nw1 = row_above[i1 - 1];
            const double se1 = row_below[i1 + 1];
            const double sw1 = row_below[i1 - 1];
            const double V_S1 = (east1 - west1) * inv_two_ds;
            const double V_SS1 = (east1 - 2.0 * V1 + west1) * inv_ds_sq;
            const double V_v1 = (north1 - south1) * inv_two_dv;
            const double V_vv1 = (north1 - 2.0 * V1 + south1) * inv_dv_sq;
            const double V_Sv1 = (ne1 - nw1 - se1 + sw1) * inv_four_ds_dv;
            out_row[i1] = V1 + dt * (half_v * stock_sq[i1] * V_SS1 +
                                     rho_xi_v * stock[i1] * V_Sv1 +
                                     half_xi2_v * V_vv1 +
                                     carry * stock[i1] * V_S1 +
                                     mean_rev * V_v1 - rate * V1);

            const double V2 = row_mid[i2];
            const double east2 = row_mid[i2 + 1];
            const double west2 = row_mid[i2 - 1];
            const double north2 = row_above[i2];
            const double south2 = row_below[i2];
            const double ne2 = row_above[i2 + 1];
            const double nw2 = row_above[i2 - 1];
            const double se2 = row_below[i2 + 1];
            const double sw2 = row_below[i2 - 1];
            const double V_S2 = (east2 - west2) * inv_two_ds;
            const double V_SS2 = (east2 - 2.0 * V2 + west2) * inv_ds_sq;
            const double V_v2 = (north2 - south2) * inv_two_dv;
            const double V_vv2 = (north2 - 2.0 * V2 + south2) * inv_dv_sq;
            const double V_Sv2 = (ne2 - nw2 - se2 + sw2) * inv_four_ds_dv;
            out_row[i2] = V2 + dt * (half_v * stock_sq[i2] * V_SS2 +
                                     rho_xi_v * stock[i2] * V_Sv2 +
                                     half_xi2_v * V_vv2 +
                                     carry * stock[i2] * V_S2 +
                                     mean_rev * V_v2 - rate * V2);

            const double V3 = row_mid[i3];
            const double east3 = row_mid[i3 + 1];
            const double west3 = row_mid[i3 - 1];
            const double north3 = row_above[i3];
            const double south3 = row_below[i3];
            const double ne3 = row_above[i3 + 1];
            const double nw3 = row_above[i3 - 1];
            const double se3 = row_below[i3 + 1];
            const double sw3 = row_below[i3 - 1];
            const double V_S3 = (east3 - west3) * inv_two_ds;
            const double V_SS3 = (east3 - 2.0 * V3 + west3) * inv_ds_sq;
            const double V_v3 = (north3 - south3) * inv_two_dv;
            const double V_vv3 = (north3 - 2.0 * V3 + south3) * inv_dv_sq;
            const double V_Sv3 = (ne3 - nw3 - se3 + sw3) * inv_four_ds_dv;
            out_row[i3] = V3 + dt * (half_v * stock_sq[i3] * V_SS3 +
                                     rho_xi_v * stock[i3] * V_Sv3 +
                                     half_xi2_v * V_vv3 +
                                     carry * stock[i3] * V_S3 +
                                     mean_rev * V_v3 - rate * V3);
        }
        // Whatever is left over when the row is not a multiple of four.
        for (; stock_i < last_stock; ++stock_i) {
            const double V = row_mid[stock_i];
            const double east = row_mid[stock_i + 1];
            const double west = row_mid[stock_i - 1];
            const double north = row_above[stock_i];
            const double south = row_below[stock_i];
            const double ne = row_above[stock_i + 1];
            const double nw = row_above[stock_i - 1];
            const double se = row_below[stock_i + 1];
            const double sw = row_below[stock_i - 1];
            const double V_S = (east - west) * inv_two_ds;
            const double V_SS = (east - 2.0 * V + west) * inv_ds_sq;
            const double V_v = (north - south) * inv_two_dv;
            const double V_vv = (north - 2.0 * V + south) * inv_dv_sq;
            const double V_Sv = (ne - nw - se + sw) * inv_four_ds_dv;
            out_row[stock_i] =
                V + dt * (half_v * stock_sq[stock_i] * V_SS +
                          rho_xi_v * stock[stock_i] * V_Sv +
                          half_xi2_v * V_vv + carry * stock[stock_i] * V_S +
                          mean_rev * V_v - rate * V);
        }
    }

    // The v=0 row does three flops per cell and is not worth unrolling. It is
    // left alone so that the diff against level 5 shows only the interior.
    const double* const row_zero = cur_base;
    const double* const row_one = cur_base + num_stock_nodes;
    double* const out_zero = next_base;
    for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
        const double V = row_zero[stock_i];
        const double V_S =
            (row_zero[stock_i + 1] - row_zero[stock_i - 1]) * inv_two_ds;
        const double V_v = (row_one[stock_i] - row_zero[stock_i]) * inv_dv;
        out_zero[stock_i] = V + dt * (carry * stock[stock_i] * V_S +
                                      kappa_theta * V_v - rate * V);
    }
}

// The two negative controls, labelled opt-ctl-order and opt-ctl-branch.

// Levels 3 and 4 come out as nulls because the baseline was already written
// correctly, leaving them nothing to win back. The way to show what they are
// worth without sabotaging an earlier rung is to write the wrong version and
// time that instead.

// Both controls start from level 5 and change exactly one structural thing,
// so whatever it costs can be attributed to that thing, and both are compared
// against level 5 for the same reason. Their arithmetic is identical to the
// rungs above, so test_opt_matches holds them to the same tolerance and a
// different answer would mean a bug rather than a technique.


// The first control is level 5 with the two loops swapped, so stock runs on
// the outside and variance on the inside. Every inner iteration now jumps a
// whole row through memory, around 16 KB on the reference grid, so a 64-byte
// cache line yields one useful double instead of eight and the prefetcher has
// no contiguous stream to follow. The arithmetic, the tables and the stencil
// are all untouched. This measurement is what gives level 3 its meaning and
// what the cache animation is based on.
static void step_ctl_order(const Grid& g, const std::vector<double>& cur,
                           std::vector<double>& next, const Config& cfg,
                           double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    const double inv_two_ds = 1.0 / (2.0 * stock_spacing);
    const double inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    const double inv_two_dv = 1.0 / (2.0 * variance_spacing);
    const double inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    const double inv_four_ds_dv =
        1.0 / (4.0 * stock_spacing * variance_spacing);
    const double inv_dv = 1.0 / variance_spacing;
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    std::vector<double> stock_sq(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        const double S = stock_i * stock_spacing;
        stock[stock_i] = S;
        stock_sq[stock_i] = S * S;
    }

    const int last_stock = num_stock_nodes - 1;
    const int last_var = num_variance_nodes - 1;
    const double* const cur_base = cur.data();
    double* const next_base = next.data();

    // The order here is wrong on purpose. The row constants can no longer be
    // hoisted either, because they now change on every inner iteration, and
    // that is the second cost of getting the traversal backwards.
    for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
        for (int var_j = 1; var_j < last_var; ++var_j) {
            const double* const row_mid = cur_base + var_j * num_stock_nodes;
            const double* const row_above = row_mid + num_stock_nodes;
            const double* const row_below = row_mid - num_stock_nodes;
            double* const out_row = next_base + var_j * num_stock_nodes;

            const double v = var_j * variance_spacing;
            const double half_v = 0.5 * v;
            const double rho_xi_v = (rho * xi) * v;
            const double half_xi2_v = ((0.5 * xi) * xi) * v;
            const double mean_rev = kappa * (theta - v);

            const double V = row_mid[stock_i];
            const double east = row_mid[stock_i + 1];
            const double west = row_mid[stock_i - 1];
            const double north = row_above[stock_i];
            const double south = row_below[stock_i];
            const double ne = row_above[stock_i + 1];
            const double nw = row_above[stock_i - 1];
            const double se = row_below[stock_i + 1];
            const double sw = row_below[stock_i - 1];
            const double V_S = (east - west) * inv_two_ds;
            const double V_SS = (east - 2.0 * V + west) * inv_ds_sq;
            const double V_v = (north - south) * inv_two_dv;
            const double V_vv = (north - 2.0 * V + south) * inv_dv_sq;
            const double V_Sv = (ne - nw - se + sw) * inv_four_ds_dv;
            out_row[stock_i] =
                V + dt * (half_v * stock_sq[stock_i] * V_SS +
                          rho_xi_v * stock[stock_i] * V_Sv +
                          half_xi2_v * V_vv + carry * stock[stock_i] * V_S +
                          mean_rev * V_v - rate * V);
        }
    }

    // This row is unchanged from level 5, because the control is only meant
    // to alter how the interior is traversed.
    const double* const row_zero = cur_base;
    const double* const row_one = cur_base + num_stock_nodes;
    double* const out_zero = next_base;
    for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
        const double V = row_zero[stock_i];
        const double V_S =
            (row_zero[stock_i + 1] - row_zero[stock_i - 1]) * inv_two_ds;
        const double V_v = (row_one[stock_i] - row_zero[stock_i]) * inv_dv;
        out_zero[stock_i] = V + dt * (carry * stock[stock_i] * V_S +
                                      kappa_theta * V_v - rate * V);
    }
}

// The second control is level 5 with the v=0 row fused back in, so one loop
// covers every row and tests `if (var_j == 0)` on every single cell. That
// branch is highly predictable, since it is false for all but one row out of
// hundreds, so a modern predictor should make it nearly free. Establishing
// that is the point of the control, because it measures whether hoisting a
// conditional out of a loop still buys anything on current hardware, and a
// cost near zero is a real answer rather than a failed experiment.
static void step_ctl_branch(const Grid& g, const std::vector<double>& cur,
                            std::vector<double>& next, const Config& cfg,
                            double dt) {
    const double rate = cfg.market.rate;
    const double div_yield = cfg.market.div_yield;
    const double kappa = cfg.heston.kappa;
    const double theta = cfg.heston.theta;
    const double xi = cfg.heston.xi;
    const double rho = cfg.heston.rho;

    const int num_stock_nodes = g.num_stock_nodes();
    const int num_variance_nodes = g.num_variance_nodes();
    const double stock_spacing = g.stock_spacing();
    const double variance_spacing = g.variance_spacing();

    const double inv_two_ds = 1.0 / (2.0 * stock_spacing);
    const double inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    const double inv_two_dv = 1.0 / (2.0 * variance_spacing);
    const double inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    const double inv_four_ds_dv =
        1.0 / (4.0 * stock_spacing * variance_spacing);
    const double inv_dv = 1.0 / variance_spacing;
    const double carry = rate - div_yield;
    const double kappa_theta = kappa * theta;

    std::vector<double> stock(static_cast<std::size_t>(num_stock_nodes));
    std::vector<double> stock_sq(static_cast<std::size_t>(num_stock_nodes));
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        const double S = stock_i * stock_spacing;
        stock[stock_i] = S;
        stock_sq[stock_i] = S * S;
    }

    const int last_stock = num_stock_nodes - 1;
    const int last_var = num_variance_nodes - 1;
    const double* const cur_base = cur.data();
    double* const next_base = next.data();

    // Fused on purpose, so the equation is chosen per cell rather than once
    // per loop.
    for (int var_j = 0; var_j < last_var; ++var_j) {
        const double* const row_mid = cur_base + var_j * num_stock_nodes;
        const double* const row_above = row_mid + num_stock_nodes;
        double* const out_row = next_base + var_j * num_stock_nodes;

        const double v = var_j * variance_spacing;
        const double half_v = 0.5 * v;
        const double rho_xi_v = (rho * xi) * v;
        const double half_xi2_v = ((0.5 * xi) * xi) * v;
        const double mean_rev = kappa * (theta - v);

        for (int stock_i = 1; stock_i < last_stock; ++stock_i) {
            if (var_j == 0) {
                // There is no row below this one, which is why the branch has
                // to sit inside the loop once the two are fused. The
                // row_below pointer cannot even be formed outside the else.
                const double V = row_mid[stock_i];
                const double V_S =
                    (row_mid[stock_i + 1] - row_mid[stock_i - 1]) * inv_two_ds;
                const double V_v =
                    (row_above[stock_i] - row_mid[stock_i]) * inv_dv;
                out_row[stock_i] = V + dt * (carry * stock[stock_i] * V_S +
                                             kappa_theta * V_v - rate * V);
            } else {
                const double* const row_below = row_mid - num_stock_nodes;
                const double V = row_mid[stock_i];
                const double east = row_mid[stock_i + 1];
                const double west = row_mid[stock_i - 1];
                const double north = row_above[stock_i];
                const double south = row_below[stock_i];
                const double ne = row_above[stock_i + 1];
                const double nw = row_above[stock_i - 1];
                const double se = row_below[stock_i + 1];
                const double sw = row_below[stock_i - 1];
                const double V_S = (east - west) * inv_two_ds;
                const double V_SS = (east - 2.0 * V + west) * inv_ds_sq;
                const double V_v = (north - south) * inv_two_dv;
                const double V_vv = (north - 2.0 * V + south) * inv_dv_sq;
                const double V_Sv = (ne - nw - se + sw) * inv_four_ds_dv;
                out_row[stock_i] =
                    V + dt * (half_v * stock_sq[stock_i] * V_SS +
                              rho_xi_v * stock[stock_i] * V_Sv +
                              half_xi2_v * V_vv +
                              carry * stock[stock_i] * V_S +
                              mean_rev * V_v - rate * V);
            }
        }
    }
}

StepKernel kernel_for_level(int level) {
    switch (level) {
        case 0: return &step_level0;
        case 1: return &step_level1;
        case 2: return &step_level2;
        case 3: return &step_level3;
        case 4: return &step_level4;
        case 5: return &step_level5;
        case 6: return &step_level6;
        // Negative controls, not ladder rungs (see the block comment above).
        case 7: return &step_ctl_order;
        case 8: return &step_ctl_branch;
        default: return nullptr;
    }
}
