// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <omp.h>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "io.h"
#include "solver_experiments.h"
#include "stencil.h"

// Updates the bottom row in float, where the variance is zero.
// Same sum as update_zero_row, in float.
static void update_zero_row_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next) {
    const float* row_zero = current.data();
    const float* row_one = current.data() + constants.num_stock_nodes;
    float* out_row = next.data();
    float variance_drift = constants.variance_drift[0];

    for (int stock_i = 1; stock_i < constants.num_stock_nodes - 1; ++stock_i) {
        float centre = row_zero[stock_i];
        float dV_dS = (row_zero[stock_i + 1] - row_zero[stock_i - 1]) * constants.inv_two_ds;
        float dV_dv = (row_one[stock_i] - centre) * constants.inv_dv;
        float stock_drift_term = constants.stock_drift * constants.stock_price[stock_i] * dV_dS;
        float variance_drift_term = variance_drift * dV_dv;
        float change = stock_drift_term + variance_drift_term - constants.rate * centre;
        out_row[stock_i] = centre + constants.dt * change;
    }
}

// Writes the known values at S = 0 and at the largest S into every row, then copies the top row.
// The edge values are worked out in double, then rounded to float.
static void apply_edges_and_top_float(const StencilConstants& constants, std::vector<float>& next, double time_to_expiry) {
    double strike_today = constants.strike * std::exp(-constants.rate * time_to_expiry);
    double zero_stock_value = 0.0;
    double max_stock_value = 0.0;
    if (constants.is_call) {
        // A call is worthless at S = 0 and worth S minus the discounted strike at the top.
        zero_stock_value = 0.0;
        max_stock_value = constants.max_stock_price * std::exp(-constants.div_yield * time_to_expiry) - strike_today;
    } else {
        // A put is worth the discounted strike at S = 0 and nothing at the top.
        zero_stock_value = strike_today;
        max_stock_value = 0.0;
    }

    int num_stock_nodes = constants.num_stock_nodes;
    int last_stock = num_stock_nodes - 1;
    for (int var_j = 0; var_j < constants.num_variance_nodes; ++var_j) {
        next[var_j * num_stock_nodes] = (float)zero_stock_value;
        next[var_j * num_stock_nodes + last_stock] = (float)max_stock_value;
    }

    // The top row copies the row below it.
    int top_row = constants.num_variance_nodes - 1;
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        next[top_row * num_stock_nodes + stock_i] = next[(top_row - 1) * num_stock_nodes + stock_i];
    }
}

// Prices the option like the omp solver with both buffers in float.
// One AVX2 register holds eight floats, and each buffer takes half the memory.
SolveResult OmpFloatSolver::solve(const Config& cfg) {
    SolveResult result;
    result.threads = omp_get_max_threads();

    // Snapshots are not supported in float.
    if (cfg.dump_every > 0) {
        throw std::runtime_error("--dump-every is not supported by the omp_float solver");
    }

    // Build the payoff and the constants in double.
    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    StencilConstants constants(grid, cfg, dt);
    FloatStencilConstants float_constants(constants);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;

    // Round the payoff to float and make a second float buffer.
    int total_cells = num_stock_nodes * num_variance_nodes;
    std::vector<float> current(total_cells, 0.0f);
    std::vector<float> next(total_cells, 0.0f);
    for (int cell = 0; cell < total_cells; ++cell) {
        current[cell] = (float)grid.current()[cell];
    }

    for (int step = 1; step <= num_timesteps; ++step) {
        double step_start = omp_get_wtime();

        // Split the interior rows between the threads.
        #pragma omp parallel for schedule(static) default(none) shared(float_constants, current, next, num_variance_nodes)
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_avx_float(float_constants, current, next, var_j);
        }

        // Boundaries on the main thread, top row last.
        double boundary_start = omp_get_wtime();
        update_zero_row_float(float_constants, current, next);
        apply_edges_and_top_float(constants, next, step * dt);

        // std::swap swaps the two vectors without copying the numbers.
        std::swap(current, next);
        double step_end = omp_get_wtime();
        result.boundary_seconds += step_end - boundary_start;
        result.seconds += step_end - step_start;
    }

    // Copy the float answer into the grid's next buffer and make it current.
    // extract_result then reads the price and Greeks in double.
    std::vector<double>& answer = grid.next();
    for (int cell = 0; cell < total_cells; ++cell) {
        answer[cell] = current[cell];
    }
    grid.swap_buffers();
    extract_result(grid, cfg, result);

    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
