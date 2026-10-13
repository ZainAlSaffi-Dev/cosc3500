// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

// OpenMP is only used here for the omp_get_wtime timer.
#include <omp.h>

#include <cstdio>
#include <vector>

#include "io.h"
#include "solver.h"
#include "stencil.h"

// Prices the option on one core with AVX2, four cells per instruction.
// The time loop is the serial one with a different row function.
SolveResult SimdSolver::solve(const Config& cfg) {
    SolveResult result;

    // Warn when the row falls back to scalar maths.
    if (!avx_row_is_vectorised()) {
        std::fprintf(stderr, "[simd] built without AVX2 and FMA, using the scalar row\n");
    }

    // Fill the current buffer with the payoff at expiry.
    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    // Work out every constant once before the time loop.
    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;

    // Read from current and write to next.
    // swap_buffers swaps the contents of these two vectors.
    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        // Start the step timer. omp_get_wtime returns seconds.
        double step_start = omp_get_wtime();

        // Update every interior row, four cells at a time.
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }

        // The boundaries. The top row goes last because it copies the row below.
        update_zero_row(constants, current, next);
        apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
        copy_top_row(constants, next);

        // Swap the buffers for the next step.
        grid.swap_buffers();

        double step_end = omp_get_wtime();
        result.seconds += step_end - step_start;

        // Write snapshots after the timer stops.
        if (cfg.dump_every > 0 && step % cfg.dump_every == 0) {
            dump_snapshot(cfg.dump_dir, step, grid);
        }
    }

    extract_result(grid, cfg, result);

    // Cells updated per second.
    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
