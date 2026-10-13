// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <omp.h>

#include <vector>

#include "io.h"
#include "solver.h"
#include "stencil.h"

// Prices the option on every core. Each thread updates its own band of rows.
// Inside each row the AVX2 function does four cells per instruction.
SolveResult OmpSolver::solve(const Config& cfg) {
    SolveResult result;

    // Record the thread count set by OMP_NUM_THREADS.
    result.threads = omp_get_max_threads();

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

        // Split the interior rows between the threads, one equal band each.
        // Threads only read current and only write their own rows of next, so no locks are needed.
        #pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }

        // The parallel for ends with every thread waiting for the others.
        // Only the main thread runs the boundaries.
        // The top row goes last because it copies the row below.
        double boundary_start = omp_get_wtime();
        update_zero_row(constants, current, next);
        apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
        copy_top_row(constants, next);

        // Swap the buffers for the next step.
        grid.swap_buffers();
        double boundary_end = omp_get_wtime();
        result.boundary_seconds += boundary_end - boundary_start;

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
