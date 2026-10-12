// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <chrono>
#include <cstdio>
#include <vector>

#include "io.h"
#include "solver.h"
#include "stencil.h"

// Prices the option on one core with AVX2, four cells per instruction.
// The time loop is the serial one, only the interior row function changes.
SolveResult SimdSolver::solve(const Config& cfg) {
    SolveResult result;

    // Without AVX2 the row falls back to scalar maths, so say so
    if (!avx_row_is_vectorised()) {
        std::fprintf(stderr, "[simd] built without AVX2 and FMA, using the scalar row\n");
    }

    // Fill the first sheet with the payoff at expiry
    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    // Work out every constant once, before the time loop
    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;

    // current is always the sheet we read, next is always the sheet we write.
    // These stay correct after swap_buffers because the swap moves the contents
    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        // auto lets the compiler work out the long clock type for us
        auto step_start = std::chrono::steady_clock::now();

        // Every interior row, four cells at a time within each row
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }

        // The boundaries, with the top row last because it copies the row below
        update_zero_row(constants, current, next);
        apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
        copy_top_row(constants, next);

        // The sheet we just wrote becomes the one we read next step
        grid.swap_buffers();

        auto step_end = std::chrono::steady_clock::now();
        result.seconds += std::chrono::duration<double>(step_end - step_start).count();

        // Snapshots are written after the timer stops so file writing is not timed
        if (cfg.dump_every > 0 && step % cfg.dump_every == 0) {
            dump_snapshot(cfg.dump_dir, step, grid);
        }
    }

    extract_result(grid, cfg, result);

    // Cells updated per second, so different grid sizes compare fairly
    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
