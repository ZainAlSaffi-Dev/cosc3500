// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <omp.h>

#include <chrono>
#include <vector>

#include "solver_experiments.h"
#include "stencil.h"

// These solvers are the measured experiments around the omp rung.
// Each one is the omp time loop with exactly one thing changed.

// The omp solver with update_row_scalar in place of update_row_avx.
// Everything else is the same as OmpSolver.
SolveResult OmpScalarSolver::solve(const Config& cfg) {
    SolveResult result;
    result.threads = omp_get_max_threads();

    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);
    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;
    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        auto step_start = std::chrono::steady_clock::now();

        // The only change, one cell at a time inside each row
        #pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes)
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_scalar(constants, current, next, var_j);
        }

        auto boundary_start = std::chrono::steady_clock::now();
        update_zero_row(constants, current, next);
        apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
        copy_top_row(constants, next);
        grid.swap_buffers();
        auto step_end = std::chrono::steady_clock::now();
        result.boundary_seconds += std::chrono::duration<double>(step_end - boundary_start).count();
        result.seconds += std::chrono::duration<double>(step_end - step_start).count();
    }

    extract_result(grid, cfg, result);
    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}

// The omp solver with the threads started once, before the time loop.
// Each step the threads split the rows, then one of them does the boundaries.
SolveResult OmpPersistentSolver::solve(const Config& cfg) {
    SolveResult result;
    result.threads = omp_get_max_threads();

    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);
    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;
    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    // Timed as a whole, because every thread is inside the loop the whole time
    auto loop_start = std::chrono::steady_clock::now();

    // Start the team once. Every thread runs the step loop below
    #pragma omp parallel default(none) shared(constants, current, next, grid, num_variance_nodes, num_timesteps, dt)
    {
        for (int step = 1; step <= num_timesteps; ++step) {
            // Share the rows out, then every thread waits at the end of the for
            #pragma omp for schedule(static)
            for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
                update_row_avx(constants, current, next, var_j);
            }

            // One thread does the boundaries and the swap, the others wait at its end
            #pragma omp single
            {
                update_zero_row(constants, current, next);
                apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
                copy_top_row(constants, next);
                grid.swap_buffers();
            }
        }
    }

    auto loop_end = std::chrono::steady_clock::now();
    result.seconds = std::chrono::duration<double>(loop_end - loop_start).count();

    extract_result(grid, cfg, result);
    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}

// The omp solver with schedule(runtime), so OMP_SCHEDULE picks the schedule.
// For example OMP_SCHEDULE=dynamic,4 hands out four rows at a time on demand.
SolveResult OmpRuntimeScheduleSolver::solve(const Config& cfg) {
    SolveResult result;
    result.threads = omp_get_max_threads();

    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);
    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;
    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        auto step_start = std::chrono::steady_clock::now();

        // The only change, the schedule comes from the environment
        #pragma omp parallel for schedule(runtime) default(none) shared(constants, current, next, num_variance_nodes)
        for (int var_j = 1; var_j < num_variance_nodes - 1; ++var_j) {
            update_row_avx(constants, current, next, var_j);
        }

        auto boundary_start = std::chrono::steady_clock::now();
        update_zero_row(constants, current, next);
        apply_stock_edges(constants, next, 0, num_variance_nodes, step * dt);
        copy_top_row(constants, next);
        grid.swap_buffers();
        auto step_end = std::chrono::steady_clock::now();
        result.boundary_seconds += std::chrono::duration<double>(step_end - boundary_start).count();
        result.seconds += std::chrono::duration<double>(step_end - step_start).count();
    }

    extract_result(grid, cfg, result);
    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}

// The omp solver with the v = 0 row and the stock edges moved inside the parallel loop.
// Only the top row copy is left serial, because it needs every other row finished first.
SolveResult OmpParallelEdgesSolver::solve(const Config& cfg) {
    SolveResult result;
    result.threads = omp_get_max_threads();

    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);
    int num_timesteps = cfg.grid.num_timesteps;
    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    StencilConstants constants(grid, cfg, dt);
    int num_stock_nodes = constants.num_stock_nodes;
    int num_variance_nodes = constants.num_variance_nodes;
    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        auto step_start = std::chrono::steady_clock::now();
        double time_to_expiry = step * dt;

        // Row 0 joins the loop, and each thread writes the stock edges of its own rows.
        // Neither reads anything another thread writes in this step
        #pragma omp parallel for schedule(static) default(none) shared(constants, current, next, num_variance_nodes, time_to_expiry)
        for (int var_j = 0; var_j < num_variance_nodes - 1; ++var_j) {
            if (var_j == 0) {
                update_zero_row(constants, current, next);
            } else {
                update_row_avx(constants, current, next, var_j);
            }
            apply_stock_edges(constants, next, var_j, var_j + 1, time_to_expiry);
        }

        // The top row copies the finished row below it, edges included
        auto boundary_start = std::chrono::steady_clock::now();
        copy_top_row(constants, next);
        grid.swap_buffers();
        auto step_end = std::chrono::steady_clock::now();
        result.boundary_seconds += std::chrono::duration<double>(step_end - boundary_start).count();
        result.seconds += std::chrono::duration<double>(step_end - step_start).count();
    }

    extract_result(grid, cfg, result);
    double total_cell_updates = (double)num_stock_nodes * num_variance_nodes * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
