// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

// OpenMP is only used here for the omp_get_wtime timer.
#include <omp.h>

#include <cmath>
#include <vector>

#include "io.h"
#include "solver.h"

// Prices the option with all the maths written out in one place.
// Every other solver must match its answers.
SolveResult BaselineSolver::solve(const Config& cfg) {
    SolveResult result;

    // Fill the current buffer with the payoff at expiry.
    Grid grid(cfg.grid, cfg.option);
    grid.init_payoff(cfg.option);

    // Short local names for the settings.
    double strike = cfg.option.strike;
    double rate = cfg.market.rate;
    double div_yield = cfg.market.div_yield;
    double kappa = cfg.heston.kappa;
    double theta = cfg.heston.theta;
    double xi = cfg.heston.xi;
    double rho = cfg.heston.rho;
    int num_timesteps = cfg.grid.num_timesteps;

    double dt = cfg.option.maturity_years / num_timesteps;
    report_stability(grid, cfg, dt, result);

    // Largest stock price on the grid, for the call's top edge.
    double max_stock_price = grid.stock_price(grid.num_stock_nodes() - 1);

    // Read from current and write to next.
    // swap_buffers swaps the contents of these two vectors.
    const std::vector<double>& current = grid.current();
    std::vector<double>& next = grid.next();

    for (int step = 1; step <= num_timesteps; ++step) {
        // Start the step timer. omp_get_wtime returns seconds.
        double step_start = omp_get_wtime();

        // Apply the full PDE to every interior cell.
        for (int var_j = 1; var_j < grid.num_variance_nodes() - 1; ++var_j) {
            for (int stock_i = 1; stock_i < grid.num_stock_nodes() - 1; ++stock_i) {
                // Stock price and variance at this cell.
                double stock_price = grid.stock_price(stock_i);
                double variance = grid.variance(var_j);

                // Read the cell and its eight neighbours.
                // East and west are one stock node away. North and south are one variance row away.
                double centre = current[grid.index(stock_i, var_j)];
                double east = current[grid.index(stock_i + 1, var_j)];
                double west = current[grid.index(stock_i - 1, var_j)];
                double north = current[grid.index(stock_i, var_j + 1)];
                double south = current[grid.index(stock_i, var_j - 1)];
                double north_east = current[grid.index(stock_i + 1, var_j + 1)];
                double north_west = current[grid.index(stock_i - 1, var_j + 1)];
                double south_east = current[grid.index(stock_i + 1, var_j - 1)];
                double south_west = current[grid.index(stock_i - 1, var_j - 1)];

                // Central differences give each derivative of the option value V.
                double dV_dS = (east - west) / (2.0 * grid.stock_spacing());
                double d2V_dS2 = (east - 2.0 * centre + west) / (grid.stock_spacing() * grid.stock_spacing());
                double dV_dv = (north - south) / (2.0 * grid.variance_spacing());
                double d2V_dv2 = (north - 2.0 * centre + south) / (grid.variance_spacing() * grid.variance_spacing());
                double d2V_dS_dv = (north_east - north_west - south_east + south_west) / (4.0 * grid.stock_spacing() * grid.variance_spacing());

                // One explicit Euler step of the Heston PDE.
                double updated = centre + dt * (0.5 * variance * stock_price * stock_price * d2V_dS2
                                                + rho * xi * variance * stock_price * d2V_dS_dv
                                                + 0.5 * xi * xi * variance * d2V_dv2
                                                + (rate - div_yield) * stock_price * dV_dS
                                                + kappa * (theta - variance) * dV_dv
                                                - rate * centre);
                next[grid.index(stock_i, var_j)] = updated;
            }
        }

        // The bottom row, where the variance is zero.
        // Every diffusion term has a factor of v, so only the drift terms are left.
        for (int stock_i = 1; stock_i < grid.num_stock_nodes() - 1; ++stock_i) {
            double stock_price = grid.stock_price(stock_i);
            double centre = current[grid.index(stock_i, 0)];
            double dV_dS = (current[grid.index(stock_i + 1, 0)] - current[grid.index(stock_i - 1, 0)]) / (2.0 * grid.stock_spacing());
            // Look one row up, because there is no row below zero.
            double dV_dv = (current[grid.index(stock_i, 1)] - centre) / grid.variance_spacing();
            double updated = centre + dt * ((rate - div_yield) * stock_price * dV_dS + kappa * theta * dV_dv - rate * centre);
            next[grid.index(stock_i, 0)] = updated;
        }

        // Time left until expiry.
        double time_to_expiry = step * dt;

        // The first and last stock columns have known values.
        double zero_stock_value = 0.0;
        double max_stock_value = 0.0;
        if (cfg.option.is_call) {
            // A call is worthless at S = 0 and acts like a forward at the top.
            zero_stock_value = 0.0;
            max_stock_value = max_stock_price * std::exp(-div_yield * time_to_expiry) - strike * std::exp(-rate * time_to_expiry);
        } else {
            // A put pays the discounted strike at S = 0 and is worthless at the top.
            zero_stock_value = strike * std::exp(-rate * time_to_expiry);
            max_stock_value = 0.0;
        }
        for (int var_j = 0; var_j < grid.num_variance_nodes(); ++var_j) {
            next[grid.index(0, var_j)] = zero_stock_value;
            next[grid.index(grid.num_stock_nodes() - 1, var_j)] = max_stock_value;
        }

        // The top row copies the row below it.
        // Must run after the row below is updated.
        int top_row = grid.num_variance_nodes() - 1;
        for (int stock_i = 0; stock_i < grid.num_stock_nodes(); ++stock_i) {
            next[grid.index(stock_i, top_row)] = next[grid.index(stock_i, top_row - 1)];
        }

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
    double total_cell_updates = (double)grid.num_stock_nodes() * grid.num_variance_nodes() * num_timesteps;
    result.cell_updates_per_sec = total_cell_updates / result.seconds;
    return result;
}
