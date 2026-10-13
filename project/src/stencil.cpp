// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "stencil.h"

#include <cmath>
#include <vector>

// Fills in every constant the update needs, once before the time loop.
// Moves the divisions and repeated multiplications out of the hot loop.
StencilConstants::StencilConstants(const Grid& grid, const Config& cfg, double timestep) {
    num_stock_nodes = grid.num_stock_nodes();
    num_variance_nodes = grid.num_variance_nodes();

    dt = timestep;
    rate = cfg.market.rate;
    div_yield = cfg.market.div_yield;
    stock_drift = cfg.market.rate - cfg.market.div_yield;

    is_call = cfg.option.is_call;
    strike = cfg.option.strike;
    max_stock_price = grid.stock_price(grid.num_stock_nodes() - 1);

    // One division each here, once per solve.
    double stock_spacing = grid.stock_spacing();
    double variance_spacing = grid.variance_spacing();
    inv_two_ds = 1.0 / (2.0 * stock_spacing);
    inv_ds_sq = 1.0 / (stock_spacing * stock_spacing);
    inv_two_dv = 1.0 / (2.0 * variance_spacing);
    inv_dv_sq = 1.0 / (variance_spacing * variance_spacing);
    inv_four_ds_dv = 1.0 / (4.0 * stock_spacing * variance_spacing);
    inv_dv = 1.0 / variance_spacing;

    // Stock price and its square for every column.
    stock_price.resize(num_stock_nodes);
    stock_price_sq.resize(num_stock_nodes);
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        double price = stock_i * stock_spacing;
        stock_price[stock_i] = price;
        stock_price_sq[stock_i] = price * price;
    }

    // The four PDE coefficients for every variance row.
    double kappa = cfg.heston.kappa;
    double theta = cfg.heston.theta;
    double xi = cfg.heston.xi;
    double rho = cfg.heston.rho;
    stock_diffusion.resize(num_variance_nodes);
    cross_diffusion.resize(num_variance_nodes);
    variance_diffusion.resize(num_variance_nodes);
    variance_drift.resize(num_variance_nodes);
    for (int var_j = 0; var_j < num_variance_nodes; ++var_j) {
        double variance = var_j * variance_spacing;
        stock_diffusion[var_j] = 0.5 * variance;
        cross_diffusion[var_j] = (rho * xi) * variance;
        variance_diffusion[var_j] = ((0.5 * xi) * xi) * variance;
        variance_drift[var_j] = kappa * (theta - variance);
    }
}

// Updates the cells first_stock to end_stock - 1 in one variance row of the next buffer.
// Each new value comes from the cell and its eight neighbours in the current buffer.
void update_cells_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j, int first_stock, int end_stock) {
    // Copy the constants into local variables first.
    // Otherwise the compiler may reload them from memory after every write.
    double dt = constants.dt;
    double rate = constants.rate;
    double stock_drift = constants.stock_drift;
    double inv_two_ds = constants.inv_two_ds;
    double inv_ds_sq = constants.inv_ds_sq;
    double inv_two_dv = constants.inv_two_dv;
    double inv_dv_sq = constants.inv_dv_sq;
    double inv_four_ds_dv = constants.inv_four_ds_dv;
    int num_stock_nodes = constants.num_stock_nodes;

    // This row's four PDE coefficients.
    double stock_diffusion = constants.stock_diffusion[var_j];
    double cross_diffusion = constants.cross_diffusion[var_j];
    double variance_diffusion = constants.variance_diffusion[var_j];
    double variance_drift = constants.variance_drift[var_j];

    // Pointers to the start of the rows to read and write.
    // The vectors still own the memory.
    const double* stock_price = constants.stock_price.data();
    const double* stock_price_sq = constants.stock_price_sq.data();
    const double* row_centre = current.data() + var_j * num_stock_nodes;
    const double* row_above = row_centre + num_stock_nodes;
    const double* row_below = row_centre - num_stock_nodes;
    double* out_row = next.data() + var_j * num_stock_nodes;

    for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {
        // Read the cell and its eight neighbours.
        // East and west are one stock node away. North and south are one variance row away.
        double centre = row_centre[stock_i];
        double east = row_centre[stock_i + 1];
        double west = row_centre[stock_i - 1];
        double north = row_above[stock_i];
        double south = row_below[stock_i];
        double north_east = row_above[stock_i + 1];
        double north_west = row_above[stock_i - 1];
        double south_east = row_below[stock_i + 1];
        double south_west = row_below[stock_i - 1];

        // Central differences give each derivative of the option value V.
        double dV_dS = (east - west) * inv_two_ds;
        double d2V_dS2 = (east - 2.0 * centre + west) * inv_ds_sq;
        double dV_dv = (north - south) * inv_two_dv;
        double d2V_dv2 = (north - 2.0 * centre + south) * inv_dv_sq;
        double d2V_dS_dv = (north_east - north_west - south_east + south_west) * inv_four_ds_dv;

        // The five terms of the Heston PDE.
        double stock_diffusion_term = stock_diffusion * stock_price_sq[stock_i] * d2V_dS2;
        double cross_term = cross_diffusion * stock_price[stock_i] * d2V_dS_dv;
        double variance_diffusion_term = variance_diffusion * d2V_dv2;
        double stock_drift_term = stock_drift * stock_price[stock_i] * dV_dS;
        double variance_drift_term = variance_drift * dV_dv;

        // Step back one dt from the current value.
        double change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - rate * centre;
        out_row[stock_i] = centre + dt * change;
    }
}

// Updates every interior cell in one variance row of the next buffer.
// The first and last columns are the stock edges and are skipped.
void update_row_scalar(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next, int var_j) {
    update_cells_scalar(constants, current, next, var_j, 1, constants.num_stock_nodes - 1);
}

// Updates the bottom row, where the variance is zero.
// Every diffusion term has a factor of v, so only the drift terms are left.
void update_zero_row(const StencilConstants& constants, const std::vector<double>& current, std::vector<double>& next) {
    const double* row_zero = current.data();
    const double* row_one = current.data() + constants.num_stock_nodes;
    double* out_row = next.data();

    // At v = 0 the variance drift is kappa times theta.
    double variance_drift = constants.variance_drift[0];

    for (int stock_i = 1; stock_i < constants.num_stock_nodes - 1; ++stock_i) {
        double centre = row_zero[stock_i];
        double dV_dS = (row_zero[stock_i + 1] - row_zero[stock_i - 1]) * constants.inv_two_ds;
        // Look one row up, because there is no row below zero.
        double dV_dv = (row_one[stock_i] - centre) * constants.inv_dv;

        double stock_drift_term = constants.stock_drift * constants.stock_price[stock_i] * dV_dS;
        double variance_drift_term = variance_drift * dV_dv;
        double change = stock_drift_term + variance_drift_term - constants.rate * centre;
        out_row[stock_i] = centre + constants.dt * change;
    }
}

// Writes the known option values in the first and last stock columns.
// Only rows row_begin to row_end - 1 are written.
void apply_stock_edges(const StencilConstants& constants, std::vector<double>& next, int row_begin, int row_end, double time_to_expiry) {
    // The strike paid at expiry, discounted back to this time.
    double strike_today = constants.strike * std::exp(-constants.rate * time_to_expiry);

    double zero_stock_value = 0.0;
    double max_stock_value = 0.0;
    if (constants.is_call) {
        // A call is worthless at S = 0 and acts like a forward at the top.
        zero_stock_value = 0.0;
        max_stock_value = constants.max_stock_price * std::exp(-constants.div_yield * time_to_expiry) - strike_today;
    } else {
        // A put pays the discounted strike at S = 0 and is worthless at the top.
        zero_stock_value = strike_today;
        max_stock_value = 0.0;
    }

    int last_stock = constants.num_stock_nodes - 1;
    for (int var_j = row_begin; var_j < row_end; ++var_j) {
        next[var_j * constants.num_stock_nodes] = zero_stock_value;
        next[var_j * constants.num_stock_nodes + last_stock] = max_stock_value;
    }
}

// Sets the top variance row equal to the row just below it.
// The value stops changing with variance at the top of the grid.
void copy_top_row(const StencilConstants& constants, std::vector<double>& next) {
    // Must run after the row below has been updated.
    int num_stock_nodes = constants.num_stock_nodes;
    int top_row = constants.num_variance_nodes - 1;
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        next[top_row * num_stock_nodes + stock_i] = next[(top_row - 1) * num_stock_nodes + stock_i];
    }
}

// Rounds every constant the update needs from double down to float.
// Both solvers start from the same double constants.
FloatStencilConstants::FloatStencilConstants(const StencilConstants& constants) {
    num_stock_nodes = constants.num_stock_nodes;
    num_variance_nodes = constants.num_variance_nodes;

    // (float) rounds a double to the nearest float.
    dt = (float)constants.dt;
    rate = (float)constants.rate;
    stock_drift = (float)constants.stock_drift;
    inv_two_ds = (float)constants.inv_two_ds;
    inv_ds_sq = (float)constants.inv_ds_sq;
    inv_two_dv = (float)constants.inv_two_dv;
    inv_dv_sq = (float)constants.inv_dv_sq;
    inv_four_ds_dv = (float)constants.inv_four_ds_dv;
    inv_dv = (float)constants.inv_dv;

    // Column tables, one float per stock column.
    stock_price.resize(num_stock_nodes);
    stock_price_sq.resize(num_stock_nodes);
    for (int stock_i = 0; stock_i < num_stock_nodes; ++stock_i) {
        stock_price[stock_i] = (float)constants.stock_price[stock_i];
        stock_price_sq[stock_i] = (float)constants.stock_price_sq[stock_i];
    }

    // Row tables, one float per variance row.
    stock_diffusion.resize(num_variance_nodes);
    cross_diffusion.resize(num_variance_nodes);
    variance_diffusion.resize(num_variance_nodes);
    variance_drift.resize(num_variance_nodes);
    for (int var_j = 0; var_j < num_variance_nodes; ++var_j) {
        stock_diffusion[var_j] = (float)constants.stock_diffusion[var_j];
        cross_diffusion[var_j] = (float)constants.cross_diffusion[var_j];
        variance_diffusion[var_j] = (float)constants.variance_diffusion[var_j];
        variance_drift[var_j] = (float)constants.variance_drift[var_j];
    }
}

// Updates the cells first_stock to end_stock - 1 of one variance row in float.
// Line for line the same sum as update_cells_scalar, in float.
void update_cells_scalar_float(const FloatStencilConstants& constants, const std::vector<float>& current, std::vector<float>& next, int var_j, int first_stock, int end_stock) {
    int num_stock_nodes = constants.num_stock_nodes;

    // This row's four PDE coefficients.
    float stock_diffusion = constants.stock_diffusion[var_j];
    float cross_diffusion = constants.cross_diffusion[var_j];
    float variance_diffusion = constants.variance_diffusion[var_j];
    float variance_drift = constants.variance_drift[var_j];

    // Pointers to the start of the rows to read and write.
    const float* row_centre = current.data() + var_j * num_stock_nodes;
    const float* row_above = row_centre + num_stock_nodes;
    const float* row_below = row_centre - num_stock_nodes;
    float* out_row = next.data() + var_j * num_stock_nodes;

    for (int stock_i = first_stock; stock_i < end_stock; ++stock_i) {
        // Read the cell and its eight neighbours.
        float centre = row_centre[stock_i];
        float east = row_centre[stock_i + 1];
        float west = row_centre[stock_i - 1];
        float north = row_above[stock_i];
        float south = row_below[stock_i];
        float north_east = row_above[stock_i + 1];
        float north_west = row_above[stock_i - 1];
        float south_east = row_below[stock_i + 1];
        float south_west = row_below[stock_i - 1];

        // Central differences give each derivative of the option value V.
        // 2.0f is a float literal. A plain 2.0 would make the sum a double.
        float dV_dS = (east - west) * constants.inv_two_ds;
        float d2V_dS2 = (east - 2.0f * centre + west) * constants.inv_ds_sq;
        float dV_dv = (north - south) * constants.inv_two_dv;
        float d2V_dv2 = (north - 2.0f * centre + south) * constants.inv_dv_sq;
        float d2V_dS_dv = (north_east - north_west - south_east + south_west) * constants.inv_four_ds_dv;

        // The five terms of the Heston PDE.
        float stock_diffusion_term = stock_diffusion * constants.stock_price_sq[stock_i] * d2V_dS2;
        float cross_term = cross_diffusion * constants.stock_price[stock_i] * d2V_dS_dv;
        float variance_diffusion_term = variance_diffusion * d2V_dv2;
        float stock_drift_term = constants.stock_drift * constants.stock_price[stock_i] * dV_dS;
        float variance_drift_term = variance_drift * dV_dv;

        // Step back one dt from the current value.
        float change = stock_diffusion_term + cross_term + variance_diffusion_term + stock_drift_term + variance_drift_term - constants.rate * centre;
        out_row[stock_i] = centre + constants.dt * change;
    }
}
