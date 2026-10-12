// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <vector>

#include "solver.h"
#include "solver_experiments.h"

// A weight above this means the quote point sits between nodes
static const double kOffNodeWarnCells = 1e-9;

// Blends the four grid values around a point that sits between nodes.
// The price and every Greek use this, so they all come off the same surface.
static double blend_cell(const Grid& grid, const std::vector<double>& values, int base_i, int base_j, double stock_weight, double variance_weight) {
    double stock_lower = 1.0 - stock_weight;
    double variance_lower = 1.0 - variance_weight;

    // The four corners around the point
    double bottom_left = values[grid.index(base_i, base_j)];
    double bottom_right = values[grid.index(base_i + 1, base_j)];
    double top_left = values[grid.index(base_i, base_j + 1)];
    double top_right = values[grid.index(base_i + 1, base_j + 1)];

    // Blend along stock on the bottom and top rows, then blend those along variance
    return variance_lower * (stock_lower * bottom_left + stock_weight * bottom_right) + variance_weight * (stock_lower * top_left + stock_weight * top_right);
}

// Reads the price, delta, gamma and vega off the finished sheet.
// Called by every solver once the time loop is done.
void Solver::extract_result(const Grid& grid, const Config& cfg, SolveResult& out) const {
    const std::vector<double>& values = grid.current();

    // Where today's stock price and variance sit on the grid
    AxisPosition spot_position = grid.stock_position(cfg.market.spot);
    AxisPosition v0_position = grid.variance_position(cfg.heston.v0);
    int base_i = spot_position.cell;
    int base_j = v0_position.cell;
    double stock_weight = spot_position.weight;
    double variance_weight = v0_position.weight;

    // Warn on stderr so stdout stays one clean CSV line
    if (stock_weight > kOffNodeWarnCells || variance_weight > kOffNodeWarnCells) {
        std::fprintf(stderr,
                     "[%s] off-node quote, spot sits %.4f cells past node %d "
                     "and v0 sits %.4f cells past node %d, so price and Greeks "
                     "are bilinearly interpolated\n",
                     name().c_str(), stock_weight, base_i, variance_weight, base_j);
    }

    // The price at the point, and the same blend one cell to each side
    out.price = blend_cell(grid, values, base_i, base_j, stock_weight, variance_weight);
    double east = blend_cell(grid, values, base_i + 1, base_j, stock_weight, variance_weight);
    double west = blend_cell(grid, values, base_i - 1, base_j, stock_weight, variance_weight);
    double north = blend_cell(grid, values, base_i, base_j + 1, stock_weight, variance_weight);
    double south = blend_cell(grid, values, base_i, base_j - 1, stock_weight, variance_weight);

    // Central differences give the Greeks
    out.delta = (east - west) / (2.0 * grid.stock_spacing());
    out.gamma = (east - 2.0 * out.price + west) / (grid.stock_spacing() * grid.stock_spacing());
    // Vega per unit of variance, not per unit of volatility
    out.vega = (north - south) / (2.0 * grid.variance_spacing());
}

// Works out the largest stable timestep and prints it next to the one we use.
// Also says whether the Feller condition holds, which is for the write up only.
void Solver::report_stability(const Grid& grid, const Config& cfg, double dt, SolveResult& out) const {
    double rate = cfg.market.rate;
    double kappa = cfg.heston.kappa;
    double theta = cfg.heston.theta;
    double xi = cfg.heston.xi;

    // The worst cell is the far corner, where the diffusion is largest
    double max_stock_price = grid.stock_price(grid.num_stock_nodes() - 1);
    double max_variance = grid.variance(grid.num_variance_nodes() - 1);
    double stock_spacing = grid.stock_spacing();
    double variance_spacing = grid.variance_spacing();
    out.dt_stable_estimate = 1.0 / (max_variance * max_stock_price * max_stock_price / (stock_spacing * stock_spacing) + xi * xi * max_variance / (variance_spacing * variance_spacing) + rate);

    // Feller holds when the variance can never reach zero
    bool feller_ok = 2.0 * kappa * theta >= xi * xi;
    const char* feller_text = "VIOLATED";
    if (feller_ok) {
        feller_text = "OK";
    }
    std::fprintf(stderr, "[%s] feller %s (2*kappa*theta=%.4f, xi^2=%.4f)  dt=%.3e  dt_stable~%.3e\n",
                 name().c_str(), feller_text, 2.0 * kappa * theta, xi * xi, dt, out.dt_stable_estimate);
}

// Creates the solver with the given name and returns it.
// The unique_ptr owns the solver and deletes it when it goes out of scope.
std::unique_ptr<Solver> make_solver(const std::string& name) {
    if (name == "baseline") {
        return std::make_unique<BaselineSolver>();
    }
    if (name == "serial") {
        return std::make_unique<SerialSolver>();
    }
    if (name == "simd") {
        return std::make_unique<SimdSolver>();
    }
    if (name == "omp") {
        return std::make_unique<OmpSolver>();
    }
    if (name == "mpi") {
        return std::make_unique<MpiSolver>();
    }
    // The measured experiments around the omp rung
    if (name == "omp_scalar") {
        return std::make_unique<OmpScalarSolver>();
    }
    if (name == "omp_persistent") {
        return std::make_unique<OmpPersistentSolver>();
    }
    if (name == "omp_runtime") {
        return std::make_unique<OmpRuntimeScheduleSolver>();
    }
    if (name == "omp_edges") {
        return std::make_unique<OmpParallelEdgesSolver>();
    }
    if (name == "omp_float") {
        return std::make_unique<OmpFloatSolver>();
    }
    throw std::runtime_error("unknown solver: " + name);
}
