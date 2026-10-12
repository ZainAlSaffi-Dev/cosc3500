// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "grid.h"

#include <algorithm>
#include <cmath>
#include <utility>

// Sets up the grid size and spacing and makes both sheets full of zeros.
// Runs once when a Grid is created, like __init__ in Python.
Grid::Grid(const GridSpec& spec, const OptionSpec& opt) {
    num_stock_nodes_ = spec.num_stock_nodes;
    num_variance_nodes_ = spec.num_variance_nodes;

    // N nodes have N - 1 gaps between them
    stock_spacing_ = spec.stock_max_multiplier * opt.strike / (num_stock_nodes_ - 1);
    variance_spacing_ = spec.variance_max / (num_variance_nodes_ - 1);

    // Each sheet is one long block of doubles, like np.zeros(n).
    // The vectors free this memory themselves when the Grid is destroyed
    std::size_t total_cells = (std::size_t)num_stock_nodes_ * num_variance_nodes_;
    current_.assign(total_cells, 0.0);
    next_.assign(total_cells, 0.0);
}

// Swaps which sheet is read and which is written.
// Only the internal pointers are swapped, so no numbers are copied.
void Grid::swap_buffers() {
    std::swap(current_, next_);
}

// Fills the current sheet with what the option pays at expiry.
// The payoff only depends on the stock price, so every row is the same.
void Grid::init_payoff(const OptionSpec& opt) {
    for (int var_j = 0; var_j < num_variance_nodes_; ++var_j) {
        for (int stock_i = 0; stock_i < num_stock_nodes_; ++stock_i) {
            double price = stock_price(stock_i);
            double payoff = 0.0;
            if (opt.is_call) {
                payoff = std::max(price - opt.strike, 0.0);
            } else {
                payoff = std::max(opt.strike - price, 0.0);
            }
            current_[index(stock_i, var_j)] = payoff;
        }
    }
}

// A point closer than this to a node counts as sitting on the node
static const double kNodeSnapCells = 1e-9;

// Finds which cell a value falls in along one axis, and how far across it.
// The cell is kept one node away from each edge so the Greek stencils fit.
static AxisPosition locate(double coord, double spacing, int num_nodes) {
    // Position measured in cells, clamped so it never runs off the grid
    double position = std::clamp(coord / spacing, 0.0, (double)(num_nodes - 1));

    // Snap to a node if we are within rounding error of one
    double nearest_node = std::round(position);
    bool on_node = std::fabs(position - nearest_node) < kNodeSnapCells;
    double lower_node = 0.0;
    if (on_node) {
        lower_node = nearest_node;
    } else {
        lower_node = std::floor(position);
    }

    // Keep one node in from each edge
    int highest_cell = std::max(1, num_nodes - 3);
    AxisPosition out;
    out.cell = std::clamp((int)lower_node, 1, highest_cell);

    // Weight is exactly 0 on a node, otherwise the fraction past the lower node
    if (on_node && out.cell == (int)lower_node) {
        out.weight = 0.0;
    } else {
        out.weight = std::clamp(position - out.cell, 0.0, 1.0);
    }
    return out;
}

// Where today's stock price sits along the stock axis.
// Used to read the price off the finished grid.
AxisPosition Grid::stock_position(double spot) const {
    return locate(spot, stock_spacing_, num_stock_nodes_);
}

// Where today's variance sits along the variance axis.
// Used to read the price off the finished grid.
AxisPosition Grid::variance_position(double v0) const {
    return locate(v0, variance_spacing_, num_variance_nodes_);
}
