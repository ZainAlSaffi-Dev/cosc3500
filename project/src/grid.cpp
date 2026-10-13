// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "grid.h"

#include <algorithm>
#include <cmath>
#include <utility>

// Sets up the grid size and spacing and makes both buffers full of zeros.
// Runs once when a Grid is created, like __init__ in Python.
// Grid:: means the function belongs to the Grid class.
Grid::Grid(const GridSpec& spec, const OptionSpec& opt) {
    num_stock_nodes_ = spec.num_stock_nodes;
    num_variance_nodes_ = spec.num_variance_nodes;

    // N nodes have N - 1 gaps between them.
    stock_spacing_ = spec.stock_max_multiplier * opt.strike / (num_stock_nodes_ - 1);
    variance_spacing_ = spec.variance_max / (num_variance_nodes_ - 1);

    // Each buffer is one block of doubles, like np.zeros(n) in Python.
    // The vectors free this memory when the Grid is destroyed.
    int total_cells = num_stock_nodes_ * num_variance_nodes_;
    current_.assign(total_cells, 0.0);
    next_.assign(total_cells, 0.0);
}

// Swaps which buffer is read and which is written.
// Only the internal pointers are swapped. No numbers are copied.
void Grid::swap_buffers() {
    std::swap(current_, next_);
}

// Fills the current buffer with what the option pays at expiry.
// The payoff depends only on the stock price, so every row is the same.
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

// A point closer than this many cells to a node counts as on the node.
static const double node_snap_tolerance = 1e-9;

// Finds which cell a value falls in along one axis, and how far across it.
// The cell is kept one node away from each edge.
static AxisPosition locate(double coord, double spacing, int num_nodes) {
    // Position in cells. 2.5 is halfway from node 2 to node 3.
    double position = coord / spacing;
    // Keep the position on the grid.
    if (position < 0.0) {
        position = 0.0;
    }
    if (position > num_nodes - 1) {
        position = num_nodes - 1;
    }

    // The node at or just below the position.
    int lower_node = (int)std::floor(position);
    // A position within rounding error of a node counts as on the node.
    double nearest_node = std::round(position);
    bool on_node = std::fabs(position - nearest_node) < node_snap_tolerance;
    if (on_node) {
        lower_node = (int)nearest_node;
    }

    // The Greeks read one node either side of the cell and the cell above it.
    // The cell must be between 1 and num_nodes - 3.
    int highest_cell = num_nodes - 3;
    if (highest_cell < 1) {
        highest_cell = 1;
    }
    AxisPosition out;
    out.cell = lower_node;
    if (out.cell < 1) {
        out.cell = 1;
    }
    if (out.cell > highest_cell) {
        out.cell = highest_cell;
    }

    // The weight is 0 on a node, otherwise the fraction past the cell.
    if (on_node && out.cell == lower_node) {
        out.weight = 0.0;
    } else {
        out.weight = position - out.cell;
        // Keep the weight between 0 and 1 when the cell was moved in from an edge.
        if (out.weight < 0.0) {
            out.weight = 0.0;
        }
        if (out.weight > 1.0) {
            out.weight = 1.0;
        }
    }
    return out;
}

// Finds where today's stock price sits along the stock axis.
// Used to read the price off the finished grid.
AxisPosition Grid::stock_position(double spot) const {
    return locate(spot, stock_spacing_, num_stock_nodes_);
}

// Finds where today's variance sits along the variance axis.
// Used to read the price off the finished grid.
AxisPosition Grid::variance_position(double v0) const {
    return locate(v0, variance_spacing_, num_variance_nodes_);
}
