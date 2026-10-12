// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_GRID_H
#define HESTON_GRID_H

#include <vector>

#include "params.h"

// Where a point such as today's stock price falls along one axis.
// cell is the node just below it, weight is how far towards the next node (0 to 1).
struct AxisPosition {
    int cell = 0;
    double weight = 0.0;
};

// The grid of option values, stock price across and variance down.
// It holds two copies of the sheet, one being read and one being written.
class Grid {
public:
    Grid(const GridSpec& spec, const OptionSpec& opt);

    // Simple getters for the grid size and spacing.
    // The const at the end means calling them cannot change the grid.
    int num_stock_nodes() const { return num_stock_nodes_; }
    int num_variance_nodes() const { return num_variance_nodes_; }
    double stock_spacing() const { return stock_spacing_; }
    double variance_spacing() const { return variance_spacing_; }

    // Stock price at column stock_i and variance at row var_j.
    double stock_price(int stock_i) const { return stock_i * stock_spacing_; }
    double variance(int var_j) const { return var_j * variance_spacing_; }

    // Position of node (stock_i, var_j) in the flat buffer.
    // Each row is laid out one after another in memory.
    int index(int stock_i, int var_j) const { return var_j * num_stock_nodes_ + stock_i; }

    // The sheet being read this step, given back read only.
    const std::vector<double>& current() const { return current_; }
    // The sheet being written this step.
    std::vector<double>& next() { return next_; }

    void swap_buffers();
    void init_payoff(const OptionSpec& opt);
    AxisPosition stock_position(double spot) const;
    AxisPosition variance_position(double v0) const;

private:
    int num_stock_nodes_ = 0;
    int num_variance_nodes_ = 0;
    // Gap between stock nodes, called ds in the maths
    double stock_spacing_ = 0.0;
    // Gap between variance nodes, called dv in the maths
    double variance_spacing_ = 0.0;
    // The two sheets
    std::vector<double> current_;
    std::vector<double> next_;
};

#endif  // HESTON_GRID_H
