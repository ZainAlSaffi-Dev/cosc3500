#ifndef HESTON_GRID_H
#define HESTON_GRID_H

#include <vector>

#include "params.h"

// Where a quote point sits along one axis. The spot and variance we are asked
// about usually fall between two nodes, so cell names the node just below and
// weight says how far across that cell the point lies, from 0 on the lower
// node to 1 on the upper one.
struct AxisPosition {
    int cell = 0;
    double weight = 0.0;
};

// The (variance x stock) sheet, held as two flat buffers. A node at column
// stock_i and row var_j lives at index var_j * num_stock_nodes + stock_i, so
// walking along the stock direction walks straight through memory. The solver
// reads the finished sheet from current(), writes the new one into next(),
// then calls swap_buffers() to exchange them.
class Grid {
public:
    Grid(const GridSpec& spec, const OptionSpec& opt);

    // The const after the brackets is a promise to the compiler that calling
    // this method cannot modify the object, and the compiler enforces it.
    // Python has no equivalent.
    int num_stock_nodes() const { return num_stock_nodes_; }
    int num_variance_nodes() const { return num_variance_nodes_; }
    double stock_spacing() const { return stock_spacing_; }
    double variance_spacing() const { return variance_spacing_; }

    double stock_price(int stock_i) const { return stock_i * stock_spacing_; }
    double variance(int var_j) const { return var_j * variance_spacing_; }

    int index(int stock_i, int var_j) const {
        return var_j * num_stock_nodes_ + stock_i;
    }

    // current() hands back a const reference, which lets the caller read the
    // buffer without being able to write to it.
    const std::vector<double>& current() const { return current_; }
    std::vector<double>& next() { return next_; }

    void swap_buffers();

    void init_payoff(const OptionSpec& opt);

    // These locate the quote point for the price readout. The cell index is
    // kept one node in from each edge so that every cell the blend and the
    // Greek stencils reach still exists. On a node-aligned grid the weight
    // comes out as exactly 0.0, so the blend reproduces the nodal value
    // bit for bit.
    AxisPosition stock_position(double spot) const;
    AxisPosition variance_position(double v0) const;

private:
    int num_stock_nodes_ = 0;
    int num_variance_nodes_ = 0;
    double stock_spacing_ = 0.0;     // the docs call this ds
    double variance_spacing_ = 0.0;  // the docs call this dv
    std::vector<double> current_;
    std::vector<double> next_;
};

#endif  // HESTON_GRID_H
