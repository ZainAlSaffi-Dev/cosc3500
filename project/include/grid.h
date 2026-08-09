#ifndef HESTON_GRID_H
#define HESTON_GRID_H

#include <vector>

#include "params.h"

// Two flat buffers over the (variance x stock) sheet.
// Layout: index = var_j * num_stock_nodes + stock_i — stock is the contiguous
// direction. (Docs shorthand for the same thing: idx = j*ns + i.)
// current() is the finished sheet (read-only), next() is being written,
// swap_buffers() after every step.
class Grid {
public:
    // Constructor allocates both buffers and computes spacings.
    Grid(const GridSpec& spec, const OptionSpec& opt);

    // Getters. 'const' after the () is a compiler-enforced promise that the
    // method does not modify the object (Python has no equivalent).
    int num_stock_nodes() const { return num_stock_nodes_; }
    int num_variance_nodes() const { return num_variance_nodes_; }
    double stock_spacing() const { return stock_spacing_; }
    double variance_spacing() const { return variance_spacing_; }

    // Coordinate values at a grid node.
    double stock_price(int stock_i) const { return stock_i * stock_spacing_; }
    double variance(int var_j) const { return var_j * variance_spacing_; }

    // Flat position of node (stock_i, var_j) in the 1-D buffers.
    int index(int stock_i, int var_j) const {
        return var_j * num_stock_nodes_ + stock_i;
    }

    // current() returns a const reference: caller can read, cannot write.
    const std::vector<double>& current() const { return current_; }
    std::vector<double>& next() { return next_; }

    // std::swap on vectors exchanges internal pointers — O(1), no copy.
    void swap_buffers();

    // Fill current() with the expiry payoff.
    void init_payoff(const OptionSpec& opt);

    // Nearest grid indices to today's (spot, v0) for the price/Greeks readout.
    int nearest_stock_index(double spot) const;
    int nearest_variance_index(double v0) const;

private:
    int num_stock_nodes_ = 0;
    int num_variance_nodes_ = 0;
    double stock_spacing_ = 0.0;     // spacing between stock nodes (docs: ds)
    double variance_spacing_ = 0.0;  // spacing between variance nodes (docs: dv)
    std::vector<double> current_;
    std::vector<double> next_;
};

#endif  // HESTON_GRID_H
