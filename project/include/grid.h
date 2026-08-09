#ifndef HESTON_GRID_H
#define HESTON_GRID_H

#include <vector>

#include "params.h"

// Two flat buffers over the (variance x stock) sheet.
// Layout: index = j*ns + i, i = stock index (contiguous), j = variance index.
// cur() is the finished sheet (read-only), next() is being written,
// swap_buffers() after every step.
class Grid {
public:
    // Constructor allocates both buffers and computes spacings.
    Grid(const GridSpec& spec, const OptionSpec& opt);

    // Getters. 'const' after () = "does not mutate self"
    // analogue; the compiler enforces it.
    int ns() const { return ns_; }
    int nv() const { return nv_; }
    double ds() const { return ds_; }
    double dv() const { return dv_; }
    double s(int i) const { return i * ds_; }
    double v(int j) const { return j * dv_; }
    int idx(int i, int j) const { return j * ns_ + i; }

    // cur() returns a const reference: caller can read, cannot write.
    const std::vector<double>& cur() const { return cur_; }
    std::vector<double>& next() { return next_; }

    // std::swap on vectors exchanges internal pointers — O(1), no copy.
    void swap_buffers();

    // Fill cur() with the expiry payoff.
    void init_payoff(const OptionSpec& opt);

    // Nearest grid indices to (spot, v0) for the price/Greeks readout.
    int nearest_i(double spot) const;
    int nearest_j(double v0) const;

private:
    int ns_ = 0, nv_ = 0;
    double ds_ = 0.0, dv_ = 0.0;
    std::vector<double> cur_;
    std::vector<double> next_;
};

#endif  // HESTON_GRID_H
