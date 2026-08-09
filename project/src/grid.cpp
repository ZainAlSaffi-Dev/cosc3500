#include "grid.h"

#include <algorithm>  // std::clamp (and std::max for the payoff)
#include <cmath>      // std::lround
#include <utility>    // std::swap

// The num_stock_nodes part is a member-initialiser list: it fills in
// the fields while the object is being created, before the constructor body
// runs — similar to assigning self.num_stock_nodes at the top of __init__ in Python.
Grid::Grid(const GridSpec& spec, const OptionSpec& opt)
    : num_stock_nodes_(spec.num_stock_nodes),
      num_variance_nodes_(spec.num_variance_nodes) {
    // N nodes spanning a range leave N-1 gaps between them — hence the -1.
    stock_spacing_ =
        spec.stock_max_multiplier * opt.strike / (num_stock_nodes_ - 1);
    variance_spacing_ = spec.variance_max / (num_variance_nodes_ - 1);
    // One flat contiguous allocation per buffer (row-major, stock contiguous):
    // a single solid block of doubles, similar to np.zeros(n) in Python —
    // never a list-of-lists. The vectors have ownership of this memory, and then grid frees it (RAII).
    const std::size_t total_cells =
        static_cast<std::size_t>(num_stock_nodes_) * num_variance_nodes_;
    current_.assign(total_cells, 0.0);
    next_.assign(total_cells, 0.0);
}

void Grid::swap_buffers() {
    // two buffers, swap pointers: avoids copying the whole grid each timestep.
    // std::swap on vectors exchanges internal pointers — O(1), no element copy.
    std::swap(current_, next_);
}

void Grid::init_payoff(const OptionSpec& opt) {
    // At expiry the option is worth exactly its payoff — no time value left —
    // and the payoff depends on the stock price only, so every variance row
    // is identical.
    //
    //
    // var_j outer, stock_i inner (row-major: the memory walk is contiguous)
    // stock_price(stock_i) is the price at that column; opt.strike is K
    // call: max(price - strike, 0.0)    
    //put: max(strike - price, 0.0)
    // opt.is_call picks which payoff to use
    // write into current_[index(stock_i, var_j)]  (the read buffer)
    for (int var_j = 0; var_j < num_variance_nodes_; ++var_j) {
        for (int stock_i = 0; stock_i < num_stock_nodes_; ++stock_i) {
            // using inline to save time
            double price = stock_price(stock_i);
            // check if std max is fast
            double call_payoff = std::max(price - opt.strike, 0.0);
            double put_payoff = std::max(opt.strike - price, 0.0);
            double payoff = opt.is_call ? call_payoff : put_payoff;
            current_[index(stock_i, var_j)] = payoff;
        }
    }
    
}

int Grid::nearest_stock_index(double spot) const {
    // Nearest column to today's spot: round to grid units, then clamp the
    // index so it stays inside [0, num_stock_nodes-1] — similar to
    // min(max(round(spot/step), 0), n-1) in Python.
    const int stock_i = static_cast<int>(std::lround(spot / stock_spacing_));
    return std::clamp(stock_i, 0, num_stock_nodes_ - 1);
}

int Grid::nearest_variance_index(double v0) const {
    // Nearest row to today's variance v0, same round-and-clamp.
    const int var_j = static_cast<int>(std::lround(v0 / variance_spacing_));
    return std::clamp(var_j, 0, num_variance_nodes_ - 1);
}
