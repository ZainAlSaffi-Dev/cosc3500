#include "grid.h"

#include <algorithm>
#include <cmath>
#include <utility>

// Everything after the colon on the next line is a member-initialiser list. It
// fills in the fields while the object is being built, before the constructor
// body starts running. This is similar to assigning self.num_stock_nodes at
// the top of __init__ in Python.
Grid::Grid(const GridSpec& spec, const OptionSpec& opt)
    : num_stock_nodes_(spec.num_stock_nodes),
      num_variance_nodes_(spec.num_variance_nodes) {
    // N nodes spanning a range leave N-1 gaps between them, hence the -1.
    stock_spacing_ =
        spec.stock_max_multiplier * opt.strike / (num_stock_nodes_ - 1);
    variance_spacing_ = spec.variance_max / (num_variance_nodes_ - 1);
    // Each buffer is one solid block of doubles rather than a list of lists,
    // much like np.zeros(n) in Python. The vectors own that memory and free
    // it when the Grid is destroyed, so nothing here is freed by hand.
    const std::size_t total_cells =
        static_cast<std::size_t>(num_stock_nodes_) * num_variance_nodes_;
    current_.assign(total_cells, 0.0);
    next_.assign(total_cells, 0.0);
}

void Grid::swap_buffers() {
    // two buffers, swap pointers: avoids copying the whole grid each timestep.
    // std::swap on two vectors only exchanges the pointers they hold inside,
    // so this costs the same no matter how large the grid is and no element
    // is ever copied.
    std::swap(current_, next_);
}

void Grid::init_payoff(const OptionSpec& opt) {
    // At expiry the option is worth exactly its payoff, because there is no
    // time left for anything else to happen. The payoff depends only on the
    // stock price, so every variance row ends up holding the same numbers.
   

    // var_j outer, stock_i inner (row-major: the memory walk is contiguous)
    // stock_price(stock_i) is the price at that column; opt.strike is K
    // call: max(price - strike, 0.0)
    // put: max(strike - price, 0.0)
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

// An unnamed namespace makes everything inside it visible only to this file,
// which is how C++ marks a private helper that is not part of the public API.
// It plays the same role as a leading underscore on a module-level function in
// Python.
namespace {

// A quote point this close to a node, measured in cells, counts as sitting on
// it. Snapping means an aligned grid produces a weight of exactly 0.0 rather
// than a leftover 1e-16 of rounding noise, which is what lets the blend
// reproduce the nodal value bit for bit.
constexpr double kNodeSnapCells = 1e-9;

AxisPosition locate(double coord, double spacing, int num_nodes) {
    // Clamping happens before any index arithmetic, so a quote point outside
    // the grid reads the nearest edge cell instead of running off the end of
    // the array.
    const double position = std::clamp(coord / spacing, 0.0,
                                       static_cast<double>(num_nodes - 1));
    const double nearest_node = std::round(position);
    const bool on_node = std::fabs(position - nearest_node) < kNodeSnapCells;
    const double lower_node = on_node ? nearest_node : std::floor(position);

    // The cell stays one node in from each edge, because the blend also reads
    // cell+1 and the gamma stencil reaches cell-1 and cell+2. The std::max
    // only matters on a grid with fewer than four nodes on an axis.
    const int highest_cell = std::max(1, num_nodes - 3);
    AxisPosition out;
    out.cell = std::clamp(static_cast<int>(lower_node), 1, highest_cell);
    // The weight is exactly 0.0 when the snap fired and the clamp then left
    // the cell alone. Otherwise it is the fraction of a cell the point sits
    // past the lower node, clamped so that the blend never extrapolates.
    out.weight = (on_node && out.cell == static_cast<int>(lower_node))
                     ? 0.0
                     : std::clamp(position - out.cell, 0.0, 1.0);
    return out;
}

}  // namespace

AxisPosition Grid::stock_position(double spot) const {
    return locate(spot, stock_spacing_, num_stock_nodes_);
}

AxisPosition Grid::variance_position(double v0) const {
    return locate(v0, variance_spacing_, num_variance_nodes_);
}
