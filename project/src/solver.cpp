// The parts of the solver that both BaselineSolver and OptSolver share, which
// means reading the price and the Greeks off the finished sheet, and the
// factory that picks a solver by name.

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <vector>

#include "io.h"
#include "solver.h"

namespace {

// Grid::stock_position already snaps anything within 1e-9 of a node to a
// weight of exactly 0.0, so a weight above this really did come from a quote
// point sitting between nodes.
constexpr double kOffNodeWarnCells = 1e-9;

// Blends the four cells surrounding the quote point. The price and every arm
// of the Greek stencils go through this same function, so they all come off
// one interpolated surface, which is what keeps the reported delta the actual
// derivative of the reported price.
double blend_cell(const Grid& g, const std::vector<double>& V, int base_i,
                  int base_j, double stock_weight, double variance_weight) {
    const double stock_lower = 1.0 - stock_weight;
    const double variance_lower = 1.0 - variance_weight;
    // This blends along the stock direction twice and then blends those two
    // results along the variance direction. When both weights are exactly 0
    // every term is either 1.0 times a nodal value or 0.0, so the nodal value
    // comes back unchanged, which is the reason the grids are node-aligned.
    return variance_lower *
               (stock_lower * V[g.index(base_i, base_j)] +
                stock_weight * V[g.index(base_i + 1, base_j)]) +
           variance_weight *
               (stock_lower * V[g.index(base_i, base_j + 1)] +
                stock_weight * V[g.index(base_i + 1, base_j + 1)]);
}

}  // namespace

void Solver::extract_result(const Grid& g, const Config& cfg,
                            SolveResult& out) const {
    const std::vector<double>& V = g.current();
    // Today's market point normally falls between grid nodes. Reading the
    // nearest cell instead would cost delta times the offset, which measured
    // as a 4.67% error on the unaligned smoke grid, so the four surrounding
    // cells get blended instead.
    const AxisPosition spot_at = g.stock_position(cfg.market.spot);
    const AxisPosition v0_at = g.variance_position(cfg.heston.v0);

    // Every config is meant to be node-aligned, so it is worth saying out loud
    // when one is not. This goes to stderr so that stdout stays a clean CSV
    // line for the benchmark scripts.
    if (spot_at.weight > kOffNodeWarnCells ||
        v0_at.weight > kOffNodeWarnCells) {
        std::fprintf(stderr,
                     "[%s] off-node quote, spot sits %.4f cells past node %d "
                     "and v0 sits %.4f cells past node %d, so price and Greeks "
                     "are bilinearly interpolated\n",
                     solver_label(cfg).c_str(), spot_at.weight, spot_at.cell,
                     v0_at.weight, v0_at.cell);
    }

    out.price = blend_cell(g, V, spot_at.cell, v0_at.cell, spot_at.weight,
                           v0_at.weight);

    // The same blend one cell to each side, so that the stencils below
    // differentiate the interpolated surface rather than the raw nodes.
    const double east = blend_cell(g, V, spot_at.cell + 1, v0_at.cell,
                                   spot_at.weight, v0_at.weight);
    const double west = blend_cell(g, V, spot_at.cell - 1, v0_at.cell,
                                   spot_at.weight, v0_at.weight);
    const double north = blend_cell(g, V, spot_at.cell, v0_at.cell + 1,
                                    spot_at.weight, v0_at.weight);
    const double south = blend_cell(g, V, spot_at.cell, v0_at.cell - 1,
                                    spot_at.weight, v0_at.weight);

    out.delta = (east - west) / (2.0 * g.stock_spacing());
    out.gamma = (east - 2.0 * out.price + west) /
                (g.stock_spacing() * g.stock_spacing());
    // This vega is measured per unit of variance. The market quotes vega per
    // unit of volatility instead, which is this number times 2*sqrt(v0) by the
    // chain rule, and the write-up says so.
    out.vega = (north - south) / (2.0 * g.variance_spacing());
}

std::unique_ptr<Solver> make_solver(const std::string& name) {
    // std::make_unique performs the allocation so that ownership is never left
    // loose. The unique_ptr it returns deletes the solver as soon as it goes
    // out of scope.
    if (name == "baseline") return std::make_unique<BaselineSolver>();
    if (name == "opt") return std::make_unique<OptSolver>();
    throw std::runtime_error("unknown solver: " + name);
}
