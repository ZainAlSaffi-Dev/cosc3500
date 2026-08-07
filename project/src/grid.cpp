#include "grid.h"

Grid::Grid(const GridSpec& spec, const OptionSpec& opt) {
    // TODO(P1): ns_/nv_ from spec; ds_ = s_max_mult*strike/(ns-1);
    // dv_ = v_max/(nv-1); cur_.assign(ns_*nv_, 0.0); same for next_.
    (void)spec;
    (void)opt;
}

void Grid::swap_buffers() {
    // TODO(P1): std::swap(cur_, next_);  // pointer exchange inside, O(1)
}

void Grid::init_payoff(const OptionSpec& opt) {
    // TODO(P2): every j-row identical at expiry:
    // call max(s(i)-K, 0), put max(K-s(i), 0). Write into cur_.
    (void)opt;
}

int Grid::nearest_i(double spot) const {
    // TODO(P2): round(spot / ds_), clamped to [0, ns_-1].
    (void)spot;
    return 0;
}

int Grid::nearest_j(double v0) const {
    // TODO(P2): round(v0 / dv_), clamped to [0, nv_-1].
    (void)v0;
    return 0;
}
