#include "black_scholes.h"

double norm_cdf(double x) {
    // TODO(P3): 0.5 * erfc(-x / sqrt(2)) — <cmath> erfc, no approximation needed.
    (void)x;
    return 0.0;
}

double bs_price(bool is_call, double spot, double strike, double rate,
                double div_yield, double sigma, double maturity_years) {
    // TODO(P3): textbook formula with continuous dividend yield q:
    // d1 = (ln(S/K) + (r - q + sigma^2/2)T) / (sigma sqrt(T)), d2 = d1 - sigma sqrt(T)
    // call = S e^{-qT} N(d1) - K e^{-rT} N(d2); put via parity or direct.
    (void)is_call;
    (void)spot;
    (void)strike;
    (void)rate;
    (void)div_yield;
    (void)sigma;
    (void)maturity_years;
    return 0.0;
}
