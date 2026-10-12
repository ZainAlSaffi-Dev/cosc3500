// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_BLACK_SCHOLES_H
#define HESTON_BLACK_SCHOLES_H

// The standard normal cumulative distribution function.
double norm_cdf(double x);

// The exact Black-Scholes price, used only by the tests.
double bs_price(bool is_call, double spot, double strike, double rate, double div_yield, double sigma, double maturity_years);

#endif  // HESTON_BLACK_SCHOLES_H
