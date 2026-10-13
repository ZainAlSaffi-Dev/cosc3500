// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "black_scholes.h"

#include <cmath>

// Chance that a standard normal value is below x.
// Uses N(x) = erfc(-x / sqrt(2)) / 2.
double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// The textbook Black-Scholes price with a continuous dividend yield.
// The tests compare the PDE solver against this with the variance fixed.
double bs_price(bool is_call, double spot, double strike, double rate, double div_yield, double sigma, double maturity_years) {
    double sqrt_maturity = std::sqrt(maturity_years);
    double d1 = (std::log(spot / strike) + (rate - div_yield + 0.5 * sigma * sigma) * maturity_years) / (sigma * sqrt_maturity);
    double d2 = d1 - sigma * sqrt_maturity;

    // Discount the stock leg by the dividends and the strike leg by the interest rate.
    double stock_leg = spot * std::exp(-div_yield * maturity_years);
    double strike_leg = strike * std::exp(-rate * maturity_years);

    if (is_call) {
        return stock_leg * norm_cdf(d1) - strike_leg * norm_cdf(d2);
    } else {
        // The put swaps the two legs.
        return strike_leg * norm_cdf(-d2) - stock_leg * norm_cdf(-d1);
    }
}
