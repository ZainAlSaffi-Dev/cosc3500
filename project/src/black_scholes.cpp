// Closed-form Black-Scholes, the known-good answer for the P3 validation
// gate. The solvers never call it, only the tests compare against it.

#include "black_scholes.h"

#include <cmath>

double norm_cdf(double x) {
    // erfc from <cmath> is the complementary error function, and the standard
    // normal CDF is an exact rescaling of it with no approximation involved,
    // since N(x) is erfc(-x / sqrt(2)) / 2.
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double bs_price(bool is_call, double spot, double strike, double rate,
                double div_yield, double sigma, double maturity_years) {
    // This is the textbook formula with a continuous dividend yield q, where
    //   d1 = (ln(S/K) + (r - q + sigma^2/2)*T) / (sigma*sqrt(T))
    //   d2 = d1 - sigma*sqrt(T)
    const double sqrt_maturity = std::sqrt(maturity_years);
    const double d1 =
        (std::log(spot / strike) +
         (rate - div_yield + 0.5 * sigma * sigma) * maturity_years) /
        (sigma * sqrt_maturity);
    const double d2 = d1 - sigma * sqrt_maturity;

    // Both legs are discounted back to today, the stock leg by the dividend
    // yield and the strike leg by the risk-free rate.
    const double stock_leg = spot * std::exp(-div_yield * maturity_years);
    const double strike_leg = strike * std::exp(-rate * maturity_years);

    if (is_call) return stock_leg * norm_cdf(d1) - strike_leg * norm_cdf(d2);
    // The put is the same formula mirrored, using N(-d) and swapping the legs.
    return strike_leg * norm_cdf(-d2) - stock_leg * norm_cdf(-d1);
}
