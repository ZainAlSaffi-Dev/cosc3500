#ifndef HESTON_BLACK_SCHOLES_H
#define HESTON_BLACK_SCHOLES_H

// Closed-form Black-Scholes for the validation gate only.
// with xi=0 and v0=theta the Heston solver must converge to this.

double norm_cdf(double x);

// sigma = volatility (sqrt of variance), q = continuous dividend yield.
double bs_price(bool is_call, double spot, double strike, double rate,
                double div_yield, double sigma, double maturity_years);

#endif  // HESTON_BLACK_SCHOLES_H
