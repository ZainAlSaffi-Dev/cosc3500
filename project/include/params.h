// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_PARAMS_H
#define HESTON_PARAMS_H

#include <string>

// The option contract being priced, with its strike, time to expiry and type.
// The numbers below are the defaults.
struct OptionSpec {
    double strike = 5250.0;
    double maturity_years = 0.25;
    bool is_call = true;
};

// The market today, the stock price, interest rate and dividend yield.
// Used to read off the price and to discount.
struct MarketSpec {
    double spot = 5200.0;
    // Risk free interest rate, compounded continuously.
    double rate = 0.045;
    // Dividend yield, compounded continuously.
    double div_yield = 0.013;
};

// The five numbers that define the Heston model.
// They control how the variance moves over time.
struct HestonParams {
    // Variance today.
    double v0 = 0.04;
    // How fast the variance pulls back towards theta.
    double kappa = 1.5;
    // Long run variance.
    double theta = 0.04;
    // Volatility of the variance.
    double xi = 0.35;
    // Correlation between stock moves and variance moves.
    double rho = -0.70;
};

// How big the grid is and how far each axis reaches.
// Config files use the keys ns, nv, nt, s_max_mult and v_max.
struct GridSpec {
    // Number of stock nodes. Neighbouring stock nodes sit next to each other in memory.
    int num_stock_nodes = 2048;
    // Number of variance nodes, one row of the grid each.
    int num_variance_nodes = 512;
    // Number of timesteps from expiry back to today.
    int num_timesteps = 2000;
    // Top of the stock axis as a multiple of the strike.
    double stock_max_multiplier = 4.0;
    // Top of the variance axis. 1.0 means 100 percent volatility.
    double variance_max = 1.0;
};

// Every setting one run needs.
// The config file fills it first, then command line flags overwrite parts of it.
struct Config {
    OptionSpec option;
    MarketSpec market;
    HestonParams heston;
    GridSpec grid;

    // Write a snapshot every this many steps. 0 means never.
    int dump_every = 0;
    // Folder for the snapshots.
    std::string dump_dir = "results";
    // How many times to repeat the solve for benchmarking.
    int bench_reps = 0;
    // Which solver runs, for example baseline, serial, simd, omp or mpi.
    std::string solver = "baseline";
};

// Reads a config file of key = value lines and builds a Config from it.
Config load_config(const std::string& path);

#endif  // HESTON_PARAMS_H
