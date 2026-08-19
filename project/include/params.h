#ifndef HESTON_PARAMS_H
#define HESTON_PARAMS_H

#include <string>

// One European option contract.
struct OptionSpec {
    double strike = 5250.0;
    double maturity_years = 0.25;
    bool is_call = true;
};

// Market state today.
struct MarketSpec {
    double spot = 5200.0;
    //interest rate, continuously compounded
    double rate = 0.045;
    //dividend yield, continuously compounded
    double div_yield = 0.013;
};

// Heston model parameters (spec §2).
struct HestonParams {
    //initial variance
    double v0 = 0.04;
    //mean-reversion speed
    double kappa = 1.5;
    //long-run variance
    double theta = 0.04;
    //vol of vol
    double xi = 0.35;
    //price/vol correlation
    double rho = -0.70;
};

// How the grid is discretised. The stock axis runs from zero up to
// stock_max_multiplier times the strike and the variance axis from zero up to
// variance_max, both with uniform spacing. Config files and the command line
// keep using the short keys ns, nv, nt, s_max_mult and v_max, which
// load_config maps onto these longer field names.
struct GridSpec {
    //stock nodes (contiguous dimension)    
    int num_stock_nodes = 2048;  // cfg/CLI key: ns
    //variance nodes
    int num_variance_nodes = 512;  // cfg/CLI key: nv
    //timesteps
    int num_timesteps = 2000;  // cfg/CLI key: nt
    // The top of the stock axis, given as a multiple of the strike. It is far
    // enough out that the option's value there is effectively known, because
    // the option is deep in or out of the money by that point.
    double stock_max_multiplier = 4.0;  // cfg/CLI key: s_max_mult
    // The top of the variance axis. Variance is volatility squared, so a
    // value of 1.0 means 100% annual volatility, which is well past anything
    // a real market does.
    double variance_max = 1.0;  // cfg/CLI key: v_max
};

// Everything one run needs. The command line overrides whatever the config
// file set.
struct Config {
    OptionSpec option;
    MarketSpec market;
    HestonParams heston;
    GridSpec grid;

    // I/O behaviour (CLI-only, not in .cfg)
    // 0 = no snapshots
    int dump_every = 0;
    //dump directory
    std::string dump_dir = "results";
    //benchmark repetitions
    int bench_reps = 0;
    //solver (baseline which is un optimised and opt which is optimised)
    std::string solver = "baseline"; // baseline | opt
    // Which rung of the optimisation ladder the opt solver runs, settable
    // only from the command line. Level 0 is the same algorithm as the
    // baseline and level 6 has every technique applied, which is why 6 is the
    // default now that the whole ladder is written. Levels 7 and 8 select the
    // negative controls and are never a default.
    int opt_level = 6;
};

// Reads a file of "key = value" lines. It throws std::runtime_error if a key
// is not recognised or a value does not parse.
Config load_config(const std::string& path);

#endif  // HESTON_PARAMS_H
