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

// Discretisation. S in [0, stock_max_multiplier*strike], v in
// [0, variance_max], uniform spacing. Config files and CLI keep the short
// keys (ns, nv, nt, s_max_mult, v_max); load_config maps them onto these
// descriptive fields.
struct GridSpec {
    //stock nodes (contiguous dimension)    
    int num_stock_nodes = 2048;  // cfg/CLI key: ns
    //variance nodes
    int num_variance_nodes = 512;  // cfg/CLI key: nv
    //timesteps
    int num_timesteps = 2000;  // cfg/CLI key: nt
    // Top of the stock axis, as a multiple of strike: S spans
    // [0, stock_max_multiplier*strike]. Far enough out that the option's
    // value there is effectively known (deep in/out of the money).
    double stock_max_multiplier = 4.0;  // cfg/CLI key: s_max_mult
    // Top of the variance axis. v is variance = volatility^2, so 1.0 means
    // 100% annual volatility — beyond any realistic market scenario.
    double variance_max = 1.0;  // cfg/CLI key: v_max
};

// Everything a run needs. CLI overrides applied on top of file values.
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
    // Optimisation-ladder level for the opt solver (PLAN §4b), CLI-only.
    // 0 = same algorithm as baseline. Stays 0 until the full ladder lands;
    // then PLAN §3's default of 6 (= all techniques) takes over.
    int opt_level = 0;
};

// Parse key = value file. Throws std::runtime_error on unknown key / bad value.
Config load_config(const std::string& path);

#endif  // HESTON_PARAMS_H
