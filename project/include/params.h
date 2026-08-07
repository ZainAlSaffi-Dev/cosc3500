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
    double rate = 0.045;       // r, continuously compounded
    double div_yield = 0.013;  // q
};

// Heston model parameters (spec §2).
struct HestonParams {
    double v0 = 0.04;     // initial variance
    double kappa = 1.5;   // mean-reversion speed
    double theta = 0.04;  // long-run variance
    double xi = 0.35;     // vol of vol
    double rho = -0.70;   // price/vol correlation
};

// Discretisation. S in [0, s_max_mult*strike], v in [0, v_max], uniform.
struct GridSpec {
    int ns = 2048;          // stock nodes (contiguous dimension)
    int nv = 512;           // variance nodes
    int nt = 2000;          // timesteps
    double s_max_mult = 4.0;
    double v_max = 1.0;
};

// Everything a run needs. CLI overrides applied on top of file values.
struct Config {
    OptionSpec option;
    MarketSpec market;
    HestonParams heston;
    GridSpec grid;

    // I/O behaviour (CLI-only, not in .cfg)
    int dump_every = 0;              // 0 = no snapshots
    std::string dump_dir = "results";
    int bench_reps = 0;              // 0 = single solve
    std::string solver = "baseline"; // baseline | opt
};

// Parse key = value file. Throws std::runtime_error on unknown key / bad value.
Config load_config(const std::string& path);

#endif  // HESTON_PARAMS_H
