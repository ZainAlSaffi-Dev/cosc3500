#include "params.h"

#include <fstream>  // std::ifstream: opens a file for reading and closes it
                    // itself at scope exit (RAII). Similar to Python's
                    // "with open(path) as f".
#include <stdexcept>
#include <string>

// File-local helpers. 'static' here means visible only inside this .cpp file,
// similar to a module-private _function in Python.

static std::string trim(const std::string& text) {
    // Strip spaces/tabs/CR from both ends, similar to text.strip().
    const std::size_t first = text.find_first_not_of(" \t\r");
    if (first == std::string::npos) return "";
    const std::size_t last = text.find_last_not_of(" \t\r");
    return text.substr(first, last - first + 1);
}

static double parse_double(const std::string& value, const std::string& where) {
    std::size_t used = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &used);  // similar to float(value)
    } catch (const std::exception&) {
        throw std::runtime_error(where + ": not a number: '" + value + "'");
    }
    // 'used' is how many characters stod consumed, so this rejects "1.5abc".
    if (used != value.size())
        throw std::runtime_error(where + ": trailing text in '" + value + "'");
    return parsed;
}

static int parse_int(const std::string& value, const std::string& where) {
    std::size_t used = 0;
    int parsed = 0;
    try {
        parsed = std::stoi(value, &used);  // similar to int(value)
    } catch (const std::exception&) {
        throw std::runtime_error(where + ": not an integer: '" + value + "'");
    }
    if (used != value.size())
        throw std::runtime_error(where + ": trailing text in '" + value + "'");
    return parsed;
}

// Each key and value pair goes into its own Config field. An unrecognised key
// is an error, so a typo in a .cfg file cannot quietly run with the defaults.
static void apply(Config& cfg, const std::string& key, const std::string& value,
                  const std::string& where) {
    if (key == "option.strike") { cfg.option.strike = parse_double(value, where); return; }
    if (key == "option.maturity_years") { cfg.option.maturity_years = parse_double(value, where); return; }
    if (key == "option.type") {
        if (value != "call" && value != "put")
            throw std::runtime_error(where + ": option.type must be call or put");
        cfg.option.is_call = (value == "call");
        return;
    }
    if (key == "market.spot") { cfg.market.spot = parse_double(value, where); return; }
    if (key == "market.rate") { cfg.market.rate = parse_double(value, where); return; }
    if (key == "market.div_yield") { cfg.market.div_yield = parse_double(value, where); return; }
    if (key == "heston.v0") { cfg.heston.v0 = parse_double(value, where); return; }
    if (key == "heston.kappa") { cfg.heston.kappa = parse_double(value, where); return; }
    if (key == "heston.theta") { cfg.heston.theta = parse_double(value, where); return; }
    if (key == "heston.xi") { cfg.heston.xi = parse_double(value, where); return; }
    if (key == "heston.rho") { cfg.heston.rho = parse_double(value, where); return; }
    // The short config keys map onto the longer GridSpec field names.
    if (key == "grid.ns") { cfg.grid.num_stock_nodes = parse_int(value, where); return; }
    if (key == "grid.nv") { cfg.grid.num_variance_nodes = parse_int(value, where); return; }
    if (key == "grid.nt") { cfg.grid.num_timesteps = parse_int(value, where); return; }
    if (key == "grid.s_max_mult") { cfg.grid.stock_max_multiplier = parse_double(value, where); return; }
    if (key == "grid.v_max") { cfg.grid.variance_max = parse_double(value, where); return; }
    throw std::runtime_error(where + ": unknown key '" + key + "'");
}

Config load_config(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open config: " + path);
    Config cfg;
    std::string line;
    int line_no = 0;
    while (std::getline(file, line)) {  // similar to: for line in file
        ++line_no;
        // Everything after '#' is a comment.
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = trim(line);
        if (line.empty()) continue;
        const std::string where = path + ":" + std::to_string(line_no);
        const std::size_t eq = line.find('=');
        if (eq == std::string::npos)
            throw std::runtime_error(where + ": expected 'key = value'");
        apply(cfg, trim(line.substr(0, eq)), trim(line.substr(eq + 1)), where);
    }
    return cfg;
}
