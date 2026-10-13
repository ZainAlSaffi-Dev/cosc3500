// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "params.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

// Stores one key and value from the file in the matching Config field.
// An unknown key is an error.
static void apply(Config& cfg, const std::string& key, const std::string& value, const std::string& where) {
    // std::stod turns text into a double and std::stoi into an int.
    if (key == "option.strike") {
        cfg.option.strike = std::stod(value);
    } else if (key == "option.maturity_years") {
        cfg.option.maturity_years = std::stod(value);
    } else if (key == "option.type") {
        if (value != "call" && value != "put") {
            throw std::runtime_error(where + ": option.type must be call or put");
        }
        cfg.option.is_call = (value == "call");
    } else if (key == "market.spot") {
        cfg.market.spot = std::stod(value);
    } else if (key == "market.rate") {
        cfg.market.rate = std::stod(value);
    } else if (key == "market.div_yield") {
        cfg.market.div_yield = std::stod(value);
    } else if (key == "heston.v0") {
        cfg.heston.v0 = std::stod(value);
    } else if (key == "heston.kappa") {
        cfg.heston.kappa = std::stod(value);
    } else if (key == "heston.theta") {
        cfg.heston.theta = std::stod(value);
    } else if (key == "heston.xi") {
        cfg.heston.xi = std::stod(value);
    } else if (key == "heston.rho") {
        cfg.heston.rho = std::stod(value);
    } else if (key == "grid.ns") {
        cfg.grid.num_stock_nodes = std::stoi(value);
    } else if (key == "grid.nv") {
        cfg.grid.num_variance_nodes = std::stoi(value);
    } else if (key == "grid.nt") {
        cfg.grid.num_timesteps = std::stoi(value);
    } else if (key == "grid.s_max_mult") {
        cfg.grid.stock_max_multiplier = std::stod(value);
    } else if (key == "grid.v_max") {
        cfg.grid.variance_max = std::stod(value);
    } else {
        // Nothing matched, so the key is a typo.
        throw std::runtime_error(where + ": unknown key '" + key + "'");
    }
}

// Reads a config file of key = value lines and builds a Config from it.
// Blank lines and lines starting with a hash are skipped.
Config load_config(const std::string& path) {
    // The file closes itself at the end of the function, like with open in Python.
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open config: " + path);
    }

    // Start from the defaults in params.h.
    Config cfg;
    std::string line;
    int line_number = 0;
    while (std::getline(file, line)) {
        line_number = line_number + 1;

        // Split the line into its first three words, like line.split() in Python.
        // A missing word stays an empty string.
        std::istringstream words(line);
        std::string key;
        std::string equals;
        std::string value;
        words >> key >> equals >> value;

        // Skip blank lines and comment lines.
        if (key.empty() || key[0] == '#') {
            continue;
        }

        // File name and line number for error messages.
        std::string where = path + ":" + std::to_string(line_number);
        if (equals != "=" || value.empty()) {
            throw std::runtime_error(where + ": expected key = value with spaces around the =");
        }

        // stod and stoi throw invalid_argument when the value is not a number.
        try {
            apply(cfg, key, value, where);
        } catch (const std::invalid_argument&) {
            throw std::runtime_error(where + ": not a number '" + value + "'");
        }
    }
    return cfg;
}
