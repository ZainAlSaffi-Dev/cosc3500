// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#include "params.h"

#include <fstream>
#include <stdexcept>
#include <string>

// The word static on these helpers means only this file can use them.

// Removes spaces, tabs and carriage returns from both ends of a string.
// Works like str.strip() in Python.
static std::string trim(const std::string& text) {
    std::size_t first = text.find_first_not_of(" \t\r");
    if (first == std::string::npos) {
        return "";
    }
    std::size_t last = text.find_last_not_of(" \t\r");
    return text.substr(first, last - first + 1);
}

// Turns text like "0.045" into a double.
// Throws an error if the text is not a clean number.
static double parse_double(const std::string& value, const std::string& where) {
    // used counts how many characters stod managed to read
    std::size_t used = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error(where + ": not a number: '" + value + "'");
    }
    // Leftover characters mean something like 1.5abc
    if (used != value.size()) {
        throw std::runtime_error(where + ": trailing text in '" + value + "'");
    }
    return parsed;
}

// Turns text like "2048" into an int.
// Throws an error if the text is not a clean whole number.
static int parse_int(const std::string& value, const std::string& where) {
    std::size_t used = 0;
    int parsed = 0;
    try {
        parsed = std::stoi(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error(where + ": not an integer: '" + value + "'");
    }
    if (used != value.size()) {
        throw std::runtime_error(where + ": trailing text in '" + value + "'");
    }
    return parsed;
}

// Stores one key and value from the file in the matching Config field.
// An unknown key is an error so a typo cannot slip through unnoticed.
static void apply(Config& cfg, const std::string& key, const std::string& value, const std::string& where) {
    if (key == "option.strike") {
        cfg.option.strike = parse_double(value, where);
    } else if (key == "option.maturity_years") {
        cfg.option.maturity_years = parse_double(value, where);
    } else if (key == "option.type") {
        if (value != "call" && value != "put") {
            throw std::runtime_error(where + ": option.type must be call or put");
        }
        cfg.option.is_call = (value == "call");
    } else if (key == "market.spot") {
        cfg.market.spot = parse_double(value, where);
    } else if (key == "market.rate") {
        cfg.market.rate = parse_double(value, where);
    } else if (key == "market.div_yield") {
        cfg.market.div_yield = parse_double(value, where);
    } else if (key == "heston.v0") {
        cfg.heston.v0 = parse_double(value, where);
    } else if (key == "heston.kappa") {
        cfg.heston.kappa = parse_double(value, where);
    } else if (key == "heston.theta") {
        cfg.heston.theta = parse_double(value, where);
    } else if (key == "heston.xi") {
        cfg.heston.xi = parse_double(value, where);
    } else if (key == "heston.rho") {
        cfg.heston.rho = parse_double(value, where);
    } else if (key == "grid.ns") {
        cfg.grid.num_stock_nodes = parse_int(value, where);
    } else if (key == "grid.nv") {
        cfg.grid.num_variance_nodes = parse_int(value, where);
    } else if (key == "grid.nt") {
        cfg.grid.num_timesteps = parse_int(value, where);
    } else if (key == "grid.s_max_mult") {
        cfg.grid.stock_max_multiplier = parse_double(value, where);
    } else if (key == "grid.v_max") {
        cfg.grid.variance_max = parse_double(value, where);
    } else {
        // Nothing matched, so the key is a typo
        throw std::runtime_error(where + ": unknown key '" + key + "'");
    }
}

// Reads a config file line by line and builds a Config from it.
// Blank lines and anything after a hash are ignored.
Config load_config(const std::string& path) {
    // The file closes itself at the end of the function, like Python's with open
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open config: " + path);
    }

    Config cfg;
    std::string line;
    int line_number = 0;
    while (std::getline(file, line)) {
        line_number = line_number + 1;

        // Cut off everything after a hash
        std::size_t hash_position = line.find('#');
        if (hash_position != std::string::npos) {
            line = line.substr(0, hash_position);
        }
        line = trim(line);
        if (line.empty()) {
            continue;
        }

        // File name and line number, used in error messages
        std::string where = path + ":" + std::to_string(line_number);

        // Split the line into the part before and after the equals sign
        std::size_t equals_position = line.find('=');
        if (equals_position == std::string::npos) {
            throw std::runtime_error(where + ": expected 'key = value'");
        }
        std::string key = trim(line.substr(0, equals_position));
        std::string value = trim(line.substr(equals_position + 1));
        apply(cfg, key, value, where);
    }
    return cfg;
}
