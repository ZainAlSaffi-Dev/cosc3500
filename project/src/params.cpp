#include "params.h"

#include <stdexcept>

Config load_config(const std::string& path) {
    // TODO(P1): open file, parse "key = value" lines, '#' comments,
    // map dotted keys (option.strike, heston.xi, grid.ns, ...) onto Config.
    // Unknown key or unparsable value -> std::runtime_error with line number.
    (void)path;
    throw std::runtime_error("load_config: not implemented (P1)");
}
