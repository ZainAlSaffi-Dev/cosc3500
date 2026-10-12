// AI declaration
// I used Claude, an AI assistant made by Anthropic, while working on this file.
// It helped me debug the code, draft pseudocode and implement it.

#ifndef HESTON_IO_H
#define HESTON_IO_H

#include <string>

#include "grid.h"
#include "params.h"
#include "solver.h"

// Saves the current sheet as a CSV file for the animation scripts.
void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid);

// Prints the one CSV line of results that the benchmark scripts read.
void print_result_csv(const SolveResult& result, const Config& cfg);

#endif  // HESTON_IO_H
