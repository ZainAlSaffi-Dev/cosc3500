#ifndef HESTON_IO_H
#define HESTON_IO_H

#include <string>

#include "grid.h"
#include "params.h"
#include "solver.h"

// Writes one sheet to results/<dir>/snap_<step>.csv, with one row per
// variance level and one column per stock node, which is the layout
// weather_map.py expects.
void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid);

// Prints the single stdout line that the benchmark scripts parse, holding
// price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec
void print_result_csv(const SolveResult& r, const Config& cfg);

// Gives the name a run identifies itself by, which is "baseline", "opt-L3" or
// one of the two control names. Keeping it in one function stops the CSV
// column and the stderr status lines from drifting apart.
std::string solver_label(const Config& cfg);

#endif  // HESTON_IO_H
