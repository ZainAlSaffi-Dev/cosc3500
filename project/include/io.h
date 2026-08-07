#ifndef HESTON_IO_H
#define HESTON_IO_H

#include <string>

#include "grid.h"
#include "params.h"
#include "solver.h"

// Snapshot one sheet as CSV: results/<dir>/snap_<step>.csv
// Row = variance level j, column = stock node i (matches weather_map.py).
void dump_snapshot(const std::string& dump_dir, int step, const Grid& grid);

// One machine-readable stdout line (CLI contract, PLAN §3):
// price,delta,gamma,vega,ns,nv,nt,solver,seconds,cell_updates_per_sec
void print_result_csv(const SolveResult& r, const Config& cfg);

#endif  // HESTON_IO_H
