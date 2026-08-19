#ifndef HESTON_SOLVER_OPT_KERNELS_H
#define HESTON_SOLVER_OPT_KERNELS_H

#include <vector>

#include "grid.h"
#include "params.h"

// One backward timestep for the optimisation ladder. Each kernel updates the
// interior and the v=0 row of `next` from `cur`, while the outer boundaries
// are handled by OptSolver::solve once the kernel returns. That outer part is
// identical for every level, so the ladder timings only ever compare kernels.
//
// StepKernel is a function pointer, meaning a variable that holds which
// function to call. OptSolver chooses it once before the time loop so that
// the hot loop never pays for a per-cell test of which level is running. This
// is similar to assigning a function to a variable in Python.
using StepKernel = void (*)(const Grid& grid,
                            const std::vector<double>& cur,
                            std::vector<double>& next, const Config& cfg,
                            double dt);

// Returns the kernel for a ladder level, or nullptr if that level has not
// been written yet.
StepKernel kernel_for_level(int level);

#endif  // HESTON_SOLVER_OPT_KERNELS_H
