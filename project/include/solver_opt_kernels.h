#ifndef HESTON_SOLVER_OPT_KERNELS_H
#define HESTON_SOLVER_OPT_KERNELS_H

#include <vector>

#include "grid.h"
#include "params.h"

// One backward timestep for the optimisation ladder (PLAN §4b). Each kernel
// updates the interior AND the v=0 degenerate row of `next` from `cur`.
// The S=0 / S=max columns and the v=vmax row are shared plumbing applied by
// OptSolver::solve after the kernel — identical for every level, so the
// ladder timings compare kernels only.
//
// StepKernel is a function pointer: a plain variable that holds WHICH
// function to call. OptSolver picks it ONCE before the time loop, so the
// hot loop never pays a per-cell "which level?" branch — similar to
// assigning a function to a variable in Python.
using StepKernel = void (*)(const Grid& grid,
                            const std::vector<double>& cur,
                            std::vector<double>& next, const Config& cfg,
                            double dt);

// The kernel for a ladder level, or nullptr while that level's kernel is
// still unwritten. OptSolver::solve throws on nullptr — before any work,
// at the edge, per the §1b exception rule.
StepKernel kernel_for_level(int level);

#endif  // HESTON_SOLVER_OPT_KERNELS_H
