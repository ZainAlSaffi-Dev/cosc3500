# Graph Report - /Users/zer0/Documents/cosc3500/matrixass  (2026-08-07)

## Corpus Check
- Corpus is ~896 words - fits in a single context window. You may not need a graph.

## Summary
- 38 nodes · 51 edges · 7 communities (4 shown, 3 thin omitted)
- Extraction: 69% EXTRACTED · 31% INFERRED · 0% AMBIGUOUS · INFERRED: 16 edges (avg confidence: 0.85)
- Token cost: 50,291 input · 4,600 output

## Community Hubs (Navigation)
- Grading Pipeline (GradeBot)
- CPU & MPI Sources + Submission
- Benchmark References (MKL/cuBLAS)
- GPU CUDA Code
- Debug CPU Job
- Debug GPU Job
- Debug MPI Job

## God Nodes (most connected - your core abstractions)
1. `COSC3500 Assignment Specification: Parallel Programming Techniques` - 10 edges
2. `matrixMultiply()` - 6 edges
3. `matrixMultiply_GPU()` - 6 edges
4. `matrixMultiply_MPI()` - 6 edges
5. `Assignment1_GradeBot` - 5 edges
6. `matrixMultiplyKernel_GPU()` - 4 edges
7. `CPU (AVX/OpenMP) Implementation` - 4 edges
8. `MPI Implementation` - 4 edges
9. `Intel MKL Reference Implementation` - 4 edges
10. `Final Submission Requirements` - 4 edges

## Surprising Connections (you probably didn't know these)
- `Final Submission Requirements` --references--> `matrixMultiply_GPU()`  [INFERRED]
  COSC3500_Assignment.pdf → matrixMultiplyGPU.cu
- `COSC3500 Assignment Specification: Parallel Programming Techniques` --references--> `matrixMultiply_GPU()`  [INFERRED]
  COSC3500_Assignment.pdf → matrixMultiplyGPU.cu
- `matrixMultiply_GPU()` --implements--> `GPU (CUDA) Implementation`  [INFERRED]
  matrixMultiplyGPU.cu → COSC3500_Assignment.pdf
- `COSC3500 Assignment Specification: Parallel Programming Techniques` --references--> `matrixMultiplyKernel_GPU()`  [INFERRED]
  COSC3500_Assignment.pdf → matrixMultiplyGPU.cu
- `COSC3500 Assignment Specification: Parallel Programming Techniques` --references--> `matrixMultiply_MPI()`  [INFERRED]
  COSC3500_Assignment.pdf → matrixMultiplyMPI.cpp

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Automated Grading Pipeline** — slurm_goslurm_cosc3500assignment_rangpurjudgementday, cosc3500_assignment_gradebot, cosc3500_assignment_benchmark_procedure, cosc3500_assignment_marking_rubric [EXTRACTED 1.00]
- **Three Hardware Target Implementations** — cosc3500_assignment_cpu_implementation, cosc3500_assignment_gpu_implementation, cosc3500_assignment_mpi_implementation [EXTRACTED 1.00]
- **Per-Section Debug Slurm Workflow** — slurm_goslurm_cosc3500assignment_rangpurdebugcpu, slurm_goslurm_cosc3500assignment_rangpurdebuggpu, slurm_goslurm_cosc3500assignment_rangpurdebugmpi, cosc3500_assignment_gradebot [INFERRED 0.85]

## Communities (7 total, 3 thin omitted)

### Community 0 - "Grading Pipeline (GradeBot)"
Cohesion: 0.25
Nodes (8): Benchmark Procedure, ERR_TOLERANCE, Assignment1_GradeBot, Rangpur Cluster, COSC3500 Assignment Specification: Parallel Programming Techniques, Benchmark Text Output Format, goslurm_COSC3500Assignment_RangpurJudgementDay script, PATH

### Community 1 - "CPU & MPI Sources + Submission"
Cohesion: 0.29
Nodes (6): Final Submission Requirements, Modifiable Files Constraint, floatType, matrixMultiply(), floatType, matrixMultiply_MPI()

### Community 2 - "Benchmark References (MKL/cuBLAS)"
Cohesion: 0.47
Nodes (6): CPU (AVX/OpenMP) Implementation, CUBLAS Reference Implementation, GPU (CUDA) Implementation, Intel MKL Reference Implementation, Marking Rubric (Table 1), MPI Implementation

### Community 3 - "GPU CUDA Code"
Cohesion: 0.40
Nodes (5): floatTypeCUDA, __global__, __host__, matrixMultiply_GPU(), matrixMultiplyKernel_GPU()

## Knowledge Gaps
- **3 isolated node(s):** `ERR_TOLERANCE`, `Benchmark Text Output Format`, `Rangpur Cluster`
  These have ≤1 connection - possible missing edges or undocumented components.
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `COSC3500 Assignment Specification: Parallel Programming Techniques` connect `Grading Pipeline (GradeBot)` to `CPU & MPI Sources + Submission`, `GPU CUDA Code`, `Debug CPU Job`, `Debug GPU Job`, `Debug MPI Job`?**
  _High betweenness centrality (0.632) - this node is a cross-community bridge._
- **Why does `Assignment1_GradeBot` connect `Grading Pipeline (GradeBot)` to `Benchmark References (MKL/cuBLAS)`?**
  _High betweenness centrality (0.195) - this node is a cross-community bridge._
- **Why does `matrixMultiply_GPU()` connect `GPU CUDA Code` to `Grading Pipeline (GradeBot)`, `CPU & MPI Sources + Submission`, `Benchmark References (MKL/cuBLAS)`?**
  _High betweenness centrality (0.166) - this node is a cross-community bridge._
- **Are the 6 inferred relationships involving `COSC3500 Assignment Specification: Parallel Programming Techniques` (e.g. with `matrixMultiply_GPU()` and `matrixMultiplyKernel_GPU()`) actually correct?**
  _`COSC3500 Assignment Specification: Parallel Programming Techniques` has 6 INFERRED edges - model-reasoned connections that need verification._
- **Are the 3 inferred relationships involving `matrixMultiply_GPU()` (e.g. with `Final Submission Requirements` and `COSC3500 Assignment Specification: Parallel Programming Techniques`) actually correct?**
  _`matrixMultiply_GPU()` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `matrixMultiply_MPI()` (e.g. with `Final Submission Requirements` and `COSC3500 Assignment Specification: Parallel Programming Techniques`) actually correct?**
  _`matrixMultiply_MPI()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **What connects `ERR_TOLERANCE`, `Benchmark Text Output Format`, `Rangpur Cluster` to the rest of the system?**
  _3 weakly-connected nodes found - possible documentation gaps or missing edges._