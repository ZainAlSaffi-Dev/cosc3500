# Graph Report - /Users/zer0/Documents/courses/cosc3500/matrixass  (2026-08-26)

## Corpus Check
- Corpus is ~42,377 words - fits in a single context window. You may not need a graph.

## Summary
- 157 nodes · 195 edges · 22 communities (11 shown, 11 thin omitted)
- Extraction: 90% EXTRACTED · 10% INFERRED · 0% AMBIGUOUS · INFERRED: 19 edges (avg confidence: 0.87)
- Token cost: 245,930 input · 0 output

## Community Hubs (Navigation)
- Assignment Spec Contract
- MPI Tutorial Collectives
- Marking and References
- Build Phases and OpenMP
- OpenMP Race Fixes
- CUDA Optimisation Techniques
- CPU Optimisation Ladder
- Rangpur and JudgementDay
- GPU Source Skeleton
- CUDA Basics
- CPU Source Skeleton
- MPI Source Skeleton
- CPU Bench Script
- GPU Bench Script
- MPI Bench Script
- Debug CPU Script
- Debug GPU Script
- Debug MPI Script
- JudgementDay Script
- Args Sweep Script
- Get Files Script

## God Nodes (most connected - your core abstractions)
1. `Assignment1_GradeBot Benchmark Harness` - 10 edges
2. `Parallel Reduction Algorithm` - 9 edges
3. `Phase 2 CPU Planar Rebuild` - 7 edges
4. `COSC3500 Assignment: Parallel Programming Techniques` - 7 edges
5. `Race Condition` - 6 edges
6. `OMP Reduction` - 6 edges
7. `OpenMP Numerical Integration Example` - 6 edges
8. `Marking Rubric (Table 1)` - 6 edges
9. `GradeBot Harness` - 5 edges
10. `GPU Ladder` - 5 edges

## Surprising Connections (you probably didn't know these)
- `GradeBot Harness` --semantically_similar_to--> `Assignment1_GradeBot`  [INFERRED] [semantically similar]
  MASTERDOC.md → docs/COSC3500_Assignment_-transcript.txt
- `Marking Scheme` --semantically_similar_to--> `Speed-Based Marking Table`  [INFERRED] [semantically similar]
  MASTERDOC.md → docs/COSC3500_Assignment_-transcript.txt
- `JudgementDay Graded Run` --semantically_similar_to--> `Judgement Day Slurm Script`  [INFERRED] [semantically similar]
  MASTERDOC.md → docs/COSC3500_Assignment_-transcript.txt
- `MPI Row Split Decomposition` --semantically_similar_to--> `Domain Decomposition with Remainder Handling`  [INFERRED] [semantically similar]
  MASTERDOC.md → docs/COSC3500_T05_MPI.mp4-transcript.txt
- `MPI_Allgatherv Reassembly` --semantically_similar_to--> `MPI Scatter and Gather`  [INFERRED] [semantically similar]
  MASTERDOC.md → docs/COSC3500_T05_MPI.mp4-transcript.txt

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Shared Numerical Integration Example** — docs_cosc3500_t04_omp_mp4_transcript_numerical_integration, docs_cosc3500_t05_mpi_mp4_transcript_numerical_integration, docs_cosc3500_t06_cuda_mp_transcript_integral_example [EXTRACTED 1.00]
- **OpenMP Race Condition Fixes Ranked** — docs_cosc3500_t04_omp_mp4_transcript_critical, docs_cosc3500_t04_omp_mp4_transcript_atomic, docs_cosc3500_t04_omp_mp4_transcript_partial_sums, docs_cosc3500_t04_omp_mp4_transcript_reduction [EXTRACTED 1.00]
- **CPU Optimisation Ladder** — masterdoc_cache_blocking, masterdoc_planar_split, masterdoc_4m_scheme, masterdoc_avx1_microkernel, masterdoc_3m_gauss_trick [EXTRACTED 1.00]
- **Three Matrix Multiply Implementations** — docs_cosc3500_assignment_cpu_implementation, docs_cosc3500_assignment_gpu_implementation, docs_cosc3500_assignment_mpi_implementation, docs_cosc3500_assignment_matrixmultiply_cpp, docs_cosc3500_assignment_matrixmultiplygpu_cu, docs_cosc3500_assignment_matrixmultiplympi_cpp [EXTRACTED 1.00]
- **GradeBot Grading Pipeline** — docs_cosc3500_assignment_gradebot, docs_cosc3500_assignment_judgementday_script, docs_cosc3500_assignment_marking_rubric, docs_cosc3500_assignment_gradebot_text_output, docs_cosc3500_assignment_chained_unitary_benchmark [EXTRACTED 1.00]
- **Final Submission Package** — docs_cosc3500_assignment_submission_format, docs_cosc3500_assignment_slurm_zip, docs_cosc3500_assignment_matrixmultiply_cpp, docs_cosc3500_assignment_matrixmultiplygpu_cu, docs_cosc3500_assignment_matrixmultiplympi_cpp [EXTRACTED 1.00]

## Communities (22 total, 11 thin omitted)

### Community 0 - "Assignment Spec Contract"
Cohesion: 0.12
Nodes (28): COSC3500 Assignment: Parallel Programming Techniques, Chained Unitary Matrix Benchmark, CPU Implementation (AVX/OpenMP), cuBLAS Reference Implementation, Debug Slurm Scripts, Three Editable Files Rule, ERR_TOLERANCE (rubric.h), Fixed matrixMultiply Function Signatures (+20 more)

### Community 1 - "MPI Tutorial Collectives"
Cohesion: 0.17
Nodes (16): MPI Allreduce, MPI Alltoall, MPI Barrier, MPI Broadcast, Cluster Job Submission, MPI Comm Size and Rank, MPI Communicator, Domain Decomposition with Remainder Handling (+8 more)

### Community 2 - "Marking and References"
Cohesion: 0.23
Nodes (13): Assignment Briefing Transcript, cuBLAS Reference Library, Assignment1_GradeBot, Speed-Based Marking Table, Intel MKL Reference, MPI Init and Finalize, MPI Tutorial Transcript, Assignment Master Plan (+5 more)

### Community 3 - "Build Phases and OpenMP"
Cohesion: 0.17
Nodes (12): Collapse Clause, Data Sharing Clauses, OpenMP Parallel Region, Schedule Clause, OMP Sections, OMP Single and Master, Work Sharing For Loop, Phase 0 Pipeline and Reconnaissance (+4 more)

### Community 4 - "OpenMP Race Fixes"
Cohesion: 0.27
Nodes (11): OMP Atomic, OMP Barrier, OMP Critical, Monte Carlo Sampling Example, OpenMP Numerical Integration Example, OpenMP Tutorial Transcript, Per-Thread Partial Sums, Race Condition (+3 more)

### Community 5 - "CUDA Optimisation Techniques"
Cohesion: 0.25
Nodes (11): Bank Conflicts, Cooperative Groups, CUDA Tutorial Transcript, Integral Under a Curve Example, Mark Harris Reduction PDF, Parallel Reduction Algorithm, Shared Memory, __syncthreads Fence (+3 more)

### Community 6 - "CPU Optimisation Ladder"
Cohesion: 0.25
Nodes (11): 3M Gauss Trick, 4M Real GEMM Scheme, AVX1 16x4 Microkernel, Cache Blocking, CPU Ladder, ERR_TOLERANCE, Iron Rule, Phase 2 CPU Planar Rebuild (+3 more)

### Community 7 - "Rangpur and JudgementDay"
Cohesion: 0.40
Nodes (6): Debug Slurm Scripts, Judgement Day Slurm Script, Rangpur Cluster, JudgementDay Graded Run, Phase 7 Dress Rehearsal, Phase 8 Freeze and Package

### Community 8 - "GPU Source Skeleton"
Cohesion: 0.40
Nodes (5): floatTypeCUDA, __global__, __host__, matrixMultiply_GPU(), matrixMultiplyKernel_GPU()

### Community 9 - "CUDA Basics"
Cohesion: 0.40
Nodes (5): CUDA Error Check Macro, cudaMalloc and cudaMemcpy, Default Stream Ordering, Global Thread Index, CUDA Kernel

## Knowledge Gaps
- **33 isolated node(s):** `COSC3500_GetAssignmentFiles.sh script`, `bench_cpu_2048.sh script`, `PATH`, `bench_gpu_2048.sh script`, `PATH` (+28 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **11 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Work Sharing For Loop` connect `Build Phases and OpenMP` to `OpenMP Race Fixes`?**
  _High betweenness centrality (0.064) - this node is a cross-community bridge._
- **Why does `OpenMP Numerical Integration Example` connect `OpenMP Race Fixes` to `MPI Tutorial Collectives`, `Build Phases and OpenMP`, `CUDA Optimisation Techniques`?**
  _High betweenness centrality (0.061) - this node is a cross-community bridge._
- **Why does `Phase 1 CPU Correct Baseline` connect `Build Phases and OpenMP` to `CPU Optimisation Ladder`?**
  _High betweenness centrality (0.057) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `Assignment1_GradeBot Benchmark Harness` (e.g. with `Debug Slurm Scripts` and `RangpurJudgementDay Grading Script`) actually correct?**
  _`Assignment1_GradeBot Benchmark Harness` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `COSC3500_GetAssignmentFiles.sh script`, `bench_cpu_2048.sh script`, `PATH` to the rest of the system?**
  _33 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Assignment Spec Contract` be split into smaller, more focused modules?**
  _Cohesion score 0.12169312169312169 - nodes in this community are weakly interconnected._