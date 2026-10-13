# Graph Report - .  (2026-10-11)

## Corpus Check
- 67 files · ~174,415 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 715 nodes · 1045 edges · 69 communities (53 shown, 16 thin omitted)
- Extraction: 89% EXTRACTED · 10% INFERRED · 0% AMBIGUOUS · INFERRED: 109 edges (avg confidence: 0.82)
- Token cost: 755,053 input · 0 output

## Community Hubs (Navigation)
- Grid Class and Stencil Core
- MPI Lecture Concepts
- Benchmark Plotting Script
- Presentation Deck Content
- Solver Class Hierarchy
- Milestone Spec Requirements
- Convergence Plot Script
- Memory Traffic Animation
- Grid Spec and Axis Lookup
- CLI Entry and Params
- Weather Map Animation
- Xi Sweep Smile Script
- Serial Optimisation Techniques
- Code Style and Ownership
- OpenMP Thread Control API
- BS Collapse Animation
- Code Shot Renderer
- Heston Spec Foundations
- OpenMP Fork-Join Basics
- Validation Checks
- SIMD Instruction Sets
- Optimisation Ladder Results
- SIMD Load Store Alignment
- OpenMP Loop Scheduling
- OpenMP Races and Overhead
- Explicit Euler Numerics
- Stability and Sweep Findings
- Benchmark Hygiene Rules
- SIMD FMA and Reductions
- Deck Timer and Budget
- Handoff and Docs Index
- Milestone 2 Scope
- CUDA Shared Memory
- OpenMP Barriers and Single
- CUDA Warp Reduction
- Cache Loop Order and AVX Plan
- Boundaries and Correctness
- AVX Intrinsics Basics
- OpenMP Critical and Atomic
- OpenMP Reduction and SIMD
- CUDA Kernel Indexing
- CUDA Sync and Divergence
- Benchmark Protocol Traps
- Slide Stills Grabber
- SIMD Masks and Blend
- OpenMP Variable Scoping
- CUDA Memory Basics
- Slurm Common Helpers
- Double Buffer Time March
- SIMD Arithmetic Costs
- OpenMP Parallel Plan
- Rangpur Job Limits
- Fused Boundary Control
- Slurm Flags Sweep
- Slurm Novec O3 Sweep
- Slurm Page Size Sweep
- Slurm Scaling Sweep
- Stencil Narration
- Slurm Controls Job
- Slurm Ladder Job
- Slurm Novec Job
- Slurm Peak Job
- Slurm Debug Serial
- Slurm Alloc Diagnostic
- Slurm Denormal Diagnostic
- Harness Control Level
- Induction Variable Level
- Grid Convergence Concept

## God Nodes (most connected - your core abstractions)
1. `Config` - 34 edges
2. `Grid` - 27 edges
3. `OptSolver::solve()` - 16 edges
4. `Solver` - 15 edges
5. `BaselineSolver::solve()` - 15 edges
6. `SolveResult` - 14 edges
7. `COSC3500 L06 MPI slides (Ian McCulloch, 2022)` - 14 edges
8. `Ablation ladder --opt-level 0..6` - 13 edges
9. `COSC3500 L06 MPI lecture transcript` - 13 edges
10. `COSC3500 T05 MPI tutorial transcript` - 12 edges

## Surprising Connections (you probably didn't know these)
- `thread_local random number generator` --semantically_similar_to--> `threadprivate`  [INFERRED] [semantically similar]
  docs/lecs/COSC3500_T04_OMP.mp4-transcript.txt → docs/lecs/COSC3500_L05_openMP_20220823.pdf
- `Reference grid fits in 30 MB L3` --conceptually_related_to--> `Benchmark protocol (rangpur compute node, pinned, 5 sweeps)`  [AMBIGUOUS]
  STUDY_GUIDE.md → PLAN.md
- `schedule(static)` --semantically_similar_to--> `Manual block decomposition (istart/iend per thread)`  [INFERRED] [semantically similar]
  docs/lecs/COSC3500_T04_OMP.mp4-transcript.txt → docs/lecs/COSC3500_L05_openMP_20220823.pdf
- `Broadcast size then data for dynamic arrays` --semantically_similar_to--> `Dynamic size buffers (MPI_Get_count, MPI_Probe)`  [INFERRED] [semantically similar]
  docs/lecs/COSC3500_L06_MPI_202-transcript.txt → docs/lecs/COSC3500_L06_MPI_20220901.pdf
- `Rank-based block decomposition (first/last index per rank)` --semantically_similar_to--> `MPI_Scatter / MPI_Scatterv`  [INFERRED] [semantically similar]
  docs/lecs/COSC3500_T05_MPI.mp4-transcript.txt → docs/lecs/COSC3500_L06_MPI_20220901.pdf

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Marking Rubric Criteria** — cosc3500_milestone_marking_rubric, cosc3500_milestone_introduction_and_background, cosc3500_milestone_optimisation, cosc3500_milestone_benchmarking, cosc3500_milestone_presentation_criterion, cosc3500_milestone_reflection_conclusion, cosc3500_milestone_code_submission, cosc3500_milestone_identity_verification, cosc3500_milestone_in_person_interview, cosc3500_milestone_submission_format [EXTRACTED 1.00]
- **Parallel Programming Technique Options for Milestone 2** — cosc3500_milestone_parallel_implementation, cosc3500_milestone_avx, cosc3500_milestone_openmp, cosc3500_milestone_mpi, cosc3500_milestone_cuda [EXTRACTED 1.00]
- **Opt-level ladder rungs 1-6** — plan_loop_invariant_hoisting, plan_strength_reduction, plan_cache_friendly_loop_order, plan_loop_splitting, plan_induction_variable_simplification, plan_loop_unrolling [EXTRACTED 1.00]
- **Correctness validation suite** — study_guide_black_scholes_collapse, study_guide_put_call_parity, study_guide_grid_convergence, study_guide_fourier_reference [EXTRACTED 1.00]
- **Measured memory effects (page size, traversal order, cache fit, roofline)** — plan_page_size_trap, plan_negative_controls, study_guide_llc_fit, study_guide_roofline, study_guide_cache_traffic_model [INFERRED 0.85]
- **Correctness gates validating the solver** — docs_initialspecdesign_bs_collapse_check, docs_initialspecdesign_put_call_parity, docs_results_test_opt_matches, docs_results_check_heston, docs_results_convergence_floor, docs_results_fourier_mc_check [EXTRACTED 1.00]
- **Benchmark provenance and validity discipline** — docs_running_build_and_validate, docs_running_write_header_provenance, docs_running_benchmark_hygiene, docs_results_r730_2_node, docs_initialspecdesign_iron_rule [INFERRED 0.85]
- **Compiler and vectorisation attribution controls** — docs_results_sweep_c_vacuous, docs_results_sweep_f_flags, docs_results_sweep_g_vectorisation, docs_results_compiler_forbidden_finding [EXTRACTED 1.00]
- **Seven-level optimisation ladder rungs** — presentation_speaker_notes_level0_harness_control, presentation_speaker_notes_level1_hoisting, presentation_presentation_strength_reduction, presentation_speaker_notes_level3_traversal_order, presentation_speaker_notes_level4_loop_splitting, presentation_speaker_notes_level5_induction_variables, presentation_speaker_notes_level6_unrolling [EXTRACTED 1.00]
- **Video recording document set** — presentation_presentation, presentation_timeline, presentation_script, presentation_shotlist, presentation_speaker_notes, presentation_presentation_deck_tex [EXTRACTED 1.00]
- **Model-to-equation derivation run** — presentation_presentation_heston_model, presentation_presentation_feynman_kac, presentation_script_curvature_argument, presentation_presentation_backwards_grid_solve [INFERRED 0.85]
- **CUDA reduction optimisation ladder v1-v4** — docs_lecs_cosc3500_t06_cuda_mp_transcript_single_block_tree_reduction, docs_lecs_cosc3500_t06_cuda_mp_transcript_multi_block_reduction, docs_lecs_cosc3500_t06_cuda_mp_transcript_shared_memory_reduction, docs_lecs_cosc3500_t06_cuda_mp_transcript_vectorized_loads, docs_lecs_cosc3500_t06_cuda_mp_transcript_cuda_reduction [EXTRACTED 1.00]
- **AVX complex dot-product pipeline** — docs_lecs_cosc3500_l04_simd_mp_transcript_complex_dot_product_demo, docs_lecs_cosc3500_l04_simd_mp_transcript_memory_alignment, docs_lecs_cosc3500_l04_simd_mp_transcript_set_pseudo_intrinsics, docs_lecs_cosc3500_l04_simd_mp_transcript_permute, docs_lecs_cosc3500_l04_simd_mp_transcript_fused_multiply_add, docs_lecs_cosc3500_l04_simd_mp_transcript_horizontal_reduction [EXTRACTED 1.00]
- **SIMD branchless select via masks** — docs_lecs_cosc3500_l04_simd_mp_transcript_comparison_masks, docs_lecs_cosc3500_l04_simd_mp_transcript_logical_ops, docs_lecs_cosc3500_l04_simd_mp_transcript_masked_branch_select, docs_lecs_cosc3500_l04_simd_mp_transcript_blend [INFERRED 0.85]
- **Ways to fix a racing shared accumulation (critical, atomic, partial sums, reduction)** — docs_lecs_cosc3500_l05_openmp_20220823_pragma_omp_critical, docs_lecs_cosc3500_t04_omp_mp4_transcript_atomic, docs_lecs_cosc3500_l05_openmp_20220823_local_partial_sum, docs_lecs_cosc3500_l05_openmp_20220823_reduction_clause, docs_lecs_cosc3500_l05_openmp_20220823_race_condition [EXTRACTED 1.00]
- **Controlling and querying OpenMP thread count** — docs_lecs_cosc3500_l05_openmp_20220823_omp_num_threads_env, docs_lecs_cosc3500_l05_openmp_20220823_omp_set_num_threads, docs_lecs_cosc3500_l05_openmp_20220823_num_threads_clause, docs_lecs_cosc3500_l05_openmp_20220823_omp_get_max_threads, docs_lecs_cosc3500_l05_openmp_20220823_omp_get_num_threads, docs_lecs_cosc3500_t04_omp_mp4_transcript_omp_set_dynamic [EXTRACTED 1.00]
- **Distributing loop iterations across threads** — docs_lecs_cosc3500_l05_openmp_20220823_pragma_omp_for, docs_lecs_cosc3500_t04_omp_mp4_transcript_schedule_clause, docs_lecs_cosc3500_t04_omp_mp4_transcript_static_schedule, docs_lecs_cosc3500_t04_omp_mp4_transcript_dynamic_schedule, docs_lecs_cosc3500_t04_omp_mp4_transcript_chunk_size, docs_lecs_cosc3500_t04_omp_mp4_transcript_collapse_clause [EXTRACTED 1.00]
- **MPI collective operations family** — docs_lecs_cosc3500_l06_mpi_20220901_mpi_barrier, docs_lecs_cosc3500_l06_mpi_20220901_mpi_bcast, docs_lecs_cosc3500_l06_mpi_20220901_mpi_scatter, docs_lecs_cosc3500_l06_mpi_20220901_mpi_gather, docs_lecs_cosc3500_l06_mpi_20220901_mpi_allgather, docs_lecs_cosc3500_l06_mpi_20220901_mpi_reduce, docs_lecs_cosc3500_l06_mpi_20220901_mpi_allreduce, docs_lecs_cosc3500_t05_mpi_mp4_transcript_mpi_alltoall [EXTRACTED 1.00]
- **Non-blocking overlap pattern (Irecv, compute/send, Wait)** — docs_lecs_cosc3500_l06_mpi_20220901_nonblocking_communication, docs_lecs_cosc3500_l06_mpi_20220901_mpi_wait_test, docs_lecs_cosc3500_l06_mpi_20220901_overlap_comm_computation, docs_lecs_cosc3500_l06_mpi_202_transcript_async_progress_pitfall [EXTRACTED 1.00]
- **Split work by rank, compute locally, reduce result workflow** — docs_lecs_cosc3500_t05_mpi_mp4_transcript_rank_based_domain_decomposition, docs_lecs_cosc3500_t05_mpi_mp4_transcript_distributed_memory_footprint, docs_lecs_cosc3500_t05_mpi_mp4_transcript_manual_master_receive_loop, docs_lecs_cosc3500_t05_mpi_mp4_transcript_reduce_allreduce_refactor, docs_lecs_cosc3500_l06_mpi_20220901_mpi_allreduce [EXTRACTED 1.00]

## Communities (69 total, 16 thin omitted)

### Community 0 - "Grid Class and Stencil Core"
Cohesion: 0.05
Nodes (59): Grid, current_, init_payoff, next_, stock_position, swap_buffers, variance_position, vector (+51 more)

### Community 1 - "MPI Lecture Concepts"
Cohesion: 0.07
Nodes (56): Collective communication, MPI communicator and rank, COSC3500 L06 MPI slides (Ian McCulloch, 2022), Hybrid distributed memory (MPI + OpenMP), Message passing ring example, MPI message tag, Message Passing Interface (MPI), MPI_Allgather (+48 more)

### Community 2 - "Benchmark Plotting Script"
Cohesion: 0.08
Nodes (28): aggregate(), check_provenance(), main(), plot_ladder(), plot_scaling(), provenance_caption(), provenance_facts(), Group rows by `key(row)` -> dict of median/min/max cell-updates/sec. (+20 more)

### Community 3 - "Presentation Deck Content"
Cohesion: 0.07
Nodes (36): config/bench.cfg, deck.html (earlier HTML recording deck with take timer), Reserved webcam outline (#pip, toggled with p), PRESENTATION.md (Milestone 1 video plan), Benchmark protocol (sbatch --exclusive, --nodelist=r730-2, median of five, CSV provenance header), xi=0 collapse to Black-Scholes validation, At -O2 the hand optimisations worth doing are the ones the compiler is not allowed to do, deck.tex (LaTeX recording deck, built by make slides) (+28 more)

### Community 4 - "Solver Class Hierarchy"
Cohesion: 0.10
Nodes (22): BaselineSolver, solve, string, OptSolver, solve, Solver, extract_result, name (+14 more)

### Community 5 - "Milestone Spec Requirements"
Cohesion: 0.10
Nodes (28): Amdahl's Law, AVX, Benchmarking (25%), Code Submission (pass/fail), Course FAQ, CPU vs GPU Comparison Approach, CUDA, Gustafson's Law (+20 more)

### Community 6 - "Convergence Plot Script"
Cohesion: 0.14
Nodes (23): check_alignment(), main(), One solve through the PLAN §3 CLI. Returns the parsed CSV result., The 'key = value' lines of a .cfg, as floats where they parse., Complain about any rung whose quote point misses a grid node.      Reads the geo, read_config(), run_solver(), black_scholes_call() (+15 more)

### Community 7 - "Memory Traffic Animation"
Cohesion: 0.13
Nodes (18): draw_buffers_scene(), draw_panel(), main(), measured_slowdown(), ModelCache, The nine cells the interior update reads, as flat indices. These are     the ten, Cell visit sequence. row_major=True is the solver's order (outer     variance, i, One memory panel: the flat array drawn the way it is laid out, so     reading le (+10 more)

### Community 8 - "Grid Spec and Axis Lookup"
Cohesion: 0.13
Nodes (17): AxisPosition, cell, weight, GridSpec, num_stock_nodes, num_timesteps, num_variance_nodes, stock_max_multiplier (+9 more)

### Community 9 - "CLI Entry and Params"
Cohesion: 0.24
Nodes (13): string, flag_value(), main(), usage(), apply(), string, load_config(), parse_double() (+5 more)

### Community 10 - "Weather Map Animation"
Cohesion: 0.21
Nodes (13): apply_config_defaults(), cir_variance_sd(), load_snapshots(), main(), Parse the same "key = value" .cfg the C++ reads (see src/params.cpp).      Every, Fill unset CLI flags from --config. Explicit flags always win., Standard deviation of the Heston variance v_T, in closed form.      The variance, The (S, v) box the market can plausibly reach, i.e. the shot worth     filming. (+5 more)

### Community 11 - "Xi Sweep Smile Script"
Cohesion: 0.21
Nodes (13): brent(), bs_call(), implied_vol(), main(), norm_cdf(), Brent's method: bisection's guaranteed bracketing plus the speed of     inverse, Volatility that makes Black-Scholes reproduce `price`, or NaN.      NaN (not an, Run ./heston once; return its single stdout CSV line. (+5 more)

### Community 12 - "Serial Optimisation Techniques"
Cohesion: 0.18
Nodes (13): Cell-updates per second metric, heston CLI contract and CSV output, Dispatch kernel once outside the time loop, Equal-work grid sweep (nt and maturity scaled per grid), Floating-point contract across the ladder, gprof profiling + chrono instrumentation, Level 5: induction-variable simplification, Level 1: loop-invariant hoisting + lookup table + CSE (+5 more)

### Community 13 - "Code Style and Ownership"
Cohesion: 0.21
Nodes (13): Explicit ownership (RAII, vector/unique_ptr), Working agreement: learning by building, Milestone 1 Build Plan (PLAN.md), Narration comments rule, OOP-shaped modern C++, memory-safe code style, Value vs reference semantics, Heston PDE Option Pricer README, One owner document per topic (+5 more)

### Community 14 - "OpenMP Thread Control API"
Cohesion: 0.17
Nodes (10): num_threads clause, omp_get_num_threads(), OMP_NUM_THREADS environment variable, omp_set_num_threads(), #pragma omp parallel for shortcut, #pragma omp parallel, Nested parallel regions, omp_set_dynamic() (+2 more)

### Community 15 - "BS Collapse Animation"
Cohesion: 0.24
Nodes (11): bs_call(), divisor_near(), main(), norm_cdf(), Largest-but-closest divisor of `total` to `target`.      The final snapshot must, Run ./heston and return (price, stderr_line). stdout is the one CSV line., Pull a single row out of a snapshot CSV.      A snapshot is nv rows (variance) b, Standard normal CDF, vectorised. Same identity as src/black_scholes.cpp:     N(x (+3 more)

### Community 16 - "Code Shot Renderer"
Cohesion: 0.30
Nodes (11): compose(), lexer_for(), main(), Path, Greedy word wrap so a caption longer than the panes folds onto     extra lines i, Panes side by side, one caption centred underneath., Return (text, mapping) where text joins the requested 1-based line     ranges wi, render_pane() (+3 more)

### Community 17 - "Heston Spec Foundations"
Cohesion: 0.20
Nodes (11): Learning debt (quiz backlog), InitialSpecDesign.md (project brief), Black-Scholes model, Collapse to Black-Scholes check (xi=0), Heston stochastic-volatility model, Heston pricing PDE, Implementation ladder: serial, serial-opt, SIMD, OpenMP, CUDA, Explicit nine-point stencil (+3 more)

### Community 18 - "OpenMP Fork-Join Basics"
Cohesion: 0.24
Nodes (11): Fork-join programming model (master thread), Ian McCulloch, COSC3500 Introduction to OpenMP (L05 slides), OpenMP, Shared-address multithreading model, Hyper-threading limits FP-bound speedup, L05 OpenMP lecture transcript, Thread oversubscription (+3 more)

### Community 19 - "Validation Checks"
Cohesion: 0.27
Nodes (10): demo.cfg 196.1054098 interpolation regression gate, Grid-refinement convergence check, Put-call parity check, Bilinear interpolated readout fix (was nearest cell, +4.67%), make check-heston (Fourier anchor for xi and rho terms), Convergence order ~2.2 then error floor ~$0.014, Fourier inversion and Monte Carlo cross-check (make mc-check), Parity is model-blind (residual unchanged across xi, rho) (+2 more)

### Community 20 - "SIMD Instruction Sets"
Cohesion: 0.22
Nodes (10): L04 SIMD Lecture Transcript, AVX (256-bit, 8 float32 / 4 float64), AVX2 (integer ops, full permute, gather), Compiler auto-vectorisation (-O3 / /O2), CPUID instruction set detection, Gather (AVX2) / scatter (AVX-512), Manual loop unrolling for vectorisation, MIMD / Multi-threading (+2 more)

### Community 21 - "Optimisation Ladder Results"
Cohesion: 0.25
Nodes (9): Reference configuration (SPX call K=5250, 2048x512), Mcell-updates/s throughput metric, Negative controls ctl-order (2.1x) and ctl-branch (1.0x), Optimisation ladder L0-L6 (sweep E, 2.56x), Config files (reference, bench, demo, smoke, unstable), heston CLI flags (--solver, --opt-level, --bench, --dump-every), Kernel chosen once via function pointer before time loop, One technique per full kernel copy design (+1 more)

### Community 22 - "SIMD Load Store Alignment"
Cohesion: 0.25
Nodes (9): Vector casts (zero-instruction reinterpretation), Aligned/unaligned load and store intrinsics, __m256 / __m256d / __m256i vector types, Aligned memory allocation (aligned_alloc to 32/64 bytes), Remainder (tail) loop for N not multiple of vector width, Data type conversion intrinsics, Keep sub-buffer offsets 256-bit aligned, Tail handling with switch fall-through for N%4 (+1 more)

### Community 23 - "OpenMP Loop Scheduling"
Cohesion: 0.28
Nodes (9): Manual block decomposition (istart/iend per thread), #pragma omp for (worksharing loop), private clause, Schedule chunk size, firstprivate clause, lastprivate clause, schedule clause, Shared loop counters race in parallel region (+1 more)

### Community 24 - "OpenMP Races and Overhead"
Cohesion: 0.25
Nodes (9): Parallelize outermost loops first, Race condition, Synchronization, Interleaved printf / IO across threads, collapse(n) clause, schedule(dynamic), Monte Carlo estimation of pi and e, Parallel overhead dominates for small problem sizes (+1 more)

### Community 25 - "Explicit Euler Numerics"
Cohesion: 0.25
Nodes (9): Phases P1-P8 with acceptance criteria, Explicit stability estimate dt_stable, Explicit Euler + central differences scheme, Heston PDE (2-D stock x variance grid), reference.cfg NaN failure and fix (nt=980000, v_max=0.64), Feynman-Kac theorem (why no randomness in code), Fourier semi-analytic + Monte Carlo external reference, Monte Carlo vs PDE route (+1 more)

### Community 26 - "Stability and Sweep Findings"
Cohesion: 0.29
Nodes (8): 2026-08-10 audit open items, Instability horror show (checkerboard blow-up), Explicit-scheme stability limit (CFL), Limitations and threats to validity, Huge-page effect (1.66x from 2 MiB THP), Sweep B problem-size scaling, Sweep C vacuous -O2 no-vectorise control, Sweep G vectorisation control (-O3 vs -fno-tree-vectorize)

### Community 27 - "Benchmark Hygiene Rules"
Cohesion: 0.29
Nodes (8): Iron rule: no parallel version benchmarked until it matches serial, Denormal hypothesis ruled out (diag_denormal.sh), rangpur r730-2 Xeon E5-2670 v3 benchmark node, Sweep F compiler-flag control (-O3, -ffast-math), test_opt_matches gate, Benchmark hygiene (--exclusive, --nodelist=r730-2, median of 5), build_and_validate() validation gate in bench_common.sh, write_header() CSV provenance stamp

### Community 28 - "SIMD FMA and Reductions"
Cohesion: 0.32
Nodes (8): Complex dot-product AVX demo, FMA3 instruction set, Fused multiply-add (_mm256_fmadd_ps), Horizontal add (_mm256_hadd_ps), Horizontal reduction of a vector accumulator, Memory-bound vs compute-bound limit, Permute intrinsics (permutevar8x32, permute4x64), set / set1 pseudo-intrinsics

### Community 29 - "Deck Timer and Budget"
Cohesion: 0.25
Nodes (8): clock() time formatter, Deck HUD take timer (per-slide budget via data-dur, red when over), show() slide navigation function, tick() HUD clock update function, Asset index (results/*.mp4, *.png, code_shots), Milestone 1 marking criteria (Intro 20, Optimisation 25, Benchmarking 25, Presentation 20, Reflection 10), Time budget at 150 words per minute (9 min target, 10 min ceiling), What was cut and why (strike slide, convergence, demo_surface, blowup, memory_buffers)

### Community 30 - "Handoff and Docs Index"
Cohesion: 0.33
Nodes (7): HANDOFF.md (next session pickup), arm64 Python venv under Rosetta, P8 10-minute video deliverable, Milestone 1 scope (serial only), Weather-map animation, RESULTS.md (measured findings), RUNNING.md (build, validate, reproduce)

### Community 31 - "Milestone 2 Scope"
Cohesion: 0.29
Nodes (7): Market-episode replay (Feb vs Mar 2020), Milestone 2 scope (SIMD, OpenMP, interview), Stress dial interactive demo, Finding: worthwhile hand optimisations are ones -O2 is forbidden to do, Roofline analysis (compute-bound, 82% of scalar peak), Level 2 strength reduction (+148.7%), Figure and animation make targets

### Community 32 - "CUDA Shared Memory"
Cohesion: 0.29
Nodes (7): 128-bit lane restriction of AVX permutes, Insert / extract 128-bit lane, Shared-memory bank conflicts (32 banks), CUDA memory hierarchy (registers, shared, global), cudaDeviceSetSharedMemConfig 8-byte banks, __shared__ memory (static and extern dynamic), v3 shared-memory reduction

### Community 33 - "OpenMP Barriers and Single"
Cohesion: 0.48
Nodes (7): #pragma omp barrier, exchange_boundaries() serial step in parallel region, Implicit barrier at end of worksharing loop, #pragma omp master, nowait clause, #pragma omp single, T04 OpenMP tutorial transcript

### Community 34 - "CUDA Warp Reduction"
Cohesion: 0.33
Nodes (7): T06 CUDA Reduction Tutorial Transcript, Cooperative groups (incl. grid-wide sync), Parallel reduction in CUDA, Mark Harris NVIDIA 'Optimizing Parallel Reduction' slides, Rangpur GPU node with module cuda/11.4, Volta independent thread scheduling (CC 7.0), Warp-level primitives (__syncwarp, shuffle)

### Community 35 - "Cache Loop Order and AVX Plan"
Cohesion: 0.43
Nodes (7): Level 3: cache-friendly loop order, Level 4: hoist conditionals / loop splitting, Milestone 2 SIMD via AVX intrinsics, Negative controls ctl-order and ctl-branch, Flat row-major grid idx = j*Ns + i, Cache traffic animation (8.5 B vs 41.2 B per cell), Roofline test: kernel is compute-bound

### Community 36 - "Boundaries and Correctness"
Cohesion: 0.33
Nodes (7): Feller condition check (2 kappa theta >= xi^2), S=0, S=Smax and v=vmax Neumann boundaries, v=0 degenerate transport boundary, Four independent correctness checks, Black-Scholes collapse validation (xi=0, v0=theta), Put-call parity validation, Nearest-cell quote-point blind spot

### Community 37 - "AVX Intrinsics Basics"
Cohesion: 0.33
Nodes (6): Accidental AVX-512 intrinsic crash pitfall, Hand-written assembly, AVX-512 (512-bit, fragmented sub-sets), immintrin.h header, Intrinsic functions, Intrinsic naming convention (_mm256_<op>_<ps|pd|epi32>)

### Community 38 - "OpenMP Critical and Atomic"
Cohesion: 0.47
Nodes (6): Thread-local partial sum then single critical combine, Named critical sections, Non-deterministic floating-point result, #pragma omp critical, Critical section in loop: 21 s parallel vs 0.5 s serial, #pragma omp atomic (update/read/write)

### Community 39 - "OpenMP Reduction and SIMD"
Cohesion: 0.40
Nodes (6): MIMD (Multiple Instruction Multiple Data), Numerical integration of 1/(1+x^2) = pi/4 (pi_serial.cpp), reduction clause, SIMD (Single Instruction Multiple Data), Integration benchmark: 561 ms serial to 71 ms on 8 threads, OpenMP simd directive (parallel for simd)

### Community 40 - "CUDA Kernel Indexing"
Cohesion: 0.33
Nodes (6): Ceil-division block count (N+T-1)/T, Default-stream kernel ordering (implicit sync), Global index blockDim.x*blockIdx.x+threadIdx.x, Grid-stride loop (stride gridDim.x*blockDim.x), v2 multi-block reduction with repeated kernel launches (~100 us), Ping-pong buffer pointer swap between launches

### Community 41 - "CUDA Sync and Divergence"
Cohesion: 0.33
Nodes (6): Read-after-write race between threads without barrier, Sequential addressing (halving stride, tid < i), v1 single-block tree reduction (~245 us), __syncthreads block barrier, __syncthreads inside conditional deadlocks, Interleaved addressing leaves warps mostly idle

### Community 42 - "Benchmark Protocol Traps"
Cohesion: 0.40
Nodes (6): Benchmark protocol (rangpur compute node, pinned, 5 sweeps), Clock-frequency ruled out (3100 MHz pinned), Milestone 1 definition of done, Page-size trap (huge pages vs 4 KiB, 1.66x), Reference grid fits in 30 MB L3, Sweep C: -fno-tree-vectorize null result

### Community 43 - "Slide Stills Grabber"
Cohesion: 0.40
Nodes (5): grab_frame(), main(), Path, Extract one representative still frame from each presentation animation.  The Be, Return the frame at `fraction` of the clip, or the last one if the     clip turn

### Community 44 - "SIMD Masks and Blend"
Cohesion: 0.40
Nodes (5): Blend intrinsic, Comparison intrinsics producing bit masks (_mm256_cmp_ps), Bitwise logical ops (and/or/xor/andnot), Branchless if/else via mask and/andnot/or, Ordered/unordered and signalling NaN comparison flags

### Community 45 - "OpenMP Variable Scoping"
Cohesion: 0.40
Nodes (5): const variable sharing compatibility (OpenMP 4.0 / GCC 9), default(none) clause, shared clause, threadprivate, Variable scoping (shared by default)

### Community 46 - "CUDA Memory Basics"
Cohesion: 0.50
Nodes (5): CUDA_CHECK error macro (do-while(false), __FILE__/__LINE__), cudaMalloc / cudaMemcpy / cudaFree, __global__ kernel (void return, output via pointer), Serial integral-under-curve baseline (-O3, ~900 us), Single-thread kernel is far slower than CPU

### Community 48 - "Double Buffer Time March"
Cohesion: 0.50
Nodes (4): Backward time-marching from expiry payoff, Two-buffer swap (read later sheet, write earlier sheet), Greeks (delta, gamma, vega), Flat row-major grid storage (Grid::index)

### Community 49 - "SIMD Arithmetic Costs"
Cohesion: 0.50
Nodes (4): Division is expensive, prefer multiply, Instruction latency vs throughput (CPI), Vector sqrt and fast approximate rsqrt, Vector add/sub/mul (_mm256_add_ps etc.)

### Community 50 - "OpenMP Parallel Plan"
Cohesion: 0.67
Nodes (4): Two buffers with pointer swap per step, Milestone 2 OpenMP on outer j loop, Amdahl's and Gustafson's laws, Parallel within a step, sequential across steps

### Community 51 - "Rangpur Job Limits"
Cohesion: 0.67
Nodes (3): Derived (not measured) full reference solve runtime, cosc3500 QOS 15-minute wall limit, Rangpur benchmark loop (make sync, remote-bench, queue, fetch)

### Community 52 - "Fused Boundary Control"
Cohesion: 0.67
Nodes (3): v=0 boundary as transport equation with one-sided slope, Fused boundary-row negative control (branch predictor absorbs it), Level 4 loop splitting (v=0 row peeled)

## Ambiguous Edges - Review These
- `Benchmark protocol (rangpur compute node, pinned, 5 sweeps)` → `Reference grid fits in 30 MB L3`  [AMBIGUOUS]
  STUDY_GUIDE.md · relation: conceptually_related_to
- `Convergence order ~2.2 then error floor ~$0.014` → `WALKTHROUGH.md (plain-English program tour)`  [AMBIGUOUS]
  docs/WALKTHROUGH.md · relation: conceptually_related_to

## Knowledge Gaps
- **151 isolated node(s):** `cell`, `weight`, `stock_position`, `variance_position`, `next_` (+146 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **16 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `Benchmark protocol (rangpur compute node, pinned, 5 sweeps)` and `Reference grid fits in 30 MB L3`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **What is the exact relationship between `Convergence order ~2.2 then error floor ~$0.014` and `WALKTHROUGH.md (plain-English program tour)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `Config` connect `Grid Class and Stencil Core` to `Grid Spec and Axis Lookup`, `CLI Entry and Params`, `Solver Class Hierarchy`?**
  _High betweenness centrality (0.018) - this node is a cross-community bridge._
- **Why does `SolveResult` connect `Grid Class and Stencil Core` to `Solver Class Hierarchy`?**
  _High betweenness centrality (0.005) - this node is a cross-community bridge._
- **Why does `Solver` connect `Solver Class Hierarchy` to `Grid Class and Stencil Core`, `CLI Entry and Params`?**
  _High betweenness centrality (0.004) - this node is a cross-community bridge._
- **Are the 13 inferred relationships involving `OptSolver::solve()` (e.g. with `current_` and `.index()`) actually correct?**
  _`OptSolver::solve()` has 13 INFERRED edges - model-reasoned connections that need verification._
- **Are the 12 inferred relationships involving `BaselineSolver::solve()` (e.g. with `current_` and `.index()`) actually correct?**
  _`BaselineSolver::solve()` has 12 INFERRED edges - model-reasoned connections that need verification._