# COSC3500 Assignment Master Plan (matrixass/)

Student 48008361. Spec at `docs/COSC3500_Assignment.pdf` (password protected, personalised, never share). The spec names no due date, but GradeBot printed 29 days until Judgement Day on 2026-08-26, which puts the grading run near 2026-09-24. Confirm the exact submission deadline on Blackboard/ECP and record it here.

Implement three matrix multiply functions, benchmarked by the provided GradeBot harness against MKL (CPU, MPI) and cuBLAS (GPU) on chained products A = A·B of random unitary N×N matrices. Goal set 2026-08-26: 7/7/7, the full 21. The original plan called CPU 7 unreachable, but that predated measurement. The planar 4M build measured 1.121 against MKL on a Haswell node at 96 percent of that node's no-FMA ceiling, and the graded Zen 2 nodes dual-issue separate mul and add at the same 32 flops per cycle as MKL's FMA path, so parity is a real target there. CPU near 1.0 in turn puts MPI 0.6 in range (half the work per node plus a pipelined allgather), and GPU 1.5 needs the full register-tile and transfer-overlap ladder. Every rung stays correctness-gated.

**Iron rule. An implementation whose GradeBot error exceeds ERR_TOLERANCE is never benchmarked, never tuned, and its numbers are never recorded. Correctness gates every rung.**

## Marking (from the personalised spec PDF, ratio = your time / reference time, lower is better)

| Marks | CPU vs 4-thread MKL | GPU vs cuBLAS | MPI (2 nodes × 4 cores) vs single-node MKL |
|---|---|---|---|
| 7 | 1.0 | 1.5 | 0.6 |
| 6 | 1.5 | 2.0 | 1.0 |
| 5 | 3.0 | 3.0 | 1.5 |
| 4 | 6.0 | 4.0 | 3.0 |
| 3 | 12.0 | 5.0 | 6.0 |
| 2 | runs, but more than twice the 3-mark threshold, no timeout | same | same |
| 1 | compiles and runs, wrong answer | same | same |
| 0 | no compile, no submission, or timeout | same | same |

## Ground truth

- Grader compiles with their Makefile. g++ `-std=c++11 -O2 -mavx -fopenmp` (AVX1 only, no FMA, no AVX2). nvcc `-O2 --gpu-architecture=sm_35`, run on A100 via PTX JIT (no tensor cores, no cp.async, 48KB usable shared memory per block, `__ldg` and warp shuffle available). JIT happens once per process and the discarded warm-up rep absorbs it.
- `rubric.h` sets `FLOATTYPE_COMPLEX64`. floatType is `std::complex<float>`, cuComplex on GPU, MKL_Complex8 as reference. ERR_TOLERANCE 1e-7 as the sum of |yours − ref|² over all N² elements. ALIGN 64, every harness CPU buffer is 64-byte aligned. rubric.h overrides any transcript on conflict.
- Graded run is `goslurm_COSC3500Assignment_RangpurJudgementDay`. N=2048, 4 threads, 2 nodes (a100-a, a100-b), one full A100 per node, 15-minute wall, all three benchmarks. CPU and GPU run once per node and the faster node counts.
- Harness behaviour, from librubric.so. It allocates a 256MiB pool of inputs, so 8 reps at N=2048 (first rep discarded as warm-up), a different A each rep, and it ping-pongs the C buffer. Cross-rep result caching is defeated by design. Static scratch and device buffers inside your functions are fine and persist across the 8 reps.
- Harness contract. `./Assignment1_GradeBot N threads runCPU runGPU runMPI [tuning ints...]`, trailing ints arrive as `args`/`argCount`. On the graded run there are no trailing ints, so defaults must be optimal and every access bounds-checked. The student ID is harvested by calling each function with N=0 and null pointers.
- Error budget context. 1e-7 summed over 4.2M elements is about 1.5e-7 RMS per element, only a few times above intrinsic float rounding noise for a K=2048 dot product. Summation order is part of correctness. Blocked or tree accumulation buys headroom, the 3M trick spends it.

## Results log convention

Every sbatch result is appended to `RESULTS.md` with date, phase and rung, git hash, job ID, node hostname, script name, reference and own matrices/second, ratio, and error. Any number quoted in any doc must trace to a `slurm-*.out`.

## Phases

Each phase is sized for one working session including queue waits. Do not start a phase before its entry condition holds.

### Phase 0. Pipeline and reconnaissance

Needs nothing.
Restore the three stubs from git, port `.env`, write `.gitignore`, `Makefile.sync`, and the four new slurm scripts (bench_cpu_2048, bench_gpu_2048, bench_mpi_2node, sweep_args). Push and submit the provided DebugCPU script with the blank stubs.
Done when the DebugCPU output is fetched locally showing the harness ran end to end. A wrong-answer 1-mark result is the expected pass. Record the node CPU model and the MKL reference throughput at N=128 in RESULTS.md, then recalibrate every ratio estimate in this document against the measured numbers. That recalibration is the highest-information ten minutes of the project.
Validated by DebugCPU.

### Phase 1. CPU correct baseline (rungs R1 and R2)

Needs Phase 0.
Replace the stub body with manual float-pair arithmetic (reinterpret_cast the complex buffers to `const float*`, compute re/im by hand, never call `std::complex` operator*), i-k-j loop order so B streams by rows, then `#pragma omp parallel for` over i with per-row accumulation.
Done when DebugCPU (N=128) passes tolerance, then one bench_cpu_2048 run is logged. Expect a ratio around 15 to 25. Never run the pre-R1 naive at N=2048, eight reps of `__mulsc3` ijk can blow the wall clock.
Validated by DebugCPU, then bench_cpu_2048.

### Phase 2. CPU planar rebuild (rungs R3 and R4)

Needs Phase 1 and a measured MKL time.
Allocate static `_mm_malloc` scratch planes (64-byte aligned, allocated on first call for the given N). Split interleaved A and B into separate real and imaginary planes, OpenMP-parallel so first-touch spreads pages across cores. Compute C as four real GEMMs (the 4M scheme) with cache blocking, starting points Kc=256, Mc=96, packed B panels. Write a hand AVX1 16×4 microkernel, two ymm rows by four broadcast columns, eight ymm accumulators, `_mm256_broadcast_ss` for B, dual-issued `_mm256_mul_ps` and `_mm256_add_ps`. Merge the planar result back to interleaved C.
Done when DebugCPU passes and bench_cpu_2048 lands at ratio 3.0 or better (the 5-mark line). Record achieved GFLOP/s against the AVX1 ceiling of 16 flops per cycle per core.
Validated by DebugCPU, then bench_cpu_2048.

### Phase 3. CPU 3M and tuning sweep (rung R5)

Needs Phase 2 and its recorded error margin.
Wire Kc/Mc/Nr into `args[0..2]` with graded-safe defaults and run sweep_args. Implement the 3M Gauss trick behind an args-selectable switch and measure its error. Ship 3M only if the error stays at or below ERR_TOLERANCE/3. Lock the best settings into compiled defaults.
Done when the best ratio is recorded (1.0 is the target, anything at or under 1.2 keeps MPI 0.6 alive) and a bench_cpu_2048 run with no extra args reproduces it. Ratios must be read on EPYC nodes (vcpu-* pins the graded silicon), Haswell r730 numbers flatter us.
Validated by sweep_args, then bench_cpu_2048.

### Phase 4. GPU correct baseline (rungs G0 and G1)

Needs Phase 1, the CPU correctness work gives a trusted reference mental model.
Write the host wrapper first. Static cached `cudaMalloc` buffers sized on first call, `cudaMemcpy` in, kernel launch, copy back, `cudaDeviceSynchronize`, every CUDA call return-checked during development (checks removed before submission along with any prints). Implement the declared `matrixMultiplyKernel_GPU` as a naive one-thread-per-element kernel first (thread x maps to column for coalescing), and keep that exact kernel defined forever since the .cuh declares it. Then add a 32×32 shared-memory tiled kernel (16KB per block) as a separate `__global__` function.
Done when DebugGPU passes tolerance at N=2048 and the tiled ratio is recorded from bench_gpu_2048. Expect around 4 to 7. DebugGPU runs on a shared GPU slice, treat its timing as smoke only.
Validated by DebugGPU, then bench_gpu_2048.

### Phase 5. GPU register blocking and transfer pipeline (rungs G2 to G4)

Needs Phase 4.
Build the 64×64 block tile kernel, 16×16 threads each computing a 4×4 cuComplex register tile, Kc=16 double-buffered shared memory at 32KB per block (inside the sm_35 48KB model). Use manual complex FMA expressions (a.x*b.x − a.y*b.y forms), never cuCmulf in the hot loop. Add float4 loads (two cuComplex per load, safe on the 64-byte aligned harness buffers) and `__ldg` with `const __restrict__` pointers. Then the memory pipeline, cache `cudaHostRegister` of the harness buffers behind a static flag, two streams, chunk A by row blocks to overlap H2D with compute and D2H of finished rows.
Done when DebugGPU passes and bench_gpu_2048 lands at 1.5, the 7-mark line. G2 and G3 alone are unlikely to get there, so G4 transfer overlap is on the critical path unless the reference turns out to absorb transfers. Record whether the cuBLAS reference time visibly includes transfers, compare its implied ms against the roughly 7ms compute-only model. That measurement decides where the remaining gap to 1.5 must come from.
Validated by DebugGPU, then bench_gpu_2048.

### Phase 6. MPI

Needs Phase 2, the MPI ratio is dominated by CPU quality.
Rank p computes C rows [r0, r1) where rows(p) = N/P + (p < N%P) and r0 = p*(N/P) + min(p, N%P), directly into `&C[r0*N]`, using a row-range helper. Define that helper non-static in matrixMultiply.cpp and forward-declare it extern in matrixMultiplyMPI.cpp, both object files link into the one GradeBot binary (verified against the grader's Makefile). Reassemble with `MPI_Allgatherv` using `MPI_IN_PLACE`, datatype `MPI_FLOAT` with all counts doubled (avoids every C++ complex-datatype portability question). The N<=0 guard runs before any MPI call. No MPI_Init or MPI_Finalize, the harness owns them. OpenMP still runs 4 threads inside each rank.
Done when DebugMPI (2 ranks, 1 node, N=128) passes, then bench_mpi_2node passes at N=2048 with ratio 0.6, the 7-mark line. The math is cpu_ratio/2 plus the allgather cost (each rank moves about 17MB), so 0.6 needs the CPU near 1.0 and the allgather hidden. Plan the hiding from the start: chunk the row range, compute chunk c while MPI_Iallgatherv of chunk c-1 is in flight, MPI_Waitall at the end. If the interconnect turns out to be InfiniBand the plain blocking allgather may already fit inside 0.6, measure before adding the pipelining complexity.
Validated by DebugMPI, then bench_mpi_2node.

### Phase 7. Dress rehearsal

Needs all three implementations passing tolerance at graded size on the bench scripts.
Run JudgementDay once, its first permitted use. Compare all three ratios against the bench numbers, the pinned a100 pair may differ from whatever free nodes the benches ran on. Triage the weakest ladder, apply at most one targeted fix, and confirm with at most one more JudgementDay run.
Done when a JudgementDay output shows all three correct with ratios recorded and no planned optimisation remains.
Validated by JudgementDay, budget two runs this phase.

### Phase 8. Freeze and package

Needs Phase 7.
Strip every print and debug artifact from all three files. Re-verify the immutables, signatures byte-identical to the headers, STUDENTID returns, N<=0 guards, no includes beyond what the provided headers pull in (`_mm_malloc` comes via immintrin.h, memset/memcpy via string.h, nothing extra is needed). Fresh rsync, `make clean && make all` on rangpur, one final JudgementDay confirmation. On rangpur, `zip slurm.zip slurm-*.out` then `zip 48008361.zip matrixMultiply.cpp matrixMultiplyGPU.cu matrixMultiplyMPI.cpp slurm.zip`. Fetch the zip and verify by listing, it must contain exactly those four items.
Done when submitted.
Validated by the final JudgementDay and a manual zip listing.

## Appendix A. CPU ladder

Work per multiply is 8N³ real flops for 4M (68.7 GFLOP at N=2048) and about 6N³ for 3M (51.5 GFLOP). The AVX1 no-FMA ceiling is 16 SP flops per cycle per core, roughly 190 GFLOP/s on 4 cores at 3GHz. Ratio estimates assume MKL cgemm at 250 to 350 GFLOP/s on an AVX2+FMA path. Phase 0 recalibrates all of this against the measured node.

| Rung | Change | Expected ratio | Risk |
|---|---|---|---|
| R0 | naive ijk with `std::complex` ops | 100 to 300+ (`__mulsc3`) | timeout at 2048, never bench |
| R1 | manual float-pair arithmetic, ikj order | 50 to 80 | none, correctness anchor |
| R2 | OpenMP over i, 4 threads | 15 to 25 | false sharing if C is not row-split |
| R3 | cache blocking on interleaved data | 6 to 10 | tile-size sensitivity |
| R4 | planar split, 4M real GEMMs, packed 16×4 AVX1 microkernel | 2.2 to 3.0 | register pressure, packing bugs |
| R5 | 3M Gauss trick plus swept Kc/Mc | 1.7 to 2.3 | error tolerance, ship only with 3× margin |

R4 is needed for 5 marks, R5 pushes toward 6. The interleaved-AVX shuffle scheme (moveldup/movehdup/addsub) stays documented as a fallback only, it caps near 60% of peak because it saturates the shuffle port.

## Appendix B. GPU ladder

cuBLAS cgemm on an A100 runs near 10 to 14 TFLOP/s effective, about 5 to 7ms compute at N=2048. Transfers are 96MiB, about 4 to 5ms at PCIe 25GB/s. Whether the reference timing includes transfers is unknown until Phase 4/5 measures it, and that answer gates rung G4.

| Rung | Change | Expected ratio | Risk |
|---|---|---|---|
| G0 | naive thread-per-element, coalesced | 15 to 50 | none, correctness anchor |
| G1 | 32×32 shared-memory tiles | 4 to 7 | bank conflicts with float2, pad or swizzle |
| G2 | 64×64 block, 4×4 register tile, double-buffered 32KB smem | 2 to 3 | register spills past 255 per thread |
| G3 | float4 loads, manual complex FMA, `__ldg` | 1.5 to 2.2 | alignment assumptions (safe, buffers are 64B aligned) |
| G4 | cached device buffers, pinned host memory, 2-stream overlap | 1.3 to 1.8 | only pays if reference absorbs transfers differently |

Device-side caching of B contents keyed on pointer identity is a guarded stretch goal, off by default, it saves about 1.3ms of 10 and risks correctness if the harness swaps buffers.

## Appendix C. MPI design

- Row split as in Phase 6, recvcounts[p] = rows(p)*N*2 floats with matching displacements.
- `MPI_FLOAT` with doubled counts is deliberate. `MPI_COMPLEX` is a Fortran type and `MPI_C_FLOAT_COMPLEX` is C99 `_Complex`, both occasionally misbehave in C++ builds. Raw floats cannot break.
- Optionally guard with `MPI_Initialized` and fall back to a plain serial call, cheap insurance.
- The reference is single-node 4-thread MKL, about 0.2s per multiply at 300 GFLOP/s. Allgather of 17MB per rank costs about 25 to 30ms on 10GbE (0.10 to 0.15 of the reference) or near nothing on InfiniBand. With CPU ratio 1.8 to 2.2 the MPI ratio lands near 1.0 to 1.25, which is 5 to 6 marks.
- A GPU-inside-MPI implementation would crush the 0.6 line since JudgementDay allocates a GPU per node during the MPI test, but the rubric words the MPI column as 2 nodes × 4 cores. Check the spec wording and, if needed, the course forum before writing a line of it. Default is CPU-only.

## Appendix D. Open uncertainties, resolved by measurement

1. Host CPU model and ISA of the graded nodes, hence MKL's real GFLOP/s and every CPU ratio estimate. Phase 0, lscpu plus the measured MKL reference line.
2. Whether the cuBLAS reference includes H2D/D2H. Resolved 2026-08-26 by run 568780. Fitting a·N² + b·N³ to register-kernel times at N 1024, 2048 and 4096 gives an N² share of about 0.7 ms at N=2048, far below what 100 MB of PCIe copies would cost, so the harness GPU buffers are effectively device-resident (managed or mapped inside librubric) and transfers are nearly free for both sides. G4 is moot, the graded ratio is decided by kernel arithmetic alone.
3. Interconnect bandwidth between the graded node pair. Resolved 2026-08-26 by run 568790 on a100-a and a100-b. Two-node time 0.1487 s against half the single-node 0.240 s says the 16.8 MB allgather costs single-digit milliseconds. Blocking MPI_Allgatherv suffices, no pipelining needed.
4. 3M error headroom under the 1e-7 budget. Phase 3, ship gate.
5. Legality of GPU use inside the MPI benchmark. Resolved 2026-08-26 by reading the spec. No explicit prohibition, but the rubric column says 2 nodes 4 cores each, and nothing guarantees the MPI job allocation exposes a GPU, so a GPU-dependent MPI risks a crash and 1 mark. Decision: CPU-only MPI. Revisit only if the CPU-only path measures above 0.6 and the course forum blesses GPU use.
6. Whether B is bit-identical across the 8 reps. Only matters for the guarded B-caching stretch.

## Appendix E. Transcript index (docs/)

- `COSC3500_Assignment_-transcript.txt` covers the briefing, the live rangpur demo, debug-vs-JudgementDay etiquette, and the note that MPI always initialises even in single-node runs (startup is the slowest part of a job).
- `COSC3500_T04_OMP.mp4-transcript.txt` covers parallel regions, work sharing, schedule/collapse, data-sharing clauses, race fixes ranked (critical, atomic, per-thread partials, reduction).
- `COSC3500_T05_MPI.mp4-transcript.txt` covers init/rank/size, send/recv, Bcast/Reduce/Allreduce/Allgather, domain decomposition with remainder handling, cluster workflow.
- `COSC3500_T06_CUDA.mp-transcript.txt` covers kernels, cudaMalloc/Memcpy, error-check macro, indexing and grid-stride loops, `__syncthreads` rules, shared memory and bank conflicts, reduction addressing schemes, vectorized loads, 128-bit per-thread fetch width.
