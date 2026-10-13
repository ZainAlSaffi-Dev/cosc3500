# Changelog

This document records every change made to the assignment, in order, with an explanation of what was done and why. It is written so you can read it top to bottom and understand each decision well enough to explain it yourself.

## Phase 0. Pipeline and reconnaissance (already complete before this session)

Phase 0 built the plumbing that every later phase depends on. Nothing in it touches the maths yet.

### What exists and why

The three submission files, `matrixMultiply.cpp`, `matrixMultiplyGPU.cu`, and `matrixMultiplyMPI.cpp`, were restored as blank stubs from the course-provided originals. Everything else in the repo belongs to the grader and is never edited. That includes the headers in `include/`, the `Makefile`, `librubric.so`, and the provided `goslurm_*` scripts.

`.gitignore` keeps secrets and build junk out of git. The most important line excludes `.env`, which holds the cluster password. A leaked password cannot be unleaked, so it is excluded from git, from rsync, and from every script.

`Makefile.sync` exists because the assignment `Makefile` only builds on rangpur. It references library paths like `/home/groups/cosc3500/shared/matmul/` that do not exist on the Mac. So the workflow is to edit locally and run remotely. `Makefile.sync` gives short commands for each step of that loop. `push` rsyncs the working tree up, `run-<script>` submits a slurm script, `queue` shows the job queue, and `fetch` pulls the output files back down. It also prepends `/opt/slurm/bin` to PATH in every remote command, because non-interactive ssh sessions do not have `sbatch` on PATH.

Four new slurm scripts were written by copying the partition names and module loads from the provided debug scripts. `bench_cpu_2048.sh`, `bench_gpu_2048.sh`, and `bench_mpi_2node.sh` run the graded problem size N=2048 on whatever node is free, because the five-minute debug scripts only run N=128 and small-N timings do not predict large-N timings. `sweep_args.sh` will later try different tuning values. All of them stay under the 15-minute wall-time cap that the course queue enforces.

### The smoke test

With the stubs still blank, the provided DebugCPU script was submitted as job 568325. The point was not to get a right answer. The point was to prove that the whole loop works end to end, meaning the code syncs up, compiles with the grader's exact flags, runs under slurm, and the output comes back down. A blank stub scores 1 mark, compiles and runs but wrong answer, and that is exactly what the output showed. The run also recorded two facts in `RESULTS.md`. The debug nodes are AMD EPYC 7542 machines, and the MKL reference does about 10117 multiplies per second at N=128. Those measured numbers replace guesses in all later planning.

## Phase 1. CPU correct baseline (this session)

### What changed

One file changed, `matrixMultiply.cpp`. The line `//WRITE YOUR CODE HERE` was replaced with about thirty lines of code. Every other byte of the file was left exactly as provided, including the trailing whitespace on some lines, because the assignment is personalised and stray formatting could be part of how the grader identifies files.

### The problem being solved

The function must compute C = A times B for two N by N matrices of complex numbers, where each complex number is a pair of 32-bit floats. The textbook definition says element (i, j) of C is the sum over k of A(i, k) times B(k, j). The matrices are stored in row-major order, meaning row 0 comes first in memory, then row 1, and so on. Each complex number is stored as two floats sitting next to each other, real part then imaginary part. So element (i, j) of a matrix lives at float positions 2*(i*N+j) and 2*(i*N+j)+1.

### Decision 1. Do the complex arithmetic by hand

The obvious way to multiply two `std::complex<float>` values is to write `a * b` and let C++ do it. That is a trap on this assignment. Under the grader's compile flags, that innocent-looking `*` becomes a call to a library function named `__mulsc3`. That function exists to handle strange edge cases involving infinities and NaN correctly, and it does so with branches and extra work on every single multiply. Called billions of times in the inner loop, it is catastrophically slow, and it also stops the compiler from vectorising anything.

The fix is to view the buffers as plain arrays of floats using `reinterpret_cast`, and write the complex multiply formula out manually. For two complex numbers (ar + ai*i) and (br + bi*i), the product is (ar*br - ai*bi) + (ar*bi + ai*br)*i. Four multiplies and two adds, no function call, no branches. This is safe here because we know the inputs are ordinary finite numbers, so the edge-case handling buys nothing.

### Decision 2. Loop order i-k-j instead of i-j-k

The textbook loop order is i, then j, then k, computing one element of C completely before moving on. The problem is what the innermost loop touches in memory. With k innermost, the loop reads B(k, j) for increasing k, which means jumping a whole row of B forward on every step. Memory is fetched in 64-byte cache lines, so each jump drags in 64 bytes and uses 8 of them. That wastes most of the memory bandwidth.

Reordering to i, then k, then j fixes it. Now the innermost loop holds A(i, k) fixed in a register and walks along row k of B and row i of C, both of which are contiguous in memory. Every byte of every cache line gets used. The result is mathematically identical, the sum over k still happens for every (i, j), just interleaved differently.

### Decision 3. Zero each row of C before accumulating

The i-k-j order accumulates into C gradually rather than computing each element in one go, so C must start at zero. The harness reuses the C buffer between repetitions and it arrives full of the previous answer, so the code cannot assume it is clean. Each row is zeroed with `memset` right before the loop that fills it.

### Decision 4. Parallelise over rows with OpenMP

The line `#pragma omp parallel for schedule(static)` above the i loop tells OpenMP to split the rows of C across the 4 threads the harness provides. Rows are the right unit of work because they are completely independent. Thread 0 writes rows 0 to 511, thread 1 writes rows 512 to 1023, and so on, and no two threads ever touch the same memory.

That independence matters for a subtle reason called false sharing. If two threads write to different variables that happen to sit in the same 64-byte cache line, the hardware forces the line to bounce between the two cores on every write, which is slow. Splitting by whole rows makes each thread's output region contiguous and separate, so the problem cannot occur.

`schedule(static)` means the rows are divided into equal chunks up front. That is the right choice when every row costs the same amount of work, which is true here, and it avoids the coordination overhead of dynamic scheduling.

### Decision 5. Swap the operands, because the harness is column-major

The first submitted version computed C = A times B in the obvious row-major reading and failed the correctness check badly, with an error of 2.6e+02 against a tolerance of 1e-7 (job 568373). An error that large is never rounding noise. Rounding error at N=128 would be around 1e-9. So the computed answer was a genuinely different matrix, which means the reference was computing something different from what the code assumed.

The investigation went like this. The grading library `librubric.so` ships with the assignment, and the `nm` tool lists which external functions a compiled library calls. It showed a call to `cgemm`, the Fortran-style MKL matrix multiply, rather than `cblas_cgemm`, the C-style one. That matters because Fortran stores matrices column by column while C code normally stores them row by row. Disassembling the library confirmed it. The reference function passes 'N' for both transpose flags, meaning a plain column-major product with no transposition.

Here is the consequence. The same buffer of numbers read column-major as matrix M looks, when read row-major, like the transpose of M. The reference computes C = A times B in column-major terms. Taking the transpose of both sides gives C transposed equals B transposed times A transposed. In row-major eyes those transposes are exactly what we see in the buffers. So to produce the bytes the grader expects, row-major code must multiply in the order B times A.

The fix was two lines. The left operand of the loop is now B and the right operand is A. The loop structure did not change at all, because which pointer plays which role is independent of the traversal order. Job 568381 validated the swap.

The lesson worth keeping is that a huge error and a tiny error mean different things. A tiny error points at accumulation order or precision. A huge error points at a layout, ordering, or indexing mismatch, and the first suspect between two libraries is always row-major versus column-major.

### What was deliberately not done yet

No SIMD intrinsics, no cache blocking, no planar data layout. Those are Phase 2. Phase 1 exists to be simple enough to trust, because every later optimisation is checked against it. The plan predicts a ratio around 15 to 25 versus MKL, which is far from the target, and that is fine. Correctness comes first and speed comes after.

### Validation

Three cluster runs tell the Phase 1 story, all logged in `RESULTS.md`.

Job 568373 ran the first version at N=128 and failed with error 2.6e+02. That failure led to the column-major discovery described in Decision 5.

Job 568381 ran the fixed version at N=128 and passed with error 2.3e-08, well under the 1e-7 tolerance. Correctness gate cleared.

Job 568383 ran the graded size N=2048 and passed with error 1.2e-08 at a ratio of 8.203 versus MKL, which the harness scores as about 3.5 marks. The plan predicted 15 to 25 for this rung, so the baseline landed ahead of schedule, probably because the compiler auto-vectorised the contiguous inner loop. One caution applies. That run landed on r730-0, an old Intel Haswell node, and the graded nodes are faster AMD EPYC machines where MKL will do better, so the graded ratio will likely be worse than 8.2.

### Where this leaves the project

Phase 1 is complete. There is a correct, moderately fast CPU baseline that every later optimisation will be checked against. Next is Phase 2, which rebuilds the data layout so the code can use the AVX vector units properly, aiming for a ratio of 3.0 or better, which is the 5-mark line.

## Phase 2. CPU planar rebuild (this session)

Again one file changed, `matrixMultiply.cpp`. The Phase 1 loops were replaced by a restructured implementation of about 250 lines. The column-major lesson from Phase 1 carries over unchanged, the code still computes B times A in row-major terms.

### The problem with the Phase 1 layout

The AVX vector unit does 8 float multiplies in one instruction, but only if the 8 floats sit next to each other in memory. Complex numbers are stored interleaved, real, imaginary, real, imaginary. Loading 8 consecutive floats therefore grabs 4 real parts and 4 imaginary parts mixed together, and the arithmetic for real and imaginary parts is different. Making that work requires shuffle instructions to sort the values, and the shuffle unit becomes the bottleneck. The plan caps that approach at about 60 percent of peak, so Phase 2 changes the data layout instead of fighting it.

### Decision 1. Split into planes, the planar layout

Before any arithmetic, the code copies each input matrix into two separate float arrays, one holding all the real parts and one holding all the imaginary parts. These are called planes. Inside a plane, 8 consecutive floats are 8 real parts of 8 neighbouring matrix elements, exactly what the vector unit wants. The copying costs one pass over the data, which is nothing next to the N cubed arithmetic.

The scratch planes are allocated once with `_mm_malloc`, which is like `malloc` but lets you demand 64-byte alignment, and kept in static variables so the harness's 8 repetitions reuse them. The split loop is parallelised over rows for a reason beyond speed. On a machine with two CPU sockets, memory pages live near whichever core first writes them, the first-touch rule. Splitting in parallel with the same row partition the compute loop uses later means each thread computes on memory that lives near its own socket.

### Decision 2. The 4M scheme

For complex values, C equals A times B expands into four real matrix products. The real part of C is Are times Bre minus Aim times Bim, and the imaginary part is Are times Bim plus Aim times Bre. So one complex multiply becomes four real multiplies on the planes, hence the name 4M. Real-only inner loops are exactly the shape the vector unit handles at full speed.

The plan describes this as four separate full passes. The code fuses them two at a time instead. The kernel that produces the real-part plane computes both of its products in one pass over the data, and likewise for the imaginary-part plane. The arithmetic is identical, but each output tile is loaded and stored once instead of twice, halving that traffic. The fused kernel needs 14 of the 16 vector registers, so it still fits.

### Decision 3. The microkernel

The innermost code computes a small tile of the output, 4 rows by 16 columns, held entirely in 8 ymm registers, the 256-bit AVX registers that each hold 8 floats. For every k it broadcasts one value from the left matrix across a register with `_mm256_broadcast_ss`, multiplies it against 16 contiguous values from the right matrix, and adds into the accumulators. Keeping the accumulators in registers for the whole k loop matters because the alternative, reading and writing C from memory every step, is exactly what made Phase 1 slow relative to peak.

The tile shape is dictated by the hardware. Wider would need more registers than exist. Narrower would not have enough independent additions in flight to cover the 3-cycle latency of the add instruction. The compile flags allow no fused multiply-add, so the kernel issues separate multiply and add instructions, which the CPU can run in the same cycle on different ports. Sixteen floats of multiply plus sixteen of add per cycle is the machine's ceiling, and this kernel shape can reach it.

The minus sign in the 4M formula is handled by a template parameter, a compile-time switch that produces one kernel variant that adds the second product and one that subtracts it, so the inner loop contains no branch.

### Decision 4. Cache blocking and packed panels

A 2048 by 2048 plane is 16MB, far bigger than any cache, so the loops are blocked. The k dimension is cut into slabs of Kc equal to 256. Within a slab, the columns of the right matrix are repacked into strips of 16, so the kernel reads its 64 bytes per step from one contiguous, aligned line instead of striding across a 16MB plane. Rows of C are processed in bands of Mc equal to 96, sized so a band's slice of the left planes stays resident in the 256KB L2 cache while every packed strip streams past it. Both numbers are the plan's starting points and get swept for better values in Phase 3.

Blocking also quietly improves accuracy. Each element of C is now built as 8 partial sums of 256 terms instead of one chain of 2048, and shorter chains accumulate less float rounding error.

### Decision 5. Pad instead of writing edge cases

A 4 by 16 kernel only works if the matrix dimensions divide by 4 and 16. Instead of writing special narrower kernels for the edges, the scratch planes are padded up to the next multiple and the padding is zeroed once at allocation. The kernel then runs full tiles everywhere. The padded rows and columns compute zeros that the final merge step simply never reads. The costs are a sliver of wasted arithmetic and memory, and the benefit is that one kernel handles every N.

### Decision 6. Test locally before spending cluster runs

Apple silicon Macs run x86 programs through Rosetta, which now supports AVX. So the exact AVX code that runs on the cluster was compiled for x86 on the Mac and checked against a naive reference computed in double precision, at sizes 1, 3, 5, 16, 33, 100, 128, and 257, the odd ones chosen to exercise the padding edges, each size run twice to exercise the static scratch reuse. All sixteen checks passed with errors at rounding-noise level. This catches logic bugs in seconds instead of burning five-minute cluster queue waits on them.

### Validation

Job 568422 ran DebugCPU at N=128 and passed with error 3.3e-08, ratio 1.506.

Job 568424 ran the graded size N=2048 and passed with error 2.3e-08 at a ratio of 1.121 versus MKL, which the harness scores as 6.7 marks. The phase gate was 3.0 and the plan's expected range for this rung was 2.2 to 3.0, so the rebuild landed well ahead of both. In absolute terms the code sustains about 141 GFLOP/s, roughly 96 percent of what the vector units can theoretically do at base clock without fused multiply-add, which means there is very little left to win from the arithmetic itself on this node.

The usual caution applies. Both runs landed on r730-0, the old Haswell node. On the graded AMD nodes MKL gets faster, but those cores can also issue our separate multiplies and adds at double the rate, so the graded ratio has to be measured, not assumed.

### Where this leaves the project

Phase 2 is complete and the CPU ladder already sits past the 6-mark line on the node it was measured on. Phase 3 wires the blocking sizes to the harness's tuning arguments, sweeps them, and decides whether the 3M trick is worth its precision cost. Given how close to the ceiling the kernel already runs, the sweep is likely to matter less than confirming the ratio on an EPYC node.

## Phase 3. CPU tuning, the 3M scheme, and the new goal (this session)

The goal changed to full marks on all three ladders, 7/7/7. The plan originally wrote off CPU 7 as unreachable, but that judgement predated any measurement. Two facts changed the picture. The planar rebuild measured 1.121 on a Haswell node while running at 96 percent of that node's ceiling, and the graded AMD nodes can issue our separate multiply and add pairs at 32 flops per cycle, the same rate as MKL's fused path, so on the graded silicon the ceiling argument no longer favours MKL. MASTERDOC now records the new targets, CPU 1.0, GPU 1.5, MPI 0.6.

### Change 1. Comment style

All comments in `matrixMultiply.cpp` were rewritten to be short and plain, with no references to phases or planning documents. Explanations live here instead.

### Change 2. Tuning arguments

The harness can pass optional integers on the command line that arrive in the function as `args`. The code now reads args[0] as the panel depth Kc, args[1] as the row band Mc, and args[2] equal to 1 as the switch for the 3M scheme. Every access is bounds checked and range checked, and when the arguments are absent, which is exactly what happens on the graded run, the compiled-in defaults are used. This exists so a single cluster job can benchmark many tuning combinations without recompiling. Mc is rounded down to a multiple of 4 because the kernel processes 4 rows at a time, and an unaligned band boundary would make two tasks compute the same tile.

### Change 3. The 3M scheme

The 4M formula uses four real products to build one complex product. A rearrangement attributed to Gauss gets away with three. Compute T1 as Lre times Rre, T2 as Lim times Rim, and T3 as the product of the two elementwise sums, (Lre plus Lim) times (Rre plus Rim). Then the real part is T1 minus T2, as before, and the imaginary part is T3 minus T1 minus T2, because expanding T3 gives all four cross terms and subtracting T1 and T2 leaves exactly the two the imaginary part needs. Three big products instead of four cuts the dominant arithmetic by a quarter.

The catch is precision. T3 multiplies values up to twice as large, and the double subtraction cancels digits, so the rounding error grows. The tolerance verdict has to come from a measurement at the graded size with the harness's own matrices, so 3M sits behind the args switch, off by default, until that number exists. The plan's rule is to ship it only if the error stays under a third of the tolerance.

The implementation adds a third left plane holding Lre plus Lim, a third packed panel holding the strip sums, and a single-product kernel that computes one T tile in registers and then adds or subtracts it into one or two output tiles. The signs are template parameters again, so each variant compiles to branch-free code.

### Change 4. Measuring on the graded silicon

Every benchmark so far landed on r730 nodes, old Intel Haswell machines, because the queue picks whatever is free. The graded nodes are AMD EPYC 7542 hosts. The partition also contains vcpu nodes with that same EPYC 7542 processor, so a new script `bench_cpu_epyc.sh` pins to vcpu-4 with `--nodelist` and runs the default configuration plus four tuning variants there. That calibrates the ratio on the silicon that decides the grade without touching the reserved a100 pair. The sweep script `sweep_args.sh` was filled in with a 9-point grid over Kc and Mc plus two 3M runs, all inside one 14-minute job.

### Change 5. The spec question about MPI

The plan flagged an open question, whether the MPI benchmark may use the GPU, since a GPU per node would make the MPI target trivial. The spec was read to settle it. Its marking table describes the MPI test as 2 nodes with 4 cores each and never mentions GPU use, and nothing guarantees the MPI job even has a GPU visible, so relying on one risks a crash and 1 mark. The decision recorded in MASTERDOC is a CPU-only MPI. The route to the 0.6 line is instead to hide the communication, each node computes its half of the rows in chunks and the finished chunks are exchanged with non-blocking gathers while later chunks are still computing. The spec also confirms the MPI file may reuse the CPU implementation, which the plan already intended.

### Validation

The sweep (job 568658, Haswell) showed the 4M grid is nearly flat between 1.09 and 1.13, so the exact blocking numbers barely matter on that node. The 3M runs were the story. 3M at Kc=384, Mc=96 hit ratio 0.881, and its error came out at 5.3e-09, six times lower than 4M's 2.3e-08 rather than higher. The likely reason is that MKL itself uses a 3M-style algorithm internally, so our 3M rounds the same way the reference rounds and the measured difference shrinks. Either way the error sits twenty times under tolerance and far under the plan's ship gate of a third of tolerance, so 3M became the default.

The EPYC calibration (job 568659, vcpu-4, the same EPYC 7542 silicon as the graded nodes) delivered the decisive numbers. The old 4M defaults already measured 1.006 there, 4M at Kc=384 hit 0.999, and 3M reached 0.806 with error 5.0e-09. The defaults were locked at Kc=384, Mc=96, 3M on.

Confirmation runs with no arguments reproduced the numbers, 0.884 on Haswell (job 568685) and 0.805 on EPYC (job 568686), both scoring 7.00, and the N=128 correctness gate passed with the new defaults (job 568709, error 3.4e-08). One job, 568684, failed because three concurrent jobs shared one build directory and a `make clean` from one deleted the binary of another between its build and run steps. The lesson recorded in RESULTS.md is to submit one job at a time per build directory.

### Where this leaves the project

The CPU ladder is done at 7 marks with a 24 percent margin on the graded silicon, pending only the dress rehearsal on the actual a100 pair. The 3M scheme, swept blocking, and locked defaults close Phase 3. GPU work was explicitly deferred to a future session. The remaining road to 21 is the GPU ladder to ratio 1.5 and the MPI ladder to 0.6, and the MPI design in MASTERDOC now assumes the CPU near parity, which this session delivered.

## Phase 4. GPU correct baseline (this session)

One file changed, `matrixMultiplyGPU.cu`. The blank stub became a host wrapper plus two kernels, about 120 lines.

### The shape of a GPU multiply

A GPU function has two halves. Host code runs on the CPU and organises things, and kernels run on the GPU itself, thousands of threads at once. The host wrapper here does four jobs. It allocates GPU memory once and keeps it in static pointers, so the harness's 8 repetitions reuse the same buffers instead of paying allocation cost every rep. It copies the two input matrices to the GPU. It launches a kernel. It copies the result back, and because that copy queues up behind the kernel on the same GPU work stream and blocks until finished, it doubles as the synchronisation point.

The column-major lesson from Phase 1 carries over unchanged. The kernels compute a plain row-major product and the host passes the operands swapped, so the bytes that come back match the Fortran-order reference.

### The two kernels

The header declares a kernel called `matrixMultiplyKernel_GPU`, so that name must exist. It became the naive version, one thread per output element, each thread walking the full k loop. Thread numbering maps the fast-moving thread index to the column, so threads sitting next to each other read and write neighbouring memory addresses. The GPU merges such neighbouring accesses into single wide transactions, which is called coalescing, and it is the first rule of GPU memory performance.

The naive kernel reads every input element from GPU main memory over and over, N reads per element. The tiled kernel fixes that with shared memory, a small fast scratchpad that the 1024 threads of one block share. The output is cut into 32 by 32 tiles. The block loads one 32 by 32 tile of each input into shared memory, every thread does 32 multiply-adds out of the scratchpad, and the block steps to the next tile along k. Each global memory element is now loaded once per tile pass instead of 32 times. A barrier called `__syncthreads` sits between loading and computing, because threads must not start reading the scratchpad before every thread has finished filling it.

Complex arithmetic is written out by hand on the float pairs, the same habit as the CPU side. The shared arrays get one extra column of padding, cheap insurance against threads landing on the same shared memory bank.

### Testing without a GPU

The Mac has no NVIDIA GPU, so the kernels cannot run locally. Instead both kernel bodies were transcribed into plain C++ with loops standing in for the thread grid, the barrier modelled by running each phase for all threads of a block before the next phase. Both emulations were checked against a double precision reference at nine sizes chosen to stress the tile edges. Both also had to agree with each other bit for bit, which they must because they add the same numbers in the same order. All passed before anything was submitted to the cluster.

### Validation

Job 568765 (DebugGPU, N=2048, shared GPU slice) passed first try with error 3.2e-08 against the 1e-7 tolerance, ratio 1.627. Job 568767 (full GPU, a100-a) recorded the real number, ratio 1.578, grade 6.82. The plan expected 4 to 7 for this rung, so the tiled kernel landed far ahead of schedule. The reference throughput implied about 15 TFLOP/s for cuBLAS, which raised a puzzle about transfer costs that Phase 5 resolved.

## Phase 5. GPU register tiling (this session)

Again only `matrixMultiplyGPU.cu` changed. A third kernel was added and made the default, with the optional args able to select the older two for comparison runs.

### Why the tiled kernel was still slow

In the 32 by 32 tiled kernel every thread computes one output element, and for each k step it does one complex multiply-add per two shared memory reads. Shared memory is fast but not free, and at that ratio the shared traffic, not the arithmetic, sets the ceiling. The fix is to make each thread compute more outputs per shared read. If a thread holds a 4 by 4 block of outputs in registers, it reads 4 values from each input tile and does 16 complex multiply-adds with them. Four times less shared traffic per unit of arithmetic moves the bottleneck back to the floating point units, which is where a matrix multiply belongs.

### The register tile design

Each block now owns a 64 by 64 tile of C with 16 by 16 threads, and each thread accumulates a 4 by 4 patch of complex results in registers. The k dimension advances in slabs of 16 through shared memory. A thread's 4 rows and 4 columns are not adjacent but strided 16 apart, so the 16 threads of a row always read 16 consecutive shared entries, which spreads them across the memory banks, and the final stores to C stay coalesced for the same reason.

Two more tricks come from the plan's ladder. Global loads go through `float4`, a 16-byte load that fetches two complex numbers in one instruction, legal because the harness aligns its buffers to 64 bytes and the even N keeps every such load aligned. Reads of the inputs also use `__ldg`, a hint that the data is read-only, which routes them through a cache path older GPUs otherwise skip. The shared tiles are double buffered. While the threads compute on the current slab, the loads for the next slab are already in flight into registers, and after the compute they are stored into the second buffer. The writers target a different buffer than the readers, so one barrier per slab is enough. Odd sizes and edges fall back to a guarded scalar path that pads with zeros.

The same emulation harness gained a transcription of this kernel, barriers and double buffering modelled faithfully, and all sizes again matched the naive kernel bit for bit before submission.

### Validation and the transfer puzzle

Job 568779 (DebugGPU) passed the correctness gate, same 3.2e-08 error, as it must since the summation order is unchanged. Job 568780 (full GPU, a100-b) recorded ratio 1.316, grade 7.00, crossing the 7-mark line with a 14 percent margin on the graded silicon.

The same job ran the naive and 32 by 32 kernels through the args switch and the default kernel at N 1024 and 4096. Splitting the times into an N-squared part and an N-cubed part showed the N-squared share is only about 0.7 ms of the 5.9 ms total. Copying 100 MB over PCIe costs 4 ms or more, so the harness's GPU buffers cannot be ordinary CPU memory, they must be managed memory already resident on the device. That resolves the plan's open question about whether the reference absorbs transfers. Transfers are nearly free for everyone, the planned pinned-memory and stream-overlap rung would buy nothing, and it was dropped. The remaining 15 or so percent gap to cuBLAS is pure kernel arithmetic, ours at about 13 TFLOP/s against their 15 on a 19.5 ceiling.

## Phase 6. MPI (this session)

Two files changed. `matrixMultiply.cpp` was refactored so the planar 3M machinery lives in a function that computes only a range of output rows, and `matrixMultiplyMPI.cpp` calls it once per rank and stitches the results together.

### The refactor

The whole body of `matrixMultiply` moved into a new function `matrixMultiplyRowRange` taking a row range [r0, r1). The public function now just calls it with the full range, so the CPU benchmark path is byte-for-byte the same arithmetic as before. The range version restricts each stage to its band. The left input is only deinterleaved for the band's rows, the result planes are only zeroed there, the compute loop only walks those row tiles, and the merge only writes those rows of C. The right input is still processed in full, because every output row needs every row of the right operand. The function is deliberately not static, meaning its name is visible to other compiled files, because the MPI file declares it extern and calls it. Both object files link into the one GradeBot binary, so this works without touching any header.

One subtlety sits at the band edges. The AVX kernel processes rows four at a time, so a band starting at row 37 begins its compute at row 36. The extra rows compute throwaway values into the scratch planes, which is harmless because each output row depends only on its own row of the left input, and the merge step never reads the throwaway rows. The local test splits matrices at deliberately unaligned boundaries to prove this.

### The MPI layer

Rank p of P gets rows so that every rank has N/P and the first N mod P ranks get one extra, the standard remainder-spreading split. Each rank computes its band directly into its own C buffer at the right offset, then one collective call, `MPI_Allgatherv`, exchanges the bands so every rank ends up holding the complete matrix. The v in the name means the piece sizes may vary, which the remainder split needs. `MPI_IN_PLACE` tells MPI the rank's own contribution is already sitting inside the destination buffer, so nothing is copied locally.

The data travels as plain floats with all counts doubled, two floats per complex number. MPI has complex datatypes, but they come from Fortran and C99 corners of the standard and occasionally misbehave in C++ builds, and raw floats cannot break. The N less or equal 0 guard runs before any MPI call, because the harness harvests the student ID by calling with N=0, possibly before MPI exists. If MPI is somehow not initialised or the rank count is unusable, the code falls back to computing the full matrix locally, which is slow but always correct.

The plan sketched a pipelined version that overlaps communication with compute in chunks. It was not built, because the measurement said the plain blocking gather already fits. Measure before adding complexity.

### Validation

The refactor was tested locally first, x86 AVX under Rosetta as in Phase 2, with three-way splits at unaligned row boundaries in both the 3M and 4M schemes, all matching a double precision reference. Job 568787 (DebugMPI, 2 ranks on one node, N=128) passed with error 3.4e-08. Job 568790 (2 nodes, N=2048) landed on a100-a and a100-b, the exact pair JudgementDay uses, and scored ratio 0.496 and 0.497 on the two ranks, grade 7.00, error 5.0e-09 identical to the single-node 3M error. The two-node time of 0.149 s against the single-node 0.240 s says the 16.8 MB gather costs only a few milliseconds on this interconnect.

### Where this leaves the project

All three ladders now measure at 7 marks on graded silicon, CPU 0.805 against the 1.0 line, GPU 1.316 against 1.5, MPI 0.496 against 0.6. Phases 0 through 6 are complete. What remains is Phase 7, one budgeted JudgementDay dress rehearsal to confirm all three numbers on the pinned pair in a single run, and Phase 8, the freeze and packaging. JudgementDay is the shared, lecturer-flagged resource, so that run waits for an explicit go.

## Phase 7. Dress rehearsal (this session)

No code changed. With the go given, the pinned JudgementDay script ran once as job 568818, its first permitted use. It executes all three benchmarks on a100-a and a100-b in a single job, exactly as the graded run will.

Every line came back at 7.00. CPU scored 0.807 and 0.809 on the two nodes, matching the vcpu-4 calibration almost exactly, which confirms the vcpu nodes really do pin the graded silicon. MPI scored 0.483 and 0.485, slightly better than the bench pair. GPU scored 1.372 and 1.383 against the bench's 1.316. That drift is on the reference side, our own throughput held between 167 and 169 matrices per second across four separate jobs while the cuBLAS figure moved between 221 and 234. The thinnest margin in the set is the GPU's 9 percent, which is enough that ordinary run-to-run variance cannot cross the 1.5 line.

The plan allowed for triage and one fix here. Nothing needed fixing, so the phase closed after a single run, leaving the second budgeted run for the freeze confirmation in Phase 8.

### Where this leaves the project

The dress rehearsal shows 21 of 21 on the graded hardware and the graded script. The only remaining work is Phase 8, stripping nothing (no debug artifacts were ever added), re-verifying the immutables, one clean rebuild and JudgementDay confirmation on rangpur, and packaging the zip for submission.

## Phase 8. Freeze and package (this session)

No code changed. The phase is verification and packaging.

The immutable checklist was walked first. All four function signatures match the headers byte for byte, the STUDENTID defines and returns are in place, the N less or equal 0 guards are the first statement of each function and the MPI one runs before any MPI call, the only includes are the three provided headers, and no file prints or writes anything. There were no debug artifacts to strip because none were ever added, every correctness check in the project ran locally in emulation instead.

A fresh rsync then went up and JudgementDay ran once more as job 568931, the second and final budgeted use, on a clean rebuild. All six lines scored 7.00 again. CPU came in at 0.810 and 0.811, MPI at 0.501 and 0.504, and GPU at 1.403 on both nodes. The GPU reference posted its fastest number of the whole project on this run, 236 matrices per second against 222 on the slowest, while our own throughput stayed inside one percent across every job. The margin under the 1.5 line at the worst observed reference speed is 6.5 percent.

The package was then built on rangpur exactly as the plan specifies. The 21 slurm output files went into slurm.zip, and 48008361.zip was made from the three source files plus slurm.zip. The listing shows exactly those four entries. The zip was fetched and unpacked locally, and the three sources inside were diffed against the working tree, byte identical.

### Where this leaves the project

The assignment is finished and frozen at a confirmed 21 of 21. The submission file is 48008361.zip, on rangpur in the build directory and locally in results/. The one remaining act is uploading it to Blackboard, and the spec deadline should be confirmed there before assuming the harness's countdown is the true date.

## Phase 9. AI declaration and line comments

No code changed. An AI declaration was added at the top of the three submission files. It states that Claude helped debug the code, draft pseudocode and implement it. A brief comment was also added above nearly every line so each file can be explained line by line in an interview. Comments sit on their own line above the code they describe, and the side comments already in the code were moved up the same way. The only side comments left are the two that came with the provided stubs.

The stub lines that carry trailing whitespace were left byte identical. Stripping the comments from the old and new files gives identical code, so the Phase 8 results still describe what gets compiled. The zip was then rebuilt locally. slurm.zip was rebuilt from the assignment directory on rangpur so it holds all 36 slurm outputs, adding 600266 and 600267, which the Phase 8 package had left out. The new results/48008361.zip holds the three commented files and that slurm.zip. The Phase 8 zip was kept as results/archive/48008361_before_declaration.zip. The new files were compiled on rangpur with the grader's flags beside the old ones. The CPU and MPI objects came out byte identical, and the GPU object has identical machine code, PTX and host code. results/zipcheck_final/ZIPCHECK.md records the checks, with both zips unpacked and a diff per source file.

