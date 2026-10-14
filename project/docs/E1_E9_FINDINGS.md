# E1 and E9 findings, and the slow test that delayed them

This file collects three results from 2026-10-14 in one place so they can go
into the deck, the README or the interview. The full record is in
`docs/RESULTS.md` M2.6 and `CHANGELOG.md` Phase 11. Every number is a median
of five reps in Mcell-updates/s on rangpur's r730 nodes, and every job passed
`make test`, `make test-mpi` and `make avx-check` before it measured
anything.

## E1. Could the compiler have written the SIMD rung?

The question is whether all the hand-written AVX2 work was needed, or whether
a compiler flag would have done it. The whole program was rebuilt with each
set of flags in place of `-O2`, and the baseline, serial and simd rungs were
timed on one thread at 2048 x 512 (jobs 659153 and 659154). The price and
Greeks were printed to 17 digits for every build to see how far the answer
moved.

| Build | baseline | serial | simd |
|---|---|---|---|
| `-O2` (shipped) | 57.6 | 158.4 | 370.2 |
| `-O3` | 112.5 | 241.5 | 377.0 |
| `-O3 -mavx2 -mfma` | 107.4 | 348.7 | 381.1 |
| `-O3 -march=native` | 106.0 | 265.5 | 455.4 |
| `-O3 -ffast-math` | 189.6 | 271.2 | 386.8 |
| `-O3 -mavx2 -mfma -ffast-math` | 381.0 | 482.2 | 378.4 |

The compiler can produce code close to the SIMD rung, but only when it is
given flags the shipped build does not use. Given AVX2 and FMA, g++
vectorises the plain scalar row four lanes wide and gets within 6% of the
hand-written row. It cannot prove the input and output rows never overlap, so
it also adds a run-time check and a scalar fallback.

The largest change comes from `-ffast-math`. It lets the compiler treat
floating point as exact arithmetic, so it may reorder sums and turn a
division into a multiplication by the reciprocal. The baseline divides by the
grid spacings inside its innermost loop, and those divisions made it slow in
Milestone 1. With `-ffast-math` the compiler computes each reciprocal once
before the loop. The disassembly shows the vectorised inner loop going from
10 divide instructions to none. The unmodified baseline then runs at 381, as
fast as the hand-written SIMD rung, and the scalar serial rung reaches 482,
1.27 times faster than it.

The shipped build still does not use `-ffast-math`. The answers moved by at
most 1.1e-14, but the flag relaxes the IEEE rules for every loop in the
program. It also lets the compiler assume no value is ever NaN or infinite,
and the tests detect an unstable run by exactly such a value. This result
also supports the Milestone 1 rule. At `-O2` the hand optimisations worth
doing are the ones the compiler is forbidden to do, and E1 shows the
compiler making the same change as soon as it is permitted.

E1 has two smaller results. The baseline gains nothing from AVX2 alone
(112 to 107), because it is limited by division and Haswell executes a
256 bit divide as two 128 bit divides. `-march=native` raises simd to 455,
because tuning for Haswell stops g++ 8.5 splitting each unaligned 256 bit
load into two 128 bit loads (the E2b finding), but it slows the
auto-vectorised serial rung to 265. The cause of that slowdown has not been
found.

## E9. Can the edge row exchange be overlapped with the update?

Across two nodes the edge row exchange is 38% of a step on the small grid
(M2.9). The `mpi` rung updates every row, then sends and receives the edge
rows and blocks until every send and receive has completed. The new
`mpi_overlap` rung updates the two edge rows first, and the master thread
posts their non-blocking sends and receives while the other threads update
the remaining rows. It calls `MPI_Waitall` only at the end of the step. Both
ran back to back in every hybrid layout (job 659155). The exchange share is
in brackets.

| Layout | 2048 x 512, mpi | overlap | 4096 x 1024, mpi | overlap |
|---|---|---|---|---|
| 1 node, 1 rank x 24 | 5731 | 5482 | 2821 | 2900 |
| 1 node, 2 ranks x 12 | 6069 (19%) | 6121 (20%) | 6298 (10%) | 6114 (10%) |
| 2 nodes, 1 rank x 24 each | 5567 (38%) | **8135 (15%)** | 8519 (28%) | 8843 (28%) |
| 2 nodes, 2 ranks x 12 each | 8672 (37%) | 8309 (35%) | 9833 (28%) | 9868 (30%) |

Overlapping helps in one layout and makes no measurable difference in the
others. With one rank on each of two nodes, the small grid runs 1.46 times
faster and the exchange falls from 38% of a step to 15%. The ranges of the
reps do not overlap (5493 to 5801 against 7884 to 8383), so the gain is
real. Everywhere else the rate moves by less than 5% and the exchange share
is unchanged.

The likely reason is how Open MPI 4.1.1 transfers a message, with the
settings read from `ompi_info` on rangpur. Between nodes it uses TCP with an
eager limit of 64 KiB, so a message up to that size is sent as soon as
`MPI_Isend` is called. Between two ranks on the same node it uses shared
memory, which has an eager limit of 4 KiB. A row is 16 KiB on the small
grid, so it is sent with the rendezvous protocol, where the data is
transferred only once the receiving process is inside an MPI call. This Open
MPI build has no progress thread, so that call is `MPI_Waitall` at the end of
the step, and the transfer does not overlap with the update. The one layout
that gains is the only one where no two ranks share a node. The explanation
fits every small grid row of the table. It does not explain the big grid's
1.04x in the same layout, and it has not been tested. Raising
`OMPI_MCA_btl_vader_eager_limit` above the row size would test it in a five
minute job.

A non-blocking call does not by itself make the transfer proceed while the
process computes. Communication overlaps with computation only when the
library can transfer the data without being called again.

## The slow test that timed out the first attempt

The first E1 and E9 jobs (657626 to 657628) all hit the 14 minute limit
before a single benchmark finished, because validation alone had grown from
8 to 14 minutes. The bench scripts pin OpenMP threads to cores
(`OMP_PROC_BIND=close`), and `make test-mpi` inherited that. Every MPI
process put its two threads on cores 0 and 1. A waiting OpenMP thread
busy-waits, so one rank's busy-waiting threads occupied the cores its
neighbour needed to finish the update and send its row. One 3-process solve
of the tiny test grid took 14.2 s pinned and 0.6 s unpinned. Letting the
threads sleep while they wait (`OMP_WAIT_POLICY=passive`) brought it to
1.2 s, which shows the busy-waiting OpenMP threads were most of the cost.
Setting `OMP_PROC_BIND=false` on the test line cut validation to about
2.5 minutes. It is a small, real example of oversubscription. It also shows
that thread binding has to be chosen for each kind of job, since the binding
that suits a benchmark slowed the tests.

## Ways to show them

- **E1 as one slide.** The slide shows grouped bars of baseline, serial and
  simd for `-O2`, `-O3 -mavx2 -mfma` and `-O3 -mavx2 -mfma -ffast-math`. The baseline bar
  rising to the height of the simd bar shows the result directly. The line
  to say is that the compiler matched the hand-written rung only by relaxing
  the IEEE rules for the whole program.
- **E9 as one slide.** The slide shows paired bars of `mpi` and
  `mpi_overlap` for the four layouts on 2048 x 512, with the exchange share written on each bar. The
  one layout with a large gain is the result to point at, and the eager
  limit explains it.
- **README.** The README gets one line for each under the experiments,
  pointing at RESULTS M2.6.
- **Interview.** E1 answers "why write intrinsics when the compiler can
  vectorise", and E9 answers "did you try to overlap the communication with
  the computation".

The data for both slides is in `results/bench_m2_exp_compiler_vector_659153.csv`,
`results/bench_m2_exp_compiler_fastmath_659154.csv` and
`results/bench_m2_exp_overlap_659155.csv` (local only, `results/` is not in
git). The 17 digit answers and the vectorisation reports are in
`results/e1_answers_*` and `results/e1_vectorised_*`.

These numbers were measured with the shipped AVX2 build of that time, which
still split each unaligned 256 bit load into two 128 bit loads. A rerun of the sweeps without the split would change the
simd and `mpi` figures, and E9 would need rerunning with it.
