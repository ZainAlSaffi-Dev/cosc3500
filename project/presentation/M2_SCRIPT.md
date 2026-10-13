# M2 SCRIPT, the word-for-word script for the Milestone 2 video

One section per take, in the order of the beats in `PLAN.md` §11c. Each
heading gives the start time and the slot. The script is budgeted at 140
words per minute, slower than the 150 used for Milestone 1, because the
Milestone 1 marker found the delivery fast. The whole script is about 1270
words, which is 9 minutes at that pace. The slots add up to 9 minutes 40
seconds, and the assessors stop at 10. The equation, grid and
Black-Scholes takes are there because the marker may not have seen the
Milestone 1 video. If a take runs long when spoken, cut words from it rather than
speaking faster.

Every number here comes from `docs/RESULTS.md`, Milestone 2, where each one
has its job id. Numbers are written as digits so they are easy to read off.

## Before recording

These come from the Milestone 1 feedback (`PLAN.md` §11a).

- No streaming overlay. The webcam sits in a bottom-right box that no
  slide uses, so it never covers a header or an axis.
- Turn on a large cursor (macOS Accessibility, Display, Pointer size) and
  point at the number being said on every figure.
- Say the named results slowly. Amdahl is AM-dahl. Gustafson is
  GUS-taf-son. Heston is HESS-ton. NUMA is NOO-mah. Fourier is FOO-ree-ay.
  AVX2 is A V X two.
- Record take by take. A take that runs over its slot is re-recorded, not
  sped up.

## Sources on the slides

Amdahl's law is in lecture 3. The task sheet names Gustafson's law, but the
lectures do not teach it. The other ideas below are not in the lecture
material either, so each slide that uses one carries a small grey source line.

- The protocol slide, denormal numbers. D. Goldberg, What every computer
  scientist should know about floating-point arithmetic, ACM Computing
  Surveys 23(1), 5 to 48, 1991.
- Both roofline slides. S. Williams, A. Waterman and D. Patterson, Roofline,
  an insightful visual performance model for multicore architectures,
  Communications of the ACM 52(4), 65 to 76, 2009.
- Strong scaling, fitting Amdahl's serial fraction to measured speedups. A. H.
  Karp and H. P. Flatt, Measuring parallel processor performance,
  Communications of the ACM 33(5), 539 to 543, 1990.
- Placement and hybrid MPI, pages living on the socket that first touched
  them. U. Drepper, What every programmer should know about memory, Red Hat,
  2007, section 5.
- Weak scaling. J. L. Gustafson, Reevaluating Amdahl's law, Communications of
  the ACM 31(5), 532 to 533, 1988.

## 0:00 Cold open (0:15)

On screen, the title slide.

This video shows how I made a Heston option pricer 101 times faster on one
node of the rangpur cluster, with the same answer on every version.

## 0:15 The problem (0:30)

On screen, `results/payoff.png`, then the weather map.

A European call is the right to buy a share at a fixed strike price on a
fixed day. On that day its value is the payoff. The question is its price
today. Black-Scholes holds the volatility constant. Heston lets the variance
move randomly too, so the price depends on the share price, the variance and
time.

## 0:45 The equation (0:30)

On screen, the equation slide.

The price V solves this equation. The second derivatives come from the
randomness in the share and its variance, the mixed one from their
correlation, and the first derivatives from the drift. Every other symbol is
a known number. On a grid each derivative becomes a difference of neighbouring
cells, and a time step adds dt times these terms to the old value.

## 1:15 The grid (0:30)

On screen, F13. Point at the payoff labels, the red block, the key, then
the red cell.

At expiry the value is known. Every cell holds the payoff, zero below the
strike and S minus K above it. So the solve starts there and steps back. The
new value at i, j comes from the red block of nine around it. The key shows
which derivative each cell feeds. No new cell depends on another new cell, so
they can all be done at once. Then the sheets swap, back to today.

## 1:45 The check against Black-Scholes (0:20)

On screen, F14.

To check the answer, I set the volatility of volatility to zero. The variance
then never moves and Heston becomes Black-Scholes, which has a closed form.
The solver is not told this. It runs the full stencil and lands 4 cents from
it.

## 2:05 The data structure (0:25)

On screen, `results/m2/code_shots/shot_data_layout.png`.

In memory the two sheets are flat arrays of doubles. A cell sits at row times
row length plus column, so the share price axis is contiguous and one step in
variance jumps a whole row. At 2048 by 512 a row is 16 kilobytes and the two
sheets are 16 megabytes.

## 2:30 How the numbers were measured (0:25)

On screen, the protocol slide.

Every number comes from one exclusive rangpur node through sbatch, the
median of 5 runs with threads pinned to cores. The tests run first in the
same job. I threw away a first round of sweeps, because tiny denormal numbers
ran slowly, and differently on each grid. Every sweep was rerun with them
flushed to zero.

## 2:55 Where the time goes (0:30)

On screen, F12, the gprof flat profile, then F11, the four -O levels.

Before changing anything I measured where the time goes. The gprof profile
puts 99.9 percent of the baseline in its one solve function, and over 98
percent of the serial and SIMD rungs in the row update, so every rung
attacks that. The compiler alone takes the baseline from 11
million updates a second at O0 to 58 at O2, and O3 doubles that by
using two-lane SSE.

## 3:25 The serial rung (0:30)

On screen, `results/m2/code_shots/shot_serial_row.png`.

The serial rung keeps the one change from the first milestone that paid.
The baseline divides by the grid spacing five times per cell. The compiler
may not turn those divisions into multiplications, because that changes the
last bits of the answer, so I did it by hand, once per solve. That gives 159
million cell updates a second, 2.8 times the baseline.

## 3:55 SIMD (0:55)

On screen, `results/m2/code_shots/shot_simd_fma.png`, then F5.

The serial rung already runs at 96 percent of the scalar peak, so
the only way up is wider instructions. An AVX2 register holds four doubles,
and every cell in a row runs the same formula on neighbours that sit next to
each other in memory, so four cells can share each instruction. Each fused
multiply-add does a multiply and an add with one rounding. By hand, with four lanes and
FMA, the rung measures 379 million updates a second, 2.4 times serial. It is not 4, because each
group of four loads nine overlapping windows, and the cells left at the end
of a row run one at a time. On the roofline this rung does 12 gigaflops,
against an AVX2 ceiling of 43 for one core.

## 4:50 OpenMP (1:05)

On screen, `results/m2/code_shots/shot_omp_pragma.png`, then F2, then F7.

OpenMP shares the rows of each step between threads with one line. Each
thread reads only the old sheet and writes only its own rows, so no locks are
needed. On 24 cores
it gives 15.8 times. An Amdahl fit gives a serial fraction of 1.2 percent,
three times the boundary work, so most of the loss is the barrier on every
step and the shared cache bandwidth. The scalar version
scales better, 18.5 times, because it moves data at half the rate. Placement
mattered as much as count. On the big grid, twelve threads on one socket
ran at 1880, the same twelve spread over both sockets at 3319, with two L3 caches
and two memory controllers. And every page of the grid lives on the socket
that first wrote it. Interleaving the pages gave the big grid 28 percent
more.

## 5:55 Hybrid MPI (0:50)

On screen, `results/m2/code_shots/shot_mpi_exchange.png`, then F8.

MPI gives each process a slab of rows. Each step, every process swaps one
edge row with each neighbour. On one
node I run one process of 12 threads per socket. Each process first touches its own copy of the grid on its own socket,
which fixes the page problem. That is 1.05 times one
24-thread process on the 2048 grid, and 2.16 times on the 4096 grid. Across
two nodes the best layout is 1.41 times one node, 70 percent efficient. The
halo is only 64 kilobytes, but each exchange costs about 45 microseconds of
latency and waiting, which is 38 percent of a step this short.

## 6:45 Same answer on every rung (0:20)

On screen, F9.

No version was timed until it matched the baseline. Every rung is within 3.3
times 10 to the minus 15 on price, delta, gamma and vega, and the OpenMP and
MPI rungs equal the SIMD rung exactly.

## 7:05 Weak scaling and grid size (0:40)

On screen, F3, then F4.

Weak scaling keeps the work per thread fixed and grows the grid with the
threads. The step time rises from 0.14 to 0.18 milliseconds, a scaled
speedup of 19 on 24 threads. Gustafson's law with the measured 8 percent
serial share predicts 22. Grid size decides the regime. Small grids lose
to the fixed cost of each step. The 2048 grid fits in the two L3 caches and
runs 36 times faster than serial. At 4096 by 1024 the sheets spill to
memory and the rate halves.

## 7:45 Bottlenecks (0:35)

On screen, F6, then F5.

This is where each step's time goes. On one thread nearly all of it is the
interior rows. At 24 threads the serial boundaries grow to 8 percent of the
step, because the main thread has to fetch rows other cores just wrote. That
overturned my guess of half a percent. Across nodes, the halo exchange takes
over. The roofline puts the large grid on the memory roof, so after threads
the limit is bandwidth, not arithmetic.

## 8:20 Tried and rejected (0:30)

On screen, F1, the right panel.

On the right, each bar is a variant against the OpenMP rung. Dynamic
scheduling lost up to 40 percent, because every row costs the same. Float
instead of double was 1.75 times faster, and 4.5 times on the big grid
because it fits back in cache. I rejected it, because on the accurate grid
it moved the price by 4 cents and gamma by 12 percent.

## 8:50 Conclusion (0:40)

On screen, the conclusion slide, the ladder numbers beside the measured
ceilings.

The ladder ends at 101 times the baseline on one node and 155 times on two.
The measured ceilings explain it. Serial sat at the scalar peak, so
SIMD was the only way up, and it reached a quarter of the AVX2 ceiling.
Threads took the large grid onto the memory roof, which measures 68
gigabytes a second.
Four of my six hypotheses held. The two that failed were the serial boundary
and the network, both costs that grow with parallelism. Next would be
overlapping the halo exchange with the interior rows.

## 9:30 Close (0:10)

On screen, F1.

Every number has a job id in the results file, and the repository has the
code and commands to reproduce it. Thank you.
