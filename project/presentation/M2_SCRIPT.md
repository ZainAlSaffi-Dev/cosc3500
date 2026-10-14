# M2 SCRIPT, the word-for-word script for the Milestone 2 video

One section per take, in the order of the recorded frames of
`presentation/deck_m2.pdf` (frames 1 to 35). The appendix after frame 35
is for the interview and is never recorded. Each heading gives the start
time and the length of the take, worked out from its word count at 140
words per minute plus a second for each slide change. The Milestone 1
marker found the delivery fast, which is why the pace is set at 140 words
per minute. The assessors stop watching at 10 minutes. If a take runs long
when spoken, cut words from it and keep the pace.

The takes are written to be spoken. They use "we" for the work, as a paper
would, contractions where they sound natural, and the technical name for
each operation, so the code on screen and the words match.

Every number here comes from `docs/RESULTS.md`, Milestone 2, where each one
has its job id, and matches the number block at the top of
`presentation/deck_m2.tex`. Numbers are written as digits so they are easy
to read off.

## Before recording

These come from the Milestone 1 feedback (`PLAN.md` §11a).

- No streaming overlay. The webcam sits in a bottom-right box that no
  slide uses, so it never covers a header or an axis.
- Turn on a large cursor and point at the number being said on every
  figure.
- Say the named results slowly. Amdahl is AM-dahl. Gustafson is
  GUS-taf-son. Heston is HESS-ton. NUMA is NOO-mah. AVX2 is A V X two.
- Record take by take. A take that runs over its slot is re-recorded at the
  same pace.

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
- Placement, where Linux allocates a page in memory. U. Drepper, What every
  programmer should know about memory, Red Hat, 2007, section 5.
- Weak scaling. J. L. Gustafson, Reevaluating Amdahl's law, Communications
  of the ACM 31(5), 532 to 533, 1988.

## The takes
## 0:00 Title (0:13)

On screen, frame 1, the title.

Hi, I'm Zain. In this video I'll show how we made a Heston option pricer
123 times faster on one rangpur node, with the same answer from every
version.

## 0:13 The project in one slide (0:17)

On screen, frame 2, the four blocks.

We're pricing a call option when the share price and its volatility both
move randomly, by solving the Heston PDE on a grid. No cell in a step
depends on another, so each step can be split up.

## 0:30 The problem (0:16)

On screen, frame 3. Point at the payoff, then at the circle on the value map.

At expiry the option is worth its payoff, on the left. Heston lets the
variance move randomly too, and half way through the solve the value above
the payoff rises with variance near the strike.

## 0:46 The equation (0:11)

On screen, frame 4. Run the cursor down the table.

The option value satisfies this PDE, and the table explains each term. On
the grid each derivative becomes a finite difference of neighbouring cells.

## 0:58 One step on the grid (0:18)

On screen, frame 5. Point at the payoff, the red block, then the red cell.

Starting from the payoff, each new value comes from the nine-point stencil
around it, in red. A step reads only the old
buffer and writes only the new one, then the two buffers swap pointers, so
no values are copied.

## 1:16 Checking the answer (0:13)

On screen, frame 6.

To validate the solver, we set the volatility of volatility to zero, which
reduces Heston to Black-Scholes, and the unchanged solver is within 4 cents
of the closed-form price.

## 1:29 The data structure (0:18)

On screen, frame 7.

Both buffers are flat arrays of doubles, with cell i, j at j times the row
length plus i. So the share price axis is contiguous in memory, and a step
in variance is a stride of one row.

## 1:47 How every number was measured (0:19)

On screen, frame 8.

Every number is the median of 5 runs on one exclusive rangpur node through
sbatch, after the tests pass in the same job. We discarded our
first sweeps because denormal numbers slowed some versions, and reran
everything with them flushed to zero.

## 2:06 Where the time goes (0:11)

On screen, frame 9.

Before changing anything we profiled the code with gprof, and one
function, the update of the interior rows, takes 99 percent of the time.

## 2:17 What the compiler does alone (0:15)

On screen, frame 10.

We first tested the compiler alone. O2 takes the baseline from 11 million
updates a second at O0 to 58, and every benchmark builds at O2, so each
version measures only our change.

## 2:32 The ladder (0:14)

On screen, frame 11. Point along the bars from left to right.

This is the ladder. Each version adds one technique and must match the
baseline before it's timed, and the number in each bar is the gain over
the version below.

## 2:46 The serial rung (0:22)

On screen, frame 12. Point at the division, the reciprocal, then the
multiplication.

The baseline divides by the grid spacing in every cell. The compiler can't
turn that into a multiplication by the reciprocal, because the rounding can
differ in the last bit, so we did it by hand, once per solve, for 159 million updates a
second, 2.8 times the baseline.

## 3:08 SIMD, four cells per instruction (0:16)

On screen, frame 13.

Serial was already at 96 percent of the scalar peak, so next we used vector
instructions. An AVX2 register holds four doubles, and neighbouring cells
use the same formula, so one instruction updates four cells.

## 3:24 SIMD code, before and after (0:18)

On screen, frame 14. Point from each serial line to the AVX2 line below it.

The serial loop is on top and the AVX2 loop below. Each scalar load
becomes a load of four doubles, the update becomes one fused multiply-add, and the counter
goes up by four. This runs at 456, 2.9 times serial.

## 3:42 The compiler was splitting the loads (0:22)

On screen, frame 15. Point at the two bars, then at the vinsertf128 lines.

Our first build ran this at 387, and the disassembly showed that g++ 8.5,
with only the AVX2 flag, splits every unaligned 256 bit load into two 128
bit loads. Two extra flags make each a single load, for 1.19 times the
speed, and every number here uses them.

## 4:04 Could the compiler have written it? (0:19)

On screen, frame 16.

We also tested whether the compiler could write this itself. With AVX2 it
vectorises the serial row to 77 percent of our version. Fast-math also
removes the divisions, but relaxes the IEEE rules for the whole program, so
we didn't use it.

## 4:23 SIMD on the roofline (0:13)

On screen, frame 17.

On the one-core roofline, serial is at the scalar peak. SIMD reaches 14
gigaflops, a third of the AVX2 ceiling and just below the memory roof of
one core.

## 4:37 OpenMP, one block of rows per thread (0:17)

On screen, frame 18.

With OpenMP each thread gets a contiguous block of rows. The threads read
the old buffer and write only their own rows, so there are no data races
and no locks. Then the main thread updates the boundaries.

## 4:54 OpenMP code (0:11)

On screen, frame 19. Point at the pragma.

The code change is one pragma above the row loop. On 24 cores this gives
6367 million updates a second, 14 times SIMD.

## 5:05 Strong scaling (0:19)

On screen, frame 20.

With the grid fixed, 24 threads give a speedup of 13.9. Fitting Amdahl's
law gives a serial fraction of 2.4 percent, while the boundary work is 0.4
percent of a step, so most of the loss is the barrier and the shared
caches.

## 5:24 Weak scaling (0:17)

On screen, frame 21.

For weak scaling each thread keeps 21 rows, so the grid grows with the
threads. The step time goes from 0.11 to 0.17 milliseconds, a scaled
speedup of 16 on 24 threads, below Gustafson's prediction of 22.

## 5:41 Grid size (0:17)

On screen, frame 22.

Small grids are dominated by the fixed cost of each step. While both
buffers fit in the two L3 caches, 24 threads run 40 times serial, and at
4096 by 1024, where they don't fit, the rate halves.

## 5:58 After threads, memory bandwidth (0:16)

On screen, frame 23.

On the node roofline the small grid runs above the memory roof, which is
the result of it fitting in L3. The big grid is on the memory roof, so it's
limited by memory bandwidth.

## 6:14 Placement (0:24)

On screen, frame 24. Point at the 4096 pair on the left, then the two pairs
on the right.

On the big grid, 12 threads gave 1894 on one socket and 3413 on both.
Linux places each page in the memory of the socket whose thread first
writes it, and our main thread fills the grid, so every page is on socket 0.
Interleaving the pages made 24 threads 1.74 times faster.

## 6:38 MPI, one block of rows per process (0:19)

On screen, frame 25.

With MPI each process owns a block of rows. Every step it sends its first
and last rows to its neighbours and receives theirs into its copied rows.
With one process per socket, each fills its own grid, so its pages are
local.

## 6:57 MPI code, the time loop (0:17)

On screen, frame 26. Point at each marked line in turn.

In the loop, rows for rank gives each process its block. It updates only
those rows and its own boundaries, then exchanges edge rows before the
swap. At the end MPI Gatherv collects the blocks onto process 0.

## 7:15 MPI code, the edge row exchange (0:12)

On screen, frame 27.

The exchange is non-blocking. MPI Irecv is posted before MPI Isend, both
return immediately, and MPI Waitall blocks until every send and receive has
completed.

## 7:26 Hybrid layouts (0:21)

On screen, frame 28.

One process per socket made the big grid 2.28 times faster than one
process of 24 threads. Two nodes give
1.21 times one node on the small grid and 1.64 on the big one, and on the
small grid the exchange takes 40 percent of each step.

## 7:48 Overlapping the exchange with the update (0:25)

On screen, frame 29. Point at the tall blue bar.

So we tried overlapping the exchange with computation, sending the edge
rows first and updating the rest while they're in transit. That only helped
with one process per node, at 1.41 times. We suspect that within a node
Open MPI waits for the receiver to call MPI Waitall before it transfers a
16 kibibyte row.

## 8:12 Where each step goes (0:14)

On screen, frame 30.

At 24 threads the serial boundaries grow to 9 percent of each step, because
the main thread reads rows other cores just wrote, and MPI adds the exchange
on top.

## 8:26 Same answer on every rung (0:13)

On screen, frame 31.

No version was timed until it matched the baseline. Every version is within
3.3 times 10 to the minus 15, and OpenMP and MPI match SIMD bit for bit.

## 8:39 Tried and rejected (0:17)

On screen, frame 32.

These are the variants we tested against OpenMP. Dynamic scheduling was
slower, because every row costs the same. Float was faster, but it moved
the price by 4 cents and gamma by 12 percent, so we rejected it.

## 8:57 Conclusion (0:18)

On screen, frame 33.

To conclude, we reached 123 times the baseline on one node and 147 on two.
Most of that is the result of SIMD and threads, and across nodes the
exchange is the main limit. Four of our six hypotheses held.

## 9:15 What surprised me, and what I would do next (0:17)

On screen, frame 34.

Looking back, I was surprised how much build and system settings changed
the results, from split loads to huge pages. Next I'd store only each
process's own rows, test our E9 explanation, and try a GPU version.

## 9:32 Close (0:11)

On screen, frame 35.

Every number here has a job id in the results file, and the repository has
the code to reproduce it. Thanks for watching.
