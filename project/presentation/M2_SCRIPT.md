# M2 SCRIPT, the word-for-word script for the Milestone 2 video

One section per take, in the order of the recorded slides of
`presentation/deck_m2.html` (slides 1 to 38). A take runs over more
than one slide where a code slide changes into the next rung's code, and the
notes under each slide (press S) hold the same take. The appendix after slide
38 is for the interview and is never recorded. Each heading gives the
start time and the length of the take, worked out from its word count at 150
words per minute plus a second for each slide in it. That is a little faster
than the 140 words per minute used before, because the animations show the
working while the take is spoken, although the Milestone 1 marker found the
delivery fast, so the numbers still need to be said slowly. The 1375
words come to 9:48, and the assessors stop watching at 10 minutes. If a
take runs long when spoken, cut words from it and keep the pace.

The takes are written to be spoken. They use "we" for the work, as a paper
would, contractions where they sound natural, and the technical name for
each operation, so the code on screen and the words match.

Every number here comes from `docs/RESULTS.md`, Milestone 2, where each one
has its job id, and matches the numbers under `n` at the top of
`presentation/deck_m2.qmd`. The test grid numbers in the opening animations,
the first step value 0.04959375 and the price 195.6355, come from a local run
of the 421 by 51 test grid with `--dump-every`, so they are not benchmark
numbers. Numbers are written as digits so they are easy to read off.

## Before recording

These come from the Milestone 1 feedback (`PLAN.md` §11a).

- No streaming overlay. The webcam sits in a bottom-right box that no
  slide uses, so it never covers a header or an axis.
- Turn on a large cursor and point at the number being said on every
  figure.
- Say the named results slowly. Amdahl is AM-dahl. Gustafson is
  GUS-taf-son. Heston is HESS-ton. NUMA is NOO-mah. AVX2 is A V X two.
- The four animations in the opening stop at the end of each step and move
  on with the next key, so press next as each new sentence starts.
- Record take by take. A take that runs over its slot is re-recorded at the
  same pace.

## Sources on the slides

Amdahl's law is in lecture 3. The task sheet names Gustafson's law, but the
lectures do not teach it. The other ideas below are not in the lecture
material either, so each slide that uses one carries a small grey source line.

- The model slide, the Heston model and its PDE. S. L. Heston, A closed-form
  solution for options with stochastic volatility with applications to bond
  and currency options, Review of Financial Studies 6(2), 327 to 343, 1993.
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

On screen, slide 1.

Hi, I'm Zain. In this video I'll show how we made a Heston option pricer
123 times faster on one rangpur node, with the same answer from every
version.

## 0:13 The project in one slide (0:16)

On screen, slide 2. The four blocks.

We're pricing a call option when the share price and its volatility both
move randomly, by solving the Heston PDE on a grid. No cell in a step
depends on another, so each step can be split up.

## 0:29 From the model to the PDE (0:25)

On screen, slide 3. Four presses, one per step of the derivation. Each press plays one step and stops on its last frame.

Heston lets the share price and its variance both move randomly, with the
variance pulled back towards 0.04 and correlated with the price at minus
0.7. The option is worth its expected discounted payoff, and Itô's lemma,
together with the rule that the discounted value can't drift, gives the
Heston PDE, each coloured term coming from one part of the model.

## 0:54 From the PDE to one line of C++ (0:24)

On screen, slide 4. Four presses, the explicit step, the differences, the reciprocals and the C++.

We step back from expiry one explicit Euler step at a time, and since that
is only stable below a set timestep, the full grid needs 980,000 steps. Each
derivative becomes a central difference over the nine-point stencil, the
divisions become reciprocals worked out once, and the whole update ends up as
these two lines of C++.

## 1:18 One cell with real numbers (0:21)

On screen, slide 5. Four presses, the cell, the differences, the terms, the new value.

Here's one real cell in the first step, at the strike, where only the slope
and curvature along the share price are nonzero. The update gives 0.04959,
exactly what the solver's own snapshot holds, and the cells above and below
get 0.074 and 0.025, because the diffusion term grows with variance.

## 1:39 One step over the whole grid (0:16)

On screen, slide 6. Point at the payoff, the red block, then the red cell. Six clicks: three cells, the rest of the row, the other rows, the swap.

Every interior cell gets this same update, reading only the old buffer and
writing only the new one, so the cells can be computed in any order. Then
the two buffers swap pointers, so no values are copied.

## 1:56 A full run (0:23)

On screen, slide 7. Three presses, the first steps, the whole run, the price.

A full run on this test grid is 56,000 of those steps, with the value
spreading out from the strike as we step back to today. It prices the option
at 195.64 on this coarse grid, and the full 2048 by 512 grid gives 196.17,
within a tenth of a cent of the Fourier price.

## 2:19 Checking the answer (0:13)

On screen, slide 8.

To validate the solver, we set the volatility of volatility to zero, where
Heston reduces to Black-Scholes, and the unchanged solver lands within 4
cents of the closed-form price.

## 2:31 The data structure (0:17)

On screen, slide 9. Three clicks: the rows flatten, the index of each cell, the stride to north.

Both buffers are flat arrays of doubles, with cell i, j at j times the row
length plus i. So the share price axis is contiguous in memory, and a step
in variance is a stride of one row.

## 2:48 How every number was measured (0:18)

On screen, slide 10.

Every number is the median of 5 runs on one exclusive rangpur node through
sbatch, after the tests pass in the same job. We discarded our
first sweeps because denormal numbers slowed some versions, and reran
everything with them flushed to zero.

## 3:06 Where the time goes (0:22)

On screen, slide 11.

Before changing anything we profiled with gprof, and one function, the
update of the interior rows, takes 99 percent of the time. The compiler alone
takes the baseline from 11 million updates a second at O0 to 58 at O2, and
every benchmark builds at O2, so each version measures only our change.

## 3:28 The ladder (0:13)

On screen, slide 12. The bars grow on their own. Point along them from left to right.

The ladder adds one technique per version, and each version must match the
baseline before it's timed, while the number in each bar is the gain over
the version below.

## 3:41 The serial rung (0:22)

On screen, slides 13 to 14. Point at a division, then press next. The code changes into the serial rung, then one more press grows the bar. On slide 14, press next for the bar.

The baseline divides by the grid spacing in every cell. The compiler can't
turn that into a multiplication by the reciprocal, because the rounding can
differ in the last bit, so we did it by hand, once per solve, for 159 million updates a
second, 2.8 times the baseline.

## 4:02 SIMD, four cells per instruction (0:15)

On screen, slide 15. Four clicks: the loads, then each group of four, then the scalar cell.

Serial was already at 96 percent of the scalar peak, so next we used vector
instructions. An AVX2 register holds four doubles, and neighbouring cells
use the same formula, so one instruction updates four cells.

## 4:17 SIMD code, before and after (0:17)

On screen, slides 16 to 17. Point at the loads, then press next. The serial loop turns into the AVX2 loop. Three more presses highlight the loads, the fused multiply-add and store, and the counter, and the last press grows the bar.

The serial loop turns into the AVX2 loop. Each scalar load
becomes a load of four doubles, the update becomes one fused multiply-add, and the counter
goes up by four. This runs at 456, 2.9 times serial.

## 4:34 The compiler was splitting the loads (0:28)

On screen, slide 18. Point at the two bars, then at the vinsertf128 lines.

Our first build ran this at 387, and the disassembly showed that g++ 8.5,
with only the AVX2 flag, splits every unaligned 256 bit load into two 128
bit loads. Two extra flags make each a single load, for 1.19 times the
speed, and every number here uses them. Given the same AVX2 flags, the
compiler's own version of the serial row only reaches 77 percent of ours.

## 5:02 SIMD on the roofline (0:13)

On screen, slide 19. The roofs draw on their own. Press next for serial, next for simd rising above it, next for the big grid.

On the one-core roofline, serial is at the scalar peak. SIMD reaches 14
gigaflops, a third of the AVX2 ceiling and just below the memory roof of
one core.

## 5:15 OpenMP, one block of rows per thread (0:16)

On screen, slide 20. Four clicks: the bands, the barrier, the boundaries, the swap.

With OpenMP each thread gets a contiguous block of rows. The threads read
the old buffer and write only their own rows, so there are no data races
and no locks. Then the main thread updates the boundaries.

## 5:31 OpenMP code (0:11)

On screen, slides 21 to 22. Press next and the pragma appears above the same loop. One more press grows the bar.

The code change is one pragma above the row loop. On 24 cores this gives
6367 million updates a second, 14 times SIMD.

## 5:42 Strong scaling (0:23)

On screen, slide 23.

With the grid fixed, 24 threads give 13.9 times, and Amdahl's law fits a
serial fraction of 2.4 percent. The boundaries are only 0.4 percent of a
one-thread step, so most of that is the fork and barrier every step and the
shared caches, while at 24 threads the boundaries grow to 9 percent.

## 6:05 Weak scaling (0:16)

On screen, slide 24.

For weak scaling each thread keeps 21 rows, so the grid grows with the
threads. The step time goes from 0.11 to 0.17 milliseconds, a scaled
speedup of 16 on 24 threads, below Gustafson's prediction of 22.

## 6:21 Grid size (0:16)

On screen, slide 25.

Small grids are dominated by the fixed cost of each step. While both
buffers fit in the two L3 caches, 24 threads run 40 times serial, and at
4096 by 1024, where they don't fit, the rate halves.

## 6:37 After threads, memory bandwidth (0:13)

On screen, slide 26. Press next for one core of simd, next for 24 threads rising, next for the big grid on the memory roof.

On the node roofline the small grid runs above the memory roof, because it
fits in L3. The big grid is on the memory roof, so it's
limited by memory bandwidth.

## 6:50 Placement (0:22)

On screen, slide 27. Point at the 4096 pair on the left, then the two pairs on the right.

On the big grid, 12 threads gave 1894 on one socket and 3413 on both.
Linux places each page in the memory of the socket whose thread first
writes it, and our main thread fills the grid, so every page is on socket 0.
Interleaving the pages made 24 threads 1.74 times faster.

## 7:13 MPI, one block of rows per process (0:18)

On screen, slide 28. Four clicks: the owned rows, the edge rows, the exchange, the wait.

With MPI each process owns a block of rows. Every step it sends its first
and last rows to its neighbours and receives theirs into its copied rows.
With one process per socket, each fills its own grid, so its pages are
local.

## 7:31 MPI code, the time loop (0:23)

On screen, slides 29 to 31. Press next and the loop changes into the MPI loop. Then one press per marked line, and a last press for the bars. On slide 31, three presses, the receive, the send, the wait.

In the loop, rows for rank gives each process its block, and it updates only
those rows and its own boundaries. The exchange is non-blocking, with MPI
Irecv posted before MPI Isend, and MPI Waitall blocks until both finish, then
at the end MPI Gatherv collects the blocks onto process 0.

## 7:54 Hybrid layouts (0:20)

On screen, slide 32.

One process per socket made the big grid 2.28 times faster than one
process of 24 threads. Two nodes give
1.21 times one node on the small grid and 1.64 on the big one, and on the
small grid the exchange takes 40 percent of each step.

## 8:14 Overlapping the exchange with the update (0:23)

On screen, slide 33. Point at the tall blue bar.

So we tried overlapping the exchange with computation, sending the edge
rows first and updating the rest while they're in transit. That only helped
with one process per node, at 1.41 times. We suspect that within a node
Open MPI waits for the receiver to call MPI Waitall before it transfers a
16 kibibyte row.

## 8:37 Same answer on every rung (0:13)

On screen, slide 34.

No version was timed until it matched the baseline. Every version is within
3.3 times 10 to the minus 15, and OpenMP and MPI match SIMD bit for bit.

## 8:50 Tried and rejected (0:15)

On screen, slide 35.

We tested these variants against OpenMP. Dynamic scheduling was slower,
because every row costs the same. Float was faster, but it moved
the price by 4 cents and gamma by 12 percent, so we rejected it.

## 9:05 Conclusion (0:17)

On screen, slide 36.

To conclude, we reached 123 times the baseline on one node and 147 on two.
Most of that is the result of SIMD and threads, and across nodes the
exchange is the main limit. Four of our six hypotheses held.

## 9:22 What surprised me, and what I would do next (0:16)

On screen, slide 37.

Looking back, I was surprised how much build and system settings changed
the results, from split loads to huge pages. Next I'd store only each
process's own rows, test our E9 explanation, and try a GPU version.

## 9:38 Close (0:10)

On screen, slide 38.

Every number here has a job id in the results file, and the repository has
the code to reproduce it. Thanks for watching.
