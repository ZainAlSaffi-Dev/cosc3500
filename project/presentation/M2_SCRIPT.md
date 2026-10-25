# M2 SCRIPT, the word-for-word script for the Milestone 2 video

One section per take, in the order of the recorded slides of
`presentation/deck_m2.html` (slides 1 to 32). There are 19 takes, and each
one teaches one idea. It says what the idea is, shows what we measured, then
says why. A take that runs over more than one slide has a line in italics
where to press next for the next slide, and the notes under each slide
(press S) hold the words said over that slide. The appendix after slide 32
is for the interview and is never recorded. It starts with the eight slides
cut from the video, the project summary, the Black-Scholes check, the
ladder, the split loads, the one-core roofline, the MPI time loop and the
overlap experiment.

Each heading gives the start time and the length of the take, worked out
from its word count at 136 words per minute plus a second for each slide in
it. 136 is the pace of two timed run-throughs on 2026-10-25, which reached
the reflection at 9:45 to 9:50 with the earlier 1385 word script. The
Milestone 1 marker found the delivery fast, so the numbers still need to be
said slowly. The 1261 words come to 9:48, plus 2 more seconds holding the AI
declaration, and the assessors stop watching at 10 minutes. If a take runs
long when spoken, cut words from it and keep the pace.

The takes are written to be spoken. They use "we" for the work, as a paper
would, contractions where they sound natural, and the technical name for
each operation, so the code on screen and the words match.

Every number here comes from `docs/RESULTS.md`, Milestone 2, where each one
has its job id, and matches the numbers under `n` at the top of
`presentation/deck_m2.qmd`. The test grid numbers in the opening, the first
step value 0.04959375 and the price 195.6355 on screen, come from a local
run of the 421 by 51 test grid with `--dump-every`, so they are not
benchmark numbers. Numbers are written as digits so they are easy to read
off.

## Before recording

These come from the Milestone 1 feedback (`PLAN.md` §11a).

- No streaming overlay. The webcam sits in a bottom-right box that no
  slide uses, so it never covers a header or an axis.
- Turn on a large cursor and point at the number being said on every
  figure.
- Say the named results slowly. Amdahl is AM-dahl. Gustafson is
  GUS-taf-son. Heston is HESS-ton. NUMA is NOO-mah. AVX2 is A V X two.
- The diagrams in the opening build one step per key press, so press next
  as each new sentence starts.
- In a take over two or three slides, the italic line says when to press
  next for the next slide. Keep talking through the change.
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
- The roofline slide. S. Williams, A. Waterman and D. Patterson, Roofline,
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

## 0:00 Title (0:12)

On screen, slide 1.

Hi, I'm Zain. This video shows how we price an option on a grid, and how we
made it 123 times faster on one rangpur node.

## 0:12 What an option is (0:29)

On screen, slide 2. Five presses, the path above the strike, the path below
it, the payoff line, the other paths, the price.

A call option is the right to buy one share at a fixed price, here 5250, in
three months. If the share ends at 5700, you buy at 5250 and the option pays
450. If it ends below 5250, you don't buy, so it pays nothing. Nobody knows
the path, so the price today is the average payoff over every path,
discounted to today.

## 0:42 From the model to the PDE (0:29)

On screen, slide 3. Four presses, the share price terms, the volatility
terms, the last two, then the steps back from expiry.

In the Heston model the share price and its volatility both move randomly,
and the option value V follows this PDE. These terms come from the share
price, these from the volatility, and the last two say the two move together
and that money later is worth less. We only know V at expiry, the payoff, so
we step back from there to today.

## 1:11 Finite differences, then the same sum in C++ (0:42)

On screen, slides 4 to 5. The PDE moves to the top. Four presses, the slope,
the curvature, north and south, the corners.

To solve it we put V on a grid of share price by variance, and replace each
derivative with a difference of neighbouring cells. At 5250 on the first
step back, the slope is 0.5, and the curvature, how fast that slope changes,
is 0.02. Every row starts as the same payoff, so the other differences are
0.

*Slide 5. Press next for the slide. Three presses, the terms, the update,
the stability limit.*

In C++ each difference is multiplied by a reciprocal worked out once. One
step of dt gives 0.04959. A larger dt would make the errors grow, so the
full grid needs 980,000 steps.

## 1:53 One step over the grid, and how it sits in memory (0:41)

On screen, slides 6 to 7. Point at the payoff, the red block, then the red
cell. Six clicks, three cells, the rest of the row, the other rows, the
swap.

Every interior cell reads nine cells of the old buffer and writes one cell
of the new one. No cell reads a value written in the same step, so the cells
can be computed in any order, and that's what lets us split them up.

*Slide 7. Press next for the slide. Three clicks, the rows flatten, the
index of each cell, the stride to north.*

Each buffer is one flat array of doubles, with cell i, j at j times the row
length plus i. So a row of 2048 share prices sits together in 16 kibibytes,
and the cell to the north is one row further on.

## 2:34 A full run, and checking it (0:18)

On screen, slide 8. Three presses, 1 and 9 days, 37 days, today.

A full run on this test grid is 56,000 steps. The full grid gives 196.17,
within a tenth of a cent of the Fourier price, and with no volatility of
volatility the solver matches Black-Scholes to 4 cents.

## 2:52 How we measured, and where the time goes (0:30)

On screen, slides 9 to 10. Four presses, one box each of the rate sum.

A cell update is one cell's new value in one step. A million cells over 2000
steps is 2.1 billion updates, and in 13.2 seconds that's 159 million a
second. Every rate is the median of 5 runs on one rangpur node that no other
job shares.

*Slide 10. Press next for the slide. Point at the tall bar of the profile.*

gprof shows one function takes 99 percent of the time, so that's what every
version speeds up.

## 3:22 The serial rung (0:29)

On screen, slides 11 to 12. Point at a division.

The baseline divides by the grid spacing in every cell of every step, and a
division costs several times the cycles of a multiplication. The compiler
can't swap in the reciprocal by itself, because the rounding can differ in
the last bit.

*Slide 12. Press next and the code changes into the serial rung. One more
press grows the bar.*

So we worked out the reciprocals once per solve, for 159 million updates a
second, 2.8 times the baseline.

## 3:51 SIMD (0:40)

On screen, slides 13 to 15. Four clicks, the loads, then each group of four,
then the scalar cell.

Serial was already at 96 percent of the scalar peak, so next we went to
vector instructions. An AVX2 register holds four doubles, and neighbouring
cells run the same formula, so one instruction updates four cells.

*Slide 14. Press next for the code. Point at the loads, then press next. The
serial loop turns into the AVX2 loop.*

*Slide 15. Three presses highlight the loads, the fused multiply-add and
store, and the counter. The last press grows the bar.*

In the code each load now reads four doubles, the update is one fused
multiply-add, and the counter steps by four. That runs at 456, 2.9 times
serial. It's below 4 times because four cells need nine overlapping loads,
and one core is already near its memory bandwidth.

## 4:31 OpenMP (0:38)

On screen, slides 16 to 18. Four clicks, the bands, the barrier, the
boundaries, the swap.

The grid has 512 rows, one per variance value. OpenMP gives each of the 24
threads a block of about 21 neighbouring rows, so thread 0 updates rows 1 to
22, thread 1 the next 22, and so on. Each thread writes only its own rows,
so no locks are needed. After a barrier, the main thread updates the
boundaries.

*Slide 17. Press next for the code, then next again and the pragma appears
above the same loop.*

*Slide 18. One more press grows the bar.*

In the code that's one pragma above the row loop, and it gives 6367 on 24
cores, 14 times SIMD.

## 5:09 Strong scaling (0:38)

On screen, slide 19. Point at the 24 thread point, then the dashed Amdahl
curve.

Strong scaling keeps the grid the same and adds threads. Ideally 24 threads
would be 24 times faster. We measured 13.9. Amdahl's law explains the gap.
Any part that runs on one thread takes the same time however many threads we
add, so it caps the speedup. Fitting the law gives a serial fraction of 2.4
percent. The boundary update, the only code on one thread, is 0.4 percent.
The rest is the threads waiting at the barrier every step and sharing each
socket's cache.

## 5:48 Weak scaling (0:37)

On screen, slide 20. Point at the step time line, then the scaled speedup
against the Gustafson line.

Weak scaling gives each thread its own 21 rows, so 24 threads solve a grid
24 times bigger. Does the step time stay flat? The step time rose from 0.11
to 0.17 milliseconds, so 24 threads do 24 times the work in 1.5 times the
time, a scaled speedup of 16. Gustafson's law predicts 22 from our 9 percent
serial share. We fall short because one thread has a whole socket's cache
and bandwidth to itself, and 24 threads share them.

## 6:24 The roofline, and why grid size matters (0:38)

On screen, slides 21 to 22. Point at the omp line where it drops.

While both buffers fit in the L3 caches, 24 threads run 40 times serial. At
4096 by 1024 they no longer fit, and the rate halves.

*Slide 22. Press next for the slide. Press next for one core of simd, next
for 24 threads rising, next for the big grid on the memory roof.*

The roofline shows why. Speed is capped by how fast the cores compute, the
flat roof, or by how fast memory delivers data, the sloped roof. The small
grid runs above the memory roof because its data comes from L3. The big grid
sits on the memory roof, so after threads the limit is memory bandwidth.

## 7:03 MPI (0:35)

On screen, slides 23 to 24. Four clicks, the owned rows, the edge rows, the
exchange, the wait.

To go past one node we need MPI, and each MPI process has its own memory. So
each process gets a block of rows. Its first row needs the row below it,
which another process owns. So every step, each process sends its edge rows
to its neighbours and receives theirs into an extra row, called a ghost row.

*Slide 24. Press next for the code. Three presses, the receive, the send,
the wait.*

In the code we post the receive, then the send, and MPI Waitall waits for
both.

## 7:38 Memory placement, and why one process per socket wins (0:43)

On screen, slides 25 to 26. Point at the pages on socket 0 bar, then the MPI
bar.

MPI helped on one node too. Each socket has its own memory, and Linux puts a
page on the socket that writes it first. In OpenMP the main thread fills the
grid, so every page lands on socket 0. With one MPI process per socket, each
fills its own copy, so its pages stay local. That made the big grid 2.28
times faster.

*Slide 26. Press next for the slide. Point at the two node bars, then the
white exchange shares.*

Two nodes give 1.21 times on the small grid and 1.64 on the big one. The
small grid gains less because the exchange is 40 percent of each step.

## 8:20 Tried and rejected (0:15)

On screen, slide 27. Point at the dynamic bars, then the float bar.

Dynamic scheduling only added overhead, because every row costs the same.
Float ran 1.73 times faster but moved the price by 4 cents and gamma by 12
percent, so we rejected it.

## 8:35 The same answer, and the conclusion (0:25)

On screen, slides 28 to 29.

No version was timed until it matched the baseline, and every version agrees
to within 3.3 times 10 to the minus 15.

*Slide 29. Press next for the slide. Point down the ladder table.*

In the end we reached 123 times the baseline on one node and 147 on two.
Most of it came from SIMD and threads, and four of our six hypotheses held.

## 9:01 What surprised me, what I would change, and what comes next (0:39)

On screen, slide 30. Point at each block in turn.

What surprised me most was that two compiler flags alone gave SIMD 1.19
times. If I started again, I'd flush denormals before the first sweep. A
denormal is a double below about 2.2 times 10 to the minus 308, and this CPU
spends tens of extra cycles on every operation that touches one. Next, each
MPI process would hold only its own rows plus one ghost row from each
neighbour. That cuts each process from 16 to 4 mebibytes on two nodes, so
the grid can grow.

## 9:40 Close and AI declaration (0:08)

On screen, slides 31 to 32.

Every number here has a job id in the results file. Thanks for watching.

*Slide 32. Press next for the AI declaration and hold it for three seconds.
Nothing is said.*
