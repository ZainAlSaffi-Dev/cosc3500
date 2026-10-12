# M2 SCRIPT, the word-for-word script for the Milestone 2 video

One section per take, in the order of the beats in `PLAN.md` §11c. Each
heading gives the start time and the slot. The script is budgeted at 140
words per minute, slower than the 150 used for Milestone 1, because the
Milestone 1 marker found the delivery fast. The whole script is about 1300
words, which is 9 minutes 20 seconds at that pace, and the assessors stop
at 10. If a take runs long when spoken, cut words from it rather than
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

## 0:00 Cold open (0:20)

On screen, F1, the final bar highlighted.

One node of the rangpur cluster prices this option 101 times faster than my
first working solver, and two nodes go further. Every version gives the same
answer to 15 digits. This video shows the four techniques that got there,
and what each one was worth.

## 0:20 The problem (1:00)

On screen, `results/payoff.png`, then the weather map.

A European call is the right to buy a share at a fixed strike price on a
fixed day. On that day its value is known, it is the payoff. The question is
its price today. The Heston model lets the volatility of the stock move
randomly as well, and turns that question into a partial differential
equation in the stock price, the variance and time. I solve it on a grid,
with a stock axis and a variance axis. At expiry every cell holds the
payoff, and each time step walks one step back towards today, updating every
cell from its eight neighbours. That is the same formula on a million cells,
every step, which is what makes it a parallel problem.

## 1:20 The data structure (0:40)

On screen, `results/m2/code_shots/shot_data_layout.png`.

The grid is two flat arrays of doubles, current and next. The index
function puts a cell at row times row length plus column, so the stock axis
is contiguous and one step in variance is a jump of a whole row. At 2048 by
512 a row is 16 kilobytes, a sheet is 8 megabytes, and the two sheets are
16. Each step reads current, writes next, then the two swap. The loop nest
runs time outside, then variance rows, then stock columns inside, so the
inner loop walks memory in order.

## 2:00 Where the time goes (0:40)

On screen, F12, the gprof flat profile, then F11, the four -O levels.

Before changing anything I measured where the time goes. The gprof profile
puts 99.9 percent of the baseline in its one solve function, and over 98
percent of the serial and SIMD rungs in the row update. The boundaries are
under half a percent at one thread. So the hot spot is the interior update,
and every rung attacks it. The compiler alone takes the baseline from 11
million updates a second at O0 to 58 at O2, and O3 doubles that by
vectorising with two-lane SSE.

## 2:40 The serial rung (0:30)

On screen, `results/m2/code_shots/shot_serial_row.png`.

The serial rung keeps the one change from the first milestone that paid.
The baseline divides by the grid spacing five times per cell. The compiler
may not turn those divisions into multiplications, because that changes the
last bits of the answer, so I did it by hand, once per solve. That gives 159
million cell updates a second, 2.8 times the baseline.

## 3:10 SIMD (1:00)

On screen, `results/m2/code_shots/shot_simd_fma.png`, then F5.

The serial rung already runs at 96 percent of the scalar peak, so
the only way up is wider instructions. An AVX2 register holds four doubles,
and every cell in a row runs the same formula on neighbours that sit next to
each other in memory, so four cells can share each instruction. Each fused
multiply-add does a multiply and an add with one rounding. The compiler's own
vectoriser at O3 gets 1.5 times with two lanes. By hand, with four lanes and
FMA, the rung measures 379 million updates a second, 2.4 times serial. It is not 4, because each
group of four loads nine overlapping windows, and the cells left at the end
of a row run one at a time. On the roofline this rung does 12 gigaflops,
against an AVX2 ceiling of 43 for one core.

## 4:10 OpenMP (1:15)

On screen, `results/m2/code_shots/shot_omp_pragma.png`, then F2, then F7.

OpenMP shares the rows of each step between threads with one line. Each
thread reads only the old sheet and writes only its own rows, so no locks are
needed, and the answer is identical to one thread, bit for bit. On 24 cores
it gives 15.8 times. Fitting Amdahl's law gives a serial fraction of 1.2
percent, three times the boundary work, so most of the loss is the fork and
barrier on every step and the shared cache bandwidth. The scalar version
scales better, 18.5 times, because it moves data at half the rate. Placement
mattered as much as count. On the big grid, twelve threads on one socket
ran at 1880, the same twelve spread over both sockets at 3319, because they get two L3 caches
and two memory controllers. And every page of the grid lives on the socket
that first wrote it. Interleaving the pages gave the big grid 28 percent
more.

## 5:25 Hybrid MPI (0:50)

On screen, `results/m2/code_shots/shot_mpi_halo.png`, then F8.

MPI gives each process a slab of rows. Each step, every process swaps one
edge row with each neighbour, posting its receive before its send. On one
node I run one process of 12 threads per socket. Each process first touches its own copy of the grid on its own socket,
which fixes the page problem. That is 1.05 times one
24-thread process on the 2048 grid, and 2.16 times on the 4096 grid. Across
two nodes the best layout is 1.41 times one node, 70 percent efficient. The
halo is only 64 kilobytes, but each exchange costs about 45 microseconds of
latency and waiting, which is 38 percent of a step this short.

## 6:15 Same answer on every rung (0:25)

On screen, F9.

No version was timed until it matched the baseline. Every rung is within 3.3
times 10 to the minus 15 on price, delta, gamma and vega, and the OpenMP and
MPI rungs must equal the SIMD rung exactly. Every benchmark job reran these
tests before timing anything.

## 6:40 Weak scaling and grid size (0:45)

On screen, F3, then F4.

Weak scaling keeps the work per thread fixed and grows the grid with the
threads. The step time rises from 0.14 to 0.18 milliseconds, a scaled
speedup of 19 on 24 threads. Grid size decides the regime. Small grids lose
to the fixed cost of each step. The 2048 grid fits in the two L3 caches and
runs 36 times faster than serial. At 4096 by 1024 the sheets spill to
memory, the rate halves, and the code sits on the measured memory roof.

## 7:25 Bottlenecks (0:35)

On screen, F6, then F5.

This is where each step's time goes. On one thread nearly all of it is the
interior rows. At 24 threads the serial boundaries grow to 8 percent of the
step, because the main thread has to fetch rows other cores just wrote. That
overturned my guess of half a percent. Across nodes, the halo exchange takes
over. The roofline puts the large grid on the memory roof, so after threads
the limit is bandwidth, not arithmetic.

## 8:00 Tried and rejected (0:30)

On screen, F1, the right panel.

Each bar on the right is a variant measured against the OpenMP rung. Dynamic
scheduling lost up to 40 percent, because every row costs the same. Float
instead of double was 1.75 times faster, and 4.5 times on the big grid
because it fits back in cache. I rejected it, because on the accurate grid
it moved the price by 4 cents and gamma by 12 percent.

## 8:30 Conclusion (0:40)

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

## 9:10 Close (0:10)

On screen, F1.

Every number has a job id in the results file, and the repository has the
code and commands to reproduce it. Thank you.
