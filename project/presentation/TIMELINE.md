# TIMELINE — the minute-by-minute plan for the video

This file turns the beat sheet in PRESENTATION.md into timestamps and
concrete assets. The beat sheet decides what to say and why it earns marks.
This file is the thing to keep next to the recording software.

The slides are `presentation/deck.pdf`, built from `deck.tex` by
`make slides`. There is one slide for each row of the table below, in the
same order. The backup slides after the closing frame are interview
material and are not recorded. `deck.html` is the earlier HTML deck. Its
take timer is still useful, but its slide order no longer matches this
table, so record over the PDF.

Every slot length is computed from the take's word count in SCRIPT.md at
150 words per minute, so a take read at a normal pace fills its slot with
nothing left over. The table ends at nine minutes forty and the assessors
stop watching at ten. If a take runs long, re-record that take rather than
letting the drift accumulate.

| start | length | what happens | what is on screen | camera |
|---|---|---|---|---|
| 0:00 | 0:20 | The opening. Say the three promises of the video. | `results/demo_timevalue.mp4` looping behind the title | face large in the corner |
| 0:20 | 0:20 | What an option is. | `results/payoff.png` | small webcam |
| 0:40 | 0:20 | Why Black-Scholes is not enough. The market prices a smile. | `results/xi_smile.mp4`, then hold its final frame | small webcam |
| 1:00 | 0:15 | Heston in one sentence, and its two equations. | a text slide | small webcam |
| 1:15 | 0:15 | Every symbol in those equations, walked in one pass. | a table slide | small webcam |
| 1:30 | 0:20 | From stochastic model to deterministic equation. There is no random number generator in the program. | a text slide | small webcam |
| 1:50 | 0:10 | The derivative notation. Capital V and lowercase v are different quantities. | a table slide | small webcam |
| 2:00 | 0:20 | Why the noise terms become second derivatives. | a text slide | small webcam |
| 2:20 | 0:15 | The equation term by term. Each part of the model contributes one term. | a table slide | small webcam |
| 2:35 | 0:30 | The full equation, known numbers separated from the one unknown, and where the strike went. | a text slide | small webcam |
| 3:05 | 0:30 | The grid and its spacing, filled in backwards from expiry. | `results/demo_timevalue.mp4` again, this time explained | small webcam |
| 3:35 | 0:20 | Validation begins. Switch the model off and it becomes one we can check exactly. | `results/bs_collapse.mp4` | small webcam |
| 3:55 | 0:15 | Three methods that share no code agree on one number. | a table slide with the three prices | small webcam |
| 4:10 | 0:20 | Why the optimisation study is a ladder rather than a before and after. | a text slide | small webcam |
| 4:30 | 0:10 | Level 0, the harness control. One sentence, keep moving. | `results/code_shots/shot_L0_anchor.png` | small webcam |
| 4:40 | 0:10 | Level 1, hoisting. The compiler had already done it. Keep moving. | `results/code_shots/shot_L1_hoisting.png` | small webcam |
| 4:50 | 0:35 | Level 2, strength reduction. The change the compiler is not allowed to make. | `results/code_shots/shot_L2_strength_reduction.png` | small webcam |
| 5:25 | 0:15 | Levels 3 and 4 in one pass. The baseline already did both, so throughput did not change. | `shot_L3_order.png` and `shot_L4_split.png` side by side | small webcam |
| 5:40 | 0:20 | Level 5, induction variables. The one small real win after level 2. | `results/code_shots/shot_L5_pointers.png` | small webcam |
| 6:00 | 0:25 | Level 6, unrolling, and the laptop regression that justifies the cluster. | `results/code_shots/shot_L6_unroll.png` | small webcam |
| 6:25 | 0:35 | The ladder, measured. Only level 2 moved, the swapped loop order shows what the traversal order is worth, and the test holds every level to the answer. | `results/bench_ladder_hugepage.png` | small webcam |
| 7:00 | 0:30 | How the benchmarks were set up so they could not mislead, and the control that misled anyway. | the protocol table, with `results/code_shots/shot_csv_header.png` spliced fullscreen | small webcam |
| 7:30 | 0:20 | How performance scales with problem size. | `results/bench_scaling.png` | small webcam |
| 7:50 | 0:25 | The bottleneck. The solver was assumed to be memory-bound and is not. | `results/roofline.png` | small webcam |
| 8:15 | 0:30 | The page-size effect. | a table slide with the page-size numbers | small webcam |
| 8:45 | 0:35 | Reflection. The four questions the spec asks, one sentence each. | a text slide | face large |
| 9:20 | 0:20 | Close by returning to the promises from the opening. | face only | face full screen |

The control that went wrong is told inside the protocol take. My
no-vectorise sweep at -O2 disabled nothing and proved nothing, because GCC
8.5 does not vectorise at -O2 at all. Re-run properly at -O3, the answer
split in two. The vectoriser does nothing to my optimised kernels, but it
nearly doubles the baseline, because division-heavy code is exactly what
it likes. That story is worth its seconds because the spec explicitly
rewards a benchmarking setup that avoids misleading results, and catching
your own control is direct evidence you look for them. The numbers behind
it are in RESULTS.md section 6.

## What was cut, and why

| asset | reason |
|---|---|
| the strike-price slide | Its point survives as one sentence in the full-equation take. The slide is in the backup section. |
| `results/convergence.png` | The convergence slide moved to the backup section. Validation is carried by the Black-Scholes collapse and the three-method agreement. |
| the level 3 and 4 slides | Levels 3 and 4 share one slide in the recording, with both code screenshots side by side. Their full per-level slides are in the backup section for the interview. Levels 0 and 1 kept their own slides with short takes. |
| `results/code_shots/shot_optmatches.png` | The speedup-cost slide moved to the backup section. Its number, 15 units in the last place under an automated test, is one sentence on the ladder slide. |
| `results/code_shots/shot_iron_rule.png` | The test-gate slide came out. Its content is one row of the protocol table and one spoken sentence in that take. |
| `results/code_shots/shot_ctl_order.png`, `results/memory_cache.mp4` | The negative-controls slide moved to the backup section. Its headline number, the 2.1x loop-order cost, is one sentence on the ladder slide. |
| `results/demo_surface.mp4` | It does the same job as the flat weather map, and the three-dimensional view is slower for a viewer to parse in the time available. |
| `results/memory_buffers.mp4` | The two-buffer idea lands in one spoken sentence, so fifteen seconds of animation would not buy fifteen seconds of marks. |
| `results/blowup.mp4` | The instability demonstration is an aside. PRESENTATION.md keeps it listed as a spare in case a recording comes in short. |
| the older `weather_map` and `surface` renders | They are earlier versions of what `demo_timevalue.mp4` now does better. |
| every `.gif` file | The GIF files exist as twins of the MP4 files for embedding in documents. The video always uses the MP4. |

The rule behind every cut is that each second on screen has to serve one of
the five marked criteria, and the spec says plainly that a shorter
high-quality presentation beats a longer low-quality one.
