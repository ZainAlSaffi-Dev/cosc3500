# Graph Report - /Users/zer0/Documents/cosc3500/project  (2026-08-07)

## Corpus Check
- Corpus is ~0 words - fits in a single context window. You may not need a graph.

## Summary
- 28 nodes · 37 edges · 4 communities
- Extraction: 92% EXTRACTED · 8% INFERRED · 0% AMBIGUOUS · INFERRED: 3 edges (avg confidence: 0.85)
- Token cost: 42,714 input · 5,200 output

## Community Hubs (Navigation)
- Milestone Deliverables
- Benchmarking & Clusters
- Parallel Techniques
- Rubric & Submission

## God Nodes (most connected - your core abstractions)
1. `Shared Marking Rubric` - 11 edges
2. `Parallel Implementation` - 7 edges
3. `Milestone 2 (Parallel Presentation)` - 6 edges
4. `Milestone 1 (Serial Presentation)` - 4 edges
5. `In-Person Interview (Milestone 2 only, pass/fail)` - 4 edges
6. `Benchmarking (25%)` - 4 edges
7. `COSC3500 Milestone 1 / Milestone 2 Specification` - 3 edges
8. `10-Minute Video Presentation` - 3 edges
9. `Identity Verification (pass/fail)` - 3 edges
10. `Code Submission (pass/fail)` - 3 edges

## Surprising Connections (you probably didn't know these)
- `Serial Implementation` --shares_data_with--> `Parallel Implementation`  [INFERRED]
  COSC3500_Milestone.pdf → COSC3500_Milestone.pdf  _Bridges community 0 → community 2_
- `Milestone 1 (Serial Presentation)` --references--> `Shared Marking Rubric`  [EXTRACTED]
  COSC3500_Milestone.pdf → COSC3500_Milestone.pdf  _Bridges community 0 → community 3_
- `Shared Marking Rubric` --references--> `Benchmarking (25%)`  [EXTRACTED]
  COSC3500_Milestone.pdf → COSC3500_Milestone.pdf  _Bridges community 3 → community 1_
- `Shared Marking Rubric` --references--> `Optimisation (25%)`  [EXTRACTED]
  COSC3500_Milestone.pdf → COSC3500_Milestone.pdf  _Bridges community 3 → community 2_

## Hyperedges (group relationships)
- **Marking Rubric Criteria** — cosc3500_milestone_marking_rubric, cosc3500_milestone_introduction_and_background, cosc3500_milestone_optimisation, cosc3500_milestone_benchmarking, cosc3500_milestone_presentation_criterion, cosc3500_milestone_reflection_conclusion, cosc3500_milestone_code_submission, cosc3500_milestone_identity_verification, cosc3500_milestone_in_person_interview, cosc3500_milestone_submission_format [EXTRACTED 1.00]
- **Parallel Programming Technique Options for Milestone 2** — cosc3500_milestone_parallel_implementation, cosc3500_milestone_avx, cosc3500_milestone_openmp, cosc3500_milestone_mpi, cosc3500_milestone_cuda [EXTRACTED 1.00]

## Communities (4 total, 0 thin omitted)

### Community 0 - "Milestone Deliverables"
Cohesion: 0.31
Nodes (9): Course FAQ, CPU vs GPU Comparison Approach, Identity Verification (pass/fail), In-Person Interview (Milestone 2 only, pass/fail), Milestone 1 (Serial Presentation), Milestone 2 (Parallel Presentation), Serial Implementation, COSC3500 Milestone 1 / Milestone 2 Specification (+1 more)

### Community 1 - "Benchmarking & Clusters"
Cohesion: 0.33
Nodes (7): Amdahl's Law, Benchmarking (25%), Code Submission (pass/fail), Gustafson's Law, Makefiles, Slurm Scripts, UQ Clusters

### Community 2 - "Parallel Techniques"
Cohesion: 0.33
Nodes (6): AVX, CUDA, MPI, OpenMP, Optimisation (25%), Parallel Implementation

### Community 3 - "Rubric & Submission"
Cohesion: 0.33
Nodes (6): H.264 Video Format, Introduction and Background (20%), Shared Marking Rubric, Presentation (20%), Reflection/Conclusion (10%), Submission Format (pass/fail)

## Knowledge Gaps
- **9 isolated node(s):** `Introduction and Background (20%)`, `Presentation (20%)`, `Reflection/Conclusion (10%)`, `AVX`, `OpenMP` (+4 more)
  These have ≤1 connection - possible missing edges or undocumented components.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Shared Marking Rubric` connect `Rubric & Submission` to `Milestone Deliverables`, `Benchmarking & Clusters`, `Parallel Techniques`?**
  _High betweenness centrality (0.693) - this node is a cross-community bridge._
- **Why does `Parallel Implementation` connect `Parallel Techniques` to `Milestone Deliverables`?**
  _High betweenness centrality (0.294) - this node is a cross-community bridge._
- **Why does `Milestone 2 (Parallel Presentation)` connect `Milestone Deliverables` to `Parallel Techniques`, `Rubric & Submission`?**
  _High betweenness centrality (0.277) - this node is a cross-community bridge._
- **What connects `Introduction and Background (20%)`, `Presentation (20%)`, `Reflection/Conclusion (10%)` to the rest of the system?**
  _9 weakly-connected nodes found - possible documentation gaps or missing edges._