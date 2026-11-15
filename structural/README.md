# Structural Track 1 research

Resumed November 15, 2026 on Iana's explicit mandate: accept a higher risk of failure for a meaningful result and a more interesting research loop. The hackathon cutoff is secondary. Primary construction target: **227 valid length-16 two-deletion codewords**, beyond the currently checked public leader of 226. No larger code has yet been found. Claude retains Track 2 and team settings.

## Research contract

Every proposed program has a hypothesis, executable method, predicted observation and falsification test in `proposals.jsonl`. Low-conflict states are retained as search evidence; they are not valid codes. Source snapshots and hashes distinguish program revisions. `fixed_verifier.py` uses literal retained-position combinations and is locked by `verifier-lock.json`; it must not evolve with the search. The frozen official arena verifier provides an additional check at n=16.

The outer proposer at this stage is Codex in the research chat, inspecting results and writing program revisions. There is no hidden paid LLM proposer or claim that a fixed menu is general autonomous program discovery. The controlled controller comparison will use an explicitly bounded program grammar and common primitives so all baselines receive the same choices.

Initial local phase: **up to three CPU-hours**, at most three concurrent single-threaded search jobs, including exact repairs and related-instance controls. Checkpoint scientific evidence at that budget; do not silently extend it or infer an unlimited paid budget. No GPU or paid model call is authorized. Time limits are experiment budgets, not hackathon-driven stopping rules.

## Current findings

- **All 65,536 descendant sets audited.** Independent Python retained-position enumeration matches the C++ bit-deletion implementation for every word; see `descendant-audit.json`.
- **The 226-word parent's neutral component has exactly two states** under feasible single-word replacements/additions. There is only one one-blocker outsider in the initial state. This excludes that operator as an escape mechanism, not larger or infeasible paths. See `neutral-component.json` and `parent-diagnostic.json`.
- **Every two-parent recombination among 32 credited public/symmetry parents has optimum at most 226.** All 496 pair unions were solved exactly by bipartite matching, with matching/vertex-cover equality and independent code checks. This is a restricted finite-family result, not a global bound.
- **All outsiders with at most two parent blockers, simultaneously:** the 411-word pool cannot contain 227 codewords. This allows an arbitrary number of parent removals; it is not a radius-two exchange result. CP-SAT reported INFEASIBLE in about 0.04 seconds. See `runs/full-two-blocker-exchanges/`.
- **Four strongest distinct public parents:** their 698-word union cannot contain 227 (about 1.9 seconds with graph clique cuts). Four symmetry images of the 226-word parent and the global three-blocker pool are separate, larger trials.
- Initial fixed-cardinality search drifted from two collisions to hundreds. Initial feasible destroy/repair drifted from 226 to ~165. Both missed their predictions. A first guided implementation had a contradictory restart trigger that prevented penalty updates; it is excluded from algorithm comparisons. Its corrected revision did update penalties but still repeatedly returned to the same basin. These failures motivate new representations; more iterations alone are not treated as progress.

## Programs and evidence

`engine.cpp` supports feasible packing and fixed-cardinality conflict representations, including worse-state paths, bit/relocation moves, guided collision penalties, and standard (1,2) independent-set swaps. `crossover.py` supports exact pair-union recombination, larger parent unions and global blocker-pool repairs, with graph clique cuts beyond single received-word slot constraints. Binaries are local generated files. Raw scratch arrays stay out of Git.

The engine's incremental collision cost is checked against a full recount during runs. Only `best-valid.json` is eligible for the independent mathematical verifier. `best-state.json` and archive entries may intentionally contain conflicts.

## Controlled research question

Can a controller use local obstruction checks, residual-conflict structure and archive diversity to decide whether to continue a plateau, cross a barrier or change representations, and obtain better valid constructions for equal compute? Compare with random restarts, greedy selection and ordinary program evolution on development instances, then freeze choices before unseen related instances and multiple seeds. Primary metric is independently verified code size at a fixed CPU budget; conflict cost and diversity are secondary diagnostics. Include controller/verification overhead. The earlier six-minute autocorrelation plateau motivates this question but does not establish it.

Program evolution and archives have substantial prior art. Sources read:

- [CPro1](https://github.com/Constructive-Codes/CPro1), especially `designs/deletion-code/progDC1.c` and its README: fixed-cardinality tabu search, exact deletion outcomes, and published constructions; some original runs took up to 48 hours. Our code was written independently after inspecting that method.
- [AlphaEvolve](https://github.com/google-deepmind/alphaevolve): program evolution with local mathematical evaluation.
- [ShinkaEvolve](https://github.com/SakanaAI/ShinkaEvolve): open-ended/sample-efficient program evolution.
- [EinsteinArena problem definition](https://einsteinarena.com/api/problems/two-deletion-code-16), downloaded public parents and discussion 278: ARW-style local search, weight-band locality and earlier restricted-repair results. Public parents remain credited by source ID/agent in `public-parents.json`.

No novelty claim is made for those ingredients. Any contribution must be a new construction or controlled evidence about the specific controller and representations.
