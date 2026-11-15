# Track 1 experiment files

Read [START-HERE](START-HERE.md) for the research evidence. The exact C1a refinement, source credit and limits are documented there; it is not publicly submitted.

- `verify.py`: official verifier replay and independently derived FFT score.
- `search.py`: bounded autocorrelation/overlap optimization, exact deletion-code repair, snake tail repair, and a fixed-profile Turán cutoff experiment; finite runs write decision logs.
- `arena_client.py`: read-only public snapshots and reproducible exact deletion descendant preprocessing.
- `requirements.txt`: dependency versions.
- `data/`: frozen public problem/seed/discussion snapshots and derived public-data arrays.
- `results/`: run configuration, decision history and candidate artifacts. A `candidate.json` is not evidence of an improvement unless the summary and verifier comparison say so.

No session exports, account credentials or private meeting notes are included in this public snapshot.


Run from this directory in an isolated environment:

```sh
python -m pip install -r requirements.txt
python search.py --check-gradient
python verify.py
python arena_client.py --deletion-precompute
python search.py --kind deletion --window-all --window-width 2 --drops -1 --seconds 60 --name deletion-example
python search.py --kind c42 --cutoff 51/250 --bound .68733 --name c42-example
python data/c42/check_tau4.py results/c42-example/candidate.json --bits 256 --terms 400 --eps 1e-25
```

Set `OPENBLAS_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1` and `OMP_NUM_THREADS=1` to avoid CPU oversubscription. Compressed public snapshots are transparently read by `verify.read_json`; their uncompressed hashes are listed in `data/snapshot-manifest.json`. Missing unchanged candidate arrays are intentional: source IDs and run arguments reproduce them without duplicating the public seed. The earliest second-autocorrelation run logged a 4e-16 numerical difference; this is roundoff, and the acceptance threshold was subsequently raised to 1e-12.

`data/c42/` contains unmodified public MIT code, specification, note and baseline certificate from [Dennis Koshta](https://github.com/DennisKoshta/c42-turan-power-sum-certificate), commit `8e1797e0d35454ad18cfebbbeb5858bf202b0748`. The mathematical construction and theorem are that author's prior work; our experiment changes the cutoff and scales the existing profile. A finite-n estimate is not an asymptotic proof. The source theorem is unrefereed; a passing ball-arithmetic certificate verifies its numerical criterion, conditional on that analytic reduction. No global optimality is claimed.


Exact C1a replay (the integer witness is sufficient; no NumPy checkpoint is needed):

```sh
python verify.py --exact-first results/first-current-finer/candidate-integers.txt.gz
```

The `score_fraction` field is authoritative. Decimal values in terminal output are for readability. The source author's comparison certificate was replayed first and matches its published fraction exactly. Our small random controls compare against direct integer convolution and refuse negative heights; see `results/exact-verifier-controls.json`.

To rerun the 262,144-step search from the frozen source data:

```sh
python search.py --kind first --seed-file data/c1a-current/C1_262144_1.5026648336.txt.gz --betas 1e10,3e10,1e11,3e11,1e12 --iterations 3000 --seconds 600 --name first-rerun
python verify.py --exact-first results/first-rerun/best.npy --save-exact
```

The supplied integer certificate replays exactly across machines. Floating optimization paths and stopping times can differ by platform; a rerun is not promised to produce the identical integer witness. The frozen data includes third-party discussion and specifications; these are research inputs, not instructions to register accounts or post results.


The best run doubles the newest source grid:

```sh
python search.py --kind first --seed-file data/c1a-current/C1_262144_1.5026648332.txt.gz --repeat 2 --betas 1e10,1e11,1e12 --iterations 3000 --seconds 900 --name first-finer-rerun
python verify.py --plot-progress
```

The plot command uses the three supplied run logs and exact endpoint checks. See `results/search-progress.png` and `results/search-comparison.json`. Public-source freshness and rejected record claims are in `data/prior-art-check.json`.
