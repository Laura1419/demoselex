# Persistent ELAI long-run driver

[observation] `analysis/elai_amy_long.py` builds once, fits in separate processes, and reports gates. The full dense build persists in `data/derived/elai_amy_dense/work_long/`. **No 200-step run was launched. The hours-long runs below are the parent's to launch after review.** No fork changes, installations, raw-data edits or git commands were made.

## Approved export correction and limits

[observation] Full density gives `1,000 × 54,370 × 25 = 1,359,250,000` ordered-joint TSV rows per K5 fit, regardless of the EM step count. Even an optimistic 40 bytes/row gives 54.37 GB of pipe traffic. A three-step full joint dump therefore is not a cheap end-to-end test. The supervisor approved full-density building/short marginal verification plus a small deterministic fixture for complete joint-FIFO, parallelism and gate verification; it explicitly withdrew fork-patch work from this lane.

[observation] `--pkk-site 1:104121000` already operates without `--pkk`. Long runs use that site export **and `--pk1`**, not a full joint dump. `--pk1` exposes per-marker **lower-cluster expected diploid copy counts**, flattened locus-major then cluster-major, one line per individual, rounded to three decimals. The consumer divides by **2**, averages over all markers for composition, and also averages alternate markers for the thinning diagnostic. It does not silently renormalize rounded values. The marginal stream is itself large: 271,850 values per individual, **271,850,000 values per K5 fit**, roughly 1.6 GB of transient text. It is drained through a FIFO as well; neither marginal nor full joint tensors are persisted.

[observation] The narrow supporting change is `dense.fit(..., full_export=True)`: its default command and joint-export behavior remain unchanged. The driver passes `full_export=False` for marginal mode, and sets **`elai_amy_dense.STEPS = args.steps` in that process** before calling `fit`; the source constant and preregistration are not edited. Calibration is disabled for these explicitly wall-capped background fits. The site posterior is captured before the upstream `--pk1` output-driven admixture update; fixture equivalence tests check this, not just the flag spelling.

[inference] Marginal compositions can support the known-structure and stability diagnostics, but **they do not establish preregistered complete ordered-joint coverage**. That component remains **NOT TESTED**, with an export-size-driven reason, not a scientific failure or a relaxed threshold. Frozen linkage-integration replays and independently refitted thinning are also NOT TESTED by this driver. Consequently these three long fits alone cannot authorize biological use: even with convergence and otherwise passing metrics, missing mandatory checks yield PARTIAL, not PASS.

## Build and persistent inputs

[observation] All commands are run from `/home/erikg/demoselex` with the existing `.venv/bin/python`. Threading is fixed to one BLAS/OpenMP thread before numerical imports. Modes write wall-clock/start times and JSON summaries; provenance uses `scripts/provenance.py` with git probes disabled before recorder construction. Builds/fits use filesystem locks, fits sharing input take read locks, and identical K/seed prefixes are refused rather than overwritten. Distinct K/seed processes share the same input read-only. Gate mode is run **after** fits finish.

```sh
.venv/bin/python -u analysis/elai_amy_long.py --mode build
```

[observation] The completed build took **92.876 s**, retained **54,370 markers**, fitted **N=1,000**, and retained **9,854 independent PC markers**. Actual sizes, rather than the preliminary ~110 MB estimate:

| Retained input | Bytes |
|---|---:|
| `panel.inp` | 163,966,772 |
| `panel.pos` | 1,476,136 |
| `genotypes.txt` | 163,956,808 |

[observation] The original builder deletes its genotype intermediate. The driver reconstructs exactly the genotype lines from `panel.inp` after its three BIMBAM headers and **retains** them in the persistent work directory. `panel.json`, `samples.json`, `markers.json`, and `independent_pcs.npy` serialize the panel, sample list, marker list and PC panel; source-named TSV summaries are retained too. No rebuilding is done in fit or gate mode. Existing builds are preserved/refused on re-build; use a fresh work directory for a new experiment.

## Exact three long-launch commands

[observation] **Timing warning from this lane's new full-density verification:** measured slope was **69.805 s/step**, with **248.459 s** estimated fixed initialization/export overhead. [inference] This projects **14,209.372 s = 3.947 h** for K5 at 200 steps, leaving only **790.628 s (~13.2 min)** below the requested 15,000-s cap. Keep the approved commands below, but parallel contention can exceed the cap; the earlier 61-s estimate is not a guarantee. At this newer slope, parallel wall is roughly **3.95 h**, if sufficient independent compute is available. These remain background jobs, not 35-minute-lane jobs.

[observation] The paths below exist; the three canonical output prefixes are unused by verification. The noncanonical verification seed 20261012 was used at full density specifically to preserve them. Launch these **three commands**, from the repository root, to run in parallel and survive terminal logout:

```sh
nohup env PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -u analysis/elai_amy_long.py --mode fit --work data/derived/elai_amy_dense/work_long --k 5 --seed 20261010 --steps 200 --wall 15000 --export marginal > data/derived/elai_amy_dense/work_long/panel_k5_seed20261010.launcher.log 2>&1 &
nohup env PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -u analysis/elai_amy_long.py --mode fit --work data/derived/elai_amy_dense/work_long --k 5 --seed 20261011 --steps 200 --wall 15000 --export marginal > data/derived/elai_amy_dense/work_long/panel_k5_seed20261011.launcher.log 2>&1 &
nohup env PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -u analysis/elai_amy_long.py --mode fit --work data/derived/elai_amy_dense/work_long --k 3 --seed 20261010 --steps 200 --wall 15000 --export marginal > data/derived/elai_amy_dense/work_long/panel_k3_seed20261010.launcher.log 2>&1 &
```

[inference] Using the task's measured **~61 s/EM step** at M=54,370, N=1,000, K=5: `200 × 61 = 12,200 s = 3 h 23 min 20 s = 3.389 h` **per K5 fit, EM only**. This is not 3.1 h: the previous lane's 55.925 s/step gives 3.107 h, explaining that approximation. For K3, quadratic scaling projects `61 × (3/5)^2 = 21.96 s/step`, hence `200 × 21.96 = 4,392 s = 1.220 h`, EM only. Initialization and posterior/export passes add fixed overhead. With adequate independent cores, the three parallel fits take approximately the maximum, **3.39 h plus overhead**, not their summed EM cost of 7.998 h. Sequential execution would take that summed cost plus overhead. The **15,000 s = 4.167 h** wall cap is per fit and includes export. Contention and late-step cost changes can invalidate these projections; completion and convergence are evaluated, not assumed.

[observation] For prefix `P=panel_k{K}_seed{S}`, expected persistent work-root files are `P.stdout.txt`, `P.run.json`, `P.summary.json`, `P.composition.tsv`, `P.thin_composition.tsv`, `P.pkk-site.tsv`, `P.snpinfo.txt`, and the redirected `P.launcher.log`. Binary files live in `output/`: `P.em.txt`, `P.admix.txt`, `P.ps21.txt`, `P.log.txt`, `P.snpinfo.txt`, and `P.pkk-site.tsv`. `output/P.pk1.txt` is a **FIFO, not a retained data table**; do not `cat` it after completion (that would block). Joint fixture mode instead leaves a FIFO named `output/P.pkk.tsv`. Prefix locks are separate; no fits share output files. The upstream success status is **1 plus `-elai: finished`**; the successful driver returns **0**. EM trace, return code, exact command and arithmetic diagnostics are in each run JSON.

[observation] After all three launcher processes exit, run:

```sh
.venv/bin/python -u analysis/elai_amy_long.py --mode gates --work data/derived/elai_amy_dense/work_long
```

[observation] This writes `gate_metrics.json` and `gates.summary.json`. It reuses dense gate formulas, applies the preregistered strict NMI inequalities, requires convergence before accepting any labels, and reports the five gates plus mandatory EM convergence. Country adjusted mutual information is the singleton-robust composition-assignment diagnostic, with no invented passing threshold. Stability reports seed, cross-K, virtual argmax and independent thinning-refit statuses separately; missing tests stay NOT TESTED. Full-versus-alternate-marker compositions are only a diagnostic, not a refit. Any numerical failure gives FAIL; otherwise missing tests give PARTIAL. Neither permits biological claims.

## Verified cheap end-to-end pipeline

[observation] The permanent `work_fixture/` is a deterministic real-input subset: first **12** selected individuals, **80** evenly spaced retained markers plus both core-gap flanks, for **82** total markers. It uses the corresponding five independent-PC rows of the full build; it is not a newly fitted PC panel. Its small N/M and only three EM steps validate orchestration/I/O/arithmetic, **not dense-scale scientific gates**. The full build remains unchanged. A second `work_fixture_marginal/` verifies the marginal driver path on identical inputs. These exact verification commands were executed:

```sh
.venv/bin/python -u analysis/elai_amy_long.py --mode build --fixture --work data/derived/elai_amy_dense/work_fixture
.venv/bin/python -u analysis/elai_amy_long.py --mode build --fixture --work data/derived/elai_amy_dense/work_fixture_marginal
.venv/bin/python -u analysis/elai_amy_long.py --mode fit --work data/derived/elai_amy_dense/work_fixture --k 5 --seed 20261010 --steps 3 --wall 900 --export joint
.venv/bin/python -u analysis/elai_amy_long.py --mode fit --work data/derived/elai_amy_dense/work_fixture_marginal --k 5 --seed 20261010 --steps 3 --wall 900 --export marginal
.venv/bin/python -u analysis/elai_amy_long.py --mode gates --work data/derived/elai_amy_dense/work_fixture
.venv/bin/python -u analysis/elai_amy_long.py --mode gates --work data/derived/elai_amy_dense/work_fixture_marginal
```

[observation] For a repeat after the real run, choose fresh names under `data/derived/elai_amy_dense/` with `--work`; the fixture build always reads the retained full `work_long/`. Existing prefixes deliberately refuse reruns. The checked-in tests can be rerun using the permanent verification fixtures without changing their fits.

[observation] Fixture primary: driver success, binary exit **1**, finished sentinel, wall **0.502 s**, completed **3/3** steps. The driver log/run JSON records this exact invocation; binary stdout independently echoes seed, dimensions and all exact updates:

```text
.../elai -g .../work_fixture/panel.inp -p 1 -pos .../work_fixture/panel.pos -C 1 -c 5 -w 0 -s 3 -mg 10 -R 20261010 -o panel_k5_seed20261010 --pkk-site 1:104121000 --pkk
## randseed = 20261010
-warning: position files contain 82 records.
-warning: file 0 has 12 individual and 82 snps
 ### Exact update 0 0  -497.117 -497.117
 ### Exact update 0 1  -452.601   44.516
 ### Exact update 0 2  -403.651   48.950
-elai: finished, for details see log: panel_k5_seed20261010.log.txt
```

[observation] The joint FIFO contained and consumed **12 individuals × 82 markers × 25 = 24,600 ordered rows**; completeness=True. Maximum ordered-sum error over all 984 individual/locus combinations was **7.771561172376096e-16**, checking the patch's sum-to-one convention. Site ordered-sum error=1.1102230246251565e-16, q symmetry=0, triangular q-sum error=2.220446049250313e-16; finite probabilities in [0,1]. The full three-step convergence trace is retained in JSON; terminal relative increment=0.10815265542939584, so convergence **FAIL**, legitimately. All outputs landed in the persistent fixture work directory.

[observation] Gates returned driver exit **0**, verdict **FAIL**. Known structure=NOT TESTED (nonconverged labels rejected); stability components=NOT TESTED, overall PARTIAL; usable input=PASS (12/12), full ordered export=PASS (24,600 rows), dated complete export=FAIL (12<100); virtual arithmetic=PASS but integration=NOT TESTED, overall PARTIAL; exclusion=PASS (246 actual retained markers audited across three fits, zero core markers); EM convergence=FAIL. Marginal fixture gates also execute, with full ordered coverage explicitly NOT TESTED. A separate isolated regression simulates convergence only in temporary test artifacts to exercise known-structure/country-AMI and seed/cross-K calculations; it does not relabel the real failed fits.

### Parallel verification

[observation] These fits were launched together against **the same built fixture**, then both were waited for:

```sh
.venv/bin/python -u analysis/elai_amy_long.py --mode fit --work data/derived/elai_amy_dense/work_fixture --k 5 --seed 20261011 --steps 3 --wall 900 --export joint > data/derived/elai_amy_dense/long_parallel_k5.log 2>&1 & A=$!
.venv/bin/python -u analysis/elai_amy_long.py --mode fit --work data/derived/elai_amy_dense/work_fixture --k 3 --seed 20261010 --steps 3 --wall 900 --export joint > data/derived/elai_amy_dense/long_parallel_k3.log 2>&1 & B=$!
wait "$A"; wait "$B"
```

[observation] Both driver exit codes=**0**, both binary return codes=**1**, both complete=True. Their actual fit intervals overlapped by **0.48499155044555664 s**. K5 consumed **24,600** joint rows; K3 consumed **8,856 = 12×82×9**. Distinct prefixes produced distinct FIFOs, logs, parameters, site exports and compositions. `long_parallelism.json` records timestamps/completeness. Reusing the completed primary prefix was separately rejected with driver exit **2** before overwriting any result (`long_collision_guard.log`).

### Marginal equivalence and regression evidence

[observation] On the same 12×82 fixture/seed/model, live binary fits produced joint+pk1, marginal+site, and standalone site-only outputs. The test sums ordered joint probabilities over k2 and reconstructs marginals from unordered q (diagonal plus half each off-diagonal unordered pair). Results (`long_equivalence.json`):

| Check | Observed maximum error/result |
|---|---:|
| joint marginal vs marginal-mode pk1/2 | 0.00024996147472021446 |
| allowed pk1 rounding error | 0.000250000001 (0.0005/2 + 1e-12) |
| ordered vs q-derived marginals | 1.1102230246251565e-16 |
| ordered sum-to-one | 6.661338147750939e-16 |
| joint driver composition vs joint-derived mean | 3.0531133177191805e-16 |
| marginal driver composition vs marginal-mode mean | 0 |
| site-only / marginal / joint-enabled site files | byte-identical |
| existing em, admixture, ps21, snpinfo, pk1 outputs with/without joint export | byte-identical |

[observation] **20 tests passed** (six new long-driver regressions plus 14 existing dense/baseline regressions), including default command equality, copy-count halving, truncated marginal rejection, gate missing-test aggregation, converged metric path, and live equivalence:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=analysis .venv/bin/python -m unittest analysis/test_elai_amy_long.py analysis/test_elai_amy_dense.py analysis/test_elai_amy.py -v
```

```text
Ran 20 tests in 0.259s
OK
```

[observation] Supporting dense-module diff: **33 added / 6 removed lines**; reconstructed original source SHA256 matches the previous dense-run manifest. New driver=295 lines; new tests=135 lines (`long_diff_summary.json`). No fork patch was made. Independent parent review remains required.

## Full-density short verification and disk evidence

[observation] A further full-density, three-step verification used a noncanonical seed, so the long-launch prefixes remain available:

```sh
.venv/bin/python -u analysis/elai_amy_long.py --mode fit --k 5 --seed 20261012 --steps 3 --wall 900 --export marginal
.venv/bin/python -u analysis/elai_amy_long.py --mode gates
```

[observation] The actual binary command ended `-C 1 -c 5 -w 0 -s 3 -mg 10 -R 20261012 -o panel_k5_seed20261012 --pkk-site 1:104121000 --pk1`, using the persistent full `work_long/panel.inp` and `panel.pos`. It **started and completed cleanly**, return code=1 with finished sentinel, fit wall **457.873 s** (driver step wall 458.163 s), no wall-stop reason. The marginal FIFO consumed **1,000 individuals, 54,370 markers, 271,850,000 values**. Composition files cover all 1,000 individuals. The site file has **25,000 ordered rows**, ordered-sum and triangular-q-sum errors both **8.881784197001252e-16**, q symmetry=0, finite probabilities in [0,1]. Rounded marginal normalization error=0.001000000000000223; this is a copy-count rounding diagnostic, not a failed joint normalization test. Full ordered coverage remains NOT TESTED, as approved.

[observation] Full-density convergence trace:

| Step | Log likelihood | Reported increment | Relative increment | Elapsed s |
|---:|---:|---:|---:|---:|
| 1 | -27536139.009 | -27536139.009 | first step | 77.055 |
| 2 | -17169750.558 | 10366388.451 | 0.37646485034128485 | 146.609 |
| 3 | -10611771.828 | 6557978.731 | 0.38194956343989533 | 216.664 |

[observation] Binary stdout contains all three updates and `-elai: finished`. JSON records completed_steps=3, cap=3, convergence=FAIL, first qualifying step=None. The default dense gate command executes successfully with verdict **PARTIAL**, because canonical analysis fits have not yet been launched; the verification seed is deliberately not treated as the preregistered primary. The fixture gate's FAIL above is the completed tiny scientific-gate test. Neither result authorizes biological interpretation.

[inference] Updated timing calculation: slope `(216.6643283 - 77.0551958)/2 = 69.8045663 s/step`; estimated initialization `77.0551958 - 69.8045663 = 7.2506295 s`; observed post-EM posterior/output overhead `457.8729141 - 216.6643283 = 241.2085857 s`. Therefore `7.2506295 + 200×69.8045663 + 241.2085857 = 14,209.3724681 s = 3.947 h`. This is an extrapolation, not a convergence promise. Marginal export has finite, material overhead; it is not free. During a full-density export snapshot at ~379 s, `ps` showed the binary at 95.7% CPU and the Python driver at 8.5% CPU. [inference] Together with complete draining of the full-scale stream, this supports that the marginal FIFO consumer is not the dominant bottleneck on this host; exact posterior passes and text generation still need budgeting. Full-density *parallel* timing was not tested; only fixture parallel I/O isolation was tested.

[observation] Independent audit (`long_independent_validation.json`) checks SHA256 of all seven retained inputs before/after the full fit, equality of `genotypes.txt` to the header-stripped `panel.inp` payload, immutable raw pgen/pvar/psam sizes/mtimes unchanged, and the three canonical prefixes unused. The retained build before fit outputs occupies **338,842,030 bytes (~323.145 MiB)**. After verification, `work_long/` including outputs occupies **680,255,214 bytes (~648.742 MiB)**; all regular files under the dense output directory total **686,813,148 bytes (~0.687 GB)** at the audit snapshot, including the earlier lane's artifacts. Retained per-fit output is **341,413,184 bytes**, chiefly the automatic legacy upper-state `ps21` table (326,221,000 bytes). A conservative forecast retaining three additional K5-sized footprints is **1,711,052,700 bytes (~1.711 GB)**, below the **2 GB** bound even including verification and earlier outputs. K3 should be smaller; this is a resource forecast, not a hard runtime guarantee. No long-run input was deleted.

[observation] `df -h /` before building (`long_df_before.txt`):

```text
Filesystem      Size  Used Avail Use% Mounted on
/dev/nvme0n1p2   14T   14T  166G  99% /
```

[observation] `df -h /` after verification (`long_df_after.txt`):

```text
Filesystem      Size  Used Avail Use% Mounted on
/dev/nvme0n1p2   14T   14T  165G  99% /
```

[observation] Provenance manifests live under `analysis/runs/*elai-amy-long*/`; slugs include the work directory, mode/prefix and a nanosecond suffix to avoid parallel manifest collisions. They record source/preregistration/binary hashes, parameters, exact commands and step-summary hashes. Large immutable raw inputs are recorded by size/mtime, not claimed newly rehashed. No files were staged by this worker; the global index was not inspected because all git commands are prohibited. Gate/doc artifacts may be regenerated later: manifests are historical snapshots, not assertions that future bytes remain unchanged.

[inference] Recommended next step: independent parent review, then launch the three documented background fits if the per-fit wall margin is acceptable; after completion run gates and respect FAIL/PARTIAL. This deliverable validates the driver and export conventions, not ancient-population trajectories or AMY biology.
