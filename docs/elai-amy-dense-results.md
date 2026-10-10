# Dense ELAI AMY revision — STOP

[observation] **Gated STOP.** This is an instrument/resource result, not an AMY biological result. No allele-frequency analysis, local-state biological summary, effective independent-state estimate or biological claim was made. No cluster labels from incomplete/non-converged EM are accepted.

## Revised inputs and failure-mode rationale

[observation] Immutable release: 52,382,872 variants × 16,339 samples. Eligible unique dated AADR-matched ancients: 5,981; fitted selection **1,000**, **1,000 carrying dates** (seed 20261010). Baseline had128; no N128 fallback. Dates span 150.0–13433.0 BP. Metadata labels never enter QC or EM.

[observation] Chr1 inclusive85,000,000–120,000,000, core103,915,000–104,327,000 excluded. **543,791 full-density biallelic SNP candidates**, **70,161 passing MAF≥0.05/call rate≥0.95 BEFORE thinning**, **54,370 retained after100-bp greedy spacing**. Baseline's thin-before-QC gave1,164. Input size 163,966,772 bytes, deleted with scratch. Hardcalls are native ALT calls (not rounded dosages); complete imputed calls are not original aDNA coverage.

[observation] Independent PC panel: **9,854 markers**, autosomes2–22 (chr1 excluded), ≥20-kb candidate spacing (121,871 candidates), cohort QC (15,258 passing), then≥100-kb retained spacing. Five eigen-PCs of standardized mean-imputed hardcall sample Gram matrix; all eigenvectors deterministic with single-thread BLAS. Baseline had145 chr2 SNPs. Per-chromosome counts: `{"10": 539, "11": 533, "12": 541, "13": 408, "14": 341, "15": 305, "16": 288, "17": 316, "18": 308, "19": 216, "2": 933, "20": 240, "21": 149, "22": 131, "3": 774, "4": 786, "5": 710, "6": 705, "7": 627, "8": 562, "9": 442}`. Builder wall 97.256s.

[inference] These prospective changes address sparse fitted LD coverage, the N128 pilot, non-convergence and inadequate independent structure validation. They were registered before fitting and motivated by diagnosed failure modes, not by revised gate values. No thresholds were relaxed.

## EM, live cost calibration and wall clocks

[observation] Fixed S=1,K=5,s=200 (exact quadratic EM), no early stop/checkpoint in this binary. Supervisor-approved live calibration uses the first five completed steps of the primary itself, not an extra fit/restart. Terminal convergence requires a completed fit and three consecutive nonnegative relative likelihood increments<1e-4. Full trace is below; labels from a partial trace are invalid.

[observation] **panel_k5_seed20261010**: wall **303.493s**, exit-15, completed fit=False, completed steps=5/200, convergence **FAIL**, final relative increment=0.355390479534299; first qualifying step=None. Reason: RESOURCE STOP: live five-step projection exceeds primary wall cap.

```sh
/home/erikg/elai-fork/src/elai -g /tmp/elai-amy-dense-46qwl_7c/panel.inp -p 1 -pos /tmp/elai-amy-dense-46qwl_7c/panel.pos -C 1 -c 5 -w 0 -s 200 -mg 10 -R 20261010 -o panel_k5_seed20261010 --pkk-site 1:104121000 --pkk
```

| Step | Log likelihood | Reported increment | Relative increment | Elapsed s |
|---:|---:|---:|---:|---:|
| 1 | -27539198.975 | -27539198.975 | N/A (first step) | 79.777 |
| 2 | -17188214.613 | 10350984.362 | 0.37586367 | 133.820 |
| 3 | -10532003.244 | 6656211.368 | 0.38725438 | 189.363 |
| 4 | -5868834.609 | 4663168.635 | 0.44276179 | 244.912 |
| 5 | -3783106.663 | 2085727.946 | 0.35539048 | 303.477 |

[observation] Live calibration: 55.925s/step (steps1–5 slope), initialization estimate 23.852s; N×M×K²=1,359,250,000; 4.1144031e-08s/unit. Projected200-step EM time **11208.857s**, cap1200s; estimated affordable steps **21** (export costs omitted, so optimistic). No configuration change/restart was made after this cost measurement.

[inference] At this measured scaling, M=54,370,N=1,000,K=5,s=200 does not fit the cap. Keeping M and200 steps would permit approximately N≤105; keeping N and200 steps would permit M≤5,717 (before export cost). Those are arithmetic resource projections, NOT evidence those settings converge or validate. User's marker/cohort targets were retained rather than silently reverting to the failed sparse setup.

| M | N | K | Projected seconds/step | Steps affordable in1200s (EM only) |
|---:|---:|---:|---:|---:|
| 5000 | 1000 | 5 | 5.143 | 233 |
| 10000 | 500 | 5 | 5.143 | 233 |
| 25000 | 200 | 5 | 5.143 | 233 |
| 50000 | 100 | 5 | 5.143 | 233 |
| 50000 | 1000 | 5 | 51.430 | 23 |
| 50000 | 1000 | 3 | 18.515 | 64 |

[observation] Parent corrected baseline arithmetic:128×1164×25=3,724,800 units,0.73/3,724,800=1.9598368e-07s/unit; M50,000,N1000,K5 is1.25e9 units (~244.98s/step at that baseline). Measured live slope is 0.210× that baseline, whose14.6s wall included export. Comparing EM slope to export-inclusive wall is only a cost heuristic.

[observation] Analysis fits run: **1/3 maximum**. SeedB and cross-K comparator NOT TESTED if absent below; no extra calibration fit, no frozen replay unless recorded. Primary-incomplete/non-converged STOP prevented spending on invalid labels.

## Numeric gates

| Gate | Numeric observations | Verdict |
|---|---|---|
| (a) Known structure | Independent PC markers=9,854; groups=513, singleton groups=349; PC silhouette=None, composition silhouette=None, group NMI=None, country NMI=None, singleton-robust country AMI=None, date-bin NMI=None, date/BP R²=None | **NOT TESTED** |
| (b) Stability | Seed=NOT TESTED, cross-K=NOT TESTED, virtual argmax=NOT TESTED, thinning refit=NOT TESTED; thinning posterior-subset diagnostic=None | **PARTIAL** |
| (c) Coverage | Usable=1000/1000 (100.0%); dated selected=1000; dated COMPLETE --pkk export=0; rows=0/1,359,250,000 | **FAIL** |
| (d) Virtual | Ordered sum error=None, q symmetry=None, triangle error=None; integration controls=0, linked/unlinked change=None | **NOT TESTED** |
| (e) Exclusion | Builder core markers=0; actual ELAI-retained core markers=0; PC chr1 overlap=0; virtual104121000 inside core=True; flanks=103785095/104328492 | Builder **PASS**; full gate **PASS** |
| (additional) EM convergence | Completed steps=5/200; final relative increment=0.355390479534299; complete fit=False; tolerance1e-4 over three terminal increments | **FAIL** |

[observation] The actual retained snpinfo was copied from the primary scratch output during EM, before cleanup, using `cp /tmp/elai-amy-dense-46qwl_7c/output/panel_k5_seed20261010.snpinfo.txt data/derived/elai_amy_dense/panel_k5_seed20261010.snpinfo.txt`. The report audit compares every retained ID/chromosome/coordinate with the selected panel. This is a pre-EM input/retention audit, not a posterior export.

[observation] `None` means unavailable, not zero or a measured failed silhouette. Known structure and its date-ordered/country-distance diagnostics cannot be computed without valid complete compositions; NOT TESTED does not imply absence of structure. Missing stability checks are NOT TESTED, overall PARTIAL, not falsely failed numerical comparisons. PARTIAL/NOT TESTED never pass the instrument gate. Coverage separates usable INPUT calls from complete EXPORT. ELAI writes snpinfo before EM: an independent exact ID/coordinate audit of that retained panel can PASS exclusion even with an incomplete EM/export; it does not authorize labels.

## Commands, provenance and resources

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -u analysis/elai_amy_dense.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=analysis .venv/bin/python -m unittest analysis/test_elai_amy_dense.py analysis/test_elai_amy.py -v
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python analysis/elai_amy_dense_report.py
```

[observation] Preregistration SHA256 `2fe2f2f576eb2c615babc3adddb9fd10c280320de06218707bd64b0bb7ea1c2c` verified unchanged, mtime before every fit. Independent artifact audit verifies counts, QC, spacing, chr1/core exclusion and scratch deletion. 14 helper regressions passed (six new dense regressions plus eight reused baseline-helper checks, `tests.log`). Required `scripts/provenance.py` RunRecorder used with git provider disabled before construction; small inputs, binary and outputs hashed, huge immutable raw inputs recorded by size/mtime rather than rehashing. No fork rebuild, install, raw edit or git command. Exact binary commands retain now-deleted scratch paths; deterministic script recreates inputs from immutable raw files.

[observation] Scratch removed=True; runner wall=400.846s. Persistent summary output=6,176,863 bytes; largest persistent artifact=2,815,458 bytes. Largest scratch input=163,966,772 bytes (<1GB); full tensor never persisted. New total output far below8GB. Final free bytes=177,412,804,608; `df -h /` before large input writes and after run:

```
Filesystem      Size  Used Avail Use% Mounted on
/dev/nvme0n1p2   14T   14T  166G  99% /
```

[inference] **STOP**. Before accepting labels, a resource budget that permits dense markers and converged EM (including full-export cost) must be secured or a new explicitly authorized/preregistered cohort/compute trade-off made. Then independent structure, complete export, linkage control, seed/cross-K and genuinely refitted thinning stability must be measured. This resource/convergence failure is not a biological negative result and does not establish ELAI fails generally. No downstream allele-frequency or AMY-history claim is authorized.
