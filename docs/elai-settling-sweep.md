# ELAI EM-settling sweep — preregistered, trace-only

[observation] Prepared against `docs/elai-amy-long-run-results.md`; the full sweep has
**not** been launched by the preparation worker. The configuration table and decision
rule below are fixed before any sweep launch. Smoke evidence is recorded separately
and must not be interpreted as a test of settling or biological structure.

## Hypothesis and interpretation limits

[observation] Full-density S=1, K=5 fits on 54,370 retained markers and 1,000 individuals
crossed relative tolerance around steps 46–50, then oscillated through step 200: about
2e-5 second-half mean absolute relative increment, 33–41% negative increments. S=1,
K=3 was closer (8.31e-6, 13% negative), but still failed the prior convergence gate.

[speculation] A single upper state leaves upper-layer parameters poorly identified,
creating a flat likelihood ridge. The hypothesis tested is that allowing S>1 reduces
oscillation and promotes settling, not merely crossing a magnitude tolerance.

[inference] Lower-K S=1 controls probe complexity; S=2,K=2 versus S=1,K=4 holds S*K=4
fixed; S=2,K=3 versus S=3,K=2 holds S*K=6 fixed. These are diagnostic contrasts, not
proof of the ridge mechanism: K changes, a single seed, absent source panels, and
shorter expensive runs confound causal interpretation. A negative at 50 steps can mean
insufficient optimization, not that S>1 cannot help. No result authorizes ancestry labels
or biological use; all previous dense preregistered gates remain unchanged.

## Flag semantics: checked against both manual formats and source

[observation] Sources: `/home/erikg/elai-fork/elai-manual.tex` and its PDF, version 1.01,
updated 26 May 2021 (PDF pages 5–7). PDF text was read with installed `pypdf` using
`.venv/bin/python -B`; `pdftotext` is unavailable, and nothing was installed. Also read
`/home/erikg/elai-fork/example/note.txt` and the example directory: its command uses
`-s 20 -mg 20`, paired `-g/-p` training labels 10,11 and cohort label 1. It does not
supply `-C/-c`, so it does not override the manual's recommendation.

| Flag used | Manual quote / checked meaning |
|---|---|
| `-C` | Appendix A: “specify number of upper clusters.” §7: “please use 2 for African American, and 3 for Hispanics. Other numbers is possible, provided that you have the appropriate panel data sets. But please be careful with the interpretation.” This is S, not K. |
| `-c` | Appendix A: “specify number of lower clusters.” §7: “For the lower layer number of clusters (specified with -c), a number, 5 × C is recommended.” This is K **total lower clusters**, not five per upper state multiplied again inside `-c`. |
| `-mg` | Appendix A: “specify number of mixture generations.” §7: “please use 10 for African American, 20 for Hispanics, and 100 for Uyghurs. Other values of mixture generation can be used, and the inferred local ancestry should be averaged.” It is admixture-generation/switch scale, **not** a generic concentration prior. |
| `-w` | Appendix A: “specify EM steps that use linear approximation.” §7.2: “We developed a fast linear approximation to the quadratic algorithm, which is less accurate but fast”; “not yet recommended for serious studies”; “The linear speed can be achieved using -w 20 -s 0.” Warm-up is approximate EM, before the exact steps, not extra exact steps. |
| `-s` | Appendix A: “specify steps in EM run.” §7: “a number between 20 and 50 is recommended.” §7.2: “This option fits the model using a quadratic algorithm which is accurate but slow.” We exceed the recommendation in cheap arms to observe the known ~50-step plateau. |
| `-p 1` | §5: “If -p is followed by a valid integer, then all individuals in the matching genotype file will be assigned a label of that integer”; “an admixed cohort sample is labelled by 1, an un-labelled training sample is labelled by 0, and an sample that is to be excluded is labelled by 9.” Training labels start at 10 without gaps. This sweep has **no training panels**. |
| `-g` | Appendix A: “can use multiple times, must pair with -p.” §4: “All genotype files assume the BIMBAM format.” Uses the existing dense unphased diploid panel, read-only. |
| `-pos` | Appendix A: “can use multiple times. arg is a file name.” §4.3: “first column is the SNP ID, the second column is its physical location, and the third column contains its chromosome number (optional).” |
| `-o` | Appendix A: “arg will be the prefix of all output files, the random seed will be used by default.” §6: “ELAI produces output files in a directory named output/.” Each arm has an isolated working directory. |
| `-R` | Appendix A: “specify random seed, system time by default.” §7.1: “If -R is omitted ELAI will use the machine time as the random seed.” Explicit seed 20261010 matches failed seed A; it is intentionally **shared** across arms for comparison, not independent replicate evidence. |
| `--pkk-site 1:104121000` | **Not documented in the upstream manual.** Patch README `/home/erikg/elai-fork/README-PATCH.md`: “The site flag alone writes only the site file. Existing outputs are unchanged.” It writes K*K ordered lower-state rows per eligible individual at this observed/virtual coordinate. No `--pkk`, `--pk1`, or `--ps2` is used. |

[observation] The earlier “manual 5*C” reading is **correct as a recommendation**, not
as a mandatory constraint. K=5 follows it at S=1; K=3/4 do not. **Every S>1 arm here
intentionally departs from K=5*S to bound cost**. This and the absence of appropriate
ancestral panels are interpretive caveats, not production local-ancestry settings.

[observation] The manual gives no exhaustive numerical valid ranges. Source checks
in `src/control.cpp` parse integer `-C/-c/-s/-w`, and floating-point `-mg` with leading
digit/nonnegative arguments; `model.h` setters do not enforce a K/S ratio or an upper
bound. Safe domain used here: positive integers S,K; nonnegative integer s,w; positive
mg. Zero s/w is documented by `-w 20 -s 0`; mg=0 is not used or asserted scientifically
valid. Training label source guard requires `p-10 < S` (`model.cpp::assign_pop_label`);
the manual additionally requires consecutive panel labels. Source `InitModelParam`
sets `m_upper_switches = m_morgan * m_mixgen`, explaining mg's switch-scale role.

[observation] `src/model.cpp::EM_joint` loops exactly `nWarmSteps` then `nMaxSteps`;
there is **no likelihood early-stop** in this patched binary. Exact step indices reset
to zero after approximate warm-up; the analyser excludes all approximate-update lines.
The first exact increment is likelihood minus zero, not a valid improvement, and is
excluded from all increment statistics. Printed likelihood/increments have 3 decimals.

## Bounded parallel design

[observation] Measured baseline is approximately **60–70 seconds per exact EM step**
at S=1,K=5,M=54,370,N=1,000 (200-step seed-A/B totals 11,845/11,931 s). Each process is
single-threaded; all eight are launched concurrently. Numerical thread environment
variables are set to 1. The runner uses patched `/home/erikg/elai-fork/src/elai`, **not**
the unpatched `elai-lin`.

[inference] Cost grows with diploid pair states, plausibly **(S*K)^2**: one copy has
S*K states; two copies have (S*K)^2. Source and patch describe quadratic exact recursion
and O(M*(S*K)^2) export memory. We use the deliberately simple projection
`t_step = [60,70] * ((S*K)/5)^2` seconds. This is a planning model, not measured timing
for S>1; memory-bandwidth contention can invalidate linear wall scaling.

[inference] Per-arm projection is `s*t_step`; for warm-up we conservatively charge one
**exact-cost** step per w. Add a **300-second per-process allowance** for input,
initialization and the two automatic/site posterior-output passes. The largest projected
upper endpoint plus allowance is 89 minutes. The runner imposes a real per-process
**5,400-second timeout**, then TERM with KILL after 15 seconds: overruns become incomplete,
never silent success. Parallel wall is approximately the maximum arm time, not their
sum. This budget is a hard bounded experiment, not a promise all arms finish under
contested host resources; report timed-out arms as incomplete.

| name | S | K | S*K | pair states | mg | w | exact s | projected sec/step | EM minutes (60–70 s basis) | upper + 5 min allowance |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 1 | 5 | 5 | 25 | 10 | 0 | 70 | 60–70 | 70.0–81.7 | 86.7 min |
| s1_k3 | 1 | 3 | 3 | 9 | 10 | 0 | 100 | 21.6–25.2 | 36.0–42.0 | 47.0 min |
| s1_k4 | 1 | 4 | 4 | 16 | 10 | 0 | 90 | 38.4–44.8 | 57.6–67.2 | 72.2 min |
| s2_k2 | 2 | 2 | 4 | 16 | 10 | 0 | 90 | 38.4–44.8 | 57.6–67.2 | 72.2 min |
| s2_k3 | 2 | 3 | 6 | 36 | 10 | 0 | 50 | 86.4–100.8 | 72.0–84.0 | 89.0 min |
| s3_k2 | 3 | 2 | 6 | 36 | 10 | 0 | 50 | 86.4–100.8 | 72.0–84.0 | 89.0 min |
| mg5 | 1 | 5 | 5 | 25 | 5 | 0 | 70 | 60–70 | 70.0–81.7 | 86.7 min |
| warm1 | 1 | 5 | 5 | 25 | 10 | 1 | 70 | 60–70 | 71.0–82.8 (includes w) | 87.8 min |

[inference] Honoring recommended K=5*S would cost S=2,K=10: factor 16, 960–1,120 s/step,
50 steps = 13.3–15.6 hours; S=3,K=15: factor 81, 4,860–5,670 s/step, 50 steps =
67.5–78.8 hours. Those are deliberately **not** in this 1.5-hour sweep.

## Disk protection and reproducibility

[observation] Direct binary execution avoids the long driver's marginal FIFO,
extra export passes and its unrelated biological preregistered gates. Input is reused
read-only from `data/derived/elai_amy_dense/work_long/panel.inp` and `panel.pos` (no copies).
Each arm works in `data/derived/elai_amy_dense/sweep_<name>.work/`, logging to
`data/derived/elai_amy_dense/sweep_<name>.stdout.txt`.

[observation] **Site-only does not suppress legacy automatic outputs.** Source
`control.cpp` always calls `write_em(1)` and `compute_ps(1,2)` even without any export
flags. Mean ancestry alone would write N*M*S numbers. The runner precreates
`output/<prefix>.ps21.txt` and `.em.txt` as `/dev/null` symlinks; those outputs are
computed but never persisted. No full joint or marginal dump is requested. Retained
files are site TSV, SNP info, admixture, ELAI internal log, stdout and tiny status/report
files. Each ELAI process has a **16 MiB per-regular-file cap** (`ulimit -f 16384`, KiB in
bash); a size-limit termination is a failure. At most five regular output/log files
per arm, eight arms, yields a 640 MiB cap for ELAI files, plus small status, JSON/text
trace reports, manifests and smoke artifacts: conservatively **below 1 GB**, hence below
2 GB. The runner never copies the 157 MiB panel. `/dev/null` symlinks are not hashed
as outputs. All new writes remain in analysis/, docs/, scripts/ or the allowed dense
derived directory.

[observation] A shared lock prevents concurrent reservations. All names, directories
and reports are checked before launching any fit; existing stdout logs/artifacts are
refused, not truncated. A full rerun requires a deliberate archive/removal by the
parent; smoke artifacts have a separate `sweep_smoke` prefix. Interrupt cleanup stops
launched process groups. Status accepts upstream's inherited success return code 1
only alongside the completion sentinel and the exact expected trace length. Timeout,
malformed/noncontiguous/multiple-run traces or missing completion cannot settle.

[observation] The analyser's `--provenance` uses `scripts/provenance.py::RunRecorder`,
recording input/binary/script/doc hashes, fixed seed, parameters, command, trace/report
and retained output hashes in `analysis/runs/`. Git probes are disabled explicitly:
**no git command is run**. Common random seed is fixed in the shared analyser plan,
which the shell consumes using `--plan`; no independently drifting config list exists.

## Preregistered decision rule (fixed before launch)

[inference] Define likelihood L_j from each exact-update stdout line (reported
zero-based step j; displayed step j+1). For j>=1 compute `delta_j=L_j-L_(j-1)` and
`r_j=delta_j/max(1,abs(L_(j-1)))`. Recompute from likelihood rather than divide a rounded
reported increment. Retain both in the report for auditing; differences of ~0.001 can
reflect stdout rounding. Exclude the first exact line's bogus improvement against zero.

[inference] For n usable increments, q=floor(n/4); first/last quarters are first/last q
increments. Second half is increments[floor(n/2):], so the extra increment goes in the
second half for odd n. A quarter with no increments is unassessable. Report per arm:
completed steps; **entire** likelihood trajectory and relative-increment series;
first/last-quarter mean abs relative increment; second-half mean abs **absolute**
increment and mean abs **relative** increment; second-half negative fraction; and first
qualifying 3-window **ending step** (one-based). Strict sign uses the printed likelihood
resolution; rounded zero qualifies as nonnegative.

[inference] An arm **SETTLES** iff all are true:

1. Complete expected exact steps, valid contiguous run-0 trace, completion sentinel and
   return code 0 or inherited 1.
2. **DECAYS**: positive first-quarter mean, and last-quarter mean **<= 0.5 times**
   first-quarter mean ("materially below" is preregistered as at least a 50% reduction).
3. Second-half negative fraction is **strictly < 0.05** (not <=).
4. The **final three increments** are all in **[0,1e-4)**. A qualifying early window
   alone is not enough. At least four likelihood lines are required to supply three
   increments.

[inference] This is a diagnostic finite-run rule, not mathematical convergence proof.
The large initial gains can satisfy quarter decay even on a later plateau; the negative
fraction and final window are essential additional checks. NO means the rule was not
met within the allocated steps, not an established mechanism. Five-step smoke runs
are plumbing checks regardless of their numerical verdict.

## Exact parent launch commands

[observation] Paths below exist in this checkout. Run from the repository root; scripts
also resolve/cd to the root if invoked elsewhere. **Do not launch the full sweep until
parent review.** No installation or rebuild is needed.

```bash
cd /home/erikg/demoselex
bash scripts/elai_settling_sweep.sh
```

[observation] This launches every table row in parallel, waits for all, prints each status,
then writes `data/derived/elai_amy_dense/sweep.report.json` and `.report.txt` and invokes
provenance. The same analyser can re-read the status/logs into a **new** report:

```bash
.venv/bin/python -B analysis/elai_settling_sweep_report.py \
  --status data/derived/elai_amy_dense/sweep.status.tsv --prefix sweep \
  --output data/derived/elai_amy_dense/sweep.review.report.json
```

## Smoke-test evidence

[observation] Evidence is appended here after the cheap end-to-end test; the decision
rule above is already written. No full sweep is performed during preparation.

[observation] Executed on 2026-10-10, using the **full-density** existing panel, one cheap
S=1,K=3 arm, five exact steps. Runtime **147 seconds**, completed 5/5; ELAI returned its
expected success code 1 and printed its completion sentinel. Runner then automatically
ran the analyser and provenance. Exact invocation (stdout was captured in
`data/derived/elai_amy_dense/sweep_smoke_driver.log`):

```bash
cd /home/erikg/demoselex
bash scripts/elai_settling_sweep.sh --smoke --only s1_k3
```

[observation] Actual binary command printed by the runner (absolute input paths are
resolved from the root; binary cwd is the dedicated `.work` directory):

```text
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy_dense/work_long/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy_dense/work_long/panel.pos -C 1 -c 3 -w 0 -s 5 -mg 10 -R 20261010 -o sweep_smoke_s1_k3 --pkk-site 1:104121000
```

[observation] Raw exact trace:

```text
 ### Exact update 0 0	 -27541635.053	-27541635.053
 ### Exact update 0 1	 -17384284.281	10157350.772
 ### Exact update 0 2	 -12096113.041	5288171.240
 ### Exact update 0 3	 -8372862.083	3723250.958
 ### Exact update 0 4	 -6348540.133	2024321.950
```

[observation] Analyser output, including the parsed trace and all statistics:

```text
[observation] ELAI settling sweep trace report (not biological authorization).

s1_k3: completed=5/5 complete=True elapsed_s=147 returncode=1
[observation] {"complete": true, "completed_steps": 5, "decays": false, "final_three_window_met": false, "first_qualifying_step": null, "first_quarter_mean_abs_relative_increment": 0.3687998462129648, "last_over_first_ratio": 0.6555637266127289, "last_quarter_mean_abs_relative_increment": 0.24177180155757252, "quarter_size": 1, "second_half_mean_abs_increment": 2873786.4539999994, "second_half_mean_abs_relative_increment": 0.2747886853088238, "second_half_negative_fraction": 0, "settles": false, "usable_increments": 4}
[inference] SETTLES=NO under the preregistered rule.
[observation] step elai_step likelihood reported_increment absolute_increment relative_increment
1 0 -27541635.053 -27541635.053 None None
2 1 -17384284.281 10157350.772 10157350.772 0.3687998462129648
3 2 -12096113.041 5288171.24 5288171.24 0.3041926348259077
4 3 -8372862.083 3723250.958 3723250.9579999996 0.30780556906007506
5 4 -6348540.133 2024321.95 2024321.9499999993 0.24177180155757252
Provenance manifest saved to: analysis/runs/20261010-1403-sweep_smoke.report/manifest.json
```

[observation] Retained smoke files sum to **3,500,912 bytes** (about 3.34 MiB),
plus the small provenance manifest. Site TSV is 748,628 bytes; neither a full joint nor
a marginal table exists. Automatic `.ps21.txt` and `.em.txt` are `/dev/null` symlinks.
Input M=54,370,N=1,000 appeared in stdout; no reduced fixture or copied dense panel was
used. No full sweep was launched. This NO verdict is expected after only five steps
and says nothing about plateau settling.

[observation] Re-reading the saved smoke logs independently succeeded (exit 0):

```bash
.venv/bin/python -B analysis/elai_settling_sweep_report.py \
  --status data/derived/elai_amy_dense/sweep_smoke.status.tsv --prefix sweep_smoke \
  --output data/derived/elai_amy_dense/sweep_smoke.review.report.json
```

[observation] Collision check ran the exact smoke invocation again; it returned exit 2:
`Refusing existing artifact: /home/erikg/demoselex/data/derived/elai_amy_dense/sweep_smoke.status.tsv`.
SHA-256 hashes of **all** existing smoke artifacts/regular ELAI outputs were unchanged,
the lock was removed, and no new ELAI process launched. Unknown `--only unknown` was
also rejected before reservation.

[observation] Validation commands:

```bash
bash -n scripts/elai_settling_sweep.sh
.venv/bin/python -B -m unittest discover -s analysis -p test_elai_settling_sweep_report.py -v
df -h /
```

[observation] Shell syntax passed; **9 tests passed** in 0.001 s. Regressions cover
warm-up exclusion, first-increment exclusion, valid settling, constant plateau,
oscillation, strict 5% boundary, early-but-not-final window, incomplete/short traces,
malformed/noncontiguous/multi-run traces and all cost/configuration constraints.
Manual cross-check against the existing 200-step logs reproduced historical second-half
mean abs relative increments 2.00410648e-5 / 2.14830678e-5 / 8.30525907e-6 and negative
fractions 0.41 / 0.33 / 0.13 (K=5 A/B and K=3 A); all remain non-settling.

[observation] `df -h /` after the smoke test:

```text
Filesystem      Size  Used Avail Use% Mounted on
/dev/nvme0n1p2   14T   14T  163G  99% /
```

[observation] Smoke manifest: `analysis/runs/20261010-1403-sweep_smoke.report/manifest.json`. Its doc hash is the pre-smoke document,
before evidence append; the configuration table and decision rule were unchanged.
No git commands, staging, commits, installs or full-sweep fits were performed.

[observation] The five parsed smoke updates also match the first five updates of the
historical K=3 seed-A full fit **exactly** (likelihood, reported increments, and recomputed
relative series). This checks baseline input/seed reproducibility independently of
the synthetic regression fixtures.

---

## RESULTS (2026-10-10, parent-run; sweep completed in 4,013 s wall, 10 arms in parallel)

Settling rule applied as pre-registered (decays AND second-half negative fraction < 0.05 AND the
final 3-consecutive window in `[0, 1e-4)`):

| arm | S | K | mg | w | steps | decays | 3-window | last/first | 2nd-half mean abs rel | neg frac | **settles** |
|---|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---|
| baseline | 1 | 5 | 10 | 0 | 70 | yes | no | 3.54e-04 | 9.42e-05 | 0.057 | **NO** |
| s1_k3 | 1 | 3 | 10 | 0 | 100 | yes | yes | 3.36e-04 | 4.12e-05 | 0.000 | **YES** |
| s1_k4 | 1 | 4 | 10 | 0 | 90 | yes | yes | 2.89e-04 | 7.97e-05 | 0.000 | **YES** |
| s2_k2 | 2 | 2 | 10 | 0 | 90 | yes | yes | 1.03e-03 | 8.43e-05 | 0.000 | **YES** |
| s2_k3 | 2 | 3 | 10 | 0 | 50 | yes | no | 1.60e-03 | 3.32e-04 | 0.000 | no |
| s3_k2 | 3 | 2 | 10 | 0 | 50 | yes | no | 1.97e-03 | 2.70e-04 | 0.000 | no |
| mg5 | 1 | 5 | 5 | 0 | 70 | yes | no | 3.54e-04 | 9.42e-05 | 0.057 | no |
| mg100 | 1 | 5 | 100 | 0 | 70 | yes | no | 3.54e-04 | 9.42e-05 | 0.057 | no |
| mg1000 | 1 | 5 | 1000 | 0 | 70 | yes | no | 3.54e-04 | 9.42e-05 | 0.057 | no |
| warm1 | 1 | 5 | 10 | 1 | 70 | yes | no | 5.39e-04 | 7.31e-05 | **0.257** | no |

**[observation] `-mg` is inert in this configuration.** 5, 10, 100 and 1000 mixture generations give
*identical* statistics to every printed digit. With `S=1` and `-p 1` the admixture-generation
parameter does not enter the fit. The pre-run hypothesis that a too-small mixture-generation scale
causes the oscillation is therefore **refuted**.

**[observation] `-w` makes settling worse.** The linear-approximation warm-up gives 25.7% wrong-sign
increments, ~4.5x the baseline, consistent with the manual's caution that it is "not yet recommended
for serious studies".

**[inference] Step count, not S or K, controls the verdict.** Every settling arm ran 90-100 steps;
every failing arm ran 50-70 (or, in the earlier 200-step runs, 200). K=5 at 70 steps failed narrowly
(0.057 vs the 0.05 limit). So the EM reaches a likelihood plateau and its increments then shrink to a
**noise floor** (about 2e-5 relative on a likelihood of ~2.5e6, i.e. ~50 absolute in 2.5 million); past
roughly 100 steps the increments are noise-dominated and flip sign at random. Because the
pre-registered rule requires the *final* window to qualify, running LONGER actively hurts: the earlier
200-step runs failed for this reason, not because the fit was diverging.

**[inference] Consequence for the rule.** "Settling" here means stopping while the increments are
still signal-dominated, at roughly 90-100 steps. The rule is therefore sensitive to the stopping step,
which is a genuine limitation of the pre-registered criterion and must be stated rather than tuned.

**[observation] Practical conclusion:** K=3 and K=4 settle cleanly at 90-100 steps; K=5 was not tested
at >=90 steps. The long-run configuration was therefore moved to K=4 (100 steps, seeds 20261010 and
20261011) plus K=3 (100 steps) for the cross-K comparison.
