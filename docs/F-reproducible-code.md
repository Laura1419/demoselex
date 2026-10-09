# F. Reproducible code

> Owner tasks: all. Status: in progress. Every code file currently in the
> repository under `scripts/`, `sim/`, `analysis/` and `lean/` is inventoried below
> with what it **consumes**, what it **produces**, the **exact command** to run it,
> **determinism/seed** notes, and **data prerequisites**. Files not yet written are
> listed as `TBD` with the task that owns them; nothing here is a claim that
> unbuilt code already runs. All commands are run from the repository root unless
> stated otherwise.

Conventions used below: `[deterministic]` = no random number generator (RNG)
touches the computation; `[seeded]` = an explicit seed is required or configured;
`[synthetic]` = needs no downloaded data; `[real-data]` = needs files acquired per
`docs/E-data-inventory.md`. Every analysis run is expected to record a manifest via
`scripts/provenance.py` (see `docs/PROVENANCE.md`); manifests live under
`analysis/runs/<timestamp>-<slug>/manifest.json`.

---

## 1. `scripts/` — infrastructure, provenance, orchestration

### `scripts/__init__.py`
- **Consumes:** nothing. **Produces:** nothing (package marker so `scripts.*` is
  importable). **Run:** imported, not executed. `[deterministic]` `[synthetic]`.

### `scripts/provenance.py`
- **Consumes:** the running process (git revision/branch, Python and package
  versions, platform), declared input files, parameters, an optional seed, and the
  command line. **Produces:** `analysis/runs/<timestamp>-<slug>/manifest.json` with
  git info, software versions, input/output SHA-256 hashes, parameters and seed.
- **Public interface:** class `RunRecorder(slug, command=None)` (context manager;
  `record_input`, `record_output`, `record_parameter`, `set_seed`, `save`) and
  `run_with_provenance(...)`. Also `verify_manifest(path)` and a `demo_run()`.
- **Run:**
  ```bash
  python3 scripts/provenance.py --help
  python3 scripts/provenance.py record --slug example --input docs/MISSION.md \
      --parameter design=demo --seed 0 --cmd "python3 example.py"
  python3 scripts/provenance.py verify analysis/runs/<ts>-<slug>/manifest.json
  python3 scripts/provenance.py demo
  ```
- **Determinism:** the harness itself is `[deterministic]`; it *records* whatever
  seed the caller supplies (Lean builds conventionally record seed `0` for
  "no randomness"). `[synthetic]` (only needs the repo). Note: `analysis/runs/` is
  gitignored (`analysis/.gitignore`), so manifests are local artifacts, not
  committed.

### `scripts/provenance_example.py`
- **Consumes:** nothing (a self-contained example). **Produces:** an example
  manifest and `analysis/runs/example/output.json` when its
  `example_analysis_function` is called. **Run:**
  `python3 scripts/provenance_example.py` (imports `provenance` from the same
  directory). `[seeded]` examples call `np.random.seed(seed)`. `[synthetic]`.

### `scripts/build_taskgraph.py`
- **Consumes:** the in-file declarative `TASKS` list. **Produces:** the WorksGood
  task graph in `.wg` (idempotent: existing ids are skipped) and publishes each
  connected component. **Run:**
  `python3 scripts/build_taskgraph.py` (add `--rebuild` to abandon and recreate).
  Requires the `wg` CLI. `[deterministic]` `[synthetic]`.
- **Note:** contains a hard-coded absolute `ROOT = "/home/erikg/demoselex"`; this
  must be updated or made relative before use on another host.

### `scripts/render_taskgraph.py`
- **Consumes:** `scripts.build_taskgraph.TASKS`. **Produces:** `docs/TASKGRAPH.md`
  (Mermaid + table) and `docs/DECISIONS.md`. **Run:**
  `python3 scripts/render_taskgraph.py`. `[deterministic]` `[synthetic]`.
  Overwrites those two docs, so it is an owner-append-only hazard elsewhere.

### `scripts/make_docs_skeleton.py`
- **Consumes:** the in-file `DOCS` dictionary of skeleton texts. **Produces:** the
  `docs/A-…J-*.md` skeletons (idempotent; existing files are kept unless
  `--force`). **Run:** `python3 scripts/make_docs_skeleton.py [--force]`.
  `[deterministic]` `[synthetic]`.

### `scripts/status.py`
- **Consumes:** `wg list --all --json`, `git log`, `wg service status`,
  `git branch --list wg/*`. **Produces:** a printed one-glance project status
  (no files). **Run:** `python3 scripts/status.py`. `[deterministic]` `[synthetic]`.

### `scripts/supervise.py`
- **Consumes:** `wg list --all --json`. **Produces:** `wg retry` side effects for
  failed tasks below a retry cap; prints retryable/exhausted/waiting counts.
  **Run:** `python3 scripts/supervise.py [--max-attempts N] [--dry-run]`.
  `[deterministic]` `[synthetic]`. Never mutates shared git state.

---

## 2. `sim/` — statistics, simulators, and validation

### `sim/trajectory_stat.py`
- **Consumes:** a pandas long table with columns
  `locus, trajectory, time, alt_count, total_count`; optional `family`
  (trajectory→family mapping) and `relatedness` (labeled family correlation
  DataFrame). **Produces:** `TrajectoryResult(slopes, summary)` — a tidy per-
  trajectory table (`locus, trajectory, slope, se, direction, n, span, family`) and
  a per-locus summary (pooled Model 0 slope, concordance, family concordance, GLS
  meta slope/SE, `z_score`, Wald `p_value`, `n_effective`, `q`, `q_df`,
  `q_p_value`, `i_squared`). Public entry point
  `trajectory_statistic(data, *, family=None, relatedness=None)`; equations in
  `docs/D-method.md` §3–§6.
- **Run (tests / validation):**
  ```bash
  python3 -m pytest -q sim/test_trajectory_stat.py
  python3 -m sim.validate_trajectory_stat --seed 31415     # writes a run manifest
  ```
- **Determinism:** `[deterministic]` in fitting (no RNG). `[synthetic]` — the
  validation script generates seeded synthetic counts; needs numpy, pandas, scipy.
- **Data prerequisites:** none for tests/validation. Real use needs the trajectory
  labels and time bins produced by Stage A (`analysis/trajectory_infer.py`, `TBD`).

### `sim/empirical_null.py`
- **Consumes:** the same long count table plus a per-locus `metadata` table
  (`locus, chromosome, position, allele_frequency, recombination_rate, local_ld,
  imputation_quality, missingness, temporal_coverage`) and a list of `candidates`.
  **Produces:** `EmpiricalNullResult(results, matches, statistics)` — candidate
  scores/percentiles/empirical P, an auditable match table, and all fitted
  statistics. Entry point
  `empirical_null(counts, metadata, candidates, *, family=None, relatedness=None,
  calipers=None, exclusion_bp=1_000_000, min_controls=20, max_controls=None)`.
- **Run:**
  ```bash
  python3 -m pytest -q sim/test_empirical_null.py sim/test_trajectory_stat.py
  python3 -m sim.validate_empirical_null --seed 251008      # writes a run manifest
  ```
- **Determinism:** `[deterministic]` matching (no RNG); the validation script is
  `[seeded]` (`--seed` required). `[synthetic]` for tests/validation.
- **Data prerequisites:** matches must be QC-derived without the tested outcome;
  a sparse real dataset may be ineligible under exact time-support matching
  (see `docs/E-data-inventory.md` for the QC fields).

### `sim/akbari_glmm.py`
- **Consumes:** a per-locus observation frame (`time_bp`, `alt_count`,
  `total_count`) plus either background genotypes (with per-column
  `marker_chromosomes` and a `target_chromosome` for LOCO) or a precomputed
  sample-by-sample GRM `K`. **Produces:** per-locus fit dict `beta_time`,
  `beta_intercept`, `tau1`, `tau2`, `se_beta_time`, `z_score`, `p_value`,
  `converged`, `n_iterations`; class helpers `compute_grm(...)` and
  `leave_one_chromosome_out(...)`; `simulate_data(...)` returns a seeded
  one-locus fixture. This is a transparent PQL/REML approximation, not the
  production PQLseqPy scan (see §9 of `docs/D-method.md`).
- **Run:**
  ```bash
  python3 -m pytest sim/test_akbari_glmm.py -q
  PYTHONPATH=sim python3 - <<'PY'
  from akbari_glmm import AkbariGLMM, simulate_data
  counts, background = simulate_data(random_state=42, selection_strength=1.2)
  print(AkbariGLMM().fit(counts, genotypes=background).get("z_score"))
  PY
  ```
- **Determinism:** fitting is `[deterministic]`; `simulate_data` requires an
  explicit `random_state` (`[seeded]`). `[synthetic]`. Needs numpy, pandas, scipy.

### `sim/quick_validation.py`
- **Consumes:** `akbari_glmm` (imported top-level, so `sim/` must be on the path).
  **Produces:** printed calibration checks only (no files). **Run:**
  `PYTHONPATH=sim python3 sim/quick_validation.py`. `[seeded]` (loops over
  `random_state=i`). `[synthetic]`. Diagnostic only; the committed regressions are
  in `sim/test_akbari_glmm.py`.

### `sim/validate_trajectory_stat.py`
- **Consumes:** `--seed`; internally generates a fixed designed trend recovery
  test (three trajectories, nine times, depth 2000, slopes ±0.45). **Produces:** a
  provenance run directory with seeded counts/slopes/summary. **Run:**
  `python3 -m sim.validate_trajectory_stat --seed 31415`. `[seeded]` `[synthetic]`.
  Requires the repo root on `sys.path` (imports `scripts.provenance` and
  `sim.trajectory_stat`).

### `sim/validate_empirical_null.py`
- **Consumes:** `--seed`; internally simulates four neutral strata (two
  chromosomes), holds out 128 candidate loci and matches independent controls.
  **Produces:** a provenance run directory with count/metadata/rank/match tables and
  calibration metrics (observed at seed 251008; reported in
  `docs/G-validation-report.md`). **Run:**
  `python3 -m sim.validate_empirical_null --seed 251008`. `[seeded]` `[synthetic]`.
  This is a matcher stress test under a stipulated simulator, **not** a real-cohort
  calibration.

### `sim/test_trajectory_stat.py`
- **Consumes:** `trajectory_stat` (top-level import). **Produces:** pytest pass/fail.
  **Run:** `python3 -m pytest -q sim/test_trajectory_stat.py`. `[seeded]` fixture
  (seed 31415), `[deterministic]` assertions. `[synthetic]`. Covers both-sign trend
  recovery, tidy schema, analytic binomial SE, duplication invariance, relatedness
  discount, zero/tied signs, heterogeneity and a turnover-only pooled trend.

### `sim/test_empirical_null.py`
- **Consumes:** `empirical_null`, `trajectory_stat`. **Produces:** pytest pass/fail.
  **Run:** `python3 -m pytest -q sim/test_empirical_null.py`. `[deterministic]`
  fixture. `[synthetic]`. Checks chromosome/exclusion/support matching, calipers,
  and rank/percentile conventions.

### `sim/test_akbari_glmm.py`
- **Consumes:** `akbari_glmm`. **Produces:** pytest pass/fail. **Run:**
  `python3 -m pytest sim/test_akbari_glmm.py -q`. `[seeded]` (fixed `random_state`s),
  `[deterministic]` assertions. `[synthetic]`. Covers GRM centering/LOCO, a 30-seed
  structured null, a 10-seed strong single-locus trend, the GLM Fisher-information
  limit, and a reference comparison against released PQLseqPy for a fixed fixture.

### `sim/scripts/simulator.py`
- **Consumes:** a scenario config dict (loaded from YAML by `load_config`). The
  forward-time diploid Wright–Fisher mosaic simulator; replacements and migration
  draw whole diploid donors from the *previous* generation. **Produces:** via
  `save_results`, an aligned genotype matrix, sample metadata, a per-locus×trajectory
  `truth_table`, and `summary.json` (the CLI wrapper adds a provenance manifest).
- **Run:** exercised through the CLI below; its functions are also imported by
  `run_simulations.py`.
- **Determinism:** `[seeded]` — `config["seed"]` is a required nonnegative int, and
  all randomness uses one local `numpy.random.Generator` (never Python `hash()`).
  `[synthetic]`. Needs numpy, pyyaml.
- **Data prerequisites:** none (simulation); this is a validation scaffold, not a
  calibrated human-demography model (no mutation, ascertainment, damage, or
  realistic recombination map).

### `sim/scripts/run_simulations.py`
- **Consumes:** one scenario YAML or a directory of them, and an output parent dir.
  **Produces:** per-scenario output directories (see `simulator.py`) plus a
  `provenance.json` and a canonical manifest under `analysis/runs/`. **Run:**
  ```bash
  python3 sim/scripts/run_simulations.py sim/scenarios sim/output/demo
  python3 sim/scripts/run_simulations.py sim/scenarios/neutral_drift.yaml sim/output/demo
  ```
  `[seeded]` (uses each config's seed; records it via `RunRecorder`). `[synthetic]`.
  Requires the repo root on `sys.path` (inserted by the script) and the `wg`-less
  provenance harness only.

### `sim/scripts/test_simulator.py`
- **Consumes:** the six `sim/scenarios/*.yaml` and the `run_simulations.py` CLI.
  **Produces:** unittest pass/fail (runs each scenario twice and compares outputs for
  determinism). **Run:**
  `python3 -m unittest discover -s sim/scripts -p 'test_simulator.py' -v`.
  `[seeded]` `[synthetic]`.

### `sim/scenarios/*.yaml` (six configs)
- `neutral_drift.yaml`, `replacement_only.yaml`, `one_population_selection.yaml`,
  `shared_selection.yaml`, `opposing_selection.yaml`,
  `selection_aligned_turnover.yaml`.
- **Consumes:** nothing (declarative). **Produces:** inputs to
  `run_simulations.py`/`simulator.py`. **Run:** passed to the CLI above.
  `[seeded]` (each pins an explicit seed; seeds 10401–10406). `[synthetic]`.
  Together they separate pure turnover from locus-specific selection.

### Generated (not committed) simulator artifacts
- `sim/results/20261008-043318-python3/manifest.json` — an example run manifest
  produced by the provenance harness during development. `sim/output/` is a
  suggested scratch output parent. `analysis/runs/` (below) is gitignored.

---

## 3. `analysis/` — data QC, trajectory inference, locus analysis, reanalysis

**Current contents (no executable code yet):**
- `analysis/.gitignore` — ignores `/analysis/runs/`.
- `analysis/runs/` — provenance manifests and validation outputs written by the
  scripts above (gitignored). Present examples include
  `20261008-0523-sim-stat-validation/` (counts/slopes/summary) and
  `20261008-0528-sim-ctrl-null-validation/` (counts/metadata/rank/match/calibration),
  plus Lean-build and recon manifests.

**Planned modules (`TBD`, `[real-data]` unless noted):**
| File | Owner task | Consumes | Produces |
| --- | --- | --- | --- |
| `analysis/qc.py` | `data-geno` | Akbari/AADR genotypes + metadata | QC-passed dosage + QC fields |
| `analysis/trajectory_infer.py` | `traj-prototype` | test-locus-excluded genome-wide genotypes + dates | trajectory labels, time support |
| `analysis/trajectory_alt.py` | `traj-alt` | same | alternative labels (kNN/GMM) |
| `analysis/traj_relatedness.py` | `traj-relatedness` | labels + background relatedness | `family` map and `R` matrix for Model 2 |
| `analysis/canonical_loci.py` | `locus-canonical` | labels + counts at known loci | per-locus replication results |
| `analysis/locus_amy.py` (`TBD` path) | `locus-amy` | AMY copy-number/haplotype tables (H1/H2) + labels | `docs/H-amy-analysis.md` inputs |
| `analysis/locus_akbari.py` (`TBD` path) | `locus-akbari` | Akbari summary statistics + replication scores | four-quadrant classification |
| `analysis/figures.py` (`TBD` path) | `locus-figures` | analysis outputs | figures in `fig/` |

Regardless of the exact filenames (to be fixed when each module is written), all
`analysis/` modules must: run from the repo root, import Stage-B statistics via
`sim/`, record a provenance manifest, and be `[seeded]` if they draw random samples.
Their data prerequisites are the `data/raw/` acquisitions in
`docs/E-data-inventory.md` (Akbari Dataverse statistics and per-chromosome BCF/PLINK2
files, AADR v66.p1 1240K/HO genotypes and `.anno` metadata) plus the AMY tables.

---

## 4. `lean/` — formalization (Lean 4 + Mathlib v4.34.1)

The canonical Lake package lives in `lean/`. There is also a repository-root Lake
entry point (`lakefile.lean`, `lake-manifest.json`, `lean-toolchain`) that builds
the same `lean/TrajectorySelection.lean`. Both pin Lean `v4.34.1` and Mathlib
`v4.34.1`; the checked-in manifests lock Mathlib and transitive revisions.

### Build
```bash
cd lean && lake exe cache get && lake build      # canonical project
# or, from the repository root:
lake exe cache get && lake build
```
`lake exe cache get` (optional for correctness, faster/memory-lighter) downloads
prebuilt Mathlib; building from source uses more time and memory. A successful build
validates import/elaboration only, **not** scientific calibration. Lean builds are
`[deterministic]` (no RNG; provenance conventionally records seed `0`) and
`[synthetic]`.

### Library modules (`lean/`)
| File | Role |
| --- | --- |
| `lean/lakefile.lean` | canonical Lake package (`TrajectorySelection` lib). |
| `lean/lean-toolchain` | pins `leanprover/lean4:v4.34.1`. |
| `lean/lake-manifest.json` | locked Mathlib/transitive revisions. |
| `lean/TrajectorySelection.lean` | root module; imports all `TrajectorySelection/*` and states the `sorry`/`axiom` policy. |
| `lean/TrajectorySelection/Defs.lean` | `TemporalTrajectory`, `AlleleFrequencyPath`, `weightedMeanTime/Frequency`, `sampleWeightedSlope`, `directionSign`, `directionalConcordance`, `ReplacementOnly`, `TrajectoryCovariance`. |
| `lean/TrajectorySelection/Statistic.lean` | concordance counts, `correlationMass`, `effectiveFamilyCount`/`effectiveTrajectoryCount`, duplication invariant, `N_eff` bounds/antitonicity, fair-sign binomial null mean/variance. |
| `lean/TrajectorySelection/Replacement.lean` | replacement vs replicated-change sign experiment; `highConcordanceMass` and the strict separation theorem. |
| `lean/TrajectorySelection/Collinearity.lean` | time–GRM alignment, `residual_energy_cosine`, `workingInformation`, VIF identities and bounds. |
| `lean/TrajectorySelection/Core.lean` | inherited scaffold types (`Trajectory`, `Direction`). |
| `lean/TrajectorySelection/Trajectory.lean` | inherited `direction_of_change`, `consistently_directed`, `rate_of_change`. |
| `lean/TrajectorySelection/FrequencyPath.lean` | inherited `FrequencyPath` + provisional `concordance_score` placeholder. |
| `lean/TrajectorySelection/Concordance.lean` | inherited `check_replication`/`statistical_significance` **placeholders**. |
| `lean/AXIOMS.md` | mandatory ledger for `sorry`/`axiom`/`native_decide` and for inherited scientific placeholders. |
| `lean/Test.lean`, `lean/TestDefs.lean`, `lean/TestStatistic.lean`, `lean/TestReplacement.lean`, `lean/TestCollinearity.lean` | in-library example checks (some `#print axioms` for trust inspection). |
| `lean/test_simple.lean` | trivial `print "hello"` smoke file. |

**Producer:** these files are consumed by `lake build`; they produce compiled
olean/IR artifacts under `lean/.lake/` (not committed). **`sorry`/`axiom` policy:**
`lean/AXIOMS.md` must be updated in the same commit as any new `sorry`, `axiom` or
`native_decide`; the inherited placeholders in `FrequencyPath.lean` and
`Concordance.lean` must not be cited as calibrated statistics.

---

## 5. `fig/` — regenerable figures

`fig/method_diagram.svg` is committed as the Stage A→B schematic maintained
alongside `docs/D-method.md`. Other principal figures are not yet present; they will
be regenerated by `analysis/figures.py` (`locus-figures`, `TBD`) from `analysis/` and
`sim/` outputs; the exact commands are `TBD` until that task exists.

---

## 6. Data prerequisites and dependency summary

- **`[synthetic]` (no downloads):** all `scripts/` utilities; `sim/trajectory_stat.py`,
  `sim/empirical_null.py`, `sim/akbari_glmm.py` and their tests; the two `validate_*`
  scripts; the demographic simulator and scenarios; the whole `lean/` project (needs
  a network only for `lake exe cache get`, which is optional).
- **Python deps:** numpy, pandas, scipy; pyyaml for the simulator; pytest for the
  test suites. The provenance harness probes optional packages and records whatever
  is installed.
- **`[real-data]` (`TBD`, owned tasks):** the Stage-A/Stage-B/reanalysis pipeline
  needs the artifacts inventoried in `docs/E-data-inventory.md` (Akbari Dataverse
  summary statistics and per-chromosome BCF or the smaller PLINK2 set; AADR v66.p1
  1240K/HO genotypes and `.anno` metadata; AMY tables H1/H2 for the case study),
  stored immutably under `data/raw/` with checksums and access dates recorded.
- **Provenance:** every run that produces a result must be wrapped so a manifest is
  written under `analysis/runs/`, recording seed, versions, inputs and outputs
  (`docs/PROVENANCE.md`).


## Third-pass prioritisation executable — 2026-10-09 (append only)

[observation] New source: analysis/prioritization_score.py; source-named regression tests: analysis/test_prioritization_score.py. Registered protocol: docs/prioritization-prereg.md. Primary output: data/derived/priority_loci.tsv (10,641 ranked variants, all component scores and alternative carriage ranks); diagnostics under data/derived/prioritization/; figure fig/prioritization_validation.png.

[observation] Reproduce from raw inputs with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python analysis/prioritization_score.py; --reuse uses the generated panel/sample/genotype cache. Defaults: 300 independent simulated neutral loci per mechanism and30 all-panel bootstraps, seed20261009. The large rebuildable genotypes.npy cache is ignored locally; no raw input is changed. Cached variant indices are verified against exact pvar allele IDs before ranking. No git command is executed by the scorer; scripts.provenance git probes are disabled explicitly.

[observation] New regressions pass7/7, measured-mixture9/9, source encoding/path/process11/11; the actual AADR ALT orientation check also passes. This is bounded prioritisation, not a 52-million-variant whole-genome scan or a phenotype-association implementation.
