# Subset-selection regime test - results (against the pre-registration)

Status: **Tier-1 scouting, negative result.** This reports the pre-registered
subset-regime test in `docs/subset-regime-prereg.md` (written before any result
was computed). All claims are labelled [observation] / [inference] / [speculation].

## What was run

[observation] Command:
`.venv/bin/python analysis/subset_regime.py --max-sites 8000 --chunk-sites 2000`.
Panel and genes/labels reused from `analysis/reanalysis_scan.py`; relatedness
matrix `data/derived/traj_relatedness.tsv` applied (so `relatedness_applied=true`).
Replication units: the **13 trajectories** in `data/derived/trajectories.tsv`
(trajectories_v3.tsv does not exist; see pre-registration (a)). Statistic: the
max-over-subsets T = max over subsets S (2 <= |S| <= 4) of the GLS meta |z|
restricted to S. Null: matched control AADR sites scored with the **same**
max-over-subsets statistic. Runtime 295 s.

[observation] Panel: 8,000 loci = 2,000 Akbari-significant candidates
(|X| >= 3.16) + 6,000 background loci (|X| < 3.16). All 2,000 candidates had a
finite empirical P (median 41 matched controls each, min 20). 62
(trajectory,bin) cells, 11,782 trajectory individuals, unmapped=0.

## Pre-registered decision rule vs observed

Pre-registered rule (from the pre-registration (d)): empirical upper-tail P vs
matched controls, then **BH-FDR at alpha = 0.05**, minimum subset size 2.

[observation]
* **Loci passing BH-FDR q <= 0.05: 0 of 2,000.** Minimum q = 0.639.
* Candidates with nominal empirical P <= 0.05: **153** (chance expectation ~100).
* Loci passing that are low-|X| (|X| < 3.16): **0** (the tested set is the
  Akbari-significant set by design, so low-|X| loci are the null pool, not tested
  for discovery; with zero passes there are also zero low-|X| passes).
* Candidates with a proper-subset argmax and a flat/opposite complement: 1,137
  of 2,000 - but **none** of them pass FDR.
* Nominal-P <= 0.05 **and** the subset signature (proper argmax, flat/opposite
  complement): 17 candidates, none survive FDR.

[observation] T is slightly higher for candidates than for matched background:
candidate mean 7.153 / median 6.153 vs background mean 6.310 / median 5.719
(mean ratio 1.13).

## Verdict against the pre-registered prediction

Pre-registered falsifier (a): "the prediction is falsified if no candidate passes
BH-FDR at alpha 0.05."

[inference] **The pre-registered prediction is FALSIFIED by falsifier (a).**
Zero of 2,000 Akbari-significant loci pass BH-FDR at alpha 0.05; the best q is
0.639. There is a mild, directionally consistent enrichment of T in the
candidate set (falsifier (b) is *not* met: candidates are higher than controls),
but it is nowhere near enough to survive correction for the max-over-subsets
multiple-testing device. There is **no discovery list**, and in particular no
low-|X|/high-replication locus.

[observation] Mechanistic reason the max-over-subsets statistic did not isolate
the target regime: the argmax subset size is **4 (the maximum allowed) for
1,736/2,000 candidates and 5,017/6,000 background loci**. The winning subset is
almost never a small coherent block; the maximum over subsets behaves mostly as
a multiple-testing inflation that lifts T for candidates *and* controls together.

[inference] This is the same structural conclusion as the first pass
(`docs/I-quadrant-first-pass.md`), reproduced under a statistic specifically
engineered to expose subset-restricted selection: the pooled-GLMM dilution story
does not, on this Tier-1 panel, translate into a calibrated subset signal that a
matched max-over-subsets null accepts.

[speculation] Possible reasons, untested here: the 13 trajectories are
over-split / partly stable, so the "subset" that selection acts on is smeared
across several correlated trajectories rather than a clean 2-4 group block; and
Tier-1 attenuation (26% of individuals, pseudo-haploid low coverage) weakens the
per-subset z before the max is taken.

## Deviations from the pre-registration

[observation] No substantive deviation. The pre-registered procedure was executed
as written (units = 13 trajectories; subsets 2 <= |S| <= 4; null = matched
controls scored with the identical max-over-subsets statistic; empirical upper-tail
P; BH-FDR alpha 0.05; minimum subset size 2). The only operational notes:
1. The panel was capped at 8,000 loci (rather than genome-wide) for the runtime
   budget; the cap and the 2,000/6,000 split are recorded in
   `data/derived/subset_regime_metrics.json`. The analysis itself is not
   exploratory on this account; the cap only limits power.
2. The figure `fig/subset_regime.png` was produced by an inline script rather than
   a committed plotting module; no pre-registered quantity was changed.

## Honest limitations

* Tier-1 panel; both statistics attenuated; trajectories partly over-split
  (see `docs/I-quadrant-first-pass.md`).
* The max-over-subsets statistic is anti-conservative by construction; only the
  matched max-over-subsets null is used, so the *comparison* is valid, but the
  absolute T values are not interpretable as genome-wide significance.
* This is a targeted test of one regime, not a genome-wide discovery scan.

## Files

* Procedure: `analysis/subset_regime.py` (writes `data/derived/subset_regime.tsv`
  incrementally per genotype chunk).
* Panel metrics: `data/derived/subset_regime_metrics.json`.
* Per-locus table: `data/derived/subset_regime.tsv` (8,000 rows).
* Figure: `fig/subset_regime.png`.
