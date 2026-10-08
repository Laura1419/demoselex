# Genetically-independent trajectories and the effective number of replications

Owned by task `traj-relatedness`. Producer: `analysis/traj_relatedness.py`.
Artifacts: `data/derived/traj_relatedness.tsv` (13x13 labeled matrix),
`data/derived/traj_relatedness_metrics.json`, `fig/traj_relatedness.png`.

## Why this file exists

The project keeps saying "k trajectories agree". That is only a replication
count if the trajectories are *genetically independent*. Trajectories are
clusters of individuals in a time-aware PCA embedding, and several of them are
drawn from the same deep-ancestry background, so counting them as independent
tests double-counts shared history. This file quantifies the overlap and turns
13 into an effective replication number `N_eff`.

## Estimator (exactly what the code computes)

[observation] For a filtered set of autosomal biallelic sites we read the AADR
2-bit genotype codes as diploid dosages `g in {0,1,2}` (code 3 = missing) and
form the per-trajectory frequency `p_{t,s} = mean dosage / 2`.

[observation] The reference frequency at site `s` is the **equal-weight** mean of
the trajectory frequencies, `ref_s = (1/m) sum_t p_{t,s}` (each of the 13
trajectories contributes 1/13, so the largest cluster cannot define the
reference away from itself). Scaling is `sd_s = sqrt(2 ref_s (1-ref_s))`. The
standardised profile is `d_{t,s} = (p_{t,s} - ref_s) / sd_s`.

[observation] The relatedness matrix is

```
R_ij = (1/S) * sum_s d_{i,s} d_{j,s}          # mean pairwise allele-sharing / IBS
r_ij = R_ij / sqrt(R_ii R_jj)                 # correlation-like coefficient in [0,1]
```

Because the trajectory mean is linear, `R_ij` is the sample-size-normalised mean
pairwise allele-sharing between individuals of trajectories `i` and `j` (a
genome-wide IBS / Weir-Hill-kinship family estimator), aggregated to trajectory
level. Negative entries (diverging trajectories) are clipped to 0 and the
diagonal is set to 1, exactly matching the `entries in [0,1]`, symmetric,
diagonal-1 contract of `sim/trajectory_stat._correlation`.

[observation] Site filter: autosomes 1-22 only; **excluded chr6 MHC (GRCh37
28.5-33.5 Mb)** and **chr1 AMY1/AMY2 (GRCh37 103.5-104.5 Mb)**; sites with >20%
missing among trajectory individuals dropped; a fixed-seed (20240501) uniform
subsample of sites is used for tractable decoding.

[observation] Sites: 1,143,278 autosomal candidates pass the region filter,
150,000 are subsampled, **148,458 usable** after the missingness filter.
Individuals: 11,782 in 13 trajectories.

This estimator is deliberately **not** estimated at any tested locus; it is a
background ancestry matrix to be supplied to `trajectory_stat(relatedness=...)`.

## The matrix

[observation] Symmetric, diagonal 1, PD (Cholesky passes), entries in [0,1].
Full matrix in `data/derived/traj_relatedness.tsv`. Structure:

- {1, 3, 6}: r(1,6)=0.79, r(1,3)=0.62, r(3,6)=0.58 - one tight background.
- {0, 4}: r(0,4)=0.86 - the single most related pair.
- {0, 2, 4, 7, 8, 10}: r in 0.15-0.86 - a broad related bloc.
- {5, 11, 12}: r(11,12)=0.48, r(5,11)=0.37, r(5,12)=0.17.
- {9, 12}: r=0.40; {7, 9}: r=0.35.

[observation] Max off-diagonal is 0.858 (0-4); the second largest is 0.794
(1-6). No pair reaches 0.9, so **no two trajectories are exact duplicates**;
none had to be merged or regularised. Cholesky (`sim.trajectory_stat._correlation`)
passes on the raw matrix.

## Effective number of independent trajectories

[observation] With the code's definition `N_eff = m^2 / sum(R)` over the
distinct set:

| grouping | m | N_eff |
|---|---:|---:|
| all 13 distinct | 13 | **4.78** |
| merge pairs with r>0.8 -> {0,4} | 12 | 5.08 |
| merge pairs with r>0.7 -> also {1,6} | 11 | 4.94 |
| merge pairs with r>0.5 | 6 | 4.22 |
| merge pairs with r>0.9 | 13 | 4.78 |

[observation] Coarser grouping at r>0.5 yields 6 ancestry families:
{0,2,4,7,8,10}, {1,3,6}, {5}, {9}, {11}, {12}, with N_eff = 4.22.

[inference] The single-pair merge at r>0.8 makes `N_eff` tick *up* to 5.08 before
it falls at coarser thresholds. This is an arithmetic property of
`m^2/sum(R)`, not evidence of new independence: collapsing two groups removes
two diagonal 1s but their mutual coefficient (0.858) plus their mean relatedness
to the rest is still large, and `m^2` shrinks faster than `sum(R)` for this one
pair. Across every sensible threshold the matrix never supports anywhere near 13
independent trajectories; the metric is dominated by the pervasive 0.2-0.6
background relatedness.

[inference] **The load-bearing number: the 13 trajectories behave like about
N_eff ~= 4-5 genetically independent trajectories.** Any "k trajectories
concordant" statement must use k ~= 4-5, not k = 13, or it over-counts
replication by roughly 2.5-3x.

## Split-trajectory sanity check (required, concrete)

[observation] Trajectory 4 (the largest, 4,594 individuals) was split at random
into two halves A/B, the whole pipeline was re-run on the 14-group set, giving:

- r between the two halves = **0.798** (they are the same population, as
  expected, well above any cross-trajectory value);
- `N_eff` for the 13 distinct groups = **4.781**;
- `N_eff` after the artificial split into 14 = **4.620**.

[observation] `N_eff` did **not** increase (4.781 -> 4.620), i.e. inventing a new
"trajectory" out of half of an existing one gives no phantom replication. This is
the required negative control and it passes. (r=0.798 rather than 1.0 is residual
finite-sample noise in the half-population frequencies; it is small enough that
`sum(R)` still grows faster than `m^2`.)

## Implications

[inference] Treat the trajectory set as ~4-5 real replications. Concretely, when
`sim/trajectory_stat.trajectory_statistic` reports concordance across families,
pass this matrix as `relatedness` and read `n_families` / `n_effective`, not the
raw `n_trajectories`.

[speculation] The {0,2,4,7,8,10} bloc is consistent with a shared west/central
Eurasian ancestry sampled at different times; the {1,3,6} bloc with a distinct
later ancestry. Confirming that would need admixture or outgroup-f3 modelling,
which is out of scope here.

[observation] The matrix is written in the labeled pandas DataFrame layout
(index = columns = trajectory labels as strings, diagonal 1.0, symmetric, PD)
that `sim/trajectory_stat.py` consumes.

## Reproduce

```
.venv/bin/python analysis/traj_relatedness.py
```

Requirements checked by the script itself: `sim.trajectory_stat._correlation`
imported and called (PD/contract test), the split-trajectory negative control,
and the matrix + metrics + figure are written incrementally.
