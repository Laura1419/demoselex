# Trajectory merge v3 - genetically-independent replication units

> Owner lane: traj-v3-merge. Producer: `analysis/trajectory_merge.py`. Status: **merged grouping produced; see verdict.**

Reproduce (single command, writes all four outputs):

```
.venv/bin/python analysis/trajectory_merge.py
```

## Why this file exists

The v2 pipeline returns 13 trajectories, but several are the *same* population history split by sampling period or admixture fraction, so counting them as 13 independent replications double-counts shared history. Two splits were diagnosed independently (`docs/traj-relatedness.md`):

* T1/T3/T6 - Avar-period Hungary split by oLowEastAsia vs oHighEastAsia ancestry proportion (pairwise r 0.58-0.79);
* T0/T4 - medieval/Avar-period Europe, r(0,4)=0.858 (T4 is the 4,594-individual blob that hit the 8,000-yr span cap).

## 1. Threshold sweep (union-find on the existing 13x13 matrix)

`[observation]` Merging all trajectory pairs whose relatedness coefficient exceeds the threshold, by union-find, on `data/derived/traj_relatedness.tsv`. `N_eff = m^2 / sum(R)` is computed on the merged groups (cluster-to-cluster relatedness = mean of member pair coefficients). Every tested threshold is shown:

| threshold | m groups | grouping | N_eff |
|---|---:|---|---:|
| r>0.9 | 13 | {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} | 4.781 |
| r>0.7 | 11 | {0,4} {1,6} {2} {3} {5} {7} {8} {9} {10} {11} {12} | 4.944 |
| r>0.5 | 6 | {0,10,2,4,7,8} {1,3,6} {5} {9} {11} {12} | 4.222 |
| r>0.4 | 4 | {0,10,2,4,7,8} {1,3,6} {5} {11,12,9} | 4.828 |
| r>0.3 | 2 | {0,10,11,12,2,4,5,7,8,9} {1,3,6} | 3.774 |
| (no merge) 13 distinct | 13 | all singleton | 4.781 |

`[observation]` The sweep is monotone in coarseness but not in N_eff: collapsing two near-duplicates can tick N_eff *up* slightly before it falls, an arithmetic property of m^2/sum(R) already documented in `docs/traj-relatedness.md`.

## 2. Is there a natural break in the relatedness values?

`[observation]` Sorted off-diagonal coefficients (28 nonzero pairs) and the largest consecutive gaps:

| rank | pair | r | gap above |
|---:|---|---:|---:|
| 1 | 2-8 | 0.6969 | 0.0967 |
| 2 | 0-8 | 0.2627 | 0.088 |
| 3 | 1-6 | 0.7935 | 0.064 |
| 4 | 5-12 | 0.1670 | 0.0638 |
| 5 | 9-11 | 0.0881 | 0.0577 |
| 6 | 8-10 | 0.4319 | 0.05 |
| 7 | 4-7 | 0.5286 | 0.0471 |
| 8 | 2-4 | 0.6248 | 0.0415 |

`[observation]` The largest gap is **0.0967** (between 2-8 at r=0.6969 and the next pair). The next largest gaps are ~0.088, ~0.064, ~0.064, ~0.058 - all comparable. **There is no sharp natural break**: the spectrum is continuous with only modest, similar-sized steps.

`[inference]` A gap-maximising rule would put the cut near r~0.7 (between the 0.794 and 0.697 pairs). But that threshold only merges {0,4} and {1,6} and would leave T3 *separate* from T1/T6 - i.e. it contradicts the independently diagnosed Avar split, whose internal T3-T6 coefficient (0.576) sits below the gap. So the gap structure does **not** single out a merge set: a purely data-driven cut and the diagnosis disagree. This is stated plainly rather than picked to give a pleasing number.

## 3. Chosen grouping (operating threshold r>0.5)

`[inference]` Operating threshold **r>0.5** is chosen on two grounds that are independent of the resulting count:

* it is the *finest tested* threshold that reunites the diagnosed T1/T3/T6 Avar cluster (needs r < 0.576, satisfied by 0.5 but not 0.7);
* coarser thresholds (0.4, 0.3) add merges that mix clearly distinct histories - at r>0.4 the steppe-like T9 joins the WHG-like T11/T12, and at r>0.3 almost everything collapses into 2 groups.

`[observation]` Chosen groups:

| group | member trajectories | name (descriptive) | n | date min-max (BP) | mean lat/lon | dominant published groups |
|---|---|---|---:|---|---|---|
| G0 | 0,2,4,7,8,10 | broad west/central-Eurasian Holocene bloc (medieval + Bronze-Age + early-farmer + steppe-adjacent material) | 10058 | 100-8850 | 47.8/15.5 | Austria_Avar (433), Belgium_HighMedieval (183), Sweden_Viking (123) |
| G1 | 1,3,6 | Avar-period Hungary, East-Asian-admixed | 679 | 150-2051 | 48.4/24.5 | Austria_Avar (116), Hungary_LateAvar-oLowEastAsia (84), Hungary_LateAvar-oHighEastAsia (44) |
| G2 | 5 | Iran/Anatolia Neolithic + Iron-Gates hunter-gatherer (heterogeneous) | 238 | 7708-14021 | 43.1/29.8 | Iran_GanjDareh_N (15), Serbia_IronGates_Mesolithic (14), Turkey_Catalhoyuk_MN (14) |
| G3 | 9 | steppe-like (Khvalynsk / Yamnaya) | 486 | 3725-8895 | 51.7/42.4 | Russia_Saratov_Eneolithic_Khvalynsk (38), Russia_Samara_Eneolithic_Ekaterinovka (37), Russia_Samara_EBA_Yamnaya (34) |
| G4 | 11 | WHG-like hunter-gatherer (Iron Gates / Narva / Spain Mesolithic) | 191 | 5746-10451 | 49.8/16.8 | Serbia_IronGates_Mesolithic (51), Latvia_EMN_Narva (17), Latvia_LateMesolithic_Kunda (11) |
| G5 | 12 | Ukrainian Neolithic (mixed WHG/EEF/steppe) | 130 | 6537-8100 | 49.9/33.4 | Ukraine_Kirovohrad_N (22), Ukraine_N (22), Ukraine_Zaporizhia_N (14) |

`[observation]` `data/derived/trajectories_v3.tsv` uses the SAME rows as `trajectories.tsv` with one added column **`merged_group`** (`G0`..`G5`, or `unassigned` for `trajectory_id == -1`). The existing `group_id` column is the AADR *published-group* annotation and is deliberately NOT reused; `merged_group` is the new merged label. All v2 columns and names are preserved, and downstream readers of `individual`, `date_BP`, `trajectory_id` are unaffected.

## 4. Recomputed relatedness on the merged groups

`[observation]` The v2 estimator (`analysis.traj_relatedness.group_accumulate` / `relatedness_from_accum`) was re-run with the merged labels, pooling all individuals per group, and written to `data/derived/traj_relatedness_v3.tsv`. Sites used: 148460; individuals: 11782.

`[observation]` The matrix passes the consumer contract `sim.trajectory_stat._correlation` (symmetric, unit diagonal, entries in [0,1], Cholesky PD): PASS.

Merged-group matrix:

| | G0 | G1 | G2 | G3 | G4 | G5 |
|---|---|---|---|---|---|---|
| G0 | 1.00 | 0.13 | 0.17 | 0.00 | 0.00 | 0.00 |
| G1 | 0.13 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| G2 | 0.17 | 0.00 | 1.00 | 0.00 | 0.00 | 0.00 |
| G3 | 0.00 | 0.00 | 0.00 | 1.00 | 0.00 | 0.06 |
| G4 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | 0.18 |
| G5 | 0.00 | 0.00 | 0.00 | 0.06 | 0.18 | 1.00 |

`[observation]` The recomputed group matrix is nearly block-diagonal: the only nonzero off-diagonals are G0-G1=0.13, G0-G2=0.17, G3-G5=0.06, G4-G5=0.18; all others are clipped to 0. Pooling many trajectories into one group's mean allele-frequency profile averages away the shared drift the 13-way matrix showed (r(0,4)=0.86), so the grouped estimator *under*-states relatedness here. This is the opposite failure mode to over-calling, and it is why the recomputed N_eff is higher than the coarse-v2 N_eff.

## 5. Final N_eff and honesty verdict

`[observation]` On the merged groups, `N_eff = m^2/sum(R)` computed from the **recomputed** v3 matrix = **5.086** (m=6). For comparison, `N_eff` on the raw 13 = 4.781, and on the merged groups using the coarse of the v2 matrix = 4.222.

`[inference]` **Split of the 6 merged groups into distinct histories vs sampling artefacts:**

* **Plausibly distinct histories (3-4):** G1 (Avar-period East-Asian-admixed), G3 (steppe-like: Khvalynsk/Yamnaya), G4 (WHG-like: Iron Gates/Narva/Spain Mesolithic), and arguably G5 (Ukrainian Neolithic).
* **Still sampling artefacts / over-merged (1-2):** **G0 is not one history.** It pools medieval/Avar-period Europe (T0/T4) with Anatolian/Iranian Bronze-Iron (T8), steppe CordedWare/Unetice material (T7), early-farmer LBK (T10) and Byzantine Anatolia (T2) into a single 'broad bloc' - exactly the period/region sampling artefact this task set out to remove. G2 (Iran Neolithic + Iron Gates + Catalhoyuk, n=238) is small and heterogeneous and may also be a sampling artefact.

`[observation]` The recomputed merged `N_eff` (5.086) is *above* the raw-13 value (4.781) - the same m^2/sum(R) arithmetic seen in `docs/traj-relatedness.md` when a tight pair ({0,4}) is collapsed. The recomputed matrix (pooled group means) and the coarse average of the v2 matrix (4.222) disagree, which is itself a warning that N_eff is not stable to the estimation route.

`[inference]` Across every route N_eff stays in the ~4-5 band and is dominated by the pervasive 0.2-0.6 background relatedness. Because G0 is demonstrably an over-merge, treating it as one independent unit **still overcounts** replication: the honest replication count is the number of *genuinely distinct* histories (~3-4), not 6 and not 13.

`[speculation]` The G0 bloc is consistent with a shared west/central Eurasian deep ancestry sampled at many times; splitting it cleanly would need admixture/outgroup-f3 modelling or a tighter within-Holocene gate, which is out of scope here.

## Verdict: is this grouping fit to be the replication unit?

`[inference]` **PARTLY.** The merged `trajectories_v3.tsv` / `traj_relatedness_v3.tsv` are a correct, reproducible collapse of the 13 v2 trajectories into 6 relatedness groups and are safe to pass as the `family`/`relatedness` inputs to `sim.trajectory_stat`. But they are **not yet a defensible replication count on their own**: the grouping fixes the diagnosed T1/T3/T6 and T0/T4 splits only to the extent that r>0.5 merges them, while the resulting G0 still bundles several distinct histories, so the merged set under-corrects. The load-bearing number remains **N_eff ~ 4-5** (consistent with `docs/traj-relatedness.md`); report locus concordance with `n_effective ~= 4` (the count of genuinely distinct histories), and treat the G0 bloc as mixed background rather than as one independent test.

## Artifacts

* code: `analysis/trajectory_merge.py`
* merged table: `data/derived/trajectories_v3.tsv` (v2 columns + `merged_group`)
* merged matrix: `data/derived/traj_relatedness_v3.tsv`
* metrics: `data/derived/traj_v3_metrics.json`
* figures: `fig/traj_v3_relatedness.png`, `fig/traj_v3_groups.png`

