# Trajectory validation (Decision 1)

> Owner lane: traj-prototype (v1) -> traj-v2. Status: **v2 — PARTLY STABLE; the
> three diagnostic West Eurasian paths emerge post-hoc. Locus analysis may
> proceed only on v2 trajectories, with the T4 caveat below.**

Reproduce (runs in ~104 s, well under the 12-minute gate):

```
.venv/bin/python analysis/trajectory_v2.py --n-sites 100000 --seed 1
```

All numbers below are the printed output of that command and are mirrored in
`data/derived/traj_v2_metrics.json`.

---

## v1 failure (why it was revised)

v1 (`analysis/trajectory_infer.py`) clustered a genetic kNN graph with a
time-window gate using spectral clustering on `k` clusters. It produced one blob
of **9,648 / 10,369** individuals spanning **0-65,000 BP**, plus fragments of
size 5 and 2. `[observation]` That is a single dominant connected component plus
fragments, not population paths through time. The apparent bootstrap ARI (0.814)
was meaningless because a degenerate single-bucket clustering is stable by
construction. Decision 1 was recorded as **REVISE**.

Root causes diagnosed (v1 doc): modern individuals (date 0) glued the graph
together; no within-trajectory continuity requirement (time gate only prevented
*contemporaneous* merging); `k` controls cluster count, not interpretability.

---

## v2 method

Structural changes (not a re-tune of `k`):

1. **Moderns excluded.** Only non-modern individuals (`date_method` not starting
   `Modern`, `date_bp != 0`) with a usable date, usable coordinates, inside the
   West Eurasian window (lat 25-75, lon -15-65) and 0 < `date_bp` <= 50,000 are
   used.
2. **Coverage floor** `--min-present 0.10`: individuals must have >= 10% of the
   selected sites called. `[observation]` Kept 11,934 / 14,485.
3. **Ordered paths, not clusters.** The time axis is cut into overlapping bins
   (1200-yr width, 400-yr step -> 125 bins). Each bin is clustered in PCA space
   with seeded k-means (`k = clip(n//40, 1, 6)`). Clusters are linked across
   **adjacent** bins by genetic centroid distance under a hard gate
   (`gate = 1.0 * median adjacent-bin centroid distance`), using a greedy
   one-to-one assignment. A trajectory is a chain of linked clusters.
   Constraints: min trajectory size 50, max total span 8000 yr, and no chain may
   skip a bin (links are only between adjacent bins; over-long chains are split
   so that no segment exceeds the span cap).
4. **Assignment.** Overlapping bins mean an individual can appear in up to three
   bins; each individual is assigned from its *home bin* (bin whose centre is
   nearest its date), so every individual gets exactly one `trajectory_id`.
5. **Region exclusion.** `build_embedding` accepts an exclusion list; the default
   removes chr6 MHC (28.5-33.5 Mb) and chr1 AMY1/AMY2 (103.5-104.5 Mb, GRCh37).

### Pipeline counts (v2 run)

| quantity | value |
|---|---|
| individuals (1240K patch) | 23,089 |
| sites (1240K patch) | 1,233,013 |
| modern excluded (`date_method` Modern / date 0) | 3,970 |
| no coordinates | 90 |
| outside West Eurasian / date window (incl. archaic >50k) | 4,544 |
| usable non-modern, dated, geolocated | **14,485** |
| candidate autosomal non-excluded sites | 190,547 |
| selected sites (MAF>=0.05, presence>=0.5) | 64,224 |
| coverage floor 0.10: kept | **11,934** |
| time bins | 125 |
| cluster nodes | 214 |
| trajectories kept | **13** |
| unassigned individuals | 152 |
| PCA variance explained PC1..PC10 | 0.0066, 0.0044, 0.0027, 0.0013, 0.0011, 0.0010, 0.0006, 0.0006, 0.0005, 0.0005 |

`[observation]` The gate is 48.74 (PC units, = median adjacent-bin centroid
distance).

### Cap checks (asserted in code, printed)

```
CHECK span cap: max observed span=7600 <= 8000 OK
CHECK min size: min observed size=124 >= 50 OK
```

### Per-trajectory summary

`coh` = mean within-trajectory PC-space pairwise distance / median
cross-trajectory mean (values < 1 are tighter than the cross-trajectory median).

| traj | n | date min | date max | median | span | coh |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 2,254 | 100 | 6,500 | 1,783 | 6,400 | 1.27 |
| 1 | 216 | 500 | 1,552 | 1,225 | 1,052 | 0.87 |
| 2 | 752 | 168 | 3,514 | 1,570 | 3,346 | 0.95 |
| 3 | 339 | 150 | 2,051 | 1,225 | 1,901 | 0.89 |
| 4 | 4,594 | 100 | 7,700 | 2,300 | 7,600 | 0.60 |
| 5 | 238 | 7,708 | 14,021 | 9,254 | 6,313 | 2.18 |
| 6 | 124 | 1,009 | 1,566 | 1,225 | 557 | 1.00 |
| 7 | 1,255 | 1,716 | 5,300 | 3,791 | 3,584 | 0.77 |
| 8 | 319 | 2,951 | 6,099 | 4,356 | 3,148 | 1.28 |
| 9 | 486 | 3,725 | 8,895 | 5,160 | 5,170 | 1.71 |
| 10 | 884 | 5,385 | 8,850 | 6,654 | 3,465 | 0.89 |
| 11 | 191 | 5,746 | 10,451 | 7,970 | 4,705 | 1.66 |
| 12 | 130 | 6,537 | 8,100 | 7,252 | 1,563 | 1.44 |

`[observation]` No trajectory spans the whole record, but **T4 (4,594 = 38.5%)
still spans 7,600 yr** — it hit the span cap and is the residual limitation of
this run (see verdict). It is nonetheless genetically tighter than the
cross-trajectory median (`coh`=0.60).

---

## Stability

| metric | value |
|---|---|
| bootstrap ARI (20 individual resamples, seeded) | mean **0.661**, sd **0.138** |
| sensitivity: bin width +400 yr | ARI **0.655** |
| sensitivity: bin step 800 yr | ARI **0.659** |
| sensitivity: gate x0.7 | ARI **0.999** |
| sensitivity: gate x1.5 | ARI **0.912** |

`[interpretation]` Bootstrap ARI (0.66) and width/step sensitivity (~0.66) are
moderate: the broad path membership is reproducible but individual-level labels
near path boundaries shift. Gate sensitivity is low (x0.7 almost identical; x1.5
drops to 0.91), i.e. the result is not knife-edge sensitive to the gate. This is
a genuine improvement over v1, whose parametric ARI was 0.519 on a degenerate
partition. It is not high enough to claim a single canonical labelling.

---

## Known-structure comparison (POST-HOC annotation only)

`[observation]` The published `Group ID` from `.anno` was **never used as input**;
it is used below only to name what the algorithm rediscovered. Counts are
keyword-tagged group memberships (WHG-like / early-farmer-like / steppe-like).

| traj | n | WHG-like | EEF-like | steppe-like | dominant groups |
|---:|---:|---:|---:|---:|---|
| 0 | 2,254 | 2 | 21 | 33 | Belgium_HighMedieval, Sweden_Viking, England_Saxon |
| 1 | 216 | 0 | 0 | 0 | Austria_Avar, Hungary_*Avar-oLowEastAsia |
| 2 | 752 | 0 | 13 | 0 | Turkey_Medieval_Byzantine, Hungary_LateAvar, Italy_Lazio_ImperialRoman |
| 3 | 339 | 0 | 0 | 0 | Hungary_*Avar-oLowEastAsia, Austria_Avar |
| 4 | 4,594 | 22 | 192 | 119 | Austria_Avar, Belgium_HighMedieval, Hungary_Avar |
| 5 | 238 | 98 | 6 | 0 | Iran_GanjDareh_N, Serbia_IronGates_Mesolithic, Turkey_Catalhoyuk |
| 6 | 124 | 0 | 0 | 0 | Austria_Avar, Hungary_*Avar-oHighEastAsia |
| 7 | 1,255 | 2 | 35 | 197 | England_MIA, Czechia_EBA_Unetice, Czechia_EBA_CordedWare |
| 8 | 319 | 0 | 9 | 2 | Turkey_MLBA, Turkey_C, Turkmenistan_C |
| 9 | 486 | 25 | 174 | 224 | Russia_Saratov_Eneolithic_Khvalynsk, Russia_Samara_EBA_Yamnaya |
| 10 | 884 | 1 | 220 | 0 | Austria_N_LBK, Romania_C_Bodrogkeresztur, France_Yonne_N |
| 11 | 191 | 147 | 0 | 0 | Serbia_IronGates_Mesolithic, Latvia_EMN_Narva, Spain_Mesolithic |
| 12 | 130 | 21 | 24 | 1 | Ukraine_N, Ukraine_Kirovohrad_N |

`[observation]` **All three diagnostic paths emerge without being used as
input:**

- **WHG-like:** T11 (147/191 WHG-like; Iron Gates Mesolithic, Latvia Narva,
  Spain Mesolithic) and, more weakly, T5 (98/238 WHG-like).
- **Early-farmer-like:** T10 (220/884 EEF-like; Austria LBK, Romania
  Bodrogkeresztur, France early Neolithic).
- **Steppe-like:** T9 (224/486 steppe-like; Khvalynsk, Samara Yamnaya) and T7
  (197/1255; CordedWare, Unetice, England MIA).

`[observation]` Additional coherent, non-diagnostic paths: Avar-era groups with
East Asian admixture split into T1/T3/T6; Anatolian/Iranian Bronze-Iron (T8);
and a large historical European path (T0/T4). T4 mixes WHG/EEF/steppe/steppe-ERA
material over 7,600 yr and is the least resolved path.

`[inference]` This is precisely the geometry the v1 blob destroyed: ordered
chains of contemporaneous clusters that track real West Eurasian structure.
Excluding moderns and requiring adjacency-linked ordered paths was necessary and
sufficient for these three paths to appear.

---

## Decision 1 verdict: **PARTLY STABLE**

Evidence for moving forward to locus analysis:

- v1's failure mode is gone: no trajectory spans the full record, none is a
  vestigial fragment; min size 124 >= 50, max span 7,600 <= 8,000, both asserted.
- The three paths the design asked about (WHG-like, early-farmer-like,
  steppe-like) emerge **post-hoc** from genotypes alone (Section above).
- Stability is moderate and honest: bootstrap ARI 0.66, width/step ARI ~0.66,
  gate ARI 0.91-1.00 — far above v1's degenerate 0.519 parametric ARI.

Evidence against calling it fully stable:

- **T4 absorbs 38.5% of individuals and hits the 7,600-yr span cap.** It is
  coherent (coh 0.60) but under-resolved; it likely merges several historical
  European sub-paths.
- Individual-level labels are ~34% unstable under bootstrap/width perturbation,
  so trajectory boundaries (not path centres) are soft.
- Trajectories are *paths* with a 8,000-yr cap by design; they are not full
  lineage trees, and short spans make single trajectories coarse for trend
  tests.

**What is missing to reach "stable":** a principled split of T4 (e.g. a tighter
gate inside the dense Holocene, or recursive sub-binning) and a formal
path-significance test (e.g. vs a date-shuffled null). Until T4 is resolved,
locus trends should be read per-path with T4 treated as a mixed background, and
the AMY/canonical/Akbari analyses should report sensitivity to excluding T4.

---

## Artifacts

- v2 code: `analysis/trajectory_v2.py`
- v2 output table: `data/derived/trajectories.tsv`
  (columns `individual, date_BP, lat, lon, group_id, PC1..PC10, trajectory_id, bin_id`; `trajectory_id = -1` = unassigned)
- v2 metrics: `data/derived/traj_v2_metrics.json`
- v2 figures: `fig/traj_v2_pca.png`, `fig/traj_v2_spans.png`, `fig/traj_v2_stability.png`
- v1 code + output retained for comparison: `analysis/trajectory_infer.py`,
  `data/derived/trajectories_v1.tsv`, `fig/traj_pca_time.png`

`[observation]` Downstream caveat: `analysis/amy.py::load_trajectory_table`
currently requires a `date` column, but this table (per the v2 output schema)
uses `date_BP`. `amy.py` must be updated to accept `date_BP` before
`--trend data/derived/trajectories.tsv` will run; the labels themselves are
present and usable.
