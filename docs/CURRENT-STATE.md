# CURRENT STATE (2026-10-13) — read this first

Single authoritative snapshot. If any older doc contradicts this one, **this one wins**
(older docs are historical records; the supersession chain is logged at the bottom).
Provenance: 109 commits at `a8e5d44`, all numbers reproducible from the named scripts.

---

## 1. The question

Does replicated directional allele-frequency change across independently inferred **local
ancestral states** add information beyond Akbari et al. 2026's population-structure-corrected
GLMM — specifically, does conditioning on the local ancestral state separate **carriage**
(the allele moves with its ancestry background: demography) from **selection** (the allele
moves *within* a fixed background)?

Motivating case: AMY1 copy-number, which sits below Akbari's significance threshold
(`|X|` = 3.795 vs the 5.45 cut) in the band the paper discarded.

## 2. Established results (stand independently of everything below)

These were derived before the local-ancestry work and are unaffected by it:

| finding | number | doc |
|---|---|---|
| AMY sits in the discarded band | `|X|` = 3.795, max over 2,204 window SNPs; **zero** of 31,939 strong variants in the window; LCT 16.5, FADS 16.2, SLC45A2 13.2 for contrast | `C-akbari-reconstruction.md` |
| Akbari's statistic reconstructed | X = Z/sqrt(CF), CF = 2.28; 5.45 = π>99%, 3.16 = FDR-50% (388,274 variants) | `akbari-X-units.md` |
| Binarised AMY short/long allele clears FDR-50% **and survives** the structure correction | X-equiv 3.64; ~97% GRM retention | `amy-binarized-results.md` |
| Signal is ~94% carriage (global-label decomposition) | within 0.151 / between 0.207 copies/kyr | `two-channel` section of `J-final-assessment.md` |
| Within-history change below detection floor at n=225 | 0.049 observed vs 0.146 floor | same |
| No genome-wide low-`|X|`/high-replication class | 4.3% vs 5% null; block-level top-100 has **0 novel** variants | `block-survey-results.md` |
| Peer review gap | R3 raised time–GRM collinearity; authors conceded signal "can be absorbed into the variance component"; **no referee asked the false-negative question**; threshold debate ran the other way (reviewers thought 5.45 too loose) | `akbari-peer-review-evaluation.md` |

## 3. The instrumentation arc — negative result, then reversal

This is the part most at risk of future confusion. Chronology:

1. **Attempt 1** (local ancestry from the reference panel): failed — insufficient coverage.
   **Attempt 2** (in-house painter): wrong model class (k-means/independent-SNP, no LD).
2. **Patched ELAI** instead: fork at `/home/erikg/elai-fork`, public at
   `github.com/ekg/ELAI`, branch `lower-state-posterior`. Adds `--pkk` (joint lower-state
   pair posterior, marginalised over upper S), `--pkk-site chrom:pos` (virtual
   missing-emission bridge site), and `--pk1` marginal export. Validated to ~1e-15 against
   exact enumeration; 12/12 existing outputs byte-identical. Known limitation:
   `--pkk-site` cannot place a site **outside** the panel's outermost markers (no bracket);
   interior gaps only.
3. **First verdict was "route closed"** — superseded. Two metrics errors along the way,
   both retracted in writing (`fda4ff7`, `9244b32`): a label-dependent seed comparison
   (invalid because of EM label switching) and the pre-correction convergence reading.
4. **The real cause of the apparent failure was configuration, not data** (`e365fca`):
   - **S=1**: with one source layer there *is no ancestry in the model at all* — lower
     clusters are haplotype clusters, not ancestral states. This also made `-mg` (mixture
     generations = switch rate) *mechanically* inert: `-mg 100` produced bit-identical
     output to `-mg 10`. The manual's guidance is `-c ≈ 5×-C`; we had `-C 1 -c 4`.
   - **35 Mb window**: 35–70× the ~0.5–1 Mb Neolithic tract scale. At that width the
     likelihood surface rewards 4-cluster structure over a collapsed single cluster by
     only **0.4%**, and 4/18 seeds collapsed to one effective cluster.
5. **Corrected configuration**: **S=2, K=10, 1 Mb windows, 100 exact EM steps, `-mg 10`**,
   two seeds. At this configuration states become confident, reproducible, and local
   (see §4).

**Honest residue**: the EM convergence *gate* remains noise-limited (noise floor ~5e-5 vs
tolerance 1e-4; identical configs pass or fail by seed). The convergence criterion is not
used as evidence anywhere in the final pipeline; confidence and cross-seed reproducibility
are the instruments that matter.

## 4. What exists now: the genome-wide local-ancestry state map

Produced by `analysis/elai_genome_scan.py` (modes: build/window/scan/assemble/compose),
from the Akbari genome-wide imputed bundle (same source as the paper's own scan), 1,000
dated ancients (same cohort as the dense build).

| | value |
|---|---|
| windows | **2,668 / 2,690** complete with both seeds (22 skipped: <150 QC markers) |
| configuration | 1 Mb windows, ≤800 markers/window uniform-spaced from 5.9M QC-passing variants |
| confidence | median max q = **0.94**; >50% of individuals confident (q>0.9) in **65%** of windows |
| cross-seed reproducibility | ARI median **0.62** (>0.5 in 72% of windows); **0.0%** of windows below chance |
| ancestry-bearing | median NMI(country \| state) = **0.18**; >0.10 in **92%** of windows (prereg threshold met) |

Artifacts:
- `data/derived/elai_genome/qc_table.csv` — per-window confidence/ARI map
- `data/derived/elai_genome/composition_table.csv` — per-window state-vs-country NMI
- `data/derived/elai_genome/windows/<tag>/w_s2_k10_seed<seed>.q.npy` — the 1,000 × 55
  unordered state-posterior matrix per window (the conditioning variable)
- interesting negatives: confident-but-not-country-predictive windows, e.g. chr3:138.5 Mb,
  chr10:74.5 Mb, chr5:170.5 Mb

**Locality evidence** (AMY region site ladder, `analysis/elai_site_ladder.py`): states are
piecewise-local — only 30–48% of individuals keep the same state between adjacent sites —
and the AMY CNV core itself is the *least* resolvable site in its region (16.8% confident
vs ~100% at flanks). Design consequence: condition on the state at **flanking markers**,
not the core. Resolution note: states are per 1 Mb window; within-window switches are
below resolution (acceptable at 0.5–1 Mb tract scale).

## 5. Not done yet (next steps, in order)

1. **Pre-register the genome-wide temporal analysis** — *before* any results:
   - conditioning variable: unordered state q at the variant's window (confidence- and
     ARI-gated per window); state defined at flanking markers where the variant sits in a
     low-confidence window
   - test: allele-frequency change over time **within** each shared state; concordance of
     that change across populations sharing the state; effective number of independent
     states as the replication ceiling (updates the N_eff ≈ 2–4 estimate)
   - detection floor: carried over from the n=225 within-history calculation (0.146)
   - decision rules fixed in advance; no threshold tuning
2. **Run it** genome-wide, then place AMY1 in the resulting distribution.
3. **Deferred**: repo-wide data-integrity audit (`docs/data-integrity-audit.md` was
   proposed, never started — motivated by the five silent defects already found: TGENO bit
   order, REF/ALT orientation, gene-span bounding, AMY exclusion coordinate, and their
   effect on LCT/FADS/SLC45A2 too).

## 6. Known limitations (carry these into any write-up)

- Agreement with Akbari reuses the same imputed genotypes — not independent validation.
- S=2 sources are *learned unsupervised* (`-p 1`), not labelled populations; states carry
  geographic structure (NMI 0.18) but are not Anatolian/WHG/Steppe calls. The S-layer
  (upper source) posterior is not exported yet — our patch exports lower pairs marginalised
  over S; a few-line fork extension would add the upper-state call if needed.
- `N_eff ≈ 4.66` (13 trajectories) / ≈ 2.0 (7 components) caps replication power.
- The 1240K capture and imputation are shared assumptions with the paper being critiqued.
- Full ordered-joint export permanently NOT TESTED at density (would be ~50 GB/fit);
  marginal + `--pkk-site` is the production mode.

## 7. Supersession chain (so nobody reads a dead verdict)

`elai-amy-long-run-results.md` → `local-ancestry-verdict.md` (route closed) →
**superseded by the CORRECTION commit `fda4ff7`** (metric invalid) →
**superseded by the REVERSAL section of `local-ancestry-verdict.md` + commit `e365fca`**
(S=1/35 Mb were the defects) → **this document** (genome-wide map at corrected
configuration). The two-channel/binarised/block-survey results in §2 were never affected.
