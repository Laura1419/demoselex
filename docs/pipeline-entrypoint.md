# Pipeline entry point — what to run, what is authoritative, what is legacy

Companion to `docs/methods-note.md` (the reasoning) and `docs/arg-metric-methods.md` (the
full development record). This file is the only thing you need to find your way around.

## The two authoritative statistics

Everything else is either infrastructure or a documented alternative. Claims should be made
with these two only.

1. **Local-background breadth and within-background frequency** (`g1k_local_bg.py`,
   `g1k_allele_analyzer.py`). For a locus: partition haplotypes into ARG-derived local
   backgrounds (most recent clade under tau), then report how many backgrounds carry the
   allele (>=5% of its carriers), the dominant background's share and size, and the allele's
   frequency inside it. **Label-free** — no population labels are used.
   Reading: one background + high within-background frequency = haplotype-confined state;
   many backgrounds = distributed across ancestries.

2. **Clade-age percentile, population-conditioned** (`g1k_pop.py`, `g1k_null.py`,
   `g1k_allele_analyzer.py`). The carrier clade's tsdate age, standardised against clades of
   the same descendant count in the same local recombination regime and compartment. Reported
   as a percentile; window-level enrichment (FDR) is secondary and fragile in low-rho regions.

Calibration (EUR, 1000G high-coverage, 1,266 haplotypes): LCT derived allele age 495 gen,
percentile 20, breadth 1, within-background frequency 1.00; ancestral allele 6,510 gen,
percentile 93, breadth 3. Window enrichments: LCT q=0.0008, SLC45A2 0.99x (flat), FADS 1.06x
(flat), AMY 1.39x (q=0.0017).

## Canonical scripts

| script | role |
|---|---|
| `analysis/g1k_phased.py` | build / date / stats: VCF -> tsinfer SampleData -> dated ARG per chromosome (10 Mb chunks) |
| `analysis/g1k_pop.py` | population-restricted statistics from dated chunks (simplify to a subset) |
| `analysis/g1k_null.py` | composite null (clade size x local rho x compartment), per-variant scoring, window FDR |
| `analysis/g1k_local_bg.py` | local backgrounds + breadth/confinement statistics |
| `analysis/g1k_allele_analyzer.py` | **per-locus analyzer for multiallelic diplotype data** (segdup/CNV loci) |
| `scripts/g1k_run_chr.sh` | one chromosome end to end (download -> build -> date -> stats), frees intermediates |

## The downstream analyzer (the point of the exercise)

Input: a CSV with one row per (locus, individual):

```
chrom,pos,sample,allele1,allele2
```

Each individual carries two alleles, each an integer label (1..10, 1..20, ...); the diplotype
is the copy-number/haplotype configuration. No dosage encoding is involved.

```
.venv/bin/python analysis/g1k_allele_analyzer.py --table diplo.csv --tau 200 --out results.csv
```

Output per (locus, allele): `n`, `freq`, `age`, `age_pct`, `n_backgrounds`, `dom_share`,
`dom_freq`, `background_n`. Plus `--selftest`, which recodes a known SNP as a 2-allele
diplotype and reproduces the expected contrast (see above).

Caveat that must be respected: the analyzer uses the ARG's local tree, and the ARG's allele
orientation is arbitrary (see below), so it is valid for *allele-label-based* questions —
where the labels come from your own genotyping, not from the ARG — and NOT for assigning
"derived" status. Since diplotype labels are external, this caveat does not apply to the
downstream application.

## Data layout

| path | content |
|---|---|
| `data/derived/g1k/chr<N>_<start>_<end>.dated.trees` | dated ARG chunks (GRCh38) |
| `data/derived/g1k/chr<N>_stats.csv` | per-variant clade age + rho (all 6,404 haplotypes) |
| `data/derived/g1k/chr<N>_EUR_stats.csv` | same, EUR-restricted (1,266 haplotypes) |
| `data/derived/g1k/composite_null.npz` | composite null cells (EUR) |
| `data/derived/g1k/composite_score_chr<N>.csv` | per-variant p_age (EUR) |
| `data/derived/g1k/composite_regional.csv` | window enrichments + FDR |
| `data/derived/g1k/samples_population.txt` | 3,202 samples, population and super-population labels |
| `docs/methods-note.md` | method, controls, limits, application plan |

## Legacy / not for claims

`arg_build.py`, `arg_merge.py`, `arg_modern.py`, `arg_traj.py`, `arg_cladenull.py`,
`arg_composite_null.py`, `arg_phased.py`, `arg_phased_null.py`, `g1k_null`'s predecessors —
superseded by the g1k pipeline above; kept for the development record.

`g1k_omega.py` — **do not use.** Mis-specified (see `arg-metric-methods.md`).

`g1k_omega2.py` — empirical branch-length statistic; mechanically fixed but not calibrated
(expectation ignores local Ne). Do not use for claims.

`g1k_relate_test.py` — Relate's actual selection test (n0 carriers given k0 lineages at
origin, analytic null). Correct per-variant on 3 of 4 controls, but (i) inherits the ARG's
arbitrary allele orientation and (ii) is not calibrated across windows. A documented
alternative, not the primary instrument.

## Known defect and its fix

**Allele orientation.** The build passed placeholder alleles (`['A','T']`) to tsinfer, so the
"derived" state is whichever allele is the global minority. For an allele near-fixed in the
analysed population this measures the complement lineage. Fix: carry REF/ALT from the VCF
into the SampleData at build time (requires re-streaming the VCFs). Until then: do not make
derived-vs-ancestral claims from the ARG's own states; external allele labels (as in the
diplotype analyzer) are unaffected.
