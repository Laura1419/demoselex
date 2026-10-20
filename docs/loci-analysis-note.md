# What we do with the loci — operational note

The 256 segmental-duplication / CNV loci, and how each one is processed into a statement
about selection versus demography. This is the "how we run it" note; the reasoning is in
`docs/methods-note.md`, the code map in `docs/pipeline-entrypoint.md`.

## 1. What a locus is, in data terms

Not a copy-number dosage. At each locus every individual carries **two alleles**, each
labelled with an integer (e.g. 1..10, or 1..20). The diplotype (a, b) *is* the copy-number /
haplotype configuration. One row per (locus, individual):

```
chrom,pos,sample,allele1,allele2
2,135851076,HG00096,1,2
2,135851076,HG00097,2,2
...
```

The allele labels come from our own genotyping, so they carry meaning independent of the
ARG (this matters: the ARG's internal allele orientation is arbitrary — see §7).

## 2. What we compute at each locus

Using the dated ARG at that locus (GRCh38), for **each distinct allele label**:

| quantity | meaning |
|---|---|
| `freq` | allele count / called haplotypes |
| `age` | tsdate time of the MRCA of that allele's carriers |
| `age_pct` | percentile of `age` among clades with the same descendant count in the same local tree (empirical, within-locus null) |
| `n_backgrounds` | how many ARG-derived local backgrounds hold ≥5% of that allele's carriers |
| `dom_share` | fraction of the allele's carriers that sit in its dominant background |
| `dom_freq` | the allele's frequency inside that dominant background |
| `background_n` | number of haplotypes in the dominant background |

A **local background** is the most recent clade, in the local tree at that position,
containing a haplotype whose root is younger than `tau` (default 200 generations). It is
derived from the genealogy alone — no population labels are used anywhere.

## 3. How to run it

```
.venv/bin/python analysis/g1k_allele_analyzer.py \
    --table diplo.csv --tau 200 --out results.csv
```

- `--table`: the diplotype CSV above; may contain many loci (`chrom`,`pos` group them).
- `--tau`: the background depth in generations. 200 is the calibrated default; breadth is
  the tau-stable quantity, so a sensitivity run at 100 and 500 is cheap and worth reporting.
- The analyzer restricts the ARG to the samples present in the table, so the genotyped and
  ARG cohorts must overlap. Samples absent from the ARG are skipped (and counted).
- `--selftest` recodes a known SNP as a two-allele diplotype; use it to confirm the install
  and to see the expected contrast.

## 4. Reading the output

Per locus, the interpretable pattern is in the **breadth** (how many backgrounds carry the
allele) crossed with **age** and **within-background frequency**:

| breadth | clade age | within-background freq | reading |
|---|---|---|---|
| 1 background | young | high (fixed in it) | candidate recent selection on that haplotype |
| 1 background | old | high | old haplotype, ancestrally confined (drift, or ancient selection) |
| many backgrounds | old | — | structure / carriage: the allele's distribution tracks ancestry turnover |
| many backgrounds | young | — | recent expansion across ancestries; check artefacts before interpreting |

Worked anchor (self-test, LCT rs4988235 recoded in EUR, tau = 200):

| allele | freq | age | age_pct | n_backgrounds | dom_freq |
|---|---|---|---|---|---|
| derived (ALT) | 0.536 | 495 gen | 20.0 | **1** | **1.00** |
| ancestral (REF) | 0.465 | 6,510 gen | 93.3 | 3 | 1.00 |

So at a locus like LCT, the derived allele is young, confined to one background and fixed
within it, while the ancestral allele is old and distributed — exactly the contrast the
method is designed to make visible.

## 5. QC per locus, before believing any row

1. **Enough called genotypes**: require a minimum number of called haplotypes; report the
   rate. Sparse loci give unstable backgrounds.
2. **Background count is meaningful**: if the locus has very few large backgrounds (say <3
   with ≥20 haplotypes), "breadth" is not informative at that locus — report it as such
   rather than as "confined".
3. **tau sensitivity**: report breadth at 100 / 200 / 500 generations; a reading that flips
   with tau is not a claim.
4. **Mapping/paralogy**: the allele labels are genotyped, so this is less of an issue than
   for SNP-based analysis, but still record each locus's mappability / duplication structure.
5. **Cohort overlap**: state how many genotyped samples were present in the ARG.

## 6. What we present

Per locus: allele label, frequency, clade age (with percentile), breadth, dominant-background
share and within-background frequency, and the number of backgrounds. Plus, in aggregate:

- the distribution of breadth across the 256 loci (how many CNV alleles are haplotype-confined
  versus ancestry-distributed),
- the joint distribution of age percentile and breadth (the two axes together separate
  regimes),
- the control rows from §4 as calibration,
- the skyline Ne(t) as the demographic frame for any absolute age quoted.

## 7. What we do not claim

- Not selection coefficients, and not dates beyond the resolution of the ARG (ages are
  substrate-dependent; compare percentiles, not absolute ages, across loci).
- Not "derived" versus "ancestral" from the ARG's own allele states: the build passed
  placeholder alleles to `tsinfer`, so its internal orientation is arbitrary. The diplotype
  labels sidestep this, but any statement about which allele is derived must come from the
  genotyping/annotation, not from the tree.
- Not a p-value from `age_pct`: it is an empirical percentile within the local tree, not a
  test statistic with a sampling distribution. Breadth and within-background frequency are
  descriptive; treat them as such.

---

## BED-level scan (available now, before the diplotype data)

`analysis/g1k_bed_scan.py loci.bed --pop EUR --breadth` summarises any BED whose 4th column
holds a gene/locus name. Per locus it reports: number of scored variants, how many fall in
the young-clade bottom 5% (k5) versus expectation, the Poisson p, the most extreme variant's
percentile and age, whether the locus is in the segdup compartment, and (with `--breadth`)
the local-background breadth of that extreme variant.

Demonstration on six gene-named loci (GRCh38, +/-20 kb, EUR):

| locus | n variants | k5/expected | enrichment | Poisson p | best variant pct | segdup | breadth |
|---|---|---|---|---|---|---|---|
| LCT | 162 | 17 / 8.1 | **2.10x** | **0.004** | 1.4 | no | 1 |
| AMY1A | 1 | 1 / 0.1 | 20x (n=1) | 0.049 | 10.9 | yes | 1 |
| HLA-B | 1,862 | 29 / 93.1 | 0.31x | 1.00 | 5.5 | no | 1 |
| CYP2D6 | 190 | 5 / 9.5 | 0.53x | 0.96 | 3.7 | yes | 1 |
| GSTM1 | 117 | 3 / 5.9 | 0.51x | 0.93 | 2.6 | yes | 1 |
| SLC24A5 | 16 | 0 / 0.8 | 0.00x | 1.00 | 73.2 | no | 3 |

Two things this shows, and they are the argument for structural genotyping:

1. The SNP-based local-ancestry metric fires where a SNP sweep exists (LCT, 2.1x, p=0.004)
   and correctly stays flat at the old-allele locus (SLC24A5).
2. It is **depleted** at exactly the loci whose biology is structural rather than SNP-level
   (CYP2D6 0.53x, GSTM1 0.51x), and near-blind at AMY because the duplication region carries
   almost no common SNPs in this panel (1 variant in AMY1A +/-20 kb). Those loci need the
   structural haplotypes - i.e. COSI-GT output - not SNP clade ages.

Note the breadth column refers to the locus's *most extreme* variant, which at a locus like
SLC24A5 is not the canonical SNP; when the diplotype data arrive, breadth should be computed
per allele label with `g1k_allele_analyzer.py`, which is the authoritative path.
