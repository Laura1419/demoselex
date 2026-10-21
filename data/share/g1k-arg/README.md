# Prerequisites: run the allele analyzer on the 1000G ARG

Everything here is what `g1k_allele_analyzer.py` needs. It takes **one input file** from you
(the allele table) and writes a per-locus, per-allele summary.

## What is in this directory

| file | what it is |
|---|---|
| `arg/chr<N>_<start>_<end>.dated.trees` | dated ancestral recombination graphs (GRCh38), one set of 10 Mb chunks per chromosome, built from the 3,202-sample 1000G high-coverage release. **This is the big part (~58 GiB, 294 files).** |
| `samples_population.txt` | the 3,202 samples **in ARG haplotype order**. Haplotype indices are `2i` and `2i+1` for the i-th row (0-based). Get this ordering wrong and every result is wrong. |
| `g1k_allele_analyzer.py` | the analyzer |
| `sv265_annotation.bed` | the 265 COSI-GT / HGSVC v3 loci with gene names and DEL/INS/COMPLEX class (optional, for reporting) |
| `genomicSuperDups_hg38.txt.gz` | segdup mask, GRCh38 (optional; used only to label loci) |
| `loci-analysis-note.md` | method, per-locus quantities, QC, what we do and do not claim |

## Input you need to produce (the only thing we need from you)

One row per sample per locus, diplotype as two integer allele labels:

```
chrom	pos	sample	allele1	allele2
1	103304997	HG00096	1	3
1	103304997	HG00097	2	2
```

- `allele1`,`allele2`: integer labels, consistent within a site. **Unordered** is fine
  (structural haplotypes are unphased). Missing = `NA` or `-1`.
- `pos` must be **GRCh38**. `chr` prefix optional (`1` and `chr1` both work).
- `sample` must be a 1000G sample ID present in `samples_population.txt`. Samples not in that
  file are skipped (and counted), so check the overlap your run reports.

## Setup (once)

```bash
python3 -m venv .venv
.venv/bin/pip install tskit numpy pandas scipy
```

## Run

```bash
.venv/bin/python g1k_allele_analyzer.py --table my_alleles.tsv --tau 200 --out my_results.csv
```

Useful variants:

```bash
# sanity check the install: recodes a known SNP as a 2-allele diplotype and prints the
# expected contrast (derived allele: young clade, one background, fixed within it)
.venv/bin/python g1k_allele_analyzer.py --selftest

# sensitivity to the background-depth dial (breadth is the tau-stable quantity)
for t in 100 200 500; do
  .venv/bin/python g1k_allele_analyzer.py --table my_alleles.tsv --tau $t --out res_$t.csv
done
```

## Output, per (locus, allele)

| column | meaning |
|---|---|
| `n` | carriers (called haplotypes) |
| `freq` | allele frequency among called haplotypes |
| `age` | tsdate age of the carriers' most recent common ancestor (generations) |
| `age_pct` | percentile of that age among clades with the same descendant count in the same local tree (empirical within-locus null) |
| `n_backgrounds` | ancestry backgrounds holding >=5% of carriers |
| `dom_share`, `dom_freq` | share of carriers in the dominant background; allele frequency inside it |
| `background_n` | haplotypes in the dominant background |

## Reading it — and one hard rule

- **Young allele** (`age` small, `age_pct` small): read `age` / `age_pct`. Recent spread on a
  haplotype background.
- **Old allele** (`age` greater than tau): read `n_backgrounds` / `dom_freq`. One background
  and high within-background frequency = confined state; many backgrounds = distributed
  across ancestries (structure / carriage).
- **Hard rule:** `n_backgrounds` is **vacuous when the carrier clade is younger than tau** —
  all carriers then merge below tau and terminate at the same node, giving breadth = 1 by
  construction. Breadth is an *old-allele* statistic. Always read `age` first.

## Notes / caveats

- The ARG's internal allele orientation is arbitrary (built with placeholder alleles), so the
  analyzer says nothing about derived-vs-ancestral. Your diplotype labels are external and
  unaffected by this; just don't ask the ARG which allele is ancestral.
- `age` is substrate-dependent (absolute ages shift with cohort/site density); compare
  percentiles across loci, not absolute ages.
- The analyzer restricts the ARG to samples present in your table, so the background and the
  null are computed within your cohort.
