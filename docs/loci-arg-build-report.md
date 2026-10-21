# Loci-only ARG pipeline — build report

*Branch `loci-arg-pipeline`, commit of 2026-10-21. This document explains what was
built, why, and how it was validated. Companion code: `analysis/g1k_loci_build.py`,
`analysis/build_input_tables.py`, `scripts/download_g1k_vcf.sh`,
`scripts/slurm_loci_arg.sh`, and the modified `analysis/g1k_allele_analyzer.py`.*

## 1. Context and goal

The diplotype analyzer (`analysis/g1k_allele_analyzer.py`) needs the dated 1000G ARG
(`data/derived/g1k/chr<N>_<start>_<end>.dated.trees`) to turn allele labels into
statements about selection vs demography. That ARG was previously built whole-genome on
another machine and is not in the repo. We rebuild it here, but only for the loci we
actually analyze: the 265 SV loci of `svs_refined.hgsvcv3.bed` (GRCh38) touch
**167 of ~2,885 ten-Mb chunks (~6% of the genome)** — so a loci-only build saves 94%
of the compute, at the cost of reduced inference context (see §3).

## 2. What was built

### 2.1 Per-locus dated chunks (`analysis/g1k_loci_build.py`)

For each 10 Mb grid chunk that a locus touches (plus the LCT validation chunk
chr2:135–136 Mb, included in every plan):

1. **stream** the 1000G high-coverage phased panel VCF (3,202 samples, GRCh38) for a
   window = chunk ± 5 Mb flank, keeping the exact site filters of `g1k_phased.py`
   (SNVs only, panel AF with min(AF,1−AF) ≥ 0.05, duplicate positions skipped);
2. **tsinfer** the window (flank = inference context);
3. carve the **central 10 Mb**, `keep_intervals().trim().simplify()` — precisely the
   operations `g1k_phased.date()` applies to a full-chromosome ARG;
4. **tsdate** `inside_outside`, µ=1.25e-8, Ne=10,000, generations — identical parameters;
5. dump under the canonical name `chr<N>_<start>_<end>.dated.trees` with positions
   relative to the chunk start, so `chunk_for()` in the analyzer works unchanged.

Engineering details that mattered in practice:

- **Resumability**: a finalised `SampleData` for a chunk is reused; a partial one is
  re-streamed. Completed chunks are skipped by existence.
- **htslib error detection**: stderr is captured during streaming; any `[E::` line
  (e.g. `bgzf_uncompress` failures) deletes the SampleData and fails the task so the
  scheduler retries, instead of silently producing a truncated ARG. This was added
  after remote streaming from EBI produced ~30% corrupted windows under 72-way
  concurrency (some "completing" with a fraction of their sites).
- **Local VCF preference**: if `data/derived/g1k/vcf/chrN…vcf.gz` exists it is used
  instead of the HTTPS URL. All 22 chromosomes were downloaded once (27 GB,
  size-verified, segmented + resumable — `scripts/download_g1k_vcf.sh`) after which
  streaming is fast and reliable.

**Faithfulness caveat**: the main build gives tsinfer whole-chromosome context; this
build gives it 20 Mb. The LCT selftest (§4) shows the difference is small at the
statistics we read.

### 2.2 Environment

A conda env at `.venv` with a self-built htslib 1.24 (the bioconda build lacks libcurl;
the pip cyvcf2 wheel bundles OpenSSL 1.0.2k from 2017 which cannot verify modern
certificates). This is what enables remote tabix access at all; the full-chromosome
downloads then made it moot for production runs.

### 2.3 The two input tables (`analysis/build_input_tables.py`)

Builds the INPUT CONTRACT files from `docs/loci-analysis-note.md`:

**`alleles_clusters.tsv`** — one row per (locus, sample), full grid, `NA` for uncalled:

| column | source |
|---|---|
| `chrom, pos, start, end` | `svs_refined.hgsvcv3.bed` (pos = interval midpoint; start/end = locus span) |
| moderns (`HG/NA…`, all 3,202) | **read-based genotypes** `geno_3202_fa/geno_3202_master.tsv.gz` (likegt, rank-1 calls) |
| ancients (1,628 samples) | `aDNA_real_1x_geno/results/*.geno1.tsv` (likegt on aDNA) |
| `allele1, allele2` | **cluster labels** from the MDS/DBSCAN haplotype clustering (`cosi_like/mds/paths_<locus>.tsv`) |
| `path1, path2` | original path keys (`sample#copy#contig`, the `:start-end` graph interval stripped — it is not a stable identity) |
| `cluster_k, silhouette, noise_frac` | per-locus clustering QC from `mds_summary.tsv` |

Label scheme: per locus, DBSCAN clusters are numbered **by descending size (1 = the
common class)**; **noise paths get private labels after k** so nothing is dropped
(see rationale below). `allele_labels_clusters.tsv` maps label → path key → cluster id.

Why cluster labels (this is the important design decision):

- With raw path-key labels, every modern's paths are private (their own assemblies):
  each label has exactly one ARG-cohort carrier, and the analyzer's `n ≥ 3` filter
  drops everything. The analysis is vacuous.
- Cluster labels define **structural haplotype classes** with many carriers, shared
  between moderns and ancients through the paths they were genotyped onto — which is
  the operational form of the "allelic haplotype sharing is the observable" argument
  in the methods note.
- Noise→NA was rejected with data: it would delete 15.2% of allele instances
  genome-wide and **78–86% at the most interesting loci** (LPA 78%, KATNAL2 86%,
  KIR2DL1 80%) — the hyper-diverse loci where DBSCAN's noise is likely real biology.
  Private labels keep those calls; they simply sit out of class-level statistics.

**`dates.tsv`** — `sample, date_bp, date_sd, source, culture` from
`old_dups_metadata-V1_dated.tsv` (date_mean_BP, date_sd_BP, date_method,
cultural_genetic_context); moderns at date 0. 5,162 dated ancients + moderns.

Known caveats (also printed by the builder): the GYP cluster has no modern graph-truth
rows when the identity table is used (per-gene CN cannot reconstruct the merged-locus
diplotype) — moot with the read-based 3202 genotyping, which genotypes the merged
locus; `--min-qv` is available but defaults to 0; the anchor `pos` is the SV-interval
midpoint (see §2.4).

### 2.4 Analyzer upgrades (`analysis/g1k_allele_analyzer.py`)

- **Ingestion hardened**: separator sniffed (TSV or CSV), `NA`/missing alleles →
  uncalled (`-1`), and `freq` divides by **called** haplotypes.
- **Multi-anchor (`--anchors K`)**: a structural class spans many local trees across a
  megabase-scale locus, and no single midpoint tree represents it. With `--anchors K`
  the analyzer samples K trees evenly across the locus span (using the table's
  `start`/`end`, clamped into the chunk) and aggregates per allele: median `age` plus
  `age_min`/`age_max`, median `age_pct`, `n_backgrounds`, `dom_share`, `dom_freq`,
  `background_n`, and `n_anchors` (how many anchors yielded `n ≥ 3`).
  Without `--anchors` (or without start/end columns) behaviour is the previous
  single-anchor mode.

## 3. Infrastructure

Chunks are independent single-threaded jobs: `scripts/slurm_loci_arg.sh` runs them as a
SLURM array on the `workers` partition (72 concurrent, 12 GB/task), reading chunk
coordinates from `data/derived/g1k/slurm_chunks.txt`. Per-chunk wall time ≈ 15 min
stream (local VCF on MooseFS) + 1–12 min tsinfer + 14–23 min tsdate ≈ 45–55 min;
159 chunks ≈ 2.5–3 h on 6 nodes.

## 4. Validation

**LCT selftest** (windowed rebuild vs the documented full-chromosome calibration,
EUR, tau=200):

| quantity | documented | this build |
|---|---|---|
| derived freq | 0.536 | 0.5355 |
| derived age | 495 gen | 516.5 |
| derived backgrounds / dom_freq | 1 / 1.00 | **1 / 1.00** |
| ancestral age / pct / backgrounds | 6,510 / 93.3 / 3 | 6,402.6 / 92.9 / **3** |

**chr22 pilot** (4 loci, all 3,202 moderns EUR-restricted to 633, `--anchors 9`):

- The full chain ran end to end: read genotypes → cluster classes → local ARG →
  per-class freq/age/backgrounds.
- Multi-anchor exposed the recombination-scrambling it was designed for: per-class
  ages vary hugely across a locus's trees (e.g. CYP2D7 allele 1: median 8,814 gen,
  range 4,203–34,078), where the single-anchor mode had returned one identical MRCA
  for every class at a locus.
- Read-based moderns reshape the spectrum vs the 233 graph-truth samples: POTEH
  changes from an apparent 96% single class to 52%/25%/5% — the graph-panel samples
  flatten structure.

## 5. Status and next steps

- 9 chunks built and validated at time of writing; the remaining 159 building on
  SLURM (ETA ~2.5–3 h from array start).
- Next: genome-wide analyzer run (256 loci; EUR and global cohorts); then the
  temporal layer (`--dates`: time-binned within-background frequencies, drift null) —
  the ancients (1,628 genotyped, 5,162 dated in the metadata) enter there.
- Interpretation reminders from `docs/loci-analysis-note.md` still apply: age first,
  breadth only when age > tau; no selection coefficients; `age_pct` is descriptive.
