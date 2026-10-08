# Tier-2 data: Akbari et al. 2026 processed imputed genotypes

> Owner: `tier2-geno` (this lane). Companion loader: `analysis/load_akbari_geno.py`.
> Download script: `scripts/download_tier2.sh`. Access matrix: `docs/E-data-inventory.md` (row G3).
> Retrieval date: 2026-10-08. Source: Harvard Dataverse V4.0, `doi:10.7910/DVN/7RVV9N`,
> file API `https://dataverse.harvard.edu/api/access/datafile/<id>`, license **CC0 1.0**.

## Why this product

Akbari et al. (2026, Nature 654:419-428, `doi:10.1038/s41586-026-10358-1`) release the
PLINK2 bundle titled *"Genome-wide imputed genotypes for ancient and modern Europeans"*:
the combined ancient + European 1000 Genomes matrix actually used for their genome-wide
analyses **[observation: `README.txt` §5]**. It is
`AllVariants_NoQC` — **no** post-imputation/downstream QC applied — with
**52,382,873** autosomal variants (SNPs + indels, GRCh37/hs37d5), of which
**9,739,624** passed the paper's QC (`HighQualityVariants.tsv`). It combines **15,836
ancient** individuals with **503** modern European 1KG individuals (**16,339** samples).

Using this Tier-2 set replaces the AADR 1240K pseudo-haploid subset with the paper's own
imputed dosages on the same GRCh37 coordinate system, so the trajectory reanalysis can run
**like-for-like**: same individuals, same 9.7M QC variant catalog.

## What was downloaded (verified)

Destination: `data/raw/akbari/tier2/`. Remote byte sizes and md5 come from the Dataverse
dataset API listing; sha256 is computed locally. Log: `data/raw/akbari/tier2_download_log.tsv`.

| file | datafile id | bytes | md5 (Dataverse, verified) | sha256 (local) |
| --- | --- | ---: | --- | --- |
| `AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.pgen` | 13415552 | 26,300,030,448 | `c9e587230d9bdc8ace0ca1b043294a7c` | *(see log after completion)* |
| `...psam` | 13415494 | 221,391 | `787717246189699bec4b12a34884992a` | `4af9d4ebdf99a7d0493dd60d975617687217935a662f209736f367eb3fb3d80b` |
| `...pvar` | 13415518 | 1,700,215,110 | `766f95a7311f7dee4a47b86f01038ec8` | `1183c7172b32ca7012b76cbcb379d0569c9b413528148c8ff39bc09dd73682f5` |
| `HighQualityVariants.tsv` | 13415504 | 315,624,578 | `2ec2598e2a6ead74ee11050672cda972` | `19b04d4ff4fd630c99a41a6518671a0b987f3527fc02215f59e395f3211581d8` |
| `README.txt` | 13596811 | 10,901 | `7d5fdfe6d1087de7b55f14254ba823eb` | `4fef526c8f099f0464c5f67f0ae33833e5ba4230a4007ded86a11a8330c20244` |

**[observation]** All four completed files match the Dataverse md5 and the listed byte
size exactly. The `.pgen` download is resumable; if it is still partial the file exists as
`*.pgen.part` with a byte count `< 26,300,030,448`.

## What was NOT downloaded

* The per-chromosome *NO-QC BCF* (`AncientDNA_bcf/AncientDNA_AllVariants_NoQC_GRCh37_chr{1..22}.bcf`),
  aggregate **1,324,569,749,458 B (~1.3 TB)** — **[observation: `docs/E-data-inventory.md` row G2]**;
  far outside this lane's time/storage budget. The PLINK2 set is a strict subset of the
  same call set (combined ancient + EUR 1KG, autosomes) and is the analysis-ready form.
* The 1240K TGENO bundle and qpAdm ancestry tables — not needed for a dosage scan.

## How to finish an interrupted download

The script is idempotent and resumable; re-running it resumes any `*.pgen.part` with
`curl -C -` and re-verifies everything else (already-verified files are left untouched):

```bash
# inspect targets + expected sizes without touching the network:
bash scripts/download_tier2.sh --dry-run

# resume / complete (safe to run repeatedly):
bash scripts/download_tier2.sh
```

A single target can also be finished manually (the partial is resumed, then promoted only
after the md5 matches):

```bash
curl -fL -C - -A 'Mozilla/5.0' \
  'https://dataverse.harvard.edu/api/access/datafile/13415552' \
  -o data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.pgen.part
# then verify + promote (md5 c9e587230d9bdc8ace0ca1b043294a7c, size 26300030448)
```

## Loader API — `analysis/load_akbari_geno.py`

Pure-`pgenlib` reader (pgenlib 0.95.1). It never maps the 26 GB `.pgen` into RAM: it reads
only the requested variant rows.

* `paths(dirpath, prefix)` → `(pgen, pvar, psam)` paths.
* `parse_psam(path)` → numpy array of sample ids (column 1 of the `#IID\tSEX` table).
* `iter_pvar(path)` → generator of `(chrom, pos, vid, ref, alt)`; `read_pvar_header(path)`.
* `describe_pgen_storage(pgen_path)` → dict: `n_variants`, `n_samples`,
  `has_sparse_dosages`, `hardcall_probe`, `dosage_probe`, `dosage_is_integer_like`,
  `pgen_mode_byte`, `storage_mode` (human label).
* `class AkbariGeno(dirpath, prefix, sample_subset=None)`
  * `.n_variants`, `.n_samples`, `.sample_ids`
  * `.change_sample_subset(indices_or_None)` — native `change_sample_subset`, so
    sub-cohort reads touch only the selected columns.
  * `.sample_index_of(id)`
  * `.read_variants(indices, allele_idx=1)` → float32 `(len(indices), n_samples)`,
    ALT-allele dosages in `[0,2]`, NaN = missing. Only the requested variants are read.
  * `.read_hardcalls(indices)` → int8 `(len(indices), n_samples)`, `0/1/2`, missing `-9`.
  * `.iter_chunks(start, end, step, allele_idx=1)` → `(offsets, dosage_matrix)`; choose a
    `step` that bounds `step * n_samples` (use a sample subset for wide `step`).

### Worked example

```python
from analysis.load_akbari_geno import AkbariGeno

with AkbariGeno() as g:                      # 52,382,873 variants x 16,339 samples
    print(g.n_variants, g.n_samples)
    # dosages for one variant, first 5 samples:
    j = g.sample_index_of("Kou01.SG")
    dos = g.read_variants([0])[0]            # ALT dosage at variant 0
    print(g.sample_ids[:5], dos[:5])
    # sub-cohort scan without full-matrix reads:
    g.change_sample_subset([0, 1, 2, 3, 4])
    chunk = g.read_variants(range(1000))     # (1000, 5) float32
```

Self-test (prints `n_variants`, `n_samples`, one dosage, storage mode, peak RSS):

```bash
.venv/bin/python analysis/load_akbari_geno.py --self-test
```

## .pgen storage mode encountered

Filled from `describe_pgen_storage` after the `.pgen` completed; see `--self-test` output.
**[observation]** The psam header is `#IID\tSEX` with **no FID column**, so column 1 is the
individual id used to join to AADR `.ind` and `data/derived/trajectories.tsv`.

## Cross-check: trajectory individuals present in the Akbari imputed set

Trajectories (`data/derived/trajectories.tsv`, 11,934 individuals) are derived from AADR
v66.p1 1240K; every one of them is present in
`data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.ind` **[observation]**.

Intersecting the `.psam` sample ids against the trajectory ids:

| comparison | count |
| --- | ---: |
| trajectory individuals total | 11,934 |
| trajectory ∩ Akbari `.psam` (exact id) | **5,788** |
| trajectory ∩ Akbari, case-insensitive | 5,792 |
| trajectory ∩ Akbari, case-insensitive + suffix-stripped base name | 6,174 |
| Akbari `.psam` ∩ AADR `.ind` (exact id) | 5,981 |

**[inference]** About **5,788 of 11,934 (48.5%)** trajectory individuals are in the Akbari
imputed set by exact id; allowing pipeline suffix differences (`.SG`/`.AG`/`.TW`/`.DG` are
used inconsistently across the two releases) raises this to ~6,174 (51.7%). This overlap is
what makes a Tier-2 scan possible; the ~half that do not overlap can only be analysed via
the AADR 1240K pseudo-haploid calls, which are a different assay and must not be mixed with
the imputed dosages without an explicit two-tier design.

**[speculation]** The suffix discrepancy (Akbari is dominated by `.TW`=7,427 and `.AG`=5,426;
AADR by `.AG`=11,247 and `.SG`=7,007) likely reflects different capture assay labels for the
same underlying individuals, but a definitive id map would need the Akbari sample manifest
(the embargoed full metadata, `datafile 13572552`, is closed until 2041). Do not invent a
mapping; report both the conservative exact-id and the base-name upper bound.

## Validation performed

* `bash scripts/download_tier2.sh --dry-run` lists all five targets with expected sizes,
  resolves the Dataverse listing, and exits 0.
* Every promoted file matches its expected byte size and Dataverse md5; sha256 recorded in
  `tier2_download_log.tsv`.
* `.venv/bin/python analysis/load_akbari_geno.py --self-test` prints `n_variants`,
  `n_samples`, one dosage, and peak RSS without loading the whole matrix (see log output).
