
## NEW AUTHORITATIVE CORRECTED FULL-PANEL EXPERIMENT — 2026-10-09T17:52:19.621129+00:00

[observation] Original registration persisted 2026-10-09T17:32:36.014282+00:00 before initial run 2026-10-09T17:39:57.291589+00:00. Original SHA256 ca49cf2edc71b4d85554a49cb8545c921766d167e62bafdf1063385b27f9772e matches retained `prereg.original.md`. The clearly marked supervisor-approved POST-RESULTS coordinate correction persisted 2026-10-09T17:46:23.590224+00:00 before CORRECTED run 2026-10-09T17:50:10.081408+00:00. Amended prereg SHA256 fef1a531027ab40cb4b0e31199556857f29af7cd23f0c3c1444b892a9c16f73f. No K, encoding, width parameter or gate threshold was retuned. Initial outputs are retained under `invalid_exclusion_run/`; **their biological exclusion gate is INVALID**, irrespective of their numerical zero intersection. This section does not rewrite previous conclusions.

### Biological exclusion repair and generalised runtime annotation assertions

[observation] Immutable hg19 refGene identifies AMY1 copies through104301314, beyond the original104255000 exclusion. Generalised derivation includes EVERY transcript of ALL core-overlapping genes, rounds outward to1kb and applies the SAME25kb buffer to CORRECTED core ends. All focal/core transcripts are now covered; runtime assertions enforce this and every feature/window interval is checked against these annotation-derived bounds. Original control exclusions also miss full focal transcripts at LCT/FADS/SLC45A2; SLC24A5 originally covers its focal gene but the generalised core additionally covers neighbouring CTXN2/MYEF2.

| locus | original_exclusion | corrected_core | corrected_exclusion | old_uncovered_focal_transcripts | corrected_uncovered_core_transcripts |
| --- | --- | --- | --- | --- | --- |
| AMY | [103915000, 104255000] | [103940000, 104302000] | [103915000, 104327000] | 5 | 0 |
| LCT | [136583646, 136633646] | [136545000, 136635000] | [136520000, 136660000] | 1 | 0 |
| FADS | [61544830, 61594830] | [61567000, 61660000] | [61542000, 61685000] | 4 | 0 |
| SLC45A2 | [33926693, 33976693] | [33944000, 33985000] | [33919000, 34010000] | 3 | 0 |
| SLC24A5 | [48401484, 48451484] | [48413000, 48496000] | [48388000, 48521000] | 0 | 0 |

[observation] ALL genes in the expanded AMY PRIMARY outer window (103415000–104827000); genes not fully excluded are explicit allowed FLANK background genes, not silently treated as covered:

| gene | start | end | transcripts | core_gene | fully_excluded |
| --- | --- | --- | --- | --- | --- |
| ACTG1P4 | 104112026 | 104114008 | 1 | True | True |
| AMY1A | 104198141 | 104301311 | 6 | True | True |
| AMY1B | 104198325 | 104301310 | 3 | True | True |
| AMY1C | 104198138 | 104301314 | 6 | True | True |
| AMY2A | 104159954 | 104168402 | 1 | True | True |
| AMY2B | 104097266 | 104122156 | 1 | True | True |
| COL11A1 | 103342023 | 103574078 | 5 | False | False |
| LOC100129138 | 104615645 | 104619693 | 1 | False | False |
| LOC101928436 | 103959602 | 104068130 | 3 | True | True |
| RNPC3 | 104068321 | 104097859 | 1 | True | True |

[observation] Equivalent all-gene and per-transcript inventories for BOTH windows of ALL five loci are saved as `*_window_genes.tsv`, `*_window_gene_transcripts.tsv`, and `*_core_transcripts.tsv`. AMY core includes ACTG1P4, AMY1A/1B/1C, AMY2A/2B, LOC101928436 and RNPC3; zero core transcripts escape the exclusion, including the3-prime end of the last AMY1 copy. Annotation source is the existing immutable hg19_refGene.txt.gz, hashed in provenance. Canonical outer windows keep their original SNP centres and200kb/1Mb spans; AMY uses the unchanged200kb/1Mb total flanks on the corrected excluded region (outer spans612kb/1.412Mb).

### Initial INVALID versus corrected authoritative result delta

[observation] Changes are material to individual assignments, not a negligible repair; ARI quantifies the delta directly. K follows the SAME predeclared selection rule, not forced to match the contaminated run:

| locus | old_K | new_K | old_new_assignment_ARI | old_features | new_features | old_size_ARI | new_size_ARI | old_group_NMI | new_group_NMI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AMY | 4 | 4 | 0.137491 | 11266 | 12133 | 0.0358703 | 0.0503559 | 0.0881703 | 0.06155 |
| LCT | 2 | 4 | 0.27971 | 10691 | 9653 | 0.0302381 | 0.0562828 | 0.0549102 | 0.0959535 |
| FADS | 2 | 2 | 0.0553207 | 12412 | 11097 | -7.87963e-05 | 0.000108175 | 0.0466704 | 0.0474203 |
| SLC45A2 | 2 | 2 | -0.00278996 | 11397 | 10943 | 0.0310407 | -0.00197335 | 0.0418045 | 0.0247979 |
| SLC24A5 | 2 | 2 | 0.883757 | 10301 | 9690 | 0.000801777 | 0.0158002 | 0.0412563 | 0.0416563 |

[inference] The initial backgrounds demonstrably contained focal-gene sequence and were invalid. Material assignment changes after removing it are consistent with contamination sensitivity, although expanding exclusions also removes additional LD markers and (for AMY) moves outer flank boundaries; the delta does not isolate a causal contribution of that one gene copy. The corrected pass/fail verdict remains the criterion below; initial statistics are not biological evidence.

### Entire panel, two distinct encodings

[observation] **ALL 16,339 psam individuals** (15,836 ancient +503 modern EUR 1KG) participate in every locus/window/K fit; no curated/reference subset and no ancestry, relative or date filter. Pgen and streamed pvar both contain **52,382,872 variants**, agreeing with the operator. The reader docstring says 52,382,873: disagreement=-1 with that docstring, not with the requested number. 32,678 nominal chromosome slots are represented by 16,339 unphased dosage rows; **32,678 phased haplotypes were NOT observed or inferred**. Genotype windows include every eligible polymorphic biallelic SNP outside the exclusions, not a hand-picked AIM panel. Indels, monomorphic sites and excluded coordinates are counted in full window audits but are not background features.

[observation] Structural encoding independently includes all 4580 supplied diplotypes/9160 H1/H2 slots and 11 labels (including non-European archive participants). 728 rows/1456 slots overlap the full imputed panel:503 modern +225 dated ancients. The remaining15,611 panel rows lack measured structural labels; they still participate in genotype clustering, but no structural label is fabricated. Structural availability therefore has NOT increased the225-person dated outcome sample.

[observation] Whole-archive structural label frequencies and full-bundle overlap (unconditional INPUT inventory, not cluster-conditioned AMY results):

| label | slots | frequency | bundle_overlap_slots |
| --- | --- | --- | --- |
| H1^a | 949 | 0.103603 | 215 |
| H2A0 | 302 | 0.0329694 | 110 |
| H2A2B2 | 141 | 0.015393 | 23 |
| H3A2 | 7 | 0.000764192 | 1 |
| H3A3B3 | 81 | 0.00884279 | 2 |
| H3^r | 4127 | 0.450546 | 626 |
| H4A2 | 228 | 0.0248908 | 12 |
| H4A2B2 | 278 | 0.0303493 | 75 |
| H5 | 2414 | 0.263537 | 259 |
| H7 | 529 | 0.0577511 | 129 |
| H9 | 104 | 0.0113537 | 4 |

[observation] All type/source/geography/component partitions are saved in `structural_partition_frequencies.tsv`. The separate published world-wide cluster-frequency table has6204 slots (ten listed classes), whereas this entire supplementary H1/H2 table has9160 slots/11 labels; these are different denominators and are not silently combined.

### Detection floor — BEFORE any outcome interpretation

[observation] Outcome-blind design calculation on225 archive-dated ancient carriers, using ONLY published baseline B frequency assumptions, the PRIMARY AMY assignments and dates:

| variant | n | frequency_assumption | time_information_kyr2 | age_span_kyr | identifiable | MDE80_ideal_per_kyr | MDE80_ideal_2kyr |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B_A_le1_published_p_0.842222 | 225 | 0.842222 | 1252.4 | 9.9565 | True | 0.0204088 | 0.0408176 |
| B_A_le2_published_p_0.742222 | 225 | 0.742222 | 1252.4 | 9.9565 | True | 0.024489 | 0.048978 |

[inference] These reported frequency/kyr floors are **optimistic independent-binomial DESIGN floors on FAILED backgrounds**, not validated HC3 floors or achieved power. Nearly all date information is retained because these unstable clusters do little validated historical conditioning. No AMY outcome regression was fitted, so a background-conditioned HC3 floor is unavailable; presenting the optimistic decrease from the previous .146 floor as a power gain would be misleading. Structural-label-specific empirical frequencies/floors and conditional slopes are deliberately not tested after gate failure.

### Gate (a): published structure / batch composition — FAILED

[observation] Published AADR group/country metadata map6490/16,339=39.72% of panel rows; all other individuals remain assigned, with explicit UNMAPPED composition. All clusters have group, country, publication, sequencing data type, suffix technology, BP range and latitude/longitude summaries in `*_composition_*.tsv`. Weighted within-group technology NMI uses1347 persons in31 groups with>=20 persons and both technologies. Suffix and publication are batch PROXIES, not verified sequencing run IDs.

| locus | small_K | primary_K | primary_silhouette | group_NMI | country_NMI | batch_proxy_NMI | gate_a |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AMY | 2 | 4 | 0.154051 | 0.06155 | 0.0107663 | 0.0556285 | False |
| LCT | 2 | 4 | 0.0746498 | 0.0959535 | 0.0305224 | 0.039669 | False |
| FADS | 4 | 2 | 0.0397222 | 0.0474203 | 0.00906635 | 0.0229734 | False |
| SLC45A2 | 2 | 2 | 0.369565 | 0.0247979 | 0.00641392 | 0.0232825 | False |
| SLC24A5 | 2 | 2 | 0.0886722 | 0.0416563 | 0.00799267 | 0.0271318 | False |

[inference] Every primary silhouette is below.10, every group NMI below.10, every country NMI below.05: known-structure gate FAILS. Batch-proxy conditional NMI passes<=.20 but cannot rescue weak known structure or unmeasured true batch. Global labels are used only for validation, never training or panel selection.

[observation] AMY primary cluster composition (full group/country/date cross-tabs archived):

| cluster | n | mapped | dated_n | min_BP | median_BP | max_BP | median_latitude | median_longitude |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1111 | 439 | 407 | 203 | 2550 | 11321 | 48.4453 | 15.2217 |
| 1 | 13325 | 5269 | 4868 | 123 | 2688.5 | 14021 | 48.1516 | 14.5 |
| 2 | 580 | 211 | 189 | 150 | 2250 | 8700 | 47.4133 | 13.65 |
| 3 | 1323 | 571 | 523 | 225 | 2844 | 10499 | 47.5138 | 16.4695 |

### K rule and gate (b): stability — FAILED

[observation] All K scores; maximal silhouette with .01 smaller-K tie rule, not outcome selection:

| locus | window | K2 | K4 | K6 | K8 | chosen_K |
| --- | --- | --- | --- | --- | --- | --- |
| AMY | 200000 | 0.212907 | 0.177915 | 0.0988791 | 0.122025 | 2 |
| AMY | 1000000 | 0.0720626 | 0.154051 | 0.0645054 | 0.0865164 | 4 |
| LCT | 200000 | 0.288979 | 0.274955 | 0.266547 | 0.242395 | 2 |
| LCT | 1000000 | 0.0607507 | 0.0746498 | 0.0662243 | 0.0550708 | 4 |
| FADS | 200000 | 0.08465 | 0.108028 | 0.102695 | 0.11553 | 4 |
| FADS | 1000000 | 0.0397222 | 0.0377496 | 0.033396 | 0.0326454 | 2 |
| SLC45A2 | 200000 | 0.331705 | 0.134085 | 0.10115 | 0.117899 | 2 |
| SLC45A2 | 1000000 | 0.369565 | 0.0554352 | 0.0234792 | 0.0420496 | 2 |
| SLC24A5 | 200000 | 0.289933 | 0.246828 | 0.185659 | 0.208793 | 2 |
| SLC24A5 | 1000000 | 0.0886722 | 0.0355027 | 0.0276249 | 0.0406526 | 2 |

[observation] ARIs (required same-K size>=.65, each adjacent K>=.60, second seed>=.80):

| locus | same_K_size_ARI | selected_K_size_ARI | neighbor_K_ARI | second_seed_ARI | exclusion_off_ARI | gate_b |
| --- | --- | --- | --- | --- | --- | --- |
| AMY | 0.0503559 | -0.00970421 | 2:-0.000438; 6:0.118486 | 0.125668 | 0.249565 | False |
| LCT | 0.0562828 | 0.00141461 | 2:0.278679; 6:0.348952 | 0.402233 | 0.335735 | False |
| FADS | 0.000108175 | 0.0191983 | 4:0.134373 | 0.386096 | 0.0343992 | False |
| SLC45A2 | -0.00197335 | -0.00197335 | 4:0.087169 | 0.0248002 | 0.965752 | False |
| SLC24A5 | 0.0158002 | 0.0158002 | 4:0.232819 | 0.436241 | 0.567876 | False |

[inference] Stability is worse, not better, than the previous painter precedent; ARI and that painter's L1/43% local agreement are different metrics and cannot be numerically equated. No retroactive choice of the smaller window or a convenient K is a valid rescue. Random projection, local LD variation and k-means fitting are possible causes, not established causes.

### Gate (c): exclusion and negative control — PASSED

[observation] Exact tested SNPs and pgen offsets, explicitly excluded:

| locus | index | chrom | pos | ID | ref | alt | excluded | eligible_SNP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LCT | 6678339 | 2 | 136608646 | 2_136608646_G_A | G | A | True | True |
| FADS | 34808594 | 11 | 61569830 | 11_61569830_C_T | C | T | True | True |
| SLC45A2 | 16783414 | 5 | 33951693 | 5_33951693_C_G | C | G | True | True |
| SLC24A5 | 42796964 | 15 | 48426484 | 15_48426484_A_G | A | G | True | True |

[observation] Canonical tested SNP and every +/-25kb site beyond the corrected gene core are absent from both feature lists. AMY entire103915000–104327000 CORRECTED core-plus-buffer is absent, not merely its top SNP. Exact zero-intersection checks, included index hashes and per-coordinate audits are saved. Software regression changes ALL forbidden dosages to2-dosage, rebuilds features, and obtains bit-identical exclusion-on sketches; exclusion-off sketches change. Independent permuted balanced pseudo-variant (199 null permutations):

| locus | excluded_sites | forbidden_features | negative_eta_squared | negative_permutation_P | gate_c |
| --- | --- | --- | --- | --- | --- |
| AMY | 2936 | 0 | 7.27094e-05 | 0.73 | True |
| LCT | 2486 | 0 | 3.38295e-05 | 0.91 | True |
| FADS | 2262 | 0 | 2.72827e-06 | 0.865 | True |
| SLC45A2 | 1660 | 0 | 0.000181779 | 0.075 | True |
| SLC24A5 | 1681 | 0 | 2.59816e-05 | 0.515 | True |

[inference] Negative controls pass, but exclusion-off ARIs are often low: it is NOT true that removing the locus changes nothing. The mandatory exclusion prevents direct circularity, while long-range LD can remain. Background sharing would not prove independent origin, even with validated clusters: recombination, migration and unphased partner chromosomes can place a single origin in several diplotype clusters.

### Gate (d): usable imputed coverage — PASSED

[observation] Every panel individual has100% observed dosage calls on every background feature, in BOTH windows for all FIVE loci. Unlike the failed raw-AIM lane, nobody is discarded for5–9% or<20% call rate. This is imputed input coverage, NOT direct ancient sequencing coverage or proof of imputation accuracy:

| locus | window | features | assigned | observed_95pct_n | observed_95pct_fraction | min_call_fraction | median_call_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| AMY | 200000 | 1221 | 16339 | 16339 | 1 | 1 | 1 |
| AMY | 1000000 | 12133 | 16339 | 16339 | 1 | 1 | 1 |
| LCT | 200000 | 622 | 16339 | 16339 | 1 | 1 | 1 |
| LCT | 1000000 | 9653 | 16339 | 16339 | 1 | 1 | 1 |
| FADS | 200000 | 848 | 16339 | 16339 | 1 | 1 | 1 |
| FADS | 1000000 | 11097 | 16339 | 16339 | 1 | 1 | 1 |
| SLC45A2 | 200000 | 1318 | 16339 | 16339 | 1 | 1 | 1 |
| SLC45A2 | 1000000 | 10943 | 16339 | 16339 | 1 | 1 | 1 |
| SLC24A5 | 200000 | 725 | 16339 | 16339 | 1 | 1 | 1 |
| SLC24A5 | 1000000 | 9690 | 16339 | 16339 | 1 | 1 | 1 |

### Structural sequence identity inventory (not a post-failure origin claim)

[observation] Sequence archive: 94 assembled haplotype contigs, 28 exact structure strings; all94 join exactly to H labels through `haplotype_assignments.txt`. The distance table uses a different graph-block nomenclature: all31 representative strings match its hap1/hap2 fields;30/31 representative contig IDs link to named labels (one unresolved). The961 ordered distance pairs comprise31 zero diagonal pairs and930 positive pairs, maximum edit distance0.624433. Stable contig identifiers (before updated interval suffix) are used ONLY for the representative-to-label join; this join is archived.

[observation] Label/structure inventory, using supplied distances rather than treating all copies of a coarse label as identical:

| label | assembled_contigs | exact_structure_strings | known_population_partitions | distance_representatives | within_label_distance_pairs | within_label_min_edit_d | within_label_max_edit_d |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1^a | 13 | 1 | 9 | 1 | 1 | 0 | 0 |
| H2A0 | 2 | 1 | 2 | 1 | 1 | 0 | 0 |
| H2A2B2 | 4 | 1 | 3 | 1 | 1 | 0 | 0 |
| H3A2 | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| H3A3B3 | 2 | 2 | 2 | 2 | 4 | 0 | 0.321152 |
| H3^r | 45 | 4 | 11 | 5 | 25 | 0 | 0.102113 |
| H4A2 | 6 | 4 | 5 | 5 | 25 | 0 | 0.299416 |
| H4A2B2 | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| H5 | 16 | 10 | 9 | 10 | 100 | 0 | 0.154757 |
| H7 | 3 | 2 | 2 | 2 | 4 | 0 | 0.118426 |
| H9 | 1 | 1 | 1 | 1 | 1 | 0 | 0 |

[inference] These observations distinguish identical graph structure (e.g. all13 assembled H1^a contigs share one string across9 known population partitions) from coarse labels that pool different strings (H5:10 strings/16 contigs). Neither structural similarity nor heterogeneity proves recurrent independent origin; mutation phylogeny and phase-linked flanking genealogy remain untested. This INPUT check is not an AMY within-background recurrence result.

### STOP verdict and canonical positive controls

[observation] **AMY interpretation STOPPED: gates(a) and(b) fail; (c)/(d) pass.** No AMY background-conditioned frequency, time slope, P value, or origin-recurrence result is generated. Canonical LCT, FADS, SLC45A2 and SLC24A5 individually fail the SAME structure/stability gates, so their gated temporal/recurrence tests also stop. Positive-control recovery is0/4 validated loci, not four demonstrated null alleles. This machinery has not shown the expected within-background temporal changes or sharing at the controls; released global selection maxima are not a substitute for that validation.

[inference] Full-panel imputation fixes the demonstrated COVERAGE bottleneck, but does not by itself validate these local diplotype backgrounds. The requested AMY within-fixed-background question remains UNRESOLVED; this run does not overturn or newly corroborate the carriage-dominated global-label result.

### Independent comparison against established global-label results

[observation] Fresh raw CN+archive dates reconstruct533 dated ancients and the supplied component means (agree at stated rounding):

| p2 | n | mean_AMY1 |
| --- | --- | --- |
| Bronze age | 46 | 6.63043 |
| CHG | 1 | 8 |
| EHG | 49 | 5.38776 |
| Early farmers | 22 | 5.45455 |
| Iron age to early modern | 333 | 7.24925 |
| Neolithic farmers | 34 | 7.14706 |
| Steppe pastoralists | 16 | 8.25 |
| WHG | 32 | 5.75 |

[observation] Reconstructed GLOBAL-label conditional slopes from the225-person cohort's actual H1/H2, with HC3 detection floor stated before slope interpretation (not using this run's clusters):

| cut | n | B_frequency | HC3_MDE80_per_kyr | within_slope | P | carriage_share |
| --- | --- | --- | --- | --- | --- | --- |
| primary_le1 | 225 | 0.842222 | 0.145764 | 0.0485754 | 0.349655 | 0.946114 |
| sensitivity_le2 | 225 | 0.742222 | 0.15238 | 0.0759689 | 0.162283 | 0.942277 |

[observation] Saved GLMM artifact check (not a refit): A<=2 |X-equivalent|=3.637970, GRM retention=96.955252%. Thus3.64 clears3.16 but not5.45; ~97% retention and~94% carriage reproduce. The supplied .049/.146 floor belongs to PRIMARY A<=1, not A<=2; A<=2 instead has within .075969 and floor .152380. Both established within slopes are below their floors.

[observation] Fresh stream of immutable RELEASED scan, using original audit windows, not the new exclusion-on backgrounds:

| locus | n | max_abs_X | above_5_45 |
| --- | --- | --- | --- |
| AMY | 2204 | 3.79496 | 0 |
| LCT | 881 | 16.5232 | 281 |
| FADS | 979 | 16.1534 | 182 |
| SLC45A2 | 1032 | 13.1496 | 221 |
| SLC24A5 | 840 | 5.99031 | 6 |

[observation] AMY2204/max |X|3.794956889/zero>5.45 reproduces exactly; canonical values agree with16.52/16.15/13.15/5.99 at supplied precision. No supplied biological number disagrees. These SNP statistics do not directly test structural copy number.

### Genome-wide cost, provenance and limits

[observation] Corrected five-locus computation after pvar cache: 89.98s of locus work, total run 93.07s, peak process RSS 835.0MiB with one BLAS/OpenMP thread. Initial (invalid-exclusion) run including a full52.38M-row coordinate stream was87.72s; each locus reads~10–13k SNP features across every16,339 individual. Sketch memory is~3.0MiB per window; a256-variant read is~16MiB before working copies. Primary feature coordinates and all assignments are saved.

[inference] Mechanical diplotype clustering is genome-wide-capable without phasing: measured average18.00s/locus suggests~5.00 CPU-hours for1000 comparable loci, excluding provenance/reporting and assuming cached coordinates/linear IO. A naive52.38M-variant-specific run would cost~29.9 CPU-years; window caching is necessary, and a separate +/-25kb exclusion is required for EVERY tested variant. This is an extrapolation, not a benchmark at that scale. Whole-panel phased haplotype analysis is NOT delivered. **Do not scale a biologically failed background model just because computation is cheap.**

[observation] Code: `analysis/haplotype_background.py`, `analysis/haplotype_background_report.py`; tests: `analysis/test_haplotype_background.py`; artifacts: `data/derived/haplotype_background/`. Commands use `.venv/bin/python`, `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`; biological main ran twice, report `--facts` freshly streams released statistics, report without that flag APPENDS this section. All software tests and full-panel artifact checks must pass separately from the failed biological gates. RunRecorder from `scripts/provenance.py` disables git probes BEFORE construction, hashes sources/preregistration/small inputs/outputs, and records large immutable input size/mtime without asserting a content hash. No raw data were edited; no git command, staging or commit was run.

## APPENDED report errata and final validation — 2026-10-09T17:55Z

[observation] Correction to the generated gate(a) explanatory sentence: **AMY silhouette .154051 and SLC45A2 .369565 PASS the .10 silhouette threshold**; the other three do not. ALL five still fail the published-group and country NMI thresholds, so the gate(a) verdict is unchanged. The table is authoritative; the earlier sentence saying every silhouette is below.10 is incorrect. The report generator now states the correct condition.

[inference] The old painter's L1/local agreement and this run's ARI are not comparable scales; neither run demonstrates the required stability. The phrase “worse, not better” above should not be read as a quantitative cross-method comparison.

[observation] Initial runs began17:36:31Z and17:39:57Z; corrected runs began17:47:59Z and17:50:10Z. The opening section's “initial run” timestamp denotes the initial-method REPETITION, not the earliest execution. Original preregistration17:32:36Z precedes both. Corrected locus statistics reproduce exactly across both corrected executions (excluding elapsed/RSS only). Nine software/artifact tests pass, including annotation-complete AMY/control intervals, all-panel assignment/exclusion checks, same forbidden-genotype perturbation invariance, floor units, K tie rule, projection determinism, fail-closed outcomes and scientific-result replication. A prematurely invoked report correctly refused incomplete biological output; its failure log is retained, and the completed report succeeds.

[observation] Authoritative repeated-run timing overlapped an independent released-statistics stream; its89.98s locus-work benchmark is conservative relative to the53-ish-second isolated corrected run, not a pure isolated performance measurement. Computational feasibility does not rescue failed validation. Unphased diplotypes collapse different cis arrangements and cannot attach H1 versus H2 to an individual flanking chromosome; sharing across clusters cannot count mutational origins.
