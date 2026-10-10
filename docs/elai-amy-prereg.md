# ELAI AMY instrument gate — preregistration

[observation] Registered 2026-10-10 before any ELAI fit in this lane. This is an instrument-validation run, not a biological analysis. Sources read: patched README-PATCH.md, all of elai-manual.tex, the eight-page PDF via pypdf, example/note.txt and example genotype/position headers. No rebuild or installation.

## Input and bounded cohort

[observation] GRCh37 chromosome 1, inclusive 95,000,000–115,000,000. Exclude **inclusive 103,915,000–104,327,000** before reading fitted genotypes; virtual missing-emission site **1:104121000**, inside that core. No genotype at that site is used. Stream pvar with zero-based variant indices. Only biallelic single-base A/C/G/T SNPs; first SNP at least 2,000 bp after the previous candidate, then cohort hardcall MAF >=0.05 and call rate >=0.95. No ancestry-informative selection. A second panel takes alternate retained markers (roughly 4 kb spacing). Counts will be recorded, not adjusted to improve gates.

[observation] Cohort: up to 128 dated ancient individuals, date BP >0, exact psam-to-AADR-1240K annotation and .ind ID matches, one representation per AADR Individual ID. Sample without replacement using numpy default_rng seed 20261010 from psam-order eligible representatives, then restore psam order (uint32 pgen subset). All selected samples are cohort -p 1, no named training panel. Labels/country/date do not enter EM or marker selection. This small random cohort is a bounded validation pilot, not the full ancient panel. Report eligible pool, fitted N, usable calls and dated exports separately.

[observation] Use pgen hardcalls from the imputed release, rather than rounding fractional ALT dosages: ELAI's documented input accepts unphased diploid allele strings, not genotype likelihoods. 0/1/2 ALT copies map REFREF/REFALT/ALTALT, -9 to ??. First header is N without '=', second M, optional IND identifiers; positions are `ID position 1`, exactly as manual/examples. Imputed call rate is not original aDNA read coverage or imputation certainty.

[observation] Independent PC validation: select chromosome 2 SNP candidates at >=200 kb spacing (stream until chromosome 3), same SNP/MAF/call-rate filters and same cohort; standardized hardcalls with cohort-mean missing imputation, deterministic SVD, first 5 PCs. These markers cannot overlap fitted chromosome 1. AADR .ind supplies published group labels; annotation supplies country and dates. Date categories are floor(BP/2000), used only for validation.

## Model and commands

[observation] Primary S=1 (`-C 1`), K=5 (`-c 5`); sensitivities K=3 and K=8. Two recorded seeds 20261010 and 20261011 for each K; alternate-marker panel K=5 seed 20261010. Exact quadratic EM `-w 0 -s 20`, within manual's 20–50 recommendation; no fast-linear-only approximation. `-mg 10` is a fixed model default-scale choice, not an ancient admixture-date estimate. K=5 follows manual's 5*C recommendation; S=1 is the requested linked-model reference-free limit. With S=1 upper ancestry is identically 1 and cannot validate population structure: gate uses lower-state composition.

[observation] Command template (working directory is disposable scratch, not ELAI repo):
```
/home/erikg/elai-fork/src/elai -g <absolute-panel.inp> -p 1 -pos <absolute-panel.pos> -C 1 -c <K> -w 0 -s 20 -mg 10 -R <seed> -o <prefix> --pk1 --pkk-site 1:104121000
```
Primary additionally `--pkk`. Stream its full-marker tensor through a FIFO into bounded diagnostics; never persist the full tensor. Save the small virtual-site TSV, fitted parameters, logs, panel manifest and sample compositions. `--pk1` is expected diploid COPY COUNTS summing to 2, divide by 2; `--pkk` ordered square sums to 1; q is repeated symmetrically, sum only its upper triangle. Upstream exits 1 on success: require finished sentinel and expected complete files, not exit status alone.

## Five gates (fixed before fitting)

[inference] All thresholds are pilot acceptance rules, not biological truth guarantees. Fail-closed: unavailable/unfinished tests fail, never count as passes.

1. **Known structure:** primary mean lower composition argmax assignments: silhouette in composition space >0.10; group NMI >0.10; country NMI >0.05; independent chr2 PC-space silhouette of assignments >0.10. At least two nonempty assignments. Also report group/country silhouettes in PC and composition spaces, date-bin NMI and date-versus-state-composition explained variance (descriptive; no minimum date correlation required). Must beat failed painter precedent, not merely match it.
2. **Stability:** ARI >=0.80 across seeds at each K; cross-K ARI >=0.60 (K3/K5 and K5/K8, each seed); alternate thinning ARI >=0.80. Same-K compositions compared after Hungarian state matching: mean per-individual L1 <=0.20, virtual-site argmax agreement >=0.80. Cross-K composition compare co-assignment/Gram similarity (report normalized Frobenius difference <=0.30); labels have no direct cross-K meaning. All planned comparisons required.
3. **Coverage:** >=95% of fitted individuals have >=95% usable hardcalls; >=100 dated ancients with complete virtual-site output. Primary all-marker export must include all N individuals and all fitted marker sites with K*K rows per site; no silent dropped individuals.
4. **Virtual sanity:** maximum ordered sum error <=1e-10; q symmetry error <=1e-12 and unique-triangle sum error <=1e-10, all finite probabilities in [0,1]. Reconfirm integration on this input with a frozen-parameter control: set fitted lower-switch probability rk at the core gap to 1 (unlinked), leaving all other parameters unchanged and no EM. Compare against a frozen rounded-parameter linked replay; require max absolute joint change >1e-6. This is instrument sensitivity, not bridge accuracy. Ordered copy symmetry also reported. No observed core genotype is inserted.
5. **Exclusion:** programmatic zero fitted markers in inclusive core (both panels and ELAI retained snpinfo), virtual coordinate inside core, same-chromosome left/right retained flanks, zero chr1/chr2 PC marker overlap.

## Decision and resource rule

[inference] If any gate fails, **STOP**, report exact numbers and thresholds and what is missing; no allele frequencies, no biological interpretation, no temporal analysis. Even all gates passing would permit only local-state count/size/composition and effective-state summary, not allele-frequency claims.

[observation] Wall-clock budget ~35 minutes. Before large writes run df -h /; starting available space 71G (14T volume at rounded 100%). New persistent output <5 GB, any written matrix <500 MB. Scratch deleted at end. If resources/time prevent complete validation, report test counts and unmeasured fields explicitly, FAIL the affected gates. No git mutations or commands; provenance git provider disabled before construction. Raw tiers are read-only.
