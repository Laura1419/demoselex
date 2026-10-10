# ELAI AMY fit results — gated STOP

[observation] **Verdict: STOP. Known-structure and stability gates FAILED.** Coverage, virtual-site arithmetic/integration, and exclusion gates passed. No allele-frequency analysis, temporal analysis, biological state interpretation or claim about AMY was performed. This is a bounded pilot on real imputed genotypes, not a full-panel validation.

## Input and instrument

[observation] Immutable source dimensions: 52,382,872 variants × 16,339 samples. Streamed 8,635,540 records, stopping at chr3. Exact AADR-1240K annotation/.ind matches, BP>0, deduplicated by annotation Individual ID: 5,981 eligible; deterministic random cohort 128 dated ancients, no modern references. 111 published groups (101 singleton groups); 31 country labels, 0 unknown countries. These sparse labels can inflate unadjusted NMI; NMI alone is not evidence of recovery.

[observation] Chr1 95,000,000–115,000,000; inclusive core 103,915,000–104,327,000 excluded. Fixed 2-kb candidate thinning gave 9,330 biallelic SNP candidates; cohort hardcall MAF≥0.05, call rate≥0.95 gave **1,164 fitted markers**, not the aspirational tens of thousands. Minimum retained MAF 0.05078125, minimum call rate 1.000000. Alternate-marker thinning 582 markers. Thinning before QC often selected rare variants; we did not revise the preregistered rule after seeing results. Independent chr2 PC panel 145 SNPs (200-kb candidate spacing), first five standardized-genotype SVD PCs. This sparse PC panel and bounded cohort are limitations.

[observation] BIMBAM unphased diploid headers N/M/IND, allele strings mapped from native ALT hardcalls (not rounded fractional dosages). The imputed data's complete calls do not establish original aDNA coverage or imputation accuracy. Input panel.inp: 466,593 bytes; panel.pos: 31,982; thin.inp: 233,944; thin.pos: 15,990. No matrix persisted over 500 MB. Full primary --pkk streamed 3,724,800 rows through a FIFO; summary checks verified every numeric sample index, fitted marker ID and ordered row sequence. Upstream leaves individual_id blank despite the manual's optional IND header; samples.tsv plus unchanged input-column order is the authoritative index-to-sample mapping. Only virtual-site tensors, small composition summaries, parameters and logs retained.

[observation] Patched binary /home/erikg/elai-fork/src/elai, no rebuild; operator identifies its public fork as github.com/ekg/ELAI, default branch lower-state-posterior. Actual binary SHA256 is recorded in provenance, without relying on an unverified branch identity. Manual-derived settings: -p 1 (cohort), -C 1 (requested linked limit; upper state is constant), -c 5 (manual 5*C), sensitivities 3/8, -w 0 -s 20 (manual recommended exact EM 20–50), -mg 10 fixed prior scale, -R recorded seeds. pk1 COPY COUNTS divide by 2; ordered pkk sums to 1; q sums to 1 only on k1≤k2. The final likelihood still improves after 20 steps (primary last increment 119.681); convergence was not established. We do not tune EM or retry toward acceptance after observing gate failure.

## Numeric gates

| Gate | Observed | Preregistered requirement | Result |
|---|---|---|---|
| (a) Known structure | composition silhouette 0.176180; independent PC silhouette **-0.057757**; group NMI 0.437534; country NMI 0.205049; 5 nonempty assignments | silhouette and PC silhouette >0.10; group NMI >0.10; country NMI >0.05 | **FAIL: PC silhouette** |
| (b) Stability | 8/8 comparisons completed; see below | seed/thinning ARI≥0.80, cross-K ARI≥0.60, same-K L1≤0.20, virtual agreement≥0.80, cross-K normalized Gram difference≤0.30 | **FAIL** |
| (c) Coverage | 128/128 usable (100%); 128 dated ancients in complete site output; 3,724,800 complete all-marker rows | ≥95% individuals with ≥95% calls; ≥100 dated exports; all-marker tensor complete | PASS |
| (d) Virtual sanity | max ordered sum error 1.11e-15; q symmetry error 0; triangular q error 8.88e-16; linked/unlinked max change **0.769908793** | ≤1e-10, ≤1e-12, ≤1e-10 respectively; finite probabilities [0,1]; change >1e-6 | PASS |
| (e) Exclusion | 0 builder/retained core markers, 0 PC overlap; site 104121000 inside core; left 103783823 < site < right 104341989 | zero excluded markers; inside core; same-chromosome flanks | PASS |

[observation] Eight source-named helper regressions pass (tests.log); independent artifact validation confirms all nine executions, registration before fitting, exact zero-core-marker exclusion, complete 128-individual exports, q formula error 0, and exactly one edited control token (independent_validation.json).

[observation] Other structure diagnostics (no minimum date correlation preregistered): date-bin NMI 0.067200, composition-predicted BP R² 0.056891; published-group/country PC silhouettes -0.126225/-0.358233, composition-space silhouettes -0.149310/-0.426456. These do not support an independent-PC structure pass. Earlier painter silhouette<0.10/group NMI<0.10/country NMI<0.05 was a failed precedent; beating its NMI values is insufficient here.

### Stability numbers

| Comparison | ARI | Matched mean L1 | Virtual argmax agreement | Result |
|---|---:|---:|---:|---|
| K=3, seeds 20261010/20261011 | 0.013642 | 0.103133 | 0.046875 | FAIL |
| K=5, seeds 20261010/20261011 | -0.021952 | 0.168408 | 0.078125 | FAIL |
| K=8, seeds 20261010/20261011 | -0.025109 | 0.209265 | 0.078125 | FAIL |
| K=5, original/alternate thinning | 0.021840 | 0.148711 | 0.148438 | FAIL |

| Cross-K (seed) | ARI | Normalized Gram difference | Result |
|---|---:|---:|---|
| K=3/5 (20261010) | -0.000130 | 0.013697 | FAIL |
| K=5/8 (20261010) | 0.019577 | 0.019757 | FAIL |
| K=3/5 (20261011) | 0.014322 | 0.017461 | FAIL |
| K=5/8 (20261011) | 0.011870 | 0.021593 | FAIL |

[observation] Composition-based Hungarian matching is applied only at the same K. Cross-K comparison uses ARI and normalized sample co-composition Gram matrices, not arbitrary label correspondence. Low Gram differences with poor ARI show why marginal similarity alone cannot pass stability. Virtual agreement is the lower marginal argmax, not an inferred phased haplotype. Prior painter's failed 43% target agreement is not rescued: this pilot gives 7.8125% seed agreement at K=5 and 14.84375% thinning agreement.

## Exclusion and virtual control demonstrated

[observation] Builder and runner assert inclusive exclusion; the runner rechecks every actual retained snpinfo record in all nine successful executions, and checks the virtual locus/coordinate and complete K*K rows for all individuals. Virtual ordered-copy symmetry max error 6.66e-16. Frozen control preserves all EM-file lines and changes only rk at the right core-gap flank to 1 (fully unlinked); original fitted file is replayed with s=0 as the linked rounded-parameter baseline. The two outputs differ by 0.769908793 maximum absolute ordered probability on this input, confirming linkage sensitivity, not biological accuracy of the missing-emission bridge. Parameters, edited control, commands and outputs are retained.

## Exact reproducible executions

[observation] Run from repository root with .venv/bin/python, OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1:
```
.venv/bin/python analysis/elai_amy_build.py
.venv/bin/python analysis/elai_amy_gate.py
.venv/bin/python analysis/elai_amy_gate.py --metrics-only
PYTHONPATH=analysis .venv/bin/python -m unittest analysis/test_elai_amy.py -v
.venv/bin/python analysis/elai_amy_report.py
```
The runner creates/deletes a scratch working directory because ELAI always writes relative output/. The exact binary commands below use persistent input paths; prefix uniquely identifies K/panel/seed. Successful upstream exit code is 1; finished sentinel, output completeness and numeric checks are also required.

[observation] frozen_linked: 1164 retained markers, 0.395 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 5 -w 0 -s 0 -mg 10 -R 20261010 -o frozen_linked --pkk-site 1:104121000 -rem /home/erikg/demoselex/data/derived/elai_amy/panel_k5_seed20261010.em.txt
```

[observation] frozen_unlinked: 1164 retained markers, 0.394 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 5 -w 0 -s 0 -mg 10 -R 20261010 -o frozen_unlinked --pkk-site 1:104121000 -rem /home/erikg/demoselex/data/derived/elai_amy/unlinked_control.em.txt
```

[observation] panel_k3_seed20261010: 1164 retained markers, 1.951 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 3 -w 0 -s 20 -mg 10 -R 20261010 -o panel_k3_seed20261010 --pk1 --pkk-site 1:104121000
```

[observation] panel_k3_seed20261011: 1164 retained markers, 1.951 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 3 -w 0 -s 20 -mg 10 -R 20261011 -o panel_k3_seed20261011 --pk1 --pkk-site 1:104121000
```

[observation] panel_k5_seed20261010: 1164 retained markers, 14.555 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 5 -w 0 -s 20 -mg 10 -R 20261010 -o panel_k5_seed20261010 --pk1 --pkk-site 1:104121000 --pkk
```

[observation] panel_k5_seed20261011: 1164 retained markers, 4.112 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 5 -w 0 -s 20 -mg 10 -R 20261011 -o panel_k5_seed20261011 --pk1 --pkk-site 1:104121000
```

[observation] panel_k8_seed20261010: 1164 retained markers, 8.943 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 8 -w 0 -s 20 -mg 10 -R 20261010 -o panel_k8_seed20261010 --pk1 --pkk-site 1:104121000
```

[observation] panel_k8_seed20261011: 1164 retained markers, 8.930 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/panel.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/panel.pos -C 1 -c 8 -w 0 -s 20 -mg 10 -R 20261011 -o panel_k8_seed20261011 --pk1 --pkk-site 1:104121000
```

[observation] thin_k5_seed20261010: 582 retained markers, 2.005 seconds, exit 1.
```sh
/home/erikg/elai-fork/src/elai -g /home/erikg/demoselex/data/derived/elai_amy/thin.inp -p 1 -pos /home/erikg/demoselex/data/derived/elai_amy/thin.pos -C 1 -c 5 -w 0 -s 20 -mg 10 -R 20261010 -o thin_k5_seed20261010 --pk1 --pkk-site 1:104121000
```

## Execution audit, resources, validation and decision

[observation] Registration hash ae4ae3e2dd0942abd40e43bbcfbd3b618e046803d0d2dde48bf77790c614d460 is unchanged across runs; its filesystem mtime precedes every recorded final fit start. Three execution-only repairs are preserved: first parser assumed the manual's six-column snpinfo description, but this binary actually writes rs,pos,chr,major,minor,maf,miss,n0,n1,n2; second control failed because parameter files are line-oriented and flattening tokens breaks -rem; a repeat-audit assertion initially required textual IDs even though README explicitly warns BIMBAM IDs may be empty. Corrected parsing/line preservation and authoritative numeric-index auditing, with **no changes to cohort, markers, model settings or gates**. Historical failure logs/metrics are explicitly named first/control/pre_control. A complete final rerun reproduced 16/16 site/composition file SHA256 hashes byte-for-byte.

[observation] Starting, pre-write and final df -h / show 71G free (rounded 100% use on shared 14T volume); final metric-time free 75,762,126,848 bytes. Persistent lane output approximately 5.46 MB (max file 0.68 MB), far below 5GB. Scratch removed and deletion recorded in scratch_cleanup.json; no raw inputs modified, no installs, no rebuild, no git mutation/staging. Final permitted read-only `git diff --cached --name-only` returned no entries (no staged files); no git command ran in elai-fork. Provenance uses scripts/provenance.py with git provider disabled before construction; small files and binary hashed, huge immutable pgen/pvar recorded by size/mtime rather than claiming fresh content hashes.

[inference] **STOP is the outcome, not a biological negative result.** Required before a new attempt: prospectively registered denser QC-then-thin markers, more representative/larger cohort and independent genome-wide PC validation, adequate EM convergence, and then renewed seed/K/thinning validation with the same fail-closed discipline. This lane does not establish that ELAI generally fails, nor that any AMY history is absent. No local-state biological summary or effective independent-state estimate is justified because gates failed. Temporal and allele-frequency analyses remain blocked.
