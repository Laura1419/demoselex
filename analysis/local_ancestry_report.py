#!/usr/bin/env python3
"""Append numeric Stage-1 results; never overwrite earlier authors' conclusions."""
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.local_ancestry import OUT


def table(frame):
    def fmt(x):
        if isinstance(x, (float, np.floating)):
            return 'NA' if not np.isfinite(x) else f'{x:.6g}'
        return str(x)
    lines = ['| '+' | '.join(map(str, frame.columns))+' |', '| '+' | '.join(['---']*len(frame.columns))+' |']
    lines += ['| '+' | '.join(fmt(x) for x in row)+' |' for row in frame.itertuples(index=False, name=None)]
    return '\n'.join(lines)


def main():
    v = json.loads((OUT/'validation.json').read_text())
    facts = json.loads((OUT/'verified_fact_counts.json').read_text())
    baseline = json.loads((OUT/'verified_global_baseline.json').read_text())
    audit = json.loads((OUT/'verified_akbari_windows.json').read_text())
    refs = pd.read_csv(OUT/'reference_counts.tsv', sep='\t')
    panel = pd.read_csv(OUT/'aim_panel.tsv', sep='\t')
    confusion = pd.read_csv(OUT/'confusion.tsv', sep='\t')
    held = pd.read_csv(OUT/'holdout_validation.tsv', sep='\t')
    confusion['unevaluable'] = [int(((held.source == s) & ~held.evaluable).sum()) for s in confusion.published_source]
    sensitivity = pd.read_csv(OUT/'sensitivity.tsv', sep='\t')
    means = pd.read_csv(OUT/'verified_CN_means.tsv', sep='\t')
    now = datetime.now(timezone.utc).isoformat()
    failed = ', '.join(k for k, val in v['gates'].items() if not val)
    summary = pd.DataFrame([dict(quantity=col, minimum=panel[col].min(), median=panel[col].median(), maximum=panel[col].max())
                            for col in ['divergence', 'Fst_equal_source', 'target_call_fraction']])
    previous = pd.DataFrame([dict(cut=k, n=baseline[k]['n'], B_frequency=baseline[k]['frequency'],
                                  within_slope=baseline[k]['beta_time'], HC3_SE=baseline[k]['HC3_se'],
                                  P=baseline[k]['p'], MDE80_per_kyr=baseline[k]['MDE80_HC3'],
                                  pooled_slope=baseline[k]['beta_pooled'], carriage_share=baseline[k]['carriage_share'])
                             for k in ['primary_le1', 'sensitivity_le2']])
    text = f'''\n## NEW COMPLETED VALIDATION / STOP VERDICT — {now}\n
### Reference panel and target cohort

[observation] Exact-label raw counts, exclusions and eligible BP ranges (years before1950) follow. Assessment=`Pass*`, all annotated family relations excluded, biological-ID duplicates collapsed using 1240K-hit priority. Eligible totals:105 references,78 training/27 held out. Labels are empirical brackets of variation, not assertions of discrete replacements.

{table(refs)}

[observation] Structural archive has{v['targets_in_archive']} dated labelled ancient diplotypes;{v['targets_mapped']} unique biological individuals match AADR,{v['targets_excluding_training']} remain after training-reference exclusion,{v['targets_global_painted']} meet genome-wide coverage criteria, and{v['targets_AMY_supported']} meet the primary AMY-vicinity support rule. Exactly383 reference+target individuals were computationally painted; only19/27 holdouts were evaluable. These are NOT the533 copy-number ancients or the previous225-person imputed cohort. The{239} supported AMY targets are diagnostic coverage counts, NOT a tested cohort. `target_mapping.tsv`, `references.tsv`, `target_composition.tsv` provide individual IDs and all denominators.

### AIM panel and painter

[observation] {v['candidate_AIMs']:,} autosomal candidates pass source divergence/call criteria;{v['AIMs']:,} markers remain after target call>=.2 and50kb greedy thinning. Full source ALT frequencies, source/target call fractions, divergence and equal-source descriptive Fst for EVERY selected marker are saved in `aim_panel.tsv`.

{table(summary)}

[inference] Method:15 unordered ancestry-pair states;2Mb windows; diploid two-source genotype probabilities and one-allele Bernoulli pseudo-haploid probabilities; emission logs tempered by.25; row-normalized transition weights exp(-4 * minimum changed chromosome slots); uniform initial pair prior and resets between chromosomes. Forward–backward posterior chromosome fractions and Viterbi pairs are saved in `primary_paint.npz`;<5 observed SNPs produces flat emissions and excluded support. Genome-wide evaluation needs100 supported windows and20% AIM call rate. AMY support needs10 observed SNPs in supported overlapping windows. This is explicitly NOT a phased ancestry-segment assignment of H1/H2.

[inference] Known approximations: source frequencies estimated from few partially missing ancient samples; SNP independence despite LD overcounts evidence and overconfident segments;50kb thinning/.25 tempering are not LD calibration; physical-window switches are not genetic-map recombination; no genotype likelihoods, phase linkage, drift calibration or source-frequency posterior integration. A historically admixed STEPPE state and non-CHG-equivalent IRAN state bracket variation but can be mutually confounded with ancestry pairs.

### Required validation: FAILED, not waived by dominant-label accuracy

[observation] Confusion matrix: rows published source, columns inferred dominant global source; unevaluable column is excluded from accuracy. There are ZERO covered IRAN holdouts, not perfect IRAN classification.

{table(confusion)}

[observation] Evaluable accuracy19/19=100%, but only19/27=70.37% of held-out individuals were evaluable. Correct-source median fractions: WHG=.957974, EHG=.514091, FARM=.742684, IRAN=NA, STEPPE=.937793. Both IRAN holdouts have AIM call rates.088445 and.052945 (<.2); WHG1 and FARM5 holdouts also fail coverage. EHG median violates the prespecified.60 predominance threshold even though its largest component is EHG. The failed gates are **{failed}**. Covered-subset classification cannot rescue missing-source and unstable-proportion validation.

[observation] Negative control selected WITHOUT ancestry results: highest-call-rate evaluable held-out WHG `hoe004.SG`; observed AIM call fraction=.962847, genome-wide WHG fraction={v['negative_control']['WHG_fraction']:.6f} (passes.80). PC panel contains{v['pc_markers']:,} disjoint,200kb-spaced polymorphic markers after missing/constant exclusion. R² of PCs1–3 predicted by global local-call proportions: {', '.join(f'{x:.6f}' for x in v['PC_R2'])}; mean={v['PC_mean_R2']:.6f} (passes.50). This is descriptive in-sample association, NOT held-out ancestry truth or proof of segment accuracy. Marker coordinates/PC scores and fractions are archived.

[observation] All nine preregistered settings follow. L1 is SUM of absolute source-fraction differences, not percent agreement. Stability gate requires EVERY alternative held-out median L1<=.25 and dominant-source agreement>=.8. It fails1Mb/2 (.779726),2Mb/2 (.544155), and4Mb/2 (.250815); dominant GLOBAL labels remain unchanged but ancestry proportions do not. Target-local dominant agreement falls to43.31% for1Mb/2 despite stable global labels.

{table(sensitivity)}

[inference] Dominant labels and PC correlations alone would misleadingly suggest success. This prototype has not demonstrated stable local ancestry proportions, much less reliable haplotype-segment backgrounds. No parameters or eligibility rules were changed to force a pass. The final-code biological rerun reproduces every finite validation quantity and sensitivity row exactly; null replaces unavailable IRAN statistics in strict JSON.

### AMY test and detection floor: deliberately NOT RUN

[observation] `amy_status.json` records STOPPED. No `amy_statistics.tsv`, local-ancestry-conditioned AMY outcome frequencies, time coefficients or P values were generated. The implemented post-gate regression/frequency/power path is unreachable in this failed biological run. **There is no valid local-conditioned detection floor to report**, because there is no validated analysis design/test.239 covered targets do not by themselves establish power after ancestry conditioning. Reporting an approximate floor based on unstable inferred fractions would violate the stop rule and invite interpretation of a broken painter.

[inference] **DOES THE AMY SIGNAL PERSIST WITHIN A FIXED LOCAL-ANCESTRY BACKGROUND? UNRESOLVED.** This local view provides no validated biological comparison and therefore DOES NOT change or corroborate the prior carriage-dominated conclusion. Power sufficiency of a valid local-background test remains unknown. Even a successful unphased dosage~local-composition regression would answer persistence conditional on estimated diploid composition, not assign either structural allele to a fixed phased chromosome ancestry. A genuine fixed-background inference requires phase/linkage validation beyond this approximation.

### Carried-forward facts independently rechecked (NOT local-painter AMY results)

[observation] Raw copy-number table + archive age metadata reproduce all supplied means on{facts['dated_CN_ancients']} dated CN ancients:

{table(means)}

[observation] Independently reconstructed the prior225 matched structural diplotypes from saved cohort/dosages, refitted GLOBAL-label HC3 models, and recomputed pooled/within channel algebra. Detection floor is stated in this baseline table BEFORE interpreting its nonsignificant within effect; units are B-frequency/kyr, forward time=-BP/1000.

{table(previous)}

[observation] Primary previous B frequency=.842222; within=.0485754/kyr versus80%-power floor=.145764/kyr (2kyr floor=.291528), P=.349655; carriage=.946114. A<=2 within=.0759689, floor=.152380, carriage=.942277. Thus the supplied.049/.146 and~94% are primary-A<=1 quantities, not the A<=2 quantities associated with X-equivalent3.64. Saved GLMM values (artifact verification, NOT fresh GLMM fitting): A<=2 |X-equivalent|={baseline['saved_GLMM']['sensitivity_le2_X_equivalent']:.6f}; primary GRM retention={baseline['saved_GLMM']['primary_GRM_retention']:.6%}. These reproduce3.64 and~97%; A<=2 clears3.16 but not5.45. Both within slopes remain below their floors; this is an underpowered baseline null, not proof of absence.

[observation] A fresh immutable released-statistics stream reproduces all canonical windows:

{table(pd.DataFrame([dict(locus=k, **value) for k, value in audit.items()]))}

[observation] Previous matched count225/288 is verified. Raw CN+AADR mapping finds{facts['CN_AADR_dated_coords']} uniquely matched dated-coordinate ancients, versus the supplied522: **disagreement={facts['CN_AADR_dated_coords']-522:+d}**, under this explicit genetic/biological/unique-prefix mapping of533 finite-CN/age ancients. Structural mapping278 is a DIFFERENT cohort, not a replacement for the522 statement. Full mapping flags are saved. Exact canonical |X| values agree with supplied rounded values; AMY2204/3.794956889/zero above5.45 reproduces. No modern/full genome local-ancestry application was attempted.

### Reproduction and limits

[observation] Commands: `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python analysis/local_ancestry.py`; same prefix with `analysis/test_aadr_orientation.py`, `analysis/local_ancestry_verify_facts.py`, `-m unittest discover -s analysis -p test_local_ancestry.py -v`, and `analysis/local_ancestry_report.py` (APPENDS). Logs and first-run snapshots are under `data/derived/local_ancestry/`; manifests under `analysis/runs/` use `scripts/provenance.py` with all git metadata probes explicitly disabled. Source/preregistration/small-output hashes and large-input size/mtime records are provided; no large-input content hash is asserted.

[observation] Orientation reproduces low-bit-first r=-.1168, high-bit-first unflipped=-.8103, high-bit-first ALT inversion=+.8103, agreement=.7223. Six software regressions cover one-versus-two-allele weighting, pair transitions, exact-enumeration forward–backward/Viterbi, missing support/chromosome reset, conditional HC3/detection floor, prereg chronology/AIM rules/confusion counts/fail-closed outcomes. Software tests do NOT validate biological ancestry. An early coding correction ensured primary AMY coverage counts use PRIMARY rather than last sensitivity windows; final-code rerun is authoritative. An initial artifact test was explicitly skipped while biological outputs were incomplete; the completed run must pass it.

[inference] Stage1 delivers reproducible machinery and an honest failed validation, NOT a validated painter or AMY answer. A separately preregistered revision should improve reference coverage, investigate EHG/STEPPE and window/transition identifiability, and validate segment accuracy with known mosaics plus phase-linked structural calls. This run must not be tuned retroactively or scaled genome-wide.
'''
    with Path('docs/local-ancestry-results.md').open('a') as f:
        f.write(text)
    print('Appended final report; validation failed and AMY stopped.')


if __name__ == '__main__':
    main()
