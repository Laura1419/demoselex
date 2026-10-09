#!/usr/bin/env python3
"""Append the marked final CN-as-SNP report and save descriptive verification tables."""
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from cn_as_snp import OUT, HAPS


def markdown(headers, rows):
    return '| ' + ' | '.join(str(x).replace('|', r'\|') for x in headers) + ' |\n|' + '|'.join(['---'] * len(headers)) + '|\n' + ''.join(
        '| ' + ' | '.join(str(x).replace('|', r'\|') for x in row) + ' |\n' for row in rows)


def main():
    state = json.loads((OUT / 'summary.json').read_text())
    if state['status'] not in ['complete', 'complete_with_fit_failures']:
        raise RuntimeError('primary run has not finished')
    diag = json.loads((OUT / 'diagnostics.json').read_text())
    extra = json.loads((OUT / 'unrestricted_null.json').read_text())
    by = {(r['locus'], r['model']): r for r in state['results']}
    unrestricted = {(r['locus'], r['model']): r for r in extra['results']}
    h9 = diag['H9_GRM_off']
    # Carry-forward checks use all released copy-number ancients, not the smaller
    # categorical-haplotype subset. Keep denominators explicit.
    copy = pd.read_csv(OUT / 'copy_number_table_for_supplement.tsv', sep='\t')
    hap = pd.read_csv(OUT / 'haplotype_structure_for_supplement.tsv', sep='\t')
    meta = pd.read_csv(OUT / 'ancient_metadata_for_supplement.tsv', sep='\t')
    ca = copy[copy.type.eq('Ancient')]
    means = ca.groupby('p2').AMY1.agg(['mean', 'count']).reset_index()
    means.to_csv(OUT / 'carried_forward_AMY1_means.tsv', sep='\t', index=False)
    a = hap[hap.type.eq('Ancient')].merge(meta[['sample', 'age']], on='sample', validate='one_to_one')
    groups = []
    for group, part in a.groupby('p2'):
        row = {'component': group, 'n': len(part), 'mean_age_BP': part.age.mean()}
        for h in HAPS:
            row[h] = (part.H1.eq(h).sum() + part.H2.eq(h).sum()) / (2 * len(part))
        groups.append(row)
    groups = pd.DataFrame(groups).sort_values('mean_age_BP', ascending=False)
    groups.to_csv(OUT / 'carried_forward_ancient_haplotype_frequencies.tsv', sep='\t', index=False)
    modern = []
    mh = hap[hap.type.eq('Modern')]
    for pop, part in mh.groupby('p1'):
        row = {'population': pop, 'n': len(part), 'region': part.p2.iloc[0]}
        for h in HAPS:
            row[h] = (part.H1.eq(h).sum() + part.H2.eq(h).sum()) / (2 * len(part))
        modern.append(row)
    modern = pd.DataFrame(modern)
    modern.to_csv(OUT / 'carried_forward_modern_haplotype_frequencies.tsv', sep='\t', index=False)

    text = '\n## New final analysis section — ' + datetime.now(timezone.utc).isoformat() + '\n\n'
    text += '[observation] **The amended ordering gate passed.** Both canonical loci are much stronger than the released lead AMY-region real SNP in the same 225 ancients. Absolute reproduction of the paper was NOT achieved and is not asserted; the supervisor explicitly authorized the ordering gate before fitting.\n\n'
    controls = []
    for n in ['LCT', 'SLC45A2', 'AMY_real_SNP']:
        r, off = by[n, 'GRM_on'], by[n, 'GRM_off']
        controls.append([n, state['references'][n]['ID'], f"{r['released_X']:.5f}",
                         f"{r['z_score']:.4f}", f"{r['X_equivalent']:.4f}",
                         f"{off['z_score']:.4f}", f"{r['empirical_P']:.3f}", f"{off['empirical_P']:.3f}"])
    text += markdown(['Control', 'Released allele (GRCh37)', 'Released X', 'Z on', 'X-equiv on', 'Z off', 'P ancestry on', 'P ancestry off'], controls)
    text += '\n[observation] Canonical GRM-on |Z| values are 8.75x (LCT) and 13.75x (SLC45A2) the lead AMY SNP, exceeding the prespecified 1.25x separation and |Z|>=2 requirements. Both fits converged. Exact released control alleles—not regional maxima substituted at different positions—were used. ALT time-BP signs are retained: positive Z means enrichment in older samples. Signed Z differs from the released canonical X signs; the exact signed paper-slope convention is not reproduced here. ALT dosage orientation follows the supplied imputed-bundle reader, not an AADR allele substitution. Only absolute values enter the gate and extrapolation.\n\n'
    text += '[observation] **Genome-wide individual-level GRM:** 1,133 polymorphic background SNPs, all 22 autosomes represented, sampled by a fixed 10,000-variant index grid then filtered for MAF>=0.05/call rate>=0.95. AMY K excludes all 75 chr1 markers, leaving 1,058 genome-wide background SNPs. LCT and SLC45A2 use the same construction rule/cohort with their own LOCO exclusion. This is not a GRM inferred from AMY and not a trajectory/group-level proxy. Marker list and actual matrices are saved. No LD pruning, paper QC recreation or paper fixed batch covariates are claimed; a sparse NoQC-bundle approximation may underrepresent structure. All control ALT dosages were integer-like (fractional counts in summary); the fitted response uses explicit ALT hardcalls, not rounded dosages.\n\n'
    text += '### Synthetic loci: the substantive GRM-on versus GRM-off comparison\n\n'
    rows = []
    for h in HAPS + ['AMY1_continuous']:
        on, off_primary = by[h, 'GRM_on'], by[h, 'GRM_off']
        off = h9 if h == 'H9' else off_primary
        counts = state.get('dosage_counts', {}).get(h)
        count = '/'.join(str(counts[str(i)]) for i in range(3)) if counts else 'quantitative, n=225'
        rows.append([h, count, f"{on['z_score']:.4f}", f"{off['z_score']:.4f}" + (' †' if h == 'H9' else ''),
                     f"{on['X_equivalent']:.4f}", f"{off['X_equivalent']:.4f}",
                     f"{on['empirical_P']:.4f}", f"{off['empirical_P']:.4f}",
                     'no/no', 'no/no'])
    text += markdown(['Locus', 'Dosage n0/n1/n2', 'Z on', 'Z off', 'X-equiv on', 'X-equiv off', 'P ancestry on', 'P ancestry off', '|X|>3.16 on/off', '|X|>5.45 on/off'], rows)
    text += '\n[observation] † H9 GRM-off **failed the primary 40-iteration fit**. The separately labelled same-model 400-iteration diagnostic converged at iteration 51 (Z=1.4125, X-equiv=0.9354, ancestry P=0.750; 199/199 diagnostic null fits succeeded). These diagnostic values do not erase the primary failure. H9 has only two alternate alleles (2/450), so numerical/inferential stability is poor. H9 GRM-on null had 165 successful and 34 failed fits out of 199 attempted; its P=0.7771 conditions on successful fits and may be biased. Every other primary converged row completed 199/199 null fits. Failed-null counts and quantiles are preserved; no failed fit is assigned Z=0.\n\n'
    text += '[observation] **No structure-absorption collapse occurs in this approximation.** H1^a retains 96.96% of its GRM-off |Z| (4.1163 versus 4.2456); H5 and H7 are unchanged, H3^r retains 99.66%, and continuous AMY1 is unchanged. H9 retains approximately 94.69% relative to its longer-iteration off diagnostic. The lack of collapse with this sparse GRM is not proof that the original scan\'s richer GRM would behave identically. GRM-off still estimates independent overdispersion tau2; it is not a no-overdispersion binomial GLM.\n\n'
    text += '### Extrapolated released-scale comparison (not measured Akbari statistics)\n\n'
    text += '[observation] The prespecified number is |X_LCT,released| × |Z_h,our| / |Z_LCT,our|, with all denominators from the same 225-person cohort/model. The table uses GRM-on; GRM-off is also saved in `statistics.tsv`. The continuous row is a Gaussian-trend analogy, not a binomial allele-scan result.\n\n'
    extrap = []
    for h in HAPS + ['AMY1_continuous']:
        r = by[h, 'GRM_on']
        sx = abs(state['references']['SLC45A2']['X']) * abs(r['z_score'] / by['SLC45A2', 'GRM_on']['z_score'])
        extrap.append([h, f"{r['AMY_to_LCT_abs_Z_ratio']:.3f}", f"{r['extrapolated_abs_X']:.3f}",
                       'yes' if r['extrapolated_abs_X'] > 3.16 else 'no',
                       'yes' if r['extrapolated_abs_X'] > 5.45 else 'no', f'{sx:.3f}'])
    text += markdown(['Locus', '|Z|/LCT |Z|', 'LCT-extrapolated |X|', '>3.16', '>5.45', 'Post-hoc SLC45A2 anchor |X|'], extrap)
    text += '\n[inference] **Uncertainty is substantial and not a numerical confidence interval.** A post-hoc alternative canonical anchor changes H1^a from 11.49 to 5.82 and H5 from 7.56 to 3.83, flipping H5\'s 5.45 classification. The loci do not share a validated universal multiplicative scaling factor. Small-cohort ascertainment, differing trajectories, locus-specific variances, approximation to PQL and sparse GRM all limit extrapolation. The ratio is internally comparable but not transportably calibrated; neither anchor estimates the paper\'s posterior pi in this cohort.\n\n'
    text += '### Permutation null: ancestry-conditioned versus unrestricted\n\n'
    text += '[observation] Primary shuffles permute individual dosage/phenotype values within the released **p2 ancestry component**, retaining each component\'s count/frequency, dates and K. Null statistics therefore need not center on zero: ancestry replacement correlated with date remains in the shuffled data. Continuous phenotype and both models are fitted anew, including variance components. An independent secondary unrestricted shuffle destroys that ancestry association. Both use seed 20261009 and 199 attempts, with plus-one two-sided empirical P. These are exchangeability nulls, not drift simulations.\n\n'
    nullrows = []
    for n in ['LCT', 'SLC45A2', 'AMY_real_SNP'] + HAPS + ['AMY1_continuous', 'background_chr3', 'background_chr4', 'background_chr6']:
        on = by[n, 'GRM_on']
        off = h9 if n == 'H9' else by[n, 'GRM_off']
        uon, uoff = unrestricted[n, 'GRM_on'], unrestricted[n, 'GRM_off']
        nullrows.append([n, f"{on['empirical_P']:.4f}", f"{off['empirical_P']:.4f}",
                         f"{uon['empirical_P']:.4f}", f"{uoff['empirical_P']:.4f}",
                         f"{uon['permutations_completed']}/{uoff['permutations_completed']}"])
    text += markdown(['Locus/control', 'Ancestry P on', 'Ancestry P off', 'Unrestricted P on', 'Unrestricted P off', 'Successful unrestricted on/off'], nullrows)
    text += '\n[inference] None of the five haplotypes has ancestry-conditioned empirical P<0.05. Canonical LCT also has ancestry-conditioned P=0.965/0.975: this null asks whether the trend exceeds what survives **within these broad ancestry labels**, not whether a known selected locus has an overall temporal trend. The null must not be advertised as recovering all established selection. GRM-on/off similarity and high conditional P can coexist: a covariance random effect and conditioning permutations on categories are different operations. No strong GRM absorption here does not establish ancestry-independent selection. Individual label exchangeability is not guaranteed under kinship or correlated drift, so these empirical P values are diagnostic rather than a demonstrated kinship-aware selection calibration; no multiple-testing correction is claimed.\n\n'
    text += '[observation] Full null distributions are `permutations.tsv`, `unrestricted_permutations.tsv` and the H9 diagnostic null TSV. `summary.json` and `unrestricted_null.json` give |Z| median/90th/95th/99th percentiles, completed counts and failures. The primary requested attempts all finished within the time budget; `complete_with_fit_failures` reflects numerical failures, not a timed-out null. The primary was rerun deterministically after correcting only that status label; biological statistics and nulls are independently compared to the preserved first-run snapshot.\n\n'
    text += '### Optional 600-person control cohort (never used in AMY ratios)\n\n'
    text += '[observation] A genotype-blind seeded sample of 600 unrelated dated ancients from Table S1 was fitted using the same marker list and LOCO construction. No synthetic AMY phenotypes are available for this cohort, so it supplies diagnostics only.\n\n'
    text += markdown(['Locus', 'n', 'Z on', 'Z off', 'X-equiv on', 'Released |X|'], [
        [n, 600, f"{next(r for r in diag['secondary_controls'] if r['locus']==n and r['model']=='GRM_on')['z_score']:.4f}",
         f"{next(r for r in diag['secondary_controls'] if r['locus']==n and r['model']=='GRM_off')['z_score']:.4f}",
         f"{next(r for r in diag['secondary_controls'] if r['locus']==n and r['model']=='GRM_on')['X_equivalent']:.4f}",
         f"{abs(state['references'][n]['X']):.4f}"] for n in ['LCT', 'SLC45A2']])
    text += '\n[observation] Absolute large-|X| values are still not reproduced, and LCT versus SLC45A2 ordering differs from the full released scan. The amended gate requires separation from the real AMY SNP, not reproducing canonical ordering against each other. No full ~16k-individual fit or original scan-equivalent cohort is claimed.\n\n'
    text += '### Frequencies in the primary 225 matched individuals\n\n'
    freq = pd.read_csv(OUT / 'frequencies.tsv', sep='\t')
    for dimension in ['ancestry', 'date_bin']:
        part = freq[freq.dimension.eq(dimension)]
        pivot = part.pivot(index='group', columns='haplotype', values='frequency')
        if dimension == 'date_bin':
            pivot = pivot.loc[sorted(pivot.index, key=lambda s: float(s.split(',')[0].lstrip('[')))]
        counts = part.groupby('group').n.first()
        text += '[observation] ' + dimension.replace('_', ' ').capitalize() + ' frequencies (alternate copies / 2n):\n\n'
        text += markdown(['Group', 'n'] + HAPS,
                         [[group, int(counts[group])] + [f'{pivot.loc[group,h]:.4f}' for h in HAPS] for group in pivot.index])
        text += '\n'
    text += '### Carried-forward facts verified from raw releases\n\n'
    text += '[observation] Mean AMY1 values refer to **533 copy-number ancients**, not merely the 288 haplotype-labelled ancients or 225 matched individuals. The seven quoted values reproduce at stated precision; CHG (n=1, mean 8.00) is an additional group omitted in the carried-forward list. The exact release label is “Iron age to early modern”, not “Iron age to modern”.\n\n'
    text += markdown(['Component', 'Copy-number n', 'Mean AMY1'], [[r.p2, int(r['count']), f"{r['mean']:.5f}"] for _, r in means.iterrows()])
    text += '\n[observation] Full 288-labelled ancient haplotype component frequencies, ordered oldest to youngest by mean date:\n\n'
    text += markdown(['Component', 'n', 'Mean BP', 'H1^a frequency', 'H5 frequency'],
                     [[r.component, int(r.n), f'{r.mean_age_BP:.1f}', f"{r['H1^a']:.4f}", f"{r['H5']:.4f}"] for _, r in groups.iterrows()])
    declines = int((np.diff(groups['H1^a']) < 0).sum())
    text += f"\n[observation] H1^a declines in {declines}/6 adjacent ancestry-component steps; H5 rises from {groups.iloc[0]['H5']:.4f} to {groups.iloc[-1]['H5']:.4f}. These are component averages, not replicated within-component trajectories.\n\n"
    pops = ['CHB', 'CHS', 'JPT', 'CDX', 'KHV', 'PEL', 'Pima', 'Quechua']
    picked = modern[modern.population.isin(pops)].set_index('population').loc[pops]
    text += '[observation] Selected modern-population checks (the broad East Asia range is evaluated on five major 1KG populations, not asserted for every tiny labelled subpopulation):\n\n'
    text += markdown(['Population', 'n', 'H5 frequency', 'H9 frequency'],
                     [[pop, int(r.n), f"{r['H5']:.4f}", f"{r['H9']:.4f}"] for pop, r in picked.iterrows()])
    ea = modern[modern.region.eq('EA')]
    text += f"\n[observation] Across *all* released EA-labelled populations, including n=1–2 groups, H5 spans {ea.H5.min():.3f}–{ea.H5.max():.3f} and H9 spans {ea.H9.min():.3f}–{ea.H9.max():.3f}; hence “0.40–0.50 / 0.00–0.03 in East Asia” is not literally true for every population. All per-population denominators are saved in the modern verification TSV.\n\n"
    text += f"[observation] Released AMY window: {state['released_AMY_window']['n_variants']:,} variants, maximum |X|={state['released_AMY_window']['max_abs_X']:.12f}, zero above 5.45, agreeing with the quoted values.\n\n"
    text += '### Explicit counterfactual answer and limitations\n\n'
    text += '[inference] **Under the directly fitted approximation, none of the five perfect-tag haplotype loci clears either transported cutoff** (largest |X-equivalent|=2.726 with GRM, 2.812 without). **Under the requested same-cohort LCT-anchored extrapolation, H1^a and H5 would clear 5.45** (11.493 and 7.556), while H7/H9/H3^r clear 3.16 but not 5.45. This second statement is an uncertain extrapolation, not a faithful rerun or a yes/no answer about the original Akbari scan. H1^a is 6.09x and H5 4.00x stronger in |Z| than the *released lead* AMY SNP in this cohort; the full 2,204-SNP window was not refitted to select a cohort-specific strongest tag. These contrasts are consistent with weak tagging, but are not proof that every SNP tag fails.\n\n'
    text += '[inference] The original-paper question remains **unresolved**: n=225 is small and biased toward later ancients; a perfect categorical haplotype tag is an assay upper bound, not a tested variant in the real scan; ancestry-conditioned nulls do not support independent temporal change at P<0.05; sparse GRM correction does not collapse the trend but cannot settle the behavior of the paper\'s full relatedness correction. Increased perfect-tag signal and compatibility with ancestry-associated frequency change coexist here; weak tagging versus structure-associated change is not definitively separated for the original scan. A counterfactual threshold exceedance is not a discovery, posterior selection probability or claim that Akbari actually assayed the copy-number allele.\n\n'
    text += '[observation] Reproducible entry points: `analysis/cn_as_snp.py`, `analysis/cn_as_snp_diagnostics.py`, `analysis/cn_as_snp_unrestricted.py`, and this report generator. All outputs are under allowed paths; raw releases unchanged. No git commands or staging were performed. The repository provenance harness was reused with only its git-metadata provider disabled under the explicit no-git rule; large raw genotype files have size/mtime identifiers rather than newly asserted SHA256 hashes.\n'
    target = 'docs/cn-as-snp-results.md'
    with open(target, 'a') as f:
        f.write(text)
    print('appended final marked section to', target)


if __name__ == '__main__':
    main()
