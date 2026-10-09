#!/usr/bin/env python3
"""Append a new results section, never replace earlier conclusions."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from amy_binarized import OUT, CUTS


def table(frame, columns):
    rows=['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']
    for _,r in frame.iterrows():
        cells=[]
        for c in columns:
            v=r[c]
            cells.append(f'{v:.6g}' if isinstance(v,(float,np.floating)) else str(v))
        rows.append('| '+' | '.join(cells)+' |')
    return '\n'.join(rows)+'\n'


def main():
    s=json.loads((OUT/'summary.json').read_text()); d=json.loads((OUT/'distribution.json').read_text())
    f=pd.read_csv(OUT/'frequencies.tsv',sep='\t'); m=pd.read_csv(OUT/'modern_frequencies.tsv',sep='\t'); checks=pd.read_csv(OUT/'haplotype_copy_check.tsv',sep='\t')
    for label,cut in CUTS.items(): checks[label+'_allele']=np.where(checks.name_AMY1_copies.le(cut),'A','B')
    lines=['\n## NEW COMPLETED ANALYSIS — '+s['finished_utc']+'\n',
           '[observation] This section is appended, not a rewrite of earlier conclusions. `docs/amy-binarized-prereg.md` was written before results (local preregistration, not an external immutable registry). Primary data extraction, matching, controls, background-GRM reconstruction and both nulls completed. '+s['status']+'. Artifacts: `data/derived/amy_binarized/`; code: `analysis/amy_binarized.py`.\n',
           '### Model choice, A/B mapping and measured-copy check\n',
           '[observation] PRIMARY A = **H1^a, one haploid AMY1 copy**; B = every longer labelled haplotype. Sensitivity A<=2 = H1^a/H2A0/H2A2B2; A<=3 also adds H3^r/H3A2/H3A3B3. The initial H integer describes AMY1, not the later AMY2A/AMY2B suffixes. This cut is a **MODELING CHOICE**; the pooled B labels are not a single mutational origin. “Ancestral” is the paper framing, not a new ancestral-sequence reconstruction.\n',
           '[observation] Copy checks use '+str(s['copy_check']['n'])+' measured diplotypes, additive design rank '+str(s['copy_check']['rank'])+'/'+str(s['copy_check']['labels'])+', residual diploid RMSE '+f"{s['copy_check']['RMSE']:.4f}"+'. Homozygote columns give measured AMY1/2 (not confounded partner means); additive coefficients estimate per-haplotype copies jointly without an intercept.\n',
           table(checks,['haplotype','name_AMY1_copies','primary_le1_allele','sensitivity_le2_allele','sensitivity_le3_allele','homozygote_n','homozygote_mean_AMY1_per_hap','additive_estimated_AMY1_per_hap','diplotype_slots']),
           '[inference] The integer mapping agrees approximately, NOT exactly, with measured CN. Common H3^r/H5 estimates are near 3/5; H1^a homozygote mean is 1.1554, not exactly 1; sparse H3A2 and H9 estimates are 3.4802/9.4564 versus names 3/9. H2A2B2 and H3A2 lack homozygotes. Assignment/measurement errors and linked partner composition limit the check; no labels were recoded after seeing CN.\n',
           '[observation] Matched 225/288 labelled, finite-dated ancients; 63 absent, no invented baseline samples. Saved cohort order, dates and diplotypes exactly reconstructed. Genome-wide background markers '+str(s['GRM_verified']['background_markers'])+' across '+str(s['GRM_verified']['chromosomes'])+' autosomes, chr1-LOCO '+str(s['GRM_verified']['AMY_LOCO_markers'])+'. Cached chr1/2/5 matrices independently rebuilt from released background genotypes and matched to 1e-12 tolerance.\n',
           table(pd.DataFrame([dict(cut=k,**v) for k,v in s['dosage_counts'].items()]),['cut','n','n0','n1','n2','B_copies','alleles','B_frequency']),
           '### Detection floor — stated BEFORE interpreting trajectories\n',
           '[observation] Forward time t=-BP/1000; effects are B-frequency/kyr. Approximate two-sided alpha=.05, 80% power MDE=(t_.975,217+normal_.8)*SE. HC3 charges seven ancestry baselines; independent-binomial ideal floor uses observed component B frequencies and two copies/person. All cuts have within time information 49.3137 kyr² versus pooled 1281.3338 kyr²: group-centering retains only 3.85% of time information.\n',
           table(pd.DataFrame([dict(cut=k,**v) for k,v in s['channels'].items()]),['cut','B_frequency','MDE80_HC3_per_kyr','MDE80_ideal_binomial_per_kyr','MDE80_HC3_2kyr','observed_within_2kyr','date_span_kyr','MDE80_HC3_full_span']),
           '[inference] Every observed within slope is BELOW both approximate floors. Primary minimum detectable change is ~.146 B-frequency/kyr (~29.2 percentage points over 2kyr), versus observed ~.0486/kyr (~9.72 points/2kyr). The null within-history result is underpowered for the observed effect, not an exclusion of selection. Full-span MDE>1 is physically impossible under a linear frequency model and warns that a sustained full-span linear effect cannot be reliably identified with this within-component design. These Wald floors are approximate, not exact minima or kinship/drift-calibrated power; they do not say the large across-component trajectory is within-history.\n',
           '### Primary trajectory: ordered ancestry components, with all counts\n']
    for scope in ['labelled_288','matched_225']:
        p=f[f.cut.eq('primary_le1')&f.scope.eq(scope)&f.dimension.eq('component')]
        lines+=['[observation] '+scope+' (oldest to youngest by actual mean BP):\n',table(p,['group','mean_BP','n','n0','n1','n2','B_copies','alleles','B_frequency'])]
        diffs=np.diff(p.B_frequency)
        lines+=['[observation] B increases in '+str(int((diffs>0).sum()))+'/'+str(len(diffs))+' adjacent steps; starts '+f'{p.B_frequency.iloc[0]:.4f}'+' and ends '+f'{p.B_frequency.iloc[-1]:.4f}'+'. These are donor/component means, not seven independent temporal replicates.\n']
    for cut in ['sensitivity_le2','sensitivity_le3']:
        lines+=['[observation] '+cut+' matched ancestry frequencies, ordered by mean BP; full labelled counterparts in frequencies.tsv:\n',table(f[f.cut.eq(cut)&f.scope.eq('matched_225')&f.dimension.eq('component')],['group','mean_BP','n','n0','n1','n2','B_copies','alleles','B_frequency'])]
    lines+=['[observation] Matched date bins:\n',table(f[f.cut.eq('primary_le1')&f.scope.eq('matched_225')&f.dimension.eq('date_bin')],['group','n','n0','n1','n2','B_copies','alleles','B_frequency']),
            '### Same-cohort GLMM, reference thresholds and both anchors\n',
            '[observation] Controls re-fitted from integer ALT hardcalls; prior ordering gate passed and values reproduce. PQL/REML is a small-data approximation, not the original optimized Akbari scan. GRM-off retains independent overdispersion; tau1 when K=0 is inactive. Positive Z means older BP enrichment; rising B has negative Z.\n',
            table(pd.DataFrame(s['controls']),['locus','model','n','z_score','X_equivalent','released_X','converged']),
            table(pd.DataFrame(s['results']),['locus','model','z_score','X_equivalent','clears_3_16','clears_5_45','LCT_extrapolated_abs_X','SLC45A2_extrapolated_abs_X']),
            '[observation] Primary |Z| on/off=4.1163/4.2456, retention=96.96% (3.04% attenuation), X-equivalent=2.7261/2.8117: below both 3.16 and 5.45. A<=2 gives 3.6380, above 3.16 but below 5.45; A<=3 gives 2.4597. Both alternative cuts have essentially 100% GRM retention. Neither primary nor sensitivities exhibit a structure-absorption collapse.\n',
            '[inference] LCT-anchored values are 10.37–15.34 and all exceed 5.45; SLC45A2 gives 4.91–7.77, changing strict classification for A<=3. Anchor-dependent extrapolations are NOT measured released X, posterior probabilities, or evidence that Akbari would call this allele. Absolute control scale remains compressed (LCT |X-equiv|3.919 vs released16.523; SLC45A2 6.157 vs13.150).\n',
            '[observation] **Primary dosage is exactly 2 minus the previous H1^a dosage**. Its primary Z magnitude and power are therefore identical to that earlier binary perfect-tag locus: pooling “all longer” has NOT created extra information relative to the H1^a-versus-rest test. It is more common than testing rare H5/H9 separately, but the earlier H1^a lane already encoded this partition. The A<=2 cut is a genuinely different binary phenotype.\n',
            '### Permutation nulls (199 requested per row)\n',
            '[observation] Whole dosages shuffled within p2 or unrestricted, dates/K fixed, variance components refit; plus-one two-sided P conditional on successful fits. Quantiles below are median/90th/95th/99th of absolute Z; full rows in permutations.tsv.\n']
    null=pd.DataFrame(s['nulls']); null['quantiles']=null.null_abs_Z_quantiles.map(lambda a:', '.join(f'{v:.4f}' for v in a))
    lines+=[table(null,['null','cut','model','successful','failed','empirical_P','quantiles']),
            '[inference] Unrestricted null asks about total trend; the ancestry null preserves composition and asks a different question. The earlier, exactly reproduced LCT fit had ancestry P=.965 on/.975 off despite known selection: this conditioning is not a universal selection detector. Neither shuffle is a calibrated demographic/kinship-aware null. Small P unrestricted plus large P conditional is compatible with substantial carriage, not proof of selection or proof of its absence.\n',
            '### Two-channel split and the effect of correction\n',
            table(pd.DataFrame([dict(cut=k,**v) for k,v in s['channels'].items()]),['cut','beta_within','within_HC3_se','within_p','beta_between','between_bootstrap_se','beta_total','beta_pooled','within_pooled_contribution','carriage_absolute_share','channel_class']),
            '[observation] Primary within=.048575, between=.032824, total=.081399 B-frequency/kyr. The requested total is an INDEX: actual pooled=.034694=.001869 within contribution+.032824 carriage. Carriage accounts for **94.61%** of absolute pooled contributions (94.23% A<=2, 99.80% A<=3). No cut has opposing channel signs.\n',
            '[inference] **Mostly carriage/composition**, not demonstrated repeated within-history selection. A group-fixed-effect correction removes the between contribution by construction and leaves a nonsignificant within signal (primary P=.350; sensitivities .162/.963), all below detection floor. It therefore removes the statistically visible pooled trend in this operational sense, NOT 94.6% of the numerical conditional slope: the conditional coefficient increases relative to pooled while losing precision. **The genome-wide GRM random-effect correction does NOT remove the signal.** These two corrections must not be conflated. Broad period/ancestry labels are not exact histories; residual identifiability and drift calibration remain limited.\n',
            '### Modern binarized replication across continents\n',
            '[observation] Available populations and full cut sensitivity; n and all alternate denominators given. Every requested exact population label is present. Values for all released populations and dosage counts are in modern_frequencies.tsv.\n']
    targets=['PEL','Pima','Quechua','Maya','Surui','Han','CHB','CHS','JPT','CEU','TSI','IBS','GBR','YRI','GWD','ESN','Mbuti']
    wide=m[m.cut.eq('primary_le1')&m.population.isin(targets)].copy()
    for cut in ['sensitivity_le2','sensitivity_le3']:
        wide=wide.merge(m[m.cut.eq(cut)][['population','B_frequency']].rename(columns={'B_frequency':cut+'_B'}),on='population',validate='one_to_one')
    wide['order']=wide.population.map({x:i for i,x in enumerate(targets)}); wide=wide.sort_values('order')
    lines+=[table(wide,['population','n','B_copies','alleles','B_frequency','sensitivity_le2_B','sensitivity_le3_B']),
            '[observation] Primary PEL .9529, Pima/Quechua 1.0, Han .9848 and CHB/CHS/JPT .9238–.9660 exceed Surui .7500 and most Europe .8485–.9019. However Maya .8421 is lower than several European/African samples, YRI .9120/ESN .9091/GWD .8673 are already high, and **Mbuti is 1.0** (28/28 alleles).\n',
            '[inference] Descriptive high-B in several agricultural populations on separate continents is present, but **not a specific or universal agricultural replication under the primary cut**: the hunter-gatherer comparator Mbuti saturates too. The higher-copy cuts sharpen some contrasts (see table) without establishing independent selection origins. Shared labels, ancestry, small Pima/Quechua denominators and postulated starch categories preclude a causal diet claim; B-commonness may erase high-copy functional distinctions.\n',
            '### Cut sensitivity and explicit first-question verdict\n',
            '[inference] **YES, a pooled directional binary AMY trend is visible** under the prespecified unrestricted-P rule, consistently across cuts. **NO, the sparse LOCO GRM does not remove it.** **YES, ancestry fixed-effect/within-history conditioning removes statistical detectability, leaving mostly carriage and an underpowered within estimate.** Primary directly fitted X-equivalent does not clear either transported threshold; the <=2 sensitivity clears 3.16 only. Thus direction and correction conclusions are robust, but cutoff classification and modern dietary specificity depend on the MODELING CHOICE. This is not a faithful original-paper counterfactual or new selection discovery.\n',
            '## SECOND QUESTION — separate descriptive released-|X| distribution\n',
            '[observation] Independently streamed '+f"{d['n_finite']:,}"+' finite released variants, '+str(d['n_missing'])+' missing X; no centering/FILTER restriction. |X| quantiles:\n',
            table(pd.DataFrame({'quantile':d['quantile_probabilities'],'abs_X':d['quantiles']}),['quantile','abs_X']),
            '[observation] Fraction strictly exceeding each canonical exact window maximum (ties excluded), alongside operator-rounded thresholds:\n',
            table(pd.DataFrame([dict(locus=k,**v) for k,v in d['thresholds'].items()]),['locus','abs_X','n_exceed','fraction_exceed']),
            '[observation] LCT16.5232 versus AMY3.79496 is a **12.7282-X gap, 4.3540-fold magnitude**; only 2 variants exceed canonical LCT versus 167,199 (1.71669%) exceeding AMY. SLC24A5 is already near the strict threshold at5.9903. Genome-wide maximum is '+f"{d['quantiles'][-1]:.9f}"+' at '+d['top_variants'][0]['ID']+', so the operator’s claim that LCT is the strongest signal is not literally true at the released single-variant level: this chr5 variant is larger than canonical LCT16.52. Top20 IDs are in distribution.json.\n',
            table(pd.DataFrame([dict(category='3.16<=|X|<=5.45',variants=d['band_variants'],variant_fraction=d['band_variants']/d['n_finite'],any_blocks=d['band_any_blocks'],max_class_blocks=d['band_max_blocks']),dict(category='|X|>5.45',variants=d['strong_variants'],variant_fraction=d['strong_variants']/d['n_finite'],any_blocks=d['strong_any_blocks'],max_class_blocks=d['strong_max_blocks'])]),['category','variants','variant_fraction','any_blocks','max_class_blocks']),
            '[observation] 250kb blocks use (chr,floor(POS/250000)), '+str(d['all_blocks'])+' represented. Any-band and any-strong sets overlap in '+str(d['overlap_blocks'])+' blocks; mutually exclusive block-MAX categories are '+str(d['band_max_blocks'])+' band-only versus '+str(d['strong_max_blocks'])+' strong. These are not independent LD units.\n',
            '[inference] The strict threshold keeps **31,939 variants in 723 blocks**, including all listed canonical windows except AMY; it discards **356,335 band variants**, including AMY’s best SNP and variants in **9,285 additional band-only blocks**. The band contains ~11.16 times as many variants as the retained strong set. There is a broad prioritization domain below the strict gate, not just the handful of spectacular LCT/FADS signals. **False-positive caveat: the discarded band is NOT all real selection**; released FDR/inflation, demography, LD and correlated tests make these follow-up candidates, not discoveries. Magnitude proximity is descriptive, not shared mechanism or causal evidence.\n',
            '### Carried-forward verification and disagreements\n',
            '[observation] Measured CN means/slopes rechecked on the 533 dated CN ancients, NOT the 225 matched or 288 labelled cohort:\n',
            table(pd.read_csv(OUT/'verified_CN.tsv',sep='\t'),['component','n','mean_AMY1','beta','p']),
            '[observation] All identified CN within-component slopes are nonsignificant, not proved zero; CHG n=1 is unidentifiable. H1^a full labelled component frequency .7667→.1283, declining 5/6 steps, and H5 .0000→.2261 reproduce; verification table follows.\n',
            table(pd.read_csv(OUT/'verified_ancient_haplotypes.tsv',sep='\t'),['component','n','mean_BP','H1a_frequency','H5_frequency']),
            table(wide,['population','n','H5_frequency','H9_frequency']),
            '[observation] Disagreements with supplied modern approximations: PEL H5=.44706, NOT .46 (rounds .45); Pima=.71429, Quechua=.50. CHB H5=.39806 and KHV H9=.03030 marginally differ from .40/.03 at precision; small EA groups are not covered by those broad ranges. PEL/Pima/Quechua H9=.11765/.10714/.16667 agree with .11–.17 after rounding. AMY released window remains 2,204 SNPs, max3.794956889, zero >5.45. Full released band356,335 is six more than prior two-channel scored moderate356,329: different sampled panel, not an inconsistency in the released counts. Full represented blocks10,672 versus prior scored panel10,634 likewise differ in ascertainment. Canonical rounded values agree at stated two decimals; exact and rounded exceedance counts differ, as explicitly tabulated.\n',
            '### Validation, provenance and limits\n',
            '[observation] Source-named tests audit dosage/complement identity, all count denominators, three cut definitions, within/pooled channel algebra, control reproduction, preregistration-before-output timestamp, null rows/P values, and genome-wide block totals. Commands and validation log are archived. Provenance manifests under analysis/runs/ use scripts/provenance.py with ONLY git metadata probes disabled. Raw inputs untouched; no git commands or staging performed. Large genotype/stat files recorded by path/size/mtime, not falsely asserted content hashes.\n',
            '[inference] Sparse GRM, n=225 late-biased sampling (174/225 below2kyr), label coarseness, copy/genotype uncertainty, approximate Wald power, uncalibrated shuffle/drift nulls and anchor-dependent scale limit interpretation. Directional pooled frequency change is measured; independent within-history selection and what the full original Akbari GRM would do remain unresolved.\n']
    with Path('docs/amy-binarized-results.md').open('a') as out: out.write('\n'.join(lines))
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    from scripts import provenance
    provenance.get_git_info=lambda: {'commit':'not inspected: no-git rule','branch':'not inspected','dirty':'unknown'}
    recorder=provenance.RunRecorder('amy-binarized-report',command=[sys.executable]+sys.argv)
    for path in [Path(__file__),OUT/'summary.json',OUT/'distribution.json',OUT/'frequencies.tsv',OUT/'modern_frequencies.tsv']:
        recorder.record_input(path)
    recorder.record_output('docs/amy-binarized-results.md')
    recorder.save()


if __name__=='__main__': main()
