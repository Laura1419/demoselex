#!/usr/bin/env python3
"""Append full-panel diagnostics without interpreting outcomes after failed gates."""
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis import haplotype_background as h
from analysis.load_akbari import iter_chunks, DEFAULT_PATH


def table(frame):
    frame = frame.copy().fillna('NA')
    def fmt(x):
        return f'{x:.6g}' if isinstance(x, (float, np.floating)) else str(x)
    columns = list(frame.columns)
    return '| '+' | '.join(columns)+' |\n| '+' | '.join(['---']*len(columns))+' |\n'+''.join(
        '| '+' | '.join(fmt(x) for x in row)+' |\n' for row in frame.itertuples(index=False, name=None))


def facts():
    meta = pd.read_csv(h.OUT/'scratch/ancient_metadata_for_supplement.tsv', sep='\t')
    cn = pd.read_csv('data/raw/amy/copy_number_table_for_supplement.tsv', sep='\t')
    ancient = cn[cn.type.eq('Ancient')].merge(meta[['sample', 'age']], on='sample', validate='one_to_one')
    ancient = ancient[np.isfinite(ancient.age) & np.isfinite(ancient.AMY1)]
    means = ancient.groupby('p2').agg(n=('sample', 'size'), mean_AMY1=('AMY1', 'mean')).reset_index()
    means.to_csv(h.OUT/'verified_global_CN_means.tsv', sep='\t', index=False)
    olddir = Path('data/derived/amy_binarized')
    cohort = pd.read_csv(olddir/'primary_cohort.tsv', sep='\t')
    old = json.loads((olddir/'summary.json').read_text())
    labels = pd.Categorical(cohort.p2).codes
    tt = -cohort.age.to_numpy(float)/1000
    baseline = {}
    for key, cut in [('primary_le1', 1), ('sensitivity_le2', 2)]:
        c1 = cohort.H1.str.extract(r'^H(\d+)')[0].astype(float)
        c2 = cohort.H2.str.extract(r'^H(\d+)')[0].astype(float)
        y = ((c1 > cut).astype(int)+(c2 > cut).astype(int)).to_numpy()/2
        fit = h.temporal_fit(y, tt, labels)
        tr = tt.copy(); mu = y.copy()
        for c in np.unique(labels):
            tr[labels == c] -= tt[labels == c].mean()
            mu[labels == c] = y[labels == c].mean()
        centered = tt-tt.mean()
        between = float(centered@mu/(centered@centered))
        pooled = float(centered@y/(centered@centered))
        within_contribution = fit['beta_per_kyr']*float(tr@tr/(centered@centered))
        carriage = abs(between)/(abs(between)+abs(within_contribution))
        np.testing.assert_allclose(fit['beta_per_kyr'], old['channels'][key]['beta_within'], atol=1e-10)
        np.testing.assert_allclose(fit['MDE80_HC3_per_kyr'], old['channels'][key]['MDE80_HC3_per_kyr'], atol=1e-10)
        baseline[key] = dict(**fit, B_frequency=float(y.mean()), pooled_slope=pooled, carriage_share=carriage)
    primary = {r['model']: r for r in old['results'] if r['locus'] == 'primary_le1'}
    le2 = [r for r in old['results'] if r['locus'] == 'sensitivity_le2' and r['model'] == 'GRM_on'][0]
    baseline['saved_GLMM'] = dict(GRM_retention=abs(primary['GRM_on']['z_score']/primary['GRM_off']['z_score']),
                                   A_le2_X_equivalent=abs(le2['X_equivalent']))
    windows = {'AMY': (1, 103500000, 104500000), 'LCT': (2, 136300000, 136900000),
               'FADS': (11, 61300000, 61700000), 'SLC45A2': (5, 33800000, 34200000),
               'SLC24A5': (15, 48200000, 48600000)}
    h.initialize_intervals()
    for name, (chrom, lo, hi) in windows.items():
        transcripts = pd.read_csv(h.OUT/f'{name}_core_transcripts.tsv', sep='\t')
        assert transcripts.chrom.eq(chrom).all()
        assert (transcripts.start.ge(lo) & transcripts.end.le(hi)).all(), 'Released audit window misses focal/core genes'
    audit = {k: dict(n=0, max_abs_X=0., above_5_45=0) for k in windows}
    for chunk in iter_chunks(chunksize=500000, usecols=['CHROM', 'POS', 'X']):
        for name, (c, lo, hi) in windows.items():
            d = chunk[chunk.CHROM.eq(c) & chunk.POS.between(lo, hi)]
            if len(d):
                audit[name]['n'] += len(d)
                audit[name]['max_abs_X'] = max(audit[name]['max_abs_X'], float(d.X.abs().max()))
                audit[name]['above_5_45'] += int(d.X.abs().gt(5.45).sum())
    assert audit['AMY']['n'] == 2204
    np.testing.assert_allclose(audit['AMY']['max_abs_X'], 3.794956889, atol=1e-9)
    h.save_json(h.OUT/'verified_comparison.json', dict(dated_CN_n=len(ancient), baseline=baseline, akbari=audit))
    print('Fresh CN means, global temporal baseline and released scan reproduced.', flush=True)


def report():
    s = json.loads((h.OUT/'summary.json').read_text())
    assert s['status'].startswith('complete')
    f = json.loads((h.OUT/'verified_comparison.json').read_text())
    now = datetime.now(timezone.utc).isoformat()
    amy = s['loci']['AMY']
    rows = []
    for name, r in s['loci'].items():
        k = str(r['chosen_K']['1000000'])
        rows.append(dict(locus=name, small_K=r['chosen_K']['200000'], primary_K=int(k),
                         primary_silhouette=r['scores']['1000000'][k],
                         group_NMI=r['composition']['group_NMI'], country_NMI=r['composition']['country_NMI'],
                         batch_proxy_NMI=r['composition']['within_group_technology_NMI'],
                         gate_a=r['gates']['a']))
    initial = json.loads((h.OUT/'invalid_exclusion_run/summary.json').read_text())
    doc = '\n## NEW AUTHORITATIVE CORRECTED FULL-PANEL EXPERIMENT — '+now+'\n\n'
    doc += '[observation] Original registration persisted '+datetime.fromtimestamp(initial['prereg_mtime_ns']/1e9, timezone.utc).isoformat()+' before initial run '+initial['started_utc']+'. Original SHA256 '+initial['prereg_sha256']+' matches retained `prereg.original.md`. The clearly marked supervisor-approved POST-RESULTS coordinate correction persisted '+datetime.fromtimestamp(s['prereg_mtime_ns']/1e9, timezone.utc).isoformat()+' before CORRECTED run '+s['started_utc']+'. Amended prereg SHA256 '+s['prereg_sha256']+'. No K, encoding, width parameter or gate threshold was retuned. Initial outputs are retained under `invalid_exclusion_run/`; **their biological exclusion gate is INVALID**, irrespective of their numerical zero intersection. This section does not rewrite previous conclusions.\n\n'
    doc += '### Biological exclusion repair and generalised runtime annotation assertions\n\n'
    interval_rows=[]
    delta_rows=[]
    for name,r in s['loci'].items():
        a=r['interval_audit']; old=initial['loci'][name]
        before=pd.read_csv(h.OUT/f'invalid_exclusion_run/{name}_assignments.tsv',sep='\t')
        after=pd.read_csv(h.OUT/f'{name}_assignments.tsv',sep='\t')
        from sklearn.metrics import adjusted_rand_score
        oldk=old['chosen_K']['1000000'];newk=r['chosen_K']['1000000']
        ari=adjusted_rand_score(before[f'window_1000000_K{oldk}'],after[f'window_1000000_K{newk}'])
        interval_rows.append(dict(locus=name,original_exclusion=str(a['original_exclusion']),corrected_core=str(a['corrected_core']),
                                  corrected_exclusion=str(a['corrected_exclusion']),old_uncovered_focal_transcripts=a['original_focal_transcripts_uncovered'],
                                  corrected_uncovered_core_transcripts=a['uncovered_core_transcripts']))
        delta_rows.append(dict(locus=name,old_K=oldk,new_K=newk,old_new_assignment_ARI=ari,
                               old_features=old['coverage']['1000000']['features'],new_features=r['coverage']['1000000']['features'],
                               old_size_ARI=old['stability']['same_K_window_ARI'],new_size_ARI=r['stability']['same_K_window_ARI'],
                               old_group_NMI=old['composition']['group_NMI'],new_group_NMI=r['composition']['group_NMI']))
    pd.DataFrame(delta_rows).to_csv(h.OUT/'coordinate_correction_deltas.tsv',sep='\t',index=False)
    doc += '[observation] Immutable hg19 refGene identifies AMY1 copies through104301314, beyond the original104255000 exclusion. Generalised derivation includes EVERY transcript of ALL core-overlapping genes, rounds outward to1kb and applies the SAME25kb buffer to CORRECTED core ends. All focal/core transcripts are now covered; runtime assertions enforce this and every feature/window interval is checked against these annotation-derived bounds. Original control exclusions also miss full focal transcripts at LCT/FADS/SLC45A2; SLC24A5 originally covers its focal gene but the generalised core additionally covers neighbouring CTXN2/MYEF2.\n\n'+table(pd.DataFrame(interval_rows))+'\n'
    doc += '[observation] ALL genes in the expanded AMY PRIMARY outer window (103415000–104827000); genes not fully excluded are explicit allowed FLANK background genes, not silently treated as covered:\n\n'+table(pd.read_csv(h.OUT/'AMY_1000000_window_genes.tsv',sep='\t'))+'\n'
    doc += '[observation] Equivalent all-gene and per-transcript inventories for BOTH windows of ALL five loci are saved as `*_window_genes.tsv`, `*_window_gene_transcripts.tsv`, and `*_core_transcripts.tsv`. AMY core includes ACTG1P4, AMY1A/1B/1C, AMY2A/2B, LOC101928436 and RNPC3; zero core transcripts escape the exclusion, including the3-prime end of the last AMY1 copy. Annotation source is the existing immutable hg19_refGene.txt.gz, hashed in provenance. Canonical outer windows keep their original SNP centres and200kb/1Mb spans; AMY uses the unchanged200kb/1Mb total flanks on the corrected excluded region (outer spans612kb/1.412Mb).\n\n'
    doc += '### Initial INVALID versus corrected authoritative result delta\n\n'
    doc += '[observation] Changes are material to individual assignments, not a negligible repair; ARI quantifies the delta directly. K follows the SAME predeclared selection rule, not forced to match the contaminated run:\n\n'+table(pd.DataFrame(delta_rows))+'\n'
    doc += '[inference] The initial backgrounds demonstrably contained focal-gene sequence and were invalid. Material assignment changes after removing it are consistent with contamination sensitivity, although expanding exclusions also removes additional LD markers and (for AMY) moves outer flank boundaries; the delta does not isolate a causal contribution of that one gene copy. The corrected pass/fail verdict remains the criterion below; initial statistics are not biological evidence.\n\n'
    doc += '### Entire panel, two distinct encodings\n\n'
    doc += '[observation] **ALL 16,339 psam individuals** (15,836 ancient +503 modern EUR 1KG) participate in every locus/window/K fit; no curated/reference subset and no ancestry, relative or date filter. Pgen and streamed pvar both contain **52,382,872 variants**, agreeing with the operator. The reader docstring says 52,382,873: disagreement=-1 with that docstring, not with the requested number. 32,678 nominal chromosome slots are represented by 16,339 unphased dosage rows; **32,678 phased haplotypes were NOT observed or inferred**. Genotype windows include every eligible polymorphic biallelic SNP outside the exclusions, not a hand-picked AIM panel. Indels, monomorphic sites and excluded coordinates are counted in full window audits but are not background features.\n\n'
    inv = s['structural_inventory']
    doc += f'[observation] Structural encoding independently includes all {inv["archive_individual_rows"]} supplied diplotypes/{inv["labelled_slots"]} H1/H2 slots and {inv["labels"]} labels (including non-European archive participants). {inv["bundle_matched_individual_rows"]} rows/{inv["bundle_labelled_slots"]} slots overlap the full imputed panel:503 modern +225 dated ancients. The remaining15,611 panel rows lack measured structural labels; they still participate in genotype clustering, but no structural label is fabricated. Structural availability therefore has NOT increased the225-person dated outcome sample.\n\n'
    doc += '[observation] Whole-archive structural label frequencies and full-bundle overlap (unconditional INPUT inventory, not cluster-conditioned AMY results):\n\n'
    doc += table(pd.read_csv(h.OUT/'structural_label_inventory.tsv', sep='\t'))+'\n'
    doc += '[observation] All type/source/geography/component partitions are saved in `structural_partition_frequencies.tsv`. The separate published world-wide cluster-frequency table has6204 slots (ten listed classes), whereas this entire supplementary H1/H2 table has9160 slots/11 labels; these are different denominators and are not silently combined.\n\n'
    doc += '### Detection floor — BEFORE any outcome interpretation\n\n'
    floor_rows = [dict(variant=name, **value) for name, value in amy['outcome_blind_design_floor'].items()]
    doc += '[observation] Outcome-blind design calculation on225 archive-dated ancient carriers, using ONLY published baseline B frequency assumptions, the PRIMARY AMY assignments and dates:\n\n'+table(pd.DataFrame(floor_rows))+'\n'
    doc += '[inference] These reported frequency/kyr floors are **optimistic independent-binomial DESIGN floors on FAILED backgrounds**, not validated HC3 floors or achieved power. Nearly all date information is retained because these unstable clusters do little validated historical conditioning. No AMY outcome regression was fitted, so a background-conditioned HC3 floor is unavailable; presenting the optimistic decrease from the previous .146 floor as a power gain would be misleading. Structural-label-specific empirical frequencies/floors and conditional slopes are deliberately not tested after gate failure.\n\n'
    doc += '### Gate (a): published structure / batch composition — FAILED\n\n'
    doc += '[observation] Published AADR group/country metadata map6490/16,339=39.72% of panel rows; all other individuals remain assigned, with explicit UNMAPPED composition. All clusters have group, country, publication, sequencing data type, suffix technology, BP range and latitude/longitude summaries in `*_composition_*.tsv`. Weighted within-group technology NMI uses1347 persons in31 groups with>=20 persons and both technologies. Suffix and publication are batch PROXIES, not verified sequencing run IDs.\n\n'+table(pd.DataFrame(rows))+'\n'
    doc += '[inference] AMY and SLC45A2 primary silhouettes exceed.10, but every group NMI is below.10 and every country NMI below.05: known-structure gate FAILS at ALL five loci. Batch-proxy conditional NMI passes<=.20 but cannot rescue weak known structure or unmeasured true batch. Global labels are used only for validation, never training or panel selection.\n\n'
    doc += '[observation] AMY primary cluster composition (full group/country/date cross-tabs archived):\n\n'+table(pd.read_csv(h.OUT/'AMY_composition_summary.tsv', sep='\t'))+'\n'
    doc += '### K rule and gate (b): stability — FAILED\n\n'
    score_rows = []
    stable_rows = []
    for name, r in s['loci'].items():
        for w, scores in r['scores'].items():
            score_rows.append(dict(locus=name, window=w, **{'K'+k: v for k, v in scores.items()}, chosen_K=r['chosen_K'][w]))
        b = r['stability']
        stable_rows.append(dict(locus=name, same_K_size_ARI=b['same_K_window_ARI'],
                                selected_K_size_ARI=b['selected_K_window_ARI'],
                                neighbor_K_ARI='; '.join(k+':'+str(round(v,6)) for k,v in b['neighbor_K_ARI'].items()),
                                second_seed_ARI=b['second_seed_ARI'], exclusion_off_ARI=b['exclusion_off_ARI'], gate_b=r['gates']['b']))
    doc += '[observation] All K scores; maximal silhouette with .01 smaller-K tie rule, not outcome selection:\n\n'+table(pd.DataFrame(score_rows))+'\n'
    doc += '[observation] ARIs (required same-K size>=.65, each adjacent K>=.60, second seed>=.80):\n\n'+table(pd.DataFrame(stable_rows))+'\n'
    doc += '[inference] This fails the requirement to demonstrate stable backgrounds; ARI and the previous painter\'s L1/43% local agreement are different metrics and cannot be ranked or numerically equated. No retroactive choice of the smaller window or a convenient K is a valid rescue. Random projection, local LD variation and k-means fitting are possible causes, not established causes.\n\n'
    doc += '### Gate (c): exclusion and negative control — PASSED\n\n'
    excl_rows = []
    for name,r in s['loci'].items():
        excl_rows.append(dict(locus=name, excluded_sites=r['exclusion']['excluded_window_variants'],
                              forbidden_features=r['exclusion']['forbidden_primary_features'],
                              negative_eta_squared=r['negative_control']['eta_squared'],
                              negative_permutation_P=r['negative_control']['p'], gate_c=r['gates']['c']))
    tested_rows=[]
    for name in list(h.LOCI)[1:]:
        centre=(h.LOCI[name][1]+h.LOCI[name][2])//2
        audit=pd.read_csv(h.OUT/f'{name}_window_variants.tsv',sep='\t')
        hit=audit[audit.pos.eq(centre)]
        assert len(hit)>=1 and hit.excluded.all()
        tested_rows.extend(dict(locus=name,**row) for row in hit.to_dict('records'))
    tested=pd.DataFrame(tested_rows)
    tested.to_csv(h.OUT/'tested_variant_exclusion_audit.tsv',sep='\t',index=False)
    doc += '[observation] Exact tested SNPs and pgen offsets, explicitly excluded:\n\n'+table(tested)+'\n'
    doc += '[observation] Canonical tested SNP and every +/-25kb site beyond the corrected gene core are absent from both feature lists. AMY entire103915000–104327000 CORRECTED core-plus-buffer is absent, not merely its top SNP. Exact zero-intersection checks, included index hashes and per-coordinate audits are saved. Software regression changes ALL forbidden dosages to2-dosage, rebuilds features, and obtains bit-identical exclusion-on sketches; exclusion-off sketches change. Independent permuted balanced pseudo-variant (199 null permutations):\n\n'+table(pd.DataFrame(excl_rows))+'\n'
    doc += '[inference] Negative controls pass, but exclusion-off ARIs are often low: it is NOT true that removing the locus changes nothing. The mandatory exclusion prevents direct circularity, while long-range LD can remain. Background sharing would not prove independent origin, even with validated clusters: recombination, migration and unphased partner chromosomes can place a single origin in several diplotype clusters.\n\n'
    doc += '### Gate (d): usable imputed coverage — PASSED\n\n'
    cov_rows=[]
    for name,r in s['loci'].items():
        for w,c in r['coverage'].items(): cov_rows.append(dict(locus=name, window=w, **c))
    doc += '[observation] Every panel individual has100% observed dosage calls on every background feature, in BOTH windows for all FIVE loci. Unlike the failed raw-AIM lane, nobody is discarded for5–9% or<20% call rate. This is imputed input coverage, NOT direct ancient sequencing coverage or proof of imputation accuracy:\n\n'+table(pd.DataFrame(cov_rows))+'\n'
    doc += '### Structural sequence identity inventory (not a post-failure origin claim)\n\n'
    seq=inv['structure_distances']
    doc += f'[observation] Sequence archive: {seq["contigs"]} assembled haplotype contigs, {seq["unique_structures"]} exact structure strings; all94 join exactly to H labels through `haplotype_assignments.txt`. The distance table uses a different graph-block nomenclature: all31 representative strings match its hap1/hap2 fields;30/31 representative contig IDs link to named labels (one unresolved). The961 ordered distance pairs comprise31 zero diagonal pairs and930 positive pairs, maximum edit distance{seq["max_edit_distance"]:.6g}. Stable contig identifiers (before updated interval suffix) are used ONLY for the representative-to-label join; this join is archived.\n\n'
    doc += '[observation] Label/structure inventory, using supplied distances rather than treating all copies of a coarse label as identical:\n\n'+table(pd.read_csv(h.OUT/'structural_label_identity_inventory.tsv',sep='\t'))+'\n'
    doc += '[inference] These observations distinguish identical graph structure (e.g. all13 assembled H1^a contigs share one string across9 known population partitions) from coarse labels that pool different strings (H5:10 strings/16 contigs). Neither structural similarity nor heterogeneity proves recurrent independent origin; mutation phylogeny and phase-linked flanking genealogy remain untested. This INPUT check is not an AMY within-background recurrence result.\n\n'
    doc += '### STOP verdict and canonical positive controls\n\n'
    doc += '[observation] **AMY interpretation STOPPED: gates(a) and(b) fail; (c)/(d) pass.** No AMY background-conditioned frequency, time slope, P value, or origin-recurrence result is generated. Canonical LCT, FADS, SLC45A2 and SLC24A5 individually fail the SAME structure/stability gates, so their gated temporal/recurrence tests also stop. Positive-control recovery is0/4 validated loci, not four demonstrated null alleles. This machinery has not shown the expected within-background temporal changes or sharing at the controls; released global selection maxima are not a substitute for that validation.\n\n'
    doc += '[inference] Full-panel imputation fixes the demonstrated COVERAGE bottleneck, but does not by itself validate these local diplotype backgrounds. The requested AMY within-fixed-background question remains UNRESOLVED; this run does not overturn or newly corroborate the carriage-dominated global-label result.\n\n'
    doc += '### Independent comparison against established global-label results\n\n'
    doc += '[observation] Fresh raw CN+archive dates reconstruct533 dated ancients and the supplied component means (agree at stated rounding):\n\n'+table(pd.read_csv(h.OUT/'verified_global_CN_means.tsv',sep='\t'))+'\n'
    baserows=[]
    for cut in ['primary_le1','sensitivity_le2']:
        b=f['baseline'][cut]
        baserows.append(dict(cut=cut,n=b['n'],B_frequency=b['B_frequency'],HC3_MDE80_per_kyr=b['MDE80_HC3_per_kyr'],
                             within_slope=b['beta_per_kyr'],P=b['p'],carriage_share=b['carriage_share']))
    doc += '[observation] Reconstructed GLOBAL-label conditional slopes from the225-person cohort\'s actual H1/H2, with HC3 detection floor stated before slope interpretation (not using this run\'s clusters):\n\n'+table(pd.DataFrame(baserows))+'\n'
    glmm=f['baseline']['saved_GLMM']
    doc += f'[observation] Saved GLMM artifact check (not a refit): A<=2 |X-equivalent|={glmm["A_le2_X_equivalent"]:.6f}, GRM retention={100*glmm["GRM_retention"]:.6f}%. Thus3.64 clears3.16 but not5.45; ~97% retention and~94% carriage reproduce. The supplied .049/.146 floor belongs to PRIMARY A<=1, not A<=2; A<=2 instead has within .075969 and floor .152380. Both established within slopes are below their floors.\n\n'
    doc += '[observation] Fresh stream of immutable RELEASED scan, using original audit windows, not the new exclusion-on backgrounds:\n\n'+table(pd.DataFrame([dict(locus=k,**v) for k,v in f['akbari'].items()]))+'\n'
    doc += '[observation] AMY2204/max |X|3.794956889/zero>5.45 reproduces exactly; canonical values agree with16.52/16.15/13.15/5.99 at supplied precision. No supplied biological number disagrees. These SNP statistics do not directly test structural copy number.\n\n'
    doc += '### Genome-wide cost, provenance and limits\n\n'
    seconds=sum(r['elapsed_seconds'] for r in s['loci'].values())
    rss=max(r['peak_RSS_MiB'] for r in s['loci'].values())
    doc += f'[observation] Corrected five-locus computation after pvar cache: {seconds:.2f}s of locus work, total run {s["elapsed_seconds"]:.2f}s, peak process RSS {rss:.1f}MiB with one BLAS/OpenMP thread. Initial (invalid-exclusion) run including a full52.38M-row coordinate stream was87.72s; each locus reads~10–13k SNP features across every16,339 individual. Sketch memory is~3.0MiB per window; a256-variant read is~16MiB before working copies. Primary feature coordinates and all assignments are saved.\n\n'
    doc += f'[inference] Mechanical diplotype clustering is genome-wide-capable without phasing: measured average{seconds/5:.2f}s/locus suggests~{seconds/5*1000/3600:.2f} CPU-hours for1000 comparable loci, excluding provenance/reporting and assuming cached coordinates/linear IO. A naive52.38M-variant-specific run would cost~{seconds/5*52382872/(3600*24*365):.1f} CPU-years; window caching is necessary, and a separate +/-25kb exclusion is required for EVERY tested variant. This is an extrapolation, not a benchmark at that scale. Whole-panel phased haplotype analysis is NOT delivered. **Do not scale a biologically failed background model just because computation is cheap.**\n\n'
    doc += '[observation] Code: `analysis/haplotype_background.py`, `analysis/haplotype_background_report.py`; tests: `analysis/test_haplotype_background.py`; artifacts: `data/derived/haplotype_background/`. Commands use `.venv/bin/python`, `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`; biological main ran twice with the original invalid exclusions and twice with corrected exclusions, report `--facts` freshly streams released statistics, report without that flag APPENDS this section. All software tests and full-panel artifact checks must pass separately from the failed biological gates. RunRecorder from `scripts/provenance.py` disables git probes BEFORE construction, hashes sources/preregistration/small inputs/outputs, and records large immutable input size/mtime without asserting a content hash. No raw data were edited; no git command, staging or commit was run.\n'
    Path('docs/haplotype-background-results.md').open('a').write(doc)
    print('Appended results with fail-closed verdict.', flush=True)


def provenance(mode):
    from scripts import provenance as p
    p.get_git_info=lambda: {'omitted': 'operator forbids ALL git commands'}
    r=p.RunRecorder('haplotype-background-'+mode)
    r.set_seed(h.SEED)
    for f in [Path(__file__), Path(h.__file__), Path('analysis/test_haplotype_background.py'),h.PRE,
              h.ARCHIVE, h.ANNO, h.ANNOTATION, Path('data/raw/amy/copy_number_table_for_supplement.tsv'),
              Path('data/derived/amy_binarized/summary.json'),Path('data/derived/amy_binarized/primary_cohort.tsv')]:
        r.record_input(f)
    st=DEFAULT_PATH.stat(); r.record_parameter(str(DEFAULT_PATH),dict(size_bytes=st.st_size,mtime_ns=st.st_mtime_ns,hash='not computed: bounded IO'))
    for f in h.OUT.glob('*'):
        if f.is_file(): r.record_output(f)
    if Path('docs/haplotype-background-results.md').exists():r.record_output('docs/haplotype-background-results.md')
    r.save()


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--facts',action='store_true');args=ap.parse_args()
    if args.facts:facts()
    else:report()
    provenance('facts' if args.facts else 'report')
