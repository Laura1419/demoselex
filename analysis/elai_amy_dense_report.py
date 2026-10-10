#!/usr/bin/env python3
"""Audit dense-run artifacts and write the numeric resource/gate report; no fitting."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/derived/elai_amy_dense'


def main():
    sys.path.insert(0,str(ROOT))
    import scripts.provenance as provenance
    provenance.get_git_info=lambda: {'disabled':'operator prohibits all git commands'}
    recorder=provenance.RunRecorder('elai-amy-dense-report')
    recorder.set_seed(20261010)
    for p in [Path(__file__), ROOT/'analysis/test_elai_amy_dense.py', ROOT/'docs/elai-amy-dense-prereg.md', *OUT.glob('*.run.json'), OUT/'gate_metrics.json', OUT/'panel.json', OUT/'tests.log']:
        recorder.record_input(p)
    panel=json.loads((OUT/'panel.json').read_text())
    result=json.loads((OUT/'gate_metrics.json').read_text())
    runs=sorted([json.loads(p.read_text()) for p in OUT.glob('*.run.json')],key=lambda r:r['started_epoch'])
    cleanup=json.loads((OUT/'scratch_cleanup.json').read_text())
    with (OUT/'markers.tsv').open() as f:
        markers=list(csv.DictReader(f,delimiter='\t'))
    with (OUT/'pc_markers.tsv').open() as f:
        pcs=list(csv.DictReader(f,delimiter='\t'))
    with (OUT/'samples.tsv').open() as f:
        samples=list(csv.DictReader(f,delimiter='\t'))
    prereg=ROOT/'docs/elai-amy-dense-prereg.md'
    assert hashlib.sha256(prereg.read_bytes()).hexdigest()==panel['prereg_sha256']
    assert all(prereg.stat().st_mtime < r['started_epoch'] for r in runs)
    assert len(runs)<=3 and len(markers)==panel['retained_markers'] and len(pcs)==panel['pc_markers']
    assert len(samples)==panel['fitted_n']>=1000
    assert all(int(r['chrom'])==1 and not 103915000<=int(r['position'])<=104327000 for r in markers)
    assert all(float(r['maf'])>=.05 and float(r['call_rate'])>=.95 for r in markers)
    assert all(2<=int(r['chrom'])<=22 for r in pcs)
    assert all(int(b['position'])-int(a['position'])>=100 for a,b in zip(markers,markers[1:]))
    assert all(a['chrom']!=b['chrom'] or int(b['position'])-int(a['position'])>=100000 for a,b in zip(pcs,pcs[1:]))
    # snpinfo is written before EM: it can establish exclusion even if EM aborts.
    retained_files=list(OUT.glob('*.snpinfo.txt'))
    retained_count=0
    for path in retained_files:
        retained=[]
        for line in path.read_text().splitlines()[1:]:
            f=line.split()
            retained.append((f[0],int(f[1]),int(f[2])))
        assert retained==[(r['id'],int(r['position']),int(r['chrom'])) for r in markers]
        assert all(not 103915000<=pos<=104327000 for _,pos,_ in retained)
        retained_count+=len(retained)
    if retained_files:
        result['gates']['exclusion'].update(retained_core_markers=0,retained_markers_checked=retained_count,
                                           status='PASS',pass_gate=True)
        (OUT/'gate_metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    assert cleanup['removed'] and not Path(cleanup['path']).exists()
    assert result['verdict']=='STOP' or all(g['pass_gate'] for g in result['gates'].values())
    persistent=sum(p.stat().st_size for p in OUT.iterdir() if p.is_file())
    maximum=max(p.stat().st_size for p in OUT.iterdir() if p.is_file())
    assert persistent<8_000_000_000 and maximum<1_000_000_000
    disk=subprocess.check_output(['df','-h','/'],text=True)
    (OUT/'final_df.txt').write_text(disk)
    audit=dict(registration_hash_verified=True,registration_before_every_fit=True,fit_count=len(runs),
               selected_core_markers=0,actual_retained_markers_checked=retained_count,actual_retained_core_markers=0 if retained_files else None,
               pc_chr1_markers=0,marker_spacing_verified=True,pc_spacing_verified=True,
               scratch_removed=True,persistent_bytes=persistent,max_persistent_file_bytes=maximum,
               final_free_bytes=shutil.disk_usage('/').free,no_git_commands=True)
    (OUT/'independent_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    lines=['# Dense ELAI AMY revision — '+result['verdict'],'',
           '[observation] **Gated '+result['verdict']+'.** This is an instrument/resource result, not an AMY biological result. No allele-frequency analysis, local-state biological summary, effective independent-state estimate or biological claim was made. No cluster labels from incomplete/non-converged EM are accepted.','',
           '## Revised inputs and failure-mode rationale','',
           f"[observation] Immutable release: {panel['source_dimensions'][0]:,} variants × {panel['source_dimensions'][1]:,} samples. Eligible unique dated AADR-matched ancients: {panel['eligible_dated_unique']:,}; fitted selection **{panel['fitted_n']:,}**, **{panel['dated_n']:,} carrying dates** (seed {panel['seed']}). Baseline had128; no N128 fallback. Dates span {min(float(s['BP']) for s in samples):.1f}–{max(float(s['BP']) for s in samples):.1f} BP. Metadata labels never enter QC or EM.",'',
           f"[observation] Chr1 inclusive85,000,000–120,000,000, core103,915,000–104,327,000 excluded. **{panel['candidates']['local']:,} full-density biallelic SNP candidates**, **{panel['pre_thinning_qc']['local']:,} passing MAF≥0.05/call rate≥0.95 BEFORE thinning**, **{panel['retained_markers']:,} retained after100-bp greedy spacing**. Baseline's thin-before-QC gave1,164. Input size {panel['input_bytes']:,} bytes, deleted with scratch. Hardcalls are native ALT calls (not rounded dosages); complete imputed calls are not original aDNA coverage.",'',
           f"[observation] Independent PC panel: **{panel['pc_markers']:,} markers**, autosomes2–22 (chr1 excluded), ≥20-kb candidate spacing ({panel['candidates']['pc']:,} candidates), cohort QC ({panel['pre_thinning_qc']['pc']:,} passing), then≥100-kb retained spacing. Five eigen-PCs of standardized mean-imputed hardcall sample Gram matrix; all eigenvectors deterministic with single-thread BLAS. Baseline had145 chr2 SNPs. Per-chromosome counts: `{json.dumps(panel['pc_chromosome_counts'],sort_keys=True)}`. Builder wall {panel['build_seconds']:.3f}s.",'',
           '[inference] These prospective changes address sparse fitted LD coverage, the N128 pilot, non-convergence and inadequate independent structure validation. They were registered before fitting and motivated by diagnosed failure modes, not by revised gate values. No thresholds were relaxed.','',
           '## EM, live cost calibration and wall clocks','',
           '[observation] Fixed S=1,K=5,s=200 (exact quadratic EM), no early stop/checkpoint in this binary. Supervisor-approved live calibration uses the first five completed steps of the primary itself, not an extra fit/restart. Terminal convergence requires a completed fit and three consecutive nonnegative relative likelihood increments<1e-4. Full trace is below; labels from a partial trace are invalid.']
    for r in runs:
        lines += ['',f"[observation] **{r['prefix']}**: wall **{r['elapsed_seconds']:.3f}s**, exit{r['returncode']}, completed fit={r['complete']}, completed steps={r['convergence']['completed_steps']}/{r['convergence']['cap']}, convergence **{r['convergence']['status']}**, final relative increment={r['convergence']['final_relative_increment']}; first qualifying step={r['convergence']['first_qualifying_step']}. Reason: {r['reason']}.",'',
                  '```sh',' '.join(r['command']),'```','',
                  '| Step | Log likelihood | Reported increment | Relative increment | Elapsed s |','|---:|---:|---:|---:|---:|']
        for t in r['convergence']['trace']:
            rel='N/A (first step)' if t['relative_increment'] is None else f"{t['relative_increment']:.8g}"
            lines.append(f"| {t['step']} | {t['likelihood']:.3f} | {t.get('reported_increment')} | {rel} | {t.get('elapsed_seconds',0):.3f} |")
        cal=r['calibration']
        if cal:
            unit=cal['seconds_per_unit']; m=panel['retained_markers']; n=panel['fitted_n']
            feasible_n=int((1200-cal['initialization_seconds'])/(200*unit*m*25))
            feasible_m=int((1200-cal['initialization_seconds'])/(200*unit*n*25))
            lines += ['',f"[observation] Live calibration: {cal['seconds_per_step']:.3f}s/step (steps1–5 slope), initialization estimate {cal['initialization_seconds']:.3f}s; N×M×K²={cal['units']:,}; {unit:.8g}s/unit. Projected200-step EM time **{cal['projected_200_step_seconds']:.3f}s**, cap1200s; estimated affordable steps **{cal['affordable_steps']}** (export costs omitted, so optimistic). No configuration change/restart was made after this cost measurement.",'',
                      f"[inference] At this measured scaling, M={m:,},N={n:,},K=5,s=200 does not fit the cap. Keeping M and200 steps would permit approximately N≤{feasible_n}; keeping N and200 steps would permit M≤{feasible_m:,} (before export cost). Those are arithmetic resource projections, NOT evidence those settings converge or validate. User's marker/cohort targets were retained rather than silently reverting to the failed sparse setup."]
            lines += ['','| M | N | K | Projected seconds/step | Steps affordable in1200s (EM only) |','|---:|---:|---:|---:|---:|']
            for mm,nn,kk in ((5000,1000,5),(10000,500,5),(25000,200,5),(50000,100,5),(50000,1000,5),(50000,1000,3)):
                sec=unit*mm*nn*kk*kk
                lines.append(f'| {mm} | {nn} | {kk} | {sec:.3f} | {int(1200//sec)} |')
            baseline=.73/(128*1164*25)
            lines += ['',f'[observation] Parent corrected baseline arithmetic:128×1164×25=3,724,800 units,0.73/3,724,800={baseline:.8g}s/unit; M50,000,N1000,K5 is1.25e9 units (~244.98s/step at that baseline). Measured live slope is {unit/baseline:.3f}× that baseline, whose14.6s wall included export. Comparing EM slope to export-inclusive wall is only a cost heuristic.']
    lines += ['',f"[observation] Analysis fits run: **{len(runs)}/3 maximum**. SeedB and cross-K comparator NOT TESTED if absent below; no extra calibration fit, no frozen replay unless recorded. Primary-incomplete/non-converged STOP prevented spending on invalid labels.",'','## Numeric gates','',
              '| Gate | Numeric observations | Verdict |','|---|---|---|']
    g=result['gates']; known=g['known_structure']; stab=g['stability']; cov=g['coverage']; virt=g['virtual_site']; exc=g['exclusion']; em=g['em_convergence']
    lines += [f"| (a) Known structure | Independent PC markers={known['pc_markers']:,}; groups={known['groups']}, singleton groups={known['singleton_groups']}; PC silhouette={known['pc_silhouette']}, composition silhouette={known['composition_silhouette']}, group NMI={known['group_nmi']}, country NMI={known['country_nmi']}, singleton-robust country AMI={known['country_ami']}, date-bin NMI={known['date_bin_nmi']}, date/BP R²={known['date_composition_r2']} | **{known['status']}** |",
              f"| (b) Stability | Seed={stab['seed']['status']}, cross-K={stab['cross_k']['status']}, virtual argmax={stab['virtual_argmax']['status']}, thinning refit={stab['thinning_refit']['status']}; thinning posterior-subset diagnostic={stab['thinning_diagnostic']} | **{stab['status']}** |",
              f"| (c) Coverage | Usable={cov['usable_individuals']}/{cov['fitted_n']} ({cov['usable_fraction']:.1%}); dated selected={cov['dated_selected']}; dated COMPLETE --pkk export={cov['dated_complete_export']}; rows={cov['full_export_rows']}/{cov['expected_rows']:,} | **{cov['status']}** |",
              f"| (d) Virtual | Ordered sum error={virt['ordered_sum_error']}, q symmetry={virt['q_symmetry_error']}, triangle error={virt['q_triangle_sum_error']}; integration controls={virt['integration_controls_completed']}, linked/unlinked change={virt['bridge_max_change']} | **{virt['status']}** |",
              f"| (e) Exclusion | Builder core markers={exc['selected_core_markers']}; actual ELAI-retained core markers={exc['retained_core_markers']}; PC chr1 overlap={exc['pc_overlap']}; virtual104121000 inside core={exc['site_inside_core']}; flanks={exc['left_flank']}/{exc['right_flank']} | Builder **PASS**; full gate **{exc['status']}** |",
              f"| (additional) EM convergence | Completed steps={em['completed_steps']}/200; final relative increment={em.get('final_relative_increment')}; complete fit={bool(runs and runs[0]['complete'])}; tolerance1e-4 over three terminal increments | **{em['status']}** |",'',
              '[observation] The actual retained snpinfo was copied from the primary scratch output during EM, before cleanup, using `cp /tmp/elai-amy-dense-46qwl_7c/output/panel_k5_seed20261010.snpinfo.txt data/derived/elai_amy_dense/panel_k5_seed20261010.snpinfo.txt`. The report audit compares every retained ID/chromosome/coordinate with the selected panel. This is a pre-EM input/retention audit, not a posterior export.','',
              '[observation] `None` means unavailable, not zero or a measured failed silhouette. Known structure and its date-ordered/country-distance diagnostics cannot be computed without valid complete compositions; NOT TESTED does not imply absence of structure. Missing stability checks are NOT TESTED, overall PARTIAL, not falsely failed numerical comparisons. PARTIAL/NOT TESTED never pass the instrument gate. Coverage separates usable INPUT calls from complete EXPORT. ELAI writes snpinfo before EM: an independent exact ID/coordinate audit of that retained panel can PASS exclusion even with an incomplete EM/export; it does not authorize labels.','',
              '## Commands, provenance and resources','',
              '```sh','OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -u analysis/elai_amy_dense.py',
              'PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=analysis .venv/bin/python -m unittest analysis/test_elai_amy_dense.py analysis/test_elai_amy.py -v',
              'PYTHONDONTWRITEBYTECODE=1 .venv/bin/python analysis/elai_amy_dense_report.py','```','',
              f"[observation] Preregistration SHA256 `{panel['prereg_sha256']}` verified unchanged, mtime before every fit. Independent artifact audit verifies counts, QC, spacing, chr1/core exclusion and scratch deletion. 14 helper regressions passed (six new dense regressions plus eight reused baseline-helper checks, `tests.log`). Required `scripts/provenance.py` RunRecorder used with git provider disabled before construction; small inputs, binary and outputs hashed, huge immutable raw inputs recorded by size/mtime rather than rehashing. No fork rebuild, install, raw edit or git command. Exact binary commands retain now-deleted scratch paths; deterministic script recreates inputs from immutable raw files.",'',
              f"[observation] Scratch removed={cleanup['removed']}; runner wall={cleanup['elapsed_seconds']:.3f}s. Persistent summary output={persistent:,} bytes; largest persistent artifact={maximum:,} bytes. Largest scratch input={panel['input_bytes']:,} bytes (<1GB); full tensor never persisted. New total output far below8GB. Final free bytes={audit['final_free_bytes']:,}; `df -h /` before large input writes and after run:",'','```',disk.strip(),'```','',
              '[inference] **STOP**. Before accepting labels, a resource budget that permits dense markers and converged EM (including full-export cost) must be secured or a new explicitly authorized/preregistered cohort/compute trade-off made. Then independent structure, complete export, linkage control, seed/cross-K and genuinely refitted thinning stability must be measured. This resource/convergence failure is not a biological negative result and does not establish ELAI fails generally. No downstream allele-frequency or AMY-history claim is authorized.']
    (ROOT/'docs/elai-amy-dense-results.md').write_text('\n'.join(lines)+'\n')
    for p in [ROOT/'docs/elai-amy-dense-results.md',OUT/'gate_metrics.json',OUT/'independent_validation.json',OUT/'final_df.txt']:
        recorder.record_output(p)
    recorder.save()
    print(json.dumps(audit,indent=2))


if __name__=='__main__':
    main()
