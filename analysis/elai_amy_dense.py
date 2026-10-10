#!/usr/bin/env python3
"""QC-then-thin dense AMY instrument run; live resource calibration, fail closed."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time

import numpy as np
from sklearn.metrics import adjusted_mutual_info_score, normalized_mutual_info_score
from elai_amy_build import cohort, excluded, SITE
from elai_amy_gate import BINARY, read_joint, joint_errors, silhouette, same_k, unlinked_parameters
from load_akbari_geno import AkbariGeno, iter_pvar, paths

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/elai_amy_dense'
PREREG = ROOT / 'docs/elai-amy-dense-prereg.md'
SEED = 20261010
REGION = (85000000, 120000000)
N = 1000
STEPS = 200
TOL = 1e-4


def save(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def select_samples():
    # Reuse the source-named baseline matcher, not its bounded random subset.
    import elai_amy_build as old
    original_rng = old.np.random.default_rng
    class AllEligible:
        def choice(self, size, *_args, **_kwargs):
            return np.arange(size)
    old.np.random.default_rng = lambda _seed: AllEligible()
    try:
        eligible, count = cohort()
    finally:
        old.np.random.default_rng = original_rng
    chosen = np.sort(original_rng(SEED).choice(count, min(N, count), replace=False))
    return [eligible[i] for i in chosen], count


def qc_mask(h):
    valid = (h >= 0) & (h <= 2)
    calls = valid.mean(1)
    af = np.divide(np.where(valid, h, 0).sum(1), 2 * valid.sum(1),
                   out=np.zeros(len(h), dtype=float), where=valid.sum(1) > 0)
    maf = np.minimum(af, 1-af)
    return (calls >= .95) & (maf >= .05), maf, calls


def retain_spaced(rows, h, mask, last, spacing):
    """Spacing state changes only AFTER passing QC, never on failed candidates."""
    kept = []
    for j in np.flatnonzero(mask):
        r = rows[j]
        if r[2] - last.get(r[1], -spacing) >= spacing:
            kept.append(j)
            last[r[1]] = r[2]
    return kept


def build(work):
    started = time.time()
    samples, eligible = select_samples()
    with (OUT / 'samples.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(samples[0]), delimiter='\t')
        w.writeheader(); w.writerows(samples)
    local, pcs, ph = [], [], []
    last_local, last_pc, last_candidate = {}, {}, {}
    candidate_counts = {'local': 0, 'pc': 0}
    qc_counts = {'local': 0, 'pc': 0}
    coverage = np.zeros(len(samples), dtype=np.int64)
    geno = work / 'genotypes.txt'
    batch = []
    def process(reader, f):
        nonlocal batch
        if not batch:
            return
        h = reader.read_hardcalls(np.array([r[0] for r in batch], dtype=np.uint32))
        mask, maf, calls = qc_mask(h)
        for kind in ('local', 'pc'):
            ix = [j for j, r in enumerate(batch) if (r[1] == 1) == (kind == 'local')]
            rows = [batch[j] for j in ix]
            submask = mask[ix]
            qc_counts[kind] += int(submask.sum())
            indices = retain_spaced(rows, h[ix], submask,
                                    last_local if kind == 'local' else last_pc, 100 if kind == 'local' else 100000)
            for j in indices:
                r = rows[j]; g = h[ix[j]]
                if kind == 'pc':
                    pcs.append(r); ph.append(g.copy())
                else:
                    local.append((*r, float(maf[ix[j]]), float(calls[ix[j]])))
                    valid = (g >= 0) & (g <= 2)
                    coverage[:] += valid
                    codes = np.array([r[4]+r[4], r[4]+r[5], r[5]+r[5], '??'])
                    f.write(r[3] + ' ' + ' '.join(codes[np.where(valid, g, 3)]) + '\n')
        batch = []
    subprocess.run(['df', '-h', '/'], check=True)
    with AkbariGeno(sample_subset=np.array([s['index'] for s in samples], dtype=np.uint32)) as reader, geno.open('w') as f:
        dimensions = [reader.n_variants, reader._raw_n_samples]
        for i, (chrom, pos, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
            if not chrom.isdigit() or not 1 <= int(chrom) <= 22:
                continue
            c = int(chrom)
            if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
                continue
            if c == 1:
                if not REGION[0] <= pos <= REGION[1] or excluded(pos):
                    continue
                candidate_counts['local'] += 1
            else:
                if pos-last_candidate.get(c, -20000) < 20000:
                    continue
                last_candidate[c] = pos
                candidate_counts['pc'] += 1
            batch.append((i, c, pos, vid, ref, alt))
            if len(batch) == 1024:
                process(reader, f)
        process(reader, f)
    m = len(local)
    with (work / 'panel.inp').open('w') as f, geno.open() as g:
        f.write(f'{len(samples)}\n{m}\nIND ' + ' '.join(s['sample'] for s in samples) + '\n')
        shutil.copyfileobj(g, f, length=1024*1024)
    geno.unlink()
    with (work / 'panel.pos').open('w') as f:
        for r in local:
            f.write(f'{r[3]} {r[2]} {r[1]}\n')
    with (OUT / 'markers.tsv').open('w') as f:
        f.write('pvar_index\tchrom\tposition\tid\tref\talt\tmaf\tcall_rate\n')
        for r in local:
            f.write('\t'.join(map(str, r))+'\n')
    with (OUT / 'pc_markers.tsv').open('w') as f:
        f.write('pvar_index\tchrom\tposition\tid\tref\talt\n')
        for r in pcs:
            f.write('\t'.join(map(str, r))+'\n')
    x = np.array(ph, dtype=float).T
    x[(x < 0) | (x > 2)] = np.nan
    mean = np.nanmean(x, axis=0)
    x = np.where(np.isfinite(x), x, mean)
    x = (x-mean)/np.sqrt(2*(mean/2)*(1-mean/2))
    values, vectors = np.linalg.eigh(x @ x.T)
    pc = vectors[:, -5:][:, ::-1] * np.sqrt(np.maximum(values[-5:][::-1], 0))
    np.savetxt(OUT / 'independent_pcs.tsv', pc, delimiter='\t')
    pos = np.array([r[2] for r in local])
    assert all(not excluded(p) for p in pos) and excluded(SITE)
    assert (pos < SITE).any() and (pos > SITE).any()
    manifest = dict(seed=SEED, source_dimensions=dimensions, region=list(REGION), eligible_dated_unique=eligible,
                    fitted_n=len(samples), dated_n=len(samples), candidates=candidate_counts, pre_thinning_qc=qc_counts,
                    retained_markers=m, local_spacing_bp=100, pc_markers=len(pcs), pc_candidate_spacing_bp=20000,
                    pc_retained_spacing_bp=100000, pc_chromosome_counts={str(c):sum(r[1]==c for r in pcs) for c in range(2,23)},
                    individual_call_rates=(coverage/m).tolist(), usable_individuals=int((coverage/m >= .95).sum()),
                    excluded_markers=0, pc_overlap=0, left_flank=int(pos[pos<SITE][-1]), right_flank=int(pos[pos>SITE][0]),
                    build_seconds=time.time()-started, prereg_sha256=hashlib.sha256(PREREG.read_bytes()).hexdigest(),
                    prereg_mtime_epoch=PREREG.stat().st_mtime, input_bytes=(work/'panel.inp').stat().st_size)
    save('panel.json', manifest)
    print('BUILD', m, 'local markers', len(pcs), 'PC markers', len(samples), 'samples', flush=True)
    return manifest, samples, local, pc


def convergence(trace, complete):
    rel = [None] + [(b['likelihood']-a['likelihood'])/max(1, abs(a['likelihood'])) for a,b in zip(trace, trace[1:])] if trace else []
    for t, r in zip(trace, rel):
        t['relative_increment'] = r
    hits = [i+1 for i in range(3, len(rel)) if all(r is not None and 0 <= r < TOL for r in rel[i-2:i+1])]
    passed = bool(complete and hits and hits[-1] == len(trace))
    return dict(status='PASS' if passed else 'FAIL', pass_gate=passed, tolerance=TOL,
                completed_steps=len(trace), cap=STEPS, final_relative_increment=rel[-1] if rel else None,
                first_qualifying_step=hits[0] if hits else None, trace=trace)


def consume(fifo, n, k, markers, result):
    comp = np.zeros((n,k)); thin = comp.copy()
    count = 0; norm = 0.; err = 0.
    try:
        with fifo.open() as f:
            assert f.readline().startswith('individual_index\t')
            for line in f:
                fields = line.split('\t')
                i, a, b = int(fields[0]), int(fields[5]), int(fields[6])
                m = (count//(k*k)) % len(markers)
                assert i == count//(len(markers)*k*k)
                assert a == count % (k*k)//k and b == count % k
                assert fields[2] == markers[m][3]
                p = float(fields[7]); q = float(fields[8])
                assert np.isfinite(p) and np.isfinite(q) and 0 <= p <= 1 and 0 <= q <= 1
                comp[i,a] += p
                if m % 2 == 0:
                    thin[i,a] += p
                norm += p; count += 1
                if count % (k*k) == 0:
                    err = max(err, abs(norm-1)); norm = 0.
        result.update(rows=count, complete=count == n*len(markers)*k*k, norm_error=err,
                      composition=(comp/len(markers)).tolist(), thin_composition=(thin/len(markers[::2])).tolist())
    except Exception as exc:
        result.update(rows=count, complete=False, error=str(exc))


def consume_marginal(fifo, n, k, markers, result):
    """Drain --pk1 lower-state copy counts, retaining only mean compositions."""
    comp, thin = [], []
    error = 0.
    try:
        with fifo.open() as f:
            for line in f:
                a = np.fromstring(line, sep=' ')
                assert a.size == len(markers)*k, 'incomplete marginal individual'
                a = a.reshape(len(markers), k)/2
                assert np.isfinite(a).all() and (a >= 0).all() and (a <= 1).all()
                error = max(error, float(np.abs(a.sum(1)-1).max()))
                comp.append(a.mean(0).tolist())
                thin.append(a[::2].mean(0).tolist())
        result.update(rows=len(comp)*len(markers)*k, complete=len(comp)==n,
                      individuals=len(comp), norm_error=error,
                      composition=comp, thin_composition=thin)
    except Exception as exc:
        result.update(complete=False, error=str(exc), individuals=len(comp))


def fit(work, panel, markers, k, seed, wall, calibration=False, full_export=True):
    prefix = f'panel_k{k}_seed{seed}'
    fifo = work/'output'/(prefix + ('.pkk.tsv' if full_export else '.pk1.txt'))
    os.mkfifo(fifo)
    streamed = {}
    consumer = consume if full_export else consume_marginal
    thread = threading.Thread(target=consumer, args=(fifo, panel['fitted_n'], k, markers, streamed), daemon=True)
    thread.start()
    cmd = [BINARY, '-g', str(work/'panel.inp'), '-p', '1', '-pos', str(work/'panel.pos'),
           '-C', '1', '-c', str(k), '-w', '0', '-s', str(STEPS), '-mg', '10', '-R', str(seed),
           '-o', prefix, '--pkk-site', '1:104121000', '--pkk' if full_export else '--pk1']
    start = time.time(); trace = []; reason = None; cal = None
    assert PREREG.stat().st_mtime < start
    with (OUT/f'{prefix}.stdout.txt').open('w') as log:
        proc = subprocess.Popen(cmd, cwd=work, stdout=log, stderr=subprocess.STDOUT)
        while proc.poll() is None:
            time.sleep(.5)
            text = (OUT/f'{prefix}.stdout.txt').read_text()
            trace = [dict(step=int(m[0])+1, likelihood=float(m[1]), reported_increment=float(m[2]))
                     for m in re.findall(r'### Exact update\s+0\s+(\d+)\s+([-\d.e+]+)\s+([-\d.e+]+)', text)]
            elapsed = time.time()-start
            for t in trace:
                # Retain timestamp when each step is first observed.
                if t['step'] not in seen_steps:
                    seen_steps[t['step']] = elapsed
                t['elapsed_seconds'] = seen_steps[t['step']]
            if calibration and len(trace) >= 5 and cal is None:
                per_step = (trace[4]['elapsed_seconds']-trace[0]['elapsed_seconds'])/4
                init = max(0., trace[0]['elapsed_seconds']-per_step)
                projection = init + STEPS*per_step
                cal = dict(completed_steps=5, seconds_per_step=per_step, initialization_seconds=init,
                           projected_200_step_seconds=projection, budget_seconds=wall,
                           units=panel['fitted_n']*len(markers)*k*k,
                           seconds_per_unit=per_step/(panel['fitted_n']*len(markers)*k*k),
                           affordable_steps=int(max(0, wall-init)//per_step))
                save('live_calibration.json', cal)
                if projection > wall:
                    reason = 'RESOURCE STOP: live five-step projection exceeds primary wall cap'
            if elapsed > wall:
                reason = 'RESOURCE STOP: hard wall cap reached'
            if reason:
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill(); proc.wait()
                break
    elapsed = time.time()-start
    complete = proc.returncode == 1 and '-elai: finished' in (OUT/f'{prefix}.stdout.txt').read_text()
    if not complete:
        # Unblock an exporter FIFO that was never opened by an aborted EM.
        try:
            fd = os.open(fifo, os.O_WRONLY | os.O_NONBLOCK); os.close(fd)
        except OSError:
            pass
    thread.join(timeout=5)
    em = convergence(trace, complete)
    result = dict(prefix=prefix, command=cmd, k=k, seed=seed, started_epoch=start, elapsed_seconds=elapsed,
                  returncode=proc.returncode, complete=complete, reason=reason, convergence=em,
                  full_export_rows=streamed.get('rows',0) if full_export else 0,
                  full_export_complete=streamed.get('complete',False) if full_export else False,
                  full_export_norm_error=streamed.get('norm_error') if full_export else None, calibration=cal)
    if not full_export:
        result.update(export_mode='marginal', marginal_export=streamed.copy())
        result['marginal_export'].pop('composition', None)
        result['marginal_export'].pop('thin_composition', None)
    if complete and streamed.get('complete'):
        np.savetxt(OUT/f'{prefix}.composition.tsv', streamed['composition'], delimiter='\t')
        np.savetxt(OUT/f'{prefix}.thin_composition.tsv', streamed['thin_composition'], delimiter='\t')
        for suffix in ('.pkk-site.tsv', '.snpinfo.txt'):
            shutil.copyfile(work/'output'/(prefix+suffix), OUT/(prefix+suffix))
        p,q = read_joint(OUT/(prefix+'.pkk-site.tsv'), panel['fitted_n'],k)
        result['joint_errors'] = joint_errors(p,q)
        result['excluded_retained_markers'] = sum(excluded(int(f[1])) for l in (OUT/(prefix+'.snpinfo.txt')).read_text().splitlines() if len(f:=l.split())>=6 and f[1].isdigit())
    save(prefix+'.run.json',result)
    return result


seen_steps = {}


def metrics(panel, samples, pc, runs):
    r = runs[0] if runs else None
    valid = bool(r and r['complete'] and r['full_export_complete'])
    missing = dict(status='NOT TESTED', pass_gate=False)
    known = dict(missing, pc_markers=panel['pc_markers'], pc_silhouette=None, composition_silhouette=None,
                 group_nmi=None, country_nmi=None, country_ami=None, date_bin_nmi=None,
                 date_composition_r2=None, groups=len(set(s['group'] for s in samples)),
                 singleton_groups=sum(sum(t['group']==g for t in samples)==1 for g in set(s['group'] for s in samples)))
    stability = dict(status='PARTIAL', pass_gate=False, seed=dict(missing), cross_k=dict(missing),
                     virtual_argmax=dict(missing), thinning_refit=dict(missing), thinning_diagnostic=None)
    virtual = dict(missing, ordered_sum_error=None, q_symmetry_error=None, q_triangle_sum_error=None,
                   integration_controls_completed=0, bridge_max_change=None)
    if valid:
        a = np.loadtxt(OUT/(r['prefix']+'.composition.tsv')); labels=a.argmax(1)
        countries=[s['country'] for s in samples]; groups=[s['group'] for s in samples]
        date=np.array([s['BP'] for s in samples]); design=np.column_stack((np.ones(len(a)),a))
        fitted=design @ np.linalg.lstsq(design,date,rcond=None)[0]
        known.update(pc_silhouette=silhouette(pc,labels), composition_silhouette=silhouette(a,labels),
                     nonempty_assignments=len(set(labels)), group_nmi=float(normalized_mutual_info_score(groups,labels)),
                     country_nmi=float(normalized_mutual_info_score(countries,labels)), country_ami=float(adjusted_mutual_info_score(countries,labels)),
                     date_bin_nmi=float(normalized_mutual_info_score(np.floor(date/2000),labels)),
                     date_composition_r2=float(1-((date-fitted)**2).sum()/((date-date.mean())**2).sum()))
        known['pass_gate']=all(known[x] is not None and known[x]>=v for x,v in [('pc_silhouette',.1),('composition_silhouette',.1),('group_nmi',.1),('country_nmi',.05)])
        known['status']='PASS' if known['pass_gate'] else 'FAIL'
        p,_=read_joint(OUT/(r['prefix']+'.pkk-site.tsv'),len(samples),5)
        thin=np.loadtxt(OUT/(r['prefix']+'.thin_composition.tsv'))
        from sklearn.metrics import adjusted_rand_score
        stability['thinning_diagnostic']=dict(ari=float(adjusted_rand_score(labels,thin.argmax(1))), mean_l1=float(np.abs(a-thin).sum(1).mean()), independent_refit=False)
        for other in runs[1:]:
            if not other['complete'] or not other['full_export_complete']:
                continue
            b=np.loadtxt(OUT/(other['prefix']+'.composition.tsv'))
            pb,_=read_joint(OUT/(other['prefix']+'.pkk-site.tsv'),len(samples),other['k'])
            if other['k']==5:
                v=same_k(a,b,p,pb)
                _,cols=__import__('scipy').optimize.linear_sum_assignment(np.abs(a[:,:,None]-b[:,None,:]).mean(0))
                v['matched_label_agreement']=float((labels==b[:,cols].argmax(1)).mean())
                v['pass_gate']=v['ari']>=.8 and v['matched_label_agreement']>=.8 and v['mean_l1']<=.2
                v['status']='PASS' if v['pass_gate'] else 'FAIL'; stability['seed']=v
                vv=dict(agreement=v['virtual_argmax_agreement'],pass_gate=v['virtual_argmax_agreement']>=.8)
                vv['status']='PASS' if vv['pass_gate'] else 'FAIL'; stability['virtual_argmax']=vv
            else:
                g,h=a@a.T,b@b.T; g/=np.linalg.norm(g); h/=np.linalg.norm(h)
                v=dict(ari=float(adjusted_rand_score(labels,b.argmax(1))),gram_difference=float(np.linalg.norm(g-h)))
                v['pass_gate']=v['ari']>=.6 and v['gram_difference']<=.3; v['status']='PASS' if v['pass_gate'] else 'FAIL'
                stability['cross_k']=v
        if any(stability[x]['status']=='FAIL' for x in ('seed','cross_k','virtual_argmax')):
            stability['status']='FAIL'
        virtual.update(r['joint_errors'])
        # No frozen replays run unless the completed primary leaves sufficient budget;
        # missing integration always remains explicitly NOT TESTED.
        virtual['arithmetic_status']='PASS' if virtual['ordered_sum_error']<=1e-10 and virtual['q_symmetry_error']<=1e-12 and virtual['q_triangle_sum_error']<=1e-10 else 'FAIL'
    coverage=dict(fitted_n=panel['fitted_n'],dated_selected=panel['dated_n'],usable_individuals=panel['usable_individuals'],
                  usable_fraction=panel['usable_individuals']/panel['fitted_n'],full_export_rows=r['full_export_rows'] if r else 0,
                  expected_rows=panel['fitted_n']*panel['retained_markers']*25,dated_complete_export=panel['dated_n'] if valid else 0,
                  status='PASS' if valid and panel['usable_individuals']/panel['fitted_n']>=.95 else 'FAIL',pass_gate=valid and panel['usable_individuals']/panel['fitted_n']>=.95)
    exclusion=dict(selected_core_markers=panel['excluded_markers'],retained_core_markers=r.get('excluded_retained_markers') if r else None,
                   site_inside_core=excluded(SITE),pc_overlap=0,left_flank=panel['left_flank'],right_flank=panel['right_flank'],
                   builder_status='PASS',status='PASS' if valid and r.get('excluded_retained_markers')==0 else 'NOT TESTED',pass_gate=valid and r.get('excluded_retained_markers')==0)
    em=r['convergence'] if r else dict(status='NOT TESTED',pass_gate=False,completed_steps=0,trace=[])
    gates=dict(known_structure=known,stability=stability,coverage=coverage,virtual_site=virtual,exclusion=exclusion,em_convergence=em)
    result=dict(verdict='PASS' if all(g['pass_gate'] for g in gates.values()) else 'STOP',gates=gates,
                fits=len(runs),disk_free_bytes=shutil.disk_usage('/').free)
    save('gate_metrics.json',result)
    return result


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    # Disable all git probes before using the required provenance harness.
    sys.path.insert(0,str(ROOT))
    import scripts.provenance as provenance
    provenance.get_git_info=lambda: {'disabled':'operator prohibits all git commands'}
    recorder=provenance.RunRecorder('elai-amy-dense')
    recorder.set_seed(SEED)
    for p in (PREREG,Path(__file__),Path(BINARY)):
        recorder.record_input(p)
    recorder.record_parameter('raw_input_stats',{str(p):dict(size=p.stat().st_size,mtime=p.stat().st_mtime) for p in paths()})
    work=Path(tempfile.mkdtemp(prefix='elai-amy-dense-',dir='/tmp')); (work/'output').mkdir()
    started=time.time(); runs=[]
    try:
        panel,samples,markers,pc=build(work)
        if panel['retained_markers']>=50000 and panel['pc_markers']>=5000 and panel['fitted_n']>=1000:
            seen_steps.clear()
            runs.append(fit(work,panel,markers,5,SEED,1200,True))
            for k,seed in ((5,SEED+1),(3,SEED)):
                if not runs[0]['complete'] or not runs[0]['convergence']['pass_gate']:
                    break
                projected=runs[0]['elapsed_seconds']*(k/5)**2
                left=1700-(time.time()-started)
                if projected>left:
                    break
                seen_steps.clear(); runs.append(fit(work,panel,markers,k,seed,min(1200,left)))
        metrics(panel,samples,pc,runs)
    except Exception as exc:
        save('execution_failure.json',dict(error=repr(exc),elapsed_seconds=time.time()-started))
        raise
    finally:
        shutil.rmtree(work)
        save('scratch_cleanup.json',dict(path=str(work),removed=not work.exists(),elapsed_seconds=time.time()-started))
        for p in OUT.iterdir():
            if p.is_file():
                recorder.record_output(p)
        recorder.save()
        subprocess.run(['df','-h','/'],check=True)


if __name__=='__main__':
    main()
