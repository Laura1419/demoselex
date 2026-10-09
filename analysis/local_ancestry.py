#!/usr/bin/env python3
"""Stage-1 unphased ancestry-pair HMM; contract: docs/local-ancestry-prereg.md."""
from __future__ import annotations
import json
import math
import re
import sys
import zipfile
from pathlib import Path
from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
from scipy.stats import t as student_t, norm
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.load_aadr import load_bundle, BUNDLE, MAIN_PREFIX

OUT = Path('data/derived/local_ancestry')
SOURCES = dict(WHG=['France_Mesolithic', 'Luxembourg_Loschbour_Mesolithic'],
               EHG=['Russia_Karelia_Mesolithic'], FARM=['Turkey_N'],
               IRAN=['Iran_GanjDareh_N'], STEPPE=['Russia_Samara_EBA_Yamnaya'])
NAMES = list(SOURCES)
PAIRS = np.array(list(combinations_with_replacement(range(5), 2)))
COMPOSITION = np.array([np.bincount(pair, minlength=5)/2 for pair in PAIRS])
SEED = 20261009


def save_json(name, obj):
    def finite(value):
        if isinstance(value, dict):
            return {k: finite(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [finite(v) for v in value]
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            return None
        return value.item() if isinstance(value, np.generic) else value
    (OUT/name).write_text(json.dumps(finite(obj), indent=2, allow_nan=False))


def is_diploid(ids):
    return np.array([x.endswith(('.DG', '.SG')) for x in ids])


def source_frequency(g, diploid):
    called = g <= 2
    weight = np.where(diploid, 2., 1.)[:, None]
    alt = np.where(called, g/2*weight, 0).sum(axis=0)
    alleles = (called*weight).sum(axis=0)
    return (alt+.5)/(alleles+1), called.mean(axis=0), called.sum(axis=0)


def transition(penalty):
    x, y = PAIRS[:, None, :], PAIRS[None, :, :]
    changes = np.minimum((x != y).sum(axis=2), (x != y[:, :, ::-1]).sum(axis=2))
    tr = np.exp(-penalty*changes)
    return tr/tr.sum(axis=1, keepdims=True)


def smooth_hmm(log_emission, penalty, viterbi=False):
    """Scaled forward/backward; independent initial prior, no cross-chromosome call."""
    n, w, k = log_emission.shape
    e = np.exp(log_emission-log_emission.max(axis=2, keepdims=True))
    tr = transition(penalty)
    f = np.empty_like(e)
    a = np.full((n, k), 1/k)
    for j in range(w):
        a = (a if j == 0 else a@tr)*e[:, j]
        a /= a.sum(axis=1, keepdims=True)
        f[:, j] = a
    b = np.ones((n, k))
    post = np.empty_like(e)
    for j in range(w-1, -1, -1):
        p = f[:, j]*b
        post[:, j] = p/p.sum(axis=1, keepdims=True)
        b = (b*e[:, j])@tr.T
        b /= b.sum(axis=1, keepdims=True)
    path = None
    if viterbi:
        back = np.zeros((n, w, k), dtype=np.uint8)
        score = log_emission[:, 0]-np.log(k)
        logtr = np.log(tr)
        for j in range(1, w):
            z = score[:, :, None]+logtr
            back[:, j] = z.argmax(axis=1)
            score = z.max(axis=1)+log_emission[:, j]
        path = np.empty((n, w), dtype=np.uint8)
        path[:, -1] = score.argmax(axis=1)
        for j in range(w-1, 0, -1):
            path[:, j-1] = back[np.arange(n), j, path[:, j]]
    return post, path


def paint(g, diploid, panel, width, penalty, viterbi=False):
    """Window likelihoods; unsupported emissions are flat, not observed ancestry."""
    p = panel[[f'p_{s}' for s in NAMES]].to_numpy()
    chrom = panel.chrom.to_numpy()
    pos = panel.pos.to_numpy()
    fractions, supports, infos, paths = [], [], [], []
    for c in sorted(set(chrom)):
        ix = np.flatnonzero(chrom == c)
        bins = pos[ix]//width
        unique, start = np.unique(bins, return_index=True)
        # Empty genomic windows are represented, so transitions never jump over gaps.
        allbins = np.arange(unique.min(), unique.max()+1)
        slots = unique-allbins[0]
        gg = g[:, ix]
        called = gg <= 2
        count = np.zeros((len(g), len(allbins)), dtype=np.int32)
        count[:, slots] = np.add.reduceat(called.astype(int), start, axis=1)
        logs = np.zeros((len(g), len(allbins), len(PAIRS)))
        for k, (a, b) in enumerate(PAIRS):
            pa, pb = p[ix, a], p[ix, b]
            probs = np.stack([(1-pa)*(1-pb), pa*(1-pb)+pb*(1-pa), pa*pb])
            vals = np.log(np.clip(probs[np.minimum(gg, 2), np.arange(len(ix))], 1e-8, 1))
            q = (pa+pb)/2
            pseudo = gg/2*np.log(q)+(1-gg/2)*np.log1p(-q)
            vals = np.where(diploid[:, None], vals, pseudo)
            vals = np.where(called, vals, 0)
            logs[:, slots, k] = np.add.reduceat(vals, start, axis=1)*.25
        supported = count >= 5
        logs[~supported] = 0
        post, path = smooth_hmm(logs, penalty, viterbi)
        fractions.append(post@COMPOSITION)
        supports.append(supported)
        infos.append(np.column_stack([np.full(len(allbins), c), allbins*width, (allbins+1)*width]))
        if path is not None:
            paths.append(path)
    return (np.concatenate(fractions, axis=1), np.concatenate(supports, axis=1),
            np.concatenate(infos), np.concatenate(paths, axis=1) if paths else None)


def summarize(fractions, supported, windows):
    outside = ~((windows[:, 0] == 1) & (windows[:, 1] < 104500000) & (windows[:, 2] > 103500000))
    valid = supported & outside[None, :]
    denominator = valid.sum(axis=1)
    global_q = (fractions*valid[:, :, None]).sum(axis=1)/np.maximum(denominator[:, None], 1)
    local = ~outside
    local_valid = supported & local[None, :]
    local_q = (fractions*local_valid[:, :, None]).sum(axis=1)/np.maximum(local_valid.sum(axis=1)[:, None], 1)
    return global_q, local_q, denominator, local_valid.sum(axis=1)


def fit_frequency(y, time, covariates):
    X = np.column_stack([np.ones(len(y)), covariates, time])
    rank = np.linalg.matrix_rank(X)
    if rank != X.shape[1] or len(y) < 30:
        return dict(identifiable=False, n=len(y), rank=int(rank), columns=X.shape[1])
    inv = np.linalg.inv(X.T@X)
    beta = inv@X.T@y
    residual = y-X@beta
    h = np.sum((X@inv)*X, axis=1)
    meat = X.T@((residual/(1-h))[:, None]**2*X)
    se = np.sqrt(np.diag(inv@meat@inv))
    df = len(y)-rank
    p = 2*student_t.sf(abs(beta[-1]/se[-1]), df)
    tr = time-X[:, :-1]@np.linalg.lstsq(X[:, :-1], time, rcond=None)[0]
    info = tr@tr
    floor = (student_t.ppf(.975, df)+norm.ppf(.8))*se[-1]
    ideal = 2.802*np.sqrt(y.mean()*(1-y.mean())/(2*info))
    return dict(identifiable=True, n=len(y), rank=int(rank), frequency=float(y.mean()),
                beta_time=float(beta[-1]), HC3_se=float(se[-1]), p=float(p),
                MDE80_HC3=float(floor), MDE80_ideal=float(ideal), MDE80_2kyr=float(2*floor),
                residual_time_information=float(info))


def extract_targets(inds, anno):
    scratch = OUT/'scratch'
    scratch.mkdir(parents=True, exist_ok=True)
    tables = {}
    with zipfile.ZipFile('data/raw/amy/amylase_diversity_project-v1.1.zip') as z:
        for name in ['haplotype_structure_for_supplement.tsv', 'ancient_metadata_for_supplement.tsv']:
            match = [n for n in z.namelist() if n.endswith('/'+name)]
            assert len(match) == 1
            path = scratch/name
            path.write_bytes(z.read(match[0]))
            tables[name] = pd.read_csv(path, sep='\t')
    h = tables['haplotype_structure_for_supplement.tsv']
    meta = tables['ancient_metadata_for_supplement.tsv']
    targets = h[h.type == 'Ancient'].merge(meta[['sample', 'age', 'p2_l0']], on='sample', validate='one_to_one')
    rows_by_bio, rows_by_prefix, exact = {}, {}, {}
    for i, ind in enumerate(inds):
        exact[ind.id] = i
        bio = anno[ind.id].individual_id if ind.id in anno else ind.id
        rows_by_bio.setdefault(bio, []).append(i)
        rows_by_prefix.setdefault(ind.id.split('.')[0], []).append(i)
    rows = []
    for _, r in targets.iterrows():
        sample = str(r['sample'])
        if sample in exact:
            candidates, method = [exact[sample]], 'exact_genetic'
        elif sample in rows_by_bio:
            candidates, method = rows_by_bio[sample], 'biological_ID'
        else:
            candidates, method = rows_by_prefix.get(sample, []), 'unique_prefix'
        # Multiple representations of ONE biological person are not ambiguous people.
        bios = {anno[inds[i].id].individual_id for i in candidates if inds[i].id in anno}
        if len(candidates) > 1 and len(bios) != 1:
            candidates = []
        if candidates:
            candidates.sort(key=lambda i: (-HITS.get(inds[i].id, 0), inds[i].id))
            i = candidates[0]
            a = anno.get(inds[i].id)
            rows.append(dict(r, aadr_index=i, genetic_id=inds[i].id, match_method=method,
                             aadr_date_bp=a.date_bp if a else None,
                             biological_id=a.individual_id if a else sample))
    matched = pd.DataFrame(rows).sort_values('sample').drop_duplicates('biological_id')
    matched.to_csv(OUT/'target_mapping.tsv', sep='\t', index=False)
    return matched, len(targets)


def record_provenance():
    from scripts import provenance as p
    p.get_git_info = lambda: {'omitted': 'operator prohibits all git commands'}
    recorder = p.RunRecorder('local-ancestry-stage1')
    recorder.set_seed(SEED)
    for path in [Path(__file__), Path('docs/local-ancestry-prereg.md')]:
        recorder.record_input(path)
    for path in Path(BUNDLE).glob('*'):
        if path.suffix in ['.geno', '.ind', '.snp', '.anno']:
            stat = path.stat()
            recorder.manifest['inputs'][str(path)] = dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    recorder.record_input('data/raw/amy/amylase_diversity_project-v1.1.zip')
    for path in OUT.glob('*'):
        if path.is_file() and path.suffix != '.log':
            recorder.record_output(path)
    recorder.save()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert Path('docs/local-ancestry-prereg.md').exists()
    gm, inds, anno = load_bundle()
    raw = pd.read_csv(Path(BUNDLE)/'v66.p1_1240K.aadr.PUB.anno', sep='\t', dtype=str).fillna('..')
    raw_rows = {r.iloc[0]: r for _, r in raw.iterrows()}
    global HITS
    HITS = {r.iloc[0]: pd.to_numeric(r.iloc[26], errors='coerce') for _, r in raw.iterrows()}
    HITS = {k: (v if np.isfinite(v) else 0) for k, v in HITS.items()}
    rng = np.random.default_rng(SEED)
    refs, counts = [], []
    for source, labels in SOURCES.items():
        eligible = []
        for i, ind in enumerate(inds):
            if ind.group not in labels or ind.id not in raw_rows:
                continue
            r = raw_rows[ind.id]
            if not r.iloc[47].startswith('Pass') or r.iloc[31] not in ['..', '', 'NA', 'n/a']:
                continue
            eligible.append(i)
        eligible.sort(key=lambda i: (-HITS[inds[i].id], inds[i].id))
        seen, unique = set(), []
        for i in eligible:
            bio = anno[inds[i].id].individual_id
            if bio not in seen:
                unique.append(i)
                seen.add(bio)
        unique.sort(key=lambda i: inds[i].id)
        shuffled = rng.permutation(unique)
        nh = min(max(2, math.ceil(len(unique)/4)), len(unique)-3)
        heldout = set(shuffled[:nh])
        for i in unique:
            refs.append(dict(aadr_index=i, genetic_id=inds[i].id, source=source, group=inds[i].group,
                             date_bp=anno[inds[i].id].date_bp, split='holdout' if i in heldout else 'train'))
        for label in labels:
            before = [i for i, x in enumerate(inds) if x.group == label]
            after = [i for i in unique if inds[i].group == label]
            dates = [anno[inds[i].id].date_bp for i in after if anno[inds[i].id].date_bp is not None]
            before_dates = [anno[inds[i].id].date_bp for i in before if inds[i].id in anno and anno[inds[i].id].date_bp is not None]
            counts.append(dict(source=source, label=label, published_n=len(before), eligible_n=len(after),
                               published_min_BP=min(before_dates, default=None), published_max_BP=max(before_dates, default=None),
                               train_n=sum(i not in heldout for i in after), holdout_n=sum(i in heldout for i in after),
                               min_BP=min(dates, default=None), max_BP=max(dates, default=None)))
    refs = pd.DataFrame(refs)
    refs.to_csv(OUT/'references.tsv', sep='\t', index=False)
    pd.DataFrame(counts).to_csv(OUT/'reference_counts.tsv', sep='\t', index=False)
    print('Reference panel:', counts, flush=True)
    targets, target_total = extract_targets(inds, anno)
    # Remove all training individuals by biological ID, even if target uses another representation.
    train_bio = {anno[r.genetic_id].individual_id for r in refs[refs.split == 'train'].itertuples()}
    targets = targets[~targets.biological_id.isin(train_bio)].reset_index(drop=True)
    snp = pd.read_csv(Path(BUNDLE)/(MAIN_PREFIX+'.snp'), sep=r'\s+', header=None,
                      names=['rsid', 'chrom', 'genetic_pos', 'pos', 'ref', 'alt'])
    candidate = snp.chrom.between(1, 22) & snp.ref.isin(list('ACGT')) & snp.alt.isin(list('ACGT')) & snp.ref.ne(snp.alt)
    freq, calls, observed = [], [], []
    for source in NAMES:
        rr = refs[(refs.source == source) & (refs.split == 'train')]
        gg = np.stack([gm.individual_genotypes(i) for i in rr.aadr_index])
        p, call, count = source_frequency(gg, is_diploid(rr.genetic_id))
        freq.append(p); calls.append(call); observed.append(count)
    freq, calls, observed = map(np.array, (freq, calls, observed))
    divergence = np.ptp(freq, axis=0)
    candidate &= (calls.min(axis=0) >= .5) & (observed.min(axis=0) >= 3) & (divergence >= .45)
    ci = np.flatnonzero(candidate)
    print('Candidate AIMs', len(ci), 'mapped targets', len(targets), flush=True)
    tg = np.stack([gm.individual_genotypes(i, ci) for i in targets.aadr_index])
    target_call = (tg <= 2).mean(axis=0)
    usable = ci[target_call >= .2]
    occupied, selected = {}, []
    import bisect
    for i in sorted(usable, key=lambda i: (-divergence[i], i)):
        c, pos = int(snp.chrom.iloc[i]), int(snp.pos.iloc[i])
        coords = occupied.setdefault(c, [])
        j = bisect.bisect_left(coords, pos)
        if (j and pos-coords[j-1] < 50000) or (j < len(coords) and coords[j]-pos < 50000):
            continue
        coords.insert(j, pos); selected.append(i)
    selected = np.array(sorted(selected, key=lambda i: (snp.chrom.iloc[i], snp.pos.iloc[i])))
    panel = snp.iloc[selected].copy()
    panel.insert(0, 'aadr_snp_index', selected)
    panel['divergence'] = divergence[selected]
    mean = freq[:, selected].mean(axis=0)
    panel['Fst_equal_source'] = freq[:, selected].var(axis=0)/(mean*(1-mean))
    panel['target_call_fraction'] = target_call[np.searchsorted(ci, selected)]
    for k, name in enumerate(NAMES):
        panel[f'p_{name}'] = freq[k, selected]
        panel[f'call_{name}'] = calls[k, selected]
    panel.to_csv(OUT/'aim_panel.tsv', sep='\t', index=False)
    allrows = sorted(set(refs.aadr_index) | set(targets.aadr_index))
    sample_lookup = {i: j for j, i in enumerate(allrows)}
    g = np.stack([gm.individual_genotypes(i, selected) for i in allrows])
    ids = [inds[i].id for i in allrows]
    diploid = is_diploid(ids)
    callrate = (g <= 2).mean(axis=1)
    held = refs[refs.split == 'holdout'].copy().reset_index(drop=True)
    hi = np.array([sample_lookup[i] for i in held.aadr_index])
    truth = np.array([NAMES.index(s) for s in held.source])
    ti = np.array([sample_lookup[i] for i in targets.aadr_index])
    print('Selected AIMs', len(panel), 'individuals', len(g), flush=True)
    q, local, nw, nl, settings = None, None, None, None, []
    # Primary first, then preregistered alternatives.
    grid = [(2000000, 4)]+[(w, p) for w in [1000000, 2000000, 4000000] for p in [2, 4, 8] if (w, p) != (2000000, 4)]
    primary_ok = None
    for width, penalty in grid:
        f, supported, windows, path = paint(g, diploid, panel, width, penalty, viterbi=width == 2000000 and penalty == 4)
        qq, ll, nn, ln = summarize(f, supported, windows)
        if q is None:
            q, local, nw, nl = qq, ll, nn, ln
            primary_windows = windows.copy()
            primary_ok = (nw >= 100) & (callrate >= .2)
            np.savez_compressed(OUT/'primary_paint.npz', fractions=f.astype('float32'), supported=supported,
                                windows=windows, viterbi_pair=path, genetic_ids=np.array(ids), pairs=PAIRS)
        valid_h = primary_ok[hi] & (nn[hi] >= 100) & (callrate[hi] >= .2)
        targetvalid = primary_ok[ti] & (nn[ti] >= 100)
        localvalid = (nl[ti] > 0) & (ln[ti] > 0)
        stat = dict(window_bp=width, switch_penalty=penalty,
                    heldout_median_L1=float(np.median(np.abs(qq[hi][valid_h]-q[hi][valid_h]).sum(axis=1))),
                    heldout_dominant_agreement=float(np.mean(qq[hi][valid_h].argmax(axis=1) == q[hi][valid_h].argmax(axis=1))),
                    target_median_L1=float(np.median(np.abs(qq[ti][targetvalid]-q[ti][targetvalid]).sum(axis=1))),
                    target_local_median_L1=float(np.median(np.abs(ll[ti][localvalid]-local[ti][localvalid]).sum(axis=1))),
                    target_local_dominant_agreement=float(np.mean(ll[ti][localvalid].argmax(axis=1) == local[ti][localvalid].argmax(axis=1))))
        settings.append(stat)
        np.savez_compressed(OUT/f'composition_{width}_{penalty}.npz', global_q=qq.astype('float32'), local_q=ll.astype('float32'),
                            supported_windows=nn, supported_local_windows=ln)
        pd.DataFrame(settings).to_csv(OUT/'sensitivity.tsv', sep='\t', index=False)
        print('Sensitivity', stat, flush=True)
    valid_h = primary_ok[hi]
    pred = q[hi].argmax(axis=1)
    confusion = np.zeros((5, 5), dtype=int)
    for a, b in zip(truth[valid_h], pred[valid_h]):
        confusion[a, b] += 1
    pd.DataFrame(confusion, index=NAMES, columns=NAMES).to_csv(OUT/'confusion.tsv', sep='\t', index_label='published_source')
    for k, source in enumerate(NAMES):
        held[f'q_{source}'] = q[hi, k]
    held['predicted'] = np.array(NAMES)[pred]
    held['evaluable'] = valid_h
    held['call_fraction'] = callrate[hi]
    held.to_csv(OUT/'holdout_validation.tsv', sep='\t', index=False)
    per_source = []
    for k, source in enumerate(NAMES):
        mask = valid_h & (truth == k)
        per_source.append(dict(source=source, evaluable=int(mask.sum()), accuracy=float(np.mean(pred[mask] == k)) if mask.any() else None,
                               median_correct_fraction=float(np.median(q[hi[mask], k])) if mask.any() else None))
    control_candidates = np.flatnonzero((truth == 0) & valid_h)
    nc = control_candidates[np.argmax(callrate[hi[control_candidates]])] if len(control_candidates) else None
    negative = dict(genetic_id=held.genetic_id.iloc[nc], source='WHG', WHG_fraction=float(q[hi[nc], 0])) if nc is not None else dict(WHG_fraction=0)
    # PCs on a disjoint, evenly-spaced marker set, not the AIM panel.
    blocked = set(selected)
    pc_sites, last = [], {}
    for i, r in enumerate(snp.itertuples(index=False)):
        if i in blocked or not 1 <= r.chrom <= 22 or (r.chrom == 1 and 103500000 <= r.pos <= 104500000):
            continue
        if r.pos-last.get(r.chrom, -200000) >= 200000:
            pc_sites.append(i); last[r.chrom] = r.pos
    if len(pc_sites) > 10000:
        pc_sites = np.array(pc_sites)[np.linspace(0, len(pc_sites)-1, 10000, dtype=int)]
    pcg = np.stack([gm.individual_genotypes(i, np.array(pc_sites)) for i in allrows]).astype(float)
    pcg[pcg > 2] = np.nan
    counts_pc = np.isfinite(pcg).sum(axis=0)
    means = np.divide(np.nansum(pcg, axis=0), counts_pc, out=np.full(pcg.shape[1], np.nan), where=counts_pc > 0)
    good = np.isfinite(means)
    pcg = pcg[:, good]; means = means[good]
    pcg = np.where(np.isnan(pcg), means, pcg)
    sd = pcg.std(axis=0)
    pcg = (pcg[:, sd > 0]-means[sd > 0])/sd[sd > 0]
    pcs = PCA(n_components=3, svd_solver='randomized', random_state=SEED).fit_transform(pcg)
    x = np.column_stack([np.ones(primary_ok.sum()), q[primary_ok, :4]])
    pp = pcs[primary_ok]
    fitted = x@np.linalg.lstsq(x, pp, rcond=None)[0]
    r2 = 1-((pp-fitted)**2).sum(axis=0)/((pp-pp.mean(axis=0))**2).sum(axis=0)
    pc_table = pd.DataFrame(dict(genetic_id=ids, PC1=pcs[:, 0], PC2=pcs[:, 1], PC3=pcs[:, 2]))
    for k, source in enumerate(NAMES): pc_table[f'q_{source}'] = q[:, k]
    pc_table.to_csv(OUT/'pc_composition.tsv', sep='\t', index=False)
    pc_table_sites = snp.iloc[pc_sites]
    pc_table_sites.to_csv(OUT/'pc_marker_panel.tsv', sep='\t', index=False)
    gates = dict(heldout_count=all(r['evaluable'] >= 2 for r in per_source),
                 overall_accuracy=bool(np.mean(pred[valid_h] == truth[valid_h]) >= .8),
                 per_source_accuracy=all(r['accuracy'] is not None and r['accuracy'] >= .6 for r in per_source),
                 source_dominance=all(r['median_correct_fraction'] is not None and r['median_correct_fraction'] >= .6 for r in per_source),
                 negative_control=negative['WHG_fraction'] >= .8,
                 PC_agreement=bool(r2.mean() >= .5),
                 stability=all(s['heldout_median_L1'] <= .25 and s['heldout_dominant_agreement'] >= .8 for s in settings))
    # Count local observed AIMs in supported overlapping windows (not interpolated gaps).
    amy_observed = np.zeros(len(g), dtype=int)
    for c, lo, end in primary_windows[(primary_windows[:, 0] == 1) & (primary_windows[:, 1] < 104500000) & (primary_windows[:, 2] > 103500000)]:
        ix = (panel.chrom.to_numpy() == c) & (panel.pos.to_numpy() >= lo) & (panel.pos.to_numpy() < end)
        count = (g[:, ix] <= 2).sum(axis=1)
        amy_observed += np.where(count >= 5, count, 0)
    local_ok = primary_ok & (amy_observed >= 10) & (nl > 0)
    target_output = targets.copy()
    target_output['global_evaluable'] = primary_ok[ti]
    target_output['AMY_evaluable'] = local_ok[ti]
    target_output['AMY_observed_AIMs'] = amy_observed[ti]
    target_output['global_supported_windows'] = nw[ti]
    target_output['AIM_call_fraction'] = callrate[ti]
    for k, source in enumerate(NAMES):
        target_output[f'global_{source}'] = q[ti, k]
        target_output[f'local_{source}'] = local[ti, k]
    target_output.to_csv(OUT/'target_composition.tsv', sep='\t', index=False)
    validation = dict(gates=gates, passed=all(gates.values()), source_validation=per_source,
                      overall_accuracy=float(np.mean(pred[valid_h] == truth[valid_h])), negative_control=negative,
                      PC_R2=r2.tolist(), PC_mean_R2=float(r2.mean()), pc_markers=pcg.shape[1],
                      targets_in_archive=target_total, targets_mapped=len(pd.read_csv(OUT/'target_mapping.tsv', sep='\t')),
                      targets_excluding_training=len(targets), targets_global_painted=int(primary_ok[ti].sum()),
                      targets_AMY_supported=int(local_ok[ti].sum()), AIMs=len(panel),
                      candidate_AIMs=len(ci), heldout_evaluable=int(valid_h.sum()))
    save_json('validation.json', validation)
    print('VALIDATION', json.dumps(validation), flush=True)
    if not validation['passed']:
        save_json('amy_status.json', dict(status='STOPPED: preregistered validation gate failed', failed_gates=[k for k, v in gates.items() if not v]))
        record_provenance()
        return
    # Only after passing all validation gates may structural outcome be constructed.
    frame = target_output[local_ok[ti] & np.isfinite(pd.to_numeric(targets.age, errors='coerce'))].copy()
    results = []
    for cut in [1, 2]:
        def long(h):
            match = re.match(r'^H(\d+)', str(h))
            if not match: raise ValueError(f'Unrecognized haplotype {h}')
            return int(match.group(1)) > cut
        labelled = frame.H1.notna() & frame.H2.notna()
        ff = frame[labelled].copy()
        yy = (ff.H1.map(long).astype(int)+ff.H2.map(long).astype(int)).to_numpy()/2
        time = -ff.age.to_numpy(float)/1000
        localcov = ff[[f'local_{s}' for s in NAMES[:4]]].to_numpy()
        globalcov = pd.get_dummies(ff.p2_l0, drop_first=True).to_numpy(float)
        for label, cov in [('local', localcov), ('global_label', globalcov), ('pooled', np.empty((len(yy), 0)))]:
            results.append(dict(cut=cut, model=label, **fit_frequency(yy, time, cov)))
        ff['B_dosage'] = yy*2
        ff['dominant_local'] = np.array(NAMES)[ff[[f'local_{s}' for s in NAMES]].to_numpy().argmax(axis=1)]
        ff['BP_bin'] = pd.cut(ff.age, [0, 2000, 4000, 6000, 8000, 10000, np.inf], right=False)
        agg = ff.groupby(['dominant_local', 'BP_bin'], observed=True).agg(n=('sample', 'size'), B_copies=('B_dosage', 'sum'), mean_BP=('age', 'mean')).reset_index()
        agg['alleles'] = agg.n*2; agg['B_frequency'] = agg.B_copies/agg.alleles
        agg.to_csv(OUT/f'amy_frequency_cut{cut}.tsv', sep='\t', index=False)
        ff.to_csv(OUT/f'amy_cohort_cut{cut}.tsv', sep='\t', index=False)
        for s, group in ff.groupby('dominant_local'):
            if len(group) >= 10 and np.ptp(group.age)/1000 >= 1:
                results.append(dict(cut=cut, model='within_dominant_'+s, **fit_frequency(group.B_dosage.to_numpy()/2, -group.age.to_numpy()/1000, np.empty((len(group), 0)))))
        # Regression sensitivity on identical primary cohort, no reselected support.
        idx = np.array([sample_lookup[i] for i in ff.aadr_index])
        for width, penalty in grid[1:]:
            cached = np.load(OUT/f'composition_{width}_{penalty}.npz')
            if np.all(cached['supported_local_windows'][idx] > 0):
                results.append(dict(cut=cut, model=f'local_{width}_{penalty}', **fit_frequency(yy, time, cached['local_q'][idx, :4])))
    pd.DataFrame(results).to_csv(OUT/'amy_statistics.tsv', sep='\t', index=False)
    save_json('amy_status.json', dict(status='completed; diploid composition conditioning, not phased haplotype test'))
    record_provenance()


if __name__ == '__main__':
    main()
