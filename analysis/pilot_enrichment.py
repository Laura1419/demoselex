#!/usr/bin/env python3
"""Within-state enrichment pilot (prereg + AMENDMENT 1, 2026-10-13).

Per (variant, state): enrichment trajectory LOR_s(b) = logit(p_s(b)) - logit(p_other(b))
over four shared date-quartile bins; step statistic STEP_s = max adjacent |dLOR|;
enrichment-at-first-detection E0_s; background expansion curve f_s(b); within-window
permutation null (dates shuffled among confident individuals); drift baseline from the
null windows; state-specificity discriminator. Deterministic (crc32-window seeds).
Writes data/derived/pilot_enrichment/{states.csv,summary.json}.
"""
import json, sys, zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, iter_pvar, paths  # noqa: E402

GENOME = ROOT / 'data/derived/elai_genome'
OUT = ROOT / 'data/derived/pilot_enrichment'
S, K = 2, 10
SEED_A = 20261020
SEED_B = 20261021
MIN_W = 30.0   # minimum CALLED genotypes (2*weighted ncp) per bin
MIN_STATE_N = 60
MIN_BIN_N = 10          # called genotypes per bin
MAX_VARIANTS = 600
CONF = 0.9
PERMS = 100             # pilot; 1000 in the full run
RNG_SEED = 20261020
KNOWN = {'LCT': 'c2_w135', 'FADS': 'c11_w61', 'SLC45A2': 'c5_w33', 'SLC24A5': 'c15_w48'}
AMY = ['c1_w103', 'c1_w104', 'c1_w105']
EPS = 1e-3


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def pilot_tags():
    qc = pd.read_csv(GENOME / 'qc_table.csv')
    gated = qc[(qc.ok) & (qc.ari > .5) & (qc.frac90 > .3)]
    exclude = set(KNOWN.values()) | set(AMY) | \
        {w['tag'] for w in json.load(open(GENOME / 'build.manifest.json'))['windows']
         if w['chrom'] == 1 and 100 <= w['bin'] <= 108}
    pool = gated[~gated.tag.isin(exclude)]
    rng = np.random.default_rng(RNG_SEED)
    picks = rng.choice(len(pool), size=24, replace=False)
    return [t for _, t in sorted(KNOWN.items(), key=lambda kv: kv[1])] + AMY + sorted(pool.tag.iloc[picks])


def collect_indices(tags, win_by_tag):
    want = {}
    for t in tags:
        w = win_by_tag[t]
        want.setdefault(w['chrom'], []).append(w)
    found = {t: [] for t in tags}
    for i, (chrom, pos, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
        c = int(chrom) if chrom.isdigit() else -1
        if c not in want:
            continue
        for w in want[c]:
            if w['first'] <= pos <= w['last'] and len(ref) == 1 and len(alt) == 1 \
                    and ref in 'ACGT' and alt in 'ACGT' and ref != alt:
                found[w['tag']].append((i, pos, vid, ref, alt))
    for t, rows in found.items():
        rows.sort()
        if len(rows) > MAX_VARIANTS:
            step = len(rows) / MAX_VARIANTS
            found[t] = [rows[int(k * step)] for k in range(MAX_VARIANTS)]
    return found


def freqs(G, members, bins, binof):
    """G (V,n) with NaN missing. Returns p (V,bins), n_called (bins,)."""
    p = np.full((G.shape[0], bins), np.nan)
    nc = np.zeros(bins, dtype=int)
    for b in range(bins):
        idx = members & (binof == b)
        if not idx.any():
            continue
        vals = G[:, idx]
        called = np.isfinite(vals)
        den = called.sum(1) * 2
        ok = den >= MIN_BIN_N
        with np.errstate(invalid='ignore', divide='ignore'):
            p[:, b] = np.where(ok, np.nansum(np.where(called, vals, 0), axis=1) / np.where(den == 0, 1, den), np.nan)
        nc[b] = int(den.mean() / 2 * idx.sum()) if idx.sum() else 0
    return p, nc


def trajectories(G, state, bins, binof, states):
    """Per state: p_s, p_other, LOR (V,bins)."""
    conf = state >= 0
    out = {}
    for s in states:
        m_s = state == s
        if m_s.sum() < MIN_STATE_N:
            continue
        p_s, _ = freqs(G, m_s, bins, binof)
        p_o, _ = freqs(G, conf & ~m_s, bins, binof)
        L = logit(p_s) - logit(p_o)
        out[s] = (p_s, p_o, L)
    return out


def step_and_e0(L):
    d = np.abs(np.diff(L, axis=1))
    step = np.nanmax(d, axis=1)
    e0 = L[:, 0]
    step = np.where(np.isfinite(step), step, np.nan)
    return step, e0


def perm_null(G, state, bins, binof, states, perms, seed):
    conf = np.where(state >= 0)[0]
    rng = np.random.default_rng(seed)
    out = {s: np.empty((G.shape[0], perms)) for s in states}
    for k in range(perms):
        shuf = conf.copy()
        rng.shuffle(shuf)
        pseudo = state.copy()
        pseudo[conf] = state[shuf]
        tr = trajectories(G, pseudo, bins, binof, states)
        for s, (_p, _o, L) in tr.items():
            out[s][:, k] = step_and_e0(L)[0]
    return out


def process_window(t, cache):
    rows = cache['rows']
    fq = GENOME / 'windows' / t / f'w_s{S}_k{K}_seed{SEED_A}.q.npy'
    if not fq.exists() or not rows:
        return []
    qc, xmap, binof, bins, PERMS_ = cache['qc'], cache['xmap'], cache['binof'], cache['bins'], cache['PERMS_']
    win_by_tag = cache['win_by_tag']
    Q = np.load(fq)
    hard = Q.argmax(1)
    conf = Q.max(1) > CONF
    state = np.where(conf, hard, -1)
    states = [int(s) for s in np.unique(state[state >= 0]) if (state == s).sum() >= MIN_STATE_N]
    fcurves = {s: [float(((state == s) & (binof == b)).sum() / max((binof == b).sum(), 1))
                   for b in range(bins)] for s in states}
    vi = np.array([r[0] for r in rows], dtype=np.uint32)
    G = cache['read'](vi).astype(float)
    G[G > 2] = np.nan
    tr = trajectories(G, state, bins, binof, states)
    seed = (RNG_SEED + zlib.crc32(t.encode())) % 2**32
    null = perm_null(G, state, bins, binof, states, PERMS_, seed)
    ari = float(qc.loc[t, 'ari'])
    out = []
    for k, (i, pos, vid, ref, alt) in enumerate(rows):
        for s in states:
            p_s, p_o, L = tr[s]
            st, e0 = step_and_e0(L[k:k + 1])
            st, e0 = float(st[0]), float(e0[0])
            if not np.isfinite(st):
                continue
            nperm = int((null[s][k] >= st).sum())
            out.append(dict(
                window=t, chrom=win_by_tag[t]['chrom'], pos=pos, vid=vid, ref=ref, alt=alt,
                akbari_X=xmap.get((win_by_tag[t]['chrom'], pos), np.nan), ari=ari,
                state=s, step=st, e0=e0, perm_p=(nperm + 1) / (PERMS_ + 1),
                L_traj=json.dumps([round(x, 3) for x in L[k]]),
                p_traj=json.dumps([round(x, 4) for x in p_s[k]]),
                f_curve=json.dumps([round(x, 3) for x in fcurves[s]]),
                n_states_conf=len(states)))
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--combine', action='store_true', help='merge states_<tag>.csv files and compute baseline/summary')
    ap.add_argument('--v2', action='store_true', help='AMENDMENT-2 statistic (top-5+OTHER, soft)')
    ap.add_argument('--combine2', action='store_true', help='merge v2 states2_<tag>.csv + baseline/summary')
    ap.add_argument('--window', default=None)
    ap.add_argument('--perms', type=int, default=PERMS)
    ap.add_argument('--out', default=str(OUT))
    a = ap.parse_args()
    if a.combine2:
        combine2(Path(a.out)); return
    perms_n = a.perms
    outdir = Path(a.out); outdir.mkdir(parents=True, exist_ok=True)
    if a.combine:
        import glob as _g
        dfs = [pd.read_csv(f) for f in sorted(_g.glob(str(outdir / 'states_c*.csv')))]
        state_rows = pd.concat(dfs, ignore_index=True)
        win_by_tag = {w['tag']: w for w in json.load(open(GENOME / 'build.manifest.json'))['windows']}
        edges = None
        df = state_rows
        if df.empty:
            print('NO ROWS'); return
        df['arm'] = np.where(df.window.isin(KNOWN.values()), 'known',
                             np.where(df.window.isin(AMY), 'amy', 'null'))
        df['step_second'] = df.groupby(['window', 'pos']).step.transform(
            lambda x: x.nlargest(2).iloc[-1] if len(x) > 1 else np.nan)
        df['grp_max'] = df.groupby(['window', 'pos']).step.transform('max')
        best = df[df.step == df.grp_max].drop_duplicates(['window', 'pos']).copy()
        base = df[df.arm == 'null'].step.dropna()
        q95, q99, q999 = base.quantile([.95, .99, .999])
        med = float(base.median())
        best['above_q99'] = best.step > q99
        best['specific'] = best.above_q99 & (best.step_second.fillna(0) < med)
        df = df.merge(best[['window', 'pos', 'state', 'above_q99', 'specific', 'step_second']],
                      on=['window', 'pos', 'state'], how='left')
        df.to_csv(outdir / 'states.csv', index=False)
        print(f'drift baseline ({len(base)} null pairs): q95={q95:.3f} q99={q99:.3f} q999={q999:.3f} median={med:.3f}', flush=True)
        summary = dict(drift_baseline=dict(q95=float(q95), q99=float(q99), q999=float(q999),
                                           median=float(med), n=int(len(base))))
        for arm in ('known', 'amy'):
            sub = best[(best.arm == arm) & best.above_q99].sort_values('step', ascending=False)
            summary[arm] = dict(
                n_pairs=int(len(best[best.arm == arm])), n_above_q99=int(len(sub)),
                n_specific=int(sub.specific.sum()) if len(sub) else 0,
                top=[dict(window=r.window, pos=int(r.pos), state=int(r.state), step=round(r.step, 3),
                          e0=round(r.e0, 3), perm_p=round(r.perm_p, 4),
                          step_second=(None if not np.isfinite(r.step_second) else round(r.step_second, 3)),
                          X=(None if not np.isfinite(r.akbari_X) else round(r.akbari_X, 2)))
                     for r in sub.head(10).itertuples()])
        both = df[np.isfinite(df.akbari_X) & (df.arm == 'null')]
        summary['x_vs_step_spearman_null'] = float(both.akbari_X.corr(both.step, method='spearman')) if len(both) > 10 else None
        json.dump(summary, open(outdir / 'summary.json', 'w'), indent=1, default=str)
        print(json.dumps(summary, indent=1, default=str))
        return
    win_by_tag = {w['tag']: w for w in json.load(open(GENOME / 'build.manifest.json'))['windows']}
    tags_all = pilot_tags()
    if a.window:
        tags_all = [a.window]
    found = collect_indices(tags_all, win_by_tag)
    samples = json.load(open(ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'))
    years = 1950 - np.array([s['BP'] for s in samples], dtype=float)
    edges = np.quantile(years, [0, .25, .5, .75, 1.0])
    binof = np.digitize(years, edges[1:-1], right=False)
    bins = 4
    print('date bin edges (yr AD):', edges.astype(int).tolist(), flush=True)

    xmap = {}
    for fn in ('akbari_top.tsv', 'akbari_signal_absXge4.tsv', 'akbari_background200k.tsv'):
        p = ROOT / 'data/derived' / fn
        if p.exists():
            dfx = pd.read_csv(p, sep='\t', usecols=['CHROM', 'POS', 'X'])
            for c, pos, x in zip(dfx.CHROM, dfx.POS, dfx.X):
                xmap[(int(c), int(pos))] = float(x)

    qc = pd.read_csv(GENOME / 'qc_table.csv').set_index('tag')
    sample_idx = np.array([s['index'] for s in samples], dtype=np.uint32)
    state_rows = []
    with AkbariGeno(sample_subset=sample_idx) as reader:
        cache = dict(win_by_tag=win_by_tag, qc=qc, xmap=xmap, binof=binof,
                     bins=bins, PERMS_=perms_n, read=reader.read_hardcalls)
        if a.v2:
            rows_out = v2_process_window(a.window, {**cache, 'rows': found.get(a.window, [])})
            df = pd.DataFrame(rows_out)
            df.to_csv(outdir / f'states2_{a.window}.csv', index=False)
            print('wrote', len(df), 'v2 rows for', a.window); return
        if a.window:
            state_rows.extend(process_window(a.window, {**cache, 'rows': found.get(a.window, [])}))
        else:
          for t in tags_all:
            state_rows.extend(process_window(t, {**cache, 'rows': found.get(t, [])}))
            print('done', t, flush=True)

    df = pd.DataFrame(state_rows)
    if a.window:
        df.to_csv(outdir / f'states_{a.window}.csv', index=False)
        print('wrote', outdir / f'states_{a.window}.csv', len(df), 'rows'); return
    if df.empty:
        print('NO ROWS'); return
    df['arm'] = np.where(df.window.isin(KNOWN.values()), 'known',
                         np.where(df.window.isin(AMY), 'amy', 'null'))
    # state-specificity: per (window,pos), best state step vs second-best
    df['step_second'] = df.groupby(['window', 'pos']).step.transform(
        lambda x: x.nlargest(2).iloc[-1] if len(x) > 1 else np.nan)
    df['grp_max'] = df.groupby(['window', 'pos']).step.transform('max')
    best = df[df.step == df.grp_max].drop_duplicates(['window', 'pos']).copy()
    base = df[df.arm == 'null'].step.dropna()
    q95, q99, q999 = base.quantile([.95, .99, .999])
    med = float(base.median())
    best['above_q99'] = best.step > q99
    best['specific'] = best.above_q99 & (best.step_second.fillna(0) < med)
    df = df.merge(best[['window', 'pos', 'state', 'above_q99', 'specific', 'step_second']],
                  on=['window', 'pos', 'state'], how='left')
    df.to_csv(outdir / 'states.csv', index=False)
    print(f'drift baseline ({len(base)} null pairs): q95={q95:.3f} q99={q99:.3f} q999={q999:.3f} median={med:.3f}', flush=True)
    summary = dict(bins=edges.astype(int).tolist(), perms=PERMS,
                   drift_baseline=dict(q95=float(q95), q99=float(q99), q999=float(q999),
                                       median=float(med), n=int(len(base))))
    for arm in ('known', 'amy'):
        sub = best[(best.arm == arm) & best.above_q99].sort_values('step', ascending=False)
        summary[arm] = dict(
            n_pairs=int(len(best[best.arm == arm])), n_above_q99=int(len(sub)),
            n_specific=int(sub.specific.sum()) if len(sub) else 0,
            top=[dict(window=r.window, pos=int(r.pos), state=int(r.state), step=round(r.step, 3),
                      e0=round(r.e0, 3), perm_p=round(r.perm_p, 4),
                      step_second=(None if not np.isfinite(r.step_second) else round(r.step_second, 3)),
                      X=(None if not np.isfinite(r.akbari_X) else round(r.akbari_X, 2)))
                 for r in sub.head(10).itertuples()])
    both = df[np.isfinite(df.akbari_X) & (df.arm == 'null')]
    summary['x_vs_step_spearman_null'] = float(both.akbari_X.corr(both.step, method='spearman')) if len(both) > 10 else None
    json.dump(summary, open(OUT / 'summary.json', 'w'), indent=1, default=str)
    print(json.dumps(summary, indent=1, default=str))


def combine2(outdir):
    import glob as _g
    dfs = [pd.read_csv(f) for f in sorted(_g.glob(str(outdir / 'states2_c*.csv')))]
    df = pd.concat(dfs, ignore_index=True)
    df['arm'] = np.where(df.window.isin(list(KNOWN.values())), 'known',
                         np.where(df.window.isin(AMY), 'amy', 'null'))
    df['step_second'] = df.groupby(['window', 'pos']).step.transform(
        lambda x: x.nlargest(2).iloc[-1] if len(x) > 1 else np.nan)
    df['grp_max'] = df.groupby(['window', 'pos']).step.transform('max')
    best = df[df.step == df.grp_max].drop_duplicates(['window', 'pos']).copy()
    base = df[df.arm == 'null'].step.dropna()
    q95, q99, q999 = base.quantile([.95, .99, .999])
    med = float(base.median())
    best['above_q99'] = best.step > q99
    best['specific'] = best.above_q99 & (best.step_second.fillna(0) < med)
    best_cols = best[['window', 'pos', 'group', 'above_q99', 'specific', 'step_second']].copy()
    df = df.merge(best_cols, on=['window', 'pos', 'group'], how='left', suffixes=('', '_b'))
    if 'arm_b' in df.columns:
        df = df.drop(columns=['arm_b'])
    df['arm'] = df['arm'].fillna('null')
    import sys as _sys; print('PRE-WRITE arms:', df.arm.value_counts().to_dict(), file=_sys.stderr, flush=True)
    df.to_csv(outdir / 'states2.csv', index=False)
    print(f'drift baseline ({len(base)} null pairs): q95={q95:.3f} q99={q99:.3f} q999={q999:.3f} median={med:.3f}', flush=True)
    summary = dict(drift_baseline=dict(q95=float(q95), q99=float(q99), q999=float(q999),
                                       median=float(med), n=int(len(base))))
    for arm in ('known', 'amy'):
        sub = best[(best.arm == arm) & best.above_q99].sort_values('step', ascending=False)
        summary[arm] = dict(
            n_pairs=int(len(best[best.arm == arm])), n_above_q99=int(len(sub)),
            n_specific=int(sub.specific.sum()) if len(sub) else 0,
            top=[dict(window=r.window, pos=int(r.pos), group=r.group_name, step=round(r.step, 3),
                      e0=(None if not np.isfinite(r.e0) else round(r.e0, 3)),
                      perm_p=round(r.perm_p, 4),
                      step_second=(None if not np.isfinite(r.step_second) else round(r.step_second, 3)),
                      X=(None if not np.isfinite(r.akbari_X) else round(r.akbari_X, 2)))
                 for r in sub.head(12).itertuples()])
    json.dump(summary, open(outdir / 'summary2.json', 'w'), indent=1, default=str)
    print(json.dumps(summary, indent=1, default=str))




# ---------------------------------------------------------------------------
# v2 (AMENDMENT 2): per-window top-5 states + OTHER, soft q-weighted assignment.
# ---------------------------------------------------------------------------

def pair_columns(k=10):
    """Column s of the unordered q matrix corresponds to pair (k1,k2), k1<=k2,
    in the fork's q(a,b) order. Rebuild that list."""
    return [(a, b) for a in range(k) for b in range(a, k)]


def soft_weights(Q, top, k=10):
    """Q (n,55) unordered-PAIR posterior; top = list of 5 pair-column ids kept as
    groups; every other pair column maps to OTHER (index len(top)). Column s IS the
    state id in hard-assignment space (the unordered pair), so grouping is at pair
    level. Returns W (n, 6) with row sums 1."""
    n = Q.shape[0]
    Gt = len(top) + 1
    grp = {int(s): gi for gi, s in enumerate(top)}
    W = np.zeros((n, Gt))
    for s in range(Q.shape[1]):
        W[:, grp.get(s, Gt - 1)] += Q[:, s]
    return W


def weighted_freqs(G, W, binof, bins, min_w=30.0):
    """G (V,n) NaN=missing; W (n,Gt) weights. P (V,Gt,bins) freq; N (Gt,bins) weighted
    called-genotype counts. freq = sum w_i g_i / (2 sum w_i called_i)."""
    V, Gt = G.shape[0], W.shape[1]
    P = np.full((V, Gt, bins), np.nan)
    N = np.zeros((Gt, bins))
    called = np.isfinite(G).astype(float)
    gv = np.nan_to_num(G, nan=0.0)
    for b in range(bins):
        mb = (binof == b)
        w = W[mb]                                # (n_b, Gt)
        c = called[:, mb]                        # (V, n_b)
        cb = called[:, mb]                       # (V, n_b) per-variant call indicator
        ncp = np.einsum('ig,ij->gj', w, cb.T)    # (Gt, V) weighted called per variant
        cnt = gv[:, mb] @ w                      # (V, Gt) weighted copy counts
        with np.errstate(invalid='ignore', divide='ignore'):
            f = cnt / (2 * ncp.T)                # (V, Gt)
        ok = (2 * ncp) >= min_w                  # (Gt, V)
        P[:, :, b] = np.where(ok.T, f, np.nan)
        N[:, b] = 2 * ncp.mean(1)
    return P, N


def v2_process_window(t, cache):
    import numpy as np, json as _json
    rows, qc, xmap, binof, bins, PERMS_ = cache['rows'], cache['qc'], cache['xmap'], cache['binof'], cache['bins'], cache['PERMS_']
    win_by_tag = cache['win_by_tag']
    fq = GENOME / 'windows' / t / f'w_s{S}_k{K}_seed{SEED_A}.q.npy'
    fq_b = GENOME / 'windows' / t / f'w_s{S}_k{K}_seed{SEED_B}.q.npy'
    if not fq.exists() or not rows:
        return []
    Q = np.load(fq)
    hard = Q.argmax(1)
    counts = np.bincount(hard, minlength=55)
    top = [int(s) for s in np.argsort(counts)[::-1][:5]]
    W = soft_weights(Q, top, k=K)
    # seed-B agreement on 6-group assignment
    ari_b = np.nan
    if fq_b.exists():
        Qb = np.load(fq_b)
        hardb = Qb.argmax(1); countsb = np.bincount(hardb, minlength=55)
        topb = [int(s) for s in np.argsort(countsb)[::-1][:5]]
        Wb = soft_weights(Qb, topb, k=K)
        ga, gb = W.argmax(1), Wb.argmax(1)
        from sklearn.metrics import adjusted_rand_score
        ari_b = float(adjusted_rand_score(ga, gb))
    fcurves = [[float(W[binof == b, g].sum() / max(W[binof == b].sum(), 1e-9)) for b in range(bins)]
               for g in range(W.shape[1])]
    vi = np.array([r[0] for r in rows], dtype=np.uint32)
    Gg = cache['read'](vi).astype(float)
    Gg[Gg > 2] = np.nan
    P, N = weighted_freqs(Gg, W, binof, bins, min_w=MIN_W)
    L = logit(P) - logit(1 - P)  # placeholder replaced below
    # enrichment LOR: group g vs pooled others (weighted), per bin
    V, Gt = P.shape[0], P.shape[1]
    LOR = np.full((V, Gt, bins), np.nan)
    for g in range(Gt):
        others = [o for o in range(Gt) if o != g]
        # pooled other frequency: weighted mean of other groups' freq by their N
        num = np.nansum(np.stack([P[:, o, :] * N[o, :][None, :] for o in others], 0), 0)
        den = np.nansum(np.stack([np.where(np.isfinite(P[:, o, :]), N[o, :][None, :], 0) for o in others], 0), 0)
        p_other = num / np.maximum(den, 1e-9)
        LOR[:, g, :] = logit(P[:, g, :]) - logit(np.clip(p_other, EPS, 1 - EPS))
    # permutation null on the SAME machinery: shuffle years (soft weights keep individuals,
    # so permute the binof vector among individuals)
    rng = np.random.default_rng((RNG_SEED + zlib.crc32(t.encode())) % 2**32)
    conf_idx = np.arange(1000)
    steps_null = np.empty((V, Gt, PERMS_))
    for kperm in range(PERMS_):
        binp = binof.copy()
        rng.shuffle(binp)
        Pp, _ = weighted_freqs(Gg, W, binp, bins, min_w=MIN_W)
        for g in range(Gt):
            others = [o for o in range(Gt) if o != g]
            num = np.nansum(np.stack([Pp[:, o, :] * N[o, :][None, :] for o in others], 0), 0)
            den = np.nansum(np.stack([np.where(np.isfinite(Pp[:, o, :]), N[o, :][None, :], 0) for o in others], 0), 0)
            po = num / np.maximum(den, 1e-9)
            Lp = logit(Pp[:, g, :]) - logit(np.clip(po, EPS, 1 - EPS))
            d = np.abs(np.diff(Lp, axis=1))
            with np.errstate(invalid='ignore'):
                steps_null[:, g, kperm] = np.nanmax(d, axis=1)
    pair_names = {s: f"{a}|{b}" for s, (a, b) in enumerate(pair_columns(K))}
    out = []
    ari = float(qc.loc[t, 'ari'])
    for k, (i, pos, vid, ref, alt) in enumerate(rows):
        for g in range(Gt):
            lor = LOR[:, g, :]
            d = np.abs(np.diff(lor[k]))
            with np.errstate(invalid='ignore'):
                st = np.nanmax(d) if np.any(np.isfinite(d)) else np.nan
            fin = np.where(np.isfinite(lor[k]))[0]
            e0 = lor[k, fin[0]] if len(fin) else np.nan
            if not np.isfinite(st):
                continue
            nperm = int(np.nansum(steps_null[k, g, :] >= st))
            out.append(dict(window=t, chrom=win_by_tag[t]['chrom'], pos=pos, vid=vid, ref=ref, alt=alt,
                            akbari_X=xmap.get((win_by_tag[t]['chrom'], pos), np.nan), ari=ari, ari_seed_b=ari_b,
                            group=g, group_name=(pair_names.get(top[g], str(top[g])) if g < len(top) else 'OTHER'),
                            step=float(st), e0=float(e0), perm_p=(nperm + 1) / (PERMS_ + 1),
                            L_traj=_json.dumps([round(x, 3) for x in lor[k]]),
                            p_traj=_json.dumps([round(x, 4) for x in P[k, g, :]]),
                            f_curve=_json.dumps([round(x, 3) for x in fcurves[g]]),
                            n_states_conf=Gt))
    return out


if __name__ == '__main__':
    main()
