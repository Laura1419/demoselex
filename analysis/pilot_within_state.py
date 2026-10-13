#!/usr/bin/env python3
"""Within-state temporal pilot (prereg: docs/pilot-within-state-prereg.md).

For each pilot window's variants: hard-assign individuals to their confident
unordered local-ancestry state (max q > 0.9, seed-A q matrix), split dates into
quartile bins per state, WLS slope of allele frequency on calendar year per
state, inverse-variance combine, weighted sign concordance. Joins Akbari X by
chrom/pos where available. Writes CSV + summary under data/derived/pilot_within_state/.

Read-only inputs; writes only under data/derived/pilot_within_state/.
"""
import csv, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, iter_pvar, paths  # noqa: E402

GENOME = ROOT / 'data/derived/elai_genome'
OUT = ROOT / 'data/derived/pilot_within_state'
S, K = 2, 10
SEED_A = 20261020
MIN_STATE_N = 60
MIN_BIN_N = 10
MAX_VARIANTS = 600
CONF = 0.9
RNG_SEED = 20261020

KNOWN = {'LCT': 'c2_w135', 'FADS': 'c11_w61', 'SLC45A2': 'c5_w33', 'SLC24A5': 'c15_w48'}
AMY = ['c1_w103', 'c1_w104', 'c1_w105']


def pilot_windows():
    qc = pd.read_csv(GENOME / 'qc_table.csv')
    gated = qc[(qc.ok) & (qc.ari > .5) & (qc.frac90 > .3)]
    exclude = set(KNOWN.values()) | set(AMY) | {w['tag'] for w in
                                                json.load(open(GENOME / 'build.manifest.json'))['windows']
                                                if w['chrom'] == 1 and 100 <= w['bin'] <= 108}
    pool = gated[~gated.tag.isin(exclude)]
    rng = np.random.default_rng(RNG_SEED)
    picks = rng.choice(len(pool), size=24, replace=False)
    return sorted(KNOWN.items(), key=lambda kv: kv[1]), AMY, sorted(pool.tag.iloc[picks])


def collect_variant_indices(tags, manifest):
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
    # uniform spacing cap per window
    for t, rows in found.items():
        rows.sort()
        if len(rows) > MAX_VARIANTS:
            step = len(rows) / MAX_VARIANTS
            rows = [rows[int(k * step)] for k in range(MAX_VARIANTS)]
        found[t] = rows
    return found


def within_state_stats(g, dates, state, n_states_space=55):
    """g: (V,n) uint8 0/1/2; dates: (n,) calendar year AD; state: (n,) int or -1."""
    V = g.shape[0]
    rows = []
    states = np.unique(state[state >= 0])
    for v in range(V):
        gv = g[v]
        called = gv <= 2
        zs, slopes, ses, ns = [], [], [], []
        for s in states:
            m = (state == s) & called
            n_s = int(m.sum())
            if n_s < MIN_STATE_N:
                continue
            t = dates[m]; y = gv[m].astype(float)
            order = np.argsort(t)
            qs = np.array_split(order, 4)
            if min(len(q) for q in qs) < MIN_BIN_N:
                continue
            xs, ps, ses_b = [], [], []
            ok = True
            for q in qs:
                nv = int((gv[m][q] <= 2).sum())
                if nv < 4:
                    ok = False; break
                p = gv[m][q][gv[m][q] <= 2].sum() / nv
                xs.append(t[q].mean()); ps.append(p)
                ses_b.append(max(np.sqrt(max(p * (1 - p), 1e-12) / nv), 1e-6))
            if not ok:
                continue
            xs, ps, w = np.array(xs), np.array(ps), 1 / np.square(np.array(ses_b))
            A = np.vstack([xs, np.ones_like(xs)]).T
            W = np.diag(w)
            coef = np.linalg.solve(A.T @ W @ A, A.T @ W @ ps)
            cov = np.linalg.inv(A.T @ W @ A)
            slope, intercept = coef
            se_slope = np.sqrt(cov[0, 0])
            zs.append(slope / se_slope); slopes.append(slope); ses.append(se_slope); ns.append(n_s)
        if not zs:
            rows.append(dict(n_states_used=0, z_comb=np.nan, concordance=np.nan,
                             slopes=json.dumps([]), zs=json.dumps([]), ns=json.dumps([])))
            continue
        w = 1 / np.square(np.array(ses))
        slope_comb = float(np.dot(w, slopes) / w.sum())
        se_comb = float(1 / np.sqrt(w.sum()))
        zc = slope_comb / se_comb
        dom = np.sign(np.median(np.array(zs)[np.abs(zs) >= 2])) if any(np.abs(np.array(zs)) >= 2) else np.sign(np.median(zs))
        conc = float(np.dot(w, np.sign(zs) == dom) / w.sum())
        rows.append(dict(n_states_used=len(zs), z_comb=float(zc), concordance=conc,
                         slopes=json.dumps([round(x, 6) for x in slopes]),
                         zs=json.dumps([round(x, 3) for x in zs]), ns=json.dumps(ns)))
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.load(open(GENOME / 'build.manifest.json'))
    global win_by_tag
    win_by_tag = {w['tag']: w for w in manifest['windows']}
    known, amy, nulls = pilot_windows()
    tags = [t for _, t in known] + amy + nulls
    missing = [t for t in tags if t not in win_by_tag]
    if missing:
        print('MISSING WINDOWS (skipped):', missing); tags = [t for t in tags if t in win_by_tag]
    samples = manifest['samples']
    dates = np.array([1950 - s['BP'] for s in json.load(open(ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'))], dtype=float)

    print('collecting variant indices...', flush=True)
    found = collect_variant_indices(tags, None)
    total = sum(len(v) for v in found.values())
    print(f'{total} variants across {len(found)} windows', flush=True)

    idx_all = sorted({i for rows in found.values() for i, *_ in rows})
    pos_of = {}
    for t, rows in found.items():
        for i, pos, vid, ref, alt in rows:
            pos_of[i] = (t, pos, vid, ref, alt)

    # X join
    xmap = {}
    for fn in ('akbari_top.tsv', 'akbari_signal_absXge4.tsv', 'akbari_background200k.tsv'):
        p = ROOT / 'data/derived' / fn
        if not p.exists():
            continue
        df = pd.read_csv(p, sep='\t', usecols=['CHROM', 'POS', 'X'])
        for c, pos, x in zip(df.CHROM, df.POS, df.X):
            xmap[(int(c), int(pos))] = float(x)

    qmat = {}
    for t in tags:
        f = GENOME / 'windows' / t / f'w_s{S}_k{K}_seed{SEED_A}.q.npy'
        if f.exists():
            qmat[t] = np.load(f)

    qc = pd.read_csv(GENOME / 'qc_table.csv').set_index('tag')
    out_rows = []
    sample_idx = np.array([s['index'] for s in json.load(open(ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'))], dtype=np.uint32)
    with AkbariGeno(sample_subset=sample_idx) as reader:
        for t in tags:
            rows = found.get(t, [])
            if not rows or t not in qmat:
                continue
            Q = qmat[t]
            hard = Q.argmax(1)
            conf = Q.max(1) > CONF
            state = np.where(conf, hard, -1)
            vi = np.array([r[0] for r in rows], dtype=np.uint32)
            g = reader.read_hardcalls(vi) if hasattr(reader, 'read_hardcalls') else None
            if g is None:
                raise RuntimeError('read_hardcalls unavailable')
            g = np.where(g > 2, 9, g)
            stats = within_state_stats(g.astype(float), dates, state)
            ari = float(qc.loc[t, 'ari']) if t in qc.index else np.nan
            for (i, pos, vid, ref, alt), st in zip(rows, stats):
                out_rows.append(dict(window=t, chrom=win_by_tag[t]['chrom'], pos=pos, vid=vid,
                                     ref=ref, alt=alt, ari=ari,
                                     akbari_X=xmap.get((win_by_tag[t]['chrom'], pos), np.nan), **st))
            print(f'{t}: {len(rows)} variants, states used median '
                  f'{np.median([s["n_states_used"] for s in stats]):.0f}', flush=True)

    df = pd.DataFrame(out_rows)
    df['arm'] = np.where(df.window.isin(KNOWN.values()), 'known',
                         np.where(df.window.isin(AMY), 'amy', 'null'))
    df.to_csv(OUT / 'variants.csv', index=False)

    summary = {}
    null = df[(df.arm == 'null') & (df.n_states_used > 0)]
    known = df[(df.arm == 'known') & (df.n_states_used > 0)]
    amy = df[(df.arm == 'amy') & (df.n_states_used > 0)]
    summary['null'] = dict(n=len(null), z_median=float(null.z_comb.median()),
                           z_q999=float(null.z_comb.quantile(.999)),
                           z_q99=float(null.z_comb.quantile(.99)),
                           abs_z_q999=float(null.z_comb.abs().quantile(.999)))
    summary['known'] = {t: dict(x_max=float(g.akbari_X.max()), z_max=float(g.z_comb.abs().max()),
                                z_at_max_x=float(g.loc[g.akbari_X.idxmax(), 'z_comb']) if np.isfinite(g.akbari_X.max()) else np.nan,
                                conc_at_max_x=float(g.loc[g.akbari_X.idxmax(), 'concordance']) if np.isfinite(g.akbari_X.max()) else np.nan)
                        for t, g in known.groupby('window')}
    summary['amy'] = {t: dict(z_max=float(g.z_comb.abs().max()),
                              z_at_max_x=float(g.loc[g.akbari_X.idxmax(), 'z_comb']) if np.isfinite(g.akbari_X.max()) else np.nan,
                              n_above_q999=int((g.z_comb.abs() > summary['null']['abs_z_q999']).sum()),
                              n_variants=len(g)) for t, g in amy.groupby('window')}
    both = df[np.isfinite(df.akbari_X) & np.isfinite(df.z_comb) & (df.n_states_used > 0)]
    summary['x_vs_z_spearman'] = float(both.akbari_X.corr(both.z_comb.abs(), method='spearman'))
    summary['x_vs_z_signed_spearman'] = float(both.akbari_X.corr(both.z_comb, method='spearman'))
    summary['join_coverage'] = dict(n=len(both), frac_with_X=float(len(both) / max(len(df[df.n_states_used > 0]), 1)))
    json.dump(summary, open(OUT / 'summary.json', 'w'), indent=1, default=str)
    print(json.dumps(summary, indent=1, default=str))


if __name__ == '__main__':
    main()
