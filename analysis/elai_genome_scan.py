#!/usr/bin/env python3
"""Genome-wide local-ancestry scan with the validated configuration.

Configuration validated on the AMY region probe (docs/local-ancestry-verdict.md,
REVERSAL section): S=2 sources, K=10 lower clusters, 1 Mb windows, 100 exact EM
steps, two seeds, unordered diploid state posterior q at the window centre via
the patched --pkk-site export.

Modes
  build     stream the Akbari genome-wide pgen once; same marker QC as the dense
            build (biallelic ACGT, call>=.95, MAF>=.05); split autosomes into
            1 Mb windows capped at ~800 markers by uniform spacing; buffer all
            kept genotypes in RAM (~2.3 GB uint8) and write per-window ELAI
            inputs under data/derived/elai_genome/windows/.
  window    run one window (both seeds), reduce the pkk-site stream to the
            per-individual unordered q matrix + diagnostics, write JSON.
  scan      run every window in parallel (xargs -P).
  assemble  aggregate per-window outputs into genome-wide arrays + QC table.
  compose   composition check: confident states vs date/country labels.

Writes only under data/derived/elai_genome/. Reads data/raw/ read-only.
"""
import argparse, json, subprocess, sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, iter_pvar, paths  # noqa: E402

GENOME = ROOT / 'data/derived/elai_genome'
BINARY = '/home/erikg/elai-fork/src/elai'
S, K, STEPS, MG = 2, 10, 100, 10
SEEDS = (20261020, 20261021)
WINDOW = 1_000_000
MAX_MARKERS = 800
MIN_MARKERS = 150
CORE = (103915000, 104327000)  # AMY1 segmental-duplication core (validated exclusion)
DENSE_SAMPLES = ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'
PLOIDY_MAP = {3: 2, 2: 1, 1: 0}  # will re-encode to ELAI's 0/1/2 with 9 = missing


def select_samples():
    return json.load(open(DENSE_SAMPLES))


def marker_qc(h):
    valid = (h >= 0) & (h <= 2)
    calls = valid.mean(1)
    af = np.divide(np.where(valid, h, 0).sum(1), 2 * valid.sum(1),
                   out=np.zeros(len(h), dtype=float), where=valid.sum(1) > 0)
    maf = np.minimum(af, 1 - af)
    return (calls >= .95) & (maf >= .05), maf


def mode_build(_args):
    samples = select_samples()
    idx = np.array([s['index'] for s in samples], dtype=np.uint32)
    n = len(samples)
    win_dir = GENOME / 'windows'
    win_dir.mkdir(parents=True, exist_ok=True)

    kept_meta = []          # (chrom, pos, vid, ref, alt)
    kept_geno = []          # uint8 rows (n,), 0/1/2, 9 = missing
    batch = []              # (variant_index, chrom, pos, vid, ref, alt)
    counts = dict(scanned=0, qc_pass=0, core_excluded=0)

    def flush(reader):
        if not batch:
            return
        h = reader.read_hardcalls(np.array([r[0] for r in batch], dtype=np.uint32))
        mask, maf = marker_qc(h)
        for j, r in enumerate(batch):
            if not mask[j]:
                continue
            counts['qc_pass'] += 1
            c, pos = r[1], r[2]
            if c == 1 and CORE[0] <= pos <= CORE[1]:
                counts['core_excluded'] += 1
                continue
            g = h[j].astype(np.uint8)
            g[~((g >= 0) & (g <= 2))] = 9
            kept_meta.append((c, pos, r[3], r[4], r[5]))
            kept_geno.append(g)
        batch.clear()

    with AkbariGeno(sample_subset=idx) as reader, open(GENOME / 'build.progress', 'w') as prog:
        nvar = reader.n_variants
        for i, (chrom, pos, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
            counts['scanned'] += 1
            if not chrom.isdigit() or not 1 <= int(chrom) <= 22:
                continue
            if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
                continue
            batch.append((i, int(chrom), pos, vid, ref, alt))
            if len(batch) >= 2048:
                flush(reader)
            if i % 500000 == 0:
                prog.write(f'{i}/{nvar} kept={len(kept_meta)}\n'); prog.flush()
        flush(reader)

    G = np.asarray(kept_geno, dtype=np.uint8)          # (M, n): 0/1/2, 9 = missing
    meta = kept_meta
    by_win = defaultdict(list)                          # (chrom, bin) -> row indices
    for j, (c, pos, *_rest) in enumerate(meta):
        if c == 1 and CORE[0] <= pos <= CORE[1]:
            continue
        by_win[(c, pos // WINDOW)].append(j)

    manifest = dict(samples=[s['sample'] for s in samples], n=n, counts=counts,
                    config=dict(S=S, K=K, steps=STEPS, mg=MG, seeds=SEEDS,
                                window=WINDOW, max_markers=MAX_MARKERS),
                    windows=[], skipped=[])
    for (c, wbin), rows in sorted(by_win.items()):
        rows.sort(key=lambda j: meta[j][1])
        span = meta[rows[-1]][1] - meta[rows[0]][1]
        spacing = max(100, span // MAX_MARKERS)
        sel, last = [], -10**12
        for j in rows:
            if meta[j][1] - last < spacing:
                continue
            sel.append(j); last = meta[j][1]
        if len(sel) < MIN_MARKERS:
            manifest['skipped'].append(dict(chrom=c, bin=wbin, markers=len(sel)))
            continue
        tag = f'c{c}_w{wbin}'
        wd = win_dir / tag
        wd.mkdir(exist_ok=True)
        center = wbin * WINDOW + WINDOW // 2
        with (wd / 'genotypes.txt').open('w') as g:
            for j in sel:
                vid, ref, alt = meta[j][2], meta[j][3], meta[j][4]
                codes = np.array([ref + ref, ref + alt, alt + alt, '??'])
                g.write(vid + ' ' + ' '.join(codes[G[j]]) + '\n')
        with (wd / 'panel.inp').open('w') as f:
            f.write(f"{n}\n{len(sel)}\nIND " + ' '.join(s['sample'] for s in samples) + '\n')
            with (wd / 'genotypes.txt').open() as src:
                f.write(src.read())
        with (wd / 'panel.pos').open('w') as f:
            for j in sel:
                f.write(f"{meta[j][2]} {meta[j][1]} {c}\n")
        manifest['windows'].append(dict(tag=tag, chrom=c, bin=wbin, markers=len(sel),
                                        center=center, first=meta[sel[0]][1], last=meta[sel[-1]][1]))
    json.dump(manifest, open(GENOME / 'build.manifest.json', 'w'))
    print(json.dumps(dict(counts=counts, windows=len(manifest['windows']),
                          skipped=len(manifest['skipped']))))


def consume_pkk_site(path, k2):
    """Reduce a pkk-site.tsv stream to per-individual unordered q (n x k2) and ordered p."""
    import pandas as pd
    d = pd.read_csv(path, sep='\t')
    if not abs(d.groupby('individual_index').p_ordered.sum().median() - 1.0) < 1e-6:
        raise RuntimeError('ordered posterior not normalised')
    n = d.individual_index.max() + 1
    # q_unordered is defined on unordered pairs {a,b}; take rows with k1<=k2 and
    # key by the (min,max) pair so the column space is exactly K(K+1)/2.
    d = d[d.k1 <= d.k2].copy()
    d['u'] = d.k1.astype(str) + '|' + d.k2.astype(str)
    q = d.groupby(['individual_index', 'u']).q_unordered.sum().unstack(fill_value=0.0)
    kcols = K * (K + 1) // 2
    Q = np.zeros((n, kcols))
    cols = {}
    for c in q.columns:
        a, b = c.split('|')
        i = int(a) * K + int(b) - int(a) * (int(a) + 1) // 2
        cols[c] = i
        Q[:, i] = q[c].to_numpy()
    hard = q.to_numpy().argmax(1)
    mx = q.to_numpy().max(1)
    ent = -(q.to_numpy() / q.to_numpy().sum(1, keepdims=True) *
            np.log(np.clip(q.to_numpy() / q.to_numpy().sum(1, keepdims=True), 1e-300, None))).sum(1)
    return Q, q.columns.tolist(), hard, mx, ent


def mode_window(args):
    m = json.load(open(GENOME / 'build.manifest.json'))
    win = {w['tag']: w for w in m['windows']}[args.window]
    wd = GENOME / 'windows' / args.window
    out = dict(tag=args.window, chrom=win['chrom'], center=win['center'], markers=win['markers'], seeds={})
    hards = {}
    for seed in SEEDS:
        prefix = f'w_s{S}_k{K}_seed{seed}'
        cmd = [BINARY, '-g', 'panel.inp', '-p', '1', '-pos', 'panel.pos', '-C', str(S), '-c', str(K),
               '-w', '0', '-s', str(STEPS), '-mg', str(MG), '-R', str(seed), '-o', prefix,
               '--pkk-site', f"{win['chrom']}:{win['center']}"]
        log = (wd / f'{prefix}.launch.log').open('w')
        proc = subprocess.run(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
        log.close()
        ok = proc.returncode == 1 and '-elai: finished' in (wd / f'{prefix}.launch.log').read_text()
        if not ok:
            out['seeds'][str(seed)] = dict(ok=False, returncode=proc.returncode)
            continue
        Q, cols, hard, mx, ent = consume_pkk_site(wd / 'output' / f'{prefix}.pkk-site.tsv', K * (K + 1) // 2)
        np.save(wd / f'{prefix}.q.npy', Q.astype(np.float32))
        hards[seed] = hard
        out['seeds'][str(seed)] = dict(ok=True, maxq_median=float(np.median(mx)),
                                       frac_gt_90=float((mx > .9).mean()),
                                       entropy_median=float(np.median(ent)))
    if len(hards) == 2:
        a, b = hards[SEEDS[0]], hards[SEEDS[1]]
        from sklearn.metrics import adjusted_rand_score
        from scipy.optimize import linear_sum_assignment
        keys = sorted(set(a.tolist()) | set(b.tolist()))
        ix = {k: i for i, k in enumerate(keys)}
        O = np.zeros((len(keys), len(keys)))
        for x, y in zip(a, b):
            O[ix[x], ix[y]] += 1
        ra, ca = linear_sum_assignment(-O)
        best = O[ra, ca].sum() / len(a)
        va = np.bincount([ix[x] for x in a], minlength=len(keys)) / len(a)
        vb = np.bincount([ix[y] for y in b], minlength=len(keys)) / len(b)
        P = np.outer(va, vb)
        rp, cp = linear_sum_assignment(-P)
        ch = P[rp, cp].sum()
        out['cross_seed'] = dict(ari=float(adjusted_rand_score(a, b)), agree=float(best),
                                 chance=float(ch), excess=float(best - ch))
    json.dump(out, open(wd / 'window.summary.json', 'w'))
    print(json.dumps(out))


def mode_scan(args):
    m = json.load(open(GENOME / 'build.manifest.json'))
    tags = [w['tag'] for w in m['windows']]
    done = {p.parent.parent.name for p in GENOME.glob('windows/*/window.summary.json')}
    todo = [t for t in tags if t not in done]
    listfile = GENOME / 'scan.list'
    listfile.write_text('\n'.join(todo) + '\n')
    cmd = (f"xargs -a {listfile} -P {args.parallel} -I{{}} "
           f"{sys.executable} -u {__file__} window --window {{}}\n")
    print(f'{len(todo)} windows to run at -P {args.parallel}')
    subprocess.run(cmd, shell=True, check=False)


def mode_assemble(_args):
    m = json.load(open(GENOME / 'build.manifest.json'))
    rows = []
    for w in m['windows']:
        f = GENOME / 'windows' / w['tag'] / 'window.summary.json'
        if not f.exists():
            continue
        d = json.load(open(f))
        s0, s1 = (d['seeds'].get(str(SEEDS[0]), {}), d['seeds'].get(str(SEEDS[1]), {}))
        rows.append(dict(tag=w['tag'], chrom=w['chrom'], center=w['center'], markers=w['markers'],
                         ok=bool(s0.get('ok') and s1.get('ok')),
                         maxq=float(np.mean([s0.get('maxq_median', np.nan), s1.get('maxq_median', np.nan)])),
                         frac90=float(np.mean([s0.get('frac_gt_90', np.nan), s1.get('frac_gt_90', np.nan)])),
                         ari=d.get('cross_seed', {}).get('ari', np.nan),
                         excess=d.get('cross_seed', {}).get('excess', np.nan)))
    import pandas as pd
    df = pd.DataFrame(rows).dropna(subset=['maxq'])
    ok = df[df.ok]
    stats = dict(windows_total=len(m['windows']), windows_done=len(df), windows_both_seeds=int(len(ok)),
                 median_maxq=float(ok.maxq.median()), frac90_median=float(ok.frac90.median()),
                 frac_windows_q90_above_half=float((ok.frac90 > .5).mean()),
                 ari_median=float(ok.ari.median()),
                 frac_windows_ari_above_0p5=float((ok.ari > .5).mean()),
                 frac_windows_excess_positive=float((ok.excess > 0).mean()))
    df.to_csv(GENOME / 'qc_table.csv', index=False)
    json.dump(stats, open(GENOME / 'qc.summary.json', 'w'), indent=1)
    print(json.dumps(stats, indent=1))


def mode_compose(_args):
    """Composition: do confident states reproduce population structure (country/date)?"""
    import pandas as pd
    from sklearn.metrics import normalized_mutual_info_score
    m = json.load(open(GENOME / 'build.manifest.json'))
    samples = m['samples']
    countries = [s['country'] for s in select_samples()]
    per_win = []
    for w in m['windows']:
        f = GENOME / 'windows' / w['tag'] / f'w_s{S}_k{K}_seed{SEEDS[0]}.q.npy'
        g = GENOME / 'windows' / w['tag'] / 'window.summary.json'
        if not (f.exists() and g.exists()):
            continue
        d = json.load(open(g))
        if d.get('cross_seed', {}).get('ari', 0) < .5:
            continue
        Q = np.load(f)
        hard = Q.argmax(1)
        conf = Q.max(1) > .9
        if conf.sum() < 50:
            continue
        nmi_c = normalized_mutual_info_score(np.array(countries)[conf], hard[conf])
        per_win.append(dict(tag=w['tag'], chrom=w['chrom'], center=w['center'], nmi_country=float(nmi_c),
                            n_conf=int(conf.sum())))
    df = pd.DataFrame(per_win)
    df.to_csv(GENOME / 'composition_table.csv', index=False)
    stats = dict(windows_scored=len(df), median_nmi_country=float(df.nmi_country.median()) if len(df) else None,
                 frac_windows_nmi_above_0p1=float((df.nmi_country > .1).mean()) if len(df) else None)
    json.dump(stats, open(GENOME / 'composition.summary.json', 'w'), indent=1)
    print(json.dumps(stats, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['build', 'window', 'scan', 'assemble', 'compose'])
    ap.add_argument('--window')
    ap.add_argument('--parallel', type=int, default=240)
    a = ap.parse_args()
    dict(build=mode_build, window=mode_window, scan=mode_scan,
         assemble=mode_assemble, compose=mode_compose)[a.mode](a)


if __name__ == '__main__':
    main()
