#!/usr/bin/env python3
"""Build the preregistered bounded ELAI panel, streaming immutable pvar/pgen."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from load_akbari_geno import AkbariGeno, iter_pvar, parse_psam, paths

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/elai_amy'
SEED = 20261010
CORE = (103915000, 104327000)
SITE = 104121000


def excluded(pos):
    return CORE[0] <= pos <= CORE[1]


def cohort():
    ann = ROOT / 'data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.anno'
    ind = ROOT / 'data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.ind'
    with ind.open() as fh:
        groups = {f[0]: f[2] for line in fh if len(f := line.split()) == 3}
    with ann.open() as fh:
        rows = list(csv.reader(fh, delimiter='\t'))[1:]
    lookup = {r[0]: r for r in rows}
    eligible, seen = [], set()
    for i, sid in enumerate(parse_psam(paths()[2])):
        r = lookup.get(sid)
        if r is None or sid not in groups:
            continue
        try:
            bp = float(r[10])
        except ValueError:
            continue
        if not np.isfinite(bp) or bp <= 0 or r[2] in seen:
            continue
        seen.add(r[2])
        eligible.append(dict(index=i, sample=sid, biological_id=r[2], BP=bp,
                             group=groups[sid], anno_group=r[14], country=r[16]))
    chosen = np.sort(np.random.default_rng(SEED).choice(len(eligible), min(128, len(eligible)), replace=False))
    return [eligible[i] for i in chosen], len(eligible)


def candidates(pvar):
    local, pcs = [], []
    last = {1: -200000, 2: -200000}
    scanned = 0
    for i, (chrom, pos, vid, ref, alt) in enumerate(iter_pvar(pvar)):
        scanned = i + 1
        if chrom not in ('1', '2'):
            if chrom.isdigit() and int(chrom) > 2:
                break
            continue
        c = int(chrom)
        if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
            continue
        spacing = 2000 if c == 1 else 200000
        if c == 1 and (not 95000000 <= pos <= 115000000 or excluded(pos)):
            continue
        if pos - last[c] < spacing:
            continue
        last[c] = pos
        (local if c == 1 else pcs).append((i, c, pos, vid, ref, alt))
    return local, pcs, scanned


def filtered(reader, variants):
    selected, hardcalls, stats = [], [], []
    for start in range(0, len(variants), 256):
        batch = variants[start:start + 256]
        hard = reader.read_hardcalls(np.array([r[0] for r in batch], dtype=np.uint32))
        for row, h in zip(batch, hard):
            valid = (h >= 0) & (h <= 2)
            call = float(valid.mean())
            af = float(h[valid].mean() / 2) if valid.any() else 0
            maf = min(af, 1 - af)
            if call >= .95 and maf >= .05:
                selected.append(row)
                hardcalls.append(h.copy())
                stats.append((maf, call))
    return selected, np.array(hardcalls, dtype=np.int8), stats


def write_panel(name, rows, h, ids):
    with (OUT / f'{name}.inp').open('w') as f, (OUT / f'{name}.pos').open('w') as p:
        f.write(f'{len(ids)}\n{len(rows)}\nIND ' + ' '.join(ids) + '\n')
        for row, calls in zip(rows, h):
            _, chrom, pos, vid, ref, alt = row
            codes = {0: ref + ref, 1: ref + alt, 2: alt + alt}
            f.write(vid + ' ' + ' '.join(codes.get(int(x), '??') for x in calls) + '\n')
            p.write(f'{vid} {pos} {chrom}\n')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    samples, eligible = cohort()
    with (OUT / 'samples.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(samples[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(samples)
    local, pcrows, scanned = candidates(paths()[1])
    with AkbariGeno(sample_subset=np.array([s['index'] for s in samples], dtype=np.uint32)) as reader:
        rows, h, stats = filtered(reader, local)
        pcrows, ph, _ = filtered(reader, pcrows)
        dimensions = [reader.n_variants, reader._raw_n_samples]
    assert all(not excluded(r[2]) for r in rows)
    assert excluded(SITE)
    pos = np.array([r[2] for r in rows])
    assert np.all(np.diff(pos) > 0) and np.any(pos < SITE) and np.any(pos > SITE)
    assert all(r[1] == 2 for r in pcrows)
    ids = [s['sample'] for s in samples]
    write_panel('panel', rows, h, ids)
    write_panel('thin', rows[::2], h[::2], ids)
    with (OUT / 'markers.tsv').open('w') as f:
        f.write('pvar_index\tchrom\tposition\tid\tref\talt\tmaf\tcall_rate\n')
        for r, st in zip(rows, stats):
            f.write('\t'.join(map(str, (*r, *st))) + '\n')
    with (OUT / 'pc_markers.tsv').open('w') as f:
        f.write('pvar_index\tchrom\tposition\tid\tref\talt\n')
        for r in pcrows:
            f.write('\t'.join(map(str, r)) + '\n')
    x = ph.astype(float).T
    x[(x < 0) | (x > 2)] = np.nan
    mean = np.nanmean(x, axis=0)
    x = np.where(np.isfinite(x), x, mean)
    freq = mean / 2
    x = (x - mean) / np.sqrt(2 * freq * (1 - freq))
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    pcs = u[:, :5] * s[:5]
    np.savetxt(OUT / 'independent_pcs.tsv', pcs, delimiter='\t', header='PC1\tPC2\tPC3\tPC4\tPC5', comments='')
    coverage = ((h >= 0) & (h <= 2)).mean(axis=0)
    manifest = dict(seed=SEED, source_dimensions=dimensions, eligible_dated_unique=eligible,
                    fitted_n=len(samples), dated_n=len(samples), candidates=len(local), markers=len(rows),
                    thin_markers=len(rows[::2]), pc_markers=len(pcrows), pvar_records_scanned=scanned,
                    min_marker_maf=min(st[0] for st in stats), min_marker_call_rate=min(st[1] for st in stats),
                    individual_call_rates=coverage.tolist(), usable_individuals=int((coverage >= .95).sum()),
                    excluded_markers=sum(excluded(r[2]) for r in rows), site_inside_core=excluded(SITE),
                    left_flank=int(pos[pos < SITE][-1]), right_flank=int(pos[pos > SITE][0]),
                    pc_overlap=0, prereg_sha256=hashlib.sha256((ROOT / 'docs/elai-amy-prereg.md').read_bytes()).hexdigest())
    (OUT / 'panel.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({k: v for k, v in manifest.items() if k != 'individual_call_rates'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
