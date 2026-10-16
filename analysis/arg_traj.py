#!/usr/bin/env python3
"""Frequency-trajectory estimator on dated ARGs (docs/arg-trajectory-plan.md, item 5).

For a variant with derived-allele carriers, compute the carrier-lineage fraction on
the dated tree sequence:

    rho(T) = (# edges present at time T whose child is in the carrier subtree)
             / (# edges present at time T)

This is the Relate-style ARG frequency proxy. Reported per variant: present frequency,
T_rise (oldest time with rho >= 50% of present rho), width (T at 25% minus T at 75% of
present rho), and the rho(T) curve. Small width + recent T_rise = punctuated rise.

Usage:
  arg_traj.py <chrom> <trees-path> <label:pos[,...]> [window-lo window-hi]
"""
import sys
from pathlib import Path

import numpy as np
import tskit

GRID = np.concatenate([np.arange(2, 200, 2.0), np.logspace(np.log10(200), np.log10(2e6), 60)])


def local_tree_edges(ts, pos):
    """Edges of the LOCAL tree at position `pos`: (child, parent, t_child, t_parent)."""
    tree = ts.at(pos)
    par = tree.parent_array
    t = ts.tables.nodes.time
    kids = np.where(par != -1)[0]
    return kids, par[kids], t[kids], t[par[kids]]


def subtree_from_par(par, root, n):
    """Nodes in the subtree of `root` given the local parent array (BFS over children)."""
    children = {}
    for v in np.where(par != -1)[0]:
        children.setdefault(par[v], []).append(v)
    seen = {root}
    stack = [root]
    while stack:
        u = stack.pop()
        for v in children.get(u, ()):
            if v not in seen:
                seen.add(v); stack.append(v)
    return seen


def trajectory(ts, mut_node, pos, grid=GRID):
    kids, pars, tc, tp = local_tree_edges(ts, pos)
    sub = subtree_from_par(ts.at(pos).parent_array, mut_node, ts.num_nodes)
    mask = np.array([k in sub for k in kids])
    rho = np.zeros(len(grid)); nlin = np.zeros(len(grid), dtype=int)
    for i, T in enumerate(grid):
        active = (tc <= T) & (T < tp)
        tot = int(active.sum()); nlin[i] = tot
        rho[i] = (active & mask).sum() / tot if tot else np.nan
    return grid, rho, nlin


def analyse(ts, pos, label=None):
    site = ts.site(position=pos)
    muts = [m for m in ts.mutations() if m.site == site.id]
    if not muts:
        return None
    mut_nodes = [m.node for m in muts]
    t = ts.tables.nodes.time
    # present frequency from genotypes at that site
    G = ts.genotype_matrix()
    row = G[site.id]
    freq = float(row.sum() / len(row))
    # trajectory: union of carrier subtrees (use the oldest mutation node's subtree)
    node = max(mut_nodes, key=lambda n: t[n])
    grid, rho, n_lin = trajectory(ts, node, pos)
    pres = np.nanmax(rho) if np.isfinite(rho).any() else np.nan
    out = dict(label=label or str(pos), pos=int(pos), freq=freq,
               origin_gen=float(t[node]), max_rho=float(pres))
    if not np.isfinite(pres) or pres <= 0:
        return out
    def cross(frac):
        """Oldest time at which rho is still at or above `frac` of its present value."""
        i2 = np.where(rho >= frac * pres)[0]
        return float(grid[i2[-1]]) if len(i2) else np.nan
    out['T_rise'] = cross(0.5)
    out['T25'] = cross(0.25); out['T75'] = cross(0.75)
    out['width'] = (out['T25'] - out['T75']) if np.isfinite(out['T25']) and np.isfinite(out['T75']) else np.nan
    out['curve'] = [(round(float(g), 1), round(float(r), 4)) for g, r in zip(grid, rho) if np.isfinite(r) and r > 0]
    return out


def main():
    chrom, path = sys.argv[1], sys.argv[2]
    targets = [(s.split(':')[0], int(s.split(':')[1])) for s in sys.argv[3].split(',')]
    ts = tskit.load(path)
    lo = None
    if len(sys.argv) > 4:
        lo = int(sys.argv[4])           # positions in the file are offset by lo after trim()
    print(f'ARG: sites {ts.num_sites} trees {ts.num_trees} samples {ts.num_samples}')
    for label, p in targets:
        use = p - lo if lo is not None else p
        r = analyse(ts, use, label)
        if r is None:
            print(f'{label}: no mutation at {p}'); continue
        print(f"{r['label']}: pos={r['pos']} freq={r['freq']:.3f} origin={r['origin_gen']:.0f} gen "
              f"max_rho={r['max_rho']:.3f} T_rise={r.get('T_rise', float('nan')):.0f} "
              f"T25={r.get('T25', float('nan')):.0f} T75={r.get('T75', float('nan')):.0f} "
              f"width={r.get('width', float('nan')):.0f} gen")


if __name__ == '__main__':
    main()
