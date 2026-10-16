#!/usr/bin/env python3
"""2D comparison: our ARG based z (clade youth) vs Akbari et al. 2026 significance.

x = our z  (negative = carrier clade younger than expected for its frequency x rho)
y = -log10(Akbari P), where P is recovered from their published X via Z = X*sqrt(CF),
    CF = 2.28 (our reconstruction) and P = 2*Phi(-|Z|).

Quadrants (the project's core contrast):
  top-left    Akbari-significant AND young clade   -> concordant (in-situ selection)
  bottom-left young clade, not Akbari-significant  -> our-specific candidates
  top-right   Akbari-significant, OLD clade        -> consistent with carriage/demography
  bottom-right neither

Top hits in each direction are labelled. Writes fig/arg_z_vs_akbari.png and .pdf
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
CF = 2.28
X_STRONG, X_FDR50 = 5.45, 3.16


def main():
    z = pd.read_csv(ROOT / 'data/derived/arggenome/genomewide_z.csv')
    z['pct'] = z.z.rank(method='min') / len(z) * 100
    tc = pd.read_csv(ROOT / 'data/derived/two_channel.tsv', sep='\t', low_memory=False,
                     usecols=['chrom', 'pos', 'akbari_X'])
    tc = tc.dropna(subset=['akbari_X'])
    m = z.merge(tc, on=['chrom', 'pos'], how='inner')
    m['Z_ak'] = m.akbari_X.abs() * np.sqrt(CF)
    m['neglog10P'] = -np.log10(np.clip(2 * stats.norm.sf(m.Z_ak), 1e-300, None))
    m['lab'] = m.chrom.astype(str) + ':' + (m.pos / 1e6).round(2).astype(str) + 'Mb'

    lo_strong = -np.log10(2 * stats.norm.sf(X_STRONG * np.sqrt(CF)))
    lo_fdr = -np.log10(2 * stats.norm.sf(X_FDR50 * np.sqrt(CF)))

    fig, (ax, axh) = plt.subplots(1, 2, figsize=(14, 7.5), width_ratios=[4, 1],
                                  gridspec_kw=dict(wspace=0.06))
    hb = ax.hexbin(m.z, m.neglog10P, gridsize=70, mincnt=1, cmap='viridis', bins='log')
    plt.colorbar(hb, ax=ax, pad=0.01, label='variants per hex (log)')
    ax.axhline(lo_strong, color='crimson', ls='--', lw=1,
               label=f"Akbari strong |X|>5.45 (P={10**-lo_strong:.1e})")
    ax.axhline(lo_fdr, color='orange', ls=':', lw=1,
               label=f"Akbari FDR-50% |X|=3.16 (P={10**-lo_fdr:.1e})")
    ax.axvline(0, color='grey', lw=0.8)

    topz = m.nsmallest(6, 'z')
    topP = m.nlargest(6, 'neglog10P')
    box_z = "top by our z (young clade)\n" + "\n".join(
        f"  {r.lab:>14s}  z={r.z:+.2f}  age={r.age:.0f}g  X={r.akbari_X:+.1f}" for r in topz.itertuples())
    box_p = "top by Akbari P\n" + "\n".join(
        f"  {r.lab:>14s}  P={10**-r.neglog10P:.1e}  age={r.age:.0f}g  z={r.z:+.2f}" for r in topP.itertuples())
    ax.text(0.012, 0.72, box_z, transform=ax.transAxes, fontsize=8, va='top', family='monospace',
            bbox=dict(boxstyle='round,pad=0.35', fc='aliceblue', ec='blue', alpha=0.92), zorder=8)
    ax.text(0.985, 0.60, box_p, transform=ax.transAxes, fontsize=8, va='top', ha='right', family='monospace',
            bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='black', alpha=0.92), zorder=8)
    known = {'LCT': (2, 136608646), 'FADS': (11, 61585000), 'SLC45A2': (5, 33951693),
             'SLC24A5': (15, 48426484), 'AMY-best': (1, 104350432)}
    for name, (c, p) in known.items():
        r = m[(m.chrom == c) & (np.abs(m.pos - p) < 200000)]
        if len(r):
            r = r.iloc[(r.pos - p).abs().argmin()]
            ax.scatter([r.z], [r.neglog10P], s=90, facecolor='none', edgecolor='red', lw=1.6, zorder=5)
            ax.annotate(name, (r.z, r.neglog10P), fontsize=9, color='red', fontweight='bold',
                        xytext=((14, -16) if name == 'LCT' else (6, -12)), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='red', lw=0.8, alpha=0.9))
    ax.set_xlabel('our ARG z  (negative = carrier clade younger than expected)', fontsize=11)
    ax.set_ylabel(r'Akbari significance  $-\log_{10}P$', fontsize=11)
    ax.set_title(f'ARG clade-youth vs Akbari et al. 2026  (n={len(m):,} shared variants;\n'
                 f'{100*len(m)/len(z):.1f}% of our {len(z):,} scored variants carry an Akbari X)', fontsize=11)
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(alpha=0.15)
    ax.set_ylim(-3, max(m.neglog10P) * 1.18)
    ax.text(0.40, 0.93, 'BOTH agree\n(concordant selection)', transform=ax.transAxes,
            fontsize=9, color='darkgreen',
            bbox=dict(boxstyle='round,pad=0.3', fc='honeydew', ec='darkgreen', alpha=0.85))
    ax.text(0.985, 0.99, 'Akbari-only\n(often OLD clade\n-> carriage/demography)', transform=ax.transAxes,
            fontsize=8.5, color='purple', ha='right',
            bbox=dict(boxstyle='round,pad=0.3', fc='lavender', ec='purple', alpha=0.85))
    ax.text(0.985, 0.20, 'young clade only\n(our-specific\ncandidates)', transform=ax.transAxes,
            fontsize=8.5, color='blue', ha='right',
            bbox=dict(boxstyle='round,pad=0.3', fc='aliceblue', ec='blue', alpha=0.85))

    axh.hist(m.z, bins=80, orientation='horizontal', color='steelblue')
    axh.axhline(0, color='grey', lw=0.8)
    axh.set_ylim(ax.get_ylim())
    axh.set_yticklabels([])
    axh.set_xlabel('count', fontsize=9)
    axh.set_title('z marginal', fontsize=9)

    out_png = ROOT / 'fig/arg_z_vs_akbari.png'
    out_pdf = ROOT / 'fig/arg_z_vs_akbari.pdf'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    print('wrote', out_png)
    print('wrote', out_pdf)
    print(f'shared variants plotted: {len(m)}')
    print('\ntop 8 by our z:')
    print(topz[['lab', 'freq', 'age', 'z', 'pct', 'akbari_X', 'neglog10P']].to_string(index=False))
    print('\ntop 8 by Akbari P:')
    print(topP[['lab', 'freq', 'age', 'z', 'pct', 'akbari_X', 'neglog10P']].to_string(index=False))


if __name__ == '__main__':
    main()
