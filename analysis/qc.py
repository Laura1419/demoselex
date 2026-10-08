"""Quality-control pass over the AADR v66.p1 1240K bundle.

Produces ``data/derived/qc_manifest.json`` and the numbers quoted in
``analysis/qc_report.md``.

What it computes
----------------
* Counts: N individuals, N sites (from the TGENO header, cross-checked against
  the .ind / .snp line counts).
* Date coverage: per-individual date (BP before 1950 CE) and lat/lon from the
  .anno table, with method categories and rows flagged as uncertain.
* Missingness + MAF on a bounded random subset of sites (default 200k, fixed
  seed) computed in a single streaming pass over the memory-mapped matrix.
* A crude pairwise relatedness screen (IBS concordance and opposite-homozygote
  rate) on the HO subset, over bounded random subsets of individuals and sites.

Everything is reproducible from a fixed seed; versions are recorded in the
manifest.
"""

from __future__ import annotations

import json
import os
import platform
import sys
from collections import Counter

import numpy as np

# Allow both `python3 analysis/qc.py` and `python3 -m analysis.qc`.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.load_aadr import (
    ANNO_NAME,
    BUNDLE,
    HO_ANNO_NAME,
    HO_PREFIX,
    MAIN_PREFIX,
    GenotypeMatrix,
    parse_anno,
    parse_ind,
)

SEED = 20261008
N_SITES_SAMPLE = 200_000          # main-set missingness / MAF subset
MIN_MAF_KEEP = 0.05               # downstream trajectory site filter (see report)
MAX_IND_MISSING_KEEP = 0.90       # individual-level missingness filter
REL_IBS_THRESHOLD = 0.90          # IBS concordance -> "related candidate"
DUP_IBS_THRESHOLD = 0.97          # IBS concordance -> "duplicate candidate"

HOREL_N_IND = 1200                # relatedness screen: individuals sampled
HOREL_N_SITE = 15_000             # relatedness screen: sites sampled

OUT_MANIFEST = "data/derived/qc_manifest.json"


def _quantiles(x: np.ndarray) -> dict:
    if x.size == 0:
        return {}
    qs = [0, 1, 5, 25, 50, 75, 95, 99, 100]
    vals = np.percentile(x, qs)
    return {f"p{q}": float(v) for q, v in zip(qs, vals)}


def date_summary(anno_path: str) -> dict:
    anno = parse_anno(anno_path)
    recs = list(anno.values())
    n = len(recs)
    has_date = sum(r.has_date for r in recs)
    has_coords = sum(r.has_coords for r in recs)
    usable = [r for r in recs if r.has_date and r.has_coords and not r.is_modern]
    dates = np.array([r.date_bp for r in usable], dtype=float)
    method_cat = Counter()
    for r in recs:
        m = r.date_method.strip()
        if m.startswith("Direct"):
            cat = "Direct"
        elif m.startswith("Context"):
            cat = "Context"
        elif m.startswith("Modern"):
            cat = "Modern"
        else:
            cat = "Other/Unlabeled"
        method_cat[cat] += 1
    # uncertainty flags
    no_date = sum(not r.has_date for r in recs)
    no_coords = sum(not r.has_coords for r in recs)
    modern = sum(r.is_modern for r in recs)
    warned = sum(("WARNING" in r.date_method or "CAUTION" in r.date_method)
                 for r in recs)
    return {
        "anno_path": anno_path,
        "n_rows": n,
        "n_has_date_bp": int(has_date),
        "n_has_coords": int(has_coords),
        "n_usable_date_and_coords_not_modern": len(usable),
        "n_missing_date": int(no_date),
        "n_missing_coords": int(no_coords),
        "n_modern_date": int(modern),
        "n_date_method_warned": int(warned),
        "date_bp_min": float(dates.min()) if dates.size else None,
        "date_bp_max": float(dates.max()) if dates.size else None,
        "date_bp_median": float(np.median(dates)) if dates.size else None,
        "date_bp_quantiles": _quantiles(dates),
        "date_method_categories": dict(method_cat.most_common()),
        "date_convention": "BP before 1950 CE (AADR 'Date mean in BP'); "
                           "Modern/Context flagged, not converted",
    }


def missingness_maf(matrix: GenotypeMatrix, anno: dict, seed: int,
                    n_sites: int) -> dict:
    rng = np.random.default_rng(seed)
    n_snp = matrix.n_snp
    take = min(n_sites, n_snp)
    sites = np.sort(rng.choice(n_snp, size=take, replace=False))

    n_ind = matrix.n_ind
    ind_missing = np.zeros(n_ind, dtype=np.float32)
    site_missing = np.zeros(take, dtype=np.int64)
    site_called = np.zeros(take, dtype=np.int64)
    site_alt = np.zeros(take, dtype=np.int64)

    for i in range(n_ind):
        codes = matrix.individual_genotypes(i, sites)
        miss = codes == 3
        ind_missing[i] = miss.mean()
        called = ~miss
        site_missing += miss
        site_called += called
        site_alt += (codes == 1).astype(np.int64) * called + \
                    (codes == 2).astype(np.int64) * 2 * called

    site_called_safe = np.where(site_called == 0, 1, site_called)
    site_miss_frac = site_missing / n_ind
    p = site_alt / (2 * site_called_safe)
    maf = np.minimum(p, 1 - p)
    maf[site_called == 0] = np.nan

    # cross-reference with anno usability
    ind = parse_ind(os.path.join(BUNDLE, MAIN_PREFIX + ".ind"))
    has_date = np.array([bool(anno.get(x.id) and anno[x.id].has_date) for x in ind])
    has_coords = np.array([bool(anno.get(x.id) and anno[x.id].has_coords) for x in ind])
    is_modern = np.array([bool(anno.get(x.id) and anno[x.id].is_modern) for x in ind])
    usable = has_date & has_coords & ~is_modern

    return {
        "seed": seed,
        "n_sites_sampled": int(take),
        "ind_missing_quantiles": _quantiles(ind_missing),
        "ind_missing_mean": float(ind_missing.mean()),
        "n_ind_missing_gt_thresh": int(
            (ind_missing > MAX_IND_MISSING_KEEP).sum()),
        "ind_missing_gt_thresh_value": MAX_IND_MISSING_KEEP,
        "site_missing_quantiles": _quantiles(site_miss_frac),
        "site_maf_quantiles": _quantiles(maf[np.isfinite(maf)]),
        "site_n_monomorphic": int(np.sum(np.isfinite(maf) & (maf == 0))),
        "site_n_maf_ge_keep": int(np.sum(np.isfinite(maf) & (maf >= MIN_MAF_KEEP))),
        "maf_keep_threshold": MIN_MAF_KEEP,
        "n_ind_with_date": int(has_date.sum()),
        "n_ind_with_coords": int(has_coords.sum()),
        "n_ind_usable": int(usable.sum()),
        "ind_missing_on_usable": _quantiles(ind_missing[usable]),
    }


def relatedness_screen(matrix: GenotypeMatrix, seed: int,
                       n_ind: int, n_site: int) -> dict:
    rng = np.random.default_rng(seed + 1)
    n_ind_total = matrix.n_ind
    n_snp = matrix.n_snp
    ind_idx = np.sort(rng.choice(n_ind_total, size=min(n_ind, n_ind_total),
                                 replace=False))
    site_idx = np.sort(rng.choice(n_snp, size=min(n_site, n_snp), replace=False))

    G = np.empty((len(ind_idx), len(site_idx)), dtype=np.uint8)
    for k, i in enumerate(ind_idx):
        G[k] = matrix.individual_genotypes(int(i), site_idx)

    missing = G == 3
    M = (~missing).astype(np.float32)
    both_nonmiss = M @ M.T
    denom = np.where(both_nonmiss == 0, np.nan, both_nonmiss)
    match = np.zeros_like(both_nonmiss)
    ibs0 = np.zeros_like(both_nonmiss)
    for c in (0, 1, 2):
        O = (G == c).astype(np.float32)
        match += O @ O.T
    O0 = (G == 0).astype(np.float32)
    O2 = (G == 2).astype(np.float32)
    ibs0 = (O0 @ O2.T) + (O2 @ O0.T)

    ibs = match / denom
    ibs0_rate = ibs0 / denom
    # ignore pairs with too little genotype overlap: rates become noise
    MIN_OVERLAP = 1000
    ibs = np.where(both_nonmiss >= MIN_OVERLAP, ibs, np.nan)
    ibs0_rate = np.where(both_nonmiss >= MIN_OVERLAP, ibs0_rate, np.nan)
    iu = np.triu_indices(len(ind_idx), k=1)
    ibs_vals = ibs[iu]
    ibs0_vals = ibs0_rate[iu]
    ibs_vals = ibs_vals[np.isfinite(ibs_vals)]
    n_overlap_ok = int(np.sum(both_nonmiss[iu] >= MIN_OVERLAP))

    rel = int(np.sum(ibs_vals >= REL_IBS_THRESHOLD))
    dup = int(np.sum(ibs_vals >= DUP_IBS_THRESHOLD))
    return {
        "note": "crude IBS concordance on sampled subsets; not F_ij/kinship",
        "n_ind_sampled": int(len(ind_idx)),
        "n_site_sampled": int(len(site_idx)),
        "min_overlap_sites": MIN_OVERLAP,
        "n_pairs_tested": int(n_overlap_ok),
        "ibs_quantiles": _quantiles(ibs_vals),
        "ibs0_quantiles": _quantiles(ibs0_vals[np.isfinite(ibs0_vals)]),
        "ibs_related_threshold": REL_IBS_THRESHOLD,
        "n_pairs_ibs_ge_related": rel,
        "ibs_dup_threshold": DUP_IBS_THRESHOLD,
        "n_pairs_ibs_ge_dup": dup,
        "n_pairs_ibs0_lt_0p005": int(np.sum(ibs0_vals[np.isfinite(ibs0_vals)] < 0.005)),
    }


def main():
    os.makedirs("data/derived", exist_ok=True)

    main_geno = os.path.join(BUNDLE, MAIN_PREFIX + ".geno")
    ho_geno = os.path.join(BUNDLE, HO_PREFIX + ".geno")
    main_mat = GenotypeMatrix(main_geno)
    ho_mat = GenotypeMatrix(ho_geno)

    anno = parse_anno(os.path.join(BUNDLE, ANNO_NAME))
    ho_anno = parse_anno(os.path.join(BUNDLE, HO_ANNO_NAME))

    print("date summary (main .anno) ...", flush=True)
    dates = date_summary(os.path.join(BUNDLE, ANNO_NAME))

    print("missingness / MAF on main set ...", flush=True)
    maf = missingness_maf(main_mat, anno, SEED, N_SITES_SAMPLE)

    print("relatedness screen on HO subset ...", flush=True)
    rel = relatedness_screen(ho_mat, SEED, HOREL_N_IND, HOREL_N_SITE)

    manifest = {
        "generated_by": "analysis/qc.py",
        "seed": SEED,
        "versions": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "main_bundle": {
            "geno": main_geno,
            "n_ind": main_mat.n_ind,
            "n_snp": main_mat.n_snp,
            "geno_bytes": os.path.getsize(main_geno),
        },
        "ho_bundle": {
            "geno": ho_geno,
            "n_ind": ho_mat.n_ind,
            "n_snp": ho_mat.n_snp,
            "geno_bytes": os.path.getsize(ho_geno),
            "n_anno_rows": len(ho_anno),
        },
        "date_summary_main": dates,
        "missingness_maf_main": maf,
        "relatedness_ho": rel,
    }

    # relatedness pairs flagged (ids) for exclusion policy -- report counts only
    with open(OUT_MANIFEST, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    print("wrote", OUT_MANIFEST, flush=True)


if __name__ == "__main__":
    main()
