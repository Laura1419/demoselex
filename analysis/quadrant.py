#!/usr/bin/env python3
"""Four-quadrant comparison: Akbari |X| vs our trajectory-replication statistic.

Quadrant definitions (stated explicitly so the figure is not over-read):
  Akbari "high"   : |X| >= 3.16  (the paper's FDR-50% threshold; 5.45 = pi>99%)
  Akbari "low"    : |X| <  3.16
  Replication high: matched-control empirical P <= 0.05 (BH-FDR reported too)
  Replication low : otherwise

Q1 high/high  consensus
Q2 high/low   likely one broad demographic event (Akbari-owned)
Q3 low/high   **candidate over-correction** - the class this project targets
Q4 low/low    weak / null

Nothing here is a discovery claim: Tier-1 data are ~26% of Akbari's individuals and
1.23M of 9.7M sites, and the empirical null is a finite-control caliper matcher.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/erikg/demoselex")
SRC = ROOT / "data/derived/replication_stats.tsv"
OUT_DOC = ROOT / "docs/I-quadrant-first-pass.md"
OUT_FIG = ROOT / "fig/quadrant_akbari_vs_replication.png"
OUT_JSON = ROOT / "data/derived/quadrant_metrics.json"

ABS_X = 3.16
GW_X = 5.45
P_CUT = 0.05


def main() -> int:
    d = pd.read_csv(SRC, sep="\t")
    d = d[d["replication_score"].notna()].copy()
    d["absX"] = d["akbari_X"].abs()
    d["akbari_high"] = d["absX"] >= ABS_X
    d["rep_high"] = d["p_empirical"] <= P_CUT

    def quad(r):
        if r["akbari_high"] and r["rep_high"]:
            return "Q1"
        if r["akbari_high"]:
            return "Q2"
        if r["rep_high"]:
            return "Q3"
        return "Q4"

    d["quadrant"] = d.apply(quad, axis=1)
    counts = d["quadrant"].value_counts().to_dict()

    q3 = d[d["quadrant"] == "Q3"].sort_values("replication_score", ascending=False)
    q1 = d[d["quadrant"] == "Q1"].sort_values("replication_score", ascending=False)

    # also the strongest replication overall, regardless of Akbari
    top_rep = d.sort_values("replication_score", ascending=False).head(25)

    metrics = {
        "n_scored": int(len(d)),
        "abs_x_threshold": ABS_X,
        "genome_wide_threshold": GW_X,
        "p_cut": P_CUT,
        "n_with_empirical_p": int(d["p_empirical"].notna().sum()),
        "quadrant_counts": {k: int(v) for k, v in counts.items()},
        "n_q3_with_fdr_lt_0_05": int((q3["bh_fdr"] < 0.05).sum()) if len(q3) else 0,
        "n_q3_with_absx_ge_5_45": int((q3["absX"] >= GW_X).sum()) if len(q3) else 0,
        "q3_top": q3.head(30)[["locus", "chrom", "pos", "akbari_X", "z_score",
                               "concordance", "n_trajectories", "n_effective",
                               "p_empirical", "bh_fdr"]].to_dict("records"),
        "q1_top": q1.head(15)[["locus", "chrom", "pos", "akbari_X", "z_score",
                               "concordance", "p_empirical", "bh_fdr"]].to_dict("records"),
    }
    OUT_JSON.write_text(json.dumps(metrics, indent=1))

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.2, 6.2))
        colors = {"Q1": "#1b7837", "Q2": "#762a83", "Q3": "#d73027", "Q4": "#bdbdbd"}
        for q in ("Q4", "Q2", "Q1", "Q3"):
            s = d[d["quadrant"] == q]
            if not len(s):
                continue
            ax.scatter(s["absX"], s["replication_score"], s=8, alpha=0.5,
                       c=colors[q], label=f"{q} (n={len(s)})", linewidths=0)
        ax.axvline(ABS_X, ls="--", c="k", lw=0.8)
        ax.axvline(GW_X, ls=":", c="k", lw=0.8)
        ax.set_xlabel("Akbari |X|   (dashed 3.16 = FDR 50%, dotted 5.45 = pi>99%)")
        ax.set_ylabel("trajectory replication score")
        ax.set_title("Akbari signal vs trajectory replication (Tier-1 first pass)")
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT_FIG, dpi=180)
        print(f"wrote {OUT_FIG}")
    except Exception as exc:  # noqa: BLE001
        print(f"figure skipped: {exc}")

    lines = [
        "# I. Akbari vs trajectory replication - first pass (Tier-1)",
        "",
        "Status: **scouting only.** ~26% of Akbari's individuals, 1.23M of 9.7M sites,",
        "partly-stable trajectories, finite-control empirical null.",
        "",
        f"- loci scored: **{len(d)}**",
        f"- with an empirical P: **{int(d['p_empirical'].notna().sum())}**",
        f"- quadrants: {counts}",
        f"- Q3 (low |X| < {ABS_X}, empirical P <= {P_CUT}): **{len(q3)}**",
        f"- of those, BH-FDR < 0.05: **{metrics['n_q3_with_fdr_lt_0_05']}**",
        "",
        "## Q3 candidates (low Akbari |X|, high replication)",
        "",
        "| locus | chrom:pos | X | z | C | n_traj | N_eff | P | FDR |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, r in q3.head(40).iterrows():
        lines.append(
            f"| {r['locus']} | {r['chrom']}:{int(r['pos'])} | {r['akbari_X']:.2f} | "
            f"{r['z_score']:.2f} | {r['concordance']:.2f} | {int(r['n_trajectories'])} | "
            f"{r['n_effective']:.2f} | {r['p_empirical']:.4f} | {r['bh_fdr']:.3g} |"
        )
    lines += [
        "",
        "## Q1 (both high) - consensus",
        "",
        "| locus | chrom:pos | X | z | C | P | FDR |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, r in q1.head(15).iterrows():
        lines.append(
            f"| {r['locus']} | {r['chrom']}:{int(r['pos'])} | {r['akbari_X']:.2f} | "
            f"{r['z_score']:.2f} | {r['concordance']:.2f} | {r['p_empirical']:.4f} | "
            f"{r['bh_fdr']:.3g} |"
        )
    lines += [
        "",
        "## Interpretation rules (do not skip)",
        "",
        "- A Q3 locus is a *candidate*, not a finding: it must survive the relatedness",
        "  discount (it did - the scan applied `traj_relatedness.tsv`), matched controls,",
        "  and ideally the Tier-2 imputed data before any claim is made.",
        "- If Q3 is dominated by ancestry-informative positions or by a few linkage",
        "  blocks, that is evidence the statistic leaks structure, not that Akbari",
        "  missed selection. Check clustering of Q3 positions by chromosome and LD.",
        "- Gene/functional annotation is NOT yet attached; positions are reported raw.",
        "",
    ]
    OUT_DOC.write_text("\n".join(lines))
    print(f"wrote {OUT_DOC}")
    print(f"quadrants: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
