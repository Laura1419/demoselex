#!/usr/bin/env python3
"""AMY (salivary amylase) case-study loader and trend tools.

This module handles the AMY structural-variation (copy-number) tables that
motivate the trajectory-replication project. AMY1 dosage is *copy number*
(structural variation), not a SNP genotype, so every consumer here works on
integer copy counts rather than allele dosages.

Inputs (all observed, verified downloads under data/raw/amy/):

  copy_number_table_for_supplement.tsv
      Columns: type, source, sample, p1, p2, AMY1, AMY2A, AMY2B.
      One row per individual. `type` in {Modern, Ancient, Archaic, GTEx, Trio}.
      AMY1 / AMY2A / AMY2B are integer copy numbers (observed calls from the
      Sudmant-lab graph-genotyping pipeline, not inferred by this module).

  amylase_diversity_project-v1.1.zip
      Versioned analysis archive. This module reads two members directly with
      `zipfile` (no extraction needed):
        .../supplementary_tables/ancient_metadata_for_supplement.tsv
            Columns: sample, p2_l0, p2_l1, site, country, region, age,
            latitude, longitude, data_source. `age` is years BP.
        .../supplementary_tables/haplotype_structure_for_supplement.tsv
            Columns: type, source, sample, p1, p2, best_matching_score, H1, H2.
            H1/H2 are the per-haplotype structural labels (observed calls).

Ancient-vs-Akbari linkage uses:
  data/raw/akbari/ENA_PRJEB106907_runs.tsv  (ENA runs; sample identities)
  data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.{ind,anno}  (dates/coords)

CLI
  python3 analysis/amy.py --describe     # N individuals + copy-number dist
  python3 analysis/amy.py --overlap      # ancient <-> AADR/Akbari overlap count
  python3 analysis/amy.py --trend PATH   # AMY1 trend per trajectory label
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
import zipfile
from collections import Counter, defaultdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AMY_DIR = os.path.join(REPO_ROOT, "data", "raw", "amy")
CN_TABLE = os.path.join(AMY_DIR, "copy_number_table_for_supplement.tsv")
AMY_ZIP = os.path.join(AMY_DIR, "amylase_diversity_project-v1.1.zip")
AKBARI_RUNS = os.path.join(REPO_ROOT, "data", "raw", "akbari", "ENA_PRJEB106907_runs.tsv")
AADR_DIR = os.path.join(REPO_ROOT, "data", "raw", "aadr-v66.p1")
AADR_ANNO = os.path.join(AADR_DIR, "v66.p1_1240K.aadr.PUB.anno")

_ZIP_ROOT = "sudmantlab-amylase_diversity_project-9332fc7/graph_genotyping/" \
            "assess_graph_genotypes/supplementary_tables/"
ANCIENT_META_MEMBER = _ZIP_ROOT + "ancient_metadata_for_supplement.tsv"
HAPLO_STRUCT_MEMBER = _ZIP_ROOT + "haplotype_structure_for_supplement.tsv"

CN_COLUMNS = ["type", "source", "sample", "p1", "p2", "AMY1", "AMY2A", "AMY2B"]
ANCIENT_TYPES = ("Ancient", "Archaic")


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_copy_number_table(path: str = CN_TABLE) -> list[dict]:
    """Load the per-individual AMY copy-number table (one row per individual)."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    missing = [c for c in CN_COLUMNS if c not in (rows[0] if rows else {})]
    if missing:
        raise ValueError(f"copy-number table missing columns: {missing}")
    return rows


def _read_zip_member(member: str, path: str = AMY_ZIP) -> list[dict]:
    with zipfile.ZipFile(path) as zf:
        raw = zf.read(member).decode("utf-8")
    return list(csv.DictReader(io.StringIO(raw), delimiter="\t"))


def load_ancient_metadata(path: str = AMY_ZIP) -> list[dict]:
    """Ancient metadata (dates/coordinates). Read from the archive in place."""
    return _read_zip_member(ANCIENT_META_MEMBER, path)


def load_haplotype_structure(path: str = AMY_ZIP) -> list[dict]:
    """Per-haplotype structural labels (H1/H2). Read from the archive in place."""
    return _read_zip_member(HAPLO_STRUCT_MEMBER, path)


def load_trajectory_table(path: str) -> list[dict]:
    """Load a trajectory table (columns: individual, date, trajectory_id).

    Fails loudly when the requested table does not exist: this module never
    invents trajectory labels for individuals that lack them.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"trajectory table not found: {path}\n"
            "AMY trends require labels produced by analysis/trajectory_infer.py. "
            "That file does not exist in this repository yet, so no trajectory "
            "labels are available. This module will not fabricate them. Run "
            "analysis/trajectory_infer.py first (see docs/D-method.md)."
        )
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    need = {"individual", "trajectory_id"}
    if not rows or not need.issubset(rows[0]):
        raise ValueError(
            f"trajectory table {path} must have columns {sorted(need)}; "
            f"got {list(rows[0]) if rows else '[]'}"
        )
    # Downstream canonical column is `date`; accept the `date_BP` spelling that
    # trajectory_v2.py emits (years before 1950 CE) as the same quantity.
    date_col = "date" if "date" in rows[0] else ("date_BP" if "date_BP" in rows[0] else None)
    if date_col is None:
        raise ValueError(
            f"trajectory table {path} must have a `date` or `date_BP` column; "
            f"got {list(rows[0])}"
        )
    if date_col != "date":
        for r in rows:
            r["date"] = r[date_col]
    return rows


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #
def describe(cn_rows: list[dict]) -> dict:
    """N individuals and AMY1/AMY2 copy-number distributions, by type."""
    by_type = defaultdict(list)
    for r in cn_rows:
        by_type[r["type"]].append(r)

    def dist(vals):
        c = Counter(vals)
        return {k: c[k] for k in sorted(c)}

    out = {"n_individuals": len(cn_rows), "types": {}}
    for t, rs in sorted(by_type.items()):
        amy1 = [int(r["AMY1"]) for r in rs]
        out["types"][t] = {
            "n": len(rs),
            "amy1_mean": round(sum(amy1) / len(amy1), 3),
            "amy1_dist": dist(amy1),
            "amy2a_dist": dist(int(r["AMY2A"]) for r in rs),
            "amy2b_dist": dist(int(r["AMY2B"]) for r in rs),
        }
    return out


def _aadr_rows(anno_path: str = AADR_ANNO) -> list[list[str]]:
    with open(anno_path, errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # header
        return [row for row in reader if len(row) >= 17]


def _aadr_token_index(rows: list[list[str]]) -> dict:
    """Map every genetic-id / individual-id / skeletal-code token -> row index."""
    index = {}
    for i, row in enumerate(rows):
        tokens = {row[0].split(".")[0]}
        for col in (2, 3):  # Individual ID, Skeletal code
            v = row[col].strip()
            if v and v not in (".", ".."):
                tokens.add(v)
        for tok in tokens:
            index.setdefault(tok, i)
    return index


def _akbari_sample_ids(runs_path: str = AKBARI_RUNS) -> set:
    ids = set()
    with open(runs_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for part in (row.get("submitted_ftp") or "").split(";"):
                m = re.search(r"([^/]+?)\.(AG|SG|DG|TW|BY|AA|EC|WGC|HO|REF)\.bam", part)
                if m:
                    ids.add(m.group(1))
    return ids


def overlap_report(cn_rows: list[dict]) -> dict:
    """Quantify ancient-AMY <-> AADR/Akbari linkability.

    Returns counts only; the identity of the channel (AADR genetic id vs
    Akbari run) is reported explicitly so a zero is unambiguous.
    """
    ancient_ids = sorted({r["sample"] for r in cn_rows if r["type"] in ANCIENT_TYPES})

    aadr_rows = _aadr_rows()
    tok2row = _aadr_token_index(aadr_rows)

    aadr_matched = {s: tok2row[s] for s in ancient_ids if s in tok2row}

    def has_date_coord(i):
        row = aadr_rows[i]
        return bool(row[9] and row[11] and row[15] and row[16])

    dated = [s for s, i in aadr_matched.items() if has_date_coord(i)]

    akbari = _akbari_sample_ids()
    akbari_aadr_rows = set()
    for x in akbari:
        base = x.split("_")[0]
        for key in (x, base):
            if key in tok2row:
                akbari_aadr_rows.add(tok2row[key])

    both = set(aadr_matched.values()) & akbari_aadr_rows

    return {
        "n_ancient_amy": len(ancient_ids),
        "n_aadr_matched": len(aadr_matched),
        "n_aadr_matched_with_date_coord": len(dated),
        "n_akbari_ena_runs": len(akbari),
        "n_akbari_linked_to_aadr": len(akbari_aadr_rows),
        "n_amy_ancient_in_akbari": len(both),
    }


# --------------------------------------------------------------------------- #
# Trend
# --------------------------------------------------------------------------- #
def _linear_slope(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    denom = sum((x - mx) ** 2 for x in xs)
    if denom == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom


def amy1_trend(traj_rows: list[dict], cn_rows: list[dict]) -> dict:
    """AMY1 copy-number trend over time for each trajectory label.

    `traj_rows` must come from load_trajectory_table (individual, date,
    trajectory_id). Individuals absent from the copy-number table are counted
    as unmatched and excluded from the slope, never imputed.
    """
    cn_by_sample = {r["sample"]: r for r in cn_rows}
    per_traj = defaultdict(list)
    unmatched = 0
    for tr in traj_rows:
        cn = cn_by_sample.get(tr["individual"])
        if cn is None:
            unmatched += 1
            continue
        per_traj[tr["trajectory_id"]].append((float(tr["date"]), int(cn["AMY1"])))

    result = {"unmatched_individuals": unmatched, "trajectories": {}}
    for tid, pts in sorted(per_traj.items()):
        pts.sort()
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        result["trajectories"][tid] = {
            "n": len(pts),
            "amy1_mean": round(sum(ys) / len(ys), 3),
            "slope_per_unit_time": _linear_slope(xs, ys),
            "first": pts[0],
            "last": pts[-1],
        }
    return result


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def make_figures(cn_rows: list[dict], report: dict, out_dir: str) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    made = []
    by_type = defaultdict(list)
    for r in cn_rows:
        by_type[r["type"]].append(int(r["AMY1"]))

    fig, ax = plt.subplots(figsize=(7, 4))
    types = sorted(by_type, key=lambda t: -len(by_type[t]))
    xs = range(0, 21)
    for t in types:
        c = Counter(by_type[t])
        ax.plot(list(xs), [c.get(k, 0) / len(by_type[t]) for k in xs], marker="o",
                ms=3, label=f"{t} (n={len(by_type[t])})")
    ax.set_xlabel("AMY1 copy number")
    ax.set_ylabel("fraction of individuals")
    ax.set_title("AMY1 copy-number distribution by cohort type")
    ax.legend()
    fig.tight_layout()
    p1 = os.path.join(out_dir, "amy1_copy_number_distribution.png")
    fig.savefig(p1, dpi=110)
    plt.close(fig)
    made.append(p1)

    fig, ax = plt.subplots(figsize=(6, 4))
    labels = ["ancient AMY\nindividuals", "linked to AADR\n(date+coords)",
              "linked to\nAkbari runs"]
    vals = [report["n_ancient_amy"], report["n_aadr_matched_with_date_coord"],
            report["n_amy_ancient_in_akbari"]]
    bars = ax.bar(labels, vals, color=["#4477aa", "#66ccee", "#cc6677"])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, str(v), ha="center", va="bottom")
    ax.set_ylabel("individuals")
    ax.set_title("Ancient AMY linkability to trajectory-label sources")
    fig.tight_layout()
    p2 = os.path.join(out_dir, "amy_overlap.png")
    fig.savefig(p2, dpi=110)
    plt.close(fig)
    made.append(p2)
    return made


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--describe", action="store_true",
                    help="print N individuals and copy-number distributions")
    ap.add_argument("--overlap", action="store_true",
                    help="print ancient AMY <-> AADR/Akbari overlap counts")
    ap.add_argument("--trend", metavar="TRAJ_TSV",
                    help="AMY1 trend per trajectory_id from a trajectory table")
    ap.add_argument("--figures", action="store_true",
                    help="write fig/amy_*.png")
    args = ap.parse_args(argv)

    cn_rows = load_copy_number_table()

    if args.describe or not (args.overlap or args.trend or args.figures):
        d = describe(cn_rows)
        print(f"N individuals (all types): {d['n_individuals']}")
        for t, s in d["types"].items():
            print(f"  {t:8s} n={s['n']:5d} AMY1 mean={s['amy1_mean']:.2f} "
                  f"dist={s['amy1_dist']}")

    if args.overlap or args.figures:
        rep = overlap_report(cn_rows)
        if args.overlap:
            print("\nAncient AMY <-> AADR / Akbari linkability")
            for k, v in rep.items():
                print(f"  {k}: {v}")
        if args.figures:
            out = os.path.join(REPO_ROOT, "fig")
            made = make_figures(cn_rows, rep, out)
            print("\nwrote figures:", *made, sep="\n  ")

    if args.trend:
        try:
            traj = load_trajectory_table(args.trend)
        except (FileNotFoundError, ValueError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        res = amy1_trend(traj, cn_rows)
        print(f"\nunmatched individuals: {res['unmatched_individuals']}")
        for tid, s in res["trajectories"].items():
            print(f"  trajectory {tid}: n={s['n']} mean AMY1={s['amy1_mean']} "
                  f"slope={s['slope_per_unit_time']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
