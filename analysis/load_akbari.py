#!/usr/bin/env python3
"""Chunked loader for the Akbari et al. 2026 genome-wide selection summary statistics.

Source file
-----------
``data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz``
(674,584,607 B, sha256 prefix ``b888d5a4``..., CC0, Akbari et al. 2026,
Nature 654:419-428, doi:10.1038/s41586-026-10358-1).

The file is a comment-prefixed (``##``) TSV with one row per variant. It has
9,739,624 data rows; we never materialise the whole table, we stream it in
chunks. Column meanings are documented in :data:`COLUMN_DOC`; they are taken
verbatim from the ``##COLUMN=<Description=...>`` header lines shipped inside
the file, which is the authoritative source of units.

The reported selection statistic (the paper's time-slope Z / "X" statistic) is
the column **``X``**: ``Selection X statistic (Z-score adjusted for inflation)``.
``S`` is the fitted per-variant selection coefficient (the time slope itself,
in the model's per-generation frequency space) and ``SE`` its standard error;
``X`` is the inflation-adjusted Z of that slope. See ``analysis/akbari_inventory.md``
for the verified interpretation and the residual ambiguity about units of ``S``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DEFAULT_PATH = REPO / "data" / "raw" / "akbari" / "Selection_Summary_Statistics_01OCT2025.tsv.gz"

N_DATA_ROWS = 9_739_624  # documented final variant set; verified by --inspect

# Columns kept in the tidy frame returned by load().
TIDY_COLUMNS = ["ID", "CHROM", "POS", "REF", "ALT", "AF", "S", "SE", "X",
                "P_X", "POSTERIOR", "FDR", "CHI2_BE", "FILTER"]

# Verified meaning of every column, quoted from the in-file ## header.
COLUMN_DOC = {
    "CHROM": ("Chromosome (GRCh37/hg19)", "categorical (1..22)"),
    "POS": ("Position on the chromosome (GRCh37/hg19)", "bp, 1-based"),
    "REF": ("Reference allele (GRCh37/hg19)", "nt string"),
    "ALT": ("Alternative allele", "nt string"),
    "ANC": ("Ancestral allele", "nt string; 'NA' when unknown"),
    "ID": ("Unique identifier for the variant, defined as CHROM_POS_REF_ALT", "string"),
    "RSID": ("Reference SNP ID", "string; 'NA' when absent"),
    "AF": ("Allele frequency of the alternative allele (ALT)", "fraction in [0,1]"),
    "S": ("Selection coefficient of the alternative allele (ALT)", "model slope"),
    "SE": ("Standard error for the selection coefficient", "same units as S"),
    "X": ("Selection X statistic (Z-score adjusted for inflation)", "dimensionless Z"),
    "P_X": ("P-value of X, assuming it follows a normal distribution", "probability"),
    "POSTERIOR": ("Posterior probability of a true selection signal, given the observed X", "probability"),
    "FDR": ("False discovery rate, expected proportion of false positives among detected signals", "probability"),
    "CHI2_BE": ("Chi-square value for batch effect (higher = more suspect)", "chi-square"),
    "FILTER": ("Status/warning/error messages associated with the variant", "string"),
}


def _reader(path, chunksize, usecols=None):
    return pd.read_csv(
        path, sep="\t", comment="#", chunksize=chunksize,
        usecols=usecols, dtype={"CHROM": "Int64", "RSID": "string", "FILTER": "string"},
        na_values=["NA"], keep_default_na=True,
    )


def iter_chunks(path=DEFAULT_PATH, chunksize=1_000_000, usecols=None):
    """Yield the summary-statistics table as successive DataFrames.

    Streaming keeps peak RAM at ~one chunk (default 1M rows, <200 MB), so the
    full 9.7M-row file never has to fit in memory.
    """
    yield from _reader(path, chunksize, usecols=usecols)


def load(path=DEFAULT_PATH, *, chrom=None, loci=None, min_abs_x=None,
         max_p_x=None, filters=None, chunksize=1_000_000):
    """Return a tidy DataFrame of variants matching the given bounds.

    Parameters
    ----------
    chrom : int or iterable of int, optional
        Restrict to one or more chromosomes (1..22).
    loci : iterable of str, optional
        Restrict to variant ``ID`` values (``CHROM_POS_REF_ALT``).
    min_abs_x : float, optional
        Keep rows with ``abs(X) >= min_abs_x`` (the selection-statistic
        threshold; pass 5.45 for the paper's genome-wide-significance cut).
    max_p_x : float, optional
        Keep rows with ``P_X <= max_p_x``.
    filters : iterable of str, optional
        Keep rows whose ``FILTER`` equals one of these values (e.g. ``["PASS"]``).
    """
    chrom_set = None
    if chrom is not None:
        chrom_set = {int(chrom)} if np.isscalar(chrom) else {int(c) for c in chrom}
    loci_set = None if loci is None else set(loci)
    filter_set = None if filters is None else set(filters)

    parts = []
    for chunk in iter_chunks(path, chunksize=chunksize):
        if chrom_set is not None:
            chunk = chunk[chunk["CHROM"].isin(chrom_set)]
        if loci_set is not None:
            chunk = chunk[chunk["ID"].isin(loci_set)]
        if min_abs_x is not None:
            chunk = chunk[chunk["X"].abs() >= min_abs_x]
        if max_p_x is not None:
            chunk = chunk[chunk["P_X"] <= max_p_x]
        if filter_set is not None:
            chunk = chunk[chunk["FILTER"].isin(filter_set)]
        if len(chunk):
            parts.append(chunk)
    if not parts:
        return pd.DataFrame(columns=TIDY_COLUMNS)
    return pd.concat(parts, ignore_index=True)[TIDY_COLUMNS]


def iter_variants(path=DEFAULT_PATH, **bounds):
    """Convenience: iterate matching rows chunk-by-chunk (low memory)."""
    for chunk in iter_chunks(path):
        for name, val in bounds.items():
            if name == "chrom":
                cs = {int(val)} if np.isscalar(val) else {int(c) for c in val}
                chunk = chunk[chunk["CHROM"].isin(cs)]
            elif name == "min_abs_x":
                chunk = chunk[chunk["X"].abs() >= val]
            elif name == "max_p_x":
                chunk = chunk[chunk["P_X"] <= val]
            elif name == "filters":
                chunk = chunk[chunk["FILTER"].isin(set(val))]
            else:
                raise ValueError(f"unknown bound {name!r}")
        if len(chunk):
            yield chunk


def inspect(path=DEFAULT_PATH, chunksize=1_000_000):
    """Print the real header/dtypes and stream the file for a row count."""
    header = []
    import gzip
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                header.append(line.rstrip("\n"))
            else:
                colnames = line.rstrip("\n").split("\t")
                break
    print("=== in-file ## header lines ===")
    for h in header:
        print(h)
    print("\n=== column names (first data line) ===")
    print("\t".join(colnames))
    print("\n=== dtypes (first chunk) ===")
    first = next(iter_chunks(path, chunksize=chunksize))
    print(first.dtypes.to_string())
    print("\n=== streaming row count ===")
    n = 0
    for ch in iter_chunks(path, chunksize=chunksize):
        n += len(ch)
    print(f"data rows: {n}")
    print(f"matches documented final variant set {N_DATA_ROWS}: {n == N_DATA_ROWS}")
    return n


def top(path=DEFAULT_PATH, n=25, by="abs_x", chunksize=1_000_000):
    """Return the top-``n`` rows by ``|X|`` (default) or signed ``X``."""
    key = (lambda d: d["X"].abs()) if by == "abs_x" else (lambda d: d["X"])
    best = None
    for chunk in iter_chunks(path, chunksize=chunksize):
        k = key(chunk)
        take = chunk.loc[k.nlargest(n).index]
        best = take if best is None else pd.concat([best, take], ignore_index=True)
        best = best.assign(_k=key(best)).nlargest(n, "_k").drop(columns="_k")
    return best.reset_index(drop=True)[TIDY_COLUMNS]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--path", type=Path, default=DEFAULT_PATH)
    ap.add_argument("--inspect", action="store_true", help="print real header + row count")
    ap.add_argument("--top", type=int, default=None, metavar="N", help="extract top-N rows by |X|")
    ap.add_argument("--chrom", type=int, default=None)
    ap.add_argument("--min-abs-x", type=float, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    args = ap.parse_args(argv)

    if args.inspect:
        inspect(args.path, chunksize=args.chunksize)
        return 0
    if args.top is not None:
        df = top(args.path, n=args.top, chunksize=args.chunksize)
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(args.out, sep="\t", index=False)
            print(f"wrote {len(df)} rows -> {args.out}")
        else:
            print(df.to_string(index=False))
        return 0
    if args.out is not None:
        df = load(args.path, chrom=args.chrom, min_abs_x=args.min_abs_x,
                  chunksize=args.chunksize)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, sep="\t", index=False)
        print(f"wrote {len(df)} rows -> {args.out}")
        return 0
    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
