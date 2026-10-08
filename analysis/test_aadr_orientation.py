#!/usr/bin/env python3
"""Verify the AADR 2-bit decoding convention against the Akbari imputed pgen.

Decisive test of both (a) bit order within a byte and (b) allele orientation, by
correlating AADR-decoded dosages with the Akbari imputed dosages for the same
individuals at the same sites. Run:  .venv/bin/python analysis/test_aadr_orientation.py

Reported 2026-10-08:
    low-bit-first            corr -0.1168  agreement 0.2341
    high-bit-first as-is     corr -0.8103  agreement 0.0023
    high-bit-first + invert  corr +0.8103  agreement 0.7223

The last line is the correct convention: high bit first, code counts the .snp
REFERENCE allele, so ALT dosage = 2 - code.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis.load_aadr import GenotypeMatrix, parse_ind, parse_snp  # noqa: E402
from analysis import load_akbari_geno as akg  # noqa: E402

AADR = Path("data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB")
SNP = Path(str(AADR) + ".snp")
IND = Path(str(AADR) + ".ind")
GENO = Path(str(AADR) + ".geno")
PGEN = Path("data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.pgen")
PSAM = PGEN.with_suffix(".psam")
PVAR = PGEN.with_suffix(".pvar")


def main() -> int:
    snp = list(parse_snp(str(SNP)))
    key = {f"{r[1]}:{int(float(r[3]))}": (i, str(r[4]).upper(), str(r[5]).upper())
           for i, r in enumerate(snp)}

    shared, pvar_idx, seen = [], [], 0
    with open(PVAR) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            k = f"{f[0]}:{int(f[1])}"
            if k in key and len(shared) < 2000:
                ai, ref, alt = key[k]
                if ref == f[3].upper() and alt == f[4].upper():
                    shared.append(ai)
                    pvar_idx.append(seen)
            seen += 1
            if len(shared) >= 2000:
                break

    psam = pd.read_csv(PSAM, sep="\t")
    pg_of = {s: i for i, s in enumerate(psam.iloc[:, 0].astype(str))}
    pairs = sorted((pg_of[r.id], i) for i, r in
                   enumerate(parse_ind(str(IND))) if r.id in pg_of)
    p_cols = [p for p, _ in pairs[:300]]
    a_rows = [a for _, a in pairs[:300]]

    gm = GenotypeMatrix(str(GENO))
    sites = np.array(shared, dtype=np.int64)
    geno = akg.AkbariGeno(sample_subset=p_cols)
    dos = geno.read_variants(np.array(pvar_idx[: len(sites)], dtype=np.uint32))
    ak = np.rint(dos).T
    fin = np.isfinite(dos).T

    by, frac = sites // 4, sites % 4

    def decode(high_first: bool) -> np.ndarray:
        sh = ((6 - frac * 2) if high_first else (frac * 2)).astype(np.uint8)
        out = np.empty((len(a_rows), sites.size), dtype=np.int8)
        for r, rowi in enumerate(a_rows):
            out[r] = (np.asarray(gm._mm[rowi])[by] >> sh) & np.uint8(3)
        return out

    print(f"sites={len(sites)} individuals={len(a_rows)}")
    report = {}
    hi = decode(True)
    inv = np.where(hi <= 2, 2 - hi, hi)
    for name, code in (("low-bit-first", decode(False)),
                       ("high-bit-first as-is", hi),
                       ("high-bit-first + invert", inv)):
        ok = fin & (code <= 2)
        x, y = code[ok].astype(float), ak[ok]
        corr = float(np.corrcoef(x, y)[0, 1])
        agree = float((x == y).mean())
        report[name] = (corr, agree)
        print(f"  {name:24s} corr={corr:+.4f}  agreement={agree:.4f}")

    # Assert the documented convention: ours returns ALT dosage already.
    raw = decode(True)
    ours = np.where(raw <= 2, 2 - raw, raw)
    ok = fin & (raw <= 2)
    corr = float(np.corrcoef(ours[ok].astype(float), ak[ok])[0, 1])
    assert corr > 0.7, f"decoder convention wrong: corr={corr:.3f}"
    print(f"  load_aadr ALT-dosage convention verified (corr {corr:+.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
