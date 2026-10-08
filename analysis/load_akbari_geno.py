#!/usr/bin/env python3
"""Chunked reader for the Akbari et al. 2026 PLINK2 imputed genotype bundle.

Source files
------------
``data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.{pgen,pvar,psam}``

Released on Harvard Dataverse V4.0 (doi:10.7910/DVN/7RVV9N, CC0) as the
"Genome-wide imputed genotypes for ancient and modern Europeans (PLINK2)" set
used for the genome-wide analyses in Akbari et al. (2026), Nature 654:419-428.
It combines 15,836 ancient individuals with 503 modern European 1000 Genomes
individuals, autosomes only, AllVariants_NoQC (no post-imputation QC applied),
52,382,873 variants on GRCh37.

Layout
------
* ``.psam`` : ``#IID\tSEX`` header, one individual per line (FID column absent,
  so the individual id in column 1 is the sample id used by AADR/trajectory files).
* ``.pvar`` : ``#CHROM\tPOS\tID\tREF\tALT`` header, one variant per line. ``ID`` is
  ``CHROM_POS_REF_ALT``, matching the selection-summary-statistics convention.
* ``.pgen`` : PLINK2 binary genotype matrix. Read through ``pgenlib.PgenReader``,
  which touches the file variant-by-variant (optionally range/list); the 26 GB
  matrix is never materialised in RAM.

The reader exposes :func:`read_variants` / :func:`iter_chunks` for dosage slices,
plus :class:`AkbariGeno` as a small stateful wrapper. :func:`describe_pgen_storage`
reports the storage mode actually encountered and whether values are integer
hard-calls or fractional dosages.
"""
from __future__ import annotations

import argparse
import resource
from pathlib import Path

import numpy as np
import pgenlib

REPO = Path(__file__).resolve().parent.parent
DEFAULT_DIR = REPO / "data" / "raw" / "akbari" / "tier2"
PREFIX = "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes"
PVAR_HEADER = ("#CHROM", "POS", "ID", "REF", "ALT")
PSAM_HEADER = ("#IID", "SEX")
MISSING_HARDCALL = -9  # pgenlib's int8 missing sentinel


def paths(dirpath=DEFAULT_DIR, prefix=PREFIX):
    d = Path(dirpath)
    return d / f"{prefix}.pgen", d / f"{prefix}.pvar", d / f"{prefix}.psam"


# --------------------------------------------------------------------------- #
# .psam / .pvar parsing
# --------------------------------------------------------------------------- #
def parse_psam(psam_path) -> np.ndarray:
    """Return the sample ids (column 1) as a numpy array of str."""
    ids = []
    with open(psam_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields:
                ids.append(fields[0])
    return np.array(ids, dtype=object)


def iter_pvar(pvar_path):
    """Yield ``(chrom, pos, vid, ref, alt)`` per variant, streaming (no full load)."""
    with open(pvar_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt = line.rstrip("\n").split("\t")[:5]
            yield chrom, int(pos), vid, ref, alt


def read_pvar_header(pvar_path):
    with open(pvar_path, "r", encoding="utf-8") as fh:
        for line in fh:
            return line.rstrip("\n").split("\t")
    return list(PVAR_HEADER)


# --------------------------------------------------------------------------- #
# pgen description
# --------------------------------------------------------------------------- #
def describe_pgen_storage(pgen_path, probe_variant=0, sample_limit=32):
    """Inspect the .pgen storage mode and dosage representation.

    Returns a dict describing: variant/sample counts, whether sparse dosages /
    phased hardcalls are present, and whether probed values look like integer
    hard-calls (0/1/2/-9) or fractional dosages.
    """
    reader = pgenlib.PgenReader(str(pgen_path).encode())
    info = {
        "n_variants": reader.get_variant_ct(),
        "n_samples": reader.get_raw_sample_ct(),
    }
    try:
        info["has_sparse_dosages"] = bool(reader.has_sparse_dosages(probe_variant))
    except Exception:
        info["has_sparse_dosages"] = None

    n = min(sample_limit, info["n_samples"])
    hard = np.empty(info["n_samples"], dtype=np.int8)
    reader.read(probe_variant, hard, allele_idx=1)
    info["hardcall_probe"] = hard[:n].tolist()

    dos = np.empty(info["n_samples"], dtype=np.float32)
    reader.read_dosages(probe_variant, dos, allele_idx=1)
    vals = dos[:n]
    finite = vals[np.isfinite(vals)]
    info["dosage_probe"] = [round(float(v), 4) for v in vals]
    info["dosage_is_integer_like"] = bool(
        finite.size and np.all(np.abs(finite - np.round(finite)) < 1e-3)
    )
    # PLINK2 pgen modes are encoded in the first two bytes of the file; a quick
    # read lets us report the mode label pgenlib's library exposes indirectly.
    with open(pgen_path, "rb") as fh:
        mode_byte = fh.read(1)[0]
    info["pgen_mode_byte"] = mode_byte
    info["storage_mode"] = _mode_label(mode_byte, info)
    reader.close()
    return info


_MODE_NAMES = {
    0x10: "PLINK1 .bed (uncompressed)",
    0x11: "vrec (variable-record, 2-bit hardcalls)",
    0x12: "vrec + sparse hardcalls",
    0x13: "vrec + dosage list",
    0x14: "vrec + sparse dosage list",
    0x15: "vrec + dosage list + sparse hardcalls",
    0x16: "vrec + dosage list + sparse hardcalls + sparse dosages",
    0x17: "vrec + dosage list + sparse hardcalls + sparse dosages + phaseinfo",
}


def _mode_label(mode_byte, info):
    label = _MODE_NAMES.get(mode_byte, f"unknown (mode byte {mode_byte:#04x})")
    if info.get("has_sparse_dosages"):
        label += " [sparse dosages present]"
    if info.get("dosage_is_integer_like"):
        label += " [observed values integer-like -> hard-call equivalent]"
    else:
        label += " [observed fractional values -> true dosages]"
    return label


# --------------------------------------------------------------------------- #
# Stateful wrapper
# --------------------------------------------------------------------------- #
class AkbariGeno:
    """Thin stateful reader over the pgen/psam/pvar bundle (low memory).

    Parameters
    ----------
    dirpath, prefix : where the three files live.
    sample_subset : optional 0-based sample indices. If given, ``n_samples``
        and all reads refer to that subset only (via the native
        ``change_sample_subset``), so only those columns are touched.
    """

    def __init__(self, dirpath=DEFAULT_DIR, prefix=PREFIX, sample_subset=None):
        pgen_path, pvar_path, psam_path = paths(dirpath, prefix)
        self.pgen_path = str(pgen_path)
        self.pvar_path = str(pvar_path)
        self.psam_path = str(psam_path)

        self.sample_ids = parse_psam(psam_path)
        self._reader = pgenlib.PgenReader(str(self.pgen_path).encode())
        self.n_variants = self._reader.get_variant_ct()
        self._raw_n_samples = self._reader.get_raw_sample_ct()

        self.sample_subset = None
        self.change_sample_subset(sample_subset)

    # -- sample selection ---------------------------------------------------
    def change_sample_subset(self, sample_subset=None):
        """Restrict reads to 0-based sample indices; ``None`` selects all."""
        if sample_subset is None:
            self._reader.change_sample_subset(None)
            self.sample_subset = None
            self.n_samples = self._raw_n_samples
        else:
            idx = np.asarray(sample_subset, dtype=np.int32)
            self._reader.change_sample_subset(idx)
            self.sample_subset = idx
            self.n_samples = idx.size

    def sample_index_of(self, sample_id):
        """0-based index (in the full sample list) of a sample id, or -1."""
        hits = np.where(self.sample_ids == sample_id)[0]
        return int(hits[0]) if hits.size else -1

    # -- dosage reads -------------------------------------------------------
    def read_variants(self, indices, allele_idx=1):
        """Return dosages for the given variant indices.

        Shape ``(len(indices), n_samples)`` (variant-major), dtype float32.
        Dosages are ALT-allele counts in ``[0, 2]``; missing is NaN. Only the
        requested variants are read from disk.
        """
        idx = np.asarray(indices, dtype=np.int32)
        out = np.empty((idx.size, self.n_samples), dtype=np.float32)
        if idx.size == 1:
            self._reader.read_dosages(int(idx[0]), out[0], allele_idx=allele_idx)
        else:
            self._reader.read_dosages_list(
                idx.tolist(), out, allele_idx=allele_idx, sample_maj=0
            )
        return out

    def read_hardcalls(self, indices, allele_idx=1):
        """Return integer hard-calls (0/1/2, missing=-9) for the given indices."""
        idx = np.asarray(indices, dtype=np.int32)
        out = np.empty((idx.size, self.n_samples), dtype=np.int8)
        if idx.size == 1:
            self._reader.read(int(idx[0]), out[0], allele_idx=allele_idx)
        else:
            self._reader.read_list(idx.tolist(), out, allele_idx=allele_idx, sample_maj=0)
        return out

    def iter_chunks(self, start=0, end=None, step=1_000_000, allele_idx=1):
        """Yield ``(variant_offsets, dosage_matrix)`` for contiguous chunks.

        Reads at most ``step`` variants at a time, so peak memory is bounded by
        ``step * n_samples`` float32 (1M variants x 16,339 samples ~ 65 GB for the
        full sample set -- choose a smaller ``step`` or a sample subset). Use the
        default ``step`` only with a sample subset, or lower it.
        """
        end = self.n_variants if end is None else min(end, self.n_variants)
        for s in range(start, end, step):
            e = min(s + step, end)
            yield np.arange(s, e), self.read_variants(np.arange(s, e), allele_idx=allele_idx)

    def close(self):
        self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #
def self_test(dirpath=DEFAULT_DIR, prefix=PREFIX, n_samples_probe=8):
    info = describe_pgen_storage(paths(dirpath, prefix)[0])
    print("=== pgen storage ===")
    print(f"n_variants        : {info['n_variants']:,}")
    print(f"n_samples         : {info['n_samples']:,}")
    print(f"storage_mode      : {info['storage_mode']}")
    print(f"dosage_is_integer_like: {info['dosage_is_integer_like']}")

    geno = AkbariGeno(dirpath, prefix)
    assert geno.n_samples == info["n_samples"], "psam/pgen sample count mismatch"
    print("\n=== first 3 pvar rows ===")
    for i, row in enumerate(iter_pvar(geno.pvar_path)):
        print(f"  idx={i} {row}")
        if i >= 2:
            break
    print("\n=== dosage at variant 0 for first samples ===")
    dos = geno.read_variants([0])[0]
    hard = geno.read_hardcalls([0])[0]
    for i in range(min(n_samples_probe, geno.n_samples)):
        print(f"  {geno.sample_ids[i]:<14} hardcall={hard[i]:>3}  dosage={dos[i]:.4f}")
    frac = dos[np.isfinite(dos)]
    print(f"\nfractional dosages among first variant: "
          f"{int(np.sum(np.abs(frac - np.round(frac)) > 1e-3))}/{frac.size}")
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"peak RSS: {peak_kb/1024:.1f} MiB")
    geno.close()
    return info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", default=str(DEFAULT_DIR))
    ap.add_argument("--prefix", default=PREFIX)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        self_test(args.dir, args.prefix)
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
