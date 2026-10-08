"""Reader for the AADR v66.p1 (Reich lab) plink-style EIGENSTRAT/TGENO bundles.

The genotype files in ``data/raw/aadr-v66.p1/`` use a 48-byte ``TGENO`` header
followed by a densely packed, **individual-major** 2-bit-per-site matrix:

    header   = b"TGENO" + n_ind + n_snp + two 32-bit hex checksums, NUL-padded to
               HEADER_SIZE (48) bytes.  Data starts at offset 48.
    data     = n_ind contiguous rows of ``bytes_per_row = ceil(n_snp / 4)`` bytes.
               Row ``i`` holds the genotypes of individual ``i`` for every site in
               .snp order, 4 sites packed per byte, high 2 bits = first site.
               Genotype codes: 0, 1, 2 = diploid genotype state in .snp
               reference/alternate orientation; 3 = missing.

This layout is verified against the on-disk file sizes (see the self-test): the
48-byte header plus ``n_ind * bytes_per_row`` equals the exact byte size.

The reader memory-maps the file, so the 7 GB matrix is never fully loaded.
Decoding is done one individual row (or a bounded site slice) at a time.

Date semantics in the .anno file (see ``DATE_CONVENTION`` below) follow the AADR
documentation: the primary time axis is ``Date mean in BP`` = mean years
**before 1950 CE** (OxCal mu for direct radiocarbon dates, midpoint/average of a
contextual range otherwise).  Rows whose date is missing, or whose method is
``Modern``, are flagged rather than guessed.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass

import numpy as np

HEADER_SIZE = 48
MAGIC = b"TGENO"
MISSING_CODE = 3

# AADR ancestor-file date convention adopted here (see docs / AADR README):
#   "Date mean in BP in years before 1950 CE [OxCal mu for a direct
#    radiocarbon date, and average of range for a contextual date]"
# So BP is relative to 1950 CE.  We keep the raw value; conversions to calBCE
# are left to callers.
DATE_CONVENTION = "BP relative to 1950 CE (AADR 'Date mean in BP')"

BUNDLE = "data/raw/aadr-v66.p1"
MAIN_PREFIX = "v66.p1_1240K.aadr.patch.PUB"
HO_PREFIX = "v66.p1_HO.aadr.patch.PUB"
ANNO_NAME = "v66.p1_1240K.aadr.PUB.anno"
HO_ANNO_NAME = "v66.p1_HO.aadr.PUB.anno"

SNP_COLUMNS = ("rsid", "chrom", "genetic_pos", "phys_pos", "ref", "alt")
IND_COLUMNS = ("id", "sex", "group")


# --------------------------------------------------------------------------- #
# Header / file sizing
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TgenoHeader:
    n_ind: int
    n_snp: int
    checksum_a: str
    checksum_b: str


def parse_header(geno_path: str) -> TgenoHeader:
    """Parse the 48-byte TGENO header of a .geno file."""
    with open(geno_path, "rb") as fh:
        raw = fh.read(HEADER_SIZE)
    if not raw.startswith(MAGIC):
        raise ValueError(f"{geno_path}: not a TGENO file (magic={raw[:8]!r})")
    fields = raw.rstrip(b"\x00").split()
    # fields = [b'TGENO', n_ind, n_snp, checksum_a, checksum_b]
    if len(fields) != 5:
        raise ValueError(f"{geno_path}: unexpected header {raw!r}")
    return TgenoHeader(
        n_ind=int(fields[1]),
        n_snp=int(fields[2]),
        checksum_a=fields[3].decode(),
        checksum_b=fields[4].decode(),
    )


def bytes_per_row(n_snp: int) -> int:
    """Number of packed bytes per individual row: ceil(n_snp / 4)."""
    return (n_snp + 3) // 4


def expected_geno_size(geno_path: str) -> int:
    h = parse_header(geno_path)
    return HEADER_SIZE + h.n_ind * bytes_per_row(h.n_snp)


# --------------------------------------------------------------------------- #
# .snp / .ind / .anno parsing
# --------------------------------------------------------------------------- #
def parse_snp(snp_path: str):
    """Return the .snp rows as a list of 6-tuples in SNP_COLUMNS order.

    For 1.2M sites this materialises ~1.2M small tuples; use :func:`iter_snp`
    if you only need to stream.  Fields: rsid(str), chrom(str), genetic_pos,
    phys_pos (float,float), ref(str), alt(str).
    """
    return list(iter_snp(snp_path))


def iter_snp(snp_path: str):
    with open(snp_path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            yield (
                parts[0],
                parts[1],
                float(parts[2]),
                float(parts[3]),
                parts[4],
                parts[5],
            )


@dataclass(frozen=True)
class Individual:
    id: str
    sex: str
    group: str


def parse_ind(ind_path: str):
    """Return the .ind rows as a list of :class:`Individual`."""
    return list(iter_ind(ind_path))


def iter_ind(ind_path: str):
    with open(ind_path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            yield Individual(id=parts[0], sex=parts[1], group=parts[2])


@dataclass
class AnnoRecord:
    genetic_id: str
    individual_id: str
    group_id: str
    date_method: str
    date_bp: float | None
    date_sd_bp: float | None
    full_date: str
    latitude: float | None
    longitude: float | None
    assessment: str

    @property
    def has_date(self) -> bool:
        return self.date_bp is not None

    @property
    def has_coords(self) -> bool:
        return self.latitude is not None and self.longitude is not None

    @property
    def is_modern(self) -> bool:
        return self.date_method.startswith("Modern")


def _to_float(value: str) -> float | None:
    value = (value or "").strip()
    if value in ("", "..", "n/a", "NA", "nan", "N/A"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_anno(anno_path: str):
    """Parse the .anno table into a dict keyed by Genetic ID (== .ind col 0).

    Column indices refer to the AADR v66 ancestor-file schema (49 columns).
    """
    COL_GENETIC = 0
    COL_INDIV = 2
    COL_GROUP = 14
    COL_METHOD = 9
    COL_DATE_BP = 10
    COL_DATE_SD = 11
    COL_FULL_DATE = 12
    COL_LAT = 17
    COL_LON = 18
    COL_ASSESS = 47

    out: dict[str, AnnoRecord] = {}
    with open(anno_path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # header
        for row in reader:
            if not row:
                continue
            out[row[COL_GENETIC]] = AnnoRecord(
                genetic_id=row[COL_GENETIC],
                individual_id=row[COL_INDIV],
                group_id=row[COL_GROUP],
                date_method=row[COL_METHOD],
                date_bp=_to_float(row[COL_DATE_BP]),
                date_sd_bp=_to_float(row[COL_DATE_SD]),
                full_date=row[COL_FULL_DATE],
                latitude=_to_float(row[COL_LAT]),
                longitude=_to_float(row[COL_LON]),
                assessment=row[COL_ASSESS],
            )
    return out


# --------------------------------------------------------------------------- #
# Genotype matrix (memory-mapped)
# --------------------------------------------------------------------------- #
class GenotypeMatrix:
    """Memory-mapped access to a packed individual-major TGENO matrix.

    Row ``i`` corresponds to the i-th line of the matching .ind file and the
    i-th individual-major row of the packed .geno data.
    """

    def __init__(self, geno_path: str):
        self.path = geno_path
        self.header = parse_header(geno_path)
        self.n_ind = self.header.n_ind
        self.n_snp = self.header.n_snp
        self.row_bytes = bytes_per_row(self.n_snp)
        self._mm = np.memmap(
            geno_path,
            dtype=np.uint8,
            mode="r",
            offset=HEADER_SIZE,
            shape=(self.n_ind, self.row_bytes),
        )

    # -- raw packed bytes ---------------------------------------------------- #
    def packed_row(self, i: int) -> np.ndarray:
        return np.asarray(self._mm[i])

    # -- decode -------------------------------------------------------------- #
    @staticmethod
    def _decode_bytes(packed: np.ndarray, sites: np.ndarray | None = None) -> np.ndarray:
        """Decode 2-bit codes from packed bytes.

        If ``sites`` is given (site indices into this row, sorted or not), only
        those sites are decoded and returned in the requested order; otherwise
        all ``4 * len(packed)`` sites are returned.
        """
        if sites is None:
            # EIGENSTRAT packs high->low: site 0 = bits 6-7.
            shifts = np.array([6, 4, 2, 0], dtype=np.uint8)
            codes = (packed[:, None] >> shifts[None, :]) & np.uint8(0x3)
            return codes.reshape(-1)
        byte_idx = sites // 4
        shift = (6 - (sites % 4) * 2).astype(np.uint8)
        return (packed[byte_idx] >> shift) & np.uint8(0x3)

    def individual_genotypes(self, i: int, sites: np.ndarray | None = None) -> np.ndarray:
        """Return uint8 ALT-allele dosages (0/1/2, 3 missing) for individual ``i``.

        ``sites`` optionally selects/orders a subset of site indices.  When
        ``sites`` is None all sites are decoded.

        ORIENTATION (verified 2026-10-08 against the Akbari imputed pgen):
        the packed 2-bit code in this release counts the **reference** allele of
        ``.snp`` column 5, i.e.  the ALT dosage is ``2 - code``.  This method
        returns ALT dosage so that it matches the convention of
        ``analysis.load_akbari_geno`` (and every other reader in this repo).
        Evidence: for 2,000 matched sites and 300 shared individuals, decoding
        high-bit-first *and* subtracting from 2 gives ALT-dosage correlation
        +0.81 / call agreement 0.72 with the imputed calls, versus -0.12 / 0.23
        for low-bit-first and -0.81 / 0.002 for high-bit-first without the
        orientation correction.  See analysis/test_aadr_orientation.py.

        A uniform per-site inversion is invariant for covariance/relatedness
        (both members' centred values flip together) and flips only the *sign*
        of any slope, so results that used the old convention keep their
        magnitudes, |z|, and concordance; only reported directions flip.
        """
        packed = self.packed_row(i)
        codes = self._decode_bytes(packed, sites)
        if sites is None:
            codes = codes[: self.n_snp]
        # missing (3) stays 3; real dosages are flipped to the ALT orientation
        alt = np.where(codes <= 2, 2 - codes, codes).astype(np.uint8)
        return alt

    # -- iteration ----------------------------------------------------------- #
    def iter_individuals(self, sites: np.ndarray | None = None):
        """Yield ``(individual_index, codes)`` one individual row at a time."""
        for i in range(self.n_ind):
            yield i, self.individual_genotypes(i, sites)

    def iter_chunks(self, ind_start: int = 0, ind_stop: int | None = None,
                    sites: np.ndarray | None = None, ind_chunk: int = 256):
        """Yield ``(indices, codes_matrix)`` in bounded individual blocks.

        ``codes_matrix`` has shape ``(len(indices), n_sites_requested)`` and
        holds uint8 codes.  Memory stays bounded by ``ind_chunk`` rows.
        """
        if ind_stop is None:
            ind_stop = self.n_ind
        for start in range(ind_start, ind_stop, ind_chunk):
            stop = min(start + ind_chunk, ind_stop)
            idx = np.arange(start, stop)
            block = np.empty(
                (stop - start, self.n_snp if sites is None else len(sites)),
                dtype=np.uint8,
            )
            for k, i in enumerate(idx):
                block[k] = self.individual_genotypes(i, sites)
            yield idx, block


# --------------------------------------------------------------------------- #
# Top-level convenience
# --------------------------------------------------------------------------- #
def load_bundle(prefix: str = MAIN_PREFIX, base: str = BUNDLE,
                anno_name: str = ANNO_NAME):
    """Return ``(matrix, ind, anno)`` for a bundle prefix under ``base``.

    ``ind`` is an ordered tuple of :class:`Individual` aligned with the matrix
    rows; ``anno`` is the Genetic-ID keyed dict from :func:`parse_anno`.
    """
    geno_path = os.path.join(base, prefix + ".geno")
    ind_path = os.path.join(base, prefix + ".ind")
    return GenotypeMatrix(geno_path), parse_ind(ind_path), parse_anno(
        os.path.join(base, anno_name)
    )


def _count_lines(path: str) -> int:
    n = 0
    with open(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n


def _self_test():
    print("AADR bundle:", BUNDLE)
    for prefix, anno in ((MAIN_PREFIX, ANNO_NAME), (HO_PREFIX, HO_ANNO_NAME)):
        geno = os.path.join(BUNDLE, prefix + ".geno")
        ind = os.path.join(BUNDLE, prefix + ".ind")
        snp = os.path.join(BUNDLE, prefix + ".snp")
        hdr = parse_header(geno)
        size = os.path.getsize(geno)
        exp = expected_geno_size(geno)
        n_ind_lines = _count_lines(ind)
        n_snp_lines = _count_lines(snp)
        print(f"\n[{prefix}]")
        print(f"  header: n_ind={hdr.n_ind} n_snp={hdr.n_snp} "
              f"checksums={hdr.checksum_a},{hdr.checksum_b}")
        print(f"  .ind lines={n_ind_lines}  .snp lines={n_snp_lines}")
        print(f"  geno size={size} expected={exp} match={size == exp}")
        assert hdr.n_ind == n_ind_lines, "header n_ind != .ind line count"
        assert hdr.n_snp == n_snp_lines, "header n_snp != .snp line count"
        assert size == exp, "geno size != header-derived size"
        # decode a small slice and sanity-check code range
        gm = GenotypeMatrix(geno)
        codes = gm.individual_genotypes(0, np.arange(0, min(1000, gm.n_snp)))
        assert codes.dtype == np.uint8 and codes.max() <= 3
        print(f"  first individual first-1000 codes unique={sorted(set(codes.tolist()))}")
    print("\nself-test OK")


if __name__ == "__main__":
    _self_test()
