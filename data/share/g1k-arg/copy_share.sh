#!/bin/bash
# Stage the allele-analyzer prerequisites onto a shared filesystem for collaborators.
#
#   ./copy_share.sh /moosefs/erikg/demoselex-arg
#
# Copies: small files (analyzer, sample order, annotations, README) + the 294 dated ARG
# chunks (~58 GiB). Safe to re-run: rsync skips what is already there.
set -u
SRC="$(cd "$(dirname "$0")" && pwd)"          # .../data/share/g1k-arg
ARGDIR="${SRC}/../../../data/derived/g1k"      # repo data/derived/g1k
DEST="${1:-/moosefs/erikg/demoselex-arg}"

echo "source small files : ${SRC}"
echo "source ARG chunks  : ${ARGDIR}"
echo "destination        : ${DEST}"

# 0. fail fast if the destination filesystem is not responding (a hung FUSE mount blocks
#    everything and is indistinguishable from a very slow copy).
if ! timeout 20 ls -d "$(dirname "${DEST}")" >/dev/null 2>&1; then
  echo "ERROR: $(dirname "${DEST}") is not responding (mount hung?). Nothing copied." >&2
  exit 1
fi

mkdir -p "${DEST}/arg" || exit 1

# 1. small files
rsync -a --info=progress2 \
  "${SRC}/README.md" \
  "${SRC}/g1k_allele_analyzer.py" \
  "${SRC}/samples_population.txt" \
  "${SRC}/sv265_annotation.bed" \
  "${SRC}/genomicSuperDups_hg38.txt.gz" \
  "${SRC}/loci-analysis-note.md" \
  "${DEST}/" || exit 1

# 2. the ARG chunks (the big part)
rsync -a --info=progress2 --ignore-existing \
  "${ARGDIR}"/chr*_*_*.dated.trees "${DEST}/arg/" || exit 1

# 3. provenance + openness
( cd "${DEST}" && sha256sum README.md g1k_allele_analyzer.py samples_population.txt \
    sv265_annotation.bed genomicSuperDups_hg38.txt.gz > SHA256SUMS.small 2>/dev/null )
( cd "${DEST}/arg" && ls -1 chr*_*_*.dated.trees | wc -l > FILE_COUNT.txt )
chmod -R a+rX "${DEST}" 2>/dev/null
find "${DEST}" -type d -exec chmod a+rx {} \; 2>/dev/null

echo
echo "done. $(ls -1 "${DEST}/arg"/chr*_*_*.dated.trees 2>/dev/null | wc -l) ARG chunks present."
echo "tell collaborators:  ${DEST}"
