#!/usr/bin/env bash
#
# download_tier2.sh -- acquire the Akbari et al. 2026 *processed imputed* genotype
# set (the PLINK2 "G3" bundle: combined ancient + EUR 1000 Genomes, autosomes,
# AllVariants, NoQC, GRCh37) from Harvard Dataverse V4.0 (CC0).
#
# doi:10.7910/DVN/7RVV9N   ->  https://dataverse.harvard.edu/api/access/datafile/<id>
#
# Why: the reanalysis wants the same individuals and the same ~9.7M variants the
# paper used, not the 1240K AADR subset. This bundle is the pgen/psam/pvar matrix.
#
# Guarantees
#   * idempotent: a file whose size AND md5 already match is left untouched.
#   * resumable: partials live in <dest>.part and resume with `curl -C -`.
#   * verifying: each finished file must match the Dataverse md5 and the
#     expected byte size; sha256 is recorded regardless. A verified file is
#     never overwritten. An unverifiable file is kept as .part, not promoted.
#   * --dry-run: resolves the Dataverse dataset listing, prints every target
#     with its expected size, touches no data file, exits 0.
#   * one log line per target -> data/raw/akbari/tier2_download_log.tsv
#     columns: timestamp  name  http_code  bytes  sha256  dest  url
#
# This script does NOT touch the ~1.3 TB per-chromosome NO-QC BCF (out of budget).
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

DATAVERSE="https://dataverse.harvard.edu/api/access/datafile"
DATASET_API="https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/7RVV9N"
DEST_DIR="data/raw/akbari/tier2"
LOG="data/raw/akbari/tier2_download_log.tsv"

# name | datafile_id | expected_bytes | md5
TARGETS=(
  "README.txt|13596811|10901|7d5fdfe6d1087de7b55f14254ba823eb"
  "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.psam|13415494|221391|787717246189699bec4b12a34884992a"
  "HighQualityVariants.tsv|13415504|315624578|2ec2598e2a6ead74ee11050672cda972"
  "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.pvar|13415518|1700215110|766f95a7311f7dee4a47b86f01038ec8"
  "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.pgen|13415552|26300030448|c9e587230d9bdc8ace0ca1b043294a7c"
)

DRY_RUN=0
case "${1:-}" in
  --dry-run) DRY_RUN=1 ;;
  "" ) ;;
  * ) echo "usage: $0 [--dry-run]" >&2; exit 2 ;;
esac

mkdir -p "$DEST_DIR"

# Resolve the dataset listing once and confirm each datafile id/size against it.
# Writes a resolved catalog (name|id|bytes|md5) to a temp file.
resolve_catalog() {
  local json listing
  json="$(mktemp)"
  if curl -fLsS --max-time 90 -A 'Mozilla/5.0' "$DATASET_API" -o "$json"; then
    listing="$json"
  else
    echo "WARN: could not fetch Dataverse dataset listing; using hard-coded expected sizes" >&2
    listing=""
  fi
  for spec in "${TARGETS[@]}"; do
    IFS='|' read -r name id bytes md5 <<<"$spec"
    if [[ -n "$listing" ]]; then
      local api
      api="$(ID="$id" .venv/bin/python - "$listing" <<'PY' 2>/dev/null || true
import json, os, sys
want = int(os.environ["ID"])
d = json.load(open(sys.argv[1]))
for f in d["data"]["latestVersion"]["files"]:
    df = f["dataFile"]
    if df["id"] == want:
        ck = df.get("md5") or (df.get("checksum") or {}).get("value", "-")
        print(f'{df["filesize"]}|{ck}')
        break
PY
)"
      if [[ -n "$api" ]]; then
        local api_bytes="${api%%|*}" api_md5="${api##*|}"
        [[ "$api_bytes" != "$bytes" ]] && echo "WARN: $name size differs: listing=$api_bytes expected=$bytes" >&2
        [[ "$api_md5" != "$md5" && "$api_md5" != "-" ]] && md5="$api_md5"
        bytes="$api_bytes"
      fi
    fi
    printf '%s|%s|%s|%s\n' "$name" "$id" "$bytes" "$md5"
  done
  rm -f "$json"
}

CATALOG="$(resolve_catalog)"

if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY RUN -- tier2 targets (no network body, no writes):"
  printf '%-64s %-12s %-12s %s\n' "NAME" "DATAFILE" "BYTES" "MD5"
  total=0
  while IFS='|' read -r name id bytes md5; do
    [[ -z "$name" ]] && continue
    printf '%-64s %-12s %-12s %s\n' "$name" "$id" "$bytes" "$md5"
    total=$(( total + bytes ))
  done <<<"$CATALOG"
  echo ""
  echo "Total expected bytes: $total ($(( total / 1000000 )) MB)"
  echo "Destination: $DEST_DIR/"
  exit 0
fi

TS() { date -u +%Y-%m-%dT%H:%M:%SZ; }
[[ -f "$LOG" ]] || printf 'timestamp\tname\thttp_code\tbytes\tsha256\tdest\turl\n' > "$LOG"

sha256_of() { sha256sum "$1" 2>/dev/null | cut -d' ' -f1; }
md5_of()    { md5sum "$1" 2>/dev/null | cut -d' ' -f1; }
size_of()   { wc -c < "$1" 2>/dev/null | tr -d ' '; }

log_line() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "$LOG"; }

overall_rc=0
while IFS='|' read -r name id bytes md5; do
  [[ -z "$name" ]] && continue
  url="$DATAVERSE/$id"
  dest="$DEST_DIR/$name"
  part="$dest.part"

  # 1) Already complete and verified? Refuse to overwrite.
  if [[ -f "$dest" && ! -f "$part" ]]; then
    have_size="$(size_of "$dest")"
    have_md5="$(md5_of "$dest")"
    if [[ "$have_size" == "$bytes" && "$have_md5" == "$md5" ]]; then
      sha="$(sha256_of "$dest")"
      echo "OK    $name (already verified: $have_size B, md5 $have_md5)"
      log_line "$(TS)" "$name" 200 "$have_size" "$sha" "$dest" "$url"
      continue
    fi
    echo "REFUSE $name: exists but fails verification (size=$have_size md5=$have_md5); not overwriting" >&2
    log_line "$(TS)" "$name" 0 "$have_size" "$(sha256_of "$dest")" "$dest" "$url"
    overall_rc=1
    continue
  fi

  # 2) Resumable fetch into .part.
  if [[ -f "$part" ]]; then
    have_part="$(size_of "$part")"
    echo "RESUME $name from $have_part / $bytes B"
  else
    echo "GET   $name ($bytes B)"
  fi
  code="$(curl -fL -C - --retry 5 --retry-delay 5 --retry-all-errors \
            -A 'Mozilla/5.0' "$url" -o "$part" -w '%{http_code}' 2>/dev/null)"
  rc=$?
  [[ -z "$code" ]] && code=0
  got="$(size_of "$part")"

  if [[ $rc -ne 0 ]]; then
    echo "ERR   $name: curl rc=$rc code=$code partial=$got B" >&2
    log_line "$(TS)" "$name" "$code" "$got" "$(sha256_of "$part")" "$part" "$url"
    overall_rc=1
    continue
  fi

  # 3) Verify size + md5, then promote .part -> dest.
  have_md5="$(md5_of "$part")"
  if [[ "$got" == "$bytes" && "$have_md5" == "$md5" ]]; then
    mv -f "$part" "$dest"
    sha="$(sha256_of "$dest")"
    echo "DONE  $name ($got B, md5 $have_md5, sha256 $sha)"
    log_line "$(TS)" "$name" "$code" "$got" "$sha" "$dest" "$url"
  else
    echo "FAIL  $name verification: size=$got/$bytes md5=$have_md5/$md5 (kept as $part)" >&2
    log_line "$(TS)" "$name" "$code" "$got" "$(sha256_of "$part")" "$part" "$url"
    overall_rc=1
  fi
done <<<"$CATALOG"

exit $overall_rc
