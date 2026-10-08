#!/usr/bin/env bash
#
# download_data.sh -- idempotent, resumable, checksum-verifying acquisition of the
# public data catalogued in docs/E-data-inventory.md.
#
# Guarantees:
#   * Never overwrites an existing final file (a verified file is left untouched).
#   * Resumable: partial downloads live in <dest>.part and are resumed with `curl -C -`.
#   * Checksum-verified: a catalogued md5 is enforced; the sha256 of every file is
#     computed and written to the run log.
#   * Size-verified: the remote byte count from the inventory is enforced.
#   * --dry-run lists every catalogued target and exits 0 without touching the network.
#
# Targets whose remote size exceeds the per-target budget (20 GiB) are recorded as
# skipped, never fetched. The embargoed full archaeological metadata (datafile
# 13572552) is probed for status only; HTTP 403 is the expected, respected result.
#
# Outcomes are appended to data/raw/download_log.tsv.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

MAX_BYTES=$((20 * 1024 * 1024 * 1024))   # 20 GiB per-target budget
LOG="data/raw/download_log.tsv"
DATAVERSE="https://dataverse.harvard.edu/api/access/datafile"
SPRINGER="https://media.springernature.com/original/springer-static/esm"

DRY_RUN=0
case "${1:-}" in
  --dry-run) DRY_RUN=1 ;;
  "" ) ;;
  * ) echo "usage: $0 [--dry-run]" >&2; exit 2 ;;
esac

# Catalog: kind | name | url | dest | expected_bytes | md5 | note
#   kind = download | meta | skip | probe
#   expected_bytes / md5 = "-" when not catalogued.
read -r -d '' CATALOG <<EOF || true
download|akbari_selection_summary_stats|$DATAVERSE/12126280|data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz|674584607|a527fd0b36c741ffe5cc480f1626f666|Akbari Dataverse V4.0 CC0 genome-wide selection statistics
download|akbari_supp_MOESM1_pdf|$SPRINGER/art%3A10.1038%2Fs41586-026-10358-1/MediaObjects/41586_2026_10358_MOESM1_ESM.pdf|data/raw/akbari/41586_2026_10358_MOESM1_ESM.pdf|16629994|-|Akbari supplementary PDF (publisher rights)
download|akbari_supp_MOESM2_pdf|$SPRINGER/art%3A10.1038%2Fs41586-026-10358-1/MediaObjects/41586_2026_10358_MOESM2_ESM.pdf|data/raw/akbari/41586_2026_10358_MOESM2_ESM.pdf|2431382|-|Akbari reporting summary PDF
download|akbari_supp_MOESM3_xlsx|$SPRINGER/art%3A10.1038%2Fs41586-026-10358-1/MediaObjects/41586_2026_10358_MOESM3_ESM.xlsx|data/raw/akbari/41586_2026_10358_MOESM3_ESM.xlsx|5365772|-|Akbari tables 1-5 XLSX
download|aadr_v66p1_1240K_geno|$DATAVERSE/13994829|data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.geno|7117276654|-|AADR V14.0/v66.p1 EIGENSTRAT 1240K genotypes
download|aadr_v66p1_1240K_ind|$DATAVERSE/13994513|data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.ind|1017460|-|AADR 1240K individual file
download|aadr_v66p1_1240K_snp|$DATAVERSE/13994514|data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.snp|77679819|-|AADR 1240K SNP file
download|aadr_v66p1_1240K_anno|$DATAVERSE/13994515|data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.anno|13443726|-|AADR 1240K annotation (dates/locations)
download|aadr_v66p1_HO_geno|$DATAVERSE/13994808|data/raw/aadr-v66.p1/v66.p1_HO.aadr.patch.PUB.geno|4029634650|-|AADR v66.p1 HO EIGENSTRAT genotypes
download|aadr_v66p1_HO_ind|$DATAVERSE/13994526|data/raw/aadr-v66.p1/v66.p1_HO.aadr.patch.PUB.ind|1176253|-|AADR HO individual file
download|aadr_v66p1_HO_snp|$DATAVERSE/13994527|data/raw/aadr-v66.p1/v66.p1_HO.aadr.patch.PUB.snp|36800253|-|AADR HO SNP file
download|aadr_v66p1_HO_anno|$DATAVERSE/13994528|data/raw/aadr-v66.p1/v66.p1_HO.aadr.PUB.anno|15465906|-|AADR HO annotation
download|amy_supp_MOESM3_xlsx|$SPRINGER/art%3A10.1038%2Fs41586-024-07911-1/MediaObjects/41586_2024_7911_MOESM3_ESM.xlsx|data/raw/amy/41586_2024_7911_MOESM3_ESM.xlsx|544662|-|AMY (Bolognini) Tables S1-S11, CC BY 4.0
download|amy_zenodo_v1.1_zip|https://zenodo.org/api/records/12753060/files/sudmantlab/amylase_diversity_project-v1.1.zip/content|data/raw/amy/amylase_diversity_project-v1.1.zip|126578650|3c841ce130eaec8c78b160cb6b971719|AMY versioned analysis archive, CC BY 4.0
download|amy_copy_number_table|https://raw.githubusercontent.com/sudmantlab/amylase_diversity_project/v1.1/graph_genotyping/assess_graph_genotypes/supplementary_tables/copy_number_table_for_supplement.tsv|data/raw/amy/copy_number_table_for_supplement.tsv|-|-|AMY copy-number table TSV
meta|ena_prjeb106907_manifest|https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB106907&result=read_run&fields=run_accession,sample_accession,fastq_ftp,fastq_bytes,submitted_ftp,submitted_bytes&limit=0|data/raw/akbari/ENA_PRJEB106907_runs.tsv|-|-|ENA run manifest (metadata only)
meta|zenodo_record_12753060|https://zenodo.org/api/records/12753060|data/raw/amy/Zenodo_12753060_record.json|-|-|Zenodo record metadata
meta|gtex_access_instructions|https://gtexportal.org/home/protectedDataAccess|data/raw/amy/GTEx_access_instructions.html|-|-|GTEx protected-access instructions (no protected data)
skip|akbari_g2_bcf_aggregate|https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/7RVV9N|data/raw/akbari/AncientDNA_AllVariants_NoQC_GRCh37_chr1.bcf|1324569749458|-|NO-QC per-chromosome BCF aggregate exceeds 20 GiB budget
skip|akbari_g3_pgen|$DATAVERSE/13415552|data/raw/akbari/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.pgen|26300030448|-|PLINK2 combined pgen exceeds 20 GiB budget
probe|akbari_m2_embargoed_metadata|$DATAVERSE/13572552|-|7346667|-|Embargoed until 2041-01-01; status probe only, expect HTTP 403
EOF

say() { printf '%s\n' "$*"; }

mkdir -p data/raw/akbari data/raw/aadr-v66.p1 data/raw/amy

if [[ $DRY_RUN -eq 1 ]]; then
  say "DRY RUN -- catalogued targets (no network, no writes):"
  printf '%-28s %-10s %-12s %s\n' "NAME" "KIND" "BYTES" "DEST/URL"
  n=0
  while IFS='|' read -r kind name url dest size md5 note; do
    [[ -z "$kind" ]] && continue
    n=$((n+1))
    target="$dest"; [[ "$dest" == "-" ]] && target="$url"
    if [[ "$kind" == "skip" && "$size" != "-" ]]; then
      target="$target  [SKIP > 20 GiB budget]"
    fi
    printf '%-28s %-10s %-12s %s\n' "$name" "$kind" "$size" "$target"
    say "    note: $note"
  done <<< "$CATALOG"
  say ""
  say "Total catalogued targets: $n"
  exit 0
fi

: > "$LOG"
printf 'timestamp_utc\tname\tkind\tstatus\thttp_code\tbytes\texpected_bytes\tsha256\tmd5\tdest\turl\tdetail\n' >> "$LOG"
TS() { date -u +%Y-%m-%dT%H:%M:%SZ; }

sha256_of() { sha256sum "$1" | cut -d' ' -f1; }
md5_of()    { md5sum "$1"    | cut -d' ' -f1; }
size_of()   { wc -c < "$1" | tr -d ' '; }

# Probe a URL for its HTTP status and a short body excerpt (for failure records).
probe_http() {
  local url="$1" code body
  code="$(curl -sS -o /tmp/dl_probe_body.$$ -w '%{http_code}' --max-time 30 "$url" 2>/dev/null)" || code="000"
  body="$(head -c 180 /tmp/dl_probe_body.$$ 2>/dev/null | tr '\t\r\n' '   ')"
  rm -f /tmp/dl_probe_body.$$
  printf '%s\t%s' "$code" "$body"
}

fail=0
while IFS='|' read -r kind name url dest size md5 note; do
  [[ -z "$kind" ]] && continue

  if [[ "$kind" == "skip" ]]; then
    say "SKIP  $name -- $note"
    printf '%s\t%s\t%s\tskipped-too-large\t-\t0\t%s\t-\t-\t%s\t%s\t%s\n' \
      "$(TS)" "$name" "$kind" "$size" "$dest" "$url" "$note" >> "$LOG"
    continue
  fi

  if [[ "$kind" == "probe" ]]; then
    read -r code body < <(probe_http "$url"; echo)
    status="probe-embargo-$code"
    [[ "$code" == "403" ]] && status="probe-embargo-403-expected"
    say "PROBE $name -> HTTP $code"
    printf '%s\t%s\t%s\t%s\t%s\t0\t%s\t-\t-\t%s\t%s\t%s\n' \
      "$(TS)" "$name" "$kind" "$status" "$code" "$size" "-" "$url" "$body" >> "$LOG"
    continue
  fi

  # kind is download or meta
  if [[ -f "$dest" ]]; then
    actual="$(size_of "$dest")"; sh="$(sha256_of "$dest")"
    ok=1
    [[ "$size" != "-" && "$actual" != "$size" ]] && ok=0
    if [[ "$md5" != "-" && "$(md5_of "$dest")" != "$md5" ]]; then ok=0; fi
    if [[ $ok -eq 1 ]]; then
      say "OK    $name already present and verified ($actual bytes)"
      printf '%s\t%s\t%s\tpresent-verified\t-\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(TS)" "$name" "$kind" "$actual" "$size" "$sh" "$md5" "$dest" "$url" "existing file left untouched" >> "$LOG"
    else
      say "WARN  $name present but failed verification; leaving untouched (never overwrite)"
      printf '%s\t%s\t%s\tpresent-mismatch\t-\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(TS)" "$name" "$kind" "$actual" "$size" "$sh" "$md5" "$dest" "$url" "size/checksum mismatch on existing file" >> "$LOG"
      fail=1
    fi
    continue
  fi

  part="$dest.part"
  say "GET   $name -> $dest"
  curl -fSL -C - --retry 3 --retry-delay 3 --connect-timeout 30 "$url" -o "$part"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    read -r code body < <(probe_http "$url"; echo)
    say "FAIL  $name (curl rc=$rc, probe HTTP $code) -- partial kept at $part"
    printf '%s\t%s\t%s\tfailed\t%s\t%s\t%s\t-\t-\t%s\t%s\tcurl rc=%s; body=%s\n' \
      "$(TS)" "$name" "$kind" "$code" "$(size_of "$part" 2>/dev/null || echo 0)" "$size" "$dest" "$url" "$rc" "$body" >> "$LOG"
    fail=1
    continue
  fi

  actual="$(size_of "$part")"
  if [[ "$size" != "-" && "$actual" != "$size" ]]; then
    say "FAIL  $name size mismatch: got $actual expected $size -- partial kept"
    printf '%s\t%s\t%s\tfailed-size\t200\t%s\t%s\t-\t-\t%s\t%s\tbyte count mismatch\n' \
      "$(TS)" "$name" "$kind" "$actual" "$size" "$dest" "$url" >> "$LOG"
    fail=1
    continue
  fi
  if [[ "$md5" != "-" && "$(md5_of "$part")" != "$md5" ]]; then
    say "FAIL  $name md5 mismatch -- partial kept"
    printf '%s\t%s\t%s\tfailed-checksum\t200\t%s\t%s\t-\t%s\t%s\t%s\tcatalogued md5 mismatch\n' \
      "$(TS)" "$name" "$kind" "$actual" "$size" "$md5" "$dest" "$url" >> "$LOG"
    fail=1
    continue
  fi

  mv -n "$part" "$dest"
  if [[ ! -f "$dest" ]]; then
    say "FAIL  $name could not finalize (dest appeared) -- not overwriting"
    printf '%s\t%s\t%s\tfailed-finalize\t200\t%s\t%s\t-\t-\t%s\t%s\tdest existed at finalize\n' \
      "$(TS)" "$name" "$kind" "$actual" "$size" "$dest" "$url" >> "$LOG"
    fail=1
    continue
  fi

  sh="$(sha256_of "$dest")"
  say "DONE  $name ($actual bytes) sha256=$sh"
  printf '%s\t%s\t%s\tdownloaded-verified\t200\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(TS)" "$name" "$kind" "$actual" "$size" "$sh" "$md5" "$dest" "$url" "$note" >> "$LOG"
done <<< "$CATALOG"

say ""
say "Run log: $LOG"
exit $fail
