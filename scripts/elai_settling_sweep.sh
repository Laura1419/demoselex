#!/usr/bin/env bash
# Trace-only, bounded parallel sweep. The parent, not this worker, launches full mode.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
PY="$ROOT/.venv/bin/python"
BINARY=/home/erikg/elai-fork/src/elai
BASE="$ROOT/data/derived/elai_amy_dense"
INPUT="$BASE/work_long/panel.inp"
POS="$BASE/work_long/panel.pos"
PREFIX=sweep
ONLY=
while (($#)); do
    case "$1" in
        --smoke) PREFIX=sweep_smoke; shift ;;
        --only)
            [[ $# -ge 2 ]] || { echo '--only needs a configuration name' >&2; exit 2; }
            ONLY=$2; shift 2 ;;
        -h|--help)
            echo 'Usage: bash scripts/elai_settling_sweep.sh [--smoke [--only NAME]]'
            echo 'Full mode runs all registered arms; smoke runs 5 exact steps on the same dense panel.'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -z "$ONLY" || "$PREFIX" == sweep_smoke ]] || { echo '--only is allowed only with --smoke' >&2; exit 2; }
for file in "$PY" "$BINARY" "$INPUT" "$POS" "$ROOT/docs/elai-settling-sweep.md"; do
    [[ -f "$file" ]] || { echo "Missing required file: $file" >&2; exit 2; }
done
[[ -x "$PY" && -x "$BINARY" ]] || { echo 'Python/binary must be executable' >&2; exit 2; }
for command in timeout setsid; do
    command -v "$command" >/dev/null || { echo "$command is required" >&2; exit 2; }
done
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PLAN=$("$PY" -B analysis/elai_settling_sweep_report.py --plan)
NAMES=(); UPPERS=(); LOWERS=(); WARMS=(); MGS=(); STEPS=(); SEEDS=()
while IFS=$'\t' read -r name upper lower warm mg steps seed; do
    [[ -z "$ONLY" || "$name" == "$ONLY" ]] || continue
    [[ "$PREFIX" != sweep_smoke ]] || steps=5
    NAMES+=("$name"); UPPERS+=("$upper"); LOWERS+=("$lower"); WARMS+=("$warm")
    MGS+=("$mg"); STEPS+=("$steps"); SEEDS+=("$seed")
done <<< "$PLAN"
((${#NAMES[@]})) || { echo "Unknown configuration: $ONLY" >&2; exit 2; }
STATUS="$BASE/${PREFIX}.status.tsv"
REPORT="$BASE/${PREFIX}.report.json"
TEXT="$BASE/${PREFIX}.report.txt"
# Serialize reservation/launch, including smoke vs full. Never remove another run's lock.
LOCK="$BASE/settling_sweep.lock"
mkdir "$LOCK" 2>/dev/null || { echo "Sweep lock exists: $LOCK (inspect before removing)" >&2; exit 2; }
PIDS=()
cleanup() {
    local pid
    for pid in "${PIDS[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then kill -TERM "$pid" 2>/dev/null || true; wait "$pid" 2>/dev/null || true; fi
    done
    rmdir "$LOCK"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
for file in "$STATUS" "$REPORT" "$TEXT"; do
    [[ ! -e "$file" && ! -L "$file" ]] || { echo "Refusing existing artifact: $file" >&2; exit 2; }
done
for name in "${NAMES[@]}"; do
    for file in "$BASE/${PREFIX}_${name}.stdout.txt" "$BASE/${PREFIX}_${name}.status.tsv" "$BASE/${PREFIX}_${name}.work"; do
        [[ ! -e "$file" && ! -L "$file" ]] || { echo "Refusing existing artifact: $file" >&2; exit 2; }
    done
done
# Reserve every stdout before any child starts. noclobber is also race-safe.
set -o noclobber
: > "$STATUS"
for name in "${NAMES[@]}"; do
    : > "$BASE/${PREFIX}_${name}.stdout.txt"
    mkdir -p "$BASE/${PREFIX}_${name}.work/output"
    # Upstream writes these automatically: discarding them prevents dense disk dumps.
    ln -s /dev/null "$BASE/${PREFIX}_${name}.work/output/${PREFIX}_${name}.ps21.txt"
    ln -s /dev/null "$BASE/${PREFIX}_${name}.work/output/${PREFIX}_${name}.em.txt"
done
for i in "${!NAMES[@]}"; do
    name=${NAMES[$i]}; prefix="${PREFIX}_${name}"
    CMD=("$BINARY" -g "$INPUT" -p 1 -pos "$POS" -C "${UPPERS[$i]}" -c "${LOWERS[$i]}"
         -w "${WARMS[$i]}" -s "${STEPS[$i]}" -mg "${MGS[$i]}" -R "${SEEDS[$i]}" -o "$prefix"
         --pkk-site 1:104121000)
    printf 'START %s: ' "$name"; printf '%q ' "${CMD[@]}"; printf '\n'
    printf '# command: ' >> "$BASE/$prefix.stdout.txt"
    printf '%q ' "${CMD[@]}" >> "$BASE/$prefix.stdout.txt"
    printf '\n' >> "$BASE/$prefix.stdout.txt"
    # Each regular ELAI output/stdout is capped at 16 MiB (bash ulimit uses KiB).
    # Each timeout/ELAI has its own process group for safe interrupt cleanup.
    (
        cd "$BASE/$prefix.work"
        ulimit -f 16384
        start=$(date +%s)
        setsid timeout --signal=TERM --kill-after=15s 5400s "${CMD[@]}" &
        child=$!
        trap 'kill -TERM -- "-$child" 2>/dev/null || true; wait "$child" 2>/dev/null || true; exit 143' TERM INT
        rc=0
        wait "$child" || rc=$?
        elapsed=$(( $(date +%s) - start ))
        printf '%s\t%s\t%s\t%s\n' "$name" "${STEPS[$i]}" "$rc" "$elapsed" > "$BASE/$prefix.status.tsv"
        exit "$rc"
    ) >> "$BASE/$prefix.stdout.txt" 2>&1 &
    PIDS+=("$!")
done
FAILED=0
for i in "${!NAMES[@]}"; do
    rc=0
    wait "${PIDS[$i]}" || rc=$?
    name=${NAMES[$i]}
    # ELAI's inherited success exit is 1: require its completion sentinel too.
    state=FAILED
    if [[ "$rc" -le 1 ]] && grep -q -- '-elai: finished, for details see log:' "$BASE/${PREFIX}_${name}.stdout.txt"; then
        state=FINISHED
    else
        FAILED=1
    fi
    if [[ -f "$BASE/${PREFIX}_${name}.status.tsv" ]]; then
        cat "$BASE/${PREFIX}_${name}.status.tsv" >> "$STATUS"
        elapsed=$(cut -f4 "$BASE/${PREFIX}_${name}.status.tsv")
    else
        elapsed=0
        printf '%s\t%s\t%s\t%s\n' "$name" "${STEPS[$i]}" "$rc" "$elapsed" >> "$STATUS"
        FAILED=1
    fi
    printf 'STATUS %s: %s exit=%s elapsed_s=%s\n' "$name" "$state" "$rc" "$elapsed"
done
# All children are reaped; avoid sending a signal to a subsequently recycled PID.
PIDS=()
rc=0
"$PY" -B analysis/elai_settling_sweep_report.py --status "$STATUS" --prefix "$PREFIX" \
    --output "$REPORT" --provenance > "$TEXT" || rc=$?
cat "$TEXT"
((rc == 0)) || FAILED=1
printf 'Report: %s\n' "$REPORT"
exit "$FAILED"
