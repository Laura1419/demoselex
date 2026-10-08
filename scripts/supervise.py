#!/usr/bin/env python3
"""Lightweight supervisor loop: retry routinely-failed tasks, report blockers.

Safe by construction: only issues `wg retry` for failed tasks whose reported
attempt count is below a cap, and never mutates shared git state.

Usage: python3 scripts/supervise.py [--max-attempts N] [--dry-run]
"""
from __future__ import annotations

import json
import re
import subprocess
import sys

ROOT = "/home/erikg/demoselex"


def sh(cmd: list[str]) -> str:
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True).stdout


def tasks() -> list[dict]:
    raw = sh(["wg", "list", "--all", "--json"])
    try:
        d = json.loads(raw)
    except Exception:
        return []
    if isinstance(d, dict):
        return d.get("tasks", [])
    return d if isinstance(d, list) else []


def main() -> int:
    cap = 6
    dry = "--dry-run" in sys.argv
    if "--max-attempts" in sys.argv:
        cap = int(sys.argv[sys.argv.index("--max-attempts") + 1])

    rows = tasks()
    retried, failed_high, waiting, done = [], [], [], []
    for r in rows:
        tid = r.get("id") or ""
        if tid.startswith("."):
            continue
        st = (r.get("status") or "").lower()
        att = r.get("attempts") or r.get("attempt_count") or 0
        if st == "failed":
            if att < cap:
                retried.append((tid, att))
            else:
                failed_high.append((tid, att))
        elif st == "waiting":
            waiting.append(tid)
        elif st == "done":
            done.append(tid)

    print(f"done={len(done)} waiting={len(waiting)} failed-retryable={len(retried)} "
          f"failed-exhausted={len(failed_high)}")

    for tid, att in retried:
        if dry:
            print(f"  would retry {tid} (attempts={att})")
        else:
            out = sh(["wg", "retry", tid])
            print(f"  retry {tid}: {out.strip().splitlines()[-1] if out.strip() else ''}")

    for tid, att in failed_high:
        print(f"  EXHAUSTED {tid} (attempts={att}) -> needs operator")
    for tid in waiting:
        print(f"  WAITING {tid} -> landing/attention needs operator")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
