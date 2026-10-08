#!/usr/bin/env python3
"""One-glance project status for the orchestrator."""
from __future__ import annotations

import json
import subprocess
from collections import Counter, defaultdict

ROOT = "/home/erikg/demoselex"


def sh(cmd: list[str]) -> str:
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True).stdout


def main() -> None:
    raw = sh(["wg", "list", "--all", "--json"])
    try:
        data = json.loads(raw)
    except Exception:
        data = None
    rows: list[dict] = []
    if isinstance(data, list):
        rows = data
    elif isinstance(data, dict):
        rows = data.get("tasks", [])

    by_status: Counter[str] = Counter()
    by_tag: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        st = (r.get("status") or r.get("state") or "?").lower()
        by_status[st] += 1
        for t in r.get("tags", []) or []:
            by_tag[t][st] += 1

    print("== task status ==")
    for st, n in sorted(by_status.items(), key=lambda kv: -kv[1]):
        print(f"  {st:14s} {n}")
    print("\n== tags ==")
    for t, c in sorted(by_tag.items()):
        print(f"  {t:16s} {dict(c)}")

    # git state
    log = sh(["git", "log", "--oneline", "-8"])
    print("\n== recent main commits ==")
    print(log.rstrip())

    # dispatcher
    print("\n== dispatcher ==")
    st = sh(["wg", "service", "status"])
    for line in st.splitlines():
        if any(k in line for k in ("Agents", "Last tick", "admission", "Dispatcher:")):
            print("  " + line.strip())

    # branches from workers
    br = sh(["git", "branch", "--list", "wg/*"])
    print("\n== worker branches ==")
    print(br.rstrip() or "  (none)")


if __name__ == "__main__":
    main()
