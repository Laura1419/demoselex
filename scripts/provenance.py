#!/usr/bin/env python3
"""Provenance recording for trajectory-replication research project.

Records git commit, software versions, parameters, and output paths for any run.
"""

import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

def get_lean_version():
    """Get the installed Lean version"""
    try:
        result = subprocess.run(['lean', '--version'], capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except Exception as e:
        return f"Error getting Lean version: {e}"

def get_lake_version():
    """Get the installed Lake version"""
    try:
        result = subprocess.run(['lake', '--version'], capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except Exception as e:
        return f"Error getting Lake version: {e}"

def record_provenance():
    """Record current environment and software versions"""
    provenance = {
        "timestamp": datetime.now().isoformat(),
        "lean_version": get_lean_version(),
        "lake_version": get_lake_version(),
        "project_sha": subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            capture_output=True, text=True, check=True
        ).stdout.strip(),
        "dependencies": {
            "mathlib": "pinned to master",
            "target": "lean 4.34.1"
        }
    }

    # Write to file
    prov_file = Path("data/derived/provenance.json")
    prov_file.parent.mkdir(parents=True, exist_ok=True)

    with open(prov_file, 'w') as f:
        json.dump(provenance, f, indent=2)

    print(f"Provenance recorded to {prov_file}")
    print(json.dumps(provenance, indent=2))

if __name__ == "__main__":
    record_provenance()