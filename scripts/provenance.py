#!/usr/bin/env python3
"""Provenance recording for trajectory-replication research project.

Records git commit, software versions, parameters, and output paths for any run.
"""

import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional


def get_git_info() -> Dict[str, Optional[str]]:
    """Get git commit hash and status."""
    try:
        commit_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], 
            text=True, 
            stderr=subprocess.DEVNULL
        ).strip()
    except subprocess.CalledProcessError:
        commit_hash = None
    
    try:
        status = subprocess.check_output(
            ["git", "status", "--porcelain"], 
            text=True, 
            stderr=subprocess.DEVNULL
        ).strip()
        is_clean = len(status) == 0
    except subprocess.CalledProcessError:
        is_clean = False
    
    return {"commit": commit_hash, "clean": is_clean}


def record_provenance(
    command: str,
    parameters: Dict[str, Any],
    inputs: Dict[str, str],
    outputs: Dict[str, str],
    seed: Optional[int] = None,
    output_dir: Optional[Path] = None
) -> Path:
    """Record provenance for a run.
    
    Args:
        command: The command that was run
        parameters: Key parameters used
        inputs: Input files with hashes or paths
        outputs: Output files/directories
        seed: Random seed if applicable
        output_dir: Where to write manifest (defaults to cwd)
        
    Returns:
        Path to the manifest file
    """
    if output_dir is None:
        output_dir = Path.cwd()
    
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    # Create a slug from the command
    slug = command.split()[0].replace("/", "-").replace(".", "-")
    run_dir = output_dir / f"{timestamp}-{slug}"
    run_dir.mkdir(parents=True, exist_ok=True)
    
    git_info = get_git_info()
    
    manifest = {
        "timestamp": datetime.now().isoformat(),
        "git": git_info,
        "command": command,
        "parameters": parameters,
        "inputs": inputs,
        "outputs": outputs,
        "seed": seed,
        "host": os.uname().nodename,
        "python_version": sys.version,
    }
    
    manifest_path = run_dir / "manifest.json"
    with manifest_path.open("w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    
    return manifest_path


def main():
    """CLI interface for recording provenance."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Record provenance for a run")
    parser.add_argument("--command", required=True, help="Command that was run")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--seed", type=int, help="Random seed")
    parser.add_argument("--parameter", nargs=2, action="append", 
                       metavar=("KEY", "VALUE"), help="Parameter key-value pairs")
    parser.add_argument("--input", nargs=2, action="append",
                       metavar=("NAME", "PATH"), help="Input files")
    parser.add_argument("--output", nargs=2, action="append",
                       metavar=("NAME", "PATH"), help="Output files")
    
    args = parser.parse_args()
    
    parameters = dict(args.parameter) if args.parameter else {}
    inputs = dict(args.input) if args.input else {}
    outputs = dict(args.output) if args.output else {}
    
    manifest_path = record_provenance(
        command=args.command,
        parameters=parameters,
        inputs=inputs,
        outputs=outputs,
        seed=args.seed,
        output_dir=Path(args.output_dir)
    )
    
    print(f"Provenance recorded: {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())