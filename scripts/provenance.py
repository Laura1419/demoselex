#!/usr/bin/env python3
"""
Provenance and run-manifest harness for the trajectory-replication project.

Records: git commit, command line, software versions, input file hashes,
parameters, random seed, and outputs for every run.

Writes to: analysis/runs/<timestamp>-<slug>/manifest.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union


class ProvenanceError(Exception):
    """Custom exception for provenance-related errors."""
    pass


def get_git_info() -> Dict[str, str]:
    """Get current git commit and branch information."""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], 
            text=True, stderr=subprocess.DEVNULL
        ).strip()
        
        # Try to get branch name
        try:
            branch = subprocess.check_output(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                text=True, stderr=subprocess.DEVNULL
            ).strip()
        except subprocess.CalledProcessError:
            branch = "detached"
        
        # Check if working tree is dirty
        try:
            subprocess.check_output(
                ["git", "diff-index", "--quiet", "HEAD", "--"],
                stderr=subprocess.DEVNULL
            )
            dirty = False
        except subprocess.CalledProcessError:
            dirty = True
            
        return {
            "commit": commit,
            "branch": branch,
            "dirty": dirty
        }
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {"commit": "unknown", "branch": "unknown", "dirty": False}


def get_software_versions() -> Dict[str, str]:
    """Get versions of key software components."""
    versions = {
        "python": platform.python_version(),
        "system": f"{platform.system()} {platform.release()}"
    }
    
    # Try to get common package versions
    try:
        import numpy
        versions["numpy"] = numpy.__version__
    except ImportError:
        versions["numpy"] = "not installed"
    
    try:
        import pandas
        versions["pandas"] = pandas.__version__
    except ImportError:
        versions["pandas"] = "not installed"
        
    try:
        import scipy
        versions["scipy"] = scipy.__version__
    except ImportError:
        versions["scipy"] = "not installed"
        
    try:
        import sklearn
        versions["sklearn"] = sklearn.__version__
    except ImportError:
        versions["sklearn"] = "not installed"
        
    return versions


def file_hash(filepath: Union[str, Path], algorithm: str = "sha256") -> str:
    """Calculate hash of a file."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    hasher = hashlib.new(algorithm)
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    
    return f"{algorithm}:{hasher.hexdigest()}"


class RunRecorder:
    """Record provenance information for a single run."""
    
    def __init__(self, slug: str, command: Optional[List[str]] = None):
        self.slug = slug
        self.command = command or sys.argv
        self.timestamp = datetime.now().isoformat()
        self.run_id = f"{datetime.now().strftime('%Y%m%d-%H%M')}-{slug}"
        
        self.manifest: Dict[str, Any] = {
            "run_id": self.run_id,
            "timestamp": self.timestamp,
            "git": get_git_info(),
            "command": self.command,
            "versions": get_software_versions(),
            "inputs": {},
            "parameters": {},
            "outputs": {},
            "seed": None
        }
        
        # Create run directory
        self.run_dir = Path("analysis") / "runs" / self.run_id
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
    def record_input(self, filepath: Union[str, Path]) -> str:
        """Record an input file and return its hash."""
        filepath = Path(filepath)
        try:
            file_hash_str = file_hash(filepath)
            self.manifest["inputs"][str(filepath)] = file_hash_str
            return file_hash_str
        except FileNotFoundError:
            # File might be created later, record path anyway
            self.manifest["inputs"][str(filepath)] = "not_found"
            return "not_found"
    
    def record_output(self, filepath: Union[str, Path]) -> str:
        """Record an output file and return its hash."""
        filepath = Path(filepath)
        if filepath.exists():
            file_hash_str = file_hash(filepath)
            self.manifest["outputs"][str(filepath)] = file_hash_str
            return file_hash_str
        else:
            raise FileNotFoundError(f"Output file does not exist: {filepath}")
    
    def record_parameter(self, key: str, value: Any) -> None:
        """Record a parameter value."""
        self.manifest["parameters"][key] = value
    
    def set_seed(self, seed: Optional[int] = None) -> int:
        """Set and record a random seed."""
        if seed is None:
            # Generate a deterministic seed based on run_id
            seed = int(hashlib.sha256(self.run_id.encode()).hexdigest()[:8], 16) % (2**31)
        
        self.manifest["seed"] = seed
        return seed
    
    def save(self) -> Path:
        """Save the manifest to disk."""
        manifest_path = self.run_dir / "manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(self.manifest, f, indent=2, default=str)
        
        print(f"Provenance manifest saved to: {manifest_path}")
        return manifest_path
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.save()
        # Don't suppress exceptions
        return False


def run_with_provenance(
    slug: str,
    inputs: List[Union[str, Path]],
    parameters: Dict[str, Any],
    command_fn: Callable,
    seed: Optional[int] = None
) -> Any:
    """
    Run a function with provenance recording.
    
    Args:
        slug: Short identifier for this run
        inputs: List of input file paths
        parameters: Dictionary of parameter values
        command_fn: Function to execute
        seed: Optional random seed (auto-generated if None)
        
    Returns:
        Result of command_fn
    """
    with RunRecorder(slug) as recorder:
        # Record inputs
        for input_path in inputs:
            recorder.record_input(input_path)
        
        # Record parameters
        for key, value in parameters.items():
            recorder.record_parameter(key, value)
        
        # Set seed
        configured_seed = recorder.set_seed(seed)
        
        # Execute the function
        result = command_fn()
        
        return result


def demo_run() -> None:
    """Run a demo to test the provenance system."""
    print("Running provenance demo...")
    
    # Create a test directory structure
    test_dir = Path("analysis") / "runs" / "test-demo"
    test_dir.mkdir(parents=True, exist_ok=True)
    
    # Create some test files
    test_input = test_dir / "test_input.txt"
    test_output = test_dir / "test_output.txt"
    
    with open(test_input, "w") as f:
        f.write("Test input data for provenance demo\n")
    
    with RunRecorder("demo-run") as recorder:
        recorder.record_input(test_input)
        recorder.record_parameter("demo_param", 42)
        recorder.record_parameter("another_param", "test_value")
        
        # Generate a simple output
        with open(test_output, "w") as f:
            f.write(f"Output from demo run at {datetime.now()}\n")
        
        recorder.record_output(test_output)
    
    print("Demo complete. Check the generated manifest.")


def verify_manifest(manifest_path: Union[str, Path]) -> bool:
    """Verify that a manifest's inputs and outputs still match their hashes."""
    manifest_path = Path(manifest_path)
    
    if not manifest_path.exists():
        print(f"Manifest not found: {manifest_path}")
        return False
    
    with open(manifest_path) as f:
        manifest = json.load(f)
    
    print(f"Verifying manifest: {manifest['run_id']}")
    print(f"Timestamp: {manifest['timestamp']}")
    print(f"Git commit: {manifest['git']['commit']}")
    
    all_verified = True
    
    # Verify inputs
    print("\nInput files:")
    for filepath, expected_hash in manifest["inputs"].items():
        if expected_hash == "not_found":
            print(f"  {filepath}: not found during recording")
            continue
            
        if Path(filepath).exists():
            current_hash = file_hash(filepath)
            if current_hash == expected_hash:
                print(f"  {filepath}: ✓")
            else:
                print(f"  {filepath}: ✗ (hash mismatch)")
                all_verified = False
        else:
            print(f"  {filepath}: ✗ (file missing)")
            all_verified = False
    
    # Verify outputs
    print("\nOutput files:")
    for filepath, expected_hash in manifest["outputs"].items():
        if Path(filepath).exists():
            current_hash = file_hash(filepath)
            if current_hash == expected_hash:
                print(f"  {filepath}: ✓")
            else:
                print(f"  {filepath}: ✗ (hash mismatch)")
                all_verified = False
        else:
            print(f"  {filepath}: ✗ (file missing)")
            all_verified = False
    
    if all_verified:
        print("\n✓ All files verified successfully")
    else:
        print("\n✗ Some files failed verification")
    
    return all_verified


def main():
    parser = argparse.ArgumentParser(description="Provenance and run-manifest harness")
    subparsers = parser.add_subparsers(dest="subcommand", help="Command to run")
    
    # Record command
    record_parser = subparsers.add_parser("record", help="Record a new run")
    record_parser.add_argument("--slug", required=True, help="Short identifier for the run")
    record_parser.add_argument("--input", action="append", dest="inputs", 
                              help="Input file to record (can be used multiple times)")
    record_parser.add_argument("--parameter", action="append", dest="parameters",
                              help="Parameter as key=value (can be used multiple times)")
    record_parser.add_argument("--seed", type=int, help="Random seed")
    record_parser.add_argument("--cmd", dest="cmd_args", help="Command being executed")
    
    # Verify command
    verify_parser = subparsers.add_parser("verify", help="Verify a manifest")
    verify_parser.add_argument("manifest", help="Path to manifest file")
    
    # Demo command
    subparsers.add_parser("demo", help="Run a demo of the provenance system")
    
    args = parser.parse_args()
    
    if args.subcommand == "demo":
        demo_run()
    
    elif args.subcommand == "verify":
        verify_manifest(args.manifest)
    
    elif args.subcommand == "record":
        cmd_list = [args.cmd_args] if args.cmd_args else None
        with RunRecorder(args.slug, cmd_list) as recorder:
            # Record inputs
            if args.inputs:
                for input_path in args.inputs:
                    recorder.record_input(input_path)
            
            # Record parameters
            if args.parameters:
                for param_str in args.parameters:
                    if "=" in param_str:
                        key, value = param_str.split("=", 1)
                        recorder.record_parameter(key, value)
                    else:
                        print(f"Warning: ignoring malformed parameter: {param_str}")
            
            # Set seed
            if args.seed is not None:
                recorder.set_seed(args.seed)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()