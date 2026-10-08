#!/usr/bin/env python3
"""
Provenance tracking for simulation runs
Records seed, versions, and execution context for reproducible results
"""

import json
import datetime
import platform
import sys
import subprocess
import hashlib
import os
from pathlib import Path
from typing import Dict, Any, Optional

class Provenance:
    """Track simulation provenance"""

    def __init__(self, seed: int, config_file: str, output_dir: str):
        self.seed = seed
        self.config_file = config_file
        self.output_dir = output_dir
        self.start_time = datetime.datetime.now().isoformat()
        self.versions = self._get_versions()
        self.system_info = self._get_system_info()

    def _get_versions(self) -> Dict[str, str]:
        """Get versions of key dependencies"""
        versions = {}

        # Python version
        versions["python"] = platform.python_version()

        # Try to get SLiM version if available
        try:
            result = subprocess.run(["slim", "-v"], capture_output=True, text=True)
            if result.returncode == 0:
                versions["slim"] = result.stdout.strip()
        except FileNotFoundError:
            versions["slim"] = "not_available"

        # Git commit hash if in git repository
        try:
            result = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                cwd=Path(__file__).parent.parent
            )
            if result.returncode == 0:
                versions["git_commit"] = result.stdout.strip()
        except Exception:
            versions["git_commit"] = "unknown"

        return versions

    def _get_system_info(self) -> Dict[str, str]:
        """Get system information"""
        return {
            "platform": platform.platform(),
            "processor": platform.processor() or "unknown",
            "hostname": platform.node()
        }

    def record_run(self, output_files: Dict[str, str], params: Dict[str, Any]) -> str:
        """Record a simulation run and return provenance file path"""
        end_time = datetime.datetime.now().isoformat()

        # Calculate file checksums
        file_checksums = {}
        for file_type, file_path in output_files.items():
            if os.path.exists(file_path):
                with open(file_path, "rb") as f:
                    file_checksums[file_type] = hashlib.sha256(f.read()).hexdigest()

        provenance_data = {
            "seed": self.seed,
            "config_file": self.config_file,
            "start_time": self.start_time,
            "end_time": end_time,
            "parameters": params,
            "output_files": output_files,
            "file_checksums": file_checksums,
            "versions": self.versions,
            "system_info": self.system_info
        }

        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

        # Write provenance file
        provenance_file = os.path.join(
            self.output_dir,
            f"provenance_{self.seed}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )

        with open(provenance_file, 'w') as f:
            json.dump(provenance_data, f, indent=2)

        return provenance_file

def main():
    """Example usage"""
    if len(sys.argv) < 4:
        print("Usage: python3 provenance.py <seed> <config_file> <output_dir>")
        sys.exit(1)

    seed = int(sys.argv[1])
    config_file = sys.argv[2]
    output_dir = sys.argv[3]

    prov = Provenance(seed, config_file, output_dir)

    # Example output files (would be created by actual simulation)
    output_files = {
        "genotype_matrix": "output/genotypes.csv",
        "truth_table": "output/truth.csv",
        "config": config_file
    }

    # Example parameters
    params = {
        "population_size": 1000,
        "migration_rate": 0.01,
        "selection_strength": 0.1
    }

    provenance_file = prov.record_run(output_files, params)
    print(f"Provenance recorded to: {provenance_file}")

if __name__ == "__main__":
    main()