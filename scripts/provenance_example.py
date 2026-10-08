#!/usr/bin/env python3
"""
Example usage of the provenance system for analysis scripts.
"""

import sys
from pathlib import Path

# Add scripts directory to path so we can import provenance
sys.path.insert(0, str(Path(__file__).parent))

from provenance import run_with_provenance, RunRecorder


def example_analysis_function(seed: int):
    """Example analysis function that uses the provenance system."""
    import numpy as np
    
    # Set the seed explicitly
    np.random.seed(seed)
    print(f"Running analysis with seed: {seed}")
    
    # Simulate some analysis work
    data = np.random.randn(100)  # This will be deterministic due to the seed
    result = {
        "mean": float(data.mean()),
        "std": float(data.std()),
        "size": len(data)
    }
    
    # Write results to a file
    output_file = Path("analysis/runs/example/output.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    import json
    with open(output_file, "w") as f:
        json.dump(result, f, indent=2)
    
    return result


def run_example_analysis():
    """Run the example analysis with full provenance tracking."""
    import json
    
    # Example input files (create dummy ones for demo)
    input_files = [
        "data/raw/example_genotypes.csv",  # Would be real data in practice
        "config/analysis_config.json"      # Would be real config in practice
    ]
    
    # Create dummy input files for the demo
    Path("data/raw").mkdir(parents=True, exist_ok=True)
    Path("config").mkdir(parents=True, exist_ok=True)
    
    with open("data/raw/example_genotypes.csv", "w") as f:
        f.write("sample_id,chr,pos,allele1,allele2\n")
        f.write("sample1,1,1000,A,T\n")
        f.write("sample2,1,1000,A,A\n")
    
    with open("config/analysis_config.json", "w") as f:
        json.dump({
            "population_threshold": 0.05,
            "min_sample_size": 10,
            "confidence_level": 0.95
        }, f, indent=2)
    
    # Use RunRecorder directly for more control over output recording
    with RunRecorder("example-analysis") as recorder:
        # Record inputs
        for input_path in input_files:
            recorder.record_input(input_path)
        
        # Record parameters
        recorder.record_parameter("population_threshold", 0.05)
        recorder.record_parameter("min_sample_size", 10)
        recorder.record_parameter("confidence_level", 0.95)
        
        # Set seed
        seed = recorder.set_seed()
        
        # Run analysis function
        result = example_analysis_function(seed)
        
        # Explicitly record the output file that was created
        output_file = Path("analysis/runs/example/output.json")
        if output_file.exists():
            recorder.record_output(output_file)
    
    print("Analysis complete!")
    print(f"Result: {result}")
    print("Provenance manifest saved in analysis/runs/")


if __name__ == "__main__":
    run_example_analysis()