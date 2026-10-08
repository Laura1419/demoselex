# F. Reproducible code

## Demographic validation simulator (`sim-scaffold`)

The executable entry point is `sim/scripts/run_simulations.py`; the six standalone,
fixed-seed YAML configurations are `sim/scenarios/*.yaml`. For example, from the
repository root:

```bash
python3 sim/scripts/run_simulations.py sim/scenarios sim/output/demo
python3 -m unittest discover -s sim/scripts -p 'test_simulator.py' -v
```

Dependencies: Python 3, NumPy, PyYAML (the provenance harness also probes
optional packages). SLiM is not installed on the development host; the
implementation is a forward-time diploid Wright–Fisher **mosaic** simulator,
not a coalescent approximation or SLiM run. Each generation samples two
fitness-weighted parents' recombinant gametes for each offspring; haplotype
switches occur between adjacent loci at the configured probability. Populations
advance synchronously. Continuous migration replaces randomly chosen diploids
with donor diploids from the *previous* generation, and replacement is a
separate, dated pulse of donor ancestry. This distinguishes locus-specific
selection from pure ancestry turnover. Founders share ancestral allele
frequencies but have population-specific drift; optional locus overrides
create known pre-turnover frequency contrasts. Each config specifies an
explicit seed; randomness uses one local `numpy.random.Generator`, never
Python's process-randomized `hash()`.

Samples from `sample_generations` are individual genotype snapshots. `years_bp`
is `(final_generation - sampled_generation) * generation_years`, so larger
values are older. Each scenario directory under the requested output parent
contains:

- `genotypes.csv`: one row per sampled diploid, `sample_id` plus ordered
  `locus_0..locus_N` allele counts in {0,1,2};
- `samples.csv`: the same sample IDs/order, trajectory, generation, and years BP;
- `truth_table.csv`: **every locus × trajectory**, signed additive effect and
  selected indicator (zero for unselected loci); no selection is inferred
  merely from a replacement pulse;
- `summary.json`: scenario type, seed, matrix dimensions;
- `provenance.json`: copy of the `scripts/provenance.py` `RunRecorder` manifest,
  with seed, git revision, package versions, config hash, parameters, and output
  hashes. The canonical manifest also lives in `analysis/runs/`.

The six configs cover drift with migration, replacement without selection,
selection in one population, shared same-direction selection, opposing signed
selection, and selection in the incoming population aligned with turnover.
This is a controlled **validation scaffold**, not calibrated human population
genetics: no mutation, ascertainment, genotyping error, realistic recombination
map, variable effective size, ancient-DNA damage, or pedigree validation.
Selection effects are additive on individual fitness and can be obscured by
finite-population drift; do not interpret one realized trajectory as proof of
statistical power. `sim-sweeps` should vary seeds, sizes and rates and check
null calibration and identifiability before applying the method to real data.
