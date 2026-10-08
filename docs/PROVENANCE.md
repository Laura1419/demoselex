# Provenance policy

Every run records: git commit, dirty flag, command line, software versions, input
file hashes, parameters, random seed, outputs. Written to
`analysis/runs/<UTC timestamp>-<slug>/manifest.json` by `scripts/provenance.py`.

Raw data under `data/raw/` is immutable. All derived data is regenerable.
