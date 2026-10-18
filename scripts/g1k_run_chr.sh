#!/bin/bash
# One chromosome end-to-end: download -> build(train) -> date -> stats -> free intermediates.
set -u
cd /home/erikg/demoselex
c=$1
log=data/derived/g1k/log_chr$c.txt
{
  echo "=== chr$c start $(date +%H:%M:%S)"
  .venv/bin/python -u analysis/g1k_phased.py download $c || exit 1
  .venv/bin/python -u analysis/g1k_phased.py build   $c || exit 1
  rm -f data/derived/g1k/chr$c.vcf.gz          # VCF no longer needed
  .venv/bin/python -u analysis/g1k_phased.py date    $c || exit 1
  .venv/bin/python -u analysis/g1k_phased.py stats   $c || exit 1
  rm -f data/derived/g1k/chr$c.trees           # undated ARG superseded by dated chunks
  echo "=== chr$c done $(date +%H:%M:%S)"
} >> "$log" 2>&1
