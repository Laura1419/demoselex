# Every locus passing the screen (COSI-GT / HGSVC v3, 265 loci)

## What the screen is, in one paragraph

For each locus we take every common SNP in its interval (GRCh38), ask how recently the
people carrying it share a common ancestor, and note whether that is unusually recent
compared with other variants with the same number of carriers in a similar recombination
context. We then count how many of the locus's variants land in the genome-wide bottom 5%
of that measure and compare with the 5% we would expect by chance. A one-sided Poisson
test gives the p, and Benjamini-Hochberg across the 264 scorable loci gives q. **Passing =
q < 0.05**, i.e. more suspiciously-young common variants around this locus than chance
explains, after correcting for having tested 264 loci.

## Result: 39 of 264 loci pass (q < 0.05)

| # | locus | chr | span kb | SV class | variants | enrichment | p | q | most extreme variant (pct) | segdup | paralogous family |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | GNAI1 | 7 | 210 | DEL | 451 | 4.12x | 1.3e-28 | 3.4e-26 | 1.46 | no | no |
| 2 | NNT | 5 | 286 | COMPLEX | 439 | 3.64x | 1.2e-21 | 1.5e-19 | 1.50 | no | no |
| 3 | MLXIPL | 7 | 211 | COMPLEX | 272 | 4.19x | 1.6e-18 | 1.4e-16 | 0.88 | no | no |
| 4 | KATNAL2 | 18 | 258 | INS | 293 | 4.03x | 2.5e-18 | 1.7e-16 | 0.80 | no | no |
| 5 | POLR2J2 | 7 | 299 | INS | 537 | 3.02x | 3.1e-17 | 1.6e-15 | 3.12 | yes | yes |
| 6 | HILPDA | 7 | 220 | DEL | 632 | 2.59x | 6.0e-14 | 2.7e-12 | 1.37 | no | no |
| 7 | BMPR2 | 2 | 217 | DEL | 281 | 3.06x | 4.3e-10 | 1.6e-08 | 0.54 | no | no |
| 8 | ADK | 10 | 278 | DEL | 444 | 2.57x | 5.0e-10 | 1.7e-08 | 0.02 | no | no |
| 9 | CTAGE6 | 7 | 488 | DEL | 569 | 2.32x | 1.3e-09 | 3.8e-08 | 1.69 | no | yes |
| 10 | ARHGEF35 | 7 | 366 | INS | 585 | 2.29x | 1.6e-09 | 4.2e-08 | 2.34 | yes | yes |
| 11 | SCFD2 | 4 | 213 | DEL | 425 | 2.26x | 4.3e-07 | 1.0e-05 | 1.18 | no | no |
| 12 | ACTA2 | 10 | 224 | DEL | 561 | 2.07x | 5.0e-07 | 1.1e-05 | 2.13 | no | no |
| 13 | SESTD1 | 2 | 216 | DEL | 426 | 2.16x | 2.3e-06 | 4.8e-05 | 1.47 | no | no |
| 14 | GEN1 | 2 | 210 | DEL | 552 | 1.99x | 2.8e-06 | 5.3e-05 | 2.17 | no | no |
| 15 | LPA | 6 | 226 | INS | 431 | 2.04x | 1.5e-05 | 2.6e-04 | 1.44 | yes | no |
| 16 | FBXO22 | 15 | 350 | DEL | 481 | 1.91x | 4.5e-05 | 7.3e-04 | 1.27 | no | no |
| 17 | SPATA31D4 | 9 | 215 | INS | 201 | 2.49x | 5.1e-05 | 7.9e-04 | 0.66 | yes | yes |
| 18 | GOLGA6L4 | 15 | 359 | INS | 388 | 2.01x | 5.9e-05 | 8.6e-04 | 2.42 | yes | yes |
| 19 | RPL24 | 3 | 219 | DEL | 330 | 2.06x | 1.1e-04 | 1.5e-03 | 0.46 | no | yes |
| 20 | PRR20A | 13 | 379 | INS | 784 | 1.63x | 1.7e-04 | 2.2e-03 | 0.01 | no | no |
| 21 | CPNE4 | 3 | 273 | DEL | 1103 | 1.52x | 1.8e-04 | 2.3e-03 | 0.89 | no | no |
| 22 | CLEC18B | 16 | 418 | INS | 901 | 1.58x | 2.1e-04 | 2.6e-03 | 0.94 | yes | no |
| 23 | ACOT1 | 14 | 253 | COMPLEX | 535 | 1.76x | 2.5e-04 | 2.9e-03 | 1.35 | yes | no |
| 24 | DPP10 | 2 | 492 | INS | 797 | 1.61x | 2.6e-04 | 2.9e-03 | 2.24 | no | no |
| 25 | AMY1A | 1 | 596 | INS | 629 | 1.69x | 2.8e-04 | 3.0e-03 | 0.07 | yes | no |
| 26 | GOLGA6D | 15 | 225 | COMPLEX | 229 | 2.18x | 3.6e-04 | 3.6e-03 | 2.42 | yes | yes |
| 27 | RGPD1 | 2 | 210 | INS | 274 | 2.04x | 4.6e-04 | 4.4e-03 | 0.01 | yes | yes |
| 28 | SLC35F4 | 14 | 215 | DEL | 800 | 1.57x | 4.7e-04 | 4.4e-03 | 2.13 | no | no |
| 29 | SPATA5 | 4 | 244 | DEL | 646 | 1.64x | 5.1e-04 | 4.7e-03 | 1.29 | no | no |
| 30 | RABGEF1 | 7 | 348 | INS | 729 | 1.59x | 6.0e-04 | 5.3e-03 | 0.46 | yes | no |
| 31 | CTNNA1 | 5 | 310 | DEL | 432 | 1.76x | 8.8e-04 | 7.5e-03 | 0.02 | no | no |
| 32 | PACSIN1 | 6 | 250 | INS | 496 | 1.69x | 1.0e-03 | 8.4e-03 | 1.03 | no | no |
| 33 | MEI4 | 6 | 211 | DEL | 703 | 1.51x | 3.0e-03 | 2.4e-02 | 0.57 | no | no |
| 34 | TTC28 | 22 | 218 | DEL | 187 | 2.03x | 3.6e-03 | 2.8e-02 | 0.42 | no | no |
| 35 | IQCJ | 3 | 210 | DEL | 465 | 1.59x | 5.1e-03 | 3.8e-02 | 0.82 | no | no |
| 36 | NRG1 | 8 | 211 | DEL | 237 | 1.86x | 5.3e-03 | 3.8e-02 | 0.66 | no | no |
| 37 | ATF6 | 1 | 407 | COMPLEX | 993 | 1.39x | 5.4e-03 | 3.8e-02 | 2.45 | no | no |
| 38 | GRID1 | 10 | 211 | DEL | 516 | 1.55x | 5.7e-03 | 4.0e-02 | 2.66 | no | no |
| 39 | NUTM2G | 9 | 476 | DEL | 816 | 1.42x | 6.5e-03 | 4.4e-02 | 0.10 | no | no |

## How to read the columns

- **enrichment** — observed vs expected count of bottom-5% variants. 1.0x = nothing.
- **p** — Poisson probability of this many by chance.
- **q** — the same p corrected for the 264 loci tested. The line is q < 0.05.
- **most extreme variant (pct)** — the single most suspicious variant in the locus, as a
  percentile; 0.07 means it is in the top 0.07% of the genome for its carrier count.
- **paralogous family** — the gene belongs to a duplicated family where mis-mapping can
  manufacture the signal; treat these as artifact-prone.
- **segdup** — the locus falls inside a segmental duplication (UCSC hg38 track).

## What this list is and is not

**Is:** 39 loci where the SNPs *around* a structurally variable region show more
recently-coalescing common variants than chance. 31 of these are not in an obviously
paralogous family and are the substantive candidates; 7 are paralogous families and should be
treated with suspicion.

**Is not:** evidence that the structural allele itself is selected. In these regions the
copies are collapsed onto one coordinate, so the structural variants are largely invisible to
a SNP-based measure - which is why the class median is 0.82x (below chance) and why this list
is a shortlist for structural analysis, not a conclusion.

Two numbers worth carrying: the class median enrichment across all 264 loci is **0.82x** with
66% below 1.0, so passing loci are a genuine tail; and the top of the list is dominated by
loci with q < 1e-7 (GNAI1, NNT, MLXIPL, KATNAL2, POLR2J2, HILPDA, BMPR2, ADK), which is a
different regime of confidence from the q ~ 0.02-0.04 entries at the bottom.

**Caveat on q itself:** variants within a locus are linked, so treating ~500 of them as
independent Poisson counts overstates certainty. Quote these as screen-level FDRs, not as
genetic p-values.
