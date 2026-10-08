# I. Akbari vs trajectory replication - first pass (Tier-1)

Status: **scouting only, and the pre-registered hypothesis is NOT supported at
this power.** ~26% of Akbari's individuals, 1.23M of 9.7M sites, partly-stable
trajectories, leave-one-out matched null for the low-|X| pool.

- loci scored: **19809**
- quadrants: {'Q4': 14229, 'Q2': 4488, 'Q3': 635, 'Q1': 457}
- low |X| < 3.16: fraction with empirical P <= 0.05 = **0.0427** (null expectation 0.05)
- high |X| >= 3.16: same fraction = **0.0924** (~2x the null rate)
- Spearman corr(|X|, replication score) = **0.086**

## Verdict for Decision 5 (first pass)

**No low-|X|/high-replication excess.** Low-Akbari loci show high replication at
a rate *at or below* the matched-null expectation, while Akbari-significant loci
are about twice as likely to show it. The two statistics therefore agree rather
than dissociate, and this first pass finds no evidence that the Akbari correction
suppresses a class of replicated signals. Per docs/MISSION.md this is the
'over-correction appears limited under this test' branch for Tier-1 data.

### Why this is not yet final

1. Tier-1: 26% of individuals, 1.23M of 9.7M sites, pseudo-haploid low-coverage
   genotypes - both statistics are attenuated.
2. Trajectories are partly stable and *over-split* (three Avar-period clusters at
   r = 0.6-0.8, and T0/T4 at r = 0.88), so the real N_eff composition is
   questionable even though N_eff ~ 4.7 is computed from a valid PD matrix.
3. The low-|X| empirical null is leave-one-out within the same pool (a scoping
   approximation), and the high-|X| null is a finite-control caliper matcher.
4. Structural reason to expect limited gain in the worst case: selection that is
   *carried in* by an expanding population produces flat within-trajectory slopes,
   so the trajectory method is conservative there too - by design, not by bug.

## Q3 rows (low |X|, empirical P <= 0.05) - null-rate, NOT a discovery list

| locus | chrom:pos | X | z | C | n_traj | N_eff | P |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2:135750794:G:A | 2:135750794 | -0.08 | -31.68 | 0.58 | 12 | 4.56 | 0.0081 |
| 4:146686561:C:A | 4:146686561 | -0.41 | 19.61 | 0.83 | 12 | 4.56 | 0.0185 |
| 14:84242268:G:A | 14:84242268 | -1.62 | 21.74 | 0.75 | 12 | 4.56 | 0.0278 |
| 2:88462352:A:G | 2:88462352 | 1.50 | -20.13 | 0.75 | 12 | 4.56 | 0.0074 |
| 10:93452141:C:T | 10:93452141 | -0.86 | 19.98 | 0.75 | 12 | 4.56 | 0.0149 |
| 6:170563478:T:C | 6:170563478 | -0.89 | 19.66 | 0.75 | 12 | 4.56 | 0.0072 |
| 16:66203027:C:A | 16:66203027 | -2.48 | 21.99 | 0.67 | 12 | 4.56 | 0.0189 |
| 2:694759:G:A | 2:694759 | -1.10 | 21.91 | 0.67 | 12 | 4.56 | 0.0081 |
| 11:2149864:T:C | 11:2149864 | -2.52 | -18.21 | 0.75 | 12 | 4.56 | 0.0085 |
| 1:205063949:C:T | 1:205063949 | -1.69 | -16.19 | 0.83 | 12 | 4.56 | 0.0133 |
| 6:32902001:G:A | 6:32902001 | 0.18 | -14.66 | 0.92 | 12 | 4.56 | 0.0091 |
| 2:40582014:T:G | 2:40582014 | -2.04 | -19.97 | 0.67 | 12 | 4.56 | 0.0098 |
| 6:31341379:C:T | 6:31341379 | -2.44 | -22.78 | 0.58 | 12 | 4.56 | 0.0169 |
| 3:133842598:A:G | 3:133842598 | 0.96 | 17.64 | 0.75 | 12 | 4.56 | 0.0052 |
| 1:247360368:T:C | 1:247360368 | -1.04 | 17.53 | 0.75 | 12 | 4.56 | 0.0118 |
| 1:177928961:C:A | 1:177928961 | 0.68 | -19.60 | 0.67 | 12 | 4.56 | 0.0043 |
| 17:7986716:T:C | 17:7986716 | 0.11 | 22.37 | 0.58 | 12 | 4.56 | 0.0156 |
| 11:10398518:C:A | 11:10398518 | -0.44 | 19.34 | 0.67 | 12 | 4.56 | 0.0139 |
| 4:87863770:T:C | 4:87863770 | 2.99 | 21.88 | 0.58 | 12 | 4.56 | 0.0061 |
| 8:134257450:A:G | 8:134257450 | -0.75 | 16.98 | 0.75 | 12 | 4.56 | 0.0435 |
| 3:102609928:A:G | 3:102609928 | -1.54 | 15.49 | 0.82 | 12 | 4.56 | 0.0097 |
| 1:157387799:C:T | 1:157387799 | 2.13 | 16.72 | 0.75 | 12 | 4.56 | 0.0102 |
| 13:44931134:A:C | 13:44931134 | 1.33 | 16.48 | 0.75 | 12 | 4.56 | 0.0286 |
| 2:19321500:G:A | 2:19321500 | -1.11 | -21.16 | 0.58 | 12 | 4.56 | 0.0177 |
| 13:101364651:A:C | 13:101364651 | 0.62 | 18.47 | 0.67 | 12 | 4.56 | 0.0286 |
| 5:113679204:T:C | 5:113679204 | -1.60 | 18.30 | 0.67 | 12 | 4.56 | 0.0132 |
| 4:102605800:G:A | 4:102605800 | 1.11 | -20.88 | 0.58 | 12 | 4.56 | 0.0125 |
| 15:59598381:A:G | 15:59598381 | 1.87 | 18.25 | 0.67 | 12 | 4.56 | 0.0263 |
| 6:29722515:A:G | 6:29722515 | -2.32 | -22.29 | 0.55 | 12 | 4.56 | 0.0114 |
| 19:2396243:C:T | 19:2396243 | 1.14 | -16.15 | 0.75 | 12 | 4.56 | 0.0200 |
| 8:2107237:A:G | 8:2107237 | 0.71 | 18.08 | 0.67 | 12 | 4.56 | 0.0080 |
| 11:27745764:C:T | 11:27745764 | 0.13 | 13.14 | 0.92 | 12 | 4.56 | 0.0066 |
| 17:63650772:C:T | 17:63650772 | -1.12 | 20.63 | 0.58 | 12 | 4.56 | 0.0417 |
| 11:2234690:A:G | 11:2234690 | -0.15 | 18.03 | 0.67 | 12 | 4.56 | 0.0169 |
| 11:128816603:G:T | 11:128816603 | 0.20 | 17.90 | 0.67 | 12 | 4.56 | 0.0167 |
| 6:17013648:A:G | 6:17013648 | 1.90 | 20.46 | 0.58 | 12 | 4.56 | 0.0339 |
| 12:104669483:T:C | 12:104669483 | 2.06 | -17.82 | 0.67 | 12 | 4.56 | 0.0070 |
| 12:28186208:T:C | 12:28186208 | 1.96 | 17.75 | 0.67 | 12 | 4.56 | 0.0125 |
| 5:164874677:C:T | 5:164874677 | 0.16 | -14.16 | 0.83 | 12 | 4.56 | 0.0141 |
| 9:9489220:T:C | 9:9489220 | 1.72 | 17.69 | 0.67 | 12 | 4.56 | 0.0102 |

## Q1 (both high) - consensus loci

| locus | chrom:pos | X | z | C | P |
| --- | --- | ---: | ---: | ---: | ---: |
| 5:33962770:G:A | 5:33962770 | 10.34 | 32.88 | 0.83 | 0.0227 |
| 2:135539967:C:T | 2:135539967 | -7.80 | -30.82 | 0.83 | 0.0104 |
| 4:38880384:A:G | 4:38880384 | 9.02 | 29.47 | 0.73 | 0.0085 |
| 6:32583299:C:T | 6:32583299 | -4.35 | 22.36 | 0.90 | 0.0116 |
| 2:136590558:C:T | 2:136590558 | -6.87 | -29.57 | 0.67 | 0.0079 |
| 6:31394700:A:G | 6:31394700 | 12.41 | -25.62 | 0.75 | 0.0120 |
| 6:31142297:C:A | 6:31142297 | -5.97 | -24.26 | 0.75 | 0.0127 |
| 6:29692562:A:G | 6:29692562 | -4.71 | -27.72 | 0.64 | 0.0435 |
| 2:238918061:C:T | 2:238918061 | -3.45 | -26.10 | 0.67 | 0.0167 |
| 6:32670576:A:G | 6:32670576 | -3.18 | -23.35 | 0.73 | 0.0083 |
| 6:32066177:T:C | 6:32066177 | -7.88 | -24.19 | 0.70 | 0.0109 |
| 6:32337686:A:G | 6:32337686 | -6.93 | -28.63 | 0.58 | 0.0123 |
| 2:135939824:A:G | 2:135939824 | -6.82 | -19.64 | 0.83 | 0.0052 |
| 10:50481447:G:A | 10:50481447 | 4.30 | -17.34 | 0.92 | 0.0065 |
| 6:31377587:G:A | 6:31377587 | -3.83 | -27.19 | 0.58 | 0.0067 |

## Interpretation rules (do not skip)

- Counting Q3 at an empirical P threshold *by construction* yields ~5% of the
  low-|X| pool; only a departure from that rate is evidence. Here the rate is
  4.3%, i.e. no departure.
- Do not read the Q3 table above as candidate loci. It is the null-rate tail.
- Gene/functional annotation is NOT attached; positions are raw.
