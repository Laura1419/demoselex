<!-- SPDX-License-Identifier: GPL-3.0-or-later -->

## 2026-10-09 — joint lower-state ELAI export

Upstream: https://github.com/haplotype/ELAI, pinned revision
`4abc4a8fcb824d20db2dde672cfaf5c92fd09853`; working fork:
`../elai-fork`, branch `lower-state-posterior`. Before patching, SHA256:

| Pinned file | SHA256 |
|---|---|
| src/diploid.h | a0e7d55bc6927311fe5e97ba1769843f0133036990b3564f3831c92d7ad3e4df |
| src/diploid.cpp | e290f442c6065056c64aef945a306a65f2ce3280caabf6e323c3e8a2a809645d |
| src/control.cpp | 1856cdf5fc21cbde2a54373b41101f308f3f0b807026400670184b2e55f1d286 |

Added `--pkk` (all markers) and `--pkk-site integerChrom:position` (one site).
The ordered K*K joint is the S-marginalised exact diploid forward/backward
posterior, **not a product of marginals**. The unordered probability is emitted
symmetrically, with unique states obtained by filtering k1<=k2. Tables are
`output/<prefix>.pkk.tsv` and `.pkk-site.tsv`; README-PATCH.md defines fields.

An absent site inserts an export-only QQ/unit-emission step. Fitted rs/rk
survival is split by physical-distance fraction, with right-flanking beta on
both steps. No EM refit, genotype leakage, posterior interpolation, or mutation
of existing outputs. Existing-site requests retain observed emission. Multiple
retained records at one coordinate are rejected with count, all marker
indices/IDs/alleles, and pointers to existing SNP-info/log files.
The full exact recursion spans computational block boundaries; memory remains
O(M*(S*K)^2). The export snapshots admixture immediately before `--pk1` because
upstream posterior-output calls themselves update admixture.

**Important upstream scale correction:** `--pk1` writes expected diploid copy
counts, not probabilities. The correct internal check is
`sum_k2 p(k1,k2) = m_pk[k1]/2`. Measured sums over K were 2.000 for all 9
individual/locus cases. `--ps2` is a normalised upper pair posterior; its rounded
sums ranged 0.999–1.001, with no factor of two.

### Validation (measured, not biological validation)

`cd ../elai-fork/src && make clean && make` exited 0 with system GSL **2.7.1**;
`../elai-fork/src/elai` is **527256 bytes**. Existing upstream compiler warnings
remain. `bash tests/run_tests.sh` exited **0**: **15 numerical error gates,
4 posterior-movement gates, 5 coordinate-rejection cases**, and **18/18
byte-identical legacy-file comparisons** against a freshly built pinned binary
with the same system-GSL Makefile. Regression includes frozen parameters and
1 warm + 2 exact EM steps.

The independent Python oracle enumerates every `(s1,k1,s2,k2)` state path for
3 diploid individuals, S=K=2, 3 linked markers (16^3 paths), and an inserted
missing site (16^4 paths). All transition/emission equations, exact input
construction, output-induced admixture updates and tolerances are documented
in the patch README and tests. Maximum absolute discrepancies:

| Check | Value |
|---|---:|
| Joint enumeration vs binary | 1.4432899320127035e-15 |
| Virtual enumeration vs binary | 1.099120794378905e-14 |
| Explicit masked marker enumeration vs binary | 1.9317880628477724e-14 |
| Ordered normalisation | 4.440892098500626e-16 |
| Unordered symmetry | 0 |
| Primary internal marginal vs m_pk/2 | 1.1102230246251565e-16 |
| Rounded CLI marginal vs --pk1/2 | 0.00020933910410692036 |
| Composed split vs original transition | 2.220446049250313e-16 |
| Original marker posterior preservation | 8.326672684688674e-15 |
| 1001-marker computational-boundary preservation | 1.1102230246251565e-16 |

Internal tolerance is 1e-12. CLI tolerance is derived, not tuned:
three-decimal copy-count rounding gives 0.0005/2 = 0.00025 plus arithmetic
allowance. Tracked source diff: **7 files, 187 insertions, 0 deletions**;
`git diff --check` passes. README/tests are additional untracked deliverables.
No files were staged; no commit, push or history modification was performed.

### Honest limitations

The bridge assumes homogeneous within-gap transition behaviour; its biological
accuracy is unvalidated. Left-flanking rather than right-flanking beta moved
the fixed-snapshot virtual posterior by **0.054451195697281352**, so the bridge
has genuine model sensitivity. Linked versus unlinked CLI gap posteriors moved
by **0.46360006631119072**, confirming the transition model is used. At position
150, observed versus masked emission changed the posterior by
**0.11145189187526616**. An explicit missing retained marker versus an absent
virtual marker differs by **2.9181859871885152e-6** at CLI level because upstream
admixture averaging changes with marker count; this is reported, not hidden.
Matched-snapshot transition composition and original-posterior preservation
pass at the errors above.

No additional phase inference, aDNA likelihood support, real ancient-data
validation, genome-build conversion or large-panel scalability validation.
Fine-state labels and inference inherit ELAI's fitted-model assumptions.
New additions are GPL-3.0-or-later compatible; upstream LICENSE is unchanged.
A reviewer must still perform independent acceptance before publication.
