const lanes = [
  {
    key: "elai-patch",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/elai-fork - a fork of haplotype/ELAI at pinned commit 4abc4a8fcb824d20db2dde672cfaf5c92fd09853 (verified: src/diploid.{h,cpp} and src/control.cpp sha256 match the pinned revision). Branch: lower-state-posterior. Goal: PATCH ELAI so it exports the JOINT LOWER-STATE diploid pair posterior, which upstream does not provide.",
      "",
      "HARD RULES",
      "- Do NOT run `git commit`, `git push`, `git checkout -b`, or any history-changing git command. The parent commits and publishes. You MAY run `git diff`, `git status`, `git log` read-only.",
      "- Write ONLY under /home/erikg/elai-fork/ (src/, tests/, docs-new/, README-PATCH.md) and /home/erikg/demoselex/docs/ + /home/erikg/demoselex/analysis/ for the write-up.",
      "- The project is GPL v3. Do NOT modify or remove LICENSE. Any new files must be GPL v3 compatible and say so.",
      "- Do NOT break existing behaviour: every existing CLI option must behave identically. The patch must be MINIMAL and reviewable.",
      "- Build with `make` in src/ (GSL 2.7.1 is installed; the Makefile now links -lgsl -lgslcblas and builds cleanly).",
      "- Stop by ~35 minutes; leave a working build.",
      "",
      "## WHY (the scientific requirement)",
      "We need allele-frequency change conditioned on a LOCAL HAPLOTYPE-SHARING STATE - the fine equivalence class in ELAI's two-layer model, i.e. the LOWER (K) clusters - shared across populations. Upstream cannot export this:",
      "- `--ps2` writes the SxS UPPER-state joint posterior p(s1,s2|data). In src/control.cpp the commented help says so, and `DipInd::compute_ps(nS,nK,pMP,comp_mode)` allocates nS*nS.",
      "- `--pk1` writes only the MARGINAL lower posterior p(k|data) into `m_pk` (nLoci x nK), via `DipInd::compute_pk(nS,nK,pMP,comp_mode)` -> `CalcAll(..., comp_mode, 11)`.",
      "- `--pk2` has phase/Viterbi wording and is not the joint pair posterior.",
      "So the JOINT lower pair p(k1,k2|data) is computed internally during the diploid forward-backward and then discarded. We need it.",
      "",
      "## WHAT TO IMPLEMENT",
      "1. Export the JOINT ORDERED lower-state pair posterior per individual and per locus: p(k1,k2 | data), marginalised over the upper states S. Note the internal diploid recursion is over (S,K) pairs - the marginalisation over S is the whole point, and it must be done on the posterior, not by combining marginals.",
      "2. Also emit the UNORDERED pair q(a,b) (a<=b) so a symmetric sharing state can be used directly.",
      "3. A NEW CLI option (do not overload existing ones; e.g. `--pkk` / `--pkk2`), documented in `--help` text and in README-PATCH.md, writing one file per individual (or a documented table) with locus, k1, k2, p_ordered, q_unordered. State the file naming and format precisely.",
      "4. VIRTUAL MISSING-EMISSION SITE: support requesting the posterior at a specified GENOMIC POSITION (e.g. `--pkk-site <chrom:pos>` or a coordinate argument) that may be a MISSING/absent marker - i.e. a virtual site with no observed genotype, bridged across an excluded gap by the linked transition model. This is required because our tested locus is EXCLUDED from the panel and the gap must be crossed. Implement it as a genuine missing-emission step in the same recursion, NOT as interpolation between flanking outputs, and NOT by multiplying marginals. Document the semantics and the assumption.",
      "5. Keep `m_pk` (the marginal) and all existing outputs byte-identical in behaviour. The new joint posterior's MARGINAL over one of the two pair indices must equal the existing marginal `m_pk` up to numerical tolerance - implement that as a TEST, because it is the strongest available internal consistency check.",
      "",
      "## TESTS (required - this is what makes it solid)",
      "Write a tests/ directory with a runnable script (e.g. tests/run_tests.sh + a python checker) covering at least:",
      " (a) EXACT ENUMERATION on a tiny case: a hand-built cohort of 2-3 diploid individuals and 2-3 linked markers with small S and K, where the full joint posterior can be computed by brute force in python (enumerate all (k1,k2,S) state paths under the SAME transition and emission model), and compare against the patched binary's output to numerical tolerance. State the model equations you enumerated in the doc; if the enumeration cannot match exactly, report the discrepancy honestly rather than loosening the tolerance to force a pass.",
      " (b) NORMALISATION: the ordered joint over (k1,k2) sums to 1 per locus per individual.",
      " (c) SYMMETRY: q(a,b) == q(b,a).",
      " (d) MARGINAL CONSISTENCY: sum_{k2} p(k1,k2) == m_pk[k1] from the existing `--pk1` output, within tolerance.",
      " (e) VIRTUAL SITE: with a marker present vs absent (masked) at the tested position, show the bridging behaves as documented, and show the output changes when the gap is crossed vs not (i.e. the virtual site is actually integrated, not ignored).",
      " (f) REGRESSION: existing `--ps2`/`--pk1` outputs are unchanged by the patch on the same input.",
      "Use the shipped example data if usable, otherwise construct a small BIMBAM-format input; document exactly how the input was built.",
      "",
      "## DELIVERABLES",
      "- The patched source (minimal diff) in src/.",
      "- README-PATCH.md at the repo root: what the patch adds, the exact CLI, the output format, the mathematical definition of the exported quantity, the virtual-site semantics and its assumptions, build instructions (system GSL), and the test instructions. It must ALSO state clearly what the patch does NOT do (e.g. it does not infer phase beyond what ELAI's model already integrates; it is not validated on real ancient data here; the fine-state assignment still inherits ELAI's model assumptions).",
      "- tests/ with the checks above and a single entry point that returns nonzero on failure.",
      "- /home/erikg/demoselex/docs/elai-patch.md: a short appended section recording the pin (repo+commit+sha256 of the three files), what was patched, the test results, and the honest limitations. Do not rewrite other authors' text.",
      "- A working `src/elai` binary and the exact build command.",
      "",
      "## Validation (report all as numbers)",
      "- `cd src && make` succeeds from clean; report the binary path and size.",
      "- `bash tests/run_tests.sh` exits 0 and prints each check's result; report the enumerated-vs-binary maximum absolute discrepancy, the normalisation error, the symmetry error and the marginal-consistency error.",
      "- Existing options produce unchanged output vs the unpatched binary on the same input (state the comparison method).",
      "- `git diff --stat` (read-only) summarises the patch size; report it.",
      "- If any test cannot pass, say so plainly and do not claim it does."
    ].join("\n")
  }
];

const results = await runs.all(lanes);
return results.map((r, i) => ({
  key: lanes[i].key,
  status: r && r.status,
  runId: r && (r.runId || r.id),
  output: r && r.output
}));
