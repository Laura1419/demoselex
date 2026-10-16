const lanes = [
  {
    key: "elai-amy-dense",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex (analysis) + /home/erikg/elai-fork (PATCHED ELAI, public at github.com/ekg/ELAI, default branch lower-state-posterior). REVISION of a failed run. Read docs/elai-amy-prereg.md and docs/elai-amy-fit-results.md FIRST - the previous attempt is the baseline you must beat and its failure modes are diagnosed.",
      "",
      "HARD RULES",
      "- No git commands; parent commits. Write ONLY under data/derived/, analysis/, docs/, fig/ in demoselex, and a scratch dir you delete. Never modify data/raw/akbari/tier2 or data/raw/aadr-v66.p1.",
      "- DISK CONTESTED (~70G free on a 100%-full shared 14T volume). Stream; never write a file over ~1 GB; total new output under 8 GB; `df -h /` before large writes; delete scratch.",
      "- .venv/bin/python. Do not install anything.",
      "- PRE-REGISTER the revised configuration in docs/elai-amy-dense-prereg.md BEFORE fitting, stating explicitly what changed vs the failed run and WHY, and that the change is motivated by the diagnosed failure modes rather than by the gate results. Label claims [observation]/[inference]/[speculation].",
      "- TIME BUDGET ~35 minutes. Bound the number of ELAI fits (see below). A gated STOP with numbers is a valid outcome; an ungated biological number is not.",
      "",
      "## THE FAILED RUN (baseline to beat)",
      "- 1,164 fitted markers (2-kb thinning gave 9,330 candidates; cohort MAF>=0.05 and call rate>=0.95 then cut to 1,164). The published template (Mies & Mathieson 2026) used 572,223 SNPs.",
      "- Cohort: 128 dated ancients drawn from 5,981 eligible AADR-matched individuals.",
      "- EM NOT CONVERGED: 'the final likelihood still improves after 20 steps (primary last increment 119.681)'.",
      "- Independent-PC panel was only 145 SNPs on chr2 - too sparse for a meaningful silhouette.",
      "- GATES: coverage PASS (128/128), virtual-site PASS, exclusion PASS; known structure FAIL (independent-PC silhouette -0.0578 vs >0.10); stability FAIL (seed ARI 0.0218, seed agreement 7.8%).",
      "- Verdict was STOP; no tuning toward acceptance was attempted, which was correct. Do NOT simply relax thresholds - change the CONFIGURATION.",
      "",
      "## WHAT TO CHANGE (all three, each addressing a diagnosed cause)",
      "1. MARKERS - the main fix. Thin AFTER QC, not before: first filter on MAF and call rate at full density, THEN select spaced markers. Target >=50,000 retained markers in the region (raise the region if needed, e.g. chr1:90,000,000-115,000,000, still excluding the core 103,915,000-104,327,000). Report the pre-thinning count, the filter, the spacing and the retained count. Currently unused: the bundle has 52,382,872 variants; there is no shortage.",
      "2. COHORT - larger but bounded by cost, and MARKERS TAKE PRIORITY. Use as many ancient individuals as fit the time budget (aim 1,000-3,000; report the exact number and how many carry dates). Do not drop to 128 again.",
      "3. EM TO CONVERGENCE - run until the relative log-likelihood improvement falls below a stated tolerance (suggest 1e-4) with a hard cap (suggest 200 steps). REPORT THE CONVERGENCE TRACE and whether it converged; a non-converged fit invalidates the cluster labels, so this is now a GATE, not a detail.",
      "4. INDEPENDENT PC PANEL - use >=5,000 genome-wide markers EXCLUDING chr1 for the structure validation, so the silhouette is meaningful. Report the count and spacing.",
      "",
      "## COST BOUNDING (mandatory - ELAI cost is O(N * M * (S*K)^2))",
      "- Run at most THREE fits: PRIMARY K=5 seed A; then K=5 seed B (stability); then K=3 seed A (cross-K) if time allows. Report wall-clock for each.",
      "- If the first fit is too slow to allow the others, SAY SO, report the primary gates, and state which stability checks did not run. Do not silently skip them.",
      "",
      "## GATES (same five; report NUMBERS and pass/fail per gate)",
      "(a) KNOWN STRUCTURE: independent-PC silhouette (>=0.10), plus composition NMI against published AADR group and country AND a date-ordered check. NOTE the previous run had 111 groups of which 101 were singletons, which makes NMI unreliable - so also report a composition measure robust to singleton groups (e.g. per-country or per-date-bin agreement, or a distance-based measure), and do not let a passing NMI carry a failing silhouette.",
      "(b) STABILITY: seed agreement, cross-K ARI, thinning sensitivity, and the virtual-site argmax agreement. Precedents that failed: ARI 0.022, seed agreement 7.8%.",
      "(c) COVERAGE: usable individuals and how many dated ancients are in the complete --pkk export.",
      "(d) VIRTUAL SITE: ordered joint normalises to 1, q symmetric, site demonstrably integrated.",
      "(e) EXCLUSION: assert no fitted marker lies in the core and the virtual site does.",
      "If a gate fails: STOP, report the numbers, and state what would be needed. No allele-frequency or biological claim either way.",
      "",
      "## DELIVERABLES",
      "- docs/elai-amy-dense-prereg.md (before fitting): revised configuration, the explicit change-vs-failure rationale, the gates, the cost bound, the decision rule.",
      "- docs/elai-amy-dense-results.md: exact commands, marker/cohort/PC counts, the EM convergence trace and verdict, wall-clock per fit, the five gate numbers, the local-state summary IF gates pass (state count, sizes, composition by label/country/date, effective number of independent states), and either the gated STOP or the pass. Do NOT extend to allele-frequency change - that is a separate pre-registered step.",
      "- analysis/ scripts, deterministic, with recorded seeds; data/derived/elai_amy_dense/ for small summaries only.",
      "",
      "## Validation",
      "- The revised configuration is registered before fitting, and each change is tied to a diagnosed cause.",
      "- EM convergence is reported as a gate with its trace, not assumed.",
      "- Every gate reports NUMBERS and an explicit pass/fail.",
      "- Wall-clock per fit and free disk after the run are reported.",
      "- If stopped on a gate, that is the stated outcome with no biological claim."
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
