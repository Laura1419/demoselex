const lanes = [
  {
    key: "elai-long-driver",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex. Build a DRIVER so a LONG ELAI run can be launched as a background job, and VERIFY it end-to-end cheaply. The previous lane had to STOP because the configuration could not fit its own wall cap; the actual cost is now known to be ~61 s per EM step at M=54,370, N=1,000, K=5, so a 200-step fit is ~3.1 hours - feasible as a background job, not inside a 35-minute lane.",
      "",
      "HARD RULES",
      "- No git commands; parent commits. Write ONLY under analysis/, docs/, data/derived/elai_amy_dense/. Never modify data/raw/**.",
      "- DISK: genotypes.txt is ~110 MB for this configuration; keep it in the persistent work dir. Total new output under 2 GB. Report `df -h /` before and after.",
      "- .venv/bin/python. Do not install anything.",
      "- TIME BUDGET ~35 minutes. Your deliverable is a VERIFIED driver plus the exact commands, NOT the long run itself - the parent launches that.",
      "- Label claims [observation]/[inference]/[speculation].",
      "",
      "## READ FIRST",
      "- analysis/elai_amy_dense.py - the module that builds the input and runs one fit. Note `build(work)`, `fit(work, panel, markers, k, seed, wall, calibration=False)`, `metrics(panel, samples, pc, runs)`, the module constants OUT/PREREG/SEED/STEPS, and that `fit` supplies `--pkk-site 1:104121000 --pkk` to the patched binary and streams the pkk output through a FIFO named `panel_k{k}_seed{seed}.pkk.tsv` inside `work/output`.",
      "- docs/elai-amy-dense-prereg.md and docs/elai-amy-dense-results.md - the configuration and the resource-gated STOP that motivates this driver.",
      "- The patched binary: /home/erikg/elai-fork/src/elai (its README-PATCH.md documents --pkk and --pkk-site).",
      "",
      "## WHAT TO WRITE: analysis/elai_amy_long.py",
      "Two modes, so that ONE build can be shared read-only by fits running IN PARALLEL with distinct prefixes:",
      "  1. `--mode build` - build the input into a PERSISTENT work directory (NOT /tmp, NOT deleted): use data/derived/elai_amy_dense/work_long/ with an output/ subdir. Save the panel dict, the marker list, the sample list and the PC panel there so later processes can load them without rebuilding. Report the retained marker count, fitted N, and the sizes of panel.inp / panel.pos / genotypes.txt.",
      "  2. `--mode fit --k K --seed S --steps N --wall SECONDS` - load the built input, then call the module's fit() with the given k/seed and the given wall cap. IMPORTANT: `fit` reads module-level state, so set the steps constant the module uses (or pass it) rather than editing the module; if the module hardcodes STEPS you may set the attribute (e.g. `elai_amy_dense.STEPS = N`) from your driver, and say in the doc exactly how you did it. Each fit must write distinct output names (the prefix already includes k and seed) so parallel fits do not collide.",
      "  3. `--mode gates` - after all fits complete, evaluate the same five gates the previous lane specified (known structure incl. independent-PC silhouette and a singleton-robust composition measure; stability; coverage; virtual-site sanity; exclusion), reporting PASS/FAIL/NOT TESTED per component and a PASS/PARTIAL/FAIL verdict. Reuse the previous lane's gate code in analysis/elai_amy_gate.py / elai_amy_dense_report.py if it fits; do not invent new thresholds - use the preregistered ones.",
      "Make every mode deterministic, log wall-clock, and write a small JSON summary per step.",
      "",
      "## VERIFY END-TO-END CHEAPLY (this is the deliverable)",
      "Run the full pipeline at a TINY step count so it completes in minutes:",
      "  python3 analysis/elai_amy_long.py --mode build",
      "  python3 analysis/elai_amy_long.py --mode fit --k 5 --seed 20261010 --steps 3 --wall 900",
      "Then confirm with evidence: the fit started and exited cleanly; the ELAI stdout shows the expected invocation; a pkk stream was produced and consumed (report the number of markers/individuals it contained and check the arithmetic the patch documents: the ordered joint sums to 1 per locus); the convergence trace was recorded; and the outputs landed in the persistent work dir. Then run `--mode gates` and report that it executes and produces a verdict (which will legitimately be a FAIL/NOT TESTED at 3 steps - that is fine and expected; you are validating the DRIVER, not the science).",
      "Also verify the PARALLELISM works: launch two tiny fits (k=5 and k=3) at the same time from the same built input and show both complete without colliding on files.",
      "",
      "## DELIVERABLES",
      "- analysis/elai_amy_long.py (the driver) in the repo.",
      "- docs/elai-amy-long-run.md: how to run it, the verified tiny-run evidence, the exact THREE LONG-LAUNCH COMMANDS to hand to the parent (K=5 seed 20261010; K=5 seed 20261011; K=3 seed 20261010; all --steps 200, --wall 15000), the expected per-fit wall-clock from the 61 s/step measurement, the total wall-clock if run in parallel, the expected output filenames, and the gate command to run afterwards. State plainly that the long run takes hours and is the parent's to launch.",
      "- Report `df -h /` after the verification and how much disk the built input occupies.",
      "",
      "## Validation",
      "- The tiny build+fit+gates cycle completes and its evidence is pasted into the doc.",
      "- Two simultaneous tiny fits both complete without file collisions.",
      "- The three long-launch commands in the doc are copy-pasteable and reference paths that exist.",
      "- The expected wall-clock is derived from the measured per-step cost, and the derivation is shown.",
      "- Nothing is deleted that the long run needs; the work dir persists."
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
