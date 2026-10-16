const lanes = [
  {
    key: "elai-amy-fit",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex (analysis repo) and /home/erikg/elai-fork (the PATCHED ELAI, public at github.com/ekg/ELAI, default branch lower-state-posterior). Apply the patched ELAI to real ancient data and GATE the result. Do NOT interpret biology unless the gates pass - that rule is the whole point of this lane.",
      "",
      "HARD RULES",
      "- Do NOT run any git command in either repo; the parent commits. Read-only git inspection (git log/show/diff) is permitted.",
      "- Write ONLY under /home/erikg/demoselex/{analysis,docs,data/derived,fig}, /home/erikg/demoselex/data/raw/derived_elai_input/, and a scratch dir you delete at the end. Never modify data/raw/akbari/tier2 or data/raw/aadr-v66.p1.",
      "- DISK IS CONTESTED (~55G free on a 100%-full 14T volume shared with unrelated projects). Stream; never write a matrix over ~500 MB; keep total new output under 5 GB; `df -h /` before large writes and report free space. Delete scratch when done.",
      "- .venv/bin/python. Do not install anything (GSL is present; the patched ELAI already builds).",
      "- PRE-REGISTER the run configuration, the gates and the decision rule in docs/elai-amy-prereg.md BEFORE fitting. Label claims [observation]/[inference]/[speculation].",
      "- Stop by ~35 minutes. A gated STOP with honest numbers is a successful outcome; an ungated biological number is not.",
      "",
      "## THE INSTRUMENT (already validated - use it, do not rebuild it)",
      "/home/erikg/elai-fork/src/elai is a patched ELAI (187 insertions, 0 deletions) that adds:",
      "  --pkk            ordered joint lower-state diploid posterior p(k1,k2|data) marginalised over the upper S states, plus unordered q(a,b)",
      "  --pkk-site chrom:position    evaluate it at a virtual MISSING-EMISSION site inside a same-chromosome gap (for an excluded locus)",
      "Validated against exact brute-force enumeration (1.4e-15 ordered; 1.1e-14 virtual), internal m_pk/2 marginal consistency 1.1e-16, and 12/12 byte-identical existing outputs. Read README-PATCH.md in that repo for the exact CLI, output format and normalisation conventions.",
      "IMPORTANT upstream convention: `--pk1` writes expected diploid COPY COUNTS summing to 2 (p(k) = --pk1[k]/2); `--ps2` is normalised (sums to 1). Always know which you are reading.",
      "Also READ the shipped manual first: /home/erikg/elai-fork/elai-manual.tex (and elai-manual.pdf) and the /home/erikg/elai-fork/example directory, and derive the correct invocation and input format from THOSE, not from guesswork. Report the exact command you use.",
      "",
      "## STEP 1 - BUILD THE ELAI INPUT from the real data",
      "- Genotypes: /home/erikg/demoselex/data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.{pgen,psam,pvar} (52,382,872 variants x 16,339 samples, imputed). Loader /home/erikg/demoselex/analysis/load_akbari_geno.py (returns ALT dosage; uint32 indices for sample-subset and batch reads; read the file first). The pvar is 1.7 GB - STREAM it to select variants by position; do not load it whole.",
      "- Region: chr1 spanning the AMY locus, e.g. chr1:95,000,000-115,000,000. EXCLUDE the AMY core chr1:103,915,000-104,327,000 from the FITTED markers (this is mandatory locus-independence: AMY1 copies run to 104,301,314 in hg19 refGene). The virtual site for --pkk-site must lie INSIDE the excluded core - use 104,121,000 unless you have a documented reason to change it, and state the choice.",
      "- Thin markers to a tractable panel within the region (state the count; aim for a few tens of thousands, well spaced). Report the marker filter (MAF, call rate) and the resulting count. Do not hand-pick ancestry-informative markers: the design is reference-free and the linked model does the work.",
      "- Cohort: state the sample choice and count. Prefer the ANCIENT panel individuals that have dates (they are the ones we can order in time), but ELAI's EM needs enough samples to fit clusters - so report both the fitted cohort size and, separately, how many of them are dated ancients that can enter the temporal analysis later. A bounded cohort is acceptable; say so.",
      "- Write the input in exactly the format the manual specifies (BIMBAM or its own format - verify against the example).",
      "",
      "## STEP 2 - FIT AND GATE (no biology before this passes)",
      "- Fit with S=1 (the linked-model limit; reference-free in the sense of requiring no named ancestral training sources) and K=5, with K=3 and K=8 as pre-registered sensitivities. Use enough EM steps and a recorded seed. Record the exact command for each.",
      "- Run --pkk for the individuals and --pkk-site chr1:104121000.",
      "- GATE (a) KNOWN STRUCTURE: do the fitted lower/upper clusters reproduce known population structure? Validate against independent AADR published group labels, country and date (from /home/erikg/demoselex/data/raw/aadr-v66.p1/*.anno and .ind) and against PC position. Report silhouette and NMI as NUMBERS. The precedent to BEAT is the earlier in-house painter: silhouette < 0.10, group NMI < 0.10, country NMI < 0.05 - which FAILED.",
      "- GATE (b) STABILITY: cluster assignment and state composition stable across K (3/5/8), across at least two seeds, and across the marker thinning. Report as numbers (e.g. ARI). Precedent that FAILED: L1 to 0.78, target-local agreement 43%.",
      "- GATE (c) COVERAGE: how many individuals have usable calls on the fitted markers, and how many dated ancients are in the `--pkk` output. Report counts.",
      "- GATE (d) VIRTUAL SITE SANITY: the ordered joint at the virtual site normalises to 1; q is symmetric; and the site is demonstrably integrated (its output changes when the gap is bridged vs not) - the patch's own tests show this, but re-confirm on THIS input.",
      "- GATE (e) EXCLUSION ASSERTION: assert programmatically that no fitted marker lies inside the excluded AMY core and that the virtual site does. Report the check.",
      "If ANY gate fails: STOP. Report the gate numbers, the failing threshold, and what would be needed. Do not proceed to allele frequencies.",
      "",
      "## DELIVERABLES",
      "- docs/elai-amy-prereg.md (written before fitting): the region, the exclusion, the virtual site, the cohort, the thinning rule, the ELAI configuration and its justification from the manual, the five gates and their thresholds, and the decision rule.",
      "- docs/elai-amy-fit-results.md: the exact commands, the input sizes, the gate numbers, and the verdict. If gates pass, state the local-state summary (state count, sizes, composition by label/country/date, and the effective number of independent states) and STOP short of allele-frequency claims - the temporal analysis is a separate pre-registered step.",
      "- analysis/ scripts to build the input and to run/validate the fit, deterministic with recorded seeds.",
      "- data/derived/elai_amy/ : the input panel description, the gate metrics and any --pkk output (small summaries only, not full posterior tensors).",
      "",
      "## Validation",
      "- Every gate reports NUMBERS, and its pass/fail is explicit.",
      "- The exact ELAI invocation is stated and reproducible from the manual's documented options.",
      "- The exclusion and virtual-site assertions are demonstrated, not asserted.",
      "- Free disk space is reported after the run.",
      "- If the run is stopped on a gate, that is stated as the outcome, with the numbers, and no biological claim is made."
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
