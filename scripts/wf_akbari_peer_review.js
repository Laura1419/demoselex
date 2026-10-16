const lanes = [
  {
    key: "akbari-peer-review",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex. Evaluate the PEER REVIEW of Akbari et al. 2026 and cross-check it against this project's findings.",
      "",
      "HARD RULES",
      "- Do NOT run any git command; parent commits. Append to docs with clearly-marked sections; never rewrite another author's conclusions.",
      "- Write ONLY under docs/, analysis/, data/derived/.",
      "- .venv/bin/python (pypdf is installed). Label claims [observation]/[inference]/[speculation]. QUOTE the reviewers and authors verbatim (with page/version) for anything you assert they said - no paraphrase presented as quotation.",
      "- Stop by ~30 minutes; incremental writes.",
      "",
      "## SOURCE (already downloaded, immutable, do not modify)",
      "data/raw/akbari/41586_2026_10358_MOESM4_ESM.pdf - the Peer Review File (61 pages), sha256 8c94cddd5cd6c851403f68a4c3168dedfe5f5d29d1a973a808c0a0960365d108. It contains ALL reviewer reports in order by version, followed by ALL author rebuttals in order by version. Extract the text with pypdf and work from the text, not from the abstract or from memory. If a page extracts badly, say so rather than guessing its content.",
      "",
      "## WHAT TO PRODUCE: docs/akbari-peer-review-evaluation.md",
      "1. INVENTORY: how many versions, how many referees per version, and the report/rebuttal structure. Note any redactions.",
      "2. CONCERN TABLE: one row per substantive reviewer concern, with columns: version, referee, the concern in a short verbatim quote, the authors' response (short verbatim quote or 'no direct response'), whether it was addressed or deferred, and your assessment of whether the response actually resolves it. Be willing to mark a response as non-responsive or as resolving only part of the concern.",
      "3. DEEP DIVE on concerns bearing on THIS PROJECT'S subject. Specifically search the reports and rebuttals for, and report verbatim what was said about:",
      "   (a) population structure / admixture / the GRM / whether the correction could remove real selection, and whether a signal could be an artefact of population turnover or migration;",
      "   (b) the choice and justification of the |X|>5.45 significance threshold, the FDR calibration, and any reviewer pushback on the threshold or on the non-normal/null behaviour of X;",
      "   (c) the AMYLASE locus (AMY1/AMY2, copy number, the Bolognini/Sudmant/Garrison work) - was it raised by a reviewer or by the authors, and how was it treated?",
      "   (d) any demand for haplotype-based or LD-based corroboration, replication across populations, or independent lines of evidence;",
      "   (e) anything about selection correlated with demography, parallel/repeated selection, or the identifiability of selection versus drift.",
      "   If a topic does not appear, say explicitly that it does not appear rather than inferring that it was considered.",
      "4. CROSS-CHECK against our results. Our established findings are:",
      "   - Akbari's AMY-region maximum is |X| = 3.794956889 over 2,204 SNPs, with ZERO above 5.45, while LCT 16.52, FADS 16.15, SLC45A2 13.15, SLC24A5 5.99 all pass. So the canonical amylase locus falls in the FDR-50%-to-pi>99% band that the headline threshold discards.",
      "   - The band 3.16<=|X|<=5.45 contains 356,335 variants (3.66% of tested) versus 31,939 above 5.45 (0.33%); at 250 kb block level, 723 blocks carry a >5.45 variant while 9,285 carry band-only evidence.",
      "   - Our within-lineage versus carriage decomposition gives AMY ~94% carriage, with the within-history slope below the detection floor at n=225; the binarised short/long allele clears FDR-50% (X-equiv 3.64) and survives the structure correction (~97% GRM retention).",
      "   - We found NO low-|X|/high-replication class overall: low-|X| loci show high replication at 4.3% versus a 5% matched-null rate.",
      "   - We found the two blind spots overlap: Akbari is blind to the ancestry-collinear direction, and a within-history statistic is blind to carried-in change; both are the same regime.",
      "   - Mies & Mathieson (2026) report their own FADS1/SLC24A5 ancestry deviations 'support either within-background change or underestimated ancestry shifts'.",
      "   For each, state whether the peer review ALREADY raised, anticipated, dismissed or ignored it, with verbatim support. This is the point of the exercise.",
      "5. WHAT THE REVIEW DID NOT SETTLE: list the questions a reviewer asked that the rebuttal left open, and any concern where the authors' response is weaker than the concern. Be specific and quote.",
      "6. HONEST ASSESSMENT: in a short closing section, state what the peer review changes about our project's framing. In particular: did the referees raise the over-correction or threshold concern at all? If they did and the authors answered it, summarise the answer; if they did not, say so plainly, because that is itself informative about the state of the field.",
      "",
      "## CONSTRAINTS",
      "- Do not manufacture a connection between a reviewer comment and our results; if the review simply does not discuss it, say so.",
      "- Distinguish what the AUTHORS asserted in rebuttal from what they DEMONSTRATED. A rebuttal claim is not evidence.",
      "- Reported numbers from the rebuttals should be quoted as the authors' statements, not adopted as facts.",
      "- Note explicitly if the file structure means some rebuttals respond to earlier versions only, so a concern raised late may have a thin response.",
      "",
      "## Validation",
      "- docs/akbari-peer-review-evaluation.md exists with the inventory, the concern table and the deep dive.",
      "- Every quoted concern and response is traceable to a page and a version in the PDF.",
      "- The absent-topic statement is explicit where a topic is not discussed.",
      "- The cross-check section names, for each of our findings listed above, whether the review raised it, with support or an explicit 'not raised'.",
      "- Any disagreement with a number quoted here is stated with the disagreeing number."
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
