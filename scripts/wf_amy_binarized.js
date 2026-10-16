const lanes = [
  {
    key: "amy-binarized",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex. Reproduce the amylase-paper style analysis on our data: collapse the AMY1 structural locus into TWO alleles - allele A = the ancestral/short state, allele B = all longer states - and ask whether the directional signal is visible in that binarised framing, and how a population-structure correction treats it.",
      "",
      "HARD RULES",
      "- Do NOT run any git command; parent commits. Append to docs with clearly-marked sections; never rewrite another author's conclusions.",
      "- Write ONLY under analysis/, docs/, fig/, data/derived/.",
      "- .venv/bin/python. Label claims [observation]/[inference]/[speculation].",
      "- PRE-REGISTER the binarisation, the statistics and the decision rule in docs/amy-binarized-prereg.md BEFORE running; say it was written before results. State that the A/B cut is a MODELING CHOICE and test its sensitivity.",
      "- Stop by ~35 minutes; write partial results incrementally.",
      "",
      "## CONTEXT YOU MUST READ FIRST",
      "- docs/cn-as-snp-results.md and docs/cn-as-snp-prereg.md - the previous counterfactual on the same 225 matched ancients (per-haplotype synthetic loci). Its conclusions: ordering gate passed but absolute |X| not reproduced (~2-4x compressed); direct X-equivalents for single haplotypes were small (H1^a 2.73, H5 -1.79); LCT-anchored extrapolation 11.49/7.56 but SLC45A2-anchored 5.82/3.83 (anchor-dependent); NO GRM absorption; and the ancestry-conditioned null is uninformative because LCT itself gets P=0.965.",
      "- docs/two-channel-results.md - the within/carriage decomposition machinery, already built.",
      "- The reason to binarise: single haplotypes are rare (H5 had dosage 145/73/7 among 225; H9 only 2 alternate alleles), so per-haplotype power is poor. Pooling ALL longer haplotypes into one allele makes B common and restores power.",
      "",
      "## WHAT TO BUILD",
      "1. docs/amy-binarized-prereg.md - specify first: the exact A/B definition, the sensitivity variants, the statistics, the null, the controls and the decision rule.",
      "2. Define A and B from the structural haplotype labels in the Bolognini/Sudmant/Garrison release (extract data/raw/amy/amylase_diversity_project-v1.1.zip to a scratch dir): graph_genotyping/assess_graph_genotypes/supplementary_tables/haplotype_structure_for_supplement.tsv gives per-individual H1/H2, and ancient_metadata_for_supplement.tsv gives age and ancestry component. The haplotype names encode copy number (H1^a = 1 AMY1, H2A0 = 2, H3^r = 3, H4A2 = 4, H5 = 5, H7 = 7, H9 = 9, etc.) - derive the mapping from the names AND check it against measured mean AMY1 per haplotype, and report the check. PRIMARY cut: A = the shortest/ancestral haplotypes, B = everything longer. State the exact cut. SENSITIVITY: at least two alternative cuts (e.g. A = H1^a only; A = all haplotypes with <=2 AMY1 copies) and report whether the conclusion depends on the cut.",
      "3. Per individual, dosage of B = number of their two haplotypes assigned to B (0/1/2). That is the binarised allele. Report dosage counts, B frequency, and B frequency by ancestry component ordered by mean date, WITH counts at every step.",
      "4. Run the same GLMM/counterfactual as the previous lane (reuse analysis/ code; sim/akbari_glmm.py; LOCO genome-wide GRM on the same 225 matched ancients) and report: Z, X-equivalent, and the LCT-anchored and SLC45A2-anchored extrapolations, against 3.16 and 5.45. Report GRM-ON and GRM-OFF for the binarised allele. Report the permutation null.",
      "5. THE QUESTION UNDER TEST - how a structure correction treats it: run the two-channel decomposition (docs/two-channel-results.md) on the binarised allele and report beta_within, beta_between and beta_total. State explicitly whether the binarised signal is (a) mostly within-history, (b) mostly carriage/composition, or (c) both - and whether a population-structure correction removes it. If the correction removes it, say so plainly and show the numbers.",
      "6. POWER / DETECTION FLOOR: with n=225 and the observed B frequency, report the smallest within-history change that would have been detectable, and state whether the observed B trajectory is above or below that floor. This is what determines whether a null result is informative or just underpowered.",
      "7. MODERN REPLICATION at the binarised level: compute B frequency in the modern populations (PEL, Pima, Quechua, Maya, Surui, Han/CHB/CHS/JPT, CEU/TSI/IBS/GBR, YRI/GWD/ESN/Mbuti) and report whether B is elevated in high-starch agricultural populations across independent continents.",
      "8. docs/amy-binarized-results.md - the A/B table, the trajectory by component with counts, the GLMM and its GRM on/off, both extrapolations, the two-channel split, the detection floor, the modern replication, the cut-sensitivity result, and an explicit verdict on: IS THE BINARISED AMYLASE SIGNAL VISIBLE, AND DOES A STRUCTURE CORRECTION REMOVE IT?",
      "",
      "## THE SECOND QUESTION (report separately, do not conflate)",
      "The operator notes that LCT is the strongest signal in the dataset and that 'you do not have to get very far from lactase before you are in the same domain as amylase'. Quantify the gap: compute the genome-wide distribution of |X| from the released statistics (analysis/load_akbari.py), and report where the canonical loci sit (LCT 16.52, FADS 16.15, SLC45A2 13.15, SLC24A5 5.99, AMY window 3.795), what fraction of tested variants exceed each, and how many variants and 250kb blocks lie in the 3.16-5.45 band versus above 5.45. State plainly what the strict threshold keeps and what it discards. This is a descriptive distributional point, NOT a claim that the discarded band is all real; state the false-positive caveat.",
      "",
      "## CARRIED-FORWARD FACTS TO VERIFY, NOT ASSUME",
      "- Mean AMY1 by ancestry component: WHG 5.75, Early farmers 5.45, EHG 5.39, Neolithic farmers 7.15, Steppe pastoralists 8.25, Bronze age 6.63, Iron age to modern 7.25. Within-component slopes are all indistinguishable from zero.",
      "- Ancient haplotype frequencies by mean date: H1^a declines 0.77 -> 0.13 (5 of 6 steps); H5 rises 0.00 -> 0.23.",
      "- Modern: H5 0.40-0.50 in East Asia, 0.46 PEL, 0.71 Pima; H9 0.11-0.17 in PEL/Pima/Quechua vs 0.00-0.03 in East Asia.",
      "- Akbari's AMY window: 2,204 SNPs, max |X| = 3.794956889, zero above 5.45.",
      "- 225/288 labellable ancients match the imputed bundle.",
      "",
      "## Validation",
      "- docs/amy-binarized-prereg.md exists before any result, contains the exact A/B cut and the sensitivity variants, and states that the cut is a modeling choice.",
      "- Dosage counts, B frequency by component, and counts at every step are reported.",
      "- GRM-ON and GRM-OFF are both reported for the binarised allele.",
      "- beta_within/beta_between/beta_total are reported and the correction's effect on the signal is stated explicitly.",
      "- The detection floor is reported BEFORE interpreting the binarised trajectory.",
      "- The |X| distributional question is answered with numbers and its false-positive caveat.",
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
