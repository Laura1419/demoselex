const lanes = [
  {
    key: "local-ancestry",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex. Build LOCAL-ANCESTRY-STRATIFIED frequency analysis for ancient DNA. The central problem: global ancestry labels ('WHG', 'steppe', 'Neolithic farmers') are categorical snapshots of a continuum - there was continuous clinal gene flow, not discrete replacement - so a within-label slope is measured inside a moving target and the within/between split depends on where tissue boundaries were drawn. Local ancestry replaces that: the unit becomes a HAPLOTYPE SEGMENT of known ancestry, and the test becomes whether an allele's frequency changes over time ON A FIXED ANCESTRY BACKGROUND. A frequency step on a fixed ancestry background is not admixture.",
      "",
      "HARD RULES",
      "- Do NOT run any git command; parent commits. Append to docs with clearly-marked sections; never rewrite another author's conclusions.",
      "- Write ONLY under analysis/, docs/, fig/, data/derived/, data/raw/annotation/.",
      "- .venv/bin/python (numpy pandas scipy sklearn matplotlib pgenlib). Do not install heavy external tools (no RFMix/Beagle); implement a documented in-house approximation.",
      "- Label claims [observation]/[inference]/[speculation]. PRE-REGISTER the method and decision rule in docs/local-ancestry-prereg.md BEFORE results.",
      "- Stop by ~35 minutes; write partial results incrementally. This is STAGE 1 (machinery + validation + AMY), not the genome-wide application.",
      "",
      "## DATA",
      "- AADR v66.p1 1240K and HO EIGENSTRAT: data/raw/aadr-v66.p1/ with .ind (id, sex, published group label) and .anno (dates, coords). Reader analysis/load_aadr.py - NOTE it now returns ALT dosage (high-bit-first plus allele-inversion convention); verify with analysis/test_aadr_orientation.py.",
      "- Akbari imputed pgen bundle: data/raw/akbari/tier2/ (52,382,872 variants x 16,339 samples; ancient + EUR 1KG). Loader analysis/load_akbari_geno.py (returns ALT dosage; uint32 indices needed for sample-subset and batch reads).",
      "- AMY structural haplotypes: data/raw/amy/ - extract data/raw/amy/amylase_diversity_project-v1.1.zip to a scratch dir. graph_genotyping/assess_graph_genotypes/supplementary_tables/haplotype_structure_for_supplement.tsv has per-individual H1/H2; ancient_metadata_for_supplement.tsv has age and ancestry component.",
      "",
      "## STAGE 1 TASKS",
      "1. docs/local-ancestry-prereg.md FIRST: the reference populations and exactly which AADR group labels define them, the AIM selection rule, the local-ancestry method and its approximations, the validation criteria, the AMY test statistic, the null, the detection floor, and the decision rule. State up front that this is an unphased, in-house approximation, not RFMix/FLARE.",
      "2. REFERENCE PANEL from AADR published group labels. Choose source populations that bracket West Eurasian variation - at minimum hunter-gatherer (WHG and/or EHG), Anatolian/early-farmer, and Iran/CHG-related (and steppe if a clean label exists). Report the exact labels used and the number of individuals per source, and their date ranges. Exclude relatives/duplicates if the labels or the relatedness screen allow.",
      "3. AIM PANEL: select autosomal SNPs with high between-source allele-frequency divergence, low missingness and adequate call rate in the targets. Report the selection rule, the number used, and the between-source Fst or divergence per selected SNP. Save the panel.",
      "4. LOCAL-ANCESTRY METHOD (implement, do not just call something): a window-based or HMM-based painting where the state is the source population and the emission uses per-window allele counts against source allele frequencies. An explicit HMM over windows with a switch penalty is preferred over independent window calls, because it produces segment-like structure. Document: window size, the emission model and why independence across SNPs is an approximation (LD makes it wrong in a known direction), the switch penalty, and how you handle missing data and low coverage. Report sensitivity to window size and switch penalty.",
      "5. VALIDATION - REQUIRED before any AMY claim: (a) individuals with an unambiguous published group label must be inferred as predominantly that source (report the confusion matrix); (b) inferred GLOBAL ancestry (genome-wide average of local calls) must correlate with the published labels and with the PC position; (c) stability across window sizes/penalties; (d) a negative control - a putatively unadmixed reference individual should show one dominant ancestry. If validation fails, say so and STOP the AMY part rather than reporting numbers from a broken painter. Report exactly how many target individuals you could paint.",
      "6. THE AMY TEST: for each target individual, estimate local ancestry in a window around chr1:103.5-104.5 Mb. Then:",
      "   (a) report AMY1 structural haplotype frequency (or the binarised long-allele dosage) BY INFERRED LOCAL ANCESTRY and BY TIME. Does the frequency change with time WITHIN a fixed local ancestry?",
      "   (b) fit dosage ~ local_ancestry + time and test whether the TIME term survives conditioning on local ancestry. That is the direct answer to 'is the rise explained by ancestry composition, or does it persist within an ancestry background?'",
      "   (c) report the detection floor for that test given n and the allele frequency, BEFORE interpreting the result.",
      "   (d) compare against the global-label result already established (docs/amy-binarized-results.md, docs/two-channel-results.md: carriage ~94% of pooled change, within-history below detection floor, 225 matched ancients).",
      "7. docs/local-ancestry-results.md - the reference panel, AIM panel, method, validation (including the confusion matrix and the negative control), the AMY analysis, the detection floor, and an explicit verdict: DOES THE AMY SIGNAL PERSIST WITHIN A FIXED LOCAL-ANCESTRY BACKGROUND? State plainly whether the local-ancestry view changes the carriage-dominated conclusion, and whether the power is sufficient to tell.",
      "",
      "## CARRIED-FORWARD FACTS TO VERIFY, NOT ASSUME",
      "- Mean AMY1 by ancestry component: WHG 5.75, Early farmers 5.45, EHG 5.39, Neolithic farmers 7.15, Steppe pastoralists 8.25, Bronze age 6.63, Iron age to modern 7.25.",
      "- Binarised AMY: A = shortest/ancestral haplotypes, B = all longer; the A<=2 cut gives X-equivalent 3.64 (clears FDR-50% 3.16, not 5.45), GRM retention ~97% (no structure absorption), carriage ~94% of pooled change, within-history slope 0.049 B-freq/kyr vs 0.146 80%-power floor.",
      "- Akbari's AMY window: 2,204 SNPs, max |X| = 3.794956889, zero above 5.45; LCT 16.52, FADS 16.15, SLC45A2 13.15, SLC24A5 5.99.",
      "- 225/288 labellable AMY ancients match the imputed bundle; 522 AMY ancients map to AADR with dates and coordinates.",
      "",
      "## WHAT SUCCESS LOOKS LIKE",
      "A validated local-ancestry painter, a confusion matrix showing it works on known individuals, and a clear answer on whether the AMY frequency change survives within a fixed local-ancestry background - with the detection floor stated so an underpowered null is not presented as a negative. A validated machinery with an underpowered final test is an acceptable Stage 1 outcome; an unvalidated painter is not.",
      "",
      "## Validation",
      "- docs/local-ancestry-prereg.md exists before any result and names the source labels, the AIM rule, the method and its approximations.",
      "- The validation section contains a confusion matrix, a stability result and a negative control, or states that validation failed and the AMY part was stopped.",
      "- The AMY test reports the time term conditioned on local ancestry WITH its detection floor.",
      "- Every claim about what local ancestry changes is supported by the numeric comparison against the global-label result.",
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
