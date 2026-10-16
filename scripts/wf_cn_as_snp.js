const lanes = [
  {
    key: "cn-as-snp",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex. Decisive counterfactual test: treat the amylase COPY-NUMBER haplotype alleles as if they were SNPs and run the SAME GLMM approach Akbari et al. used, to ask whether those alleles would have cleared their significance threshold had they been assayed as ordinary biallelic variants.",
      "",
      "HARD RULES",
      "- Do NOT run any git command; parent commits. Append to docs with clearly-marked sections; never rewrite another author's conclusions.",
      "- Write ONLY under analysis/, docs/, fig/, data/derived/.",
      "- .venv/bin/python. Label claims [observation]/[inference]/[speculation].",
      "- PRE-REGISTER the procedure and the decision rule in docs/cn-as-snp-prereg.md BEFORE running the counterfactual; say it was written before results.",
      "- Stop by ~35 minutes; write partial results incrementally.",
      "",
      "## WHY THIS IS DECISIVE",
      "Locus chr1:103.5-104.5 Mb (AMY1/AMY2) shows max |X| = 3.795 over 2,204 SNPs in Akbari's released statistics, versus 16.52 (LCT), 16.15 (FADS), 13.15 (SLC45A2), 5.99 (SLC24A5) - i.e. AMY falls in the discarded band between FDR-50% (3.16) and pi>99% (5.45). Two competing explanations: (a) the SNP signal is a WEAK TAG for the copy-number variant, so a real signal is being missed by the assay; (b) the underlying frequency change is modest or structure-explained. Building the haplotype allele directly as a synthetic biallelic locus and running the same GLMM separates these, because the synthetic locus is a PERFECT tag by construction.",
      "",
      "## PROCEDURE",
      "1. docs/cn-as-snp-prereg.md - specify first: the synthetic locus definition, the GLMM, the statistic, the null, the controls, the decision rule, and the interpretation of each outcome.",
      "2. VALIDATION FIRST - do not proceed to AMY until this passes. Reproduce Akbari's released |X| for at least two canonical loci using genome-wide genotypes for the same individuals: LCT/MCM6 (chr2:136,608,646) and one other (FADS chr11 ~61.59 Mb or SLC45A2 chr5:33,951,693). Use the imputed bundle data/raw/akbari/tier2/ (analysis/load_akbari_geno.py, which returns ALT dosage and requires uint32 indices for sample-subset and batch reads; verify with analysis/test_aadr_orientation.py for the AADR reader if you use AADR instead) and the GLMM reimplementation sim/akbari_glmm.py (read it first; reuse rather than rewrite). If your pipeline cannot roughly reproduce those large |X| values, say so and STOP - the whole counterfactual depends on the GLMM being faithful.",
      "3. Build synthetic biallelic loci for the AMY haplotypes. Haplotype labels are in the Bolognini/Sudmant/Garrison release (extract data/raw/amy/amylase_diversity_project-v1.1.zip to a scratch dir): graph_genotyping/assess_graph_genotypes/supplementary_tables/haplotype_structure_for_supplement.tsv gives per-individual H1 and H2, and ancient_metadata_for_supplement.tsv gives each ancient sample's age and ancestry component. For a target haplotype h, code each individual as a dosage in {0,1,2} = number of their two haplotypes equal to h. That IS a perfect-tag biallelic locus. Build these for H1^a (mean AMY1 4.30), H5 (9.01), H7 (11.16), H9 (13.95), and H3^r (6.49).",
      "4. ALSO build the continuous trait locus: mean AMY1 copy number per individual, as a quantitative phenotype, and run the analogous time-trend test.",
      "5. The GLMM must be run with an individual-level relatedness/GRM term built from GENOME-WIDE genotypes for those same individuals (NOT from the AMY locus itself - that would be circular). The AMY ancients map to AADR, so use the pgen/AADR genotypes of the matched individuals; state the match rate and how many individuals you could actually use. If the GRM cannot be built for enough of them, run the model without it AND with it, and report both - the difference is exactly the 'does structure absorb it' question.",
      "6. Null: permute the haplotype labels across individuals (preserving the date distribution and, where possible, the ancestry-component composition) many times; report the null distribution of the statistic and an empirical P for the observed value. Also report the same for a set of control variants.",
      "7. Convert the statistic to the paper's scale where possible: their X = Z/sqrt(CF) with CF = 2.28 (so X = Z/1.51). Report both Z and the X-equivalent, and compare against 3.16 (FDR ~50%) and 5.45 (pi>99%).",
      "",
      "## WHAT TO REPORT (docs/cn-as-snp-results.md)",
      "- Whether the validation step reproduced the canonical loci; if not, that is the headline.",
      "- For each AMY haplotype: dosage counts, frequency by ancestry component and by date bin, the GLMM statistic, the X-equivalent, the empirical P under the label permutation null, and whether it would have cleared 3.16 and/or 5.45.",
      "- The explicit answer to: WOULD the AMY copy-number alleles have passed Akbari's GLMM if they had been assayed as SNPs? With numbers, and with the honest caveat that the synthetic locus is a perfect tag (an upper bound), that n is small (288 labellable ancients), and that passing a threshold in a counterfactual is not the same as the real scan having tested the variant.",
      "- The separate question of whether the signal survives a relatedness/GRM term: report the statistic with and without it. If the signal disappears when the GRM is included, that is the structure-absorption case and must be stated plainly.",
      "",
      "## CARRIED-FORWARD FACTS TO VERIFY, NOT ASSUME",
      "- Mean AMY1 by ancestry component: WHG 5.75, Early farmers 5.45, EHG 5.39, Neolithic farmers 7.15, Steppe pastoralists 8.25, Bronze age 6.63, Iron age to modern 7.25.",
      "- Ancient haplotype frequencies ordered by mean date show H1^a declining (0.77 -> 0.13, 5 of 6 steps) and H5 rising (0.00 -> 0.23).",
      "- Modern: H5 is 0.40-0.50 in East Asia and 0.46/0.71 in PEL/Pima; H9 (highest copy, 13.95) is 0.11-0.17 in PEL/Pima/Quechua versus 0.00-0.03 in East Asia.",
      "- Akbari's AMY window: 2,204 SNPs, max |X| = 3.794956889, zero above 5.45.",
      "",
      "## Validation",
      "- docs/cn-as-snp-prereg.md exists before any counterfactual result.",
      "- The canonical-locus reproduction step is reported with numbers, or its failure is the headline.",
      "- Every AMY haplotype has an X-equivalent and an empirical P, with the GRM-on and GRM-off versions both reported.",
      "- The GRM is built from genome-wide data, never from the AMY locus.",
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
