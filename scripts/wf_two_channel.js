const lanes = [
  {
    key: "two-channel",
    agent: "worker",
    model: "openai-codex/gpt-6.1-sol",
    task: [
      "Repo: /home/erikg/demoselex. Decompose each locus into TWO channels and report BOTH; do not use one to cancel the other.",
      "",
      "HARD RULES",
      "- Do NOT run any git command; parent commits. Append to docs with clearly-marked sections; never rewrite another author's conclusions.",
      "- Write ONLY under analysis/, docs/, fig/, data/derived/, data/raw/annotation/.",
      "- .venv/bin/python. Label claims [observation]/[inference]/[speculation].",
      "- PRE-REGISTER the decomposition and the ordered-series test (docs/two-channel-prereg.md) BEFORE computing, and say it was written before results.",
      "- Stop by ~35 minutes; write partial results incrementally. Never run git.",
      "",
      "## THE POINT",
      "The previous passes built only the WITHIN-deme channel and treated between-deme (carriage) change as nuisance. That produced a conservatism posture ('selection explains ~2%'). The operator's requirement is INCLUSIVE: report both channels, plus the ordering evidence across the donor series, and rank on the combination. A series of distinct incoming populations that carry the allele at progressively HIGHER frequency is potential selection and must be capturable, not explained away.",
      "",
      "## WHAT TO BUILD",
      "1. docs/two-channel-prereg.md - specify before computing:",
      "   (a) the WITHIN channel: within-deme slope of frequency (or copy number) on time, with measured group baselines and within-group time centring (reuse analysis/mixture_null.py / analysis/prioritization_score.py; read docs/mixture-null-prereg.md). Call it beta_within.",
      "   (b) the BETWEEN / CARRIAGE channel: the change predicted by ancestry turnover alone, beta_between. Estimate it from measured donor means and the ancestry proportions per time slice, e.g. fit observed ~ sum_k w_k(t)*mu_k with mu_k measured from component members, and read off the time trend of the composition-weighted mean. NO per-locus free donor frequency.",
      "   (c) the total: beta_total = beta_within + beta_between, reported alongside both components. Nothing is zeroed.",
      "   (d) the ORDERED-SERIES test: order the donor components by mean date; ask whether their measured means progress monotonically (or largely monotonically - allow partial credit, e.g. Spearman or the number of concordant successive steps) and how unlikely that ordering is under a drift null that accounts for donor non-independence. State how you use a relatedness/phylogenetic correlation among donors (analysis/traj_relatedness.py, data/derived/traj_relatedness*.tsv) and how you discount nested donors. This is the operator's original argument and it must be a real test with a real null.",
      "   (e) the ranking rule for prioritisation: inclusive, combining the two channels and the ordering evidence, with the false-positive cost stated. And the decision rule for what would falsify the claim that carriage-with-ordering is selection evidence.",
      "2. analysis/two_channel.py - implements it deterministically with recorded seeds. Applies to both (i) SNP frequency genome-wide and (ii) AMY1 copy number.",
      "3. Genome-wide run: use the merged prioritization panels at data/derived/prio_panels/{shard,moderate}*_out/priority_loci.tsv (538,262 scored variants; all 31,939 Akbari-strong, all 356,329 moderate, 150k controls) plus the Akbari statistics (analysis/load_akbari.py; |X|>5.45 = pi>99%, |X|=3.16 = FDR 50%, X is NOT centred). Report beta_within, beta_between, beta_total, and the ordering statistic per locus; and the counts where the two channels DISAGREE (large between, small within = carriage-class; large within = within-deme class).",
      "4. AMY case study - this is the motivating positive example and must be treated as such, not explained away:",
      "   - Verify independently that Akbari's released statistics do NOT detect the locus: the operator's parent measured max |X| = 3.795 over 2,204 SNPs in chr1:103.5-104.5 Mb (vs 16.5 LCT, 16.1 FADS, 13.2 SLC45A2, 6.0 SLC24A5) against a 5.45 threshold. Recompute it and say if it differs.",
      "   - Compute both channels for AMY1 copy number using the Bolognini/Sudmant/Garrison release in data/raw/amy/ (and the analysis archive inside amylase_diversity_project-v1.1.zip, extract to a scratch dir). Known values to verify: mean AMY1 by ancestry component WHG 5.75, Early farmers 5.45, EHG 5.39, Neolithic farmers 7.15, Steppe pastoralists 8.25, Bronze age 6.63, Iron age to modern 7.25; within-component slopes are all indistinguishable from zero. So beta_within ~ 0 and beta_between is large.",
      "   - Then run the ORDERED-SERIES test on those components and report whether the progression is unlikely under drift given donor non-independence. State the donor independence caveat honestly (West Eurasian components are partly nested: steppe carries CHG/EHG ancestry).",
      "   - Also report the global replication context: PEL 9.52 / Pima 10.29 vs Surui 5.25, and that high-copy haplotypes are shared (H5 0.44 East Asia, 0.18 Europe, 0.46 PEL, 0.71 Pima) so carriage is a live explanation that the ordering test must confront rather than ignore.",
      "5. docs/two-channel-results.md - both channels with numbers, the disagreement classes, the AMY case, the false-positive cost, and an explicit statement of what this does and does not establish. Understand that a large beta_between with a significant ordered progression is a CANDIDATE selection signal, not a dismissal, and say so - while also stating when carriage without ordering is not evidence.",
      "",
      "## CARRIED-FORWARD LIMITS (state them, do not rediscover)",
      "- The score/label reuse: agreement with Akbari is NOT independent validation (same imputed genotypes).",
      "- LD: report block-level (250 kb) results, not only variant-level.",
      "- Rank instability: earlier top-N bootstrap retention was 43.7%.",
      "- A calibrated LD-honest block drift-null was NOT achievable with the existing simulation and remains unavailable; do not paper over it.",
      "",
      "## Validation",
      "- docs/two-channel-prereg.md exists with (a)-(e) before any result appears.",
      "- data/derived/two_channel.tsv has per-locus beta_within, beta_between, beta_total and the ordering statistic.",
      "- The AMY section reports both channels and the ordered-series result with numbers.",
      "- Every genome-wide claim reports the control-locus pass fraction.",
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
