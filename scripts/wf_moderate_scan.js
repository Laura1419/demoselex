const lanes = [0,1,2,3].map(k => ({
  key: `moderate-${k}`,
  agent: "worker",
  task: [
    `Repo: /home/erikg/demoselex. PURE COMPUTE lane - do not interpret science, do not analyse.`,
    ``,
    `Run exactly this command and let it finish (expect a few minutes):`,
    `  .venv/bin/python analysis/prioritization_score.py --panel-file data/derived/prio_panels/moderate${k}.tsv --out-dir data/derived/prio_panels/moderate${k}_out --bootstrap 0 --simulations 0`,
    ``,
    `DO NOT run any git command. Do not modify any source file. Do not write outside data/derived/prio_panels/moderate${k}_out (the script writes there itself).`,
    `Use .venv/bin/python. If it fails, report the FULL traceback verbatim and stop - do not debug the script.`,
    ``,
    `When it finishes report: (1) input row count of data/derived/prio_panels/moderate${k}.tsv, (2) output row count of data/derived/prio_panels/moderate${k}_out/priority_loci.tsv, (3) number of output rows with akbari_abs_X > 5.45 (expect 0), (4) wall-clock seconds.`,
    `Write one line with those numbers to analysis/prioritization-moderate${k}.log. That is the entire deliverable.`
  ].join("\n")
}));

const results = await runs.all(lanes);
return results.map((r, i) => ({
  key: `moderate-${i}`,
  status: r && r.status,
  runId: r && (r.runId || r.id),
  output: r && r.output
}));
