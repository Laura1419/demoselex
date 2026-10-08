#!/usr/bin/env python3
"""Build the WorksGood task graph for the trajectory-replication project.

Declarative and idempotent: re-running skips tasks that already exist.
Usage:  python3 scripts/build_taskgraph.py [--rebuild]
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

ROOT = "/home/erikg/demoselex"


@dataclass
class T:
    id: str
    title: str
    deps: list[str] = field(default_factory=list)
    deliverable: str = ""
    desc: str = ""
    tags: list[str] = field(default_factory=list)
    hours: float = 2.0
    validation: str = ""


TASKS: list[T] = [
    # ---------------- WS0 recon / infrastructure ----------------
    T("recon-akbari", "Reconstruct Akbari et al. 2026 model and available data",
      [], "docs/C-akbari-reconstruction.md",
      "Read the full Nature 2026 paper, all Supplementary Information sections (S1-S6), the "
      "supplementary tables, the code repository, and the ENA data-release docs. Document: sample "
      "size; geographic scope; temporal range; genotype/imputation pipeline; variants tested; "
      "individual clustering; GRM construction; leave-one-chromosome-out; the GLMM likelihood; "
      "definition and units of the selection coefficient; locus-specific variance handling; "
      "calibration; filters; released per-SNP statistics; public genotype/dosage data; sample dates; "
      "coordinates/region labels; individual ancestry. Explain mathematically what is downweighted "
      "when the calendar-time vector aligns with major GRM eigenvectors.",
      ["recon", "akbari"], 6.0,
      "Doc exists with every bullet answered or explicitly marked 'not released'; likelihood written "
      "out; a worked derivation of the time-GRM collinearity effect."),
    T("recon-lit", "Literature map of closest existing methods",
      [], "docs/B-literature-map.md",
      "Systematically survey: ancient-DNA allele-frequency time-series selection; admixture-aware "
      "selection inference; population-graph methods; PolyGraph; GRoSS; CLES/CLUES and ancestry-"
      "stratified CLUES; local-ancestry selection; Simon & Coop drift/migration/selection "
      "decomposition; Bayenv/LFMM; convergent/parallel adaptation; phylogenetic comparative "
      "methods; trajectory inference; continuity/replacement detection without archaeological "
      "labels. For EACH: (1) unit of analysis; (2) how structure is represented; (3) what is "
      "nuisance; (4) what counts as selection evidence; (5) can selection differ among populations; "
      "(6) can it detect repeated directionality across histories; (7) adapt to this problem. "
      "End with an explicit novelty statement: what is solved, what remains open.",
      ["recon", "literature"], 6.0,
      "Table with all seven columns per method; explicit novelty statement."),
    T("recon-data", "Data inventory and access plan",
      [], "docs/E-data-inventory.md",
      "Enumerate exact files for: Akbari per-SNP/summary statistics and supplementary tables; "
      "public genotype/dosage data (ENA PRJEB106907); AADR / Reich lab releases (v66 / 1240k / "
      "ho); sample dates and coordinates; AMY structural-variation data (Bolognini/Sudmant/"
      "Garrison). For each: accession, URL, format, size, license, access requirement, "
      "download command. Flag anything not openly available and state the fallback.",
      ["recon", "data"], 4.0,
      "Every source has accession + URL + command + access status; a fallback is stated for each "
      "closed source."),
    T("repro-harness", "Provenance and run-manifest harness",
      [], "scripts/provenance.py",
      "A small library + CLI that, for any run, records: git commit, command line, software "
      "versions, input file hashes, parameters, random seed, outputs. Writes "
      "analysis/runs/<timestamp>-<slug>/manifest.json. Wire it so pipeline steps call it. Keep raw "
      "data immutable; derived data regenerable.",
      ["infra", "reproducibility"], 3.0,
      "Demo run produces a manifest containing git commit, versions, input hashes, seed; "
      "documented in docs/PROVENANCE.md."),
    T("taskgraph-render", "Render the task graph and decision map",
      [], "docs/TASKGRAPH.md",
      "Render the WG graph to a stable document (mermaid + table) and write docs/DECISIONS.md with "
      "the five §20 decision points as explicit go/no-go records. Keep regenerable from "
      "scripts/render_taskgraph.py.",
      ["infra", "planning"], 1.5,
      "docs/TASKGRAPH.md lists every task with deps; docs/DECISIONS.md has a row per decision "
      "point with status and evidence link."),

    # ---------------- WS1 Lean formalization ----------------
    T("lean-setup", "Bootstrap the Lean 4 + Mathlib project",
      [], "lean/lakefile.lean",
      "Create a working `lake` project using the installed toolchain (Lean 4.34.1). Pin Mathlib. "
      "Add TrajectorySelection.lean root importing submodules. Create lean/AXIOMS.md with the "
      "honest accounting policy: every `sorry`, `axiom`, and `native_decide` documented.",
      ["lean", "infra"], 3.0,
      "`lake build` succeeds; root imports resolve; AXIOMS.md exists with policy."),
    T("lean-defs", "Core definitions: trajectory, frequency path, concordance",
      ["lean-setup"], "lean/TrajectorySelection/Defs.lean",
      "Define: a temporal population trajectory as an indexed family with a time map; an "
      "allele-frequency path along a trajectory; sample-size-weighted slope; direction sign; "
      "directional concordance (fraction in majority direction); the replacement-only model "
      "(pooled path = mixture of two constant-frequency sources); and a relatedness/covariance "
      "structure over trajectories. Keep definitions minimal and general enough for the theorems.",
      ["lean", "formalization"], 6.0,
      "Definitions compile with no `sorry`; each definition has a docstring and a plain-language "
      "comment."),
    T("lean-replacement", "Replacement vs replicated change; null distribution theorem",
      ["lean-defs"], "lean/TrajectorySelection/Replacement.lean",
      "Formalize Scenario A and Scenario B and prove (as far as feasible) that under a "
      "replacement-only model (each trajectory constant in time) the concordance statistic is "
      "determined by random signs and therefore matches the exchangeable null, whereas under "
      "replicated directional change with fixed positive drift the probability of high concordance "
      "strictly exceeds the null. State assumptions explicitly.",
      ["lean", "formalization", "decision-2"], 10.0,
      "Main theorem statement compiles; proof either complete or `sorry` documented in AXIOMS.md "
      "with the missing mathematical steps spelled out."),
    T("lean-collinearity", "Time-GRM collinearity implies GLMM identifiability loss",
      ["lean-defs", "recon-akbari"], "lean/TrajectorySelection/Collinearity.lean",
      "Model the Akbari GLMM as allele frequency = f(time) + population-identity random effect. "
      "Prove the linear-algebra statement: if the time vector is (near-)parallel to a dominant "
      "eigenvector of the GRM, the residual (within-population) variance available to identify the "
      "time slope collapses; quantify via the projection 1 - cos^2(theta) or an explicit "
      "variance-inflation factor. This is the precise meaning of 'what is downweighted'.",
      ["lean", "formalization", "akbari"], 10.0,
      "A theorem relating cos similarity between time vector and GRM eigenvector to the residual "
      "variance/inflation factor; proof or documented sorry."),
    T("lean-statistic", "Concordance statistic and effective number of independent trajectories",
      ["lean-defs"], "lean/TrajectorySelection/Statistic.lean",
      "Define the concordance statistic C_l = (#trajectories in majority direction)/(#informative) "
      "and the relatedness-discounted effective count N_eff. Prove basic properties: monotonicity, "
      "range, behaviour under duplication of a trajectory (N_eff must not increase when a "
      "near-identical copy is added), and the exchangeable null mean/variance of C_l.",
      ["lean", "formalization"], 8.0,
      "Properties proven or documented; the duplication-invariance property is explicit."),
    T("lean-power", "Power separation between replacement and replication",
      ["lean-replacement", "lean-statistic"], "lean/TrajectorySelection/Power.lean",
      "Formalize the comparison statistic and prove, under stated conditions, that the expected "
      "concordance (or its N_eff-weighted form) under replicated selection exceeds that under "
      "replacement, while a matched-null calibration keeps the rejection probability bounded. "
      "Identify the exact conditions under which the separation fails.",
      ["lean", "formalization"], 8.0,
      "Statement + proof or documented sorry; a clearly stated 'when it fails' lemma."),
    T("lean-review", "Lean proof-quality and axiom audit",
      ["lean-replacement", "lean-collinearity", "lean-statistic", "lean-power"],
      "lean/AXIOMS.md",
      "Independently audit every Lean file: no hidden axioms, no vacuous hypotheses, every `sorry` "
      "accounted for and justified. Confirm each theorem's hypotheses are non-empty (witnessed by "
      "at least one example instance). Produce a table theorem -> status -> axioms used.",
      ["lean", "review"], 4.0,
      "Table complete; any vacuous-hypothesis or unsound-assumption finding reported and fixed or "
      "escalated."),

    # ---------------- WS2 Simulation ----------------
    T("sim-scaffold", "Demographic simulator scaffold and scenarios",
      ["recon-akbari"], "sim/scenarios/",
      "Build a simulator that produces dated individuals sampled along population trajectories with "
      "genome-wide genotypes, configurable migration, replacement, drift, and selection. Scenarios: "
      "neutral drift; replacement only; one-population selection; shared selection; opposing "
      "selection; selection aligned with turnover. Use SLiM if practical, otherwise a coalescent/"
      "mosaic simulator. Every scenario is a config file with a fixed seed.",
      ["sim", "demography"], 8.0,
      "Each scenario runs end-to-end and emits a dated genotype matrix + a truth table (which "
      "trajectory was selected, effect size)."),
    T("sim-glmm", "Reimplement the Akbari-style GLMM statistic",
      ["recon-akbari"], "sim/akbari_glmm.py",
      "Faithful reimplementation from the reconstruction: allele frequency ~ time with population "
      "random effects / GRM, leave-one-chromosome-out, locus-specific variance, producing a Z or "
      "posterior probability of selection. Validate against any released statistics or the paper's "
      "figures. Must run on simulated data.",
      ["sim", "akbari"], 8.0,
      "Unit tests: on a known no-selection simulation the statistic is calibrated; on a strong "
      "single-locus sweep it is significant; interface documented."),
    T("sim-stat", "Implement the trajectory replication statistic",
      ["lean-statistic"], "sim/trajectory_stat.py",
      "Model 0 -> Model 2: per-trajectory weighted slopes, concordance, relatedness-discounted "
      "meta-analysis (GLS / phylogenetic-style), heterogeneity. Deterministic given seed. Returns a "
      "tidy table (locus, trajectory, slope, se, direction, n, span).",
      ["sim", "statistic"], 6.0,
      "Recovers the correct sign under a designed simulation with |slope| > noise; emits the tidy "
      "table with all required columns."),
    T("sim-sweeps", "FPR, power, bias, and robustness sweeps (Decision 3)",
      ["sim-scaffold", "sim-glmm", "sim-stat"], "docs/G-validation-report.md",
      "For every scenario sweep effect size, sampling unevenness, missing genotype rate, trajectory-"
      "inference error, and rate of ancestral turnover. Measure FPR, power, effect-size bias for "
      "both statistics. The central result: under which forms of true selection does the Akbari "
      "GLMM lose power, and does the trajectory method recover it without unacceptable FPR under "
      "migration alone?",
      ["sim", "validation", "decision-3"], 12.0,
      "Power/FPR tables + figures; explicit answer to the central question; Decision 3 recorded in "
      "docs/DECISIONS.md."),
    T("sim-ctrl-null", "Matched-control empirical null",
      ["sim-stat"], "sim/empirical_null.py",
      "For each candidate locus, select matched control loci on allele frequency, recombination "
      "rate, local LD, chromosome, imputation quality, missingness, temporal coverage. Run the "
      "identical trajectory analysis on controls. Return an empirical percentile / P for the "
      "candidate. Validate the matcher on simulated null data.",
      ["sim", "validation"], 6.0,
      "Empirical-null calibration verified on a simulated neutral locus (uniform-ish P values)."),

    # ---------------- WS3 Data ----------------
    T("data-access", "Acquire Akbari per-SNP statistics and supplementary tables",
      ["recon-data"], "data/raw/akbari/",
      "Download the supplementary tables, per-locus statistics, and any released per-variant "
      "statistics, plus the code repository. Record accession, URL, date, size, sha256 in the "
      "provenance manifest. Never modify raw files.",
      ["data", "akbari"], 3.0,
      "Files present with hashes recorded; a loader reads them into a tidy table."),
    T("data-meta", "Sample metadata: dates, coordinates, ancestry",
      ["recon-data"], "data/derived/sample_meta.parquet",
      "Assemble individual-level metadata (date, latitude/longitude, region, published ancestry "
      "labels for annotation only, coverage) for the ancient individuals available. Document "
      "every join key, source, and any imputed/approximate date.",
      ["data", "metadata"], 5.0,
      "Tidy table with no duplicate individual ids; provenance for every column; count reconciled "
      "against the paper's N."),
    T("data-geno", "Genotype/dosage acquisition and QC",
      ["recon-data", "data-access"], "analysis/qc.py",
      "Acquire permitted genotype or dosage data (AADR 1240k / ho, or the released capture "
      "dataset if accessible). Run QC: missingness per individual and site, contamination proxies "
      "if present, relatedness, sex checks where available, MAF filters, and an explicit "
      "chromosome/locus exclusion policy. Emit a machine-readable QC manifest.",
      ["data", "qc"], 10.0,
      "QC manifest with pass/fail counts; documented exclusion policy; dataset shape matches "
      "metadata."),

    # ---------------- WS4 Trajectory inference ----------------
    T("traj-prototype", "Prototype trajectory inference: PCA + time-aware clustering",
      ["data-meta", "data-geno"], "analysis/trajectory_infer.py",
      "Using genome-wide variants excluding each candidate locus (and preferably the whole "
      "candidate chromosome), infer a modest number of temporal population trajectories from "
      "genetic similarity + dates (+ optionally geography). Impose no archaeological names. Start "
      "with PCA plus time-constrained clustering or a forward-in-time nearest-neighbour graph.",
      ["trajectory"], 10.0,
      "Emits trajectories as (id, ordered sample list); documented parameters; deterministic."),
    T("traj-alt", "Alternative inference: temporal kNN graph / GMM in time windows",
      ["data-meta", "data-geno"], "analysis/trajectory_alt.py",
      "Implement at least one independent method (time-window Gaussian mixtures or a "
      "genetic-distance forward-time graph). The purpose is robustness: if two methods disagree "
      "wildly, Decision 1 is at risk.",
      ["trajectory"], 8.0,
      "Second method produces trajectories in the same schema as traj-prototype."),
    T("traj-validate", "Trajectory validation and stability (Decision 1)",
      ["traj-prototype", "traj-alt"], "docs/traj-validation.md",
      "Visualize in PCA/ancestry space, geography, and time. Compare with known ancient structure. "
      "Check that major known turnovers emerge naturally. Bootstrap/resample to test robustness. "
      "Use archaeological labels only as post-hoc annotation. If trajectories are unstable or "
      "uninterpretable, STOP and improve before any locus analysis.",
      ["trajectory", "validation", "decision-1"], 8.0,
      "Stability metrics + figures; explicit Decision-1 verdict recorded in docs/DECISIONS.md."),
    T("traj-relatedness", "Trajectory genetic relatedness and effective N",
      ["traj-validate", "lean-statistic"], "analysis/traj_relatedness.py",
      "Estimate pairwise genetic relatedness among inferred trajectories from genome-wide data. "
      "Derive an effective number of independent trajectories (or the covariance matrix for GLS). "
      "Five near-identical trajectories must not count as five. Keep the correction interpretable.",
      ["trajectory", "statistic"], 6.0,
      "Relatedness matrix + N_eff per locus; sanity check on duplicated trajectories."),

    # ---------------- WS5 Locus analysis ----------------
    T("locus-canonical", "Canonical-locus validation (~20-50 known loci) (Decision 4a)",
      ["traj-validate", "sim-stat"], "analysis/canonical_loci.py",
      "Analyze LCT, FADS1/FADS2, pigmentation, immune, metabolism, pathogen-response loci and "
      "strongest Akbari signals. For each: Akbari statistic, trajectory-specific slopes, "
      "concordance, empirical-null percentile, heterogeneity. Determine whether behaviour is "
      "biologically interpretable.",
      ["locus", "validation"], 8.0,
      "Table for every locus; interpretation notes; Decision-4 first verdict."),
    T("locus-amy", "AMY1/AMY2 motivating case study",
      ["traj-validate", "data-access"], "docs/H-amy-analysis.md",
      "Obtain ancient AMY structural-variation data/code (Bolognini/Sudmant/Garrison). Relate the "
      "ancient individuals to the Akbari sample set. Because AMY is copy number, use dosage / "
      "structural-haplotype representation. Test: (1) does AMY copy number rise along multiple "
      "independent trajectories? (2) are rises confined to replacements? (3) do several structural "
      "haplotypes rise? (4) do independent mutational backgrounds converge? (5) does the evidence "
      "survive excluding close trajectories? Report whatever the data show - do not assume AMY "
      "validates the method.",
      ["locus", "amy"], 10.0,
      "Doc answers all five questions with evidence; explicit statement of whether AMY supports or "
      "fails to support the method."),
    T("locus-akbari", "Akbari reanalysis and four-quadrant classification (Decision 5)",
      ["traj-relatedness", "sim-sweeps", "locus-canonical", "sim-ctrl-null"],
      "docs/I-akbari-reanalysis.md",
      "Apply to all significant and near-significant Akbari loci; extend genome-wide if feasible. "
      "Per locus: Akbari Z, trajectory-replication statistic, #informative trajectories, N_eff, "
      "common directional effect, heterogeneity, matched-null P/FDR. Scatter plot and classify "
      "into high/low x high/low quadrants. Investigate low-Akbari/high-replication in depth.",
      ["locus", "akbari", "decision-5"], 14.0,
      "Scatter + quadrant table; low-Akbari/high-replication loci individually characterized; "
      "Decision 5 verdict recorded."),
    T("locus-figures", "Principal regenerable figures",
      ["locus-akbari", "locus-amy"], "fig/",
      "One script per figure, all reading derived data, producing vector output: trajectories in "
      "PCA/space/time, simulation power/FPR, canonical-locus forest plots, AMY trajectories, and "
      "the Akbari-vs-replication scatter. Every figure regenerable from code.",
      ["figures"], 6.0,
      "`make figures` regenerates all; figure list in README matches files present."),

    # ---------------- WS6 Synthesis ----------------
    T("synth-assessment", "Final assessment (J): supported / partly / unsupported",
      ["locus-akbari"], "docs/J-final-assessment.md",
      "Using the decision log and simulation + empirical results, state clearly whether the initial "
      "criticism is strongly supported, partly supported, unsupported, or dependent on particular "
      "forms of selection. Distinguish observation / inference / speculation throughout. Do not "
      "force the project toward the original hypothesis.",
      ["synthesis"], 5.0,
      "Doc with an explicit one-sentence verdict and the evidence chain for it."),
    T("synth-paper", "Draft manuscript",
      ["synth-assessment"], "docs/paper/",
      "Assemble the living document into a manuscript draft: methods, results, figures, "
      "limitations, data/code availability. One-paragraph plain-language method summary as "
      "specified in the plan.",
      ["synthesis", "writing"], 10.0,
      "Complete draft with all sections; every figure referenced exists in fig/."),
    T("repro-verify", "Reproducibility verification from a clean checkout",
      ["synth-paper"], "docs/repro-check.md",
      "From a fresh clone, attempt to regenerate the principal results using the documented "
      "commands. Record every failure and ambiguity. Fix or document. Confirm no result depends on "
      "unrecorded state.",
      ["reproducibility", "review"], 6.0,
      "Checklist of every principal figure/table with a clean-run status; failures enumerated."),
]


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=check)


def existing_ids() -> set[str]:
    out = run(["wg", "list", "--all", "--json"], check=False)
    try:
        data = json.loads(out.stdout)
    except Exception:
        return set()
    ids = set()
    if isinstance(data, list):
        for row in data:
            if isinstance(row, dict):
                ids.add(row.get("id") or row.get("task_id"))
    elif isinstance(data, dict):
        for row in data.get("tasks", []):
            ids.add(row.get("id") or row.get("task_id"))
    return {i for i in ids if i}


def main() -> int:
    rebuild = "--rebuild" in sys.argv
    have = existing_ids()
    if rebuild:
        for t in TASKS:
            run(["wg", "abandon", t.id], check=False)
        have = set()

    created: list[T] = []
    for t in TASKS:
        if t.id in have:
            print(f"skip (exists): {t.id}")
            continue
        desc = t.desc
        if t.validation:
            desc += f"\n\n## Validation\n{t.validation}"
        desc += "\n\n## Exit criteria\n"
        desc += ("- Deliverable written at the declared path and committed to git.\n"
                 "- Provenance recorded (see scripts/provenance.py) with seed and versions.\n"
                 "- Relevant docs section updated; docs/JOURNAL.md appended.\n"
                 "- No un-accounted `sorry`/axiom (Lean) and no unseeded randomness.\n")
        # positional TITLE must come first: --after is variadic and would eat it
        cmd = ["wg", "add", t.title, "--id", t.id, "-d", desc,
               "--deliverable", t.deliverable, "--hours", str(t.hours)]
        for tg in t.tags:
            cmd += ["-t", tg]
        for d in t.deps:
            cmd += ["--after", d]
        r = run(cmd, check=False)
        if r.returncode != 0:
            print(f"FAILED {t.id}: {r.stderr.strip() or r.stdout.strip()}")
        else:
            print(f"added: {t.id}")
        created.append(t)

    # publish whole graph (WCC from a sink reaches all connected tasks)
    sinks = [t.id for t in TASKS if not any(t.id in o.deps for o in TASKS)]
    for s in sinks:
        r = run(["wg", "publish", s, "--wcc"], check=False)
        print(f"publish {s}: {r.stdout.strip()[-200:] or r.stderr.strip()[-200:]}")
    
    # Add provenance demonstration task
    provenance_task = T(
        id="provenance-demo",
        title="Demonstrate provenance system",
        desc="""Run the provenance system demo to verify it works correctly.

## Validation
- Demo run produces a manifest at analysis/runs/*/manifest.json
- Manifest contains git commit, versions, input hashes, and seed
- Verify command: python3 scripts/provenance.py demo && python3 scripts/provenance.py verify analysis/runs/*/manifest.json
""",
        tags=["infra", "reproducibility"],
        hours=0.5
    )
    
    if provenance_task.id not in known_ids:
        wg(provenance_task.id, provenance_task.title, 
           desc=provenance_task.desc, 
           tags=provenance_task.tags, 
           hours=provenance_task.hours)
        publish_tasks.append(provenance_task.id)
        print(f"  + {provenance_task.id}")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
