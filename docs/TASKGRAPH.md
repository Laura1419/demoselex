# Task Graph

Dependency graph and task listing for the trajectory-replication research project.

## Mermaid Flowchart

```mermaid
%% Task Graph for trajectory-replication research project
graph TD
    subgraph recon
        recon_akbari[Reconstruct Akbari et al. 2026 model and]
        recon_lit[Literature map of closest existing metho]
        recon_data[Data inventory and access plan]
    end
    subgraph infra
        repro_harness[Provenance and run-manifest harness]
        taskgraph_render[Render the task graph and decision map]
    end
    subgraph lean
        lean_setup[Bootstrap the Lean 4 + Mathlib project]
        lean_defs[Core definitions: trajectory, frequency ]
        lean_replacement[Replacement vs replicated change; null d]
        lean_collinearity[Time-GRM collinearity implies GLMM ident]
        lean_statistic[Concordance statistic and effective numb]
        lean_power[Power separation between replacement and]
        lean_review[Lean proof-quality and axiom audit]
    end
    subgraph sim
        sim_scaffold[Demographic simulator scaffold and scena]
        sim_glmm[Reimplement the Akbari-style GLMM statis]
        sim_stat[Implement the trajectory replication sta]
        sim_sweeps[FPR, power, bias, and robustness sweeps ]
        sim_ctrl_null[Matched-control empirical null]
    end
    subgraph data
        data_access[Acquire Akbari per-SNP statistics and su]
        data_meta[Sample metadata: dates, coordinates, anc]
        data_geno[Genotype/dosage acquisition and QC]
    end
    subgraph trajectory
        traj_prototype[Prototype trajectory inference: PCA + ti]
        traj_alt[Alternative inference: temporal kNN grap]
        traj_validate[Trajectory validation and stability (Dec]
        traj_relatedness[Trajectory genetic relatedness and effec]
    end
    subgraph locus
        locus_canonical[Canonical-locus validation (~20-50 known]
        locus_amy[AMY1/AMY2 motivating case study]
        locus_akbari[Akbari reanalysis and four-quadrant clas]
    end
    subgraph figures
        locus_figures[Principal regenerable figures]
    end
    subgraph synthesis
        synth_assessment[Final assessment (J): supported / partly]
        synth_paper[Draft manuscript]
    end
    subgraph reproducibility
        repro_verify[Reproducibility verification from a clea]
    end

    %% Dependencies
    lean_setup --> lean_defs
    lean_defs --> lean_replacement
    lean_defs --> lean_collinearity
    recon_akbari --> lean_collinearity
    lean_defs --> lean_statistic
    lean_replacement --> lean_power
    lean_statistic --> lean_power
    lean_replacement --> lean_review
    lean_collinearity --> lean_review
    lean_statistic --> lean_review
    lean_power --> lean_review
    recon_akbari --> sim_scaffold
    recon_akbari --> sim_glmm
    lean_statistic --> sim_stat
    sim_scaffold --> sim_sweeps
    sim_glmm --> sim_sweeps
    sim_stat --> sim_sweeps
    sim_stat --> sim_ctrl_null
    recon_data --> data_access
    recon_data --> data_meta
    recon_data --> data_geno
    data_access --> data_geno
    data_meta --> traj_prototype
    data_geno --> traj_prototype
    data_meta --> traj_alt
    data_geno --> traj_alt
    traj_prototype --> traj_validate
    traj_alt --> traj_validate
    traj_validate --> traj_relatedness
    lean_statistic --> traj_relatedness
    traj_validate --> locus_canonical
    sim_stat --> locus_canonical
    traj_validate --> locus_amy
    data_access --> locus_amy
    traj_relatedness --> locus_akbari
    sim_sweeps --> locus_akbari
    locus_canonical --> locus_akbari
    sim_ctrl_null --> locus_akbari
    locus_akbari --> locus_figures
    locus_amy --> locus_figures
    locus_akbari --> synth_assessment
    synth_assessment --> synth_paper
    synth_paper --> repro_verify
```

## Task Table

| ID | Title | Dependencies | Tags | Hours | Deliverable |
|----|-------|--------------|------|-------|-------------|
| `recon-akbari` | Reconstruct Akbari et al. 2026 model and available data | — | `recon`, `akbari` | 6.0 | `docs/C-akbari-reconstruction.md` |
| `recon-lit` | Literature map of closest existing methods | — | `recon`, `literature` | 6.0 | `docs/B-literature-map.md` |
| `recon-data` | Data inventory and access plan | — | `recon`, `data` | 4.0 | `docs/E-data-inventory.md` |
| `repro-harness` | Provenance and run-manifest harness | — | `infra`, `reproducibility` | 3.0 | `scripts/provenance.py` |
| `taskgraph-render` | Render the task graph and decision map | — | `infra`, `planning` | 1.5 | `docs/TASKGRAPH.md` |
| `lean-setup` | Bootstrap the Lean 4 + Mathlib project | — | `lean`, `infra` | 3.0 | `lean/lakefile.lean` |
| `lean-defs` | Core definitions: trajectory, frequency path, concordance | `lean-setup` | `lean`, `formalization` | 6.0 | `lean/TrajectorySelection/Defs.lean` |
| `lean-replacement` | Replacement vs replicated change; null distribution theorem | `lean-defs` | `lean`, `formalization`, `decision-2` | 10.0 | `lean/TrajectorySelection/Replacement.lean` |
| `lean-collinearity` | Time-GRM collinearity implies GLMM identifiability loss | `lean-defs`, `recon-akbari` | `lean`, `formalization`, `akbari` | 10.0 | `lean/TrajectorySelection/Collinearity.lean` |
| `lean-statistic` | Concordance statistic and effective number of independent trajectories | `lean-defs` | `lean`, `formalization` | 8.0 | `lean/TrajectorySelection/Statistic.lean` |
| `lean-power` | Power separation between replacement and replication | `lean-replacement`, `lean-statistic` | `lean`, `formalization` | 8.0 | `lean/TrajectorySelection/Power.lean` |
| `lean-review` | Lean proof-quality and axiom audit | `lean-replacement`, `lean-collinearity`, `lean-statistic`, `lean-power` | `lean`, `review` | 4.0 | `lean/AXIOMS.md` |
| `sim-scaffold` | Demographic simulator scaffold and scenarios | `recon-akbari` | `sim`, `demography` | 8.0 | `sim/scenarios/` |
| `sim-glmm` | Reimplement the Akbari-style GLMM statistic | `recon-akbari` | `sim`, `akbari` | 8.0 | `sim/akbari_glmm.py` |
| `sim-stat` | Implement the trajectory replication statistic | `lean-statistic` | `sim`, `statistic` | 6.0 | `sim/trajectory_stat.py` |
| `sim-sweeps` | FPR, power, bias, and robustness sweeps (Decision 3) | `sim-scaffold`, `sim-glmm`, `sim-stat` | `sim`, `validation`, `decision-3` | 12.0 | `docs/G-validation-report.md` |
| `sim-ctrl-null` | Matched-control empirical null | `sim-stat` | `sim`, `validation` | 6.0 | `sim/empirical_null.py` |
| `data-access` | Acquire Akbari per-SNP statistics and supplementary tables | `recon-data` | `data`, `akbari` | 3.0 | `data/raw/akbari/` |
| `data-meta` | Sample metadata: dates, coordinates, ancestry | `recon-data` | `data`, `metadata` | 5.0 | `data/derived/sample_meta.parquet` |
| `data-geno` | Genotype/dosage acquisition and QC | `recon-data`, `data-access` | `data`, `qc` | 10.0 | `analysis/qc.py` |
| `traj-prototype` | Prototype trajectory inference: PCA + time-aware clustering | `data-meta`, `data-geno` | `trajectory` | 10.0 | `analysis/trajectory_infer.py` |
| `traj-alt` | Alternative inference: temporal kNN graph / GMM in time windows | `data-meta`, `data-geno` | `trajectory` | 8.0 | `analysis/trajectory_alt.py` |
| `traj-validate` | Trajectory validation and stability (Decision 1) | `traj-prototype`, `traj-alt` | `trajectory`, `validation`, `decision-1` | 8.0 | `docs/traj-validation.md` |
| `traj-relatedness` | Trajectory genetic relatedness and effective N | `traj-validate`, `lean-statistic` | `trajectory`, `statistic` | 6.0 | `analysis/traj_relatedness.py` |
| `locus-canonical` | Canonical-locus validation (~20-50 known loci) (Decision 4a) | `traj-validate`, `sim-stat` | `locus`, `validation` | 8.0 | `analysis/canonical_loci.py` |
| `locus-amy` | AMY1/AMY2 motivating case study | `traj-validate`, `data-access` | `locus`, `amy` | 10.0 | `docs/H-amy-analysis.md` |
| `locus-akbari` | Akbari reanalysis and four-quadrant classification (Decision 5) | `traj-relatedness`, `sim-sweeps`, `locus-canonical`, `sim-ctrl-null` | `locus`, `akbari`, `decision-5` | 14.0 | `docs/I-akbari-reanalysis.md` |
| `locus-figures` | Principal regenerable figures | `locus-akbari`, `locus-amy` | `figures` | 6.0 | `fig/` |
| `synth-assessment` | Final assessment (J): supported / partly / unsupported | `locus-akbari` | `synthesis` | 5.0 | `docs/J-final-assessment.md` |
| `synth-paper` | Draft manuscript | `synth-assessment` | `synthesis`, `writing` | 10.0 | `docs/paper/` |
| `repro-verify` | Reproducibility verification from a clean checkout | `synth-paper` | `reproducibility`, `review` | 6.0 | `docs/repro-check.md` |

## Legend

- **ID**: Task identifier used in WG commands
- **Title**: Human-readable description
- **Dependencies**: Tasks that must complete before this task can start
- **Tags**: Categorization labels
- **Hours**: Estimated effort (person-hours)
- **Deliverable**: Expected output file/path

*Generated by `scripts/render_taskgraph.py`*