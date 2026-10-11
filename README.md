# csmGmmCorr reproduction archive

This repository reproduces the simulations, data analyses, figures, and tables for
*Correlated csmGmm: Composite Null Testing and Posterior Ranking in Genome-Wide
Studies* and its Supplementary Material. The statistical method is provided by
the separate `csmGmmCorr` R package (version 1.3.0); this repository contains only
the reproducibility workflow and the comparison methods used in the paper.

## Complete run

On an LSF login node with internet access, a fresh clone can be started with:

```bash
git clone https://github.com/ZLiu20/csmGmmCorr_reproduce.git
cd csmGmmCorr_reproduce
mkdir -p logs
nohup bash run_all.sh > logs/run_all.log 2>&1 &
```

`run_all.sh` installs dependencies, compiles two native helper kernels, submits all
simulation arrays, downloads public summary statistics, schedules the real-data
analyses through LSF dependencies, and finally runs `make_all.sh`. Finished
replicates are detected by their CSV files, so the same command can resume an
interrupted run without regenerating completed work.

Useful modes are:

```bash
bash run_all.sh status
bash run_all.sh final
bash run_all.sh make
DRYRUN=1 bash run_all.sh
```

Simulation seeds depend only on scenario, design, dimension, threshold, replicate,
and data role. Every method evaluated in the same cell therefore receives the same
generated data. Floating-point reductions can differ slightly across machines or
thread schedules; the numerical tolerances used by the validation scripts account
for such changes.

## Setup and execution

Install the R dependencies and compile the native helpers once:

```bash
Rscript setup.R
```

A local package archive can be supplied as the first argument or through
`CSMGMMCORR_SRC`. Otherwise, setup installs `ZLiu20/csmGmmCorr` from GitHub. Package
versions are recorded in `logs/versions.txt`.

Every simulation array has an R entry point named `<array>_sim.R` and an LSF file
named `run_<array>.lsf`. Submit from the repository root, for example:

```bash
bsub < Fig2/run_main_ours_pl_pen.lsf
```

Each array task covers ten consecutive replicates of a single `(design, K, t)`
cell. The default experiment uses 100,000 features. For a local smoke run, set
`SIM_N=2000`; smoke output is isolated under `output_smoke/`. `SIM_MAX_REPS` limits
the number of replicates executed by one task, and `MAX_ITER` overrides the fitting
iteration cap.

## Figures and tables

The major products are organized as follows.

| Paper item | Folder | Main purpose |
|---|---|---|
| Figure 1 | `Fig1` | Exact illustration of posterior ranking |
| Figure 2 | `Fig2` | Main and symmetric simulations for all three targets |
| Table 1 | `Tab1` | Design summary |
| Figure 3 | `Fig3` | Correlation sensitivity |
| Table 2 | `Tab2` | Alternative estimators and data-generating mechanisms |
| Tables 3 and Figures 4--5 | `Fig4` | Real-data and holdout analyses |
| Supplementary Figures 1--12 | `SuppFig1` | Additional simulation panels |
| Supplementary Tables 2--26 | `SuppTab*` | Counts, sensitivity analyses, and diagnostics |
| Package validation | `Validation` | Comparison with the independent `csmGmm` implementation |

Once all required jobs finish, run:

```bash
bash make_all.sh
```

The script rebuilds cached summaries, runs every plotting/table program in isolation,
and writes one log per program under `logs/make/`. A failure is reported but does not
prevent the remaining products from being attempted.

## Simulation design

The primary experiments use `n=1000` subjects, 100,000 variants, dimensions
`K=4,6`, thresholds `t=1,2`, and correlated outcomes or mediator errors. The main
latent dependence is AR(1) with correlation 0.2. Additional studies use symmetric
statistics, quantitative outcomes, exchangeable correlation 0.4, stronger
correlation, linked genotypes, and independent test data. Penalized and maximum-
likelihood correlated csmGmm fits are compared with target-appropriate methods:
QCH and adaFilter for pleiotropy, repfdr and directional partial conjunction for
replication, and DEI, MediationTest, and a parameter-free pathway safeguard for
mediation.

Approximate run-time ceilings per array task are 16 hours for the heaviest Figure 2
and component-sensitivity fits, 12 hours for the remaining simulation arrays, four
hours for the weight-floor study, and up to 23 hours for real-data fits after data
preparation. Actual elapsed time depends strongly on cluster concurrency.

## Public data

The applications use public summary statistics from the Global Lipids Genetics
Consortium, ten ancestry-by-sex LDL-C strata, the Neale Laboratory UK Biobank GWAS,
GTEx liver eQTL resources, and 1000 Genomes reference data. The preparation entry
points are:

```bash
bash Data/download_data.sh
bsub < Data/run_prepare.lsf
bsub < Fig4/run_realdata.lsf
```

See `Data/README.md` for file-level sources, checks, and expected outputs. Raw data
are not committed to this repository.

## Repository map

- `SupportingCode/sim_data.R`: designs, seeds, subject-level generation, and truth.
- `SupportingCode/sim_tasks.R`: all simulation arrays and LSF task grids.
- `SupportingCode/run_sim_task.R`: fitting, decisions, performance, and diagnostics.
- `SupportingCode/comparators/`: target-specific comparison methods.
- `SupportingCode/incongruity.R`: decision- and score-level incongruity counts.
- `SupportingCode/summarize.R`: output aggregation and table helpers.
- `SupportingCode/realdata.R`: real-data fitting and inference.
- `SupportingCode/glgc_prepare.R`: GLGC harmonization utilities.
- `SupportingCode/floor_study.R`: the mixture-weight floor experiment.
- `SupportingCode/src/`: native kernels for large-scale adjustment and hashing.
- `Data/`: downloads and preparation.
- `Validation/`: numerical comparisons with Sun et al.'s `csmGmm` package.
- `tests/`: syntax, preparation, and smoke checks.

The reference software environment is R 4.5.1 with `csmGmmCorr` 1.3.0. The
computations described in the paper were run under LSF on the Seadragon cluster at
MD Anderson Cancer Center.

## License

The repository is distributed under the MIT License. The bundled
`SupportingCode/packages/mediation.test` source retains its own license notices.
