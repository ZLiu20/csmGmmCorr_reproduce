# csmGmmCorr_reproduce

Code to reproduce the simulations, the data analyses, and every figure and table of the paper
"Correlated csmGmm: Composite Null Testing with Correlated Statistics" by Zili Liu and Ryan Sun,
and of its Supplementary Material.

The method is implemented in the R package **csmGmmCorr** (version 1.3.0), which is in its own
repository. This repository holds only the code that uses the package. It follows the layout of
`csmGmm_reproduce` (Sun, McCaw and Lin, 2025).

## Overview

- Run all commands from the repository root, or from any folder inside it. File locations are
  resolved with the `here` package.
- Set up once (for example on the login node of the cluster):

  ```
  Rscript setup.R
  ```

  This installs csmGmmCorr (from GitHub, or from a local copy with
  `Rscript setup.R /path/to/csmGmmCorr`) and the packages of the comparison methods. It also
  compiles two small kernels in `SupportingCode/src` and checks them against the package.
- The code is written for a cluster with the LSF scheduler.
- There is one folder per figure or table. When a figure or table has no folder of its own,
  its code is in the preceding folder (see the table below). For example, Supplementary
  Figures 1-12 are all made in `SuppFig1`, and Figure 5 is made in `Fig4`.
- Running the computations in a folder:
  - Each simulation array of a folder has a script `<array>_sim.R` and a batch file
    `run_<array>.lsf`.
  - Submit every `.lsf` file of the folder from the repository root, e.g.
    `bsub < Fig2/run_main_ours_pl_pen.lsf`, and wait until all jobs have finished.
  - The comment `# run 1-N` in each `.lsf` file gives the number of tasks. Queue names, memory
    units and the R module name may need to be adapted to your cluster.
- Making the figures and tables: each folder has scripts that start with `plot_` or `make_`.
  Run them with `Rscript` once the jobs have finished. Their PDF and LaTeX output is written
  into the folder.
- Simulation output:
  - Each replicate writes one CSV file to `<folder>/output/<study>/<method type>/`, with one row
    per target, threshold t and method.
  - A resubmitted task skips replicates that are already on disk.
  - Every replicate regenerates its data from fixed seeds (`sim_seed()` in
    `SupportingCode/sim_data.R`), so no data files are stored.
- Testing: setting `SIM_N=2000` (features per data set) writes to `output_smoke/` instead and is
  for testing only. `bash tests/test_all.sh submit` runs one replicate of every array in this
  way.
- Run times:
  - Each simulation task runs 10 replicates of 100,000 features.
  - With 4 threads, the penalized fit of one replicate took a median of about 0.5 minutes at
    K = 4 and 1.5 minutes at K = 6; maximum likelihood took about three times as long. The
    generation of the data takes under a minute (pleiotropy, replication) to about 3 minutes
    (mediation) per replicate.
  - The `-W` values in the `.lsf` files are upper bounds per task. The total time of a folder
    depends on how many tasks run at once (at most 200 per array, set by `%200`).

## Figures and tables

| Paper | Folder | Arrays (tasks) | Script |
|---|---|---|---|
| Figure 1 (Example 1, Section 4.1) | Fig1 | none (exact computation, seconds) | plot_fig1.R |
| Figure 2 | Fig2 | 14 arrays (12,000 tasks): penalized fit (500 replicates), ML fit, sign-symmetric statistics, QCH/adaFilter, DEI/MediationTest/JS-PC | plot_fig2.R |
| Table 1 | Tab1 | uses Fig2 | make_tab1.R |
| Figure 3 | Fig3 | 9 arrays (2,400 tasks): correlation study, penalized and ML fits, csmGmm | plot_fig3.R |
| Table 2 | Tab2 | 16 arrays (2,440 tasks): csmGmm comparison, quantitative outcomes, exchangeable correlation, linkage disequilibrium, independent test set, null-correlation oracle; uses Fig2 and Fig3 | make_tab2.R |
| Table 3, Figures 4-5 | Fig4 | `run_realdata.lsf` (13 fits); requires `Data/` | make_tab3.R, plot_fig4.R, plot_fig5.R, holdout_validation.R |
| Supplementary Figures 1-12 | SuppFig1 | uses Fig2 | plot_suppfig1_12.R |
| Supplementary Table 1 | none | roadmap, text only | |
| Supplementary Table 2 | SuppTab2 | none (exact counts) | make_supptab2.R |
| Supplementary Table 3 | SuppTab3 | 18 arrays (720 tasks): penalty constants | make_supptab3.R |
| Supplementary Table 4 | SuppTab4 | `run_SuppTab4.lsf` (10 tasks): weight floors | make_supptab4.R |
| Supplementary Table 5 | SuppTab5 | none (design table) | make_supptab5.R |
| Supplementary Tables 6-8 | SuppTab6 | uses Fig2, Fig3 | make_supptab6_8.R |
| Supplementary Table 9 | SuppTab9 | 6 arrays (240 tasks): numbers of magnitude components | make_supptab9.R |
| Supplementary Tables 10-11 | SuppTab10 | uses Tab2, Fig2, Fig3 | make_supptab10_11.R |
| Supplementary Tables 12-15 | SuppTab12 | 1 array (480 tasks): repfdr and directional PC; uses Fig2, Fig3, Tab2 | make_supptab12_15.R |
| Supplementary Tables 16-18 | SuppTab16 | uses Fig2 | make_supptab16_18.R |
| Supplementary Tables 19-23 | SuppTab19 | uses Fig4 | make_supptab19_23.R |
| Supplementary Tables 24-26 | SuppTab24 | uses Fig4 LDL mediation (`run_mediation_ldl.lsf`, `run_mediation_ldl_all_t.lsf`, `run_mediation_ldl_compare.lsf`) | make_supptab24_26.R |
| Package checks (Supplementary Section 2.7) | Validation | none | validate_csmgmm.R, validate2_csmgmm.R, diag_incon.R |

Upper bounds on the wall time per task (`-W`):

| Folder | Upper bound per task |
|---|---|
| Fig2 | 16 hours |
| Fig3, Tab2, SuppTab3, SuppTab12 | 12 hours |
| SuppTab9 | 16 hours |
| SuppTab4 | 4 hours |
| Fig4 | 23 hours, after the data preparation (8 hours) |

All other scripts finish in minutes.

## Data (applications)

The applications use public summary statistics:

- the Global Lipids Genetics Consortium (GLGC, 2021) trans-ancestry lipid results;
- the GLGC LDL-cholesterol results of ten ancestry-by-sex strata;
- the round-2 UK Biobank GWAS of the Neale laboratory.

Prepare them as follows (see `Data/README.md`):

```
bash Data/download_data.sh          # about 25 GB
bsub < Data/run_prepare.lsf         # GLGC replication, GLGC lipids, UK Biobank
bsub < Fig4/run_realdata.lsf        # 13 fits
```

## Repository structure

```
SupportingCode/      code sourced by all scripts (load_all.R)
  sim_data.R           data generation, designs and seeds of all simulation studies
  sim_tasks.R          the simulation arrays: folder, study, method, designs, replicates, resources
  run_sim_task.R       one task: data, fits, decisions, FDP, power and incongruity counts
  comparators/         QCH, adaFilter, repfdr, directional PC, DEI, MediationTest, JS-PC
  incongruity.R        incongruity counts (score and decision level)
  floor_study.R        weight-floor study (Supplementary Section 2.6)
  realdata.R           applications: loading, fitting, inference, clumping
  glgc_prepare.R       preparation of the GLGC statistics
  summarize.R          reading and summarizing the simulation output, table helpers
  plot_style.R         common figure style (colors and line types fixed per method)
  src/                 C++ kernels: monotone adjustment and incongruity for millions of
                       variants; BLAKE2b hash for the GLGC fitting sample
  packages/            mediation.test (Chambaz and Miles), used by the MediationTest comparator
Data/                download and preparation of the summary statistics
Fig1, Fig2, ...      one folder per figure or table (see above)
Validation/          comparison of csmGmmCorr with Sun et al.'s csmGmm package
tests/               syntax check and smoke test (tests/test_all.sh)
```

## Software

The paper was produced with:

- R 4.5.1 and csmGmmCorr 1.3.0;
- data.table, ggplot2, cowplot and here;
- qch, repfdr and mediation.test for the comparison methods. adaFilter uses the
  implementation in `SupportingCode/comparators/adafilter.R` when the package is not installed.

The computations ran on the Seadragon cluster of MD Anderson Cancer Center (LSF).

## License

MIT (see LICENSE). `SupportingCode/packages/mediation.test` keeps its own MIT license.
