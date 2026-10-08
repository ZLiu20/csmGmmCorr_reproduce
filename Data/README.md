# Data

The applications use public GWAS summary statistics. They are not included because of their
size (about 25 GB).

## Download

```
bash Data/download_data.sh
```

This writes the files below to `Data/raw/`. You can also download them by hand to the same
paths.

| Analysis | Source | Files |
|---|---|---|
| GLGC lipid pleiotropy (Supplementary Section 8.6) | GLGC 2021 trans-ancestry results, https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/ | `with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_{HDL,LDL,TC,logTG,nonHDL}_INV_ALL_with_N_1.gz` |
| GLGC replication (Section 6) | GLGC 2021 sex- and ancestry-specific results, https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/sex_and_ancestry_specific_summary_stats/ | `LDL_INV_{AFR_HRC_1KGP3_others,EAS_1KGP3,EUR_HRC_1KGP3_others,HIS_1KGP3,SAS_HRC_1KGP3_others}_{FEMALE,MALE}.meta.singlevar.results.gz` |
| UK Biobank pleiotropy (Section 6) | Neale lab round 2, https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/ | `{21001,4080,30870,30760,30750}_irnt.gwas.imputed_v3.both_sexes.tsv.bgz` (BMI, SBP, TG, HDL, HbA1c) |

## Preparation

```
bsub < Data/run_prepare.lsf
```

This runs three jobs. Each needs up to 8 hours and 64 GB of memory.

| Job | Command | Output in `Data/prepared/` | Expected variants (fitting) |
|---|---|---|---|
| 1 | `Rscript Data/prepare_glgc.R replication` | `glgc_replication/` | 4,473,732 (258,705) |
| 2 | `Rscript Data/prepare_glgc.R pleiotropy` | `glgc_pleiotropy/` | 9,917,860 (264,749) |
| 3 | `Rscript Data/prepare_ukb.R` | `ukb_pleiotropy/` | 9,370,133 (264,709) |

The GLGC steps are described in Supplementary Section 8.2 and implemented in
`SupportingCode/glgc_prepare.R`. The coordinates of each analysis are listed in
`Data/glgc_*_manifest.tsv`, and the steps are:

- variants are matched on the build-37 position and both alleles (ALT is the effect allele);
- z = beta / SE;
- minor allele frequency at least 0.01;
- N at least half of the largest N of the file;
- variants present in every coordinate;
- statistics capped at ±12;
- one fitting variant per 10-kb window, chosen by a hash of the variant key.

The UK Biobank steps are described in Supplementary Section 8.1:

- autosomal variants that are not low-confidence for any trait;
- minor allele frequency at least 0.01;
- the t-statistic;
- one random fitting variant per 10-kb window, with a fixed seed.

## Format of the prepared files

Each folder in `Data/prepared/` contains:

- `chr<k>.tsv.gz`: one row per variant. Columns:
  - `variant_id`: `chr<chr>:<pos>:<ref>:<alt>` for GLGC, `<chr>:<pos>:<ref>:<alt>` for UK Biobank;
  - `rsid` (GLGC only);
  - `chr` and `pos_b37`;
  - `ref` and `alt` (GLGC only);
  - one z-statistic column per coordinate, named as in the manifest or BMI, SBP, TG, HDL, HbA1c.
- `chromosome_manifest.tsv`: `chr`, `path` and `n_rows`.
- `fit_sample.tsv.gz`: the fitting variants, in the same format. The test set of an analysis is
  every variant that is not in its fitting sample.
- GLGC only: `source_qc.tsv` (counts removed by each filter, per source) and
  `prepare_summary.tsv`.

The replication files contain ten strata. The analysis uses EUR, EAS and AFR (female and
male). The HIS and SAS strata are used only for the held-out validation
(`Fig4/holdout_validation.R`).
