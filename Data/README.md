# Data

The applications use public GWAS summary statistics. The raw files are not
distributed with this repository because their combined size is approximately
25 GB.

Run `bash Data/download_data.sh` to create `Data/raw/`, or place the files at
the paths listed in the two GLGC manifests. Then submit `Data/run_prepare.lsf`.
The preparation enforces build-37 position and allele matching, minor allele
frequency at least 0.01, source-specific sample-size filtering, intersection
across coordinates, clipping at absolute Z equal to 12, and one fitting variant
per 10-kb window. Prepared chromosome files and the fitting sample are written
below `Data/prepared/`.

The mediation application uses three liver-expressed genes and LDL cholesterol.
Its inputs are the GTEx v10 liver cis-eQTL catalogue, GENCODE v26 annotation,
GTEx v8 variant lookup and liver MASHR prediction model, plus GLGC LDL summary
statistics. Run `Data/download_mediation_ldl.sh` and then submit
`Data/run_prepare_mediation_ldl.lsf`.
