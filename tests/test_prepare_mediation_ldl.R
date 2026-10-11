suppressPackageStartupMessages(library(data.table))

root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[1L]
} else ".")
workspace <- tempfile("mediation_ldl_fixture_")
source_dir <- file.path(workspace, "source")
glgc_dir <- file.path(workspace, "glgc")
output_dir <- file.path(workspace, "output")
dir.create(source_dir, recursive = TRUE)
dir.create(glgc_dir, recursive = TRUE)
on.exit(unlink(workspace, recursive = TRUE), add = TRUE)

write_gzip_table <- function(x, path, sep = "\t") {
  data.table::fwrite(x, path, sep = sep, quote = FALSE)
}

genes <- sprintf("ENSG%02d", 1:7)
model_variant <- sprintf("chr1_%d_A_G_b38", seq(1000, 7000, by = 1000))
weights <- data.table(
  varID = model_variant,
  gene = paste0(genes, ".1"),
  eff_allele = rep(c("G", "A"), length.out = length(genes)),
  weight = c(0.5, -0.3, 0.4, 0.7, -0.2, 0.6, 0.3)
)
extra <- data.table(gene = paste0(genes, ".1"),
                    gene_name = paste0("GENE", seq_along(genes)))
database <- file.path(source_dir, "mashr_liver.db")
sql <- c(
  "CREATE TABLE weights (varID TEXT, gene TEXT, eff_allele TEXT, weight REAL);",
  "CREATE TABLE extra (gene TEXT, gene_name TEXT);",
  vapply(seq_len(nrow(weights)), function(i) {
    sprintf("INSERT INTO weights VALUES ('%s','%s','%s',%.17g);",
            weights$varID[i], weights$gene[i], weights$eff_allele[i],
            weights$weight[i])
  }, character(1)),
  vapply(seq_len(nrow(extra)), function(i) {
    sprintf("INSERT INTO extra VALUES ('%s','%s');",
            extra$gene[i], extra$gene_name[i])
  }, character(1))
)
sql_file <- file.path(workspace, "fixture.sql")
writeLines(sql, sql_file)
if (system2("sqlite3", shQuote(database), stdin = sql_file) != 0L) {
  stop("Could not create the SQLite fixture")
}

covariance <- data.table(
  GENE = paste0(genes, ".1"), RSID1 = model_variant,
  RSID2 = model_variant, VALUE = 1
)
write_gzip_table(covariance, file.path(source_dir, "mashr_liver.txt.gz"), sep = " ")

lookup <- data.table(
  variant_id = model_variant,
  variant_id_b37 = sprintf("1_%d_A_G_b37", seq(1100, 7100, by = 1000))
)
write_gzip_table(lookup, file.path(source_dir, "GTEx_v8_lookup_table.txt.gz"))

make_glgc <- function(scale, swap_second = FALSE) {
  ref <- rep("A", length(genes))
  alt <- rep("G", length(genes))
  if (swap_second) {
    ref[2L] <- "G"
    alt[2L] <- "A"
  }
  data.table(
    CHROM = 1L, POS_b37 = seq(1100, 7100, by = 1000),
    rsID = paste0("rs", seq_along(genes)), REF = ref, ALT = alt,
    EFFECT_SIZE = scale * c(0.02, -0.03, 0.04, -0.01, 0.05, -0.02, 0.03),
    SE = rep(0.01, length(genes)), POOLED_ALT_AF = 0.3, N = 10000
  )
}
write_gzip_table(
  make_glgc(1, FALSE),
  file.path(glgc_dir,
            "LDL_INV_EUR_HRC_1KGP3_others_FEMALE.meta.singlevar.results.gz")
)
write_gzip_table(
  make_glgc(0.8, TRUE),
  file.path(glgc_dir,
            "LDL_INV_EUR_HRC_1KGP3_others_MALE.meta.singlevar.results.gz")
)

gtf <- data.table(
  chr = c(rep("chr1", 5), rep("chr2", 2)), source = "fixture",
  feature = "gene", start = c(5000, 7000, 25000, 12000, 6000, 1000, 4000),
  end = c(9000, 8000, 30000, 15000, 7000, 2000, 6000), score = ".",
  strand = c("+", "-", "+", "+", "+", "+", "+"), frame = ".",
  attribute = sprintf('gene_id "%s.1"; gene_type "protein_coding";', genes)
)
fwrite(gtf, file.path(source_dir, "gencode.v26.GRCh38.genes.gtf"),
       sep = "\t", col.names = FALSE, quote = FALSE)

qtl_rows <- rbindlist(list(
  data.table(chromosome = 1L, position = 5000L,
             gene_id = paste0(genes[c(1, 2, 3, 6)], ".1")),
  data.table(chromosome = 1L, position = 14500L,
             gene_id = paste0(genes[c(4, 5, 6, 2)], ".1")),
  data.table(chromosome = 2L, position = 5000L,
             gene_id = paste0(genes[c(5, 4, 7, 1)], ".1"))
))
qtl_rows[, `:=`(
  variant_id = sprintf("chr%d_%d_A_T", chromosome, position),
  ref = "A", alt = "T", maf = 0.2,
  beta = seq_len(.N) / 20, se = 0.05, type = "SNP"
)]
setcolorder(qtl_rows, c("gene_id", "chromosome", "position", "ref", "alt",
                        "variant_id", "maf", "beta", "se", "type"))
write_gzip_table(qtl_rows, file.path(source_dir, "QTD000266.all.tsv.gz"))

script <- file.path(root, "Data", "prepare_mediation_ldl.R")
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(script, source_dir, glgc_dir),
  env = paste0("MEDLDL_OUT=", output_dir)
)
if (status != 0L) stop("prepare_mediation_ldl.R failed")

twas <- fread(file.path(output_dir, "mediation", "ldl_twas_liver.tsv.gz"))
stopifnot(nrow(twas) == length(genes), setequal(twas$gene, genes),
          all(is.finite(twas$z_twas)))

nearest <- fread(file.path(output_dir, "mediation",
                           "mediation_ldl_liver_nearest.tsv.gz"))
noshare <- fread(file.path(output_dir, "mediation",
                           "mediation_ldl_liver_noshare.tsv.gz"))
required <- as.vector(rbind(paste0("alpha", 1:3), paste0("beta", 1:3)))
stopifnot(nrow(nearest) == 3L, nrow(noshare) == 3L,
          all(required %in% names(nearest)), all(required %in% names(noshare)),
          all(complete.cases(nearest[, ..required])),
          all(complete.cases(noshare[, ..required])))

# Every model contains one variant, so no-share triples must have zero overlap.
stopifnot(all(noshare$shared == 0L))

# The selected feature in each chromosome/window is nearest the 5-kb midpoint.
key <- nearest[, paste(chr, pos_b38, sep = ":")]
stopifnot(setequal(key, c("1:5000", "1:14500", "2:5000")))

# With a one-variant model, S-PrediXcan is sign(weight) times the aligned LDL z.
model_z <- fread(file.path(output_dir, "prepared", "mediation_ldl",
                           "ldl_model_variants.tsv.gz"))
expected <- sign(weights$weight) * model_z$z_alt[match(weights$varID,
                                                       model_z$variant)]
observed <- twas$z_twas[match(genes, twas$gene)]
stopifnot(isTRUE(all.equal(observed, expected, tolerance = 1e-8)))

cat("MEDIATION_PREPARATION_FIXTURE_OK\n")

# Exercise the public real-data workflows on deterministic synthetic inputs.
# This is an interface and artifact contract, not a scientific calibration test.
source(file.path(root, "SupportingCode", "realdata.R"), local = FALSE)

analysis_root <- file.path(workspace, "analysis_root")
prepared_root <- file.path(analysis_root, "Data", "prepared", "glgc_pleiotropy")
mediation_root <- file.path(analysis_root, "Data", "mediation")
dir.create(prepared_root, recursive = TRUE)
dir.create(mediation_root, recursive = TRUE)

set.seed(20261011)
statistic_cols <- RD_ANALYSES$glgc_pleiotropy$cols
sigma <- matrix(0.25, nrow = length(statistic_cols), ncol = length(statistic_cols))
diag(sigma) <- 1
chol_sigma <- chol(sigma)
make_z <- function(n) {
  z <- matrix(stats::rnorm(n * length(statistic_cols)), nrow = n) %*% chol_sigma
  colnames(z) <- statistic_cols
  z
}

fit_z <- make_z(160L)
fit_table <- data.table(
  variant_id = sprintf("fit_%04d", seq_len(nrow(fit_z))),
  chr = 22L,
  pos_b37 = seq_len(nrow(fit_z)) * 100L
)
fit_table <- cbind(fit_table, as.data.table(fit_z))
write_gzip_table(fit_table, file.path(prepared_root, "fit_sample.tsv.gz"))

test_z <- make_z(180L)
test_table <- data.table(
  variant_id = sprintf("test_%04d", seq_len(nrow(test_z))),
  chr = 22L,
  pos_b37 = seq_len(nrow(test_z)) * 1000L
)
test_table <- cbind(test_table, as.data.table(test_z))
partition_name <- "chr22.tsv.gz"
write_gzip_table(test_table, file.path(prepared_root, partition_name))
fwrite(
  data.table(chr = 22L, path = partition_name, n_rows = nrow(test_table)),
  file.path(prepared_root, "chromosome_manifest.tsv"), sep = "\t"
)

copied <- file.copy(
  file.path(output_dir, "mediation", "mediation_ldl_liver_nearest.tsv.gz"),
  file.path(mediation_root, "mediation_ldl_liver_noshare.tsv.gz")
)
stopifnot(copied)

old_root <- getOption("csmgmm.reproduce.root")
old_rerun <- Sys.getenv("RERUN", unset = NA_character_)
options(csmgmm.reproduce.root = analysis_root)
Sys.setenv(RERUN = "1")
on.exit({
  options(csmgmm.reproduce.root = old_root)
  if (is.na(old_rerun)) Sys.unsetenv("RERUN") else Sys.setenv(RERUN = old_rerun)
}, add = TRUE)

real_out_root <- file.path(workspace, "realdata_contract")
real_out <- run_realdata(
  "glgc_pleiotropy", variant = "indep", out_root = real_out_root,
  cores = 1L, q = 0.1, ts = 1L, smoke = TRUE
)
real_required <- c(
  "fit.rds", "fit_runs.tsv", "sigma.tsv", "iqr_sd.tsv",
  "support_weights.tsv", "plot_sample.tsv.gz", "summary.tsv",
  "discoveries_t1.tsv.gz", "DONE"
)
stopifnot(all(file.exists(file.path(real_out, real_required))))
real_summary <- fread(file.path(real_out, "summary.tsv"))
stopifnot(
  nrow(real_summary) == 1L,
  real_summary$analysis == "glgc_pleiotropy",
  real_summary$target == "pleiotropy",
  real_summary$K == 5L,
  real_summary$t == 1L,
  isTRUE(real_summary$converged)
)

med_out_root <- file.path(workspace, "mediation_contract")
med_out <- run_mediation_ldl(
  "noshare", variant = "indep", out_root = med_out_root,
  cores = 1L, q = 0.1, ts = 1:2
)
med_required <- c(
  "fit.rds", "sigma.tsv", "iqr_sd.tsv", "support_weights.tsv",
  "summary.tsv", "discoveries.tsv.gz", "DONE"
)
stopifnot(all(file.exists(file.path(med_out, med_required))))
med_summary <- fread(file.path(med_out, "summary.tsv"))
stopifnot(
  nrow(med_summary) == 2L,
  med_summary$set[1L] == "noshare",
  identical(as.integer(med_summary$t), 1:2),
  all(med_summary$converged)
)

all_t_path <- run_mediation_all_t(
  1L, out_root = med_out_root, data_root = mediation_root,
  cores = 1L, q = 0.1
)
all_t <- fread(all_t_path)
stopifnot(nrow(all_t) == 5L, identical(as.integer(all_t$t), 1:5))

cat("REALDATA_OUTPUT_CONTRACTS_OK general_rows=", nrow(real_summary),
    " mediation_rows=", nrow(med_summary), " all_t_rows=", nrow(all_t), "\n", sep = "")
