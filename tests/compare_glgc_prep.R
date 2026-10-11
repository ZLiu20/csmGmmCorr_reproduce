args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 3L) args[3L] else ".")
source(file.path(root, "SupportingCode", "glgc_prepare.R"))

read_output <- function(path) {
  if (grepl("[.]gz$", path)) {
    data.table::fread(cmd = paste("gzip -cd", shQuote(path)))
  } else {
    data.table::fread(path)
  }
}

compare_prepared <- function(raw_root, prepared) {
  raw_root <- normalizePath(raw_root)
  prepared <- normalizePath(prepared)
  manifest <- if (grepl("replication", prepared)) {
    file.path(root, "Data", "glgc_replication_manifest.tsv")
  } else {
    file.path(root, "Data", "glgc_pleiotropy_manifest.tsv")
  }
  scratch <- tempfile("glgc_compare_")
  prepare_glgc(manifest, raw_root, scratch, overwrite = TRUE)
  files <- c("fit_sample.tsv.gz", "chromosome_manifest.tsv",
             "source_qc.tsv", "prepare_summary.tsv")
  for (name in files) {
    observed <- read_output(file.path(scratch, name))
    reference <- read_output(file.path(prepared, name))
    stopifnot(identical(names(observed), names(reference)),
              nrow(observed) == nrow(reference))
  }
  cat("GLGC_PREPARATION_EQUIVALENCE_OK\n")
}

fixture_test <- function() {
  workspace <- tempfile("glgc_fixture_")
  raw <- file.path(workspace, "raw")
  output1 <- file.path(workspace, "prepared_one")
  output2 <- file.path(workspace, "prepared_two")
  dir.create(raw, recursive = TRUE)

  common <- data.table::data.table(
    rsID = c("rs100", "rs100_lowN", "rs150", "rs12000", "rs22000",
             "bad_position", "rare", "small_sample", "source1_only"),
    CHROM = "1", POS_b37 = c(100L, 100L, 150L, 12000L, 22000L,
                              -1L, 30000L, 31000L, 40000L),
    REF = "A", ALT = "G",
    BETA = c(.2, .9, -.1, .4, -.2, .1, .2, .2, .3),
    SE = c(.1, .1, .05, .1, NA, .1, .1, .1, .1),
    P = c(.05, .01, .05, .001, .01, .5, .05, .05, .01),
    AF = c(.2, .2, .3, .4, .25, .2, .001, .2, .2),
    N = c(1000, 900, 1000, 1000, 1000, 1000, 1000, 400, 1000)
  )
  second <- data.table::data.table(
    rsID = c("other100", "other150", "other12000", "other22000", "source2_only"),
    CHROM = "1", POS_b37 = c(100L, 150L, 12000L, 22000L, 50000L),
    REF = "A", ALT = "G", BETA = c(-.3, .2, 9, .3, .4),
    SE = c(.1, .1, 1, .1, .1), P = c(.01, .05, 1e-10, .01, .01),
    AF = c(.2, .3, .4, .25, .2), N = 1000
  )
  data.table::fwrite(common, file.path(raw, "one.tsv.gz"), sep = "\t", na = "NA")
  data.table::fwrite(second, file.path(raw, "two.tsv.gz"), sep = "\t", na = "NA")
  manifest <- data.table::data.table(
    label = c("trait_one", "trait_two"),
    relative_path = c("one.tsv.gz", "two.tsv.gz"),
    beta_col = "BETA", se_col = "SE", p_col = "P",
    af_col = "AF", n_col = "N"
  )
  manifest_path <- file.path(workspace, "manifest.tsv")
  data.table::fwrite(manifest, manifest_path, sep = "\t")

  .glgc_hash <<- function(keys) sprintf("%016d", seq_along(keys))
  first <- prepare_glgc(manifest_path, raw, output1, z_cap = 5,
                        overwrite = TRUE)
  second_run <- prepare_glgc(manifest_path, raw, output2, z_cap = 5,
                             overwrite = TRUE)
  stopifnot(first$common == 4L, first$fit == 3L,
            identical(first, second_run))

  chromosome <- read_output(file.path(output1, "chr1.tsv.gz"))
  fit <- read_output(file.path(output1, "fit_sample.tsv.gz"))
  qc <- data.table::fread(file.path(output1, "source_qc.tsv"))
  summary <- data.table::fread(file.path(output1, "prepare_summary.tsv"))
  stopifnot(
    nrow(chromosome) == 4L,
    nrow(fit) == 3L,
    chromosome[pos_b37 == 100L, rsid] == "rs100",
    abs(chromosome[pos_b37 == 100L, trait_one] - 2) < 1e-10,
    chromosome[pos_b37 == 12000L, trait_two] == 5,
    chromosome[pos_b37 == 22000L, trait_one] < 0,
    qc[label == "trait_one", duplicate_key_rows] == 1L,
    qc[label == "trait_one", maf_fail] == 1L,
    qc[label == "trait_one", n_fail] == 1L,
    summary$strict_common == 4L,
    file.exists(file.path(output1, "PREPARATION_COMPLETE"))
  )
  for (name in c("fit_sample.tsv.gz", "chromosome_manifest.tsv",
                 "source_qc.tsv", "prepare_summary.tsv", "chr1.tsv.gz")) {
    left <- read_output(file.path(output1, name))
    right <- read_output(file.path(output2, name))
    stopifnot(identical(left, right))
  }
  cat("GLGC_PREPARATION_FIXTURE_OK\n")
}

if (length(args) >= 2L) compare_prepared(args[1L], args[2L]) else fixture_test()
