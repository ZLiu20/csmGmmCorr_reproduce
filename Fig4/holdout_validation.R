# Replication in the held-out Hispanic and South Asian strata (main-text Section 6;
# Supplementary Section S5.4 and Table S18). For the lead variant of each locus (greedy clumping,
# +-500 kb), the direction d is the sign of the sum of its six analysed statistics and
# S = d (z_HIS,F + z_HIS,M + z_SAS,F + z_SAS,M) / 2 combines the four held-out strata.
#   Rscript Fig4/holdout_validation.R      (after the replication fits)
here::i_am("Fig4/holdout_validation.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
a <- RD_ANALYSES$glgc_replication; ho <- a$holdout
pdir <- here::here("Data", "prepared", a$prep)
od <- function(va) here::here("Fig4", "output", paste0("glgc_replication_", va))
disc <- function(va) { d <- fread(file.path(od(va), "discoveries_t1.tsv.gz")); d[rej_raw == TRUE] }
leads <- function(d) d[clump_leads(as.character(d$chr), d$pos, d$lfdr, rep(TRUE, nrow(d)))]
near <- function(x, y, win = 5e5) vapply(seq_len(nrow(x)), function(i) any(y$chr == x$chr[i] & abs(y$pos - x$pos[i]) <= win), TRUE)
dm <- disc("main"); di <- disc("indep")
lm <- leads(dm); li <- leads(di)
li[, near_main := near(li, dm)]
# all prepared variants with the analysed and held-out statistics; test set = not in the fitting sample
man <- fread(file.path(pdir, "chromosome_manifest.tsv"))
fs <- fread(cmd = paste("gzip -cd", shQuote(file.path(pdir, "fit_sample.tsv.gz"))), select = "variant_id")
all <- rbindlist(lapply(file.path(pdir, basename(man$path)), function(f)
  fread(cmd = paste("gzip -cd", shQuote(f)), select = c("variant_id", a$cols, ho))))
test <- all[!(variant_id %chin% fs$variant_id)]
set.seed(20261007); rnd <- test[sort(sample(.N, 50000L))]
score <- function(ids, lab) {
  x <- all[match(ids, variant_id)]
  d <- sign(rowSums(as.matrix(x[, a$cols, with = FALSE])))
  H <- as.matrix(x[, ho, with = FALSE]); S <- d * rowSums(H) / 2
  data.table(set = lab, n = length(ids), same_direction = mean(S > 0), S_gt_1.645 = mean(S > 1.645), S_gt_3.09 = mean(S > 3.09),
             strata_same_sign = mean(rowSums(sign(H) == d)))
}
res <- rbind(score(lm$id, "Correlated csmGmm"), score(li[near_main == TRUE]$id, "Sigma = I, near a correlated csmGmm discovery"),
             score(li[near_main == FALSE]$id, "Sigma = I only"), score(rnd$variant_id, "Random variants"))
print(res, digits = 3)
fwrite(res, here::here("Fig4", "output", "holdout_validation.csv"))
