# External support for the mediator--outcome link of the cis mediation discoveries, in the spirit of
# the held-out strata of the replication application. ILCCO-based S-PrediXcan statistics were used
# for fitting; here we look up the S-PrediXcan statistics for SCC in the Million Veteran Program
# (MVP), an independent GWAS, for lung and whole blood (the tissues available in MVP).
# For each discovered pair and tissue, the MVP statistic is signed by the fitted (ILCCO) statistic.
# Groups: rejected by both fits, by the fit with estimated Sigma only, by Sigma = I only, by neither.
# Inputs: Fig4/output/mediation_scc/cis15_{main,indep,robust}/fit.rds,
#         Data/mediation/mvp/mvp_scc__PM__Lung.csv and mvp_scc__PM__Blood.csv
#   Rscript Fig4/mediation_scc_mvp.R
here::i_am("Fig4/mediation_scc_mvp.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages({ library(data.table); library(csmGmmCorr) })
cols <- c("alpha_lung", "beta_lung", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver")
od <- here::here("Fig4", "output", "mediation_scc"); q <- 0.1
d <- fread(here::here("Data", "mediation", "mediation_scc_cis.tsv.gz"))
d <- d[complete.cases(d[, ..cols])]
Z <- capz(as.matrix(d[, ..cols]), 12)
for (v in c("main", "indep", "robust")) {
  f <- file.path(od, paste0("cis15_", v), "fit.rds"); if (!file.exists(f)) next
  fit <- readRDS(f)
  for (t in 1:2) {
    lf <- csm_lfdr(fit, Z, "mediation", t)
    d[, (sprintf("lfdr_%s_t%d", v, t)) := lf]
    d[, (sprintf("rej_%s_t%d", v, t)) := csm_stepup(lf, q)$rejected]
  }
}
read_mvp <- function(f) {
  x <- fread(f, select = c("gene_name", "zscore"))[!is.na(zscore)]
  x <- x[!(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE))]
  setnames(x, c("gene", "z"))[]
}
mvpdir <- Sys.getenv("MVP_DIR", here::here("Data", "mediation", "mvp"))
mvp <- list(lung = read_mvp(file.path(mvpdir, "mvp_scc__PM__Lung.csv")),
            blood = read_mvp(file.path(mvpdir, "mvp_scc__PM__Blood.csv")))
long <- rbindlist(lapply(names(mvp), function(tis) {
  x <- merge(d[, c("snp", "gene", paste0("beta_", tis), grep("^rej_", names(d), value = TRUE)), with = FALSE],
             mvp[[tis]], by = "gene")
  setnames(x, paste0("beta_", tis), "beta_fit")
  x[, tissue := tis][abs(beta_fit) > 0][, signed_mvp := z * sign(beta_fit)][]
}))
summ <- function(other) {
  out <- list()
  for (t in 1:2) {
    a <- long[[sprintf("rej_main_t%d", t)]]; b <- long[[sprintf("rej_%s_t%d", other, t)]]
    g <- fifelse(a & b, "both", fifelse(a, "estimated only", fifelse(b, paste(other, "only"), "neither")))
    out[[t]] <- long[, .(pairs = uniqueN(paste(snp, gene)), pair_tissues = .N,
                         same_sign = round(mean(signed_mvp > 0), 2), one_sided_05 = round(mean(signed_mvp > qnorm(0.95)), 2)),
                     by = .(group = g)][, `:=`(t = t, comparison = paste("estimated vs", other))][]
  }
  rbindlist(out)
}
res <- rbindlist(lapply(intersect(c("indep", "robust"), sub("rej_(.*)_t1", "\\1", grep("_t1$", grep("^rej_", names(d), value = TRUE), value = TRUE))), summ))
setcolorder(res, c("comparison", "t", "group"))
options(width = 200); print(res[order(comparison, t, group)])
fwrite(res, file.path(od, "cis15_mvp_support.tsv"), sep = "\t")
fwrite(long, file.path(od, "cis15_mvp_pairs.tsv"), sep = "\t")
