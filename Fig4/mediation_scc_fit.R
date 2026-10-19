# Correlated csmGmm for the lung SCC mediation data (Data/prepare_mediation_scc.R).
# Feature = SNP-gene pair; G = 3 tissue pathways (lung, whole blood, liver); K = 6.
# A pair is a discovery when at least one tissue pathway has both links nonzero (mediation, t = 1);
# t = 2 is also reported. Settings as in the main text (structure B, a_mu = 0.5, b_mu = 15,
# lambda_mu = 30, lambda_w = 10, two starts, q = 0.1, statistics capped at +-12).
# Jobs: pair set x fit (estimated Sigma, Sigma = I). Pair sets: all = every gene for the 15 SNPs;
# top3 = every gene for the 3 SNPs ranked highest in ILCCO; cis15 / cis3 = only genes whose TSS is
# within 1 Mb of the SNP (Data/prepare_mediation_scc_cis.R), for the 15 SNPs / the top 3 SNPs.
# Jobs 9-11 (cis15 only), as in the replication and pleiotropy applications: Sigma fixed at the
# projected robust correlation, and a_mu = 1.5 with estimated Sigma and with Sigma = I.
#   Rscript Fig4/mediation_scc_fit.R <job>      (jobs 1-4: all, top3; 5-8: cis15, cis3; 9-11: cis15 sensitivity)
here::i_am("Fig4/mediation_scc_fit.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages({ library(data.table); library(csmGmmCorr) })
MEDSCC_COLS <- c("alpha_lung", "beta_lung", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver")
JOBS <- rbind(expand.grid(variant = c("main", "indep"), snps = c("all", "top3", "cis15", "cis3"), stringsAsFactors = FALSE),
              data.frame(variant = c("robust", "main_mu15", "indep_mu15"), snps = "cis15"))
job <- as.integer(commandArgs(trailingOnly = TRUE)[1]); stopifnot(job %in% seq_len(nrow(JOBS)))
variant <- JOBS$variant[job]; snpset <- JOBS$snps[job]
cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")); q <- 0.1; ts <- 1:2
a_mu <- if (grepl("_mu15$", variant)) 1.5 else 0.5

dd <- here::here("Data", "mediation")
d <- freadgz(file.path(dd, if (startsWith(snpset, "cis")) "mediation_scc_cis.tsv.gz" else "mediation_scc.tsv.gz"))
if (snpset %in% c("top3", "cis3")) {
  top <- head(fread(file.path(dd, "mediation_scc_snp_rank.tsv"))[order(rank_ILCCO)]$snp, 3)
  d <- d[snp %in% top]
}
d <- d[stats::complete.cases(d[, MEDSCC_COLS, with = FALSE])]
CAP <- as.numeric(Sys.getenv("CAP", "12"))      # 8.1 = the cap of the mediation analysis of Sun, McCaw and Lin (2025)
Z <- capz(as.matrix(d[, MEDSCC_COLS, with = FALSE]), CAP); K <- ncol(Z)
vlab <- paste0(snpset, "_", variant, if (CAP != 12) paste0("_cap", CAP) else "")
out <- here::here("Fig4", "output", "mediation_scc", vlab); dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, "DONE")) && Sys.getenv("RERUN", "0") != "1") { say("already done: %s", out); quit(save = "no") }
say("mediation SCC, %s: %d SNP-gene pairs, %d SNPs, %d genes, K = %d", vlab, nrow(Z), uniqueN(d$snp), uniqueN(d$gene), K)

ref <- reference_cors(Z)
fs <- switch(sub("_mu15$", "", variant), main = FALSE, indep = TRUE, robust = project_cor(ref$robust, 0.05))
spec <- csm_spec(K, M = csm_structure(K, 2L, 1L), a_mu = a_mu, b_mu = 15)
ctl <- csm_control(max_iter = 5000L, penalty = 30, weight_penalty = 10, fix_sigma = fs, cores = cores)
t0 <- proc.time()[3]
fit <- csm_fit(Z, spec, starts = c("null", "shared"), control = ctl)
fit_min <- (proc.time()[3] - t0) / 60
print(fit)
saveRDS(fit, file.path(out, "fit.rds"))
Sig <- fit$params$Sigma; ij <- which(upper.tri(Sig), arr.ind = TRUE)
fwrite(data.table(fit = vlab, stat1 = MEDSCC_COLS[ij[, 1]], stat2 = MEDSCC_COLS[ij[, 2]], fitted = Sig[ij],
                  robust = ref$robust[ij], pearson = ref$pearson[ij]), file.path(out, "sigma.tsv"), sep = "\t")
fwrite(data.table(fit = vlab, stat = MEDSCC_COLS, iqr_sd = ref$iqr_sd), file.path(out, "iqr_sd.tsv"), sep = "\t")

rows <- list(); disc <- list()
for (t in ts) {
  if (!fit$converged) { say("fit not converged"); break }
  lf <- csm_lfdr(fit, Z, "mediation", t, cores = cores)
  su <- csm_stepup(lf, q); ad <- csm_adjust(Z, lf, "mediation", t, q)
  ex <- if (su$n_rejected > 0) csm_extrapolation(fit, Z, "mediation", t, rejected = su$rejected, cutoff = su$cutoff)$flag else rep(NA, nrow(Z))
  r <- su$rejected
  rows[[length(rows) + 1L]] <- data.table(fit = vlab, t = t, n_pairs = nrow(Z), n_snps = uniqueN(d$snp), n_genes = uniqueN(d$gene),
    converged = fit$converged, loglik = fit$loglik, fit_minutes = round(fit_min, 2),
    null_weight = sum(fit$params$varpi[fit$spec$grp_size == 0]),
    R_raw = sum(r), R_adj = sum(ad$rejected), genes_raw = uniqueN(d$gene[r]), snps_raw = uniqueN(d$snp[r]),
    extrap_flags = sum(ex, na.rm = TRUE))
  s <- r | ad$rejected
  disc[[length(disc) + 1L]] <- data.table(t = t, d[s, .(snp, locus_gene, gene)], round(Z[s, , drop = FALSE], 3),
                                          lfdr = lf[s], rej = r[s], rej_adj = ad$rejected[s], extrap = ex[s])
}
if (length(rows)) { fwrite(rbindlist(rows), file.path(out, "summary.tsv"), sep = "\t"); print(rbindlist(rows)) }
if (length(disc)) fwrite(rbindlist(disc), file.path(out, "discoveries.tsv.gz"), sep = "\t")
writeLines("DONE", file.path(out, "DONE"))
