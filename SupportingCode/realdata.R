# Applications (main-text Section 6; Supplementary Section 8).
#   glgc_replication  replication of LDL-cholesterol association across six GLGC ancestry-by-sex
#                     strata (K = 6); the Hispanic and South Asian strata are held out
#   ukb_pleiotropy    pleiotropy across five metabolic-syndrome traits, UK Biobank (K = 5)
#   glgc_pleiotropy   pleiotropy across five GLGC lipid traits (K = 5; Supplementary Section 8.6)
# Inputs: Data/prepared/<analysis>/ (Data/prepare_glgc.R, Data/prepare_ukb.R):
#   fit_sample.tsv.gz (one variant per occupied 10-kb window; used to fit the model),
#   chr<k>.tsv.gz and chromosome_manifest.tsv (all variants; the test set is every variant
#   not in the fitting sample).
RD_ANALYSES <- list(
  glgc_replication = list(target = "replication", prep = "glgc_replication",
                          cols = c("EUR_female", "EUR_male", "EAS_female", "EAS_male", "AFR_female", "AFR_male"),
                          holdout = c("HIS_female", "HIS_male", "SAS_female", "SAS_male")),
  ukb_pleiotropy   = list(target = "pleiotropy", prep = "ukb_pleiotropy",
                          cols = c("BMI", "SBP", "TG", "HDL", "HbA1c")),
  glgc_pleiotropy  = list(target = "pleiotropy", prep = "glgc_pleiotropy",
                          cols = c("HDL", "LDL", "TC", "logTG", "nonHDL")))

# The 13 fits of the manuscript: analysis, Sigma variant, lower bound a_mu of the magnitudes.
#   main   Sigma estimated (correlated csmGmm)
#   indep  Sigma fixed at I (csmGmm of Sun et al. with the same structure and penalty)
#   robust Sigma fixed at the projected robust correlation of the fitting statistics
RD_JOBS <- data.frame(
  analysis = c(rep("glgc_replication", 5), rep("ukb_pleiotropy", 4), rep("glgc_pleiotropy", 4)),
  variant  = c("main", "indep", "robust", "main", "indep", "main", "indep", "robust", "main",
               "main", "indep", "robust", "main"),
  a_mu     = c(0.5, 0.5, 0.5, 1.5, 1.5, 0.5, 0.5, 0.5, 1.5, 0.5, 0.5, 0.5, 1.5),
  stringsAsFactors = FALSE)
rd_label <- function(variant, a_mu) paste0(variant, if (a_mu != 0.5) sprintf("_amu%g", 10 * a_mu) else "")

freadgz <- function(f, ...) {
  if (grepl("\\.(gz|bgz)$", f)) data.table::fread(cmd = paste("gzip -cd", shQuote(f)), ...) else data.table::fread(f, ...)
}
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n"); flush.console() }
capz <- function(Z, cap) { storage.mode(Z) <- "double"; Z[Z > cap] <- cap; Z[Z < -cap] <- -cap; Z }

# Fitting sample and test set. extra: additional columns kept in info (e.g. held-out strata).
rd_load <- function(an, cap = 12, smoke = FALSE, extra = NULL) {
  a <- RD_ANALYSES[[an]]
  pdir <- here::here("Data", "prepared", a$prep)
  if (!file.exists(file.path(pdir, "chromosome_manifest.tsv"))) stop("prepared data not found: ", pdir)
  fs <- freadgz(file.path(pdir, "fit_sample.tsv.gz"), select = c("variant_id", a$cols))
  if (smoke) { set.seed(1); fs <- fs[sort(sample(.N, min(.N, 20000L)))] }
  Zfit <- capz(as.matrix(fs[, a$cols, with = FALSE]), cap)
  man <- data.table::fread(file.path(pdir, "chromosome_manifest.tsv"))
  if (smoke) man <- man[as.character(chr) %in% "22"]
  parts <- lapply(seq_len(nrow(man)), function(i) {
    p <- file.path(pdir, basename(man$path[i]))
    d <- freadgz(p, select = c("variant_id", "chr", "pos_b37", a$cols, extra))
    if (nrow(d) != man$n_rows[i]) stop("row count mismatch in ", p)
    d[!(variant_id %chin% fs$variant_id)]
  })
  d <- data.table::rbindlist(parts)
  Z <- capz(as.matrix(d[, a$cols, with = FALSE]), cap)
  info <- data.table::data.table(id = d$variant_id, chr = as.character(d$chr), pos = d$pos_b37)
  if (length(extra)) info <- cbind(info, d[, extra, with = FALSE])
  list(Zfit = Zfit, Z = Z, info = info, n_fit_excluded = sum(man$n_rows) - nrow(d))
}

# Reference correlations of the fitting statistics: Pearson; Pearson of the features with all
# |z| < 2; robust (IQR-based, Gnanadesikan-Kettenring) correlation.
robust_cor <- function(Z) {
  s <- function(v) (stats::quantile(v, 0.75, names = FALSE) - stats::quantile(v, 0.25, names = FALSE)) / 1.349
  K <- ncol(Z); R <- diag(K)
  for (i in seq_len(K - 1L)) for (j in (i + 1L):K) {
    a <- Z[, i] / s(Z[, i]); b <- Z[, j] / s(Z[, j])
    R[i, j] <- R[j, i] <- (s(a + b)^2 - s(a - b)^2) / (s(a + b)^2 + s(a - b)^2)
  }
  dimnames(R) <- list(colnames(Z), colnames(Z)); R
}
project_cor <- function(R, kappa = 0.05, margin = 1e-3) {
  for (it in 1:100) {
    e <- eigen((R + t(R)) / 2, symmetric = TRUE)
    if (min(e$values) >= kappa + margin / 2) break
    R <- e$vectors %*% diag(pmax(e$values, kappa + margin)) %*% t(e$vectors)
    R <- stats::cov2cor(R)
  }
  R <- (R + t(R)) / 2; diag(R) <- 1; R
}
reference_cors <- function(Z) {
  keep <- rowSums(abs(Z) < 2) == ncol(Z)
  list(pearson = stats::cor(Z), pearson_abs_lt2 = stats::cor(Z[keep, , drop = FALSE]), robust = robust_cor(Z),
       iqr_sd = apply(Z, 2, function(v) (stats::quantile(v, .75) - stats::quantile(v, .25)) / 1.349))
}

# Monotone adjustment and decision-level incongruity for millions of variants (compiled
# kernels, SupportingCode/src/rd_kernels.cpp; equal to csmGmmCorr::csm_adjust(all = TRUE) and
# incon_decision()).
kmode <- function(target) if (target == "replication") 1L else 0L
rd_majorant  <- function(Z, s, target, t) .Call("rd_majorant", Z, as.double(s), kmode(target), as.integer(t), 32L)
rd_incon_dec <- function(Z, rej, target, t) .Call("rd_incon_dec", Z, as.logical(rej), kmode(target), as.integer(t), 32L)

# Greedy distance clumping: leads in order of increasing score; a variant within win bp of a
# lead on the same chromosome joins that lead's locus. Returns the indices of the leads.
clump_leads <- function(chr, pos, score, sel, win = 5e5) {
  idx <- which(sel); if (!length(idx)) return(integer(0))
  idx <- idx[order(score[idx], pos[idx])]
  leads <- list(); keep <- integer(0)
  for (i in idx) {
    c0 <- chr[i]; v <- leads[[c0]]
    if (!is.null(v)) {
      k <- findInterval(pos[i], v)
      if ((k >= 1L && pos[i] - v[k] <= win) || (k < length(v) && v[k + 1L] - pos[i] <= win)) next
      leads[[c0]] <- append(v, pos[i], after = k)
    } else leads[[c0]] <- pos[i]
    keep <- c(keep, i)
  }
  keep
}
clump <- function(chr, pos, score, sel, win = 5e5) length(clump_leads(chr, pos, score, sel, win))

# One fit and its inference at t = 1, 2. Output in <out_dir>/<analysis>_<label>/:
#   fit.rds, fit_runs.tsv, sigma.tsv, iqr_sd.tsv, support_weights.tsv, summary.tsv,
#   discoveries_t<t>.tsv.gz (rejected, adjusted-rejected and incongruous variants with their
#   statistics), plot_sample.tsv.gz (200,000 random test variants), DONE.
run_realdata <- function(an, variant, a_mu = 0.5, out_root = here::here("Fig4", "output"),
                         cores = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "8")), q = 0.1, cap = 12,
                         ts = 1:2, smoke = Sys.getenv("SMOKE", "0") == "1") {
  suppressMessages({ library(data.table); library(csmGmmCorr) })
  load_kernels()
  stopifnot(an %in% names(RD_ANALYSES), variant %in% c("main", "indep", "robust"))
  a <- RD_ANALYSES[[an]]; target <- a$target; vlab <- rd_label(variant, a_mu)
  maxit <- if (smoke) 300L else 5000L
  out <- file.path(if (smoke) sub("output$", "output_smoke", out_root) else out_root, paste0(an, "_", vlab))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  say("analysis %s, variant %s, a_mu = %g, csmGmmCorr %s, cores %d, smoke %s", an, variant, a_mu,
      as.character(utils::packageVersion("csmGmmCorr")), cores, smoke)
  dat <- rd_load(an, cap, smoke)
  Zfit <- dat$Zfit; Z <- dat$Z; info <- dat$info; K <- ncol(Z)
  say("fit n = %d, test n = %d, K = %d", nrow(Zfit), nrow(Z), K)
  ref <- reference_cors(Zfit)
  fs <- switch(variant, main = FALSE, indep = TRUE, robust = project_cor(ref$robust, 0.05))
  spec <- csm_spec(K, M = csm_structure(K, 2L, 1L), a_mu = a_mu, b_mu = 15)
  ctl <- csm_control(max_iter = maxit, penalty = 30, weight_penalty = 10, fix_sigma = fs, cores = cores)
  t0 <- proc.time()[3]
  fit <- csm_fit(Zfit, spec, starts = c("null", "shared"), control = ctl)
  fit_min <- (proc.time()[3] - t0) / 60
  print(fit)
  saveRDS(fit, file.path(out, "fit.rds"))
  fwrite(cbind(analysis = an, variant = vlab, fit$runs), file.path(out, "fit_runs.tsv"), sep = "\t")
  Sig <- if (is.null(fit$params)) matrix(NA, K, K) else fit$params$Sigma
  dimnames(Sig) <- list(a$cols, a$cols)
  long <- function(M, what) { ij <- which(upper.tri(M), arr.ind = TRUE)
    data.table(analysis = an, variant = vlab, what = what, stat1 = a$cols[ij[, 1]], stat2 = a$cols[ij[, 2]], value = M[ij]) }
  fwrite(rbind(long(Sig, "fitted_Sigma"), long(ref$pearson, "pearson"), long(ref$pearson_abs_lt2, "pearson_abs_lt2"),
               long(ref$robust, "robust"), long(if (is.matrix(fs)) fs else diag(K), "fixed_Sigma_input")),
         file.path(out, "sigma.tsv"), sep = "\t")
  fwrite(data.table(analysis = an, stat = a$cols, iqr_sd = ref$iqr_sd), file.path(out, "iqr_sd.tsv"), sep = "\t")
  if (!is.null(fit$params)) {
    sp <- fit$spec; mu <- fit$params$mu; colnames(mu) <- paste0("mu_", a$cols)
    fwrite(data.table(analysis = an, variant = vlab, support = apply(sp$grp_A * 1L, 1, paste, collapse = ""),
                      size = sp$grp_size, component = sp$grp_m, weight = fit$params$varpi, mu),
           file.path(out, "support_weights.tsv"), sep = "\t")
  }
  set.seed(20261006); samp <- sort(sample(nrow(Z), min(nrow(Z), 200000L)))
  plot_dt <- data.table(info[samp], Z[samp, , drop = FALSE])
  rows <- list()
  for (t in ts) {
    if (!fit$converged) { say("fit not converged: no discoveries"); break }
    t1 <- proc.time()[3]
    lf <- csm_lfdr(fit, Z, target, t, cores = cores)
    su <- csm_stepup(lf, q)
    adj <- rd_majorant(Z, lf, target, t)
    if (any(adj < lf - 1e-12)) stop("adjusted lfdr below the unadjusted lfdr")
    sa <- csm_stepup(pmin(1, adj), q)
    dec_raw <- rd_incon_dec(Z, su$rejected, target, t); dec_adj <- rd_incon_dec(Z, sa$rejected, target, t)
    inc_raw <- incon_score(Z, lf, target, t, 5000L); inc_adj <- incon_score(Z, adj, target, t, 5000L)
    ex <- if (su$n_rejected > 0) csm_extrapolation(fit, Z, target, t, rejected = su$rejected, cutoff = su$cutoff)$flag else rep(NA, nrow(Z))
    rows[[length(rows) + 1L]] <- data.table(
      analysis = an, variant = vlab, a_mu = a_mu, target = target, K = K, t = t, q = q, n_fit = nrow(Zfit), n_test = nrow(Z),
      converged = fit$converged, selected_start = fit$selected, iterations = fit$runs$iterations[fit$runs$start == fit$selected],
      loglik = fit$loglik, fit_minutes = round(fit_min, 2),
      R_raw = su$n_rejected, R_adj = sa$n_rejected, R_both = sum(su$rejected & sa$rejected),
      frac_raw = su$n_rejected / nrow(Z), cutoff_raw = su$cutoff, avg_lfdr_raw = su$average,
      loci_raw = clump(info$chr, info$pos, lf, su$rejected), loci_adj = clump(info$chr, info$pos, adj, sa$rejected),
      n_adjusted_up = sum(adj > lf + 1e-12), incon_score_raw = inc_raw, incon_score_adj = inc_adj,
      incon_dec_raw = sum(dec_raw), incon_dec_adj = sum(dec_adj), extrap_flags = sum(ex, na.rm = TRUE),
      minutes_inference = round((proc.time()[3] - t1) / 60, 2))
    print(rows[[length(rows)]])
    sel <- su$rejected | sa$rejected | dec_raw
    disc <- data.table(info[sel], Z[sel, , drop = FALSE], lfdr = lf[sel], lfdr_adj = adj[sel],
                       rej_raw = su$rejected[sel], rej_adj = sa$rejected[sel], incon_dec = dec_raw[sel], extrap = ex[sel])
    fwrite(disc, file.path(out, sprintf("discoveries_t%d.tsv.gz", t)), sep = "\t")
    plot_dt[, (paste0("lfdr_t", t)) := lf[samp]]
    plot_dt[, (paste0("lfdr_adj_t", t)) := adj[samp]]
    rm(lf, adj, dec_raw, dec_adj, ex); gc(verbose = FALSE)
  }
  fwrite(plot_dt, file.path(out, "plot_sample.tsv.gz"), sep = "\t")
  if (length(rows)) fwrite(rbindlist(rows), file.path(out, "summary.tsv"), sep = "\t")
  writeLines("DONE", file.path(out, "DONE"))
  say("DONE %s %s", an, vlab)
  invisible(out)
}

# Combine the outputs of all fits (used by the plotting and table scripts).
rd_collect <- function(out_root = here::here("Fig4", "output")) {
  dirs <- list.dirs(out_root, recursive = FALSE)
  # analysis and fit label are taken from the folder name <analysis>_<label>
  rd <- function(f) data.table::rbindlist(lapply(dirs[file.exists(file.path(dirs, f))], function(d) {
    x <- data.table::fread(file.path(d, f)); b <- basename(d)
    an <- names(RD_ANALYSES)[startsWith(b, paste0(names(RD_ANALYSES), "_"))]
    if (length(an) != 1L) return(NULL)
    x[, analysis := an]
    if ("variant" %in% names(x)) x[, variant := substring(b, nchar(an) + 2L)]
    x
  }), fill = TRUE)
  list(summary = rd("summary.tsv"), sigma = rd("sigma.tsv"), fit_runs = rd("fit_runs.tsv"), iqr_sd = rd("iqr_sd.tsv"))
}

# ---- Mediation: lung squamous cell carcinoma, one analysis per sentinel SNP --------------------
# Design of Sun, McCaw and Lin (2025, Section 6): for each of 15 lung-cancer sentinel SNPs, the
# features are genes; alpha_g = SNP -> expression in tissue g (GTEx v8), beta_g = expression in
# tissue g -> squamous cell carcinoma (S-PrediXcan). Three tissues give K = 6. The primary tissue
# set is MED_TISSUES[[1]]; the others are sensitivity analyses.
MED_FILE <- function() here::here("Data", "mediation", "mediation_scc_bysnp.tsv.gz")
MED_SNPS <- c("rs71658797", "rs6920364", "rs11780471", "rs55781567", "rs56113850", "rs13080835", "rs7705526",
              "rs4236709", "rs885518", "rs11591710", "rs1056562", "rs77468143", "rs41309931", "rs116822326", "rs7953330")
MED_TISSUES <- list(blood_liver_adipose = c("blood", "liver", "adipose"),
                    blood_adipose_aorta = c("blood", "adipose", "aorta"),
                    blood_liver_aorta   = c("blood", "liver", "aorta"),
                    blood = "blood", liver = "liver", adipose = "adipose", aorta = "aorta", coronary = "coronary")
.med_grid <- function(tissues, a_mu) expand.grid(snp = MED_SNPS, variant = c("main", "indep", "robust"), tissues = tissues,
                                                 a_mu = a_mu, stringsAsFactors = FALSE)
MED_COMBOS <- names(MED_TISSUES)[1:3]; MED_SINGLE <- names(MED_TISSUES)[4:8]
# jobs 1-135: three tissues, a_mu = 0.5; 136-405: three tissues, a_mu = 1.5 and 2;
# 406-855: one tissue (K = 2, the design of Sun et al.), a_mu = 0.5 and 1.5
MED_JOBS <- rbind(.med_grid(MED_COMBOS, 0.5), .med_grid(MED_COMBOS, c(1.5, 2)), .med_grid(MED_SINGLE, c(0.5, 1.5)))

run_mediation <- function(snp, variant, tissues = names(MED_TISSUES)[1], a_mu = 0.5, out_root = here::here("Fig4", "output", "mediation"),
                          cores = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")), q = 0.1, cap = 12, ts = 1:2) {
  suppressMessages({ library(data.table); library(csmGmmCorr) })
  stopifnot(snp %in% MED_SNPS, variant %in% c("main", "indep", "robust"), tissues %in% names(MED_TISSUES))
  cols <- as.vector(outer(c("alpha_", "beta_"), MED_TISSUES[[tissues]], paste0))
  s0 <- snp; d <- freadgz(MED_FILE()); d <- d[d[["snp"]] == s0]
  d <- d[stats::complete.cases(d[, cols, with = FALSE])]
  Z <- capz(as.matrix(d[, cols, with = FALSE]), cap); K <- ncol(Z)
  vlab <- rd_label(variant, a_mu)
  out <- file.path(out_root, tissues, paste0(snp, "_", vlab)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  say("mediation %s (%s), %s, a_mu = %g, tissues %s: %d genes, K = %d", snp, d$snp_gene[1], variant, a_mu, tissues, nrow(Z), K)
  ref <- reference_cors(Z)
  fs <- switch(variant, main = FALSE, indep = TRUE, robust = project_cor(ref$robust, 0.05))
  spec <- csm_spec(K, M = csm_structure(K, 2L, 1L), a_mu = a_mu, b_mu = 15)
  ctl <- csm_control(max_iter = 5000L, penalty = 30, weight_penalty = 10, fix_sigma = fs, cores = cores)
  t0 <- proc.time()[3]
  fit <- csm_fit(Z, spec, starts = c("null", "shared"), control = ctl)
  fit_min <- (proc.time()[3] - t0) / 60
  variant <- vlab
  Sig <- if (is.null(fit$params)) matrix(NA, K, K) else fit$params$Sigma
  ij <- which(upper.tri(Sig), arr.ind = TRUE)
  fwrite(data.table(snp = snp, variant = variant, tissues = tissues, stat1 = cols[ij[, 1]], stat2 = cols[ij[, 2]],
                    fitted = Sig[ij], robust = ref$robust[ij], pearson = ref$pearson[ij]), file.path(out, "sigma.tsv"), sep = "\t")
  rows <- list(); disc <- list()
  ts <- ts[ts < K]                                   # K = 2: t = 1 only
  for (t in ts) {
    if (!fit$converged) break
    lf <- csm_lfdr(fit, Z, "mediation", t, cores = cores)
    su <- csm_stepup(lf, q)
    ad <- csm_adjust(Z, lf, "mediation", t, q)
    ex <- if (su$n_rejected > 0) csm_extrapolation(fit, Z, "mediation", t, rejected = su$rejected, cutoff = su$cutoff)$flag else rep(NA, nrow(Z))
    rows[[length(rows) + 1L]] <- data.table(snp = snp, snp_gene = d$snp_gene[1], variant = variant, a_mu = a_mu, tissues = tissues, t = t,
      n_genes = nrow(Z), converged = fit$converged, loglik = fit$loglik, fit_minutes = round(fit_min, 2),
      null_weight = sum(fit$params$varpi[fit$spec$grp_size == 0]),
      R_raw = su$n_rejected, R_adj = sum(ad$rejected), extrap_flags = sum(ex, na.rm = TRUE),
      genes = paste(d$gene[su$rejected][order(lf[su$rejected])], collapse = ","))
    if (any(su$rejected)) disc[[length(disc) + 1L]] <- data.table(t = t, gene = d$gene[su$rejected], Z[su$rejected, , drop = FALSE],
                                                                    lfdr = lf[su$rejected], extrap = ex[su$rejected])
  }
  if (length(rows)) fwrite(rbindlist(rows), file.path(out, "summary.tsv"), sep = "\t")
  if (length(disc)) fwrite(rbindlist(disc), file.path(out, "discoveries.tsv"), sep = "\t")
  writeLines("DONE", file.path(out, "DONE"))
  invisible(out)
}

# ---- Mediation with two mediators per feature (cis design; the model of main-text Section 2.1) ------
# Feature = gene (its cis genetic variation is the exposure); M1 = expression in whole blood,
# M2 = expression in liver; outcome = lung squamous cell carcinoma; K = 4, coordinates
# (alpha_blood, beta_blood, alpha_liver, beta_liver). Data/prepare_mediation_cis.R.
MEDCIS_JOBS <- expand.grid(variant = c("main", "indep", "robust"), a_mu = c(0.5, 1.5), stringsAsFactors = FALSE)
run_mediation_cis <- function(variant, a_mu = 0.5, out_root = here::here("Fig4", "output", "mediation_cis"),
                              cores = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")), q = 0.1, cap = 12, ts = 1:3) {
  suppressMessages({ library(data.table); library(csmGmmCorr) })
  stopifnot(variant %in% c("main", "indep", "robust"))
  cols <- c("alpha_blood", "beta_blood", "alpha_liver", "beta_liver")
  d <- freadgz(here::here("Data", "mediation", "mediation_cis_blood_liver.tsv.gz"))
  d <- d[stats::complete.cases(d[, cols, with = FALSE])]
  Z <- capz(as.matrix(d[, cols, with = FALSE]), cap); K <- ncol(Z)
  vlab <- rd_label(variant, a_mu)
  out <- file.path(out_root, vlab); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  say("cis mediation (blood, liver), %s: %d genes, K = %d", vlab, nrow(Z), K)
  ref <- reference_cors(Z)
  fs <- switch(variant, main = FALSE, indep = TRUE, robust = project_cor(ref$robust, 0.05))
  spec <- csm_spec(K, M = csm_structure(K, 2L, 1L), a_mu = a_mu, b_mu = 15)
  ctl <- csm_control(max_iter = 5000L, penalty = 30, weight_penalty = 10, fix_sigma = fs, cores = cores)
  t0 <- proc.time()[3]
  fit <- csm_fit(Z, spec, starts = c("null", "shared"), control = ctl)
  fit_min <- (proc.time()[3] - t0) / 60
  saveRDS(fit, file.path(out, "fit.rds"))
  Sig <- fit$params$Sigma; ij <- which(upper.tri(Sig), arr.ind = TRUE)
  fwrite(data.table(variant = vlab, stat1 = cols[ij[, 1]], stat2 = cols[ij[, 2]], fitted = Sig[ij], robust = ref$robust[ij],
                    pearson = ref$pearson[ij]), file.path(out, "sigma.tsv"), sep = "\t")
  sp <- fit$spec; mu <- fit$params$mu; colnames(mu) <- paste0("mu_", cols)
  fwrite(data.table(support = apply(sp$grp_A * 1L, 1, paste, collapse = ""), size = sp$grp_size, component = sp$grp_m,
                    weight = fit$params$varpi, mu), file.path(out, "support_weights.tsv"), sep = "\t")
  rows <- list(); disc <- list()
  for (t in ts) {
    if (!fit$converged) break
    lf <- csm_lfdr(fit, Z, "mediation", t, cores = cores)
    su <- csm_stepup(lf, q); ad <- csm_adjust(Z, lf, "mediation", t, q)
    ex <- if (su$n_rejected > 0) csm_extrapolation(fit, Z, "mediation", t, rejected = su$rejected, cutoff = su$cutoff)$flag else rep(NA, nrow(Z))
    rows[[length(rows) + 1L]] <- data.table(variant = vlab, a_mu = a_mu, t = t, n_genes = nrow(Z), converged = fit$converged,
      loglik = fit$loglik, fit_minutes = round(fit_min, 2), null_weight = sum(fit$params$varpi[sp$grp_size == 0]),
      R_raw = su$n_rejected, R_adj = sum(ad$rejected), extrap_flags = sum(ex, na.rm = TRUE))
    disc[[length(disc) + 1L]] <- data.table(t = t, d[su$rejected | ad$rejected, .(gene_id, gene_name, chr, start_b38)],
                                            Z[su$rejected | ad$rejected, , drop = FALSE], lfdr = lf[su$rejected | ad$rejected],
                                            rej = su$rejected[su$rejected | ad$rejected], rej_adj = ad$rejected[su$rejected | ad$rejected])
  }
  fwrite(rbindlist(rows), file.path(out, "summary.tsv"), sep = "\t")
  fwrite(rbindlist(disc), file.path(out, "discoveries.tsv"), sep = "\t")
  print(rbindlist(rows))
  writeLines("DONE", file.path(out, "DONE"))
  invisible(out)
}
