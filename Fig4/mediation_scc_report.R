# Lung SCC mediation application (main-text Section 6.3 and Table 3; Supplementary Section S5.7):
# numbers, tables and figures from the cis15 fits. Run after
#   Fig4/mediation_scc_fit.R (jobs 5, 6, 9-11, and jobs 5-6 with CAP=8.1 and CAP=6),
#   Fig4/mediation_scc_compare.R cis15 and Fig4/mediation_scc_mvp.R.
# Outputs in Fig4/output/mediation_scc/report/:
#   cis15_lfdr_all.tsv   lfdr of every cis pair for the three fits at t = 1, 2
#   overlap.tsv          pairs declared by each fit and shared by each pair of fits
#   methods.tsv          Supplementary Table S22 (pairs, strong pairs, weak-alpha pairs, incongruity)
#   sensitivity.tsv      Supplementary Table S23 (caps and a_mu)
#   mvp_support.tsv      Supplementary Table S25 and Figure 4I (MVP support, t = 1)
# and the figures Fig4/fig_mediation_lfdr.pdf (Supplementary Figure S13) and
# Fig4/fig_mediation_loci.pdf (Supplementary Figure S14).
#   Rscript Fig4/mediation_scc_report.R
here::i_am("Fig4/mediation_scc_report.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages({ library(data.table); library(csmGmmCorr); library(ggplot2); library(cowplot) })
MEDSCC_COLS <- c("alpha_lung", "beta_lung", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver")
q <- 0.1
od <- here::here("Fig4", "output", "mediation_scc"); rp <- file.path(od, "report"); dir.create(rp, showWarnings = FALSE)
d <- freadgz(here::here("Data", "mediation", "mediation_scc_cis.tsv.gz"))
d <- d[stats::complete.cases(d[, MEDSCC_COLS, with = FALSE])]
Z <- capz(as.matrix(d[, MEDSCC_COLS, with = FALSE]), 12); A <- abs(Z)
strong <- (A[, 1] >= 4 & A[, 2] >= 4) | (A[, 3] >= 4 & A[, 4] >= 4) | (A[, 5] >= 4 & A[, 6] >= 4)
weak_alpha <- pmax(A[, 1], A[, 3], A[, 5]) < 2
FITS <- c(main = "cis15_main", robust = "cis15_robust", indep = "cis15_indep")

# lfdr and rejections of the three fits
L <- d[, .(snp, locus_gene, chr, pos_b37, gene, tss_dist_kb)]
for (v in names(FITS)) {
  fit <- readRDS(file.path(od, FITS[[v]], "fit.rds"))
  for (t in 1:2) {
    lf <- csm_lfdr(fit, Z, "mediation", t)
    L[, (sprintf("lfdr_%s_t%d", v, t)) := lf]
    L[, (sprintf("rej_%s_t%d", v, t)) := csm_stepup(lf, q)$rejected]
  }
}
fwrite(L, file.path(rp, "cis15_lfdr_all.tsv"), sep = "\t")
rej <- function(v, t) L[[sprintf("rej_%s_t%d", v, t)]]
ov <- rbindlist(lapply(1:2, function(t) data.table(t = t,
  main = sum(rej("main", t)), robust = sum(rej("robust", t)), indep = sum(rej("indep", t)),
  main_and_robust = sum(rej("main", t) & rej("robust", t)), robust_and_indep = sum(rej("robust", t) & rej("indep", t)),
  main_and_indep = sum(rej("main", t) & rej("indep", t)),
  indep_same_pairs_t1_t2 = identical(rej("indep", 1), rej("indep", 2)))))
fwrite(ov, file.path(rp, "overlap.tsv"), sep = "\t"); print(ov)

# Table S22: five methods, with incongruity counts under the order of Section 5.2
COMP <- list(DEI = run_dei_mediation_R1, MediationTest = run_mediation_test_R1, "JS-PC" = run_mediation_safeguard_R1)
comparable <- function(Z) {     # pairs (i, j) in the same sign orthant with |z_i| <= |z_j| coordinatewise
  rp_ <- .incon_repr(Z, "mediation", 1L); W <- rp_$w; n <- 0L
  for (g in unique(rp_$group)) {
    ix <- which(rp_$group == g); if (length(ix) < 2L) next
    for (a in ix) for (b in ix) if (a < b && (all(W[a, ] <= W[b, ]) || all(W[b, ] <= W[a, ]))) n <- n + 1L
  }
  n
}
meth <- rbindlist(lapply(1:2, function(t) {
  res <- list()
  for (v in c("main", "indep")) {
    lf <- L[[sprintf("lfdr_%s_t%d", v, t)]]; r <- rej(v, t)
    res[[c(main = "Correlated csmGmm", indep = "csmGmm, Sigma = I")[[v]]]] <- list(r = r, s = lf)
  }
  for (m in names(COMP)) {
    o <- COMP[[m]](Z, t, q)
    res[[m]] <- list(r = as.logical(o$rejected), s = if (!is.null(o$score)) o$score else NULL)
  }
  rbindlist(lapply(names(res), function(m) data.table(t = t, method = m, pairs = sum(res[[m]]$r),
    strong = sum(res[[m]]$r & strong), weak_alpha = sum(res[[m]]$r & weak_alpha),
    incon_score = incon_score(Z, res[[m]]$s, "mediation", t), incon_dec = incon_decision(Z, res[[m]]$r, "mediation", t))))
}))
fwrite(meth, file.path(rp, "methods.tsv"), sep = "\t"); print(meth)
cat(sprintf("strong pairs: %d; comparable pairs under the absolute-value order: %d\n", sum(strong), comparable(Z)))

# Table S23: sensitivity (folders written by mediation_scc_fit.R)
sens_rows <- data.table(setting = c("Cap 12, a_mu 0.5", "Cap 12, a_mu 1.5", "Cap 8.1, a_mu 0.5", "Cap 6, a_mu 0.5"),
                        sfx_main = c("", "_mu15", "_cap8.1", "_cap6"))
rd_sum <- function(f) { p <- file.path(od, f, "summary.tsv"); if (file.exists(p)) fread(p) else NULL }
sens <- rbindlist(lapply(seq_len(nrow(sens_rows)), function(i) {
  s <- sens_rows$sfx_main[i]
  get <- function(v) { x <- rd_sum(paste0("cis15_", v, s)); if (is.null(x)) c(NA, NA) else x[order(t)]$R_raw }
  data.table(setting = sens_rows$setting[i], main_t1 = get("main")[1], main_t2 = get("main")[2],
             robust_t1 = get("robust")[1], robust_t2 = get("robust")[2], indep_t1 = get("indep")[1], indep_t2 = get("indep")[2],
             adj_changed = { x <- rd_sum(paste0("cis15_main", s)); if (is.null(x)) NA else any(x$R_adj != x$R_raw) },
             extrap_flags = { x <- rd_sum(paste0("cis15_main", s)); if (is.null(x)) NA else sum(x$extrap_flags) })
}))
fwrite(sens, file.path(rp, "sensitivity.tsv"), sep = "\t"); print(sens)

# Table S25 and Figure 4I: MVP support at t = 1 (pair-tissue statistics, Fig4/mediation_scc_mvp.R)
mvp <- fread(file.path(od, "cis15_mvp_pairs.tsv"))
mvp <- merge(mvp[, .(snp, gene, tissue, signed_mvp)], L[, .(snp, gene, m = rej_main_t1, i = rej_indep_t1)], by = c("snp", "gene"))
mvp[, group := fifelse(m & i, "both", fifelse(m, "correlated only", fifelse(i, "Sigma = I only", "neither")))]
ms <- mvp[, .(pairs = uniqueN(paste(snp, gene)), statistics = .N, same_sign = mean(signed_mvp > 0), gt_1.645 = mean(signed_mvp > qnorm(0.95))), by = group]
ms <- rbind(ms, mvp[m == TRUE, .(group = "correlated csmGmm (all)", pairs = uniqueN(paste(snp, gene)), statistics = .N,
                                 same_sign = mean(signed_mvp > 0), gt_1.645 = mean(signed_mvp > qnorm(0.95)))])
fwrite(ms, file.path(rp, "mvp_support.tsv"), sep = "\t"); print(ms)

# Figures S13 and S14
CAPL <- 20
L[, cat := factor(fifelse(rej_main_t1 & rej_indep_t1, "Both fits", fifelse(rej_main_t1, sprintf("Estimated %s only", SIGMA),
                  fifelse(rej_indep_t1, sprintf("%s = I only", SIGMA), "Neither"))),
                  levels = c("Both fits", sprintf("Estimated %s only", SIGMA), sprintf("%s = I only", SIGMA), "Neither"))]
L[, `:=`(y = pmin(-log10(pmax(lfdr_main_t1, 1e-300)), CAPL), x = pmin(-log10(pmax(lfdr_indep_t1, 1e-300)), CAPL))]
COL <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#b5b4ad"), levels(L$cat)); SHP <- setNames(c(16, 17, 15, 16), levels(L$cat))
cut_m <- -log10(max(L[rej_main_t1 == TRUE]$lfdr_main_t1)); cut_i <- -log10(max(L[rej_indep_t1 == TRUE]$lfdr_indep_t1))
sl <- scales::trans_new("slog", function(x) asinh(x), function(x) sinh(x))
brk <- c(0, 0.5, 1, 2, 5, 10, 20)
lab <- c("MPZL3", "TNF", "AGER", "AGPAT1", "C4A", "PSMA4", "HLA-C", "RAD52")
pA <- ggplot(L[order(cat != "Neither")], aes(x, y, color = cat, shape = cat)) +
  geom_abline(intercept = 0, slope = 1, color = "#d9d8d3") +
  geom_hline(yintercept = cut_m, linetype = 2, color = COL[[2]]) + geom_vline(xintercept = cut_i, linetype = 2, color = COL[[3]]) +
  geom_point(size = 2.4) +
  geom_text(data = L[gene %in% lab], aes(label = gene), color = "black", size = 2.6, hjust = -0.15, vjust = -0.3, show.legend = FALSE) +
  scale_color_manual(values = COL, labels = sprintf("%s (%d)", levels(L$cat), as.integer(table(L$cat))), name = NULL) +
  scale_shape_manual(values = SHP, labels = sprintf("%s (%d)", levels(L$cat), as.integer(table(L$cat))), name = NULL) +
  scale_x_continuous(trans = sl, breaks = brk) + scale_y_continuous(trans = sl, breaks = brk) +
  xlab(sprintf("-log10 lfdr, %s = I", SIGMA)) + ylab(sprintf("-log10 lfdr, estimated %s", SIGMA)) +
  theme_cowplot(11) + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
ggsave(here::here("Fig4", "fig_mediation_lfdr.pdf"), pA, width = 4.2, height = 4.0, device = PDF_DEVICE)
loci <- unique(L[, .(snp, locus_gene, chr, pos_b37)])[order(chr, pos_b37)]
L[, panel := factor(sprintf("%s\n%s, chr%d", snp, locus_gene, chr), levels = sprintf("%s\n%s, chr%d", loci$snp, loci$locus_gene, loci$chr))]
top <- L[(rej_main_t1 & lfdr_main_t1 < 1e-2) | (rej_indep_t1 & !rej_main_t1 & lfdr_indep_t1 < 1e-3)][order(lfdr_main_t1)][, head(.SD, 3), by = panel]
pB <- ggplot(L[order(cat != "Neither")], aes(tss_dist_kb / 1000, y, color = cat, shape = cat)) +
  geom_vline(xintercept = 0, color = "#d9d8d3") + geom_hline(yintercept = cut_m, linetype = 2, color = COL[[2]], linewidth = 0.4) +
  geom_point(size = 1.8) + geom_text(data = top, aes(label = gene), color = "black", size = 2, hjust = -0.15, show.legend = FALSE) +
  facet_wrap(~panel, nrow = 3) + scale_color_manual(values = COL, name = NULL) + scale_shape_manual(values = SHP, name = NULL) +
  coord_cartesian(xlim = c(-1.05, 1.05), ylim = c(-0.5, CAPL + 1)) +
  xlab("TSS minus SNP (Mb)") + ylab("-log10 lfdr") + theme_cowplot(9) +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(size = 7))
ggsave(here::here("Fig4", "fig_mediation_loci.pdf"), pB, width = 7.2, height = 5.4, device = PDF_DEVICE)
say("wrote %s and Fig4/fig_mediation_{lfdr,loci}.pdf", rp)
