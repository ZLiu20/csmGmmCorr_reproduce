# Table 3 and the numbers reported in main-text Section 6.
#   Rscript Fig4/make_tab3.R      (after all 13 fits of Fig4/run_realdata.lsf)
here::i_am("Fig4/make_tab3.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
X <- rd_collect(); S <- X$summary; X$iqr_sd <- unique(X$iqr_sd[, .(analysis, stat, iqr_sd)])
fwrite(S, here::here("Fig4", "output", "summary_all.tsv"), sep = "\t")
fwrite(X$sigma, here::here("Fig4", "output", "sigma_all.tsv"), sep = "\t")
fwrite(X$fit_runs, here::here("Fig4", "output", "fit_runs_all.tsv"), sep = "\t")
fwrite(X$iqr_sd, here::here("Fig4", "output", "iqr_sd_all.tsv"), sep = "\t")
v <- function(an, va, t0, col) { x <- S[analysis == an & variant == va & t == t0][[col]]; if (length(x)) x[1] else NA }
L <- c("\\begin{table}[!htbp]", "\\centering",
  paste0("\\caption{Discoveries at $q=0.10$: variants and loci of the correlated csmGmm, loci after its monotone adjustment, ",
         "and loci of the fits with $\\bSigma$ fixed at a robust correlation estimate and at $I$ (csmGmm). Loci: greedy clumping ",
         "in $\\pm500$~kb windows.}\\label{tab:realdata}"),
  "\\small", "\\begin{tabular}{lc rr r rr}", "\\toprule",
  " & & \\multicolumn{3}{c}{Correlated csmGmm} & Robust $\\bSigma$ & $\\bSigma=I$\\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-6}\\cmidrule(lr){7-7}", "Analysis & $t$ & Variants & Loci & Adjusted loci & Loci & Loci\\\\", "\\midrule")
for (an in c("glgc_replication", "ukb_pleiotropy")) for (t0 in 1:2) {
  lab <- if (t0 == 1) c(glgc_replication = "Replication, GLGC", ukb_pleiotropy = "Pleiotropy, UK Biobank")[[an]] else
    sprintf("\\quad %s tested", fcomma(v(an, "main", 1, "n_test")))
  L <- c(L, sprintf("%s & %d & %s & %s & %s & %s & %s\\\\", lab, t0, fcomma(v(an, "main", t0, "R_raw")), fcomma(v(an, "main", t0, "loci_raw")),
                    fcomma(v(an, "main", t0, "loci_adj")), fcomma(v(an, "robust", t0, "loci_raw")), fcomma(v(an, "indep", t0, "loci_raw"))))
}
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("Fig4", "Tab3.tex"))

# ---- numbers in the text ----
rdisc <- function(an, va, t0) fread(here::here("Fig4", "output", paste0(an, "_", va), sprintf("discoveries_t%d.tsv.gz", t0)))
leads_of <- function(d) { d <- d[rej_raw == TRUE]; d[clump_leads(as.character(d$chr), d$pos, d$lfdr, rep(TRUE, nrow(d)))] }
near <- function(a, b, win = 5e5) {   # is each variant of a within win of some variant of b (same chromosome)?
  out <- logical(nrow(a))
  for (ch in unique(a$chr)) {
    i <- which(a$chr == ch); p <- sort(b$pos[b$chr == ch])
    if (!length(p)) next
    k <- findInterval(a$pos[i], p)
    dl <- ifelse(k >= 1, a$pos[i] - p[pmax(k, 1)], Inf); dr <- ifelse(k < length(p), p[pmin(k + 1, length(p))] - a$pos[i], Inf)
    out[i] <- pmin(dl, dr) <= win
  }
  out
}
tr <- RD_ANALYSES$ukb_pleiotropy$cols
m2 <- rdisc("ukb_pleiotropy", "main", 2); i2 <- rdisc("ukb_pleiotropy", "indep", 2)
lm <- leads_of(m2); li <- leads_of(i2)
only <- li[!near(li, m2[rej_raw == TRUE])]
mx <- function(d) apply(abs(as.matrix(d[, tr, with = FALSE])), 1, max)
nb <- function(d) rowSums(abs(as.matrix(d[, tr, with = FALSE])) > 3.29)
out <- c(sprintf("UKB pleiotropy, t = 2: Sigma = I only loci: %d; median largest |z| of their leads %.1f vs %.1f for the correlated csmGmm loci; median number of traits with |z| > 3.29 among their leads: %g",
                 nrow(only), median(mx(only)), median(mx(lm)), median(nb(only))),
         sprintf("UKB fitted weight of the all-zero configuration (estimated Sigma): %s",
                 { w <- fread(here::here("Fig4", "output", "ukb_pleiotropy_main", "support_weights.tsv")); sprintf("%.3f", sum(w[size == 0]$weight)) }),
         sprintf("Extrapolation flags (estimated Sigma): UKB t = 1: %d, t = 2: %d of %s", v("ukb_pleiotropy", "main", 1, "extrap_flags"),
                 v("ukb_pleiotropy", "main", 2, "extrap_flags"), fcomma(v("ukb_pleiotropy", "main", 2, "R_raw"))),
         sprintf("Monotone adjustment, t = 1: replication keeps %.0f%% of variants and %.0f%% of loci; pleiotropy %.0f%% and %.0f%%",
                 100 * v("glgc_replication", "main", 1, "R_adj") / v("glgc_replication", "main", 1, "R_raw"),
                 100 * v("glgc_replication", "main", 1, "loci_adj") / v("glgc_replication", "main", 1, "loci_raw"),
                 100 * v("ukb_pleiotropy", "main", 1, "R_adj") / v("ukb_pleiotropy", "main", 1, "R_raw"),
                 100 * v("ukb_pleiotropy", "main", 1, "loci_adj") / v("ukb_pleiotropy", "main", 1, "loci_raw")),
         sprintf("Replication loci, a_mu = 1.5: estimated Sigma %s, Sigma = I %s", fcomma(v("glgc_replication", "main_amu15", 1, "loci_raw")),
                 fcomma(v("glgc_replication", "indep_amu15", 1, "loci_raw"))),
         sprintf("IQR-based SD of the fitting statistics: %s",
                 paste(sprintf("%s %.2f", X$iqr_sd$stat, X$iqr_sd$iqr_sd), collapse = "; ")))
writeLines(out); writeLines(out, here::here("Fig4", "Fig4_numbers.txt"))
