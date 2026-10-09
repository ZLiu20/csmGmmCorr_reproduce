# Supplementary Tables 24-26: the three-mediator LDL mediation application
# (Supplementary Section "Mediation with summary statistics").
#   24 comparison of methods (t = 1, 2): discoveries, loci, discoveries with a strong pathway,
#      discoveries without any exposure-mediator evidence
#   25 top loci of the correlated csmGmm (a_mu = 3, t = 1) with the number of methods that reject
#      the lead variant, in the format of Table 1 of Sun, McCaw and Lin (2025)
#   26 diagnostics across fits: data sets (marginal, one variant per gene triple, conditional
#      alpha), Sigma, a_mu and t
# Requires Fig4/mediation_ldl_compare.R and Fig4/mediation_ldl_all_t.R to have been run.
#   Rscript SuppTab24/make_supptab24_26.R
here::i_am("SuppTab24/make_supptab24_26.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
od <- here::here("Fig4", "output", "mediation_ldl")
M <- fread(file.path(od, "compare", "methods.tsv"))
TOP <- fread(file.path(od, "compare", "top_loci_t1.tsv"))
d <- freadgz(here::here("Data", "mediation", "mediation_ldl_liver_noshare.tsv.gz"))
A <- abs(as.matrix(d[, MEDLDL_COLS, with = FALSE]))
n_strong <- sum((A[, 1] >= 4 & A[, 2] >= 4) | (A[, 3] >= 4 & A[, 4] >= 4) | (A[, 5] >= 4 & A[, 6] >= 4), na.rm = TRUE)
METH <- c("Correlated csmGmm (a_mu 3)" = "Correlated csmGmm, $a_\\mu=3$", "csmGmm, Sigma = I (a_mu 3)" = "csmGmm ($\\bSigma=I$), $a_\\mu=3$",
          "Correlated csmGmm (a_mu 0.5)" = "Correlated csmGmm, $a_\\mu=0.5$", "csmGmm, Sigma = I (a_mu 0.5)" = "csmGmm ($\\bSigma=I$), $a_\\mu=0.5$",
          DEI = "DEI", MediationTest = "MediationTest", "JS-PC" = "JS-PC")

# Table 24
L <- c("\\begin{table}[htbp]", "\\centering",
  paste0("\\caption{Mediation of LDL cholesterol by the liver expression of three nearby genes ($K=6$): comparison of methods ",
    "at $q=0.1$. $R$: discoveries (variants); loci: greedy clumping in $\\pm500$~kb windows; strong: discoveries with ",
    "$|\\alpha_g|\\ge4$ and $|\\beta_g|\\ge4$ for some $g$ (", fcomma(n_strong), " variants in the data); ",
    "$\\max_g|\\alpha_g|<2$: discoveries without exposure--mediator evidence.}"),
  "\\label{tab:med-ldl-methods}", "\\small", "\\begin{tabular}{l c rrrr}", "\\toprule",
  "Method & $t$ & $R$ & Loci & Strong & $\\max_g|\\alpha_g|<2$\\\\", "\\midrule")
for (t0 in 1:2) {
  for (m in names(METH)) {
    r <- M[t == t0 & method == m]; if (!nrow(r)) next
    L <- c(L, sprintf("%s & %d & %s & %s & %s & %s\\\\", METH[[m]], t0, fcomma(r$rejected), fcomma(r$loci),
                      fcomma(r$rejected_strong_pair), fcomma(r$rejected_all_alpha_lt2)))
  }
  if (t0 == 1) L <- c(L, "\\midrule")
}
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab24", "SuppTab24.tex"))

# Table 25
mcols <- intersect(names(METH), names(TOP)); nm <- length(mcols)
z1 <- function(x) sprintf("%.1f", x)
L <- c("\\begin{table}[htbp]", "\\centering",
  paste0("\\caption{Top 20 loci of the correlated csmGmm ($a_\\mu=3$, $t=1$), in the format of Table~1 of ",
    "\\citet{Sun2025}: lead variant, the three genes (by distance of the transcription start site), the $z$-statistics ",
    "of the three exposure--mediator ($\\alpha$) and mediator--outcome ($\\beta$) effects, the lfdr, and the number of the ",
    nm, " methods of Table~\\ref{tab:med-ldl-methods} that reject the lead variant at $q=0.1$ (numRej). ",
    "$^{*}$: strong pathway. $z$-statistics are capped at $\\pm12$.}"),
  "\\label{tab:med-ldl-top}", "\\scriptsize\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{r r l rr rr rr r c}", "\\toprule",
  "Rank & Chr & Genes $(M_1/M_2/M_3)$ & $\\alpha_1$ & $\\beta_1$ & $\\alpha_2$ & $\\beta_2$ & $\\alpha_3$ & $\\beta_3$ & lfdr & numRej\\\\",
  "\\midrule")
T20 <- head(TOP[order(rank)], 20)
T20[, numRej := rowSums(as.matrix(T20[, mcols, with = FALSE]))]
for (i in seq_len(nrow(T20))) {
  r <- T20[i]
  L <- c(L, sprintf("%d%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %d\\\\", r$rank, if (isTRUE(r$strong_pair)) "$^{*}$" else "",
    r$chr, tex_escape(paste(r$gene1, r$gene2, r$gene3, sep = "/")), z1(r$alpha1), z1(r$beta1), z1(r$alpha2), z1(r$beta2),
    z1(r$alpha3), z1(r$beta3), formatC(r$lfdr, format = "e", digits = 1), r$numRej))
}
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab24", "SuppTab25.tex"))

# Table 26
dirs <- list.dirs(od, recursive = FALSE)
X <- rbindlist(lapply(file.path(dirs, "all_t.tsv"), function(f) if (file.exists(f)) fread(f)), fill = TRUE)
FITS <- c(noshare_main = "Marginal $\\alpha$, estimated $\\bSigma$, $a_\\mu=0.5$",
          noshare_main_amu15 = "Marginal $\\alpha$, estimated $\\bSigma$, $a_\\mu=1.5$",
          noshare_main_amu30 = "Marginal $\\alpha$, estimated $\\bSigma$, $a_\\mu=3$",
          noshare_indep_amu30 = "Marginal $\\alpha$, $\\bSigma=I$, $a_\\mu=3$",
          noshare_robust_amu30 = "Marginal $\\alpha$, robust $\\bSigma$, $a_\\mu=3$",
          noshare_dedup_main_amu15 = "One variant per gene triple, $a_\\mu=1.5$",
          noshare_cond_main = "Conditional $\\alpha$, estimated $\\bSigma$, $a_\\mu=0.5$",
          noshare_cond_main_amu15 = "Conditional $\\alpha$, estimated $\\bSigma$, $a_\\mu=1.5$")
L <- c("{\\scriptsize\\setlength{\\tabcolsep}{3pt}", "\\begin{longtable}{l c rr rr rr}",
  paste0("\\caption{Mediation of LDL cholesterol: diagnostics of the csmGmm fits for $t=1,\\ldots,5$ ($t\\ge3$: the mediation ",
    "null equals the pleiotropy null). $n$: variants; $R$: discoveries; Strong: variants with a strong pathway and more than ",
    "$t$ coordinates with $|z|\\ge4$, in the data and among the discoveries; $\\max|\\alpha|<2$: discoveries without ",
    "exposure--mediator evidence; loci and loci containing a strong discovery.}\\label{tab:med-ldl-diag}\\\\"),
  "\\toprule", "Fit & $t$ & $R$ & Strong (data/disc.) & $\\max|\\alpha|<2$ & \\% & Loci & Strong loci\\\\", "\\midrule",
  "\\endfirsthead", "\\toprule", "Fit & $t$ & $R$ & Strong (data/disc.) & $\\max|\\alpha|<2$ & \\% & Loci & Strong loci\\\\",
  "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
for (f in names(FITS)) {
  x <- X[fit == f][order(t)]; if (!nrow(x)) next
  for (i in seq_len(nrow(x))) {
    r <- x[i]
    L <- c(L, sprintf("%s & %d & %s & %s/%s & %s & %s & %s & %s\\\\", if (i == 1) paste0(FITS[[f]], " ($n=", fcomma(r$n), "$)") else "",
      r$t, fcomma(r$R_raw), fcomma(r$strong_target_data), fcomma(r$strong_target_rej), fcomma(r$all_alpha_lt2),
      if (r$R_raw > 0) sprintf("%.0f", 100 * r$all_alpha_lt2 / r$R_raw) else "--", fcomma(r$loci), fcomma(r$loci_with_strong)))
  }
  L <- c(L, "\\midrule")
}
L[length(L)] <- "\\end{longtable}}"
write_tex(L, here::here("SuppTab24", "SuppTab26.tex"))
