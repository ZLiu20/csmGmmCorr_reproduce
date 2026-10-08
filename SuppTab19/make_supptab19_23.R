# Supplementary Tables 19-23 (details of the applications), from the output of the Fig4 folder:
#   19 coordinates of the GLGC analyses (fixed), 20 all fits, 21 fitted and empirical
#   correlations, 22 sensitivity to a_mu, 23 replication in held-out strata
#   (Fig4/holdout_validation.R must have been run).
#   Rscript SuppTab19/make_supptab19_23.R
here::i_am("SuppTab19/make_supptab19_23.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
X <- rd_collect(); S <- X$summary; G <- X$sigma; FR <- X$fit_runs
AN <- c(ukb_pleiotropy = "Pleiotropy (UKB)", glgc_replication = "Replication (GLGC)", glgc_pleiotropy = "Pleiotropy (GLGC lipids)")
FIT <- c(main = "Estimated $\\bSigma$", robust = "Robust $\\bSigma$", indep = "$\\bSigma=I$", main_amu15 = "Estimated, $a_\\mu=1.5$",
         indep_amu15 = "$\\bSigma=I$, $a_\\mu=1.5$")

# Table 19
L <- c("\\begin{table}[htbp]", "\\centering", "\\caption{Coordinates of the GLGC analyses, in model order.}", "\\label{sup:tab:glgc-coordinates}",
  "\\small", "\\begin{tabular}{rll}", "\\toprule", "Coordinate & Pleiotropy ($K=5$) & Replication ($K=6$, LDL cholesterol)\\\\", "\\midrule")
pl <- c(HDL = "HDL cholesterol", LDL = "LDL cholesterol", TC = "Total cholesterol", logTG = "Log triglycerides", nonHDL = "Non-HDL cholesterol")
rp <- c(EUR = "European", EAS = "East Asian", AFR = "African or admixed African")
pc <- RD_ANALYSES$glgc_pleiotropy$cols; rc <- RD_ANALYSES$glgc_replication$cols
for (k in 1:6) L <- c(L, sprintf("%d & %s & %s, %s\\\\", k, if (k <= 5) pl[[pc[k]]] else "---", rp[[sub("_.*", "", rc[k])]], sub(".*_", "", rc[k])))
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab19", "SuppTab19.tex"))

# Table 20
H <- "Target & Fit & $t$ & $R$ & Adj. & Loci & Adj. loci & Score & Decision & Iter. & Min.\\\\"
L <- c("{\\footnotesize\\setlength{\\tabcolsep}{3.5pt}", "\\begin{longtable}{ll c rrrr rrr r}",
  paste0("\\caption{Genetic association analyses: all fits (UK Biobank pleiotropy, GLGC replication, and the GLGC lipid ",
    "pleiotropy analysis of Section~\\ref{sec:glgc-pleiotropy}). $R$: discoveries (variants); Adj.: after the monotone ",
    "adjustment; loci: greedy clumping in $\\pm500$~kb windows; Score and Decision: incongruity counts ",
    "(Section~\\ref{sec:incongruity}; score level among the 5{,}000 smallest unadjusted lfdr); Iter.: ECM iterations of the ",
    "selected start; Min.: fitting time in minutes with eight threads.}\\label{tab:realdata-all}\\\\"),
  "\\toprule", H, "\\midrule", "\\endfirsthead", "\\toprule", H, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
for (an in names(AN)) {
  first <- TRUE
  for (va in names(FIT)) for (t0 in 1:2) {
    r <- S[analysis == an & variant == va & t == t0]
    if (!nrow(r)) next
    L <- c(L, paste0(paste(c(if (first) AN[[an]] else "", if (t0 == 1) FIT[[va]] else "", t0, fcomma(r$R_raw), fcomma(r$R_adj),
                             fcomma(r$loci_raw), fcomma(r$loci_adj), fcomma(r$incon_score_raw), fcomma(r$incon_dec_raw),
                             if (t0 == 1) fcomma(r$iterations) else "", if (t0 == 1) sprintf("%.1f", r$fit_minutes) else ""),
                           collapse = " & "), "\\\\")); first <- FALSE
  }
  if (an != "glgc_pleiotropy") L <- c(L, "\\midrule")
}
write_tex(c(L, "\\end{longtable}}"), here::here("SuppTab19", "SuppTab20.tex"))

# Table 21
nm <- c(BMI = "BMI", SBP = "SBP", TG = "TG", HbA1c = "HbA1c", HDL = "HDL", LDL = "LDL", TC = "TC", logTG = "log TG", nonHDL = "non-HDL",
        EUR_female = "EUR F", EUR_male = "EUR M", EAS_female = "EAS F", EAS_male = "EAS M", AFR_female = "AFR F", AFR_male = "AFR M")
H <- "Target & Statistic 1 & Statistic 2 & Fitted & Pearson & $|z|<2$ & Robust\\\\"
L <- c("{\\footnotesize", "\\begin{longtable}{lll rrrr}",
  paste0("\\caption{Fitted correlation matrix (estimated $\\bSigma$, $a_\\mu=0.5$) and correlations of the fitting statistics: ",
    "Pearson, Pearson among features with all $|z|<2$, and the robust estimate used for the fit with fixed $\\bSigma$ ",
    "(Section~\\ref{sec:realdata-fits}).}\\label{tab:realdata-sigma}\\\\"),
  "\\toprule", H, "\\midrule", "\\endfirsthead", "\\toprule", H, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
for (an in names(AN)) {
  w <- dcast(G[analysis == an & variant == "main"], stat1 + stat2 ~ what, value.var = "value")
  cols <- RD_ANALYSES[[an]]$cols; w <- w[order(match(stat1, cols), match(stat2, cols))]
  for (i in seq_len(nrow(w))) L <- c(L, sprintf("%s & %s & %s & %.3f & %.3f & %.3f & %.3f\\\\", if (i == 1) AN[[an]] else "",
                                                nm[[w$stat1[i]]], nm[[w$stat2[i]]], w$fitted_Sigma[i], w$pearson[i], w$pearson_abs_lt2[i], w$robust[i]))
  if (an != "glgc_pleiotropy") L <- c(L, "\\midrule")
}
write_tex(c(L, "\\end{longtable}}"), here::here("SuppTab19", "SuppTab21.tex"))

# Table 22: log-likelihood of the fitting statistics (average log-likelihood of the selected start times n_fit)
ll <- function(an, va) {
  s <- S[analysis == an & variant == va][1]; f <- FR[analysis == an & variant == va]
  f[is.na(start) | start == "", start := "null"]
  f[start == s$selected_start]$loglik[1] * s$n_fit
}
L <- c("\\begin{table}[!htbp]", "\\centering",
  paste0("\\caption{Sensitivity to the lower bound $a_\\mu$ of the nonzero magnitudes (estimated $\\bSigma$, $t=1$). Discoveries are ",
    "loci. $\\Delta\\ell$: log-likelihood of the fitting statistics relative to $a_\\mu=0.5$ (whose feasible set contains those ",
    "of the larger bounds).}\\label{tab:realdata-amu}"),
  "\\small", "\\begin{tabular}{l rr rr}", "\\toprule", " & \\multicolumn{2}{c}{$a_\\mu=0.5$} & \\multicolumn{2}{c}{$a_\\mu=1.5$}\\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}", "Target & Loci & $\\Delta\\ell$ & Loci & $\\Delta\\ell$\\\\", "\\midrule")
for (an in names(AN)) {
  base <- ll(an, "main"); r15 <- S[analysis == an & variant == "main_amu15" & t == 1]
  L <- c(L, sprintf("%s & %s & 0 & %s & %s\\\\", AN[[an]], fcomma(S[analysis == an & variant == "main" & t == 1]$loci_raw),
                    if (nrow(r15)) fcomma(r15$loci_raw) else "--", if (nrow(r15)) fcomma(round(ll(an, "main_amu15") - base)) else "--"))
}
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab19", "SuppTab22.tex"))

# Table 23
hv <- fread(here::here("Fig4", "output", "holdout_validation.csv"))
lab <- c("Correlated csmGmm" = "Correlated csmGmm", "Sigma = I, near a correlated csmGmm discovery" = "$\\bSigma=I$, near a correlated csmGmm discovery",
         "Sigma = I only" = "$\\bSigma=I$ only", "Random variants" = "Random variants")
L <- c("\\begin{table}[!htbp]", "\\centering",
  paste0("\\caption{Replication in held-out strata. For the lead variant of each locus (greedy clumping, $\\pm500$~kb), the direction ",
    "$d$ is the sign of the sum of its six analysed statistics, and $S=d\\,(z_{\\mathrm{HIS},F}+z_{\\mathrm{HIS},M}+z_{\\mathrm{SAS},F}+",
    "z_{\\mathrm{SAS},M})/2$ combines the four Hispanic and South Asian strata, which were not used in the analysis. Columns: share ",
    "with $S>0$ (same direction), with $S>1.645$ and with $S>3.09$ (one-sided $p<0.05$ and $p<0.001$), and the mean number of ",
    "the four held-out strata with the same sign. ``$\\bSigma=I$ only'': loci of the csmGmm with $\\bSigma=I$ more than 500~kb from ",
    "every discovery of the correlated csmGmm. Random variants: 50{,}000 tested variants, with $d$ defined in the same way.}\\label{tab:realdata-holdout}"),
  "\\small", "\\begin{tabular}{lr cccc}", "\\toprule", "Loci & $n$ & $S>0$ & $S>1.645$ & $S>3.09$ & Strata with same sign\\\\", "\\midrule")
for (i in seq_len(nrow(hv))) L <- c(L, sprintf("%s & %s & %.2f & %.2f & %s & %.2f\\\\", lab[[hv$set[i]]], fcomma(hv$n[i]), hv$same_direction[i],
                                               hv$S_gt_1.645[i], if (hv$S_gt_3.09[i] < 0.01) sprintf("%.3f", hv$S_gt_3.09[i]) else sprintf("%.2f", hv$S_gt_3.09[i]),
                                               hv$strata_same_sign[i]))
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab19", "SuppTab23.tex"))
