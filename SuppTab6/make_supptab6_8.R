# Supplementary Tables S26-S28 (folder numbering 6-8): design-level results of the main simulation (Table 6), with
# sign-symmetric statistics (Table 7) and of the correlation study (Table 8).
# Uses the output of Fig2 and Fig3:  Rscript SuppTab6/make_supptab6_8.R
here::i_am("SuppTab6/make_supptab6_8.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))

MNAME <- c(penalized = "Correlated csmGmm", penalized_adjusted = "Correlated csmGmm, adjusted",
           ML = "Correlated csmGmm (ML)", PCPathSafeguard = "JS-PC")
dtab <- function(D, caption, label, meths, file) {
  H <- "Target & $K$ & $t$ & Method & Design & $x$ & FDP (s.e.) & Power & $\\bar R$\\\\"
  L <- c("\\begin{longtable}{lcclrrrrr}", sprintf("\\caption{%s}\\label{%s}\\\\", caption, label), "\\toprule", H,
         "\\midrule", "\\endfirsthead", "\\toprule", H, "\\midrule", "\\endhead")
  D <- D[method %in% meths]
  D[, mo := match(method, meths)]
  setorder(D, target, K, t, mo, axis, x)
  for (i in seq_len(nrow(D))) {
    r <- D[i]; m <- if (r$method %in% names(MNAME)) MNAME[[r$method]] else r$method
    L <- c(L, sprintf("%s & %d & %d & %s & %s & %g & %.3f (%.3f) & %.3f & %.1f\\\\", TARGET_NAME[[r$target]], r$K, r$t,
                      tex_escape(m), tex_escape(r$design), r$x, r$FDP, r$FDP_se, r$power, r$R))
  }
  write_tex(c(L, "\\bottomrule", "\\end{longtable}"), file)
}
D <- design_means(read_sim("Fig2", "main"))
cm <- D[method %in% c("QCH", "adaFilter", "DEI", "MediationTest", "PCPathSafeguard")]
dtab(D, sprintf(paste0("Design-level results of the main simulation (500 replicates per design for the correlated csmGmm and ",
                       "its adjustment, at least %d for the other methods)."), min(c(D[method == "ML"]$n, cm$n))),
     "tab:sim-design", c("penalized", "penalized_adjusted", "ML", "QCH", "adaFilter", "DEI", "MediationTest", "PCPathSafeguard"),
     here::here("SuppTab6", "SuppTab6.tex"))
S <- design_means(read_sim("Fig2", "sym"))
dtab(S, paste0("Design-level results with symmetric statistics (logistic intercept 0; all twelve designs of each of the ",
               "three targets; 200 replicates per design)."), "tab:sim-sym", c("penalized", "ML"), here::here("SuppTab6", "SuppTab7.tex"))
R <- design_means(read_sim("Fig3", "rho"), by = c("target", "K", "t", "method", "design", "x"))
R[, axis := "rho"]
dtab(R, "Results of the correlation study (effect\\_05 setting of each target; 200 replicates per value of $\\rho$).",
     "tab:sim-rho", c("penalized", "ML"), here::here("SuppTab6", "SuppTab8.tex"))
