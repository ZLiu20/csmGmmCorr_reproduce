# Table 1: simulation results averaged over the twelve designs of each target.
# Uses the output of the Fig2 folder:  Rscript Tab1/make_tab1.R   (writes Tab1/Tab1.tex)
here::i_am("Tab1/make_tab1.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))

D <- design_means(read_sim("Fig2", "main"))
O <- overall_means(D)
g <- function(tg, K0, t0, m, col) { v <- O[target == tg & K == K0 & t == t0 & method == m][[col]]; if (length(v)) v else NA }
P <- list("penalized", c("FDP", "FDPmax", "power"))
panels <- list(
  list("pleiotropy", "(a) Pleiotropy", list(P, list("QCH", c("FDP", "power")), list("adaFilter", c("FDP", "power"))), c("QCH", "adaFilter")),
  list("replication", "(b) Replication", list(P), character(0)),
  list("mediation", "(c) Mediation", list(P, list("DEI", c("FDP", "power")), list("MediationTest", c("FDP", "power")),
                                         list("PCPathSafeguard", c("FDP", "power"))), c("DEI", "MediationTest", "JS-PC")))
se_pen <- max(D[method == "penalized"]$FDP_se)
se_cmp <- max(D[method %in% c("QCH", "adaFilter", "DEI", "MediationTest", "PCPathSafeguard")]$FDP_se)
L <- c("\\begin{table}[!htbp]", "\\centering",
  sprintf(paste0("\\caption{Simulation results averaged over the twelve designs of each target ($q=0.10$): mean FDP, largest ",
    "design-level mean FDP (Max) and mean power of the correlated csmGmm (500 replicates per design), and mean FDP and power ",
    "of the comparison procedures (200 replicates; Section~\\ref{sec:sim-designs}). Monte Carlo standard errors of design-level ",
    "mean FDPs are at most %.3f for the correlated csmGmm and %.3f for the comparators.}\\label{tab:simulation-overview}"), se_pen, se_cmp),
  "\\small", "\\setlength{\\tabcolsep}{4pt}")
for (pn in panels) {
  tg <- pn[[1]]; spec <- pn[[3]]; nms <- pn[[4]]
  ncol <- 2 + sum(sapply(spec, function(s) length(s[[2]])))
  L <- c(L, sprintf("\\begin{tabular}{cc%s}", paste(sapply(spec, function(s) strrep("c", length(s[[2]]))), collapse = " ")),
         sprintf("\\multicolumn{%d}{l}{\\textit{%s}}\\\\", ncol, pn[[2]]), "\\toprule",
         paste0(" & & \\multicolumn{3}{c}{Correlated csmGmm}", paste0(sprintf(" & \\multicolumn{2}{c}{%s}", nms), collapse = ""), "\\\\"),
         paste0("\\cmidrule(lr){3-5}", paste0(sprintf("\\cmidrule(lr){%d-%d}", 6 + 2 * (seq_along(nms) - 1), 7 + 2 * (seq_along(nms) - 1)), collapse = "")))
  sub <- c("$K$", "$t$", unlist(lapply(spec, function(s) c(FDP = "FDP", FDPmax = "Max", power = "Power")[s[[2]]])))
  L <- c(L, paste0(paste(sub, collapse = " & "), "\\\\"), "\\midrule")
  for (K0 in c(4, 6)) for (t0 in 1:2)
    L <- c(L, paste0(paste(c(K0, t0, unlist(lapply(spec, function(s) sapply(s[[2]], function(cc) f3(g(tg, K0, t0, s[[1]], cc)))))),
                           collapse = " & "), "\\\\"))
  L <- c(L, "\\bottomrule", "\\end{tabular}", if (tg != "mediation") "\\par\\medskip")
}
L <- c(L, "\\end{table}")
write_tex(L, here::here("Tab1", "Tab1.tex"))
fwrite(O, here::here("Tab1", "Tab1_overall_means.csv"))
