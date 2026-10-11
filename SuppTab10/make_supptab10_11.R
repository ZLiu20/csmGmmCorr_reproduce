# Supplementary Tables 10-11: comparison with the csmGmm of Sun, McCaw and Lin (Sigma = I),
# the plug-in Sigma and the penalized fit with Sigma = I, in the main simulation (designs 1-9;
# arrays main_comp_*_csm of the Tab2 folder) and in the correlation study (rho_comp_*_csm of the
# Fig3 folder). The correlated csmGmm is taken from Fig2 and Fig3.
#   Rscript SuppTab10/make_supptab10_11.R
here::i_am("SuppTab10/make_supptab10_11.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
TG <- c("pleiotropy", "replication", "mediation")
cm <- c("csmGmm", "csmGmm_cor", "csmGmm_pen", "penalized")

d <- rbind(read_sim("Tab2", "main"), read_sim("Fig2", "main")[method == "penalized"], fill = TRUE)
d <- d[method %in% cm & !design %in% c("effect_06", "effect_07", "effect_08")]
gd <- d[, .(F = mean(FDP), P = mean(power, na.rm = TRUE), n = .N), by = .(target, K, t, method, design)]
go <- gd[, .(F = mean(F), Fx = max(F), P = mean(P)), by = .(target, K, t, method)]
fwrite(go, here::here("SuppTab10", "SuppTab10.csv"))
H <- "Target & $K$ & $t$ & FDP & Max & Power & FDP & Max & Power & FDP & Max & Power & FDP & Max & Power\\\\"
L <- c("\\begin{table}[!htbp]", "\\centering",
  sprintf(paste0("\\caption{Comparison with the csmGmm of \\citet{Sun2025} (designs 1--9, at least %d replicates per design, 500 ",
    "for the correlated csmGmm; averages over designs, ``Max'' the largest design-level mean FDP). csmGmm: $\\bSigma=I$, maximum ",
    "likelihood; plug-in: $\\bSigma$ fixed at the projected sample correlation of $\\mathbf Z$, maximum likelihood; csmGmm ",
    "(penalized): $\\bSigma=I$ with the penalty of Equation~\\eqref{main-eq:penalized}; correlated csmGmm: $\\bSigma$ estimated ",
    "with the penalty. All use the same magnitude structure, starting values and stopping rule.}\\label{tab:csm-comparison}"),
    min(gd[method != "penalized"]$n)),
  "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{lcc ccc ccc ccc ccc}", "\\toprule",
  " & & & \\multicolumn{3}{c}{csmGmm} & \\multicolumn{3}{c}{Plug-in $\\bSigma$} & \\multicolumn{3}{c}{csmGmm (penalized)} & \\multicolumn{3}{c}{Correlated csmGmm}\\\\",
  "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}\\cmidrule(lr){10-12}\\cmidrule(lr){13-15}", H, "\\midrule")
for (tg in TG) {
  first <- TRUE
  for (K0 in c(4, 6)) for (t0 in 1:2) {
    cells <- unlist(lapply(cm, function(m) { r <- go[target == tg & K == K0 & t == t0 & method == m]
      if (!nrow(r)) rep("--", 3) else c(f3(r$F), f3(r$Fx), f3(r$P)) }))
    L <- c(L, paste0(paste(c(if (first) TARGET_NAME[[tg]] else "", K0, t0, cells), collapse = " & "), "\\\\")); first <- FALSE
  }
  L <- c(L, if (tg != "mediation") "\\midrule" else "\\bottomrule")
}
write_tex(c(L, "\\end{tabular}", "\\end{table}"), here::here("SuppTab10", "SuppTab10.tex"))

d2 <- read_sim("Fig3", "rho")[method %in% cm]
g2 <- d2[, .(F = mean(FDP), n = .N), by = .(target, K, t, x, method)]
gg <- g2[, .(F = mean(F), Fx = max(F)), by = .(target, x, method)]
fwrite(gg, here::here("SuppTab10", "SuppTab11.csv"))
H <- "Target & $\\rho$ & Mean & Max & Mean & Max & Mean & Max & Mean & Max\\\\"
L <- c("\\begin{table}[!htbp]", "\\centering",
  sprintf(paste0("\\caption{Comparison with the csmGmm in the correlation study (fifth effect-size design of each target; FDP ",
    "averaged over $(K,t)\\in\\{4,6\\}\\times\\{1,2\\}$, ``Max'' the largest of the four; at least %d replicates per setting). ",
    "Fits as in Table~\\ref{tab:csm-comparison}.}\\label{tab:csm-rho}"), min(g2$n)),
  "\\footnotesize", "\\begin{tabular}{lc cc cc cc cc}", "\\toprule",
  " & & \\multicolumn{2}{c}{csmGmm} & \\multicolumn{2}{c}{Plug-in $\\bSigma$} & \\multicolumn{2}{c}{csmGmm (penalized)} & \\multicolumn{2}{c}{Correlated csmGmm}\\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}", H, "\\midrule")
for (tg in TG) {
  xs <- sort(unique(gg[target == tg]$x))
  for (i in seq_along(xs)) {
    cells <- unlist(lapply(cm, function(m) { r <- gg[target == tg & x == xs[i] & method == m]
      if (!nrow(r)) c("--", "--") else c(f3(r$F), f3(r$Fx)) }))
    L <- c(L, paste0(paste(c(if (i == 1) TARGET_NAME[[tg]] else "", format(xs[i]), cells), collapse = " & "), "\\\\"))
  }
  L <- c(L, if (tg != "mediation") "\\midrule" else "\\bottomrule")
}
write_tex(c(L, "\\end{tabular}", "\\end{table}"), here::here("SuppTab10", "SuppTab11.tex"))
