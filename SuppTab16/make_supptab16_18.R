# Supplementary Tables 16-18: incongruous results (averaged over designs, Table 16; by design,
# Table 17) and the cost of the monotone adjustment (Table 18), main simulation.
# Uses the output of Fig2:  Rscript SuppTab16/make_supptab16_18.R
here::i_am("SuppTab16/make_supptab16_18.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
TG <- c("pleiotropy", "replication", "mediation")
fm <- function(x) ifelse(is.na(x), "--", ifelse(abs(x) < 100, sprintf("%.2f", x), sprintf("%.1f", x)))
fs <- function(x) ifelse(is.na(x), "--", ifelse(x < 10, sprintf("%.2f", x), sprintf("%.1f", x)))
c1 <- c(pleiotropy = "QCH", mediation = "DEI"); c2 <- c(pleiotropy = "adaFilter", mediation = "MediationTest")

d <- read_sim("Fig2", "main")
D <- design_means(d); O <- overall_means(D)
og <- function(tg, K0, t0, m, col) { v <- O[target == tg & K == K0 & t == t0 & method == m][[col]]; if (length(v)) v else NA }
dg <- function(tg, K0, t0, m, de, col) { v <- D[target == tg & K == K0 & t == t0 & method == m & design == de][[col]]; if (length(v)) v else NA }
sh <- "\\shortstack{QCH\\\\or DEI} & \\shortstack{adaFilter or\\\\MediationTest}"

# Table 16
L <- c("\\begin{table}[!htbp]", "\\centering",
  paste0("\\caption{Incongruous results, averaged over the twelve designs of each target (500 replicates per design for the ",
    "correlated csmGmm, 200 for the other methods). Score level: the count of \\citet{Sun2025}, i.e.\\ the number of features ",
    "whose score exceeds that of a feature they dominate in the prespecified order, among the 5{,}000 smallest scores below 0.99 ",
    "(lfdr for the csmGmm fits and QCH; ranking $p$-value for the other methods). Decision level: the number of features not ",
    "rejected although they dominate a rejected feature. The monotone adjustment has no score-level incongruities by ",
    "construction. Comparators: QCH and adaFilter for pleiotropy, DEI and MediationTest for mediation (scores).}\\label{tab:incongruity-overall}"),
  "\\small", "\\begin{tabular}{lcc cc c c cc}", "\\toprule",
  " & & & \\multicolumn{2}{c}{Correlated csmGmm} & Adjusted & ML & \\multicolumn{2}{c}{Comparators}\\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-6}\\cmidrule(lr){7-7}\\cmidrule(lr){8-9}",
  paste0("Target & $K$ & $t$ & Score & Decision & Decision & Score & ", sh, "\\\\"), "\\midrule")
for (tg in TG) {
  first <- TRUE
  for (K0 in c(4, 6)) for (t0 in 1:2) {
    row <- c(if (first) TARGET_NAME[[tg]] else "", K0, t0, fs(og(tg, K0, t0, "penalized", "incon_score")),
             fm(og(tg, K0, t0, "penalized", "incon_dec")), fm(og(tg, K0, t0, "penalized_adjusted", "incon_dec")),
             fs(og(tg, K0, t0, "ML", "incon_score")),
             if (tg %in% names(c1)) fs(og(tg, K0, t0, c1[[tg]], "incon_score")) else "--",
             if (tg %in% names(c2)) fs(og(tg, K0, t0, c2[[tg]], "incon_score")) else "--")
    L <- c(L, paste0(paste(row, collapse = " & "), "\\\\")); first <- FALSE
  }
  L <- c(L, if (tg != "mediation") "\\midrule" else "\\bottomrule")
}
write_tex(c(L, "\\end{tabular}", "\\end{table}"), here::here("SuppTab16", "SuppTab16.tex"))

# Table 17
H <- paste0("Target & $K$ & $t$ & Design & Score & Decision & Decision & Score & ", sh, "\\\\")
L <- c("{\\footnotesize\\setlength{\\tabcolsep}{4pt}", "\\begin{longtable}{lccl rr r r rr}",
  paste0("\\caption{Incongruous results by design (means per replicate; 500 replicates per design for the correlated csmGmm, ",
    "200 for the other methods). Score level: count of \\citet{Sun2025} computed from unadjusted scores; decision level: ",
    "features not rejected although they dominate a rejected feature. Comparators: QCH and adaFilter for pleiotropy, DEI and ",
    "MediationTest for mediation (scores); JS-PC has zero counts throughout.}\\label{tab:incongruity-design}\\\\"),
  "\\toprule", " & & & & \\multicolumn{2}{c}{Correlated csmGmm} & Adjusted & ML & \\multicolumn{2}{c}{Comparators (score)}\\\\",
  H, "\\midrule", "\\endfirsthead", "\\toprule", H, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
for (tg in TG) for (K0 in c(4, 6)) for (t0 in 1:2) {
  des <- unique(D[target == tg & method == "penalized"]$design)
  des <- des[order(!startsWith(des, "effect"), des)]
  for (i in seq_along(des)) {
    de <- des[i]
    row <- c(if (i == 1) TARGET_NAME[[tg]] else "", if (i == 1) K0 else "", if (i == 1) t0 else "", tex_escape(de),
             fs(dg(tg, K0, t0, "penalized", de, "incon_score")), fm(dg(tg, K0, t0, "penalized", de, "incon_dec")),
             fm(dg(tg, K0, t0, "penalized_adjusted", de, "incon_dec")), fs(dg(tg, K0, t0, "ML", de, "incon_score")),
             if (tg %in% names(c1)) fs(dg(tg, K0, t0, c1[[tg]], de, "incon_score")) else "--",
             if (tg %in% names(c2)) fs(dg(tg, K0, t0, c2[[tg]], de, "incon_score")) else "--")
    L <- c(L, paste0(paste(row, collapse = " & "), "\\\\"))
  }
  if (!(tg == "mediation" && K0 == 6 && t0 == 2)) L <- c(L, "\\midrule")
}
write_tex(c(L, "\\end{longtable}}"), here::here("SuppTab16", "SuppTab17.tex"))

# Table 18: paired by replicate
key <- c("scen", "design", "K", "t", "rep", "target")
a <- d[method == "penalized", c(key, "R", "power", "FDP"), with = FALSE]
b <- d[method == "penalized_adjusted", c(key, "R", "power", "FDP"), with = FALSE]
m <- merge(a, b, by = key, suffixes = c("", "_a"))
m[, `:=`(lost = R - R_a, pl = power - power_a, dF = FDP_a - FDP)]
L <- c("\\begin{table}[!htbp]", "\\centering",
  sprintf(paste0("\\caption{Cost of the monotone adjustment for the correlated csmGmm in the main simulation (all twelve designs; ",
    "%s replicates per row). Columns: percentage of replicates with fewer rejections after adjustment, mean number of ",
    "rejections before adjustment, mean and maximum number of rejections removed, mean and maximum loss of power, and mean ",
    "change of FDP.}\\label{tab:adjustment-cost}"), fcomma(m[target == "pleiotropy" & K == 4 & t == 1, .N])),
  "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{lcc rrrr rr r}", "\\toprule",
  paste0("Target & $K$ & $t$ & \\% fewer & $\\bar R$ & \\shortstack{Mean\\\\lost} & \\shortstack{Max\\\\lost} & ",
         "\\shortstack{Mean power\\\\loss} & \\shortstack{Max power\\\\loss} & \\shortstack{Mean\\\\$\\Delta$FDP}\\\\"), "\\midrule")
for (tg in c("mediation", "pleiotropy", "replication")) {
  first <- TRUE
  for (K0 in c(4, 6)) for (t0 in 1:2) {
    s <- m[target == tg & K == K0 & t == t0]
    L <- c(L, paste0(paste(c(if (first) TARGET_NAME[[tg]] else "", K0, t0, sprintf("%.1f", 100 * mean(s$lost > 0)),
                             sprintf("%.1f", mean(s$R)), sprintf("%.3f", mean(s$lost)), max(s$lost),
                             sprintf("%.5f", mean(s$pl, na.rm = TRUE)), sprintf("%.3f", max(s$pl, na.rm = TRUE)),
                             sprintf("%+.4f", mean(s$dF))), collapse = " & "), "\\\\")); first <- FALSE
  }
}
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab16", "SuppTab18.tex"))
