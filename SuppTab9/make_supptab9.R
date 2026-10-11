# Supplementary Table 9: sensitivity to the number of magnitude components (structures A = 1/1,
# B = 2/1 (main), C = 2/2 for supports of size one / larger; replicates 1-50 of designs 1, 5, 7;
# B is taken from the Fig2 output).  Rscript SuppTab9/make_supptab9.R
here::i_am("SuppTab9/make_supptab9.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))

ms <- c("penalized_A", "penalized", "penalized_C")
d <- rbind(read_sim("SuppTab9", "main"), three_designs(read_sim("Fig2", "main"))[method == "penalized"], fill = TRUE)
d <- d[method %in% ms & rep <= 50]
g <- d[, .(F = mean(FDP), se = sd(FDP) / sqrt(.N), P = mean(power, na.rm = TRUE)), by = .(target, K, t, design, method)]
fwrite(g, here::here("SuppTab9", "SuppTab9.csv"))
w <- dcast(g, target + K + t + design ~ method, value.var = c("F", "se", "P")); setorder(w, target, K, t, design)
H <- "Target & $K$ & $t$ & Design & A & B & C & A & B & C\\\\"
L <- c("{\\footnotesize\\setlength{\\tabcolsep}{4pt}", "\\begin{longtable}{lccl ccc ccc}",
  paste0("\\caption{Sensitivity to the number of magnitude components (penalized fit; replicates 1--50 of three designs per ",
    "target). Structures A, B and C have one/one, two/one and two/two components for supports of size one/larger; B is used ",
    "in the main text. Entries are mean FDP (Monte Carlo s.e.) and mean power.}\\label{tab:structure-sensitivity}\\\\"),
  "\\toprule", " & & & & \\multicolumn{3}{c}{FDP} & \\multicolumn{3}{c}{Power}\\\\", "\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
  H, "\\midrule", "\\endfirsthead", "\\toprule", H, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
for (i in seq_len(nrow(w))) {
  r <- w[i]
  L <- c(L, paste0(paste(c(TARGET_NAME[[r$target]], r$K, r$t, tex_escape(r$design),
                           sapply(ms, function(m) sprintf("%.3f (%.3f)", r[[paste0("F_", m)]], r[[paste0("se_", m)]])),
                           sapply(ms, function(m) f3(r[[paste0("P_", m)]]))), collapse = " & "), "\\\\"))
}
write_tex(c(L, "\\end{longtable}}"), here::here("SuppTab9", "SuppTab9.tex"))
