# Supplementary Table 3: sensitivity to the penalty constants lambda_mu and lambda_w
# (replicates 1-50 of designs 1, 5, 7; the main setting lambda_mu = 30, lambda_w = 10 is taken
# from the Fig2 output, replicates 1-50).  Rscript SuppTab3/make_supptab3.R
here::i_am("SuppTab3/make_supptab3.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))

ps <- c("penalized_L0A10", "penalized_L10A10", "penalized", "penalized_L100A10", "penalized_L30A0", "penalized_L30A3", "penalized_L30A30")
d <- rbind(read_sim("SuppTab3", "main"), three_designs(read_sim("Fig2", "main"))[method == "penalized"], fill = TRUE)
d <- d[method %in% ps & rep <= 50]
g <- d[, .(F = mean(FDP), se = sd(FDP) / sqrt(.N), P = mean(power, na.rm = TRUE)), by = .(target, K, t, design, method)]
fwrite(g, here::here("SuppTab3", "SuppTab3.csv"))
pw <- dcast(g, target + K + t + design ~ method, value.var = "P")
pmax_ <- max(abs(as.matrix(pw[, ps, with = FALSE]) - pw$penalized), na.rm = TRUE)
semax <- max(g$se)
fd <- dcast(g, target + K + t + design ~ method, value.var = "F"); setorder(fd, target, K, t, design)
H <- "Target & $K$ & $t$ & Design & 0 & 10 & 30 & 100 & 0 & 3 & 30\\\\"
L <- c("{\\footnotesize\\setlength{\\tabcolsep}{4pt}", "\\begin{longtable}{lccl cccc ccc}",
  sprintf(paste0("\\caption{Sensitivity to the penalty constants (replicates 1--50 of three designs per target). Entries are mean ",
    "FDP. Left: $\\lambda_\\mu\\in\\{0,10,30,100\\}$ with $\\lambda_w=10$; right: $\\lambda_w\\in\\{0,3,30\\}$ with $\\lambda_\\mu=30$. ",
    "The column $\\lambda_\\mu=30$ is the setting used in the main text ($\\lambda_\\mu=30$, $\\lambda_w=10$). Monte Carlo standard ",
    "errors are at most %.3f; mean power differs from that of $\\lambda_\\mu=30,\\lambda_w=10$ by at most %.3f.}\\label{tab:penalty-choice}\\\\"),
    semax, pmax_),
  "\\toprule", " & & & & \\multicolumn{4}{c}{$\\lambda_\\mu$ ($\\lambda_w=10$)} & \\multicolumn{3}{c}{$\\lambda_w$ ($\\lambda_\\mu=30$)}\\\\",
  "\\cmidrule(lr){5-8}\\cmidrule(lr){9-11}", H, "\\midrule", "\\endfirsthead", "\\toprule", H, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
for (i in seq_len(nrow(fd))) {
  r <- fd[i]
  L <- c(L, paste0(paste(c(TARGET_NAME[[r$target]], r$K, r$t, tex_escape(r$design), sapply(ps, function(m) f3(r[[m]]))), collapse = " & "), "\\\\"))
}
write_tex(c(L, "\\end{longtable}}"), here::here("SuppTab3", "SuppTab3.tex"))
