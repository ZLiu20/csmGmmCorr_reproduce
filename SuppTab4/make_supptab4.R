# Supplementary Table 4: mean FDP (Monte Carlo s.e.) and mean power over 10 replicates for the
# oracle posterior and for exchangeable-restriction fits with floors 1/64 and 1e-4.
#   Rscript SuppTab4/make_supptab4.R
here::i_am("SuppTab4/make_supptab4.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
d <- rbindlist(lapply(list.files(here::here("SuppTab4", "output"), pattern = "^floor_rep.*csv$", full.names = TRUE), fread))
g <- d[, .(FDP = mean(FDP), se = sd(FDP) / sqrt(.N), power = mean(power), mass2to4 = mean(omega2to4), n = .N), by = .(fit, target)]
print(g); fwrite(g, here::here("SuppTab4", "SuppTab4.csv"))
cell <- function(tg, f) { r <- g[target == tg & fit == f]; sprintf("%.3f (%.3f), %.3f", r$FDP, r$se, r$power) }
L <- c("\\begin{table}[htbp]", "\\centering",
  sprintf(paste0("\\caption{Mean FDP (Monte Carlo standard error) and mean power over %d replicates for the oracle posterior ",
    "and for exchangeable-restriction fits with two floors. The design is described in the text.}\\label{tab:floor-bias}"), max(g$n)),
  "\\begin{tabular}{@{}lccc@{}}", "\\toprule", "Target & Oracle & $a_\\omega=1/64$ & $a_\\omega=10^{-4}$\\\\", "\\midrule")
for (tg in c("pleiotropy", "replication", "mediation"))
  L <- c(L, sprintf("%s & %s & %s & %s\\\\", TARGET_NAME[[tg]], cell(tg, "oracle"), cell(tg, "floor_0.015625"), cell(tg, "floor_1e-04")))
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab4", "SuppTab4.tex"))
cat(sprintf("fitted mass on sizes 2-4 with floor 1/64: %s (3/64 = %.4f)\n",
            paste(sprintf("%.4f", range(d[fit == "floor_0.015625"]$omega2to4)), collapse = "-"), 3 / 64))
