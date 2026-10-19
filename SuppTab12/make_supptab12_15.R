# Supplementary Tables S8-S11 (robustness studies; folder numbering 12-15):
#   12  mean FDP at t = 2 under alternative estimators and data-generating processes
#   13  mean power for the settings of Table 2 and Supplementary Table S8
#   14  linkage disequilibrium study: variant-level, locus-level and union-of-effects FDP
#   15  replication comparators (repfdr, directional partial conjunction) by design
# Arrays of this folder: main_comp_rp (Table 15). Uses the output of Fig2, Fig3 and Tab2.
#   Rscript SuppTab12/make_supptab12_15.R
here::i_am("SuppTab12/make_supptab12_15.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
TG <- c("pleiotropy", "replication", "mediation")

# ---- Table 12 ----
T <- sources_means("FDP")
cap12 <- paste0("Mean FDP at $q=0.10$ and $t=2$ under alternative estimators and data-generating processes, averaged ",
  "over the first and fifth effect-size designs and one signal-proportion design of each target ($\\tau_1=0.02$ for ",
  "pleiotropy and replication, $0.04$ for mediation); the last block uses the fifth effect-size design with $\\rho=0.8$. ",
  "Replicates: 500 (first row), 200 (estimation rows and $\\rho=0.8$; 50 for the fixed null correlation), 100 (independent ",
  "test set and other data-generating rows; 200 for symmetric statistics). Power is in Table~\\ref{tab:sources-power}.")
write_tex(sources_tex(T, 2, "tab:sources-t2", cap12), here::here("SuppTab12", "SuppTab12.tex"))

# ---- Table 13 ----
P <- sources_means("power")
fwrite(P, here::here("SuppTab12", "SuppTab13_sources_power.csv"))
L <- c("{\\footnotesize", "\\begin{longtable}{l cccc cccc cccc}",
  "\\caption{Mean power for the settings of main-text Table~\\ref{main-tab:sources} and Table~\\ref{tab:sources-t2}.}\\label{tab:sources-power}\\\\",
  "\\toprule", " & \\multicolumn{4}{c}{Pleiotropy} & \\multicolumn{4}{c}{Replication} & \\multicolumn{4}{c}{Mediation}\\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-13}",
  "$(K,t)$ & (4,1) & (4,2) & (6,1) & (6,2) & (4,1) & (4,2) & (6,1) & (6,2) & (4,1) & (4,2) & (6,1) & (6,2)\\\\",
  "\\midrule", "\\endfirsthead", "\\bottomrule", "\\endfoot")
for (r in SOURCE_ROWS) {
  if (is.na(r[1])) { L <- c(L, sprintf("\\multicolumn{13}{l}{%s}\\\\", r[2])); next }
  cells <- unlist(lapply(TG, function(tg) unlist(lapply(c(4, 6), function(K0) sapply(1:2, function(t0) {
    v <- P[target == tg & K == K0 & t == t0][[r[1]]]; f3(if (length(v)) v else NA) })))))
  L <- c(L, paste0(r[2], " & ", paste(cells, collapse = " & "), "\\\\"))
}
write_tex(c(L, "\\end{longtable}}"), here::here("SuppTab12", "SuppTab13.tex"))

# ---- Table 14: LD study ----
ld <- read_sim("Tab2", "ld")[method == "penalized" & target %in% c("pleiotropy", "replication")]
G <- ld[, .(FDP = mean(FDP), FDP_loc = mean(FDP_loc), FDP_union1 = mean(FDP_union1), FDP_union2 = mean(FDP_union2),
            FDP_union5 = mean(FDP_union5), power = mean(power, na.rm = TRUE), power_union2 = mean(power_union2, na.rm = TRUE), n = .N),
        by = .(target, K, t, design)]
setorder(G, target, K, t, design)
fwrite(G, here::here("SuppTab12", "SuppTab14_ld.csv"))
L <- c("\\begin{table}[!htbp]", "\\centering",
  paste0("\\caption{Linkage disequilibrium study (latent haplotype correlation 0.5 between neighboring variants; ",
         "100 replicates per design). FDP: variant level, with the truth defined by the causal variants; locus level: a ",
         "discovery counts as true if a causal alternative lies within five positions; union: the configuration of a ",
         "variant is the signed union of the effects within two positions.}\\label{tab:ld}"),
  "\\small", "\\begin{tabular}{lccl cccc}", "\\toprule",
  "Target & $K$ & $t$ & Design & FDP & Locus-level FDP & Union FDP & Power\\\\", "\\midrule")
for (i in seq_len(nrow(G))) with(G[i], L <<- c(L, sprintf("%s & %d & %d & %s & %s & %s & %s & %s\\\\", TARGET_NAME[[target]], K, t,
                                                       tex_escape(design), f3(FDP), f3(FDP_loc), f3(FDP_union2), f3(power))))
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab12", "SuppTab14.tex"))

# ---- Table 15: replication comparators ----
rp <- read_sim("SuppTab12", "main")
pen <- read_sim("Fig2", "main")[target == "replication" & method == "penalized"]
ran <- rp[method == "repfdr", .(ran = 100 * mean(is.na(error) | error == "")), by = .(K, t, design)]
rf <- rp[method == "repfdr" & (is.na(error) | error == ""), .(rf_FDP = mean(FDP), rf_pow = mean(power, na.rm = TRUE)), by = .(K, t, design)]
dp <- rp[method == "DirectionalPC", .(dp_FDP = mean(FDP), dp_pow = mean(power, na.rm = TRUE)), by = .(K, t, design)]
cs <- pen[, .(cs_FDP = mean(FDP), cs_pow = mean(power, na.rm = TRUE)), by = .(K, t, design)]
G <- Reduce(function(a, b) merge(a, b, by = c("K", "t", "design"), all.x = TRUE), list(ran, rf, dp, cs))
setorder(G, K, t, design)
fwrite(G, here::here("SuppTab12", "SuppTab15_rep_comparators.csv"))
H <- "$K$ & $t$ & Design & Ran (\\%) & FDP & Power & FDP & Power & FDP & Power\\\\"
L <- c("{\\footnotesize\\setlength{\\tabcolsep}{4pt}", "\\begin{longtable}{cc l r cc cc cc}",
  paste0("\\caption{Replication comparators by design. Ran: percentage of replicates in which repfdr returned results (it ",
         "stops when no coordinate has an estimated null proportion below one); its FDP and power are averaged over those ",
         "replicates only. Directional PC: one-sided partial-conjunction $p$-values for more than $t$ effects in a common ",
         "direction with BH, which ignores effects of the opposite sign and is therefore not valid for the replication null. ",
         "500 replicates for the correlated csmGmm, 200 for the comparators.}\\label{tab:rep-comparators}\\\\"),
  "\\toprule", "$K$ & $t$ & Design & Ran (\\%) & \\multicolumn{2}{c}{repfdr} & \\multicolumn{2}{c}{Directional PC} & \\multicolumn{2}{c}{Correlated csmGmm}\\\\",
  "\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}", " & & & & FDP & Power & FDP & Power & FDP & Power\\\\",
  "\\midrule", "\\endfirsthead", "\\toprule", H, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot")
prev <- NULL
for (i in seq_len(nrow(G))) {
  r <- G[i]; key <- paste(r$K, r$t)
  if (!is.null(prev) && key != prev) L <- c(L, "\\midrule")
  prev <- key
  L <- c(L, sprintf("%d & %d & %s & %.0f & %s & %s & %s & %s & %s & %s\\\\", r$K, r$t, tex_escape(r$design), r$ran,
                    f3(r$rf_FDP), f3(r$rf_pow), f3(r$dp_FDP), f3(r$dp_pow), f3(r$cs_FDP), f3(r$cs_pow)))
}
write_tex(c(L, "\\end{longtable}}"), here::here("SuppTab12", "SuppTab15.tex"))
