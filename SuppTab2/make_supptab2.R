# Supplementary Table 2: represented and alternative-state counts at t = 2 for the full model
# (S = K) and for a restriction to supports of size at most S = 4 at K = 10 (exact counts).
#   Rscript SuppTab2/make_supptab2.R
here::i_am("SuppTab2/make_supptab2.R")
source(here::here("SupportingCode", "load_all.R"))
counts <- function(K, S, t = 2) {
  s <- 0:S; nconf <- 2^s * choose(K, s)
  pair <- vapply(s, function(x) complete_pair_fraction(K, x), 0)   # sim_data.R
  c(K = K, S = S, all = sum(nconf), pleiotropy = sum(nconf[s > t]), replication = sum(2 * choose(K, s)[s > t]),
    mediation = sum((nconf * pair)[s > t]))
}
tab <- rbind(counts(4, 4), counts(10, 4), counts(10, 10))
print(tab)
fmt <- function(x) formatC(x, format = "d", big.mark = ",")
L <- c("\\begin{table}[htbp]", "\\centering",
  "\\caption{Represented and alternative-state counts at $t=2$ for the full model ($S=K$) and, for illustration, for a restriction to supports of size at most $S=4$ at $K=10$.}",
  "\\label{tab:counts}", "\\begin{tabular}{@{}rrrrrr@{}}", "\\toprule", "$K$ & $S$ & All states & Pleiotropy & Replication & Mediation\\\\", "\\midrule")
for (i in seq_len(nrow(tab))) L <- c(L, paste0(paste(fmt(tab[i, ]), collapse = " & "), "\\\\"))
write_tex(c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}"), here::here("SuppTab2", "SuppTab2.tex"))
