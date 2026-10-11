here::i_am("SuppTab2/make_supptab2.R")
source(here::here("SupportingCode", "load_all.R"))
counts <- function(K, S, t = 2L) {
  H <- as.matrix(expand.grid(rep(list(-1L:1L), K)))
  support <- rowSums(abs(H))
  represented <- support <= S
  same_direction <- support == abs(rowSums(sign(H)))
  paired <- matrix(seq_len(K), ncol = 2L, byrow = TRUE)
  complete_path <- apply(H, 1L, function(h) {
    any(h[paired[, 1L]] != 0L & h[paired[, 2L]] != 0L)
  })
  eligible <- represented & support > t
  c(K = K, S = S, states = sum(represented),
    pleiotropy = sum(eligible),
    replication = sum(eligible & same_direction),
    mediation = sum(eligible & complete_path))
}
tab <- rbind(counts(4,4), counts(10,4), counts(10,10))
lines <- c("\\begin{table}[!htbp]", "\\centering",
  "\\caption{Represented and alternative-state counts at $t=2$ for the full and restricted models.}\\label{tab:counts}",
  "\\begin{tabular}{rrrrrr}\\toprule", "$K$ & $S$ & All states & Pleiotropy & Replication & Mediation\\\\\\midrule")
for (i in seq_len(nrow(tab))) lines <- c(lines, paste(formatC(tab[i,], format="d", big.mark=","), collapse=" & "), "\\\\")
write_tex(c(lines, "\\bottomrule\\end{tabular}", "\\end{table}"), here::here("SuppTab2", "SuppTab2.tex"))
