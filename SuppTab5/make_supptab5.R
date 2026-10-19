# Supplementary Table S2: simulation designs (from design_table() in SupportingCode/sim_data.R).
#   Rscript SuppTab5/make_supptab5.R
here::i_am("SuppTab5/make_supptab5.R")
source(here::here("SupportingCode", "load_all.R"))
v <- function(scen, rows, col) paste(format(design_table("main", scen)[rows, col], nsmall = 2), collapse = ",")
row <- function(lab, rows, colpl, colmed, pre) sprintf("%s & $%s=%s$ & $%s=%s$ & $%s=%s$\\\\", lab, pre[1], v("pl", rows, colpl),
                                                      pre[1], v("rp", rows, colpl), pre[2], v("med", rows, colmed))
L <- c("\\begin{table}[htbp]", "\\centering",
  "\\caption{Simulation designs. $\\beta$: common non-null coefficient; $u$: lower end of the range $[u,u+0.18]$ of $|\\alpha_{jg}|$, with $|\\beta_{jg}|$ on $[0.14+u,0.24+u]$.}\\label{tab:sim-designs}",
  "\\footnotesize", "\\begin{tabular}{@{}lccc@{}}", "\\toprule", "Designs & Pleiotropy & Replication & Mediation\\\\", "\\midrule",
  row("1--5 (effect size)", 1:5, "effect", "a_min", c("\\beta", "u")),
  sprintf("6--9 (signal proportion) & $\\tau_1=%s$ & $\\tau_1=%s$ & $\\tau_1=%s$\\\\", v("pl", 6:9, "tau1"), v("rp", 6:9, "tau1"), v("med", 6:9, "tau1")),
  row("10--12 (effect size)", 10:12, "effect", "a_min", c("\\beta", "u")), "\\bottomrule", "\\end{tabular}", "\\end{table}")
write_tex(L, here::here("SuppTab5", "SuppTab5.tex"))
print(design_table("main", "pl")); print(design_table("main", "rp")); print(design_table("main", "med"))
