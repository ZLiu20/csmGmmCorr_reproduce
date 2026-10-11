# Supplementary Figures 1-12: one figure per target, K and t (pleiotropy K = 4, t = 1, 2;
# K = 6, t = 1, 2; then replication and mediation). Panels: mean FDP (A, B), power (C, D) and
# incongruous results per replicate (E, F) against the effect size and the signal proportion.
# Uses the output of the Fig2 folder:  Rscript SuppFig1/plot_suppfig1_12.R
here::i_am("SuppFig1/plot_suppfig1_12.R")
suppressMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))

D <- design_means(read_sim("Fig2", "main"))
S <- design_means(read_sim("Fig2", "sym"))
METH <- list(pleiotropy = c("penalized", "ML", "QCH", "adaFilter"), replication = c("penalized", "ML"),
             mediation = c("penalized", "ML", "DEI", "MediationTest", "PCPathSafeguard"))
k <- 0
for (tg in names(METH)) for (K0 in c(4, 6)) for (t0 in 1:2) {
  k <- k + 1
  z <- D[target == tg & K == K0 & t == t0]
  sy <- S[target == tg & K == K0 & t == t0 & method == "penalized"]
  lab <- sprintf("%s, K = %d, t = %d", TARGET_NAME[[tg]], K0, t0)
  top <- fdp_top(c(z[method %in% METH[[tg]]]$FDP, sy$FDP))
  inc_top <- max(1, 1.1 * max(z[method %in% METH[[tg]]]$incon_score, na.rm = TRUE))
  panels <- list()
  for (row in c("FDP", "power", "incon_score")) for (ax in c("effect", "tau1")) {
    xl <- if (ax == "effect") EFFECT_XLAB[[tg]] else expression("Signal proportion " * tau[1])
    yl <- switch(row, FDP = sprintf("FDP (%s)", lab), power = sprintf("Power (%s)", lab), incon_score = "Incongruous results per replicate")
    lim <- switch(row, FDP = c(0, top), power = c(0, 1), incon_score = c(0, inc_top))
    panels[[length(panels) + 1]] <- sim_panel(z[axis == ax], row, METH[[tg]], xl, yl, lim, hline = row == "FDP",
                                              sym = if (row == "FDP") sy[axis == ax] else NULL)
  }
  p <- combine_panels(panels, nrow = 3, legend = legend_plot(c(METH[[tg]], "sym")), rel = 0.1, ncol_legend = 2)
  f <- here::here("SuppFig1", sprintf("SuppFig%d_%s_K%d_t%d.pdf", k, tg, K0, t0))
  ggsave(f, p, width = 14, height = 18); message("wrote ", f)
}
