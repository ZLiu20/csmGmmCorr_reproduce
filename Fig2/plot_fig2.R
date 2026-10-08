# Figure 2: FDP and power in the main simulation, K = 6, t = 1 (rows: pleiotropy, replication,
# mediation; columns: FDP against effect size and signal proportion, power against both).
# Run after all arrays of this folder have finished:  Rscript Fig2/plot_fig2.R
here::i_am("Fig2/plot_fig2.R")
suppressMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))

D <- design_means(read_sim("Fig2", "main"))
S <- design_means(read_sim("Fig2", "sym"))
fwrite(D, here::here("Fig2", "output", "design_main.csv")); fwrite(S, here::here("Fig2", "output", "design_sym.csv"))
METH <- list(pleiotropy = c("penalized", "QCH", "adaFilter"), replication = "penalized",
             mediation = c("penalized", "DEI", "MediationTest", "PCPathSafeguard"))
K0 <- 6; t0 <- 1
panels <- list()
for (tg in names(METH)) {
  z <- D[target == tg & K == K0 & t == t0]
  sy <- S[target == tg & K == K0 & t == t0 & method == "penalized"]
  top <- fdp_top(c(z[method %in% METH[[tg]]]$FDP, sy$FDP))
  for (ax in c("effect", "tau1")) {
    xl <- if (ax == "effect") EFFECT_XLAB[[tg]] else expression("Signal proportion " * tau[1])
    panels[[paste(tg, "FDP", ax)]] <- sim_panel(z[axis == ax], "FDP", METH[[tg]], xl, sprintf("FDP (%s)", TARGET_NAME[[tg]]),
                                                c(0, top), hline = TRUE, sym = sy[axis == ax])
  }
  for (ax in c("effect", "tau1")) {
    xl <- if (ax == "effect") EFFECT_XLAB[[tg]] else expression("Signal proportion " * tau[1])
    panels[[paste(tg, "power", ax)]] <- sim_panel(z[axis == ax], "power", METH[[tg]], xl, sprintf("Power (%s)", TARGET_NAME[[tg]]), c(0, 1))
  }
}
leg <- legend_plot(c("penalized", "QCH", "adaFilter", "DEI", "MediationTest", "PCPathSafeguard", "sym"))
p <- combine_panels(panels, nrow = 3, legend = leg, rel = 0.1, ncol_legend = 4)
ggsave(here::here("Fig2", "Fig2.pdf"), p, width = 24, height = 17)
message("wrote Fig2/Fig2.pdf")
