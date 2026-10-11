here::i_am("Fig2/plot_fig2.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))
D <- design_means(read_sim("Fig2", "main"))
S <- design_means(read_sim("Fig2", "sym"))
fwrite(D, here::here("Fig2", "output", "design_main.csv"))
fwrite(S, here::here("Fig2", "output", "design_sym.csv"))
methods <- list(pleiotropy = c("penalized", "QCH", "adaFilter"),
                replication = "penalized",
                mediation = c("penalized", "DEI", "MediationTest", "PCPathSafeguard"))
panels <- list()
for (target in names(methods)) {
  z <- D[target == ..target & K == 6L & t == 1L]
  sym <- S[target == ..target & K == 6L & t == 1L & method == "penalized"]
  for (axis in c("effect", "tau1")) {
    xlab <- if (axis == "effect") EFFECT_XLAB[[target]] else expression("Signal proportion " * tau[1])
    panels[[paste(target, "FDP", axis)]] <- sim_panel(z[axis == ..axis], "FDP", methods[[target]], xlab,
      sprintf("FDP (%s)", TARGET_NAME[[target]]), c(0, fdp_top(z$FDP)), TRUE, sym[axis == ..axis])
    panels[[paste(target, "power", axis)]] <- sim_panel(z[axis == ..axis], "power", methods[[target]], xlab,
      sprintf("Power (%s)", TARGET_NAME[[target]]), c(0, 1), FALSE, sym[axis == ..axis])
  }
}
legend <- legend_plot(c("penalized", "QCH", "adaFilter", "DEI", "MediationTest", "PCPathSafeguard", "sym"))
plot <- combine_panels(panels, nrow = 3L, legend = legend, rel = .1, ncol_legend = 4L)
ggsave(here::here("Fig2", "Fig2.pdf"), plot, width = 24, height = 17, device = PDF_DEVICE)
