here::i_am("Fig3/plot_fig3.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))
R <- design_means(read_sim("Fig3", "rho"), by = c("target", "K", "t", "method", "design", "x"))
fwrite(R, here::here("Fig3", "output", "design_rho.csv"))
R <- R[method == "penalized"]
R[, Setting := factor(sprintf("K = %d, t = %d", K, t),
  levels = c("K = 4, t = 1", "K = 4, t = 2", "K = 6, t = 1", "K = 6, t = 2"))]
colors <- grDevices::hcl(seq(15, 375, length.out = 5L)[1:4], 100, 65)
panels <- list()
for (target in c("pleiotropy", "replication", "mediation")) for (metric in c("FDP", "power")) {
  z <- R[target == ..target]
  p <- ggplot(z, aes(x, .data[[metric]], group = Setting)) +
    geom_line(aes(color = Setting, linetype = Setting), linewidth = 1.2) +
    scale_color_manual(values = colors) + scale_linetype_manual(values = 1:4) +
    scale_x_continuous(breaks = c(0, .2, .4, .6, .8)) +
    xlab(expression("Correlation " * rho)) +
    ylab(sprintf("%s (%s)", ifelse(metric == "FDP", "FDP", "Power"), TARGET_NAME[[target]])) + theme_paper()
  if (metric == "FDP") p <- p + fdp_line() + coord_cartesian(ylim = c(0, .5)) else p <- p + coord_cartesian(ylim = c(0, 1))
  panels[[length(panels) + 1L]] <- p
}
legend <- panels[[1L]]
plot <- combine_panels(panels, nrow = 3L, legend = legend, rel = .08, ncol_legend = 4L)
ggsave(here::here("Fig3", "Fig3.pdf"), plot, width = 14, height = 18, device = PDF_DEVICE)
