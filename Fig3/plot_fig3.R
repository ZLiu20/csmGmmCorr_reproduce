# Figure 3: FDP and power of the correlated csmGmm against the correlation rho of the outcomes
# (pleiotropy, replication) or mediator errors (mediation); fifth effect-size design.
# Run after all arrays of this folder have finished:  Rscript Fig3/plot_fig3.R
here::i_am("Fig3/plot_fig3.R")
suppressMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))

R <- design_means(read_sim("Fig3", "rho"), by = c("target", "K", "t", "method", "design", "x"))
fwrite(R, here::here("Fig3", "output", "design_rho.csv"))
R <- R[method == "penalized"]
R[, Setting := factor(sprintf("K = %d, t = %d", K, t), levels = c("K = 4, t = 1", "K = 4, t = 2", "K = 6, t = 1", "K = 6, t = 2"))]
cols <- grDevices::hcl(h = seq(15, 375, length = 5)[1:4], c = 100, l = 65)
panels <- list()
for (tg in c("pleiotropy", "replication", "mediation")) for (what in c("FDP", "power")) {
  z <- as.data.frame(R[target == tg])
  p <- ggplot(z, aes(x = x, y = .data[[what]], group = Setting)) +
    geom_line(aes(color = Setting, linetype = Setting), linewidth = 1.2) +
    scale_color_manual(values = cols) + scale_linetype_manual(values = 1:4) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
    xlab(expression("Correlation " * rho)) +
    ylab(sprintf("%s (%s)", if (what == "FDP") "FDP" else "Power", TARGET_NAME[[tg]])) + theme_paper()
  p <- if (what == "FDP") p + fdp_line() + coord_cartesian(ylim = c(0, fdp_top(z$FDP))) else p + coord_cartesian(ylim = c(0, 1))
  panels[[length(panels) + 1]] <- p
}
lg <- ggplot(data.frame(x = rep(1:2, 4), y = 1, Setting = factor(rep(levels(R$Setting), each = 2), levels = levels(R$Setting))),
             aes(x, y, group = Setting)) + geom_line(aes(color = Setting, linetype = Setting), linewidth = 1.2) +
  scale_color_manual(values = cols) + scale_linetype_manual(values = 1:4) + theme_paper()
p <- combine_panels(panels, nrow = 3, legend = lg, rel = 0.08, ncol_legend = 4)
ggsave(here::here("Fig3", "Fig3.pdf"), p, width = 14, height = 18)
message("wrote Fig3/Fig3.pdf")
