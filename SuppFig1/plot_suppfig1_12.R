here::i_am("SuppFig1/plot_suppfig1_12.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))
d <- design_means(read_sim("Fig2", "main"))
panels <- list()
for (target in c("pleiotropy", "replication", "mediation")) {
  for (K in c(4L, 6L)) for (threshold in c(1L, 2L)) {
    z <- d[target == ..target & K == ..K & t == threshold & axis == "effect"]
    if (!nrow(z)) next
    panels[[length(panels) + 1L]] <- ggplot(z, aes(x, FDP, color = method, linetype = method)) +
      geom_line(linewidth = 1) + fdp_line() + coord_cartesian(ylim = c(0, fdp_top(z$FDP))) +
      xlab(EFFECT_XLAB[[target]]) + ylab("FDP") +
      ggtitle(sprintf("%s, K=%d, t=%d", TARGET_NAME[[target]], K, threshold)) + theme_paper(12)
  }
}
p <- cowplot::plot_grid(plotlist = lapply(panels, function(x) x + theme(legend.position = "none")),
                        ncol = 3L, labels = LETTERS[seq_along(panels)])
ggsave(here::here("SuppFig1", "SuppFig1_12.pdf"), p, width = 18, height = 20, device = PDF_DEVICE)
