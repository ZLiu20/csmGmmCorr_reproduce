here::i_am("Fig4/plot_fig5.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(here::here("SupportingCode", "load_all.R"))
path <- here::here("Fig4", "output", "mediation_ldl", "combined_summary.tsv")
if (!file.exists(path)) stop("Run Fig4/summarize_mediation_ldl.R first.")
d <- fread(path)
p <- ggplot(d, aes(factor(t), R_raw, color = variant, group = variant)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
  xlab("Threshold t") + ylab("Mediation discoveries") + theme_paper()
ggsave(here::here("Fig4", "Fig5.pdf"), p, width = 10, height = 7, device = PDF_DEVICE)
