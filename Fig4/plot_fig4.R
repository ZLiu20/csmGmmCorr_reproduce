here::i_am("Fig4/plot_fig4.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))
x <- rd_collect(here::here("Fig4", "output"))
s <- x$summary
if (!nrow(s)) stop("No real-data summaries were found.")
p1 <- ggplot(s, aes(factor(t), R_raw, fill = variant)) + geom_col(position = "dodge") +
  facet_wrap(~analysis, scales = "free_y") + xlab("Threshold t") + ylab("Discoveries") + theme_paper(13)
p2 <- ggplot(s, aes(factor(t), R_adjusted / pmax(R_raw, 1) * 100, fill = variant)) + geom_col(position = "dodge") +
  facet_wrap(~analysis) + xlab("Threshold t") + ylab("Adjusted/raw discoveries (%)") + theme_paper(13)
p <- plot_grid(p1 + theme(legend.position="bottom"), p2 + theme(legend.position="bottom"), nrow=2, labels=c("A","B"))
ggsave(here::here("Fig4", "Fig4.pdf"), p, width = 14, height = 12, device = PDF_DEVICE)
