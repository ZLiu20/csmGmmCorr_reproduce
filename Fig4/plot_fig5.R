# Figure 5: discoveries of the correlated csmGmm at t = 1 and q = 0.10, -log10(lfdr) against
# genomic position (values above 20 shown at 20), in the style of Sun, McCaw and Lin (Figure 4).
# (A) replication: all discoveries, by whether they are also rejected at t = 2;
# (B) pleiotropy: the lead variant of each locus, by whether it is also rejected at t = 2.
#   Rscript Fig4/plot_fig5.R
here::i_am("Fig4/plot_fig5.R")
suppressMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))
CAP <- 20
od <- function(an) here::here("Fig4", "output", paste0(an, "_main"))
rd <- function(an, t0) fread(file.path(od(an), sprintf("discoveries_t%d.tsv.gz", t0)), select = c("id", "chr", "pos", "lfdr", "rej_raw"))
offsets <- function(...) {
  a <- rbindlist(lapply(list(...), function(d) d[, .(chr = as.integer(chr), pos)]))
  mx <- a[, .(mx = max(pos)), by = chr][match(1:22, chr)]$mx; mx[is.na(mx)] <- 0
  list(off = c(0, cumsum(mx)[-22]), mx = mx)
}
panel <- function(d, o, cats, title, legtitle) {
  d <- copy(d); d[, chr := as.integer(chr)]
  d[, x := pos + o$off[chr]]; d[, y := -log10(pmax(lfdr, 10^-CAP))]
  d[, cat := factor(cat, levels = cats)]
  labs <- sprintf("%s (%s)", cats, formatC(as.vector(table(d$cat)[cats]), format = "d", big.mark = ","))
  brk <- o$off + o$mx / 2; lab <- ifelse(1:22 %in% c(9, 11, 13, 15, 16, 17, 19, 20, 21), "", as.character(1:22))
  ggplot(d, aes(x, y, color = cat, shape = cat)) + geom_point(size = if (nrow(d) > 1e5) 0.6 else 1.6) +
    scale_color_manual(values = setNames(grDevices::hcl(c(15, 255), 100, 65), cats), labels = labs, name = legtitle) +
    scale_shape_manual(values = setNames(c(16, 17), cats), labels = labs, name = legtitle) +
    scale_x_continuous(breaks = brk, labels = lab, limits = c(0, sum(o$mx) * 1.005), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, CAP + 0.8)) + xlab("Chr") + ylab(expression(-log[10](lfdr))) + ggtitle(title) +
    guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4))) +
    theme_paper(16) + theme(plot.title = element_text(size = 18, face = "plain"), legend.key.size = unit(1.5, "line"))
}
# replication
r1 <- rd("glgc_replication", 1)[rej_raw == TRUE]; r2 <- rd("glgc_replication", 2)[rej_raw == TRUE]
r1[, cat := ifelse(id %chin% r2$id, "Also rejected at t = 2", "Rejected at t = 1 only")]
ps <- fread(file.path(od("glgc_replication"), "plot_sample.tsv.gz"), select = c("chr", "pos"))
pA <- panel(r1, offsets(r1, ps), c("Rejected at t = 1 only", "Also rejected at t = 2"),
            "Replication of LDL association across six ancestry-by-sex strata, GLGC", "Replication")
# pleiotropy: leads of the loci at t = 1
p1 <- rd("ukb_pleiotropy", 1)[rej_raw == TRUE]; p2 <- rd("ukb_pleiotropy", 2)[rej_raw == TRUE]
p1 <- p1[clump_leads(as.character(p1$chr), p1$pos, p1$lfdr, rep(TRUE, nrow(p1)))]
p1[, cat := ifelse(id %chin% p2$id, "Locus lead also rejected at t = 2", "Locus lead rejected at t = 1 only")]
psu <- fread(file.path(od("ukb_pleiotropy"), "plot_sample.tsv.gz"), select = c("chr", "pos"))
pB <- panel(p1, offsets(p1, psu), c("Locus lead rejected at t = 1 only", "Locus lead also rejected at t = 2"),
            "Pleiotropy across five metabolic-syndrome traits, UK Biobank (lead variant of each locus)", "Pleiotropy")
p <- plot_grid(pA, pB, labels = c("A", "B"), ncol = 1, label_size = 22)
ggsave(here::here("Fig4", "Fig5.pdf"), p, width = 16, height = 11.5)
message(sprintf("wrote Fig4/Fig5.pdf (%d replication discoveries, %d pleiotropy loci)", nrow(r1), nrow(p1)))
