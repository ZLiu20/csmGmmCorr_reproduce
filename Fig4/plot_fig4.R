# Figure 4: applications. (A, B) loci (variants) discovered at q = 0.10 by the correlated
# csmGmm, the fit with Sigma fixed at the robust correlation and the csmGmm (Sigma = I);
# (C) discoveries kept by the monotone adjustment; (D) fitted against empirical correlations;
# (E, F) fitted Sigma (upper triangle) and robust correlation (lower triangle).
#   Rscript Fig4/plot_fig4.R      (after Fig4/make_tab3.R or directly after the fits)
here::i_am("Fig4/plot_fig4.R")
suppressMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))
X <- rd_collect(); S <- X$summary; G <- X$sigma
v <- function(an, va, t0, col) S[analysis == an & variant == va & t == t0][[col]][1]
FITS <- data.table(variant = c("main", "robust", "indep"),
                   label = c(sprintf("Correlated csmGmm (estimated %s)", SIGMA), sprintf("Correlated csmGmm (%s fixed at robust correlation)", SIGMA),
                             sprintf("csmGmm (%s = I)", SIGMA)),
                   color = c(grDevices::hcl(15, 100, 65), grDevices::hcl(15, 60, 82), "#5A5A5A"))
FITS[, label := factor(label, levels = label)]
bars <- function(an, title, variants_in_millions) {
  d <- S[analysis == an & variant %in% FITS$variant, .(variant, t, loci = loci_raw, R = R_raw)]
  d <- merge(d, FITS, by = "variant")
  d[, txt := if (variants_in_millions) sprintf("%s\n(%.2fM)", format(loci, big.mark = ","), R / 1e6) else
               sprintf("%s\n(%s)", format(loci, big.mark = ","), format(R, big.mark = ","))]
  ggplot(d, aes(x = factor(sprintf("t = %d", t)), y = loci, fill = label)) +
    geom_col(position = position_dodge(0.8), width = 0.76, color = "white") +
    geom_text(aes(label = txt), position = position_dodge(0.8), vjust = -0.2, size = 3.6, lineheight = 0.9) +
    scale_fill_manual(values = setNames(FITS$color, FITS$label), name = "Fit") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    xlab(NULL) + ylab("Loci discovered") + ggtitle(title) + theme_paper(15) + theme(plot.title = element_text(size = 18, face = "plain"))
}
pA <- bars("glgc_replication", "Replication, GLGC (K = 6)", FALSE)
pB <- bars("ukb_pleiotropy", "Pleiotropy, UK Biobank (K = 5)", TRUE)
adj <- rbindlist(lapply(c("glgc_replication", "ukb_pleiotropy"), function(an) data.table(
  analysis = c(glgc_replication = "Replication", ukb_pleiotropy = "Pleiotropy")[[an]], t = sprintf("t = %d", 1:2),
  pct = sapply(1:2, function(t0) 100 * v(an, "main", t0, "R_adj") / v(an, "main", t0, "R_raw")))))
adj[, analysis := factor(analysis, levels = c("Replication", "Pleiotropy"))]
pC <- ggplot(adj, aes(analysis, pct, fill = t)) + geom_col(position = position_dodge(0.8), width = 0.76) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_dodge(0.8), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("#3C3C3C", "#A8A8A8"), name = NULL) + coord_cartesian(ylim = c(0, 110)) +
  xlab(NULL) + ylab("Adjusted / unadjusted\ndiscoveries (%)") + ggtitle("Monotone adjustment") +
  theme_paper(15) + theme(plot.title = element_text(size = 18, face = "plain"), legend.position = c(0.85, 0.8), legend.key.size = unit(1.2, "line"))
gw <- dcast(G[variant == "main" & analysis %in% c("glgc_replication", "ukb_pleiotropy")], analysis + stat1 + stat2 ~ what, value.var = "value")
gw[, Analysis := c(glgc_replication = "Replication, GLGC (K = 6)", ukb_pleiotropy = "Pleiotropy, UK Biobank (K = 5)")[analysis]]
mcol <- c("Replication, GLGC (K = 6)" = grDevices::hcl(240, 100, 55), "Pleiotropy, UK Biobank (K = 5)" = grDevices::hcl(105, 100, 55))
gw[, open_shape := ifelse(analysis == "glgc_replication", 0, 1)]
pD <- ggplot(gw) + geom_abline(intercept = 0, slope = 1, color = "grey", linetype = 2) +
  geom_point(aes(pearson, fitted_Sigma, color = Analysis), shape = gw$open_shape, size = 3.5, stroke = 1, show.legend = FALSE) +
  geom_point(aes(robust, fitted_Sigma, color = Analysis, shape = Analysis), size = 3.5) +
  scale_color_manual(values = mcol, name = "filled: robust; open: Pearson") +
  scale_shape_manual(values = c("Replication, GLGC (K = 6)" = 15, "Pleiotropy, UK Biobank (K = 5)" = 16), name = "filled: robust; open: Pearson") +
  xlab("Empirical correlation of the statistics") + ylab(sprintf("Fitted %s (estimated)", SIGMA)) + theme_paper(15) +
  theme(legend.position = c(0.02, 0.88), legend.justification = c(0, 1),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12), legend.key.size = unit(1.2, "line"))
heat <- function(an, labs, title) {
  cols <- RD_ANALYSES[[an]]$cols; n <- length(cols); M <- diag(n)
  g <- G[analysis == an & variant == "main"]
  for (i in seq_len(nrow(g))) {
    a <- match(g$stat1[i], cols); b <- match(g$stat2[i], cols); lo <- min(a, b); hi <- max(a, b)
    if (g$what[i] == "fitted_Sigma") M[lo, hi] <- g$value[i]
    if (g$what[i] == "robust") M[hi, lo] <- g$value[i]
  }
  d <- data.table(row = rep(seq_len(n), n), col = rep(seq_len(n), each = n), value = as.vector(M))[row != col]
  ggplot(d, aes(col, row, fill = value)) + geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value), color = abs(value) > 0.6), size = 4.2, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), guide = "none") +
    scale_x_continuous(breaks = seq_len(n), labels = labs, expand = c(0, 0)) +
    scale_y_reverse(breaks = seq_len(n), labels = labs, expand = c(0, 0)) + coord_equal() +
    xlab(NULL) + ylab(NULL) + ggtitle(title) + theme_cowplot(14) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(size = 14, face = "plain"))
}
pE <- heat("glgc_replication", c("EUR F", "EUR M", "EAS F", "EAS M", "AFR F", "AFR M"), sprintf("Replication: fitted %s (upper)\nvs robust correlation (lower)", SIGMA))
pF <- heat("ukb_pleiotropy", c("BMI", "SBP", "TG", "HDL", "HbA1c"), sprintf("Pleiotropy: fitted %s (upper)\nvs robust correlation (lower)", SIGMA))
g <- plot_grid(pA + theme(legend.position = "none"), pB + theme(legend.position = "none"), pC, pD, pE, pF,
               labels = LETTERS[1:6], nrow = 2, label_size = 22)
lg <- legend_grob(pA + theme(legend.position = "bottom", legend.text = element_text(size = 15)))
ggsave(here::here("Fig4", "Fig4.pdf"), plot_grid(g, lg, ncol = 1, rel_heights = c(1, 0.07)), width = 20, height = 13, device = PDF_DEVICE)
message("wrote Fig4/Fig4.pdf")
