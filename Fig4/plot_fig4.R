# Figure 4: the three applications, one column each (replication, pleiotropy, mediation).
# (A-C) loci (variants) or variant-gene pairs discovered at q = 0.10 by the correlated csmGmm,
# the fit with Sigma fixed at the robust correlation and the csmGmm (Sigma = I);
# (D-F) fitted Sigma (upper triangle) and robust correlation (lower triangle);
# (G) fitted against empirical correlations; (H) discoveries kept by the monotone adjustment;
# (I) support in external statistics (held-out GLGC strata; MVP for mediation).
#   Rscript Fig4/plot_fig4.R   (after the fits, Fig4/holdout_validation.R, the mediation fits,
#                               Fig4/mediation_scc_mvp.R and Fig4/mediation_scc_report.R)
here::i_am("Fig4/plot_fig4.R")
suppressMessages({ library(data.table); library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "load_all.R"))
X <- rd_collect(); S <- X$summary; G <- X$sigma
# lung SCC mediation (cis pairs): summaries of the three fits and Sigma of the estimated fit
MD <- here::here("Fig4", "output", "mediation_scc")
MEDSCC_COLS <- c("alpha_lung", "beta_lung", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver")
# The mediation data are not public (README); without them the mediation panels are left empty.
has_med <- all(file.exists(c(file.path(MD, paste0("cis15_", c("main", "robust", "indep")), "summary.tsv"),
                             file.path(MD, "cis15_main", "sigma.tsv"), file.path(MD, "report", "mvp_support.tsv"))))
if (has_med) {
  ms <- rbindlist(lapply(c("main", "robust", "indep"), function(va) fread(file.path(MD, paste0("cis15_", va), "summary.tsv"))[, variant := va]))
  S <- rbind(S, ms[, .(analysis = "mediation_scc", variant, t, R_raw, R_adj, loci_raw = R_raw, loci_adj = R_adj)], fill = TRUE)
  mg <- fread(file.path(MD, "cis15_main", "sigma.tsv"))
  G <- rbind(G, melt(mg[, .(stat1, stat2, fitted_Sigma = fitted, robust, pearson)], id.vars = c("stat1", "stat2"),
                     variable.name = "what", value.name = "value")[, `:=`(analysis = "mediation_scc", variant = "main", what = as.character(what))], fill = TRUE)
}
RD_ANALYSES$mediation_scc <- list(cols = MEDSCC_COLS)
no_med <- function() ggdraw() + draw_label("Mediation data not available\n(see README)", size = 14, colour = "grey40")
v <- function(an, va, t0, col) S[analysis == an & variant == va & t == t0][[col]][1]
FITS <- data.table(variant = c("main", "robust", "indep"),
                   label = c(sprintf("Correlated csmGmm (estimated %s)", SIGMA), sprintf("Correlated csmGmm (%s fixed at robust correlation)", SIGMA),
                             sprintf("csmGmm (%s = I)", SIGMA)),
                   color = c(grDevices::hcl(15, 100, 65), grDevices::hcl(15, 60, 82), "#5A5A5A"))
FITS[, label := factor(label, levels = label)]
bars <- function(an, title, variants_in_millions, pairs = FALSE) {
  d <- S[analysis == an & variant %in% FITS$variant, .(variant, t, loci = loci_raw, R = R_raw)]
  d <- merge(d, FITS, by = "variant")
  d[, txt := if (pairs) format(loci) else if (variants_in_millions) sprintf("%s\n(%.2fM)", format(loci, big.mark = ","), R / 1e6) else
               sprintf("%s\n(%s)", format(loci, big.mark = ","), format(R, big.mark = ","))]
  ggplot(d, aes(x = factor(sprintf("t = %d", t)), y = loci, fill = label)) +
    geom_col(position = position_dodge(0.8), width = 0.76, color = "white") +
    geom_text(aes(label = txt), position = position_dodge(0.8), vjust = -0.2, size = 3.6, lineheight = 0.9) +
    scale_fill_manual(values = setNames(FITS$color, FITS$label), name = "Fit") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    xlab(NULL) + ylab(if (pairs) "Pairs discovered" else "Loci discovered") + ggtitle(title) + theme_paper(15) + theme(plot.title = element_text(size = 18, face = "plain"))
}
pA <- bars("glgc_replication", "Replication, GLGC (K = 6)", FALSE)
pB <- bars("ukb_pleiotropy", "Pleiotropy, UK Biobank (K = 5)", TRUE)
pC <- if (has_med) bars("mediation_scc", "Mediation, lung SCC (K = 6)", FALSE, pairs = TRUE) else no_med()
AN3 <- c(glgc_replication = "Replication", ukb_pleiotropy = "Pleiotropy", mediation_scc = "Mediation")
if (!has_med) AN3 <- AN3[1:2]
adj <- rbindlist(lapply(names(AN3), function(an) data.table(
  analysis = AN3[[an]], t = sprintf("t = %d", 1:2),
  pct = sapply(1:2, function(t0) 100 * v(an, "main", t0, "R_adj") / v(an, "main", t0, "R_raw")))))
adj[, analysis := factor(analysis, levels = AN3)]
pH <- ggplot(adj, aes(analysis, pct, fill = t)) + geom_col(position = position_dodge(0.8), width = 0.76) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_dodge(0.8), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("#3C3C3C", "#A8A8A8"), name = NULL) + coord_cartesian(ylim = c(0, 118)) +
  xlab(NULL) + ylab("Adjusted / unadjusted\ndiscoveries (%)") + ggtitle("Monotone adjustment") +
  theme_paper(15) + theme(plot.title = element_text(size = 18, face = "plain"), legend.position = c(0.05, 0.95), legend.direction = "horizontal", legend.key.size = unit(1.2, "line"))
TT <- c(glgc_replication = "Replication, GLGC (K = 6)", ukb_pleiotropy = "Pleiotropy, UK Biobank (K = 5)", mediation_scc = "Mediation, lung SCC (K = 6)")
gw <- dcast(G[variant == "main" & analysis %in% names(TT) & what %in% c("fitted_Sigma", "robust", "pearson")], analysis + stat1 + stat2 ~ what, value.var = "value")
gw[, Analysis := factor(TT[analysis], levels = TT)]
mcol <- setNames(c(grDevices::hcl(240, 100, 55), grDevices::hcl(105, 100, 55), grDevices::hcl(285, 70, 45)), TT)
gw[, open_shape := c(glgc_replication = 0, ukb_pleiotropy = 1, mediation_scc = 2)[analysis]]
pG <- ggplot(gw) + geom_abline(intercept = 0, slope = 1, color = "grey", linetype = 2) +
  geom_point(aes(pearson, fitted_Sigma, color = Analysis), shape = gw$open_shape, size = 3.5, stroke = 1, show.legend = FALSE) +
  geom_point(aes(robust, fitted_Sigma, color = Analysis, shape = Analysis), size = 3.5) +
  scale_color_manual(values = mcol, name = "filled: robust; open: Pearson") +
  scale_shape_manual(values = setNames(c(15, 16, 17), TT), name = "filled: robust; open: Pearson") +
  xlab("Empirical correlation of the statistics") + ylab(sprintf("Fitted %s (estimated)", SIGMA)) +
  ggtitle("Fitted and empirical correlations") + theme_paper(15) + theme(plot.title = element_text(size = 18, face = "plain")) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12), legend.key.size = unit(1.2, "line"))
heat <- function(an, labs, title, blocks = FALSE) {
  cols <- RD_ANALYSES[[an]]$cols; n <- length(cols); M <- diag(n)
  g <- G[analysis == an & variant == "main"]
  for (i in seq_len(nrow(g))) {
    a <- match(g$stat1[i], cols); b <- match(g$stat2[i], cols); lo <- min(a, b); hi <- max(a, b)
    if (g$what[i] == "fitted_Sigma") M[lo, hi] <- g$value[i]
    if (g$what[i] == "robust") M[hi, lo] <- g$value[i]
  }
  d <- data.table(row = rep(seq_len(n), n), col = rep(seq_len(n), each = n), value = as.vector(M))[row != col]
  ggplot(d, aes(col, row, fill = value)) + geom_tile(color = "white") +
    geom_text(aes(label = sub("^-0.00$", "0.00", sprintf("%.2f", value)), color = abs(value) > 0.6), size = 4.2, show.legend = FALSE) +
    (if (blocks) list(geom_hline(yintercept = n / 2 + 0.5), geom_vline(xintercept = n / 2 + 0.5)) else NULL) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), guide = "none") +
    scale_x_continuous(breaks = seq_len(n), labels = labs, expand = c(0, 0)) +
    scale_y_reverse(breaks = seq_len(n), labels = labs, expand = c(0, 0)) + coord_equal() +
    xlab(NULL) + ylab(NULL) + ggtitle(title) + theme_cowplot(14) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(size = 14, face = "plain"))
}
pD <- heat("glgc_replication", c("EUR F", "EUR M", "EAS F", "EAS M", "AFR F", "AFR M"), sprintf("Replication: fitted %s (upper)\nvs robust correlation (lower)", SIGMA))
pE <- heat("ukb_pleiotropy", c("BMI", "SBP", "TG", "HDL", "HbA1c"), sprintf("Pleiotropy: fitted %s (upper)\nvs robust correlation (lower)", SIGMA))
# mediation: alpha coordinates first, then beta, so that the two blocks are visible
RD_ANALYSES$mediation_scc$cols <- c("alpha_lung", "alpha_blood", "alpha_liver", "beta_lung", "beta_blood", "beta_liver")
GR <- if (isTRUE(capabilities("cairo"))) c("\u03b1", "\u03b2") else c("alpha", "beta")
pF <- if (!has_med) no_med() else heat("mediation_scc", c(paste(GR[1], c("Lung", "Blood", "Liver")), paste(GR[2], c("Lung", "Blood", "Liver"))),
           sprintf("Mediation: fitted %s (upper)\nvs robust correlation (lower)", SIGMA), blocks = TRUE)
# I: external support. Replication: summed held-out strata at lead variants (holdout_validation.R);
# mediation: MVP S-PrediXcan statistics of the pairs (mediation_scc_report.R), t = 1.
ho <- fread(here::here("Fig4", "output", "holdout_validation.csv"))
mv <- if (has_med) fread(file.path(MD, "report", "mvp_support.tsv")) else NULL
GRP <- c("Correlated csmGmm", sprintf("%s = I only", SIGMA), "Background")
ext <- rbind(data.table(app = "Replication\n(held-out strata)", group = GRP,
                        pct = 100 * ho[match(c("Correlated csmGmm", "Sigma = I only", "Random variants"), set)]$S_gt_1.645),
             if (has_med) data.table(app = "Mediation\n(MVP)", group = GRP,
                        pct = 100 * mv[match(c("correlated csmGmm (all)", "Sigma = I only", "neither"), group)]$gt_1.645))
ext[, `:=`(group = factor(group, levels = GRP), app = factor(app, levels = unique(app)))]
pI <- ggplot(ext, aes(app, pct, fill = group)) + geom_col(position = position_dodge(0.8), width = 0.76) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_dodge(0.8), vjust = -0.3, size = 4) +
  scale_fill_manual(values = setNames(c(FITS$color[1], FITS$color[3], "#C8C8C8"), GRP), name = NULL) +
  coord_cartesian(ylim = c(0, 70)) + xlab(NULL) + ylab("Signed statistics beyond\none-sided 0.05 (%)") + ggtitle("External support") +
  theme_paper(15) + theme(plot.title = element_text(size = 18, face = "plain"), legend.position = c(0.55, 0.95),
                          legend.text = element_text(size = 12), legend.key.size = unit(1.2, "line"))
g <- plot_grid(pA + theme(legend.position = "none"), pB + theme(legend.position = "none"), pC + theme(legend.position = "none"),
               pD, pE, pF, pG, pH, pI, labels = LETTERS[1:9], nrow = 3, label_size = 22)
lg <- legend_grob(pA + theme(legend.position = "bottom", legend.text = element_text(size = 15)))
ggsave(here::here("Fig4", "Fig4.pdf"), plot_grid(g, lg, ncol = 1, rel_heights = c(1, 0.045)), width = 21, height = 19.5, device = PDF_DEVICE)
message("wrote Fig4/Fig4.pdf")
