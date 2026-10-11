# Common look of the figures (as in Sun, McCaw and Lin, 2025): ggplot2 with theme_cowplot,
# lines of width 1.2, hue colors, a grey dashed line at q = 0.1. Each method keeps the same
# color and line type in every figure of the paper and the supplement.
hue <- function(h) grDevices::hcl(h = h, c = 100, l = 65)
METHOD_STYLE <- data.frame(
  method   = c("penalized", "ML", "QCH", "adaFilter", "DEI", "MediationTest", "PCPathSafeguard", "sym"),
  label    = c("Correlated csmGmm", "Correlated csmGmm (ML)", "QCH", "adaFilter", "DEI", "MediationTest",
               "JS-PC", "Correlated csmGmm, symmetric statistics"),
  color    = hue(c(15, 195, 105, 285, 60, 240, 330, 150)),
  linetype = c("solid", "22", "13", "1343", "44", "121412", "6222", "solid"),
  stringsAsFactors = FALSE)
method_scales <- function(methods) {
  s <- METHOD_STYLE[match(methods, METHOD_STYLE$method), ]
  list(ggplot2::scale_color_manual(values = stats::setNames(s$color, s$label), breaks = s$label, name = "Method"),
       ggplot2::scale_linetype_manual(values = stats::setNames(s$linetype, s$label), breaks = s$label, name = "Method"))
}
TARGET_NAME <- c(pleiotropy = "Pleiotropy", replication = "Replication", mediation = "Mediation")
EFFECT_XLAB <- c(pleiotropy = "Effect Magnitude", replication = "Effect Magnitude", mediation = "Min Effect Magnitude")

theme_paper <- function(base = 16) {
  cowplot::theme_cowplot() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 20), axis.text = ggplot2::element_text(size = base),
                   legend.title = ggplot2::element_text(size = 20), legend.text = ggplot2::element_text(size = 18),
                   legend.key.size = ggplot2::unit(3, "line"))
}
fdp_line <- function() ggplot2::geom_hline(yintercept = 0.1, linetype = 2, color = "grey")
# Upper limit of the FDP axis: at least 0.5, otherwise the next multiple of 0.1 above the data.
fdp_top <- function(x) max(0.5, ceiling(max(x, na.rm = TRUE) * 1.05 * 10) / 10)

# One panel: y against the design value x for each method (penalized fit with sign-symmetric
# statistics added as a green line with points when sym is given).
sim_panel <- function(z, y, methods, xlab, ylab, ylim, hline = FALSE, sym = NULL) {
  z <- as.data.frame(z); z <- z[z$method %in% methods, ]
  z$Method <- METHOD_STYLE$label[match(z$method, METHOD_STYLE$method)]
  meth <- methods
  p <- ggplot2::ggplot(z, ggplot2::aes(x = x, y = .data[[y]], group = Method)) +
    ggplot2::geom_line(ggplot2::aes(color = Method, linetype = Method), linewidth = 1.2)
  if (!is.null(sym) && nrow(sym)) {
    sym <- as.data.frame(sym)
    sym$Method <- METHOD_STYLE$label[METHOD_STYLE$method == "sym"]
    p <- p + ggplot2::geom_line(data = sym, ggplot2::aes(color = Method, linetype = Method), linewidth = 1.2) +
      ggplot2::geom_point(data = sym, ggplot2::aes(color = Method), size = 3, show.legend = FALSE)
    meth <- c(methods, "sym")
  }
  if (hline) p <- p + fdp_line()
  p + method_scales(meth) + ggplot2::scale_x_continuous(breaks = sort(unique(z$x))) +
    ggplot2::coord_cartesian(ylim = ylim) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + theme_paper()
}

# Legend of a plot as a grob (works with ggplot2 3.4 and 3.5, where the legend box of a
# bottom legend is named guide-box-bottom).
legend_grob <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  i <- grep("^guide-box", g$layout$name)
  i <- i[!vapply(g$grobs[i], function(x) inherits(x, "zeroGrob"), TRUE)]
  g$grobs[[i[1]]]
}
# Combine panels with labels and one horizontal legend at the bottom. legend: index of the
# panel whose legend is shown, or a ggplot that only serves to draw the legend.
combine_panels <- function(panels, nrow, legend = 1, rel = 0.08, label_size = 22, ncol_legend = 4) {
  g <- cowplot::plot_grid(plotlist = lapply(panels, function(p) p + ggplot2::theme(legend.position = "none")),
                          labels = LETTERS[seq_along(panels)], nrow = nrow, label_size = label_size)
  lp <- if (inherits(legend, "ggplot")) legend else panels[[legend]]
  lg <- legend_grob(lp + ggplot2::guides(color = ggplot2::guide_legend(ncol = ncol_legend, byrow = TRUE),
                                         linetype = ggplot2::guide_legend(ncol = ncol_legend, byrow = TRUE)) +
                      ggplot2::theme(legend.position = "bottom", legend.justification = "center"))
  cowplot::plot_grid(g, lg, ncol = 1, rel_heights = c(1, rel))
}
# A plot whose only purpose is a legend with the given methods.
legend_plot <- function(methods, points = "sym") {
  z <- data.frame(x = rep(1:2, length(methods)), y = 1, method = rep(methods, each = 2))
  z$Method <- METHOD_STYLE$label[match(z$method, METHOD_STYLE$method)]
  ggplot2::ggplot(z, ggplot2::aes(x, y, group = Method)) +
    ggplot2::geom_line(ggplot2::aes(color = Method, linetype = Method), linewidth = 1.2) +
    method_scales(methods) + theme_paper()
}
# PDF device that can draw the Greek letters used in some labels.
PDF_DEVICE <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
SIGMA <- if (isTRUE(capabilities("cairo"))) "Σ" else "Sigma"
