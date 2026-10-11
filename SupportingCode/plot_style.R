.hue <- function(h) grDevices::hcl(h = h, c = 100, l = 65)

METHOD_STYLE <- data.frame(
  method = c("penalized", "ML", "QCH", "adaFilter", "DEI", "MediationTest",
             "PCPathSafeguard", "sym", "Correlated csmGmm",
             "Correlated csmGmm (ML)", "QCH", "adaFilter", "DEI",
             "MediationTest", "S-PC", "Correlated csmGmm, symmetric statistics"),
  label = c("Correlated csmGmm", "Correlated csmGmm (ML)", "QCH", "adaFilter",
            "DEI", "MediationTest", "S-PC", "Correlated csmGmm, symmetric statistics",
            "Correlated csmGmm", "Correlated csmGmm (ML)", "QCH", "adaFilter",
            "DEI", "MediationTest", "S-PC", "Correlated csmGmm, symmetric statistics"),
  color = .hue(c(15, 195, 105, 285, 60, 240, 330, 150,
                 15, 195, 105, 285, 60, 240, 330, 150)),
  linetype = c("solid", "22", "13", "1343", "44", "121412", "6222", "solid",
               "solid", "22", "13", "1343", "44", "121412", "6222", "solid"),
  stringsAsFactors = FALSE
)

method_scales <- function(methods) {
  style <- METHOD_STYLE[match(methods, METHOD_STYLE$method), , drop = FALSE]
  list(
    ggplot2::scale_color_manual(values = stats::setNames(style$color, style$label),
                                breaks = style$label, name = "Method"),
    ggplot2::scale_linetype_manual(values = stats::setNames(style$linetype, style$label),
                                   breaks = style$label, name = "Method")
  )
}

TARGET_NAME <- c(pleiotropy = "Pleiotropy", replication = "Replication", mediation = "Mediation")
EFFECT_XLAB <- c(pleiotropy = "Effect Magnitude", replication = "Effect Magnitude",
                 mediation = "Min Effect Magnitude")

theme_paper <- function(base = 16) {
  cowplot::theme_cowplot() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 20),
      axis.text = ggplot2::element_text(size = base),
      legend.title = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 18),
      legend.key.size = ggplot2::unit(3, "line")
    )
}

fdp_line <- function() ggplot2::geom_hline(yintercept = 0.1, linetype = 2, color = "grey")

fdp_top <- function(x) max(0.5, ceiling(max(x, na.rm = TRUE) * 10.5) / 10)

sim_panel <- function(z, y, methods, xlab, ylab, ylim, hline = FALSE, sym = NULL) {
  z <- as.data.frame(z)
  z <- z[z$method %in% methods, , drop = FALSE]
  z$Method <- METHOD_STYLE$label[match(z$method, METHOD_STYLE$method)]
  used <- methods
  panel <- ggplot2::ggplot(z, ggplot2::aes(x = x, y = .data[[y]], group = Method)) +
    ggplot2::geom_line(ggplot2::aes(color = Method, linetype = Method), linewidth = 1.2)
  if (!is.null(sym) && nrow(sym)) {
    sym <- as.data.frame(sym)
    sym$Method <- METHOD_STYLE$label[match(sym$method, METHOD_STYLE$method)]
    panel <- panel + ggplot2::geom_line(data = sym,
      ggplot2::aes(color = Method, linetype = Method), linewidth = 1.2) +
      ggplot2::geom_point(data = sym, ggplot2::aes(color = Method), size = 3,
                          show.legend = FALSE)
    used <- c(methods, "sym")
  }
  if (hline) panel <- panel + fdp_line()
  panel + method_scales(used) +
    ggplot2::scale_x_continuous(breaks = sort(unique(z$x))) +
    ggplot2::coord_cartesian(ylim = ylim) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
    theme_paper()
}

legend_grob <- function(plot) {
  grob <- ggplot2::ggplotGrob(plot)
  index <- grep("^guide-box", grob$layout$name)
  index <- index[!vapply(grob$grobs[index], function(x) inherits(x, "zeroGrob"), logical(1))]
  grob$grobs[[index[1L]]]
}

combine_panels <- function(panels, nrow, legend = 1L, rel = 0.08,
                           label_size = 22, ncol_legend = 4L) {
  body <- cowplot::plot_grid(
    plotlist = lapply(panels, function(p) p + ggplot2::theme(legend.position = "none")),
    labels = LETTERS[seq_along(panels)], nrow = nrow, label_size = label_size
  )
  source <- if (inherits(legend, "ggplot")) legend else panels[[legend]]
  guide <- legend_grob(source + ggplot2::guides(
    color = ggplot2::guide_legend(ncol = ncol_legend, byrow = TRUE),
    linetype = ggplot2::guide_legend(ncol = ncol_legend, byrow = TRUE)
  ) + ggplot2::theme(legend.position = "bottom", legend.justification = "center"))
  cowplot::plot_grid(body, guide, ncol = 1, rel_heights = c(1, rel))
}

legend_plot <- function(methods, points = "sym") {
  z <- data.frame(x = rep(1:2, length(methods)), y = 1,
                  method = rep(methods, each = 2L))
  z$Method <- METHOD_STYLE$label[match(z$method, METHOD_STYLE$method)]
  ggplot2::ggplot(z, ggplot2::aes(x, y, group = Method)) +
    ggplot2::geom_line(ggplot2::aes(color = Method, linetype = Method), linewidth = 1.2) +
    method_scales(methods) + theme_paper()
}

PDF_DEVICE <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
SIGMA <- if (isTRUE(capabilities("cairo"))) "Σ" else "Sigma"
