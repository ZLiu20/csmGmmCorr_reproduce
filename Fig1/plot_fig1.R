# Figure 1 and the numerical statements of Example 1 and Section 4.1 (exact computations; no
# simulation). Model of Example 1: K = 4, one magnitude component per support with all
# magnitudes 2, support-size probabilities (0.6, 0.2, 0.1, 0.06, 0.04), spread uniformly
# within sizes and signs, and equicorrelation 0.3.
#   Rscript Fig1/plot_fig1.R      (writes Fig1/Fig1.pdf and Fig1/Fig1_numbers.txt)
here::i_am("Fig1/plot_fig1.R")
suppressMessages({ library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "plot_style.R"))

K <- 4
H <- as.matrix(expand.grid(rep(list(c(-1, 0, 1)), K)))
S <- rowSums(H != 0)
# exact lfdr: magnitude of each configuration given by mag(size); weights probs[size]/(choose(K,s) 2^s)
lfdr <- function(Z, rho, target, t = 1, probs = c(0.6, 0.2, 0.1, 0.06, 0.04), mag = function(s) rep(2, length(s))) {
  Z <- matrix(Z, ncol = K)
  Sig <- matrix(rho, K, K); diag(Sig) <- 1; Om <- solve(Sig)
  M <- H * mag(S)
  w <- probs[S + 1] / choose(K, S) / 2^S
  ll <- sapply(seq_len(nrow(H)), function(r) {
    D <- sweep(Z, 2, M[r, ]); -0.5 * rowSums((D %*% Om) * D) + log(w[r]) })
  ll <- matrix(ll, nrow = nrow(Z))
  p <- exp(ll - apply(ll, 1, max))
  null <- if (target == "pleiotropy") S <= t else !(S > t & abs(rowSums(H)) == S)
  rowSums(p[, null, drop = FALSE]) / rowSums(p)
}
x <- seq(0, 14, length.out = 561)
ray <- function(base, dirn) t(sapply(x, function(s) base + s * dirn))

out <- c(
  sprintf("Example 1: pleiotropy lfdr (t = 1) at (12,0,0,0): rho = 0.3: %.4f; Sigma = I: %.2f",
          lfdr(c(12, 0, 0, 0), 0.3, "pleiotropy"), lfdr(c(12, 0, 0, 0), 0, "pleiotropy")),
  sprintf("Example 1: lfdr(12,3,0,0) = %.4f", lfdr(c(12, 3, 0, 0), 0.3, "pleiotropy")))
B <- lfdr(ray(c(0, 1.5, 1.5, 1.5), c(1, 0, 0, 0)), 0.3, "pleiotropy")
inc <- diff(B) > 0
# turning points: sign changes of the slope on the grid, refined by root finding
slope <- function(z) (lfdr(c(z + 1e-6, 1.5, 1.5, 1.5), 0.3, "pleiotropy") - lfdr(c(z - 1e-6, 1.5, 1.5, 1.5), 0.3, "pleiotropy")) / 2e-6
ch <- which(diff(inc) != 0)
tp <- sapply(ch, function(i) uniroot(slope, c(x[i], x[i + 2]))$root)
bounds <- c(if (inc[1]) 0, tp, if (inc[length(inc)]) max(x))
iv <- matrix(bounds, ncol = 2, byrow = TRUE)
out <- c(out, paste("Example 1: the lfdr along (z1,1.5,1.5,1.5) increases on",
                    paste(sprintf("(%.2f, %.2f)", iv[, 1], iv[, 2]), collapse = " and ")),
         sprintf("Proposition 2: replication lfdr along (z1,1.5,1.5,1.5) at z1 = 4 and 12: %.2f, %.3f",
                 lfdr(c(4, 1.5, 1.5, 1.5), 0.3, "replication"), lfdr(c(12, 1.5, 1.5, 1.5), 0.3, "replication")))
# Section 4.1: incongruity without correlation (Sigma = I), magnitudes 4 for supports of sizes
# one, three and four and 1 for size two
mag2 <- function(s) ifelse(s == 2, 1, 4)
out <- c(out, sprintf("Section 4.1: Sigma = I, lfdr along (z1,1,0,0) at z1 = 2 and 4: %.2f, %.2f",
                      lfdr(c(2, 1, 0, 0), 0, "pleiotropy", mag = mag2), lfdr(c(4, 1, 0, 0), 0, "pleiotropy", mag = mag2)))
writeLines(out); writeLines(out, here::here("Fig1", "Fig1_numbers.txt"))

c_cor <- grDevices::hcl(15, 100, 65); c_ind <- grDevices::hcl(195, 100, 65); c_diag <- grDevices::hcl(105, 100, 55)
lv <- c("rho = 0.3", "Sigma = I", "rho = 0.3, along (s,s,s,s), s = z1/2")
labs_ <- c(expression(rho == 0.3), expression(Sigma == I), expression(paste(rho == 0.3, ", along (s, s, s, s), ", s == z[1] / 2)))
mk <- function(df, ylab, title, shade = NULL) {
  df$setting <- factor(df$setting, levels = lv)
  p <- ggplot(df, aes(x, y, color = setting, linetype = setting))
  if (!is.null(shade)) p <- p + annotate("rect", xmin = shade[, 1], xmax = shade[, 2], ymin = -Inf, ymax = Inf, fill = "grey90")
  p + geom_line(linewidth = 1.2) +
    scale_color_manual(values = setNames(c(c_cor, c_ind, c_diag), lv), labels = function(b) labs_[match(b, lv)], name = NULL) +
    scale_linetype_manual(values = setNames(c("solid", "22", "12"), lv), labels = function(b) labs_[match(b, lv)], name = NULL) +
    coord_cartesian(xlim = c(0, 14), ylim = c(0, 1)) + xlab(expression(z[1])) + ylab(ylab) + ggtitle(title) +
    theme_cowplot() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
                            legend.text = element_text(size = 12), legend.position = c(0.55, 0.5), plot.title = element_text(size = 16, face = "plain"))
}
A <- rbind(data.frame(x = x, y = lfdr(ray(c(0, 0, 0, 0), c(1, 0, 0, 0)), 0.3, "pleiotropy"), setting = lv[1]),
           data.frame(x = x, y = lfdr(ray(c(0, 0, 0, 0), c(1, 0, 0, 0)), 0, "pleiotropy"), setting = lv[2]))
Bd <- rbind(data.frame(x = x, y = B, setting = lv[1]),
            data.frame(x = x, y = lfdr(ray(c(0, 1.5, 1.5, 1.5), c(1, 0, 0, 0)), 0, "pleiotropy"), setting = lv[2]))
C <- rbind(data.frame(x = x, y = lfdr(ray(c(0, 1.5, 1.5, 1.5), c(1, 0, 0, 0)), 0.3, "replication"), setting = lv[1]),
           data.frame(x = x, y = lfdr(ray(c(0, 1.5, 1.5, 1.5), c(1, 0, 0, 0)), 0, "replication"), setting = lv[2]),
           data.frame(x = x, y = lfdr(cbind(x / 2, x / 2, x / 2, x / 2), 0.3, "replication"), setting = lv[3]))
p <- plot_grid(mk(A, "Pleiotropy lfdr (t = 1)", "A. Along (z1, 0, 0, 0)"),
               mk(Bd, "Pleiotropy lfdr (t = 1)", "B. Along (z1, 1.5, 1.5, 1.5)", shade = iv) + theme(legend.position = "none"),
               mk(C, "Replication lfdr (t = 1)", "C. Replication; lines 1-2 along (z1, 1.5, 1.5, 1.5)"), nrow = 1)
ggsave(here::here("Fig1", "Fig1.pdf"), p, width = 18, height = 5.6)
message("wrote Fig1/Fig1.pdf")
