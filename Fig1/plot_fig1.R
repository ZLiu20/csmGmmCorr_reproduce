here::i_am("Fig1/plot_fig1.R")
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
source(here::here("SupportingCode", "plot_style.R"))

K <- 4L
H <- as.matrix(expand.grid(rep(list(c(-1, 0, 1)), K)))
S <- rowSums(H != 0)

exact_lfdr <- function(z, rho, target, t = 1L,
                       probs = c(.6, .2, .1, .06, .04),
                       magnitude = function(s) rep(2, length(s))) {
  z <- matrix(z, ncol = K)
  Sigma <- matrix(rho, K, K); diag(Sigma) <- 1
  precision <- solve(Sigma)
  means <- H * magnitude(S)
  weight <- probs[S + 1L] / choose(K, S) / 2^S
  logd <- vapply(seq_len(nrow(H)), function(i) {
    residual <- sweep(z, 2L, means[i, ], "-")
    -0.5 * rowSums((residual %*% precision) * residual) + log(weight[i])
  }, numeric(nrow(z)))
  logd <- sweep(logd, 1L, apply(logd, 1L, max), "-")
  density <- exp(logd)
  null <- if (target == "pleiotropy") S <= t else !(S > t & abs(rowSums(H)) == S)
  rowSums(density[, null, drop = FALSE]) / rowSums(density)
}

x <- seq(0, 14, length.out = 561L)
ray <- function(base, direction) t(vapply(x, function(s) base + s * direction, numeric(K)))
numbers <- c(
  sprintf("Example 1: pleiotropy lfdr at (12,0,0,0): rho=0.3 %.4f; Sigma=I %.2f",
          exact_lfdr(c(12,0,0,0), .3, "pleiotropy"), exact_lfdr(c(12,0,0,0), 0, "pleiotropy")),
  sprintf("Example 1: lfdr(12,3,0,0) %.4f", exact_lfdr(c(12,3,0,0), .3, "pleiotropy")),
  sprintf("Proposition 2: replication lfdr at z1=4 and 12: %.2f, %.3f",
          exact_lfdr(c(4,1.5,1.5,1.5), .3, "replication"),
          exact_lfdr(c(12,1.5,1.5,1.5), .3, "replication"))
)
writeLines(numbers, here::here("Fig1", "Fig1_numbers.txt"))

settings <- c("rho = 0.3", "Sigma = I", "rho = 0.3, along (s,s,s,s), s = z1/2")
make_curve <- function(z, target, label) {
  data.frame(x = x, y = exact_lfdr(z, label, target),
             setting = if (label == .3) settings[1L] else settings[2L])
}
A <- rbind(make_curve(ray(c(0,0,0,0), c(1,0,0,0)), "pleiotropy", .3),
           make_curve(ray(c(0,0,0,0), c(1,0,0,0)), "pleiotropy", 0))
B <- rbind(make_curve(ray(c(0,1.5,1.5,1.5), c(1,0,0,0)), "pleiotropy", .3),
           make_curve(ray(c(0,1.5,1.5,1.5), c(1,0,0,0)), "pleiotropy", 0))
C <- rbind(make_curve(ray(c(0,1.5,1.5,1.5), c(1,0,0,0)), "replication", .3),
           make_curve(ray(c(0,1.5,1.5,1.5), c(1,0,0,0)), "replication", 0),
           transform(data.frame(x=x, y=exact_lfdr(cbind(x/2,x/2,x/2,x/2), .3, "replication")), setting=settings[3L]))

draw <- function(d, title) ggplot(d, aes(x, y, color = setting, linetype = setting)) +
  geom_line(linewidth = 1.15) + coord_cartesian(xlim = c(0,14), ylim = c(0,1)) +
  scale_color_manual(values = setNames(c("#D85C41", "#367B91", "#54A24B"), settings)) +
  scale_linetype_manual(values = setNames(c("solid", "22", "12"), settings)) +
  xlab(expression(z[1])) + ylab("lfdr") + ggtitle(title) + theme_cowplot(16)
p <- plot_grid(draw(A, "Pleiotropy: axis ray") + theme(legend.position="none"),
               draw(B, "Pleiotropy: shifted ray") + theme(legend.position="none"),
               draw(C, "Replication") + theme(legend.position=c(.35,.22)), nrow = 1,
               labels = c("A", "B", "C"))
ggsave(here::here("Fig1", "Fig1.pdf"), p, width = 18, height = 5.6, device = PDF_DEVICE)
