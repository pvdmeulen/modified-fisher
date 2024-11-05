# TESTING MODIFIED FISHER EXACT TEST ==========================================

library(dplyr)
library(ggplot2)

source("R/modified_fisher.R")

u5_12_v7_11 <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1)

u13_41_v6_47 <- modified_fisher_exact_test(u = 13, m = 41, v = 6, n = 47, odds_ratio = 1)

u37_65_v28_71 <- modified_fisher_exact_test(u = 37, m = 65, v = 28, n = 71, odds_ratio = 1)


require(ggplot2)
require(latex2exp)


u <- 37
m <- 65

v <- 28
n <- 71

odds_ratio <- 1
alpha <- 0.05


plot_subtitle <- paste0("$m = ", m, ", n = ", n, ", H_0",
                        " : ", odds_ratio, "$")

ggplot(u37_65_v28_71$local.size.data, aes(x = pi1, y = size)) +
  geom_line(aes(colour = method), linewidth = 1) +
  #scale_y_continuous(limits = c(-4, 4*alpha*100)) +
  geom_hline(yintercept = alpha*100, linetype = 2, colour = "grey30",
             linewidth = 0.4) +
  coord_cartesian(ylim = c(0, 1.05*alpha*100)) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Size of Modified Fisher Exact Test",
    subtitle = latex2exp::TeX(paste0("\\textit{$m=", m, ",\\,n=", n, ",\\,H_{0}=",
                                     odds_ratio, "$}")),
    y = "Size (%)",
    x = bquote(pi[.("1")])
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(family = "STIX Two Text", face = "italic"),
    axis.title.x = element_text(family = "STIX Two Text", face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

ggsave(last_plot(), width = 160, height = 120, units = "mm", dpi = "retina",
       filename = paste0("Example/size_plot_or", odds_ratio, "_m", m, "_n", n,
                         ".png"))
