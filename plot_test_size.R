# TESTING MODIFIED FISHER EXACT TEST ==========================================

library(dplyr)
library(ggplot2)
library(stringr)

library(ggplot2)
library(latex2exp)

source("R/modified_fisher.R")

alpha <- 0.05

# Test u = 5, v = 7, m = 12, n = 11 ===========================================

testing_set <- tibble(
  u = rep(5, 5),
  v = rep(7, 5),
  m = rep(12, 5),
  n = rep(11, 5),
  or = c(0.1, 0.5, 1, 2, 10)
  )

testing_set <- bind_rows(
  testing_set %>% mutate(method = "trust"),
  testing_set %>% mutate(method = "zoom")
)

# Store test results in one column, and local size data in another for plots:

testing_set <- testing_set %>%
  rowwise() %>%
  mutate(
    results = list(modified_fisher_exact_test(u = u, m = m, v = v, n = n,
                                              odds_ratio = or, alpha = alpha,
                                              method = method)),
    localsize = list(results[["local.size.data"]])
  ) %>%
  ungroup()

# Create plots:

for(chosen_or in testing_set$or){

  df <- testing_set %>%
    filter(or == chosen_or)

  u <- df %>% pull(u) %>% unique()
  v <- df %>% pull(v) %>% unique()
  m <- df %>% pull(m) %>% unique()
  n <- df %>% pull(n) %>% unique()

  plotdata <- df %>%
    pull(localsize) %>%
    bind_rows() %>%
    as_tibble()

  # Plot:

  plot_subtitle <- paste0("$m = ", m, ", n = ", n, ", H_0",
                          " : ", chosen_or, "$")

  file_name <- paste0("Example/size_plot_or", chosen_or, "_m", m, "_n",
                      n, ".png")

  plot <- ggplot(data = plotdata, aes(x = pi1, y = size)) +
    geom_line(aes(colour = method), linewidth = 1) +
    #scale_y_continuous(limits = c(-4, 4*alpha*100)) +
    geom_hline(yintercept = alpha*100, linetype = 2, colour = "grey30",
               linewidth = 0.4) +
    coord_cartesian(ylim = c(0, 1.05*alpha*100)) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = "Size of Modified Fisher Exact Test",
      subtitle = latex2exp::TeX(paste0("\\textit{$m=", m, ",\\,n=", n,
                                       ",\\,H_{0}=", chosen_or, "$}")),
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

  ggsave(plot, width = 160, height = 120, units = "mm", dpi = "retina",
         filename = file_name)

}

# Test u = 5, v = 7, m = 12, n = 11 ===========================================

# Takes a while!

testing_set2 <- tibble(
  u = rep(13, 3),
  v = rep(6, 3),
  m = rep(41, 3),
  n = rep(47, 3),
  or = c(0.5, 1, 2)
)

# Store test results in one column, and local size data in another for plots:

testing_set2 <- testing_set2 %>%
  rowwise() %>%
  mutate(
    results = list(modified_fisher_exact_test(u = u, m = m, v = v, n = n,
                                              odds_ratio = or, alpha = alpha)),
    localsize = list(results[["local.size.data"]])
  ) %>%
  ungroup()

# Create plots:

for(chosen_or in testing_set2$or){

  df <- testing_set2 %>%
    filter(or == chosen_or)

  u <- df %>% pull(u)
  v <- df %>% pull(v)
  m <- df %>% pull(m)
  n <- df %>% pull(n)

  plotdata <- df %>%
    pull(localsize) %>%
    as.data.frame()

  # Plot:

  plot_subtitle <- paste0("$m = ", m, ", n = ", n, ", H_0",
                          " : ", chosen_or, "$")

  file_name <- paste0("Example/size_plot_or", chosen_or, "_m", m, "_n",
                      n, ".png")

  plot <- ggplot(data = plotdata, aes(x = pi1, y = size)) +
    geom_line(aes(colour = method), linewidth = 1) +
    #scale_y_continuous(limits = c(-4, 4*alpha*100)) +
    geom_hline(yintercept = alpha*100, linetype = 2, colour = "grey30",
               linewidth = 0.4) +
    coord_cartesian(ylim = c(0, 1.05*alpha*100)) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = "Size of Modified Fisher Exact Test",
      subtitle = latex2exp::TeX(paste0("\\textit{$m=", m, ",\\,n=", n,
                                       ",\\,H_{0}=", chosen_or, "$}")),
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

  ggsave(plot, width = 160, height = 120, units = "mm", dpi = "retina",
         filename = file_name)

}

# Test u = 37, v = 28, m = 65, n = 71 ===========================================

# Takes even longer...

testing_set3 <- tibble(
  u = rep(37, 3),
  v = rep(28, 3),
  m = rep(65, 3),
  n = rep(71, 3),
  or = c(0.5, 1, 2)
)

# Store test results in one column, and local size data in another for plots:

testing_set3 <- testing_set3 %>%
  rowwise() %>%
  mutate(
    results = list(modified_fisher_exact_test(u = u, m = m, v = v, n = n,
                                              odds_ratio = or, alpha = alpha)),
    localsize = list(results[["local.size.data"]])
  ) %>%
  ungroup()

# Create plots:

for(chosen_or in testing_set3$or){

  df <- testing_set3 %>%
    filter(or == chosen_or)

  u <- df %>% pull(u)
  v <- df %>% pull(v)
  m <- df %>% pull(m)
  n <- df %>% pull(n)

  plotdata <- df %>%
    pull(localsize) %>%
    as.data.frame()

  # Plot:

  plot_subtitle <- paste0("$m = ", m, ", n = ", n, ", H_0",
                          " : ", chosen_or, "$")

  file_name <- paste0("Example/size_plot_or", chosen_or, "_m", m, "_n",
                      n, ".png")

  plot <- ggplot(data = plotdata, aes(x = pi1, y = size)) +
    geom_line(aes(colour = method), linewidth = 1) +
    #scale_y_continuous(limits = c(-4, 4*alpha*100)) +
    geom_hline(yintercept = alpha*100, linetype = 2, colour = "grey30",
               linewidth = 0.4) +
    coord_cartesian(ylim = c(0, 1.05*alpha*100)) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = "Size of Modified Fisher Exact Test",
      subtitle = latex2exp::TeX(paste0("\\textit{$m=", m, ",\\,n=", n,
                                       ",\\,H_{0}=", chosen_or, "$}")),
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

  ggsave(plot, width = 160, height = 120, units = "mm", dpi = "retina",
         filename = file_name)

}
