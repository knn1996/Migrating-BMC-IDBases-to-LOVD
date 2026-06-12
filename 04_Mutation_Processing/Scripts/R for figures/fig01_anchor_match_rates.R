library(dplyr)
library(ggplot2)

infile <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step2_RefCheck/gene_confidence.csv"

df <- read.csv(infile, stringsAsFactors = FALSE)

dat <- df %>%
  filter(status %in% c("ok", "below_threshold"),
         !is.na(wilson_ci_low), n_in_test > 0) %>%
  mutate(pass = ifelse(status == "ok",
                       "\u226590% (accepted, Track A)",
                       "<90% (fell through to Track B)")) %>%
  arrange(match_rate, n_in_test) %>%
  mutate(gene = factor(gene, levels = gene))

ann <- dat %>%
  filter(status == "below_threshold") %>%
  mutate(lab = sprintf("p = %.0e", binomial_p))

p <- ggplot(dat, aes(x = match_rate * 100, y = gene, color = pass)) +
  geom_vline(xintercept = 25, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey30") +
  geom_errorbarh(aes(xmin = wilson_ci_low * 100, xmax = wilson_ci_high * 100),
                 height = 0, linewidth = 0.4) +
  geom_point(aes(size = n_in_test)) +
  geom_text(data = ann, aes(x = wilson_ci_low * 100 - 2, label = lab),
            hjust = 1, size = 2.4, show.legend = FALSE) +
  scale_color_manual(values = c("\u226590% (accepted, Track A)" = "#2C6E91",
                                "<90% (fell through to Track B)" = "#C0392B")) +
  scale_size_continuous(trans = "log10", range = c(0.5, 4.5),
                        breaks = c(1, 10, 100, 500),
                        name = "Anchorable\nsubstitutions (n)") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0.10, 0.02))) +
  labs(x = "Anchor REF-base match rate (%)", y = NULL, color = NULL,
       caption = paste("Per-gene match rate with 95% Wilson CI; point size = number of anchorable substitutions tested.",
                       "Dashed = 90% acceptance threshold; dotted = 25% chance baseline.")) +
  guides(color = guide_legend(override.aes = list(size = 2.5), order = 1),
         size = guide_legend(order = 2)) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top",
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7, color = "grey40"),
        axis.text.y = element_text(size = 5.5))

ggsave("C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures/fig01_anchor_match_rates.pdf", p, width = 6.5, height = 11, units = "in")
ggsave("C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures/fig01_anchor_match_rates.png", p, width = 6.5, height = 11, units = "in", dpi = 300)