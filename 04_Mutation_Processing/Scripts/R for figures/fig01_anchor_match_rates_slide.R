library(dplyr)
library(ggplot2)

infile <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step2_RefCheck/gene_confidence.csv"

df <- read.csv(infile, stringsAsFactors = FALSE)

dat <- df %>%
  filter(status %in% c("ok", "below_threshold"),
         !is.na(wilson_ci_low), n_in_test > 0) %>%
  mutate(pass = ifelse(status == "ok",
                       "\u226590% (accepted, Track A)",
                       "<90% (fell through to Track B)"))

n_show_pass <- 6

sel <- bind_rows(
  filter(dat, status == "below_threshold"),
  slice_max(filter(dat, status == "ok"), order_by = n_in_test,
            n = n_show_pass, with_ties = FALSE)
) %>%
  arrange(match_rate, n_in_test) %>%
  mutate(gene = factor(gene, levels = gene))

ann <- sel %>%
  filter(status == "below_threshold") %>%
  mutate(lab = sprintf("p = %.0e", binomial_p))

p <- ggplot(sel, aes(x = match_rate * 100, y = gene, color = pass)) +
  geom_vline(xintercept = 25, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey30") +
  geom_errorbarh(aes(xmin = wilson_ci_low * 100, xmax = wilson_ci_high * 100),
                 height = 0, linewidth = 0.6) +
  geom_point(aes(size = n_in_test)) +
  geom_text(data = ann, aes(x = wilson_ci_low * 100 - 2, label = lab),
            hjust = 1, size = 3.6, show.legend = FALSE) +
  scale_color_manual(values = c("\u226590% (accepted, Track A)" = "#2C6E91",
                                "<90% (fell through to Track B)" = "#C0392B")) +
  scale_size_continuous(trans = "log10", range = c(2.5, 9),
                        breaks = c(1, 10, 100, 500),
                        name = "Anchorable\nsubstitutions (n)") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0.12, 0.03))) +
  labs(x = "Anchor REF-base match rate (%)", y = NULL, color = NULL,
       caption = paste("95% Wilson CI; point size = anchorable substitutions.")) +
  guides(color = guide_legend(override.aes = list(size = 3.5), order = 1),
         size = guide_legend(order = 2)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 9, color = "grey40"),
        axis.text.y = element_text(size = 13))

ggsave("C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures/fig01_anchor_match_rates_slide.pdf", p, width = 10, height = 5.6, units = "in")
ggsave("C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures/fig01_anchor_match_rates_slide.png", p, width = 10, height = 5.6, units = "in", dpi = 300)
