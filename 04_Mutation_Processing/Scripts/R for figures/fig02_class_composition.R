library(dplyr)
library(ggplot2)

dir <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step8_Merging"

ded <- read.delim(file.path(dir, "dedup_merged_variants.tsv"), stringsAsFactors = FALSE)

lev <- c("substitution", "deletion", "duplication", "insertion")

comp <- ded %>%
  filter(mut_type %in% lev) %>%
  mutate(patient_count = as.integer(patient_count)) %>%
  group_by(mut_type) %>%
  summarise(n = n(), recur = mean(patient_count > 1) * 100, .groups = "drop") %>%
  mutate(prop = n / sum(n),
         mut_type = factor(mut_type, levels = lev)) %>%
  arrange(mut_type) %>%
  mutate(xlab = sprintf("%s\n(%.0f%% in \u22652 patients)", mut_type, recur))

comp$xlab <- factor(comp$xlab, levels = comp$xlab)

p <- ggplot(comp, aes(x = xlab, y = prop * 100)) +
  geom_col(width = 0.62, fill = "#2C6E91") +
  geom_text(aes(label = sprintf("%.1f%%\n(n = %d)", prop * 100, n)),
            vjust = -0.3, size = 2.9, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 75), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Share of distinct variants (%)",
       caption = paste("Composition of the 2,240 distinct variants (post-normalization classes).",
                       "Sub-labels: share of each class seen 
in >1 patient, declining monotonically.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7, color = "grey40"))

ggsave("C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures/fig02_class_composition.pdf", p, width = 5.5, height = 4, units = "in")
ggsave("C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures/fig02_class_composition.png", p, width = 5.5, height = 4, units = "in", dpi = 300)