library(ggplot2)
library(dplyr)
library(patchwork)

DEDUP_TSV <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step8_Merging/dedup_merged_variants.tsv"

class_cols <- c(substitution = "#2C6E91", deletion = "#5A93AC",
                duplication = "#9CB7C5", insertion = "#C0392B")

d <- read.delim(DEDUP_TSV, stringsAsFactors = FALSE)

n_distinct <- nrow(d)
n_pairs    <- sum(d$patient_count, na.rm = TRUE)

bins <- c(0, 1, 2, 3, 5, 10, 20, 50, Inf)
blab <- c("1", "2", "3", "4-5", "6-10", "11-20", "21-50", ">50")
d$bin <- cut(d$patient_count, breaks = bins, labels = blab, right = TRUE)
binct <- d %>% count(bin, .drop = FALSE)
pct1 <- round(binct$n[binct$bin == "1"] / n_distinct * 100, 1)

pA <- ggplot(binct, aes(x = bin, y = n)) +
  geom_col(fill = "#2C6E91", width = 0.78) +
  geom_text(aes(label = n), vjust = -0.35, size = 2.7, color = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Patients per variant", y = "Distinct variants",
       title = "A  Recurrence distribution") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

top <- d %>%
  arrange(desc(patient_count)) %>%
  slice(1:10) %>%
  mutate(lab = paste0(gene, "  ", dedup_key),
         lab = factor(lab, levels = rev(lab)),
         mut_type = factor(mut_type, levels = names(class_cols)))

pB <- ggplot(top, aes(x = patient_count, y = lab, fill = mut_type)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = patient_count), hjust = -0.25, size = 2.7, color = "grey20") +
  scale_fill_manual(values = class_cols, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = "Patients", y = NULL, title = "B  Most recurrent variants") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7.5),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

p <- pA + pB + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(
    caption = paste0("Patient recurrence across ", format(n_distinct, big.mark = ","),
                     " distinct variants (", format(n_pairs, big.mark = ","),
                     " patient-variant pairs). ", pct1,
                     "% are private (one patient); a small tail of founder/recurrent\n",
                     "mutations dominates, led by SBDS and AIRE. Panel B coloured by variant class."),
    theme = theme(plot.caption = element_text(hjust = 0.5, size = 7, color = "grey40")))

outdir <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures"
ggsave(file.path(outdir, "fig05_recurrence.pdf"), p, width = 8, height = 4, units = "in")
ggsave(file.path(outdir, "fig05_recurrence.png"), p, width = 8, height = 4, units = "in", dpi = 300)
