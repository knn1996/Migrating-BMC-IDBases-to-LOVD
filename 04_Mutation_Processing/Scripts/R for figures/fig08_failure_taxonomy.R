library(ggplot2)
library(dplyr)

DISP_TSV <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step8_Merging/unresolved_disposition.tsv"

disp_cols <- c(RESCUABLE_AUTO = "#2C6E91",
               RESCUABLE_MANUAL = "#5A93AC",
               UNRESCUABLE_DATA_ERROR = "#C0392B",
               UNRESCUABLE = "#7B241C",
               UNKNOWN = "#B0B0B0")
disp_labs <- c(RESCUABLE_AUTO = "Rescuable (automated)",
               RESCUABLE_MANUAL = "Rescuable (manual)",
               UNRESCUABLE_DATA_ERROR = "Unrescuable (data error)",
               UNRESCUABLE = "Unrescuable",
               UNKNOWN = "Unknown")

d <- read.delim(DISP_TSV, stringsAsFactors = FALSE)
n_total <- nrow(d)
n_resc  <- sum(grepl("^RESCUABLE", d$disposition))
pct_resc <- round(n_resc / n_total * 100, 1)

agg <- d %>%
  count(category, disposition, name = "n") %>%
  mutate(disposition = factor(disposition, levels = names(disp_cols)))

cat_order <- agg %>% group_by(category) %>% summarise(t = sum(n)) %>% arrange(t) %>% pull(category)
agg$category <- factor(agg$category, levels = cat_order)

cat_tot <- agg %>% group_by(category) %>% summarise(t = sum(n), .groups = "drop")

p <- ggplot(agg, aes(x = n, y = category, fill = disposition)) +
  geom_col(width = 0.7) +
  geom_text(data = cat_tot, aes(x = t, y = category, label = t, fill = NULL),
            hjust = -0.25, size = 2.8, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = disp_cols, labels = disp_labs, name = NULL, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Distinct unresolved variants", y = NULL,
       title = "Unresolved-variant failure taxonomy",
       caption = paste0(format(n_total, big.mark = ","),
                        " distinct unresolved variants by Mutalyzer error category and disposition. ",
                        pct_resc, "% are rescuable\n",
                        "(transcript-selector or notation issues); the remainder are genuine ",
                        "reference/data errors in the source curation.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 11, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7, color = "grey40"),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.margin = margin(5, 12, 8, 5, "mm")) +
  guides(fill = guide_legend(nrow = 2))

outdir <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures"
ggsave(file.path(outdir, "fig08_failure_taxonomy.pdf"), p, width = 6.5, height = 4, units = "in")
ggsave(file.path(outdir, "fig08_failure_taxonomy.png"), p, width = 6.5, height = 4, units = "in", dpi = 300)
