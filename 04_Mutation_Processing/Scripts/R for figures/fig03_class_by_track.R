library(dplyr)
library(ggplot2)

dir <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step8_Merging"
ded <- read.delim(file.path(dir, "dedup_merged_variants.tsv"), stringsAsFactors = FALSE)

track_map <- c(NG_IDRefseq = "Track A", NM_MANE = "Track B", NM_IDRefseq = "Track C")
class_lev <- c("substitution", "deletion", "duplication", "insertion")
class_cols <- c(substitution = "blue", deletion = "#5A93AC",
                duplication = "#9CB7C5", insertion = "#C0392B")

df <- ded %>%
  filter(mut_type %in% class_lev) %>%
  mutate(track = factor(track_map[source_track], levels = c("Track A", "Track B", "Track C")),
         mut_type = factor(mut_type, levels = class_lev)) %>%
  count(track, mut_type, name = "n")

totals <- df %>% group_by(track) %>% summarise(total = sum(n), .groups = "drop")

df <- df %>% left_join(totals, by = "track")

p <- ggplot(df, aes(x = track, y = n, fill = mut_type)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(n >= 30, n, "")),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "white", fontface = "bold") +
  geom_text(data = totals, aes(x = track, y = total, label = total, fill = NULL),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = class_cols, name = "Variant class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Distinct variants",
       caption = paste0("Variant-class composition by resolution track. ",
              "Track A: offset-validated (NG); Track B: MANE Select (NM);\n",
              "Track C: IDRefSeq NM. Total = ", sum(totals$total), ".")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7, color = "grey40"),
        plot.margin = margin(5, 5, 10, 5, "mm"),
        legend.position = "right")

outdir <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures"
ggsave(file.path(outdir, "fig03_class_by_track.pdf"), p, width = 5.5, height = 4, units = "in")
ggsave(file.path(outdir, "fig03_class_by_track.png"), p, width = 5.5, height = 4, units = "in", dpi = 300)
