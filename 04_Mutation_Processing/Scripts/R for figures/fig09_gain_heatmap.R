library(ggplot2)

history_tsv <- file.path("..", "..", "Output", "Step10_Audit", "history", "audit_history.tsv")
out_png <- file.path("..", "..", "..", "docs", "figures", "fig09_gain_heatmap.png")
if (!file.exists(history_tsv)) {
  stop("audit_history.tsv not found at ", normalizePath(history_tsv, mustWork = FALSE),
       "; run this script from '04_Mutation_Processing/Scripts/R for figures'.")
}

h <- read.delim(history_tsv, stringsAsFactors = FALSE)
measures <- data.frame(
  column = c("in_final_pct", "protein_good_pct", "pair_consistency_pct", "excess_pairs", "pairs_coherence_pct"),
  label = c("Recovery", "Protein matches source", "Pair consistency", "Excess pairs", "c./p. coherence"),
  sign = c(1, 1, 1, -1, 1),
  stringsAsFactors = FALSE
)
plateau_runs <- 11:13
plateau_limit <- 0.03

cells <- do.call(rbind, lapply(seq_len(nrow(measures)), function(i) {
  v <- as.numeric(h[[measures$column[i]]])
  raw <- diff(v)
  total <- (v[length(v)] - v[1]) * measures$sign[i]
  data.frame(measure = measures$label[i], run = seq_along(raw), raw = raw,
             share = raw * measures$sign[i] / total * 100)
}))

plateau_raw <- cells$raw[cells$run %in% plateau_runs]
if (any(abs(plateau_raw) > plateau_limit + 1e-9)) {
  stop("plateau condition fails: a measure moved by more than ", plateau_limit,
       " in Runs ", min(plateau_runs), "-", max(plateau_runs))
}

cells$label <- ifelse(abs(cells$raw) < 1e-9, "0",
                ifelse(abs(cells$share) < 1, ifelse(cells$share > 0, "<1%", "−<1%"),
                       sub("^-", "−", sprintf("%d%%", as.integer(round(cells$share))))))
cells$fill <- pmax(pmin(cells$share, 100), -100)
cells$ink <- ifelse(abs(cells$fill) >= 55, "white", "#1f1f1f")
cells$measure <- factor(cells$measure, levels = rev(measures$label))

n_rows <- nrow(measures)
p <- ggplot(cells, aes(run, measure)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label, colour = ink), size = 3.1) +
  annotate("rect", xmin = min(plateau_runs) - 0.5, xmax = max(plateau_runs) + 0.5,
           ymin = 0.5, ymax = n_rows + 0.5, fill = NA, colour = "grey25", linewidth = 0.5) +
  annotate("text", x = mean(plateau_runs), y = n_rows + 0.75,
           label = paste0("plateau: no measure moved > ", plateau_limit, " points"),
           size = 3, colour = "grey25") +
  scale_fill_gradient2(low = "#B5541B", mid = "#f2f1ee", high = "#1F5F8B", midpoint = 0,
                       limits = c(-100, 100), breaks = c(-100, -50, 0, 50, 100),
                       labels = c("≤−100%", "−50%", "0", "50%", "≥100%"),
                       name = "Share of total\nchange") +
  scale_colour_identity() +
  scale_x_continuous(breaks = seq_len(max(cells$run)), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0, 0.6))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Fix run", y = NULL,
    title = "Where each measure's improvement came from",
    subtitle = "Share of the total Run 0 → 13 change contributed by each fix run",
    caption = paste0("Negative shares are losses: Run 1 traded pair consistency and excess pairs for recovery.\n",
                     "Excess pairs: a decrease counts as improvement. Shares beyond ±100% saturate the colour scale.")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, colour = "grey30", size = 8.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

ggsave(out_png, p, width = 9, height = 3.8, dpi = 300, bg = "white")
