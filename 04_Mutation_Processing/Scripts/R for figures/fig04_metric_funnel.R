library(ggplot2)

DEDUP_TSV <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step8_Merging/dedup_merged_variants.tsv"
DISP_TSV  <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step8_Merging/unresolved_disposition.tsv"

SOURCE_RECORDS    <- 5030
TOTAL_INDIVIDUALS <- 6259

d <- read.delim(DEDUP_TSV, stringsAsFactors = FALSE)
u <- read.delim(DISP_TSV, stringsAsFactors = FALSE)

resolved    <- nrow(d)
unresolved  <- nrow(u)
pairs       <- sum(d$patient_count, na.rm = TRUE)
individuals <- length(unique(trimws(unlist(strsplit(paste(d$accession_list, collapse = ";"), ";")))))
distinct    <- resolved + unresolved
pct_res     <- resolved / distinct * 100
rescuable   <- sum(grepl("^RESCUABLE", u$disposition))
pct_ceiling <- (resolved + rescuable) / distinct * 100

fmt <- function(x) format(x, big.mark = ",")

df <- data.frame(
  group = c("main", "resolved", "unresolved"),
  n = c(SOURCE_RECORDS, resolved, unresolved)
)
order_lab <- c("Source records\n(IDbases extraction)",
               "Distinct variants\nresolved / unresolved")
df$y <- factor(c(2, 1, 1), levels = c(1, 2), labels = rev(order_lab))
df$group <- factor(df$group, levels = c("main", "resolved", "unresolved"))

fills <- c(main = "#5A93AC", resolved = "#2C6E91", unresolved = "#C0392B")

lab_df <- data.frame(
  y = factor(c(2, 1), levels = c(1, 2), labels = rev(order_lab)),
  n = c(SOURCE_RECORDS, distinct),
  txt = c(fmt(SOURCE_RECORDS),
          paste0(fmt(resolved), " resolved  /  ", fmt(unresolved), " unresolved"))
)

xmax <- SOURCE_RECORDS * 1.3

p <- ggplot(df, aes(x = n, y = y, fill = group)) +
  geom_col(width = 0.55) +
  geom_text(data = lab_df, aes(x = n, y = y, label = txt),
            hjust = -0.05, size = 3, fontface = "bold",
            color = "grey20", inherit.aes = FALSE) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0))) +
  labs(x = "Count", y = NULL,
       title = "Pipeline data reduction and resolution",
       caption = paste0(fmt(SOURCE_RECORDS), " patient-variant records collapse to ",
                        fmt(distinct), " distinct variants; ",
                        sprintf("%.1f%%", pct_res), " resolved to GRCh38 HGVS ",
                        "(addressable ceiling ", sprintf("%.1f%%", pct_ceiling),
                        " incl. rescuable failures).\n",
                        "Resolved variants cover ", fmt(pairs),
                        " patient-variant pairs across ", fmt(individuals),
                        " individuals (of ", fmt(TOTAL_INDIVIDUALS), " in IDbases).")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9, color = "grey20", lineheight = 0.95),
        plot.title = element_text(size = 11, face = "bold"),
        plot.caption = element_text(hjust = 0.5, size = 7, color = "grey40"),
        plot.margin = margin(5, 10, 8, 5, "mm"))

outdir <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures"
ggsave(file.path(outdir, "fig04_metric_funnel.pdf"), p, width = 6.5, height = 3, units = "in")
ggsave(file.path(outdir, "fig04_metric_funnel.png"), p, width = 6.5, height = 3, units = "in", dpi = 300)
