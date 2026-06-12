library(ggplot2)
library(dplyr)

LF_TSV <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/Step8_Merging/lovd_flat_with_patients.tsv"

grp_cols <- c(Variant = "#2C6E91", Reference = "#9CB7C5", `Individual / Phenotype` = "#C0392B")

spec <- tibble::tribble(
  ~field,                          ~group,
  "VariantOnTranscript/DNA",       "Variant",
  "VariantOnGenome/DNA",           "Variant",
  "transcriptid",                  "Variant",
  "VariantOnTranscript/RNA",       "Variant",
  "VariantOnTranscript/Protein",   "Variant",
  "Reference/Location",            "Reference",
  "Reference/Authors",             "Reference",
  "Reference/Title",               "Reference",
  "Reference/PubMed_ID",           "Reference",
  "Individual/Gender",             "Individual / Phenotype",
  "Phenotype/Disease",             "Individual / Phenotype",
  "Individual/Origin/Geographic",  "Individual / Phenotype",
  "Phenotype/Additional",          "Individual / Phenotype",
  "Individual/Remarks",            "Individual / Phenotype",
  "Individual/Age",                "Individual / Phenotype"
)

d <- read.delim(LF_TSV, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
n <- nrow(d)

spec$pct <- sapply(spec$field, function(c) {
  v <- d[[c]]
  100 * sum(!is.na(v) & trimws(as.character(v)) != "") / n
})

spec <- spec %>%
  arrange(pct) %>%
  mutate(field = factor(field, levels = field),
         group = factor(group, levels = names(grp_cols)))

p <- ggplot(spec, aes(x = pct, y = field, fill = group)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            hjust = -0.2, size = 2.7, color = "grey20") +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 108), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Records with field populated (%)", y = NULL,
       title = "Completeness of migrated LOVD fields",
       caption = paste0("Field population across ", format(n, big.mark = ","),
                        " patient-variant records. Variant-level and reference fields are near-complete;\n",
                        "individual and phenotype annotation is sparser, reflecting heterogeneous curation in the source IDbases.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 11, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7, color = "grey40"),
        legend.position = "bottom",
        plot.margin = margin(5, 12, 8, 5, "mm"))

outdir <- "C:/Users/BornLoser/Desktop/Assignment/Thesis/04_Mutation_Processing/Output/figures"
ggsave(file.path(outdir, "fig07_field_completeness.pdf"), p, width = 6.5, height = 4.2, units = "in")
ggsave(file.path(outdir, "fig07_field_completeness.png"), p, width = 6.5, height = 4.2, units = "in", dpi = 300)
