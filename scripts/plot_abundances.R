# Plot relative abundances from a simulation manifest
# Usage:
#   Rscript scripts/plot_abundances.R manifest.tsv output.pdf

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_abundances.R manifest.tsv output.pdf", call. = FALSE)
}

manifest_path <- args[[1]]
out_path <- args[[2]]

library(ggplot2)
library(readr)
library(gridExtra)
d <- read_tsv(manifest_path, show_col_types = FALSE)

if (!all(c("species", "relative_abundance", "sample") %in% names(d))) {
  stop("Manifest must contain columns: sample, species, relative_abundance", call. = FALSE)
}

# Aggregate per species across genomes within each sample
agg <- aggregate(relative_abundance ~ sample + species, d, sum)

plots <- list()
for (s in unique(agg$sample)) {
  sub <- subset(agg, sample == s)
  sub <- sub[order(sub$relative_abundance, decreasing = TRUE), ]
  sub$species <- factor(sub$species, levels = sub$species)
  p <- ggplot(sub, aes(x = species, y = relative_abundance)) +
    geom_bar(stat = "identity", fill = "#4C78A8") +
    labs(title = paste("Relative abundance -", s),
         x = "Species",
         y = "Relative abundance") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(face = "bold"))
  plots[[s]] <- p
}

if (length(plots) == 1) {
  ggsave(out_path, plots[[1]], width = 10, height = 6)
} else {
  # If multiple samples, save one page per sample
  ggsave(out_path, width = 10, height = 6, plot = marrangeGrob(grobs = plots, nrow = 1, ncol = 1))
}

message("Wrote plot to ", out_path)
