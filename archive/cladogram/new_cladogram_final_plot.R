# ════════════════════════════════════════════════════════════════
# LIBRARIES
# ═══════════════════════════════════════════════════════════════library(ape)
library(phangorn)
library(ggdendro)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(grid)
library(ggplot2)
library(readr)
library(dplyr)
library(magrittr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)

# ════════════════════════════════════════════════════════════════
# INPUT DIRECTORY
# ════════════════════════════════════════════════════════════════
input_dir <- "/home/aditri04/lu2025-17-24/BIMB60/TCGA/methylation_betas"

# ════════════════════════════════════════════════════════════════
# READ FUNCTION
# ════════════════════════════════════════════════════════════════
read_beta <- function(file, sample_name) {
  read_tsv(file,
           col_names = c("cpg_id", "beta", "chr", "start", "end"),
           col_types = "cncii") %>%
    mutate(sample = sample_name)
}

# ════════════════════════════════════════════════════════════════
# SAFE WILCOXON
# ════════════════════════════════════════════════════════════════
safe_wilcox <- function(x) {
  if (length(unique(x)) < 2) return(1)
  wilcox.test(x, mu = 0, exact = FALSE)$p.value
}

# ════════════════════════════════════════════════════════════════
# MAIN FUNCTION (PHYLOGENETIC EVENTS)
# ════════════════════════════════════════════════════════════════
run_phylo_events <- function(sample, k = 4) {

  cat("Processing:", sample, "\n")
  flush.console()

  tumor_file  <- file.path(input_dir, paste0(sample, "_T.tsv"))
  normal_file <- file.path(input_dir, paste0(sample, "_N.tsv"))

  if (!file.exists(normal_file) || !file.exists(tumor_file)) return(NULL)

  tumor  <- read_beta(tumor_file, "tumor")
  normal <- read_beta(normal_file, "normal")

  betas <- inner_join(normal, tumor,
                      by = c("cpg_id", "chr", "start", "end"),
                      suffix = c("_normal", "_tumor")) %>%
    filter(!is.na(beta_normal), !is.na(beta_tumor))

  # 5Mb windows
  windows <- betas %>%
    mutate(window = paste0(chr, ":", floor(start / 5e6) * 5e6)) %>%
    group_by(window, chr) %>%
    summarise(
      normal_beta = mean(beta_normal),
      tumor_beta  = mean(beta_tumor),
      .groups = "drop"
    ) %>%
    mutate(delta_beta = tumor_beta - normal_beta)

  if (nrow(windows) < 10) return(NULL)

  # CLUSTERING (phylogeny)
  mat <- windows %>%
    select(window, normal_beta, tumor_beta) %>%
    column_to_rownames("window") %>%
    as.matrix()

  hc <- as.hclust(upgma(dist(mat)))
  cluster_ids <- cutree(hc, k = k)
  windows$subclone <- LETTERS[cluster_ids]

  # THRESHOLDS
  db_sd <- sd(windows$delta_beta)
  hi_cut <- db_sd
  lo_cut <- -db_sd

  # EXTRACT EVENTS
  events <- windows %>%
    group_by(subclone, chr) %>%
    summarise(mean_db = mean(delta_beta), .groups = "drop") %>%
    mutate(
      state = case_when(
        mean_db >  hi_cut ~ "hypermethylation",
        mean_db <  lo_cut ~ "hypomethylation",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(state)) %>%
    mutate(
      sample = sample,
      cancer = str_extract(sample, "^[A-Z]+")
    )

  return(events)
}

# ════════════════════════════════════════════════════════════════
# RUN ALL SAMPLES
# ════════════════════════════════════════════════════════════════
tumor_files <- list.files(input_dir, pattern = "_T\\.tsv$")
sample_ids <- sub("_T\\.tsv$", "", tumor_files)

all_events <- map_dfr(sample_ids, run_phylo_events)

# ════════════════════════════════════════════════════════════════
# SAMPLE-LEVEL COUNTS (FOR STATISTICS)
# ════════════════════════════════════════════════════════════════
sample_counts <- all_events %>%
  group_by(cancer, sample, chr, state) %>%
  summarise(n = n(), .groups = "drop")

# ════════════════════════════════════════════════════════════════
# AGGREGATE FOR PLOTTING
# ════════════════════════════════════════════════════════════════
plot_data <- sample_counts %>%
  group_by(cancer, chr, state) %>%
  summarise(mean_count = mean(n), .groups = "drop")

# ════════════════════════════════════════════════════════════════
# SIGNIFICANCE TEST
# ════════════════════════════════════════════════════════════════
stats <- sample_counts %>%
  group_by(cancer, chr, state) %>%
  summarise(p = safe_wilcox(n), .groups = "drop") %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

plot_data <- left_join(plot_data, stats,
                       by = c("cancer", "chr", "state"))

# ════════════════════════════════════════════════════════════════
# FIX CHROMOSOME ORDER
# ════════════════════════════════════════════════════════════════
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
plot_data$chr <- factor(plot_data$chr, levels = chr_levels)

# ════════════════════════════════════════════════════════════════
# FINAL PLOT
# ════════════════════════════════════════════════════════════════
p <- ggplot(plot_data,
            aes(x = chr, y = mean_count, fill = state)) +
  
  geom_col(position = "stack", width = 0.8, color = "black") +
  
  geom_text(
    data = subset(plot_data, signif != ""),
    aes(label = signif),
    position = position_stack(vjust = 1.05),
    size = 3
  ) +
  
  facet_wrap(~ cancer, ncol = 5, scales = "free_y") +
  
  scale_x_discrete(labels = function(x) gsub("chr", "", x)) +
  
  scale_fill_manual(values = c(
    "hypermethylation" = "#D73027",
    "hypomethylation"  = "#4575B4"
  )) +
  
  labs(
    title = "Phylogenetic methylation event frequencies across chromosomes",
    x = "Chromosome",
    y = "Average number of events per sample",
    fill = NULL
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "top"
  )

print(p)

# SAVE (IMPORTANT: separate line)
ggsave("phylo_methylation_final.pdf", p, width = 14, height = 8)
ggsave("phylo_methylation_final.png", p, width = 14, height = 8, dpi = 300)
