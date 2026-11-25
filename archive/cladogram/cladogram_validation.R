# ════════════════════════════════════════════════════════════════
# LIBRARIES (unchanged)
library(ape)
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
# PATHS (UNCHANGED)
input_dir <- "/home/aditri04/lu2025-17-24/BIMB60/TCGA/methylation_betas"

# ════════════════════════════════════════════════════════════════
# READ FUNCTION (UNCHANGED)
read_beta <- function(file, sample_name) {
  read_tsv(file, col_names = c("cpg_id", "beta", "chr", "start", "end"),
           col_types = "cncii") %>%
    mutate(sample = sample_name)
}

# ════════════════════════════════════════════════════════════════
# MAIN PIPELINE FUNCTION (UNCHANGED)
run_cladogram <- function(sample, k = 4) {

  tumor_file  <- file.path(input_dir, paste0(sample, "_T.tsv"))
  normal_file <- file.path(input_dir, paste0(sample, "_N.tsv"))

  if (!file.exists(normal_file) || !file.exists(tumor_file)) return(NULL)

  tumor  <- read_beta(tumor_file, "tumor")
  normal <- read_beta(normal_file, "normal")

  betas <- inner_join(normal, tumor,
                      by = c("cpg_id", "chr", "start", "end"),
                      suffix = c("_normal", "_tumor")) %>%
    filter(!is.na(beta_normal), !is.na(beta_tumor))

  window_size <- 5e6

  windows <- betas %>%
    mutate(window = paste0(chr, ":", floor(start / window_size) * window_size)) %>%
    group_by(window, chr) %>%
    summarise(
      normal_beta = mean(beta_normal),
      tumor_beta  = mean(beta_tumor),
      n_cpg       = n(),
      .groups     = "drop"
    ) %>%
    filter(n_cpg >= 20) %>%
    mutate(delta_beta = tumor_beta - normal_beta)

  mat <- windows %>%
    select(window, normal_beta, tumor_beta) %>%
    column_to_rownames("window") %>%
    as.matrix()

  hc <- as.hclust(upgma(dist(mat)))

  cluster_ids <- cutree(hc, k = k)
  windows$subclone <- LETTERS[cluster_ids]

  trunk_cluster <- windows %>%
    group_by(subclone) %>%
    summarise(var_db = var(delta_beta), .groups = "drop") %>%
    slice_min(var_db) %>%
    pull(subclone)

  windows <- windows %>%
    mutate(clonal_type = if_else(subclone == trunk_cluster, "trunk", "private"))

  db_sd <- sd(windows$delta_beta)

  windows <- windows %>%
    mutate(
      state = case_when(
        delta_beta >  db_sd ~ "hyper",
        delta_beta < -db_sd ~ "hypo",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(state)) %>%
    mutate(sample = sample)

  return(windows)
}

# ════════════════════════════════════════════════════════════════
# RUN PIPELINE (UNCHANGED)
tumor_files <- list.files(input_dir, pattern = "_T\\.tsv$")
sample_ids <- sub("_T\\.tsv$", "", tumor_files)

all_data <- purrr::map_dfr(sample_ids, run_cladogram)

# ════════════════════════════════════════════════════════════════
# 🔬 VALIDATION (ADDED — DOES NOT MODIFY MAIN ANALYSIS)

set.seed(42)
n_iter <- 20

validation_results <- purrr::map_dfr(1:n_iter, function(i) {

  sampled_ids <- sample(sample_ids, size = round(length(sample_ids) * 0.7))

  resampled_data <- purrr::map_dfr(sampled_ids, run_cladogram)

  resampled_data$cancer <- stringr::str_extract(resampled_data$sample, "^[A-Z]+")

  tmp <- resampled_data %>%
    group_by(cancer, sample, chr, clonal_type, state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cancer, chr, clonal_type, state) %>%
    summarise(mean_events = mean(n), .groups = "drop") %>%
    group_by(cancer) %>%
    mutate(
      total_events = sum(mean_events),
      proportion   = mean_events / total_events
    ) %>%
    ungroup() %>%
    mutate(iteration = i)

  return(tmp)
})

# ════════════════════════════════════════════════════════════════
# MAIN ANALYSIS (UNCHANGED)

all_data$cancer <- stringr::str_extract(all_data$sample, "^[A-Z]+")

plot_data <- all_data %>%
  group_by(cancer, sample, chr, clonal_type, state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cancer, chr, clonal_type, state) %>%
  summarise(mean_events = mean(n), .groups = "drop")

plot_data <- plot_data %>%
  group_by(cancer) %>%
  mutate(
    total_events = sum(mean_events),
    proportion   = mean_events / total_events
  ) %>%
  ungroup()

plot_data$value_abs <- if_else(plot_data$state == "hyper",
                              plot_data$mean_events,
                              -plot_data$mean_events)

plot_data$value_prop <- if_else(plot_data$state == "hyper",
                               plot_data$proportion,
                               -plot_data$proportion)

plot_data$chr <- factor(plot_data$chr,
                        levels = paste0("chr", c(1:22, "X", "Y")))

plot_data$chr_label <- gsub("chr", "", plot_data$chr)

plot_data$chr_label <- factor(
  plot_data$chr_label,
  levels = c(as.character(1:22), "X", "Y")
)

cols <- c(trunk = "#F8766D", private = "#619CFF")

# ABSOLUTE PLOT (unchanged)
p_abs <- ggplot(plot_data,
                aes(x = chr_label, y = value_abs, fill = clonal_type)) +
  geom_col(width = 0.85) +
  facet_wrap(~ cancer, ncol = 4) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "Chromosomal distribution of trunk and private methylation events",
    x = "Chromosome",
    y = "Average number of methylation events (hyper = increased, hypo = decreased)"
  ) +
  theme_bw(base_size = 11)

# PROPORTION PLOT (unchanged)
p_prop <- ggplot(plot_data,
                 aes(x = chr_label, y = value_prop, fill = clonal_type)) +
  geom_col(width = 0.85) +
  facet_wrap(~ cancer, ncol = 4) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "Relative chromosomal distribution of trunk and private methylation events",
    x = "Chromosome",
    y = "Proportion of methylation events (hyper = increased, hypo = decreased)"
  ) +
  theme_bw(base_size = 11)

# SAVE (unchanged)
ggsave("methylation_absolute.png", p_abs, width = 18, height = 10, dpi = 300)
ggsave("methylation_absolute.pdf", p_abs, width = 18, height = 10)

ggsave("methylation_proportion.png", p_prop, width = 18, height = 10, dpi = 300)
ggsave("methylation_proportion.pdf", p_prop, width = 18, height = 10)

# ════════════════════════════════════════════════════════════════
# 🔬 VALIDATION RESULT (ADDED)

validation_summary <- validation_results %>%
  group_by(cancer, chr, clonal_type, state) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

comparison <- plot_data %>%
  select(cancer, chr, clonal_type, state, proportion) %>%
  left_join(validation_summary,
            by = c("cancer", "chr", "clonal_type", "state"))

correlation <- cor(comparison$proportion, comparison$mean_prop, use = "complete.obs")

print(paste("VALIDATION CORRELATION:", round(correlation, 3)))

cat("DONE — FINAL PIPELINE + VALIDATION COMPLETE\n")

