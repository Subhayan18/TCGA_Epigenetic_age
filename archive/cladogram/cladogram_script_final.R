# ════════════════════════════════════════════════════════════════
# LIBRARIES
# ════════════════════════════════════════════════════════════════
library(ape)
library(phangorn)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)

# ════════════════════════════════════════════════════════════════
# INPUT DIRECTORY (UNCHANGED)
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
# MAIN FUNCTION (PER SAMPLE)
# ════════════════════════════════════════════════════════════════
process_sample <- function(sample_id) {

  cat("Processing:", sample_id, "\n")
  flush.console()

  tumor_file  <- file.path(input_dir, paste0(sample_id, "_T.tsv"))
  normal_file <- file.path(input_dir, paste0(sample_id, "_N.tsv"))

  if (!file.exists(tumor_file) || !file.exists(normal_file)) return(NULL)

  tumor  <- read_beta(tumor_file, "tumor")
  normal <- read_beta(normal_file, "normal")

  betas <- inner_join(normal, tumor,
                      by = c("cpg_id", "chr", "start", "end"),
                      suffix = c("_normal", "_tumor")) %>%
    filter(!is.na(beta_normal), !is.na(beta_tumor))

  # ─────────────────────────────────────────────
  # WINDOWING
  # ─────────────────────────────────────────────
  windows <- betas %>%
    mutate(window = paste0(chr, ":", floor(start / 5e6) * 5e6)) %>%
    group_by(window, chr) %>%
    summarise(
      normal_beta = mean(beta_normal),
      tumor_beta  = mean(beta_tumor),
      n_cpg = n(),
      .groups = "drop"
    ) %>%
    filter(n_cpg >= 20) %>%
    mutate(delta_beta = tumor_beta - normal_beta)

  if (nrow(windows) == 0) return(NULL)

  # ─────────────────────────────────────────────
  # TREE (UNCHANGED CORE)
  # ─────────────────────────────────────────────
  mat <- windows %>%
    select(window, normal_beta, tumor_beta) %>%
    column_to_rownames("window") %>%
    as.matrix()

  hc <- as.hclust(upgma(dist(mat)))

  k <- 4
  cluster_ids <- cutree(hc, k = k)
  windows$subclone <- LETTERS[cluster_ids]

  # ─────────────────────────────────────────────
  # TRUNK (UNCHANGED)
  # ─────────────────────────────────────────────
  trunk_cluster <- windows %>%
    group_by(subclone) %>%
    summarise(var_db = var(delta_beta), .groups = "drop") %>%
    slice_min(var_db) %>%
    pull(subclone)

  windows <- windows %>%
    mutate(clonal_type = if_else(subclone == trunk_cluster,
                                 "Trunk", "Branch"))

  # ─────────────────────────────────────────────
  # EVENT DEFINITION (ALL EVENTS, NO FILTERING)
  # ─────────────────────────────────────────────
  events <- windows %>%
    mutate(
      state = case_when(
        delta_beta > 0  ~ "Hypermethylation",
        delta_beta < 0  ~ "Hypomethylation",
        TRUE            ~ NA_character_
      )
    ) %>%
    filter(!is.na(state))

  # ─────────────────────────────────────────────
  # CHROMOSOME SUMMARY (PER SAMPLE)
  # ─────────────────────────────────────────────
  chr_summary <- events %>%
    group_by(chr, state, clonal_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(sample = sample_id)

  return(chr_summary)
}

# ════════════════════════════════════════════════════════════════
# RUN ALL SAMPLES
# ════════════════════════════════════════════════════════════════
tumor_files <- list.files(input_dir, pattern = "_T\\.tsv$")
sample_ids  <- sub("_T\\.tsv$", "", tumor_files)

all_data <- map_dfr(sample_ids, process_sample)

# ════════════════════════════════════════════════════════════════
# ADD CANCER TYPE
# ════════════════════════════════════════════════════════════════
all_data$cancer <- str_extract(all_data$sample, "^[A-Z]+")

# FIX chromosome order
all_data$chr <- factor(
  all_data$chr,
  levels = paste0("chr", c(1:22, "X", "Y"))
)

all_data$chr_label <- gsub("chr", "", all_data$chr)

# ════════════════════════════════════════════════════════════════
# AGGREGATE TO CANCER LEVEL
# ════════════════════════════════════════════════════════════════
cancer_data <- all_data %>%
  group_by(cancer, chr, chr_label, state) %>%
  summarise(n = sum(n), .groups = "drop")

# ════════════════════════════════════════════════════════════════
# PLOT 1 — COUNTS
# ════════════════════════════════════════════════════════════════
p_counts <- ggplot(cancer_data,
                   aes(x = chr_label, y = n, fill = state)) +

  geom_col(width = 0.75) +

  facet_wrap(~ cancer, ncol = 4, scales = "free_y") +

  scale_fill_manual(values = c(
    "Hypermethylation" = "#D7302  "Hypomethylation"  = "#4575B4"
  ) +

  lbs(
    title = "Chromosome-level distribution of methylation alterations across cancers",
    x = "Chromosome",
    y = "Number of methylation events"
  ) 

  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text  = element_text(face = "bold")
    panel.grid  = element_blank(),
    legend.position = "top"
  )

# ═══════════════════════════════════════════════════════════════
# PLOT 2 — PROPORTIONS
# ════════════════════════════════════════════════════════════════
p_prop <- ggplot(cancer_data,
                 aes(x = chr_label, y = n, fill = state)) +
  geom_col(position = "fill", width = 0.75) +

  facet_wrap(~ cancer, ncol = 4) +

  scale_fill_manual(values = c(
    "Hypermethylation" = "#D73027",
    "Hypomethylation"  = "#4575B4"
  )) +

  labs
    title = "Chromosome-level distribution of methylation alterations across cancers",
    x = "Chromosome",
    y = "Proportion of methylation events"
  ) +

  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text  = element_text(face = "bold"),
    panel.grid  = element_blank(),
    legend.position = "top"
  )

# ════════════════════════════════════════════════════════
