# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
library(ape)
library(phangorn)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)

# ══════════════════════════════════════════════════════════════════════════════
# INPUT
input_dir <- "/home/aditri04/lu2025-17-24/BIMB60/TCGA/methylation_betas"

read_beta <- function(file, sample_name) {
  read_tsv(file,
           col_names = c("cpg_id", "beta", "chr", "start", "end"),
           col_types = "cncii") %>%
    mutate(sample = sample_name)
}

# ══════════════════════════════════════════════════════════════════════════════
# PROCESS SINGLE SAMPLE (UNCHANGED CORE LOGIC)
process_sample <- function(sample, k = 4) {

  tumor_file  <- file.path(input_dir, paste0(sample, "_T.tsv"))
  normal_file <- file.path(input_dir, paste0(sample, "_N.tsv"))

  if (!file.exists(tumor_file) || !file.exists(normal_file)) return(NULL)

  tumor  <- read_beta(tumor_file, "tumor")
  normal <- read_beta(normal_file, "normal")

  betas <- inner_join(normal, tumor,
                      by = c("cpg_id","chr","start","end"),
                      suffix = c("_normal","_tumor")) %>%
    filter(!is.na(beta_normal), !is.na(beta_tumor))

  # STEP 2 (unchanged)
  window_size <- 5e6

  windows <- betas %>%
    mutate(window = paste0(chr, ":", floor(start/window_size)*window_size)) %>%
    group_by(window, chr) %>%
    summarise(
      normal_beta = mean(beta_normal),
      tumor_beta  = mean(beta_tumor),
      n_cpg = n(),
      .groups="drop"
    ) %>%
    filter(n_cpg >= 20) %>%
    mutate(delta_beta = tumor_beta - normal_beta)

  # STEP 3–4 (unchanged)
  mat <- windows %>%
    select(window, normal_beta, tumor_beta) %>%
    column_to_rownames("window") %>%
    as.matrix()

  tree <- upgma(dist(mat))
  hc   <- as.hclust(tree)

  # STEP 5 (unchanged)
  cluster_ids <- cutree(hc, k = k)
  LETTERS_map <- setNames(LETTERS[1:k], 1:k)

  windows <- windows %>%
    left_join(tibble(window = names(cluster_ids),
                     subclone = LETTERS_map[as.character(cluster_ids)]),
              by="window")

  # STEP 6 (unchanged)
  trunk_cluster <- windows %>%
    group_by(subclone) %>%
    summarise(var_db = var(delta_beta), .groups="drop") %>%
    slice_min(var_db, n=1) %>%
    pull(subclone)

  windows <- windows %>%
    mutate(clonal_type = if_else(subclone == trunk_cluster, "trunk", "private"))

  # STEP 7 (MODIFIED: ALL EVENTS, NOT TOP 3)
  db_sd <- sd(windows$delta_beta)

  windows <- windows %>%
    mutate(
      state = case_when(
        delta_beta >  db_sd  ~ "hyper",
        delta_beta < -db_sd  ~ "hypo",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(state))

  # add cancer type
  cancer <- str_extract(sample, "^[A-Z]+")

  windows %>%
    mutate(sample = sample,
           cancer = cancer)
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN ALL SAMPLES
tumor_files <- list.files(input_dir, pattern="_T\\.tsv$", full.names=TRUE)
sample_ids  <- sub("_T\\.tsv$", "", basename(tumor_files))

all_data <- map_dfr(sample_ids, process_sample)

# ══════════════════════════════════════════════════════════════════════════════
# AGGREGATE TO CANCER LEVEL (KEY STEP)

plot_data <- all_data %>%
  group_by(cancer, sample, chr, clonal_type, state) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(cancer, chr, clonal_type, state) %>%
  summarise(mean_events = mean(n), .groups="drop")

plot_data <- plot_data %>%
  group_by(cancer) %>%
  mutate(
    total_events = sum(mean_events),
    proportion   = mean_events / total_events
  ) %>%
  ungroup()

# SIGN DIRECTION
# Absolute
plot_data$value_abs <- if_else(plot_data$state == "hyper",
                              plot_data$mean_events,
                              -plot_data$mean_events)

# Proportion
plot_data$value_prop <- if_else(plot_data$state == "hyper",
                               plot_data$proportion,
                               -plot_data$proportion)
# FIX CHROMOSOME ORDER (IMPORTANT)

plot_data$chr <- factor(
  plot_data$chr,
  levels = paste0("chr", c(1:22, "X", "Y"))
)

# CLEAN LABELS (REMOVE "chr")
plot_data$chr_label <- gsub("chr", "", plot_data$chr)

plot_data$chr_label <- factor(
  plot_data$chr_label,
  levels = c(as.character(1:22), "X", "Y")
)
# ══════════════════════════════════════════════════════════════════════════════
# FINAL PLOT

p_abs <- ggplot(plot_data,
                aes(x = chr, y = value_abs, fill = clonal_type)) +

  geom_col(width = 0.8) +

  facet_wrap(~ cancer, ncol = 4) +

  scale_fill_manual(values = c(
    trunk   = "#E41A1C",
    private = "#377EB8"
  )) +

  labs(
    title = "Chromosome-level burden of trunk and private methylation alterations across cancers",
    x = "Chromosome",
    y = "Average number of methylation events"
  ) +

  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold"),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  )
p_prop <- ggplot(plot_data,
                 aes(x = chr, y = value_prop, fill = clonal_type)) +

  geom_col(width = 0.8) +

  facet_wrap(~ cancer, ncol = 4) +

  scale_fill_manual(values = c(
    trunk   = "#E41A1C",
    private = "#377EB8"
  )) +

  labs(
    title = "Relative chromosomal enrichment of trunk and private methylation alterations across cancers",
    x = "Chromosome",
    y = "Proportion of methylation events"
  ) +

  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold"),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  )
# ══════════════════════════════════════════════════════════════════════════════
# SAVE
ggsave("methylation_events_absolute.pdf", p_abs, width = 18, height = 10)
ggsave("methylation_events_absolute.png", p_abs, width = 18, height = 10, dpi = 300)

ggsave("methylation_events_proportion.pdf", p_prop, width = 18, height = 10)
ggsave("methylation_events_proportion.png", p_prop, width = 18, height = 10, dpi = 300)
cat("DONE — both plots created and saved\n")
