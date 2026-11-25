# ════════════════════════════════════════════════════════════════
# LIBRARIES
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
library(openxlsx)

# ════════════════════════════════════════════════════════════════
input_dir <- "/home/aditri04/lu2025-17-24/BIMB60/TCGA/methylation_betas"
# ════════════════════════════════════════════════════════════════
# READ FUNCTION 
read_beta <- function(file, sample_name) {
  read_tsv(file, col_names = c("cpg_id", "beta", "chr", "start", "end"),
           col_types = "cncii") %>%
    mutate(sample = sample_name)
}

# ════════════════════════════════════════════════════════════════
# COMPLEXITY SCORE FUNCTION
#
# Implements the formula from Martin & Hirschberg adapted for phylogenetic trees:
#   Total complexity = { Σᵢ rᵢ × 2^(i−1) / (K × (K+1)) }^(1/(N−1))
#
#   rᵢ = number of leaf (end) nodes first appearing at generation i
#   K  = total number of edges in the tree
#   N  = total number of nodes in the tree
#
# Complexity is computed on the UPGMA hclust object produced per sample.

compute_tree_complexity <- function(tree_input) {
  # Accepts either a hclust or phylo object directly
  tree     <- if (inherits(tree_input, "phylo")) tree_input else as.phylo(tree_input)
  n_leaves <- length(tree$tip.label)
  N         <- n_leaves + tree$Nnode   # total nodes (leaves + internal)
  K         <- nrow(tree$edge)         # total edges (= N-1 for a tree)
  root_node <- n_leaves + 1L

  # BFS to assign generation depth (0 = root) to every node
  node_gen            <- rep(NA_integer_, N)
  node_gen[root_node] <- 0L
  queue               <- root_node

  while (length(queue) > 0L) {
    cur   <- queue[1L]
    queue <- queue[-1L]
    kids  <- tree$edge[tree$edge[, 1L] == cur, 2L]
    for (kid in kids) {
      node_gen[kid] <- node_gen[cur] + 1L
      queue <- c(queue, kid)
    }
  }

  max_gen    <- max(node_gen, na.rm = TRUE)
  leaf_nodes <- seq_len(n_leaves)   # tip nodes are indexed 1..n_leaves in ape

  # r_i = count of leaf nodes at generation i  (i = 1 … max_gen)
  # Generation 0 is the root (never a leaf), contributing 0 by definition.
  r <- vapply(seq_len(max_gen), function(i) {
    sum(which(node_gen == i) %in% leaf_nodes)
  }, integer(1L))

  # Σᵢ rᵢ × 2^(i−1) for i = 1 … max_gen
  raw_score  <- sum(r * 2^(seq_len(max_gen) - 1L))
  complexity <- (raw_score / (K * (K + 1L)))^(1 / (N - 1L))

  tibble(
    n_nodes       = N,
    n_edges       = K,
    n_leaves      = n_leaves,
    n_generations = max_gen,
    raw_score     = raw_score,
    complexity    = complexity
  )
}

# ════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# Now returns a list: $windows (existing output) + $complexity (new)

run_cladogram <- function(sample, k = 4) {

  tumor_file  <- file.path(input_dir, paste0(sample, "_T.tsv"))
  normal_file <- file.path(input_dir, paste0(sample, "_N.tsv"))

  if (!file.exists(normal_file) || !file.exists(tumor_file)) return(NULL)

  tumor  <- read_beta(tumor_file,  "tumor")
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

  # TREE + CLUSTERING 
  mat <- windows %>%
    select(window, normal_beta, tumor_beta) %>%
    column_to_rownames("window") %>%
    as.matrix()

  hc <- as.hclust(upgma(dist(mat)))

  cluster_ids      <- cutree(hc, k = k)
  windows$subclone <- LETTERS[cluster_ids]

  # ── COMPLEXITY SCORE ──────────────────────────────────────────
  # Cut the original UPGMA tree at k=4 branches using keep.tip().
  # One representative window per subclone (closest to cluster centre)
  # is kept. The pruned phylo object is passed directly to
  # compute_tree_complexity — no as.hclust() conversion needed,
  # which previously broke ultrametricity after pruning.
  complexity_row <- tryCatch({
    rep_tips <- vapply(LETTERS[seq_len(k)], function(cl) {
      idx  <- which(windows$subclone == cl)
      if (length(idx) == 0L) stop("empty subclone: ", cl)
      if (length(idx) == 1L) return(idx)
      dmat <- as.matrix(dist(mat[idx, , drop = FALSE]))
      idx[which.min(rowMeans(dmat))]
    }, integer(1L))

    tip_names   <- rownames(mat)[rep_tips]
    pruned_tree <- keep.tip(as.phylo(hc), tip_names)

    # Pass pruned phylo directly — avoids unreliable as.hclust() on
    # a non-ultrametric pruned tree
    compute_tree_complexity(pruned_tree) %>% mutate(sample = sample, .before = 1)
  },
  error = function(e) {
    warning("Complexity computation failed for sample: ", sample, " — ", conditionMessage(e))
    tibble(sample = sample,
           n_nodes = NA_integer_, n_edges = NA_integer_,
           n_leaves = NA_integer_, n_generations = NA_integer_,
           raw_score = NA_real_, complexity = NA_real_)
  })
  # ─────────────────────────────────────────────────────────────

  trunk_cluster <- windows %>%
    group_by(subclone) %>%
    summarise(var_db = var(delta_beta), .groups = "drop") %>%
    slice_min(var_db) %>%
    pull(subclone)

  windows <- windows %>%
    mutate(clonal_type = if_else(subclone == trunk_cluster, "trunk", "private"))

  # METHYLATION STATE 
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

  list(windows = windows, complexity = complexity_row)
}

# ════════════════════════════════════════════════════════════════
# RUN PIPELINE

tumor_files <- list.files(input_dir, pattern = "_T\\.tsv$")
sample_ids  <- sub("_T\\.tsv$", "", tumor_files)

# Run pipeline — returns a list of lists; NULL entries (missing files) are dropped
all_results <- purrr::map(sample_ids, run_cladogram) %>%
  purrr::compact()

# Extract windows and complexity tables separately
all_data        <- purrr::map_dfr(all_results, "windows")
complexity_data <- purrr::map_dfr(all_results, "complexity")

# ════════════════════════════════════════════════════════════════
# ADD CANCER TYPE

all_data$cancer        <- stringr::str_extract(all_data$sample,        "^[A-Z]+")
complexity_data$cancer <- stringr::str_extract(complexity_data$sample, "^[A-Z]+")

# ════════════════════════════════════════════════════════════════
# AGGREGATION 

plot_data <- all_data %>%
  group_by(cancer, sample, chr, clonal_type, state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cancer, chr, clonal_type, state) %>%
  summarise(mean_events = mean(n), .groups = "drop")

# ════════════════════════════════════════════════════════════════
# PROPORTIONS 

plot_data <- plot_data %>%
  group_by(cancer) %>%
  mutate(
    total_events = sum(mean_events),
    proportion   = mean_events / total_events
  ) %>%
  ungroup()

# ════════════════════════════════════════════════════════════════
# SIGNED VALUES (HYPER UP, HYPO DOWN)

plot_data$value_abs <- if_else(plot_data$state == "hyper",
                               plot_data$mean_events,
                               -plot_data$mean_events)

plot_data$value_prop <- if_else(plot_data$state == "hyper",
                                plot_data$proportion,
                                -plot_data$proportion)

# ════════════════════════════════════════════════════════════════
# CLEAN CHROMOSOME LABELS 

plot_data$chr <- factor(plot_data$chr,
                        levels = paste0("chr", c(1:22, "X", "Y")))

plot_data$chr_label <- gsub("chr", "", plot_data$chr)

plot_data$chr_label <- factor(
  plot_data$chr_label,
  levels = c(as.character(1:22), "X", "Y")
)

# ════════════════════════════════════════════════════════════════
# COLORS 

cols <- c(
  trunk   = "#F8766D",
  private = "#619CFF"
)

# ════════════════════════════════════════════════════════════════
# ABSOLUTE PLOT

p_abs <- ggplot(plot_data,
                aes(x = chr_label, y = value_abs, fill = clonal_type)) +
  geom_col(width = 0.85) +
  facet_wrap(~ cancer, ncol = 4) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "Chromosomal distribution of trunk and private methylation events",
    x = "Chromosome",
    y = "Average number of methylation events (↑ hyper, ↓ hypo)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold", size = 10),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ════════════════════════════════════════════════════════════════
# PROPORTION PLOT

p_prop <- ggplot(plot_data,
                 aes(x = chr_label, y = value_prop, fill = clonal_type)) +
  geom_col(width = 0.85) +
  facet_wrap(~ cancer, ncol = 4) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "Relative chromosomal distribution of trunk and private methylation events",
    x = "Chromosome",
    y = "Proportion of methylation events (↑ hyper, ↓ hypo)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = element_text(face = "bold", size = 10),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ════════════════════════════════════════════════════════════════
# SAVE (PNG + PDF)

ggsave("methylation_absolute.png", p_abs, width = 18, height = 10, dpi = 300)
ggsave("methylation_absolute.pdf", p_abs, width = 18, height = 10)

ggsave("methylation_proportion.png", p_prop, width = 18, height = 10, dpi = 300)
ggsave("methylation_proportion.pdf", p_prop, width = 18, height = 10)

# ════════════════════════════════════════════════════════════════
# COMPLEXITY BOXPLOT
#
# One-sample Wilcoxon test per cancer type vs. the overall median complexity,
# with Benjamini-Hochberg correction across cancer types.
# This avoids false-positive inflation from pairwise comparisons while still
# giving a meaningful significance label (above/below average complexity).

global_median_complexity <- median(complexity_data$complexity, na.rm = TRUE)

complexity_stats <- complexity_data %>%
  group_by(cancer) %>%
  summarise(
    median_c = median(complexity, na.rm = TRUE),
    p_raw    = tryCatch(
      wilcox.test(complexity, mu = global_median_complexity, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj  = p.adjust(p_raw, method = "BH"),
    sig    = case_when(
      is.na(p_adj)  ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

# Sort cancer types by median complexity (descending = highest at top)
cancer_order <- complexity_stats %>%
  arrange(median_c) %>%
  pull(cancer)

complexity_data$cancer  <- factor(complexity_data$cancer,  levels = cancer_order)
complexity_stats$cancer <- factor(complexity_stats$cancer, levels = cancer_order)

# All significance labels at a single fixed x position (just past the global max)
x_label_pos <- max(complexity_data$complexity, na.rm = TRUE) +
               0.02 * diff(range(complexity_data$complexity, na.rm = TRUE))

p_complexity <- ggplot(complexity_data,
                       aes(x = complexity, y = cancer)) +
  geom_vline(xintercept = global_median_complexity,
             linetype = "dashed", colour = "#CC0000", linewidth = 0.6) +
  geom_boxplot(fill      = "#6BAED6",
               colour    = "black",
               linewidth = 0.45,
               outlier.shape = 16,
               outlier.size  = 1.2) +
  geom_text(data = complexity_stats,
            aes(y = cancer, label = sig),
            x = x_label_pos, hjust = 0, size = 3.5, inherit.aes = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title   = "Phylogenetic Tree Complexity by Cancer Type",
    x       = "Complexity Score",
    y       = "Cancer Type",
    caption = "Red dashed line = cross-cancer median complexity\nSignificance vs. median: one-sample Wilcoxon, BH-adjusted"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.caption  = element_text(size = 8, colour = "grey40", hjust = 0),
    axis.text.y   = element_text(size = 10),
    axis.text.x   = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave("complexity_boxplot.png", p_complexity, width = 8, height = 6, dpi = 300)
ggsave("complexity_boxplot.pdf", p_complexity, width = 8, height = 6)

# ════════════════════════════════════════════════════════════════
# EXPORT MULTI-SHEET EXCEL (SUPERVISOR-FRIENDLY, CONSISTENT COLORS)

# ---- 1. CLEAN EVENT TABLE (LONG FORMAT)
events_sheet <- all_data %>%
  mutate(chr_label = gsub("chr", "", chr)) %>%
  select(
    cancer,
    sample,
    chr_label,
    window,
    delta_beta,
    subclone,
    clonal_type,
    state
  ) %>%
  arrange(cancer, sample, chr_label)

# ---- 2. CHROMOSOME SUMMARY
chr_summary <- events_sheet %>%
  group_by(cancer, chr_label, clonal_type, state) %>%
  summarise(n_events = n(), .groups = "drop")

# ---- 3. SUBCLONE SUMMARY
subclone_summary <- events_sheet %>%
  group_by(cancer, sample, subclone, clonal_type) %>%
  summarise(n_events = n(), .groups = "drop")

# ---- 4. SAMPLE SUMMARY
sample_summary <- events_sheet %>%
  group_by(cancer, sample, clonal_type, state) %>%
  summarise(n_events = n(), .groups = "drop")

# ════════════════════════════════════════════════════════════════
# COMPLEXITY SHEETS
#
# Sheet 5 — complexity_raw:   one row per sample, all tree metrics + complexity score
# Sheet 6 — complexity_cancer: one row per cancer type, summary statistics

# ---- 5. COMPLEXITY RAW (per sample)
complexity_raw <- complexity_data %>%
  select(cancer, sample, n_nodes, n_edges, n_leaves, n_generations, raw_score, complexity) %>%
  arrange(cancer, sample)

# ---- 6. COMPLEXITY SUMMARY (per cancer type)
complexity_cancer <- complexity_raw %>%
  group_by(cancer) %>%
  summarise(
    n_samples          = n(),
    mean_complexity    = mean(complexity,    na.rm = TRUE),
    sd_complexity      = sd(complexity,      na.rm = TRUE),
    median_complexity  = median(complexity,  na.rm = TRUE),
    min_complexity     = min(complexity,     na.rm = TRUE),
    max_complexity     = max(complexity,     na.rm = TRUE),
    mean_n_nodes       = mean(n_nodes,       na.rm = TRUE),
    mean_n_edges       = mean(n_edges,       na.rm = TRUE),
    mean_n_leaves      = mean(n_leaves,      na.rm = TRUE),
    mean_n_generations = mean(n_generations, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_complexity))

# ════════════════════════════════════════════════════════════════
# ---- CREATE WORKBOOK

wb <- createWorkbook()

# ── Style definitions ─────────────────────────────────────────
headerStyle    <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9",
                              border = "Bottom", borderColour = "#808080")
headerStylePurple <- createStyle(textDecoration = "bold", fgFill = "#E8D5F5",
                                 border = "Bottom", borderColour = "#6A0DAD")
numStyle2dp    <- createStyle(numFmt = "0.00")
numStyle4dp    <- createStyle(numFmt = "0.0000")

# ========================
# EVENTS SHEET
# ========================
addWorksheet(wb, "events")
writeData(wb, "events", events_sheet, withFilter = TRUE)

addStyle(wb, "events", headerStyle,
         rows = 1, cols = 1:ncol(events_sheet), gridExpand = TRUE)

setColWidths(wb, "events", cols = 1:ncol(events_sheet), widths = "auto")

# Color by clonal type (matches plots)
clonal_col <- which(names(events_sheet) == "clonal_type")

conditionalFormatting(
  wb, "events",
  cols = clonal_col,
  rows = 2:(nrow(events_sheet) + 1),
  rule = '=="trunk"',
  style = createStyle(bgFill = "#FEE0D2")
)

conditionalFormatting(
  wb, "events",
  cols = clonal_col,
  rows = 2:(nrow(events_sheet) + 1),
  rule = '=="private"',
  style = createStyle(bgFill = "#DEEBF7")
)

# ========================
# OTHER SUMMARY SHEETS
# ========================
addWorksheet(wb, "chromosome_summary")
writeData(wb, "chromosome_summary", chr_summary, withFilter = TRUE)
setColWidths(wb, "chromosome_summary", cols = 1:ncol(chr_summary), widths = "auto")

addWorksheet(wb, "subclone_summary")
writeData(wb, "subclone_summary", subclone_summary, withFilter = TRUE)
setColWidths(wb, "subclone_summary", cols = 1:ncol(subclone_summary), widths = "auto")

addWorksheet(wb, "sample_summary")
writeData(wb, "sample_summary", sample_summary, withFilter = TRUE)
setColWidths(wb, "sample_summary", cols = 1:ncol(sample_summary), widths = "auto")

# ========================
# COMPLEXITY RAW SHEET
# ========================
addWorksheet(wb, "complexity_raw")
writeData(wb, "complexity_raw", complexity_raw, withFilter = TRUE)

addStyle(wb, "complexity_raw", headerStylePurple,
         rows = 1, cols = 1:ncol(complexity_raw), gridExpand = TRUE)

setColWidths(wb, "complexity_raw", cols = 1:ncol(complexity_raw), widths = "auto")

# Format complexity column to 4 dp for readability
complexity_col_idx <- which(names(complexity_raw) == "complexity")
raw_score_col_idx  <- which(names(complexity_raw) == "raw_score")

addStyle(wb, "complexity_raw", numStyle4dp,
         rows = 2:(nrow(complexity_raw) + 1),
         cols = complexity_col_idx, gridExpand = TRUE)

addStyle(wb, "complexity_raw", numStyle2dp,
         rows = 2:(nrow(complexity_raw) + 1),
         cols = raw_score_col_idx, gridExpand = TRUE)

# Manual purple-intensity fill on complexity column (avoids colourScale corruption)
# Rescale complexity to 0-1 within the sheet, then map to white→purple hex fills
complexity_vals <- complexity_raw$complexity
c_min <- min(complexity_vals, na.rm = TRUE)
c_max <- max(complexity_vals, na.rm = TRUE)
c_range <- if (c_max == c_min) 1 else c_max - c_min

for (i in seq_along(complexity_vals)) {
  if (!is.na(complexity_vals[i])) {
    scale_val <- (complexity_vals[i] - c_min) / c_range  # 0 (low) → 1 (high)
    r_hex <- format(as.hexmode(round(255 - scale_val * (255 - 106))), width = 2, upper.case = TRUE)
    g_hex <- format(as.hexmode(round(255 - scale_val * (255 - 13))),  width = 2, upper.case = TRUE)
    b_hex <- format(as.hexmode(round(255 - scale_val * (255 - 173))), width = 2, upper.case = TRUE)
    fill_hex <- paste0(r_hex, g_hex, b_hex)
    addStyle(wb, "complexity_raw",
             style = createStyle(fgFill = fill_hex, numFmt = "0.0000"),
             rows = i + 1, cols = complexity_col_idx)
  }
}

# ========================
# COMPLEXITY CANCER SUMMARY SHEET
# ========================
addWorksheet(wb, "complexity_cancer")
writeData(wb, "complexity_cancer", complexity_cancer, withFilter = TRUE)

addStyle(wb, "complexity_cancer", headerStylePurple,
         rows = 1, cols = 1:ncol(complexity_cancer), gridExpand = TRUE)

setColWidths(wb, "complexity_cancer", cols = 1:ncol(complexity_cancer), widths = "auto")

# Format all numeric columns to 4 dp
num_cols_cancer <- which(sapply(complexity_cancer, is.numeric))
addStyle(wb, "complexity_cancer", numStyle4dp,
         rows = 2:(nrow(complexity_cancer) + 1),
         cols = num_cols_cancer, gridExpand = TRUE)

# Manual purple-intensity fill on mean_complexity column
mean_complexity_col <- which(names(complexity_cancer) == "mean_complexity")
mc_vals  <- complexity_cancer$mean_complexity
mc_min   <- min(mc_vals, na.rm = TRUE)
mc_max   <- max(mc_vals, na.rm = TRUE)
mc_range <- if (mc_max == mc_min) 1 else mc_max - mc_min

for (i in seq_along(mc_vals)) {
  if (!is.na(mc_vals[i])) {
    scale_val <- (mc_vals[i] - mc_min) / mc_range
    r_hex <- format(as.hexmode(round(255 - scale_val * (255 - 106))), width = 2, upper.case = TRUE)
    g_hex <- format(as.hexmode(round(255 - scale_val * (255 - 13))),  width = 2, upper.case = TRUE)
    b_hex <- format(as.hexmode(round(255 - scale_val * (255 - 173))), width = 2, upper.case = TRUE)
    fill_hex <- paste0(r_hex, g_hex, b_hex)
    addStyle(wb, "complexity_cancer",
             style = createStyle(fgFill = fill_hex, numFmt = "0.0000"),
             rows = i + 1, cols = mean_complexity_col)
  }
}

# ---- SAVE FILE
saveWorkbook(wb, "methylation_events.xlsx", overwrite = TRUE)

cat("Excel file saved: methylation_events.xlsx\n")
cat("  Sheets: events | chromosome_summary | subclone_summary | sample_summary |\n")
cat("          complexity_raw | complexity_cancer\n")
cat("DONE — FINAL PIPELINE + PLOTS + COMPLEXITY COMPLETE\n")
