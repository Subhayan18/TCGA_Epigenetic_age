# ══════════════════════════════════════════════════════════════════════════════
#loading required libraries
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
# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: READ METHYLATION FILE & LOAD DATA

base_dir <- "/lunarc/nobackup/projects/mehtyl_mrd_dat/Meth_TCGA_GSE213478/cladogram_results"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

read_beta <- function(file, sample_name) {
  read_tsv(file, col_names = c("cpg_id", "beta", "chr", "start", "end"),
           col_types = "cncii") %>%
    mutate(sample = sample_name)
}

run_cladogram <- function(sample, k = 4) {
 tryCatch({

  message("Processing sample: ", sample)

  tumor_file  <- paste0(sample, "_T.tsv")
  normal_file <- paste0(sample, "_N.tsv")

  if (!file.exists(normal_file) | !file.exists(tumor_file)) {
    message("Skipping (missing file): ", sample)
    return(NULL)
  }

    #CREATE SAMPLE SPECIFIC FOLDER
    sample_dir <- file.path(base_dir, sample)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

    #LOAD AND MERGE TUMOR AND NORMAL SAMPLES
    tumor  <- read_beta(tumor_file, "tumor")
    normal <- read_beta(normal_file, "normal")

    betas <- inner_join(normal, tumor,
                        by     = c("cpg_id", "chr", "start", "end"),
                        suffix = c("_normal", "_tumor")) %>%
    filter(!is.na(beta_normal), !is.na(beta_tumor))
    message("CpGs after NA filtering: ", nrow(betas))

    # ══════════════════════════════════════════════════════════════════════════════
    #STEP 2: CREATE CHROMOSOMAL WINDOWS (5Mb)
    window_size <- 5e6

    windows <- betas %>%
    mutate(window = paste0(chr, ":", floor(start / window_size) * window_size)) %>%
    group_by(window, chr) %>%
    summarise(
        normal_beta = mean(beta_normal, na.rm = TRUE),
        tumor_beta  = mean(beta_tumor,  na.rm = TRUE),
        n_cpg       = n(),
        .groups     = "drop"
    ) %>%
    filter(n_cpg >= 20) %>%
    mutate(delta_beta = tumor_beta - normal_beta)

    message("Genomic windows after filtering: ", nrow(windows))
    message("delta_beta range: ",
        round(min(windows$delta_beta), 3), " to ",
        round(max(windows$delta_beta), 3))

    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 3: BUILD DISTANCE MATRIX & TREE
    mat <- windows %>%
        select(window, normal_beta, tumor_beta) %>%
        column_to_rownames("window") %>%
        as.matrix()

    dist_mat <- dist(mat, method = "euclidean")

    tree <- upgma(dist_mat)
    hc   <- as.hclust(tree)
    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 4: ASSIGN SUBCLONES
    cluster_ids <- cutree(hc, k = k)
    LETTERS_map <- setNames(LETTERS[seq_len(k)], seq_len(k))

    windows <- windows %>%
        left_join(
        tibble(window   = names(cluster_ids),
                subclone = LETTERS_map[as.character(cluster_ids)]),
        by = "window"
    )

    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 5: IDENTIFY TRUNK CLUSTER
    trunk_cluster <- windows %>%
        group_by(subclone) %>%
        summarise(var_db = var(delta_beta), .groups = "drop") %>%
        slice_min(var_db, n = 1) %>%
        pull(subclone)

    windows <- windows %>%
        mutate(clonal_type = if_else(subclone == trunk_cluster, "trunk", "private"))

    message("Trunk cluster: ", trunk_cluster)
    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 6: DETECT METHYLATION EVENTS (TOP #3 PER CLUSTER)
    db_sd <- sd(windows$delta_beta) 
    hi_cut <-  db_sd
    lo_cut <- -db_sd

    message("Thresholds — hyper: >", round(hi_cut, 3),
            "  hypo: <", round(lo_cut, 3))

    meth_annotation <- windows %>%
        group_by(subclone, chr) %>%
        summarise(mean_db = mean(delta_beta), .groups = "drop") %>%
        mutate(
            state = case_when(
         mean_db >  hi_cut ~ "hyper",   # strong gain of methylation
         mean_db <  lo_cut ~ "hypo",    # strong loss of methylation
         TRUE               ~ NA_character_
            )
        ) %>%
        filter(!is.na(state)) %>% 
        mutate(strength = abs(mean_db)) %>% 
        arrange(subclone, desc(strength)) %>%  #arrange subclones by strongest signals
        group_by(subclone) %>% 
        slice_head(n = 3) %>% # prevents long and unreadable tables
        mutate(event_label = paste0(sub("chr", "", chr), " ", state)) %>%
        summarise(meth_label = paste(event_label, collapse = ", "), .groups = "drop")

    windows <- windows %>%
        left_join(meth_annotation, by = "subclone")
    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 7: BUILD CLUSTER SUMMMARY TABLE
    cluster_summary <- windows %>%
        group_by(subclone, clonal_type) %>%
        summarise(
            n_windows   = n(),
            mean_delta  = round(mean(delta_beta), 4),
            sd_delta    = round(sd(delta_beta),   4),
            meth_events = first(meth_label),
            .groups     = "drop"
        ) %>%
        arrange(subclone) %>%
        mutate(
            meth_events = replace_na(meth_events, "stable"),
            clonal_type = if_else(clonal_type == "trunk",
                                paste0(clonal_type, " *"), clonal_type)
        ) %>%
        rename(
            Subclone      = subclone,
            Type          = clonal_type,
            Windows       = n_windows,
            `Mean dB`     = mean_delta,
            `SD dB`       = sd_delta,
            `Meth events` = meth_events
        )

    print(cluster_summary)
    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 8: DENDROGRAM DATA
    ddata <- dendro_data(hc, type = "rectangle")

    tip_meta <- ddata$labels %>%
        rename(window = label) %>%
        left_join(windows, by = "window")

    seg_meta <- ddata$segments %>%
        left_join(tip_meta %>% select(x, subclone),
                by = c("xend" = "x")) %>%
    mutate(
        subclone = replace_na(subclone, "unassigned"),
        is_trunk = subclone == "unassigned"
    )
    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 9: COLORS AND GEOMETRY
    subclone_colors <- setNames(
        colorRampPalette(brewer.pal(8, "Set2"))(k),
        LETTERS[seq_len(k)]
    )

    cluster_ranges <- tip_meta %>%
        group_by(subclone) %>%
        summarise(
            x_min = min(x) - 0.4,
            x_max = max(x) + 0.4,
            .groups = "drop"
        )

    y_max <- max(ddata$segments$y)
    y_bar <- -y_max * 0.04
    y_lab <- -y_max * 0.10
    y_evt <- -y_max * 0.17
    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 10: PLOT
    message("Plotting sample: ", sample)
    p_tree <- ggplot() +

    # ---- TRUNK ----
        geom_segment(
            data = seg_meta %>% filter(is_trunk),
            aes(x = x, y = y, xend = xend, yend = yend),
            color = "#C0392B",
            linewidth = 2.2,
            lineend = "round"
        ) +

    # ---- SUBCLONE BRANCHES ----
        geom_segment(
            data = seg_meta %>% filter(!is_trunk),
            aes(x = x, y = y, xend = xend, yend = yend, color = subclone),
            linewidth = 0.8,
            lineend = "round"
        ) +

    # ---- TIPS ----
        geom_point(
            data = tip_meta,
            aes(x = x, y = y, fill = subclone),
            shape = 21,
            color = "white",
            size = 1.5
        ) +

    # ---- CLUSTER BARS ----
        geom_segment(
            data = cluster_ranges,
            aes(x = x_min, xend = x_max,
                y = y_bar, yend = y_bar,
                color = subclone),
            linewidth = 3,
            lineend = "butt"
        ) +

    # ---- SUBCLONE LABELS ----
        geom_text(
            data = cluster_ranges %>%
                left_join(cluster_summary %>% select(Subclone, `Meth events`),
                        by = c("subclone" = "Subclone")),
            aes(x = (x_min + x_max) / 2,
                y = y_lab,
                label = subclone,
                color = subclone),
            size = 5,
            fontface = "bold"
        ) +

    # ---- EVENT LABELS ----
        geom_text(
            data = cluster_ranges %>%
                left_join(cluster_summary %>% select(Subclone, `Meth events`),
                        by = c("subclone" = "Subclone")),
            aes(x = (x_min + x_max) / 2,
                y = y_evt,
                label = `Meth events`),
            size = 2.8,
            color = "grey40",
            fontface = "italic"
        ) +

    # ---- COLORS ----
        scale_color_manual(values = subclone_colors,
                            breaks = LETTERS[seq_len(k)]) +
        scale_fill_manual(values = subclone_colors,
                            breaks = LETTERS[seq_len(k)]) +

    # ---- AXES ----
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
        scale_y_continuous(
            expand = expansion(mult = c(0.05, 0.1)),
            limits = c(y_evt - y_max * 0.05, y_max * 1.05)
        ) +

        coord_flip() +

    # ---- THEME ----
        theme_classic(base_size = 12) +
        theme(
            plot.background  = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            axis.line        = element_blank(),
            axis.text        = element_blank(),
            axis.ticks       = element_blank(),
            axis.title       = element_blank(),
            panel.grid       = element_blank(),
            legend.position  = "none",
            plot.title       = element_text(size = 14, face = "bold", color = "grey20"),
            plot.subtitle    = element_text(size = 9,  color = "grey50"),
            plot.margin      = margin(20, 20, 50, 20)
        ) +

    # AUTOMATED TITLE
        labs(
            title = paste0("Epigenetic cladogram - ", sample),
            subtitle = paste0(
                k, " subclones  |  UPGMA  |  5 Mb windows  |  ",
                "red = shared trunk (common ancestor)"
            )
        )
    # ══════════════════════════════════════════════════════════════════════════════
    # STEP 11: BUILD TABLE GROB
    table_grob <- tableGrob(
        cluster_summary,
        rows  = NULL,
        theme = ttheme_minimal(
            core = list(
                fg_params = list(cex = 0.80),
                bg_params = list(
                    fill = rep(c("white", "#f0f0f0"), length.out = nrow(cluster_summary)),
                    col  = NA
                )
            ),
            colhead = list(
                fg_params = list(cex = 0.85, fontface = "bold"),
                bg_params = list(fill = "grey88", col = NA)
            ),
            padding = unit(c(3, 6), "mm")
        )
    )
    p_final <- ggdraw() +
        draw_plot(p_tree,
                        x = 0, y = 0, width = 1, height = 1) +
        draw_grob(table_grob,
                        x      = 0.70,   # left edge of table (70% across canvas)
                        y      = 0.62,   # bottom edge of table (62% up canvas)
                        width  = 0.27,   # table takes 27% of canvas width
                        height = 0.32)   # table takes 32% of canvas height
    # ==========================================
    #  STEP 12: SAVE OUTPUTS
    ggsave(
    file.path(sample_dir, paste0(sample, ".pdf")),
    p_final,
    width = 20, height = 10
    )

  write.csv(
    cluster_summary,
    file.path(sample_dir, paste0(sample, "_cluster_summary.csv")),
    row.names = FALSE
    )
 message("Finished sample: ", sample)
    return(TRUE)

  }, error = function(e) {

    message("❌ ERROR in sample: ", sample)
    message("Reason: ", e$message)

    return(NULL)

  })

}

# ==========================================
# STEP 12: RUN PIPELINE
# ==========================================
sample_ids <- sub("_T.tsv$", "", list.files(pattern = "_T.tsv$"))

# test first
sample_ids <- sample_ids[1:3]

purrr::map(sample_ids, run_cladogram)
