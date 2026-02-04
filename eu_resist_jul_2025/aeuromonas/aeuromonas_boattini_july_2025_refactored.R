# aeuromonas_boattini_july_2025_refactored.R
# Refactored from original exploratory script to a modular, reusable pipeline.
# Outputs: cleaned MIC/ZONE tables, resistance summaries, heatmaps and bubble plots in ../out/

library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pheatmap)

# --- CONFIG ---------------------------------------------------------------
CFG <- list(
  raw_xlsx = "../data/2.AEuROMONAS_dataset_matteo_only.xlsx",
  lookup_csv = "../data/Look_up_AEuROMONAS_Complete_dataset.csv",
  out_dir = "../out",
  mic_col_pattern = " MIC",
  zone_col_pattern = " inhib zone diam",
  min_samples = 10,
  heatmap_width = 7,
  heatmap_height = 7,
  bubble_width = 8,
  bubble_height = 8
)

# Ensure output dir exists
if (!dir.exists(CFG$out_dir)) dir.create(CFG$out_dir, recursive = TRUE)

# --- Helpers --------------------------------------------------------------
load_lookup <- function(path = CFG$lookup_csv) {
  dt <- tryCatch(fread(path), error = function(e) stop("lookup not found: ", path))
  return(dt)
}

safe_read_raw <- function(path = CFG$raw_xlsx, sheet = "1. Complete dataset") {
  df <- read_xlsx(path = path, sheet = sheet, col_types = "text", trim_ws = TRUE)
  return(as.data.frame(df))
}

gather_mic <- function(df, mic_pattern = CFG$mic_col_pattern) {
  mic_cols <- grep(mic_pattern, names(df), value = TRUE)
  res <- df %>% select(`Code event`, `Species identification`, all_of(mic_cols)) %>%
    gather(Abx, MIC, all_of(mic_cols)) %>%
    mutate(Abx = str_replace(Abx, fixed(mic_pattern), ""))

  # common token cleaning
  res <- res %>% mutate(MIC = ifelse(MIC %in% c("≤1", "<=1"), "0.1", MIC))
  res$MIC <- suppressWarnings(as.numeric(res$MIC))
  return(res)
}

gather_zone <- function(df, zone_pattern = CFG$zone_col_pattern) {
  zone_cols <- grep(zone_pattern, names(df), value = TRUE)
  res <- df %>% select(`Code event`, `Species identification`, all_of(zone_cols)) %>%
    gather(Abx, ZONE, all_of(zone_cols)) %>%
    mutate(Abx = str_replace(Abx, fixed(zone_pattern), ""))

  res$ZONE <- suppressWarnings(as.numeric(res$ZONE))
  return(res)
}

flag_resistance_mic <- function(mic_df, lookup) {
  lk <- lookup %>% as.data.frame()
  mic_df <- mic_df %>% left_join(lk %>% select(`Species identification`, Abx, EUCAST_mic_big, CLSI_mic_bigeq_m100, CASFM_mic_big),
                                 by = c("Species identification" = "Species identification", "Abx" = "Abx"))
  mic_df <- mic_df %>%
    mutate(EUCAST_mic_Resist = ifelse(!is.na(MIC) & !is.na(EUCAST_mic_big) & MIC > EUCAST_mic_big, 1, ifelse(!is.na(MIC), 0, NA)),
           CLIST_mic_m100_Resist = ifelse(!is.na(MIC) & !is.na(CLSI_mic_bigeq_m100) & MIC >= CLSI_mic_bigeq_m100, 1, ifelse(!is.na(MIC), 0, NA)),
           CASFM_mic_Resist = ifelse(!is.na(MIC) & !is.na(CASFM_mic_big) & MIC > CASFM_mic_big, 1, ifelse(!is.na(MIC), 0, NA)))
  return(mic_df)
}

flag_resistance_zone <- function(zone_df, lookup) {
  lk <- lookup %>% as.data.frame()
  zone_df <- zone_df %>% left_join(lk %>% select(`Species identification`, Abx, EUCAST_disc_lower, CLSI_disc_lower_m100, CASFM_disc_lower),
                                   by = c("Species identification" = "Species identification", "Abx" = "Abx"))
  zone_df <- zone_df %>%
    mutate(EUCAST_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(EUCAST_disc_lower) & ZONE < EUCAST_disc_lower, 1, ifelse(!is.na(ZONE), 0, NA)),
           CLIST_Diam_m100_Resist = ifelse(!is.na(ZONE) & !is.na(CLSI_disc_lower_m100) & ZONE < CLSI_disc_lower_m100, 1, ifelse(!is.na(ZONE), 0, NA)),
           CASFM_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(CASFM_disc_lower) & ZONE < CASFM_disc_lower, 1, ifelse(!is.na(ZONE), 0, NA)))
  return(zone_df)
}

make_resistance_summary <- function(df, flag_col, group_cols = c("Species identification", "Abx"), min_samples = CFG$min_samples) {
  temp <- df %>% select(all_of(group_cols), all_of(flag_col)) %>% drop_na() %>%
    group_by_at(group_cols) %>%
    summarise(n_samples = n(), resistance_rate = mean(.data[[flag_col]], na.rm = TRUE) * 100) %>%
    ungroup() %>% filter(n_samples > min_samples)
  return(temp)
}

save_heatmap_and_bubble <- function(resistance_summary, filename_prefix, out_dir = CFG$out_dir, min_samples = CFG$min_samples) {
  heatmap_data <- resistance_summary %>% select(`Species identification`, Abx, resistance_rate) %>%
    pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

  if (nrow(heatmap_data) > 1) {
    heatmap_matrix <- as.matrix(heatmap_data[, -1])
    rownames(heatmap_matrix) <- heatmap_data$`Species identification`
    heatmap_matrix[is.na(heatmap_matrix)] <- -99
    display_matrix <- round(heatmap_matrix, 0)

    ph <- pheatmap(heatmap_matrix,
                   color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   na_col = "white",
                   number_color = "black",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   display_numbers = display_matrix,
                   main = paste0(filename_prefix, " \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >", min_samples, " Samples] "))

    svg(file.path(out_dir, paste0("dendo_", filename_prefix, ".svg")), width = CFG$heatmap_width, height = CFG$heatmap_height)
    grid::grid.newpage(); grid::grid.draw(ph$gtable)
    dev.off()
  }

  # Bubble
  bplot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
    geom_point(aes(size = n_samples, color = resistance_rate), alpha = 0.9) +
    scale_size(range = c(1, 20)) +
    scale_color_gradient(low = "lightgray", high = "midnightblue") +
    labs(title = paste0(filename_prefix, " \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >", min_samples, " Samples]"),
         x = "Antibiotic", y = "Species", size = "Sample Size", color = "Proportion Resistant") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  ggsave(file = file.path(out_dir, paste0("bubble_", filename_prefix, ".svg")), plot = bplot, width = CFG$bubble_width, height = CFG$bubble_height)
}

# --- Main pipeline --------------------------------------------------------
main <- function() {
  message("Loading lookup and raw dataset...")
  lookup <- load_lookup(CFG$lookup_csv)
  raw <- safe_read_raw(CFG$raw_xlsx)

  message("Processing MIC data...")
  mic_df <- gather_mic(raw)
  mic_df_flagged <- flag_resistance_mic(mic_df, lookup)
  fwrite(as.data.table(mic_df_flagged), file.path(CFG$out_dir, "MIC_data_refactored.csv"))

  message("Creating MIC summaries and plots (EUCAST/CLSI/CASFM)...")
  mic_eu_summary <- make_resistance_summary(mic_df_flagged, "EUCAST_mic_Resist")
  save_heatmap_and_bubble(mic_eu_summary, "10plus_mic_eucast")

  mic_clsi_summary <- make_resistance_summary(mic_df_flagged, "CLIST_mic_m100_Resist")
  save_heatmap_and_bubble(mic_clsi_summary, "10plus_mic_clsi_m100")

  mic_casfm_summary <- make_resistance_summary(mic_df_flagged, "CASFM_mic_Resist")
  save_heatmap_and_bubble(mic_casfm_summary, "10plus_mic_casfm")

  message("Processing zone inhibition data...")
  zone_df <- gather_zone(raw)
  zone_df_flagged <- flag_resistance_zone(zone_df, lookup)
  fwrite(as.data.table(zone_df_flagged), file.path(CFG$out_dir, "ZONE_data_refactored.csv"))

  zone_eu_summary <- make_resistance_summary(zone_df_flagged, "EUCAST_Diam_Resist")
  save_heatmap_and_bubble(zone_eu_summary, "10plus_zone_inhib_eucast")

  zone_casfm_summary <- make_resistance_summary(zone_df_flagged, "CASFM_Diam_Resist")
  save_heatmap_and_bubble(zone_casfm_summary, "10plus_zone_inhibition_casfm")

  message("Writing summary aggregates (collapsed) ...")
  # collapsed summaries similar to original scripts: per-Abx aggregates
  mic_summary_collapsed <- mic_df_flagged %>% filter(!is.na(MIC)) %>% group_by(Abx) %>%
    summarise(mean_MIC = mean(MIC, na.rm = TRUE), median_MIC = median(MIC, na.rm = TRUE), n = n())
  fwrite(as.data.table(mic_summary_collapsed), file.path(CFG$out_dir, "MIC_summary_collapsed.csv"))

  zone_summary_collapsed <- zone_df_flagged %>% filter(!is.na(ZONE)) %>% group_by(Abx) %>%
    summarise(mean_ZONE = mean(ZONE, na.rm = TRUE), median_ZONE = median(ZONE, na.rm = TRUE), n = n())
  fwrite(as.data.table(zone_summary_collapsed), file.path(CFG$out_dir, "ZONE_summary_collapsed.csv"))

  message("Done. Outputs in: ", CFG$out_dir)
}

# Run when executed directly
if (identical(environment(), globalenv())) {
  tryCatch(
    main(),
    error = function(e) message("Pipeline failed: ", e$message)
  )
}
