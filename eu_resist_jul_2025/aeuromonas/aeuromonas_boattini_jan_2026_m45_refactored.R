
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pheatmap)

# --- CONFIG ----------------------------------------------------------------
CFG <- list(
  raw_xlsx = "../data/2.AEuROMONAS_dataset_matteo_only.xlsx",
  lookup_csv = "../data/Look_up_AEuROMONAS_Complete_dataset.csv",
  lookup_m45_xlsx = "../data/lookup_CLSI_m45_rev_matteo.xlsx", # optional extra lookup for m45
  out_dir = "../out",
  mic_pattern = " MIC",
  zone_pattern = " inhib zone diam",
  min_samples = 10,
  heatmap_w = 7,
  heatmap_h = 7,
  bubble_w = 8,
  bubble_h = 8
)

if (!dir.exists(CFG$out_dir)) dir.create(CFG$out_dir, recursive = TRUE)

# --- Helpers ----------------------------------------------------------------
load_lookup <- function(lookup_path = CFG$lookup_csv, lookup_m45_path = CFG$lookup_m45_xlsx) {
  lk <- tryCatch(fread(lookup_path), error = function(e) stop("lookup CSV missing: ", lookup_path))

  # if an m45 lookup workbook exists, join its m45 columns (Abx, CLSI_mic_bigeq_m45, CLSI_disc_lower_m45)
  if (!is.null(lookup_m45_path) && file.exists(lookup_m45_path)) {
    m45 <- read_xlsx(path = lookup_m45_path, col_types = "text", trim_ws = TRUE)
    m45 <- m45 %>% select(Abx, CLSI_mic_bigeq_m45, CLSI_disc_lower_m45) %>% distinct()
    # ensure numeric columns where possible
    m45$CLSI_mic_bigeq_m45 <- suppressWarnings(as.numeric(m45$CLSI_mic_bigeq_m45))
    m45$CLSI_disc_lower_m45 <- suppressWarnings(as.numeric(m45$CLSI_disc_lower_m45))
    lk <- lk %>% left_join(m45, by = c("Abx" = "Abx"))
  }
  return(lk)
}

safe_read_raw <- function(path = CFG$raw_xlsx, sheet = "1. Complete dataset") {
  df <- read_xlsx(path = path, sheet = sheet, col_types = "text", trim_ws = TRUE)
  return(as.data.frame(df))
}

gather_mic <- function(df, mic_pat = CFG$mic_pattern) {
  mic_cols <- grep(mic_pat, names(df), value = TRUE)
  res <- df %>% select(`Code event`, `Species identification`, all_of(mic_cols)) %>%
    gather(Abx, MIC, all_of(mic_cols)) %>%
    mutate(Abx = str_replace(Abx, fixed(mic_pat), ""))

  # standard token fixes
  res <- res %>% mutate(MIC = ifelse(MIC %in% c("≤1", "<=1"), "0.1", MIC))
  res$MIC <- suppressWarnings(as.numeric(res$MIC))
  return(res)
}

gather_zone <- function(df, zone_pat = CFG$zone_pattern) {
  zone_cols <- grep(zone_pat, names(df), value = TRUE)
  res <- df %>% select(`Code event`, `Species identification`, all_of(zone_cols)) %>%
    gather(Abx, ZONE, all_of(zone_cols)) %>%
    mutate(Abx = str_replace(Abx, fixed(zone_pat), ""))

  res$ZONE <- suppressWarnings(as.numeric(res$ZONE))
  return(res)
}

flag_resistance_mic <- function(mic_df, lookup) {
  lk <- lookup %>% as.data.frame()
  mic_df <- mic_df %>% left_join(lk %>% select(`Species identification`, Abx, EUCAST_mic_big, CLSI_mic_bigeq_m100, CLSI_mic_bigeq_m45, CASFM_mic_big),
                                 by = c("Species identification" = "Species identification", "Abx" = "Abx"))

  mic_df <- mic_df %>%
    mutate(
      EUCAST_mic_Resist = ifelse(!is.na(MIC) & !is.na(EUCAST_mic_big) & MIC > EUCAST_mic_big, 1, ifelse(!is.na(MIC), 0, NA)),
      CLIST_mic_m100_Resist = ifelse(!is.na(MIC) & !is.na(CLSI_mic_bigeq_m100) & MIC >= CLSI_mic_bigeq_m100, 1, ifelse(!is.na(MIC), 0, NA)),
      CLIST_mic_m45_Resist = ifelse(!is.na(MIC) & !is.na(CLSI_mic_bigeq_m45) & MIC >= CLSI_mic_bigeq_m45, 1, ifelse(!is.na(MIC), 0, NA)),
      CASFM_mic_Resist = ifelse(!is.na(MIC) & !is.na(CASFM_mic_big) & MIC > CASFM_mic_big, 1, ifelse(!is.na(MIC), 0, NA))
    )
  return(mic_df)
}

flag_resistance_zone <- function(zone_df, lookup) {
  lk <- lookup %>% as.data.frame()
  zone_df <- zone_df %>% left_join(lk %>% select(`Species identification`, Abx, EUCAST_disc_lower, CLSI_disc_lower_m100, CLSI_disc_lower_m45, CASFM_disc_lower),
                                   by = c("Species identification" = "Species identification", "Abx" = "Abx"))
  zone_df <- zone_df %>%
    mutate(
      EUCAST_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(EUCAST_disc_lower) & ZONE < EUCAST_disc_lower, 1, ifelse(!is.na(ZONE), 0, NA)),
      CLIST_Diam_m100_Resist = ifelse(!is.na(ZONE) & !is.na(CLSI_disc_lower_m100) & ZONE < CLSI_disc_lower_m100, 1, ifelse(!is.na(ZONE), 0, NA)),
      CLIST_Diam_m45_Resist = ifelse(!is.na(ZONE) & !is.na(CLSI_disc_lower_m45) & ZONE < CLSI_disc_lower_m45, 1, ifelse(!is.na(ZONE), 0, NA)),
      CASFM_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(CASFM_disc_lower) & ZONE < CASFM_disc_lower, 1, ifelse(!is.na(ZONE), 0, NA))
    )
  return(zone_df)
}

make_resistance_summary <- function(df, flag_col, group_cols = c("Species identification", "Abx"), min_samples = CFG$min_samples) {
  temp <- df %>% select(all_of(group_cols), all_of(flag_col)) %>% drop_na() %>%
    group_by_at(group_cols) %>%
    summarise(n_samples = n(), resistance_rate = mean(.data[[flag_col]], na.rm = TRUE) * 100) %>%
    ungroup() %>% filter(n_samples > min_samples)
  return(temp)
}

save_heatmap_and_bubble <- function(resistance_summary, prefix, out_dir = CFG$out_dir, min_samples = CFG$min_samples) {
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
                   main = paste0(prefix, " \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >", min_samples, " Samples] "))

    svg(file.path(out_dir, paste0("dendo_", prefix, ".svg")), width = CFG$heatmap_w, height = CFG$heatmap_h)
    grid::grid.newpage(); grid::grid.draw(ph$gtable)
    dev.off()
  }

  bplot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
    geom_point(aes(size = n_samples, color = resistance_rate), alpha = 0.9) +
    scale_size(range = c(1, 20)) +
    scale_color_gradient(low = "lightgray", high = "midnightblue") +
    labs(title = paste0(prefix, " \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >", min_samples, " Samples]"),
         x = "Antibiotic", y = "Species", size = "Sample Size", color = "Proportion Resistant") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  ggsave(file = file.path(out_dir, paste0("bubble_", prefix, ".svg")), plot = bplot, width = CFG$bubble_w, height = CFG$bubble_h)
}

# --- Main pipeline ---------------------------------------------------------
main <- function() {
  message("Loading lookup (including optional m45) and raw dataset...")
  lookup <- load_lookup(CFG$lookup_csv, CFG$lookup_m45_xlsx)
  raw <- safe_read_raw(CFG$raw_xlsx)

  # MIC processing
  message("Processing MIC data and applying EUCAST/CLSI m45/m100/CASFM flags...")
  mic_df <- gather_mic(raw)
  mic_flagged <- flag_resistance_mic(mic_df, lookup)
  fwrite(as.data.table(mic_flagged), file.path(CFG$out_dir, "MIC_data_jan_2026_m45_refactored.csv"))

  # per-species MIC summaries and plots
  mic_eu <- make_resistance_summary(mic_flagged, "EUCAST_mic_Resist")
  save_heatmap_and_bubble(mic_eu, "10plus_mic_eucast_m45")

  mic_clsi_m45 <- make_resistance_summary(mic_flagged, "CLIST_mic_m45_Resist")
  save_heatmap_and_bubble(mic_clsi_m45, "10plus_mic_clsi_m45")

  mic_casfm <- make_resistance_summary(mic_flagged, "CASFM_mic_Resist")
  save_heatmap_and_bubble(mic_casfm, "10plus_mic_casfm")

  # Zone processing
  message("Processing zone inhibition data and applying diameter cutoffs (m45 where present)...")
  zone_df <- gather_zone(raw)
  zone_flagged <- flag_resistance_zone(zone_df, lookup)
  fwrite(as.data.table(zone_flagged), file.path(CFG$out_dir, "ZONE_data_jan_2026_m45_refactored.csv"))

  zone_eu <- make_resistance_summary(zone_flagged, "EUCAST_Diam_Resist")
  save_heatmap_and_bubble(zone_eu, "10plus_zone_inhib_eucast_m45")

  zone_clsi_m45 <- make_resistance_summary(zone_flagged, "CLIST_Diam_m45_Resist")
  save_heatmap_and_bubble(zone_clsi_m45, "10plus_zone_clsi_m45")

  zone_casfm <- make_resistance_summary(zone_flagged, "CASFM_Diam_Resist")
  save_heatmap_and_bubble(zone_casfm, "10plus_zone_casfm")

  # collapsed summaries per-antibiotic
  message("Writing collapsed per-antibiotic summaries...")
  mic_collapsed <- mic_flagged %>% filter(!is.na(MIC)) %>% group_by(Abx) %>% summarise(mean_MIC = mean(MIC, na.rm = TRUE), median_MIC = median(MIC, na.rm = TRUE), n = n())
  fwrite(as.data.table(mic_collapsed), file.path(CFG$out_dir, "MIC_summary_collapsed_jan_2026_m45.csv"))

  zone_collapsed <- zone_flagged %>% filter(!is.na(ZONE)) %>% group_by(Abx) %>% summarise(mean_ZONE = mean(ZONE, na.rm = TRUE), median_ZONE = median(ZONE, na.rm = TRUE), n = n())
  fwrite(as.data.table(zone_collapsed), file.path(CFG$out_dir, "ZONE_summary_collapsed_jan_2026_m45.csv"))

  message("Pipeline complete. Outputs written to: ", CFG$out_dir)
}

# run when executed directly
if (identical(environment(), globalenv())) {
  tryCatch(
    main(),
    error = function(e) message("Pipeline failed: ", e$message)
  )
}
