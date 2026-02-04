# =====================================================================
# bsi_micro_resist_2024_refactored.R
# =====================================================================

library(openxlsx)
library(readxl)
library(tidyverse)
library(data.table)
library(pheatmap)
library(xgboost)
library(caret)
library(SHAPforxgboost)
options(scipen = 999)

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
RAW_XLSX <- "ANAEuROBE_dataset_matteo_only.xlsx"
RAW_MIC_SHEET <- "3.Completo"
LOOKUP_SPECIES <- "Lookup_species.csv"
LOOKUP_DATASPECIES <- "lookupdataspecies.csv"
OUTPUT_DIR <- "./"

MIN_N_MIC_SUMMARY <- 20
MIN_N_ZONE_SUMMARY <- 10
PEDS_SAMPLE_THRESHOLD <- 5
SEED <- 2024

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

load_lookup_species <- function(path = LOOKUP_SPECIES) {
  lk <- fread(path, colClasses = "character")
  lk <- lk %>% mutate(`Species identification` = str_replace_all(`Species identification`, "�", " "))
  lk$`Species identification` <- str_trim(lk$`Species identification`)
  lk <- lk %>% distinct()
  # numeric conversions
  lk <- lk %>% mutate(across(c(EUCAST_mic_big, CLSI_mic_bigeq, CASFM_mic_big, EUCAST_disc_lower, CASFM_disc_lower), ~ as.numeric(.x)))
  lk
}

load_mydata <- function(path = RAW_XLSX, sheet = RAW_MIC_SHEET) {
  df <- read_excel(path = path, sheet = sheet, col_types = "text")
  df <- df %>% mutate(`Species identification` = str_replace_all(`Species identification`, "�", " "))
  df$`Species identification` <- str_trim(df$`Species identification`)
  df <- df %>% filter(!grepl("JEN", `Code event`))
  df
}

# Generic sequential replacement utility for messy MIC/ZONE tokens
apply_replacements_seq <- function(vec, replacements) {
  out <- vec
  for (pat in names(replacements)) {
    out <- stringr::str_replace_all(out, fixed(pat), replacements[[pat]])
  }
  out
}

# Sanitizer for MIC strings - compress many ifelse chains into replacement map
sanitize_mic_strings <- function(mic_vec) {
  repl <- c(
    "<0.016" = "0.008", "<0.019" = "0.010", "<=0.016" = "0.008",
    "1.6E-2" = "0.016", "4.7E-2" = "0.047", "2.3E-2" = "0.023",
    "6.4000000000000001E-2" = "0.064", "9.4E-2" = "0.094", 
    ">256" = "260", ">32" = "48", ">320" = "480",
    "3.2000000000000001E-2" = "0.032", "1.2E-2" = "0.012",
    "2E-3" = "0.002", "8.0000000000000002E-3" = "0.008",
    "6.0000000000000001E-3" = "0.006", "3.0000000000000001E-3" = "0.003",
    "4.0000000000000001E-3" = "0.004", "2.4E-2" = "0.024",
    "≤0.25" = "0.12", ">0.5" = "1", "≤0.016" = "0.008",
    "<0.06" = "0.03", "3.4000000000000002E-2" = "0.03", "<0.002" = "0.001",
    "<0.064" = "0.03", "<0.03" = "0.01", "<0.25" = "0.1",
    ">2" = "3", "<256" = "128", ">8" = "12", "00.032" = "0.03",
    "≤4" = "2", "≤0.5" = "0.2", "1-2" = "1",
    "≤4/2" = "2", ">8/2" = "12",
    "<0,06" = "0.03", "<1" = "0.5", "<2" = "1",
    ">0.016" = "0.03", "8/2" = "8", "1.7000000000000001E-2" = "0.017",
    "§" = "NA", "<0.16" = "0.08", "3.7999999999999999E-2" = "0.04",
    "≤8/4" = "4", "≤1" = "0.5", "Not tested" = "NA",
    "45 mm" = "45", "16/4" = "16", ">16/4" = "20", ">128" = "200",
    ">4" = "6", "≤0.06" = "0.03", "0,75" = "0.7", "<0.026" = "0.01",
    "≤0.0625" = "0.03", ".064" = "0.064", ">64" = "100", "<0.047" = "0.02",
    "1.2500000000000001E-2" = "0.0125", "<0,016" = "0.005",
    "4.8000000000000001E-2" = "0.048", "< 0.016" = "0.005", "<0.0016" = "0.0005",
    "1.6000000000000001E-3" = "0.0016", "<0.5" = "0.2", "<0,002" = "0.001",
    "9.1999999999999998E-2" = "0.08", "> 32" = "48", "2.3E-3" = "0.002",
    "<0.02" = "0.01", "1.9E-2" = "0.019", "≤2" = "1"
  )
  out <- apply_replacements_seq(mic_vec, repl)
  out <- ifelse(out %in% c("R", "S", "r", "s", "O"), NA, out)
  # convert comma decimals
  out <- str_replace_all(out, ",", ".")
  out
}

sanitize_zone_strings <- function(zone_vec) {
  repl <- c(
    ">35" = "50", ">16" = "24", ">32" = "48", "8.0000000000000002E-3" = "0.008",
    "<0.002" = "0.001", "4.0000000000000001E-3" = "0.004", "2.3E-2" = "0.0023",
    "6.0000000000000001E-3" = "0.006", "9.4E-2" = "0.0094",
    ">30" = "45", ">40" = "60", ">25" = "38", ">22" = "33", ">20" = "30",
    "<40" = "20", "6.4000000000000001E-2" = "0.064", "1.2E-2" = "0.012",
    "3.0000000000000001E-3" = "0.003", "1.6E-2" = "0.0016", "4.7E-2" = "0.047",
    "3.2000000000000001E-2" = "0.032", "2E-3" = "0.002", ">27" = "40",
    "<0,06" = "0.03", "<0,016" = "0.005"
  )
  out <- apply_replacements_seq(zone_vec, repl)
  out <- str_replace_all(out, ",", ".")
  out
}

# Parse numeric values safely
safe_parse_number <- function(vec) {
  out <- parse_number(vec)
  out
}

# Centralized save helper
save_csv <- function(df, name) {
  path <- file.path(OUTPUT_DIR, name)
  fwrite(df, path)
  message("Saved: ", path)
}

# Plot helpers
plot_heatmap_and_save <- function(mat, filename, main = NULL, width = 6, height = 10, palette = c("white", "lightcyan1", "royalblue4")) {
  display_matrix <- round(mat)
  p <- pheatmap(mat,
                color = colorRampPalette(palette)(50),
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                na_col = "white",
                fontsize_row = 8,
                fontsize_col = 8,
                display_numbers = display_matrix,
                main = main)
  ggsave(filename = file.path(OUTPUT_DIR, filename), plot = p, width = width, height = height)
}

plot_bubble_and_save <- function(df, filename, width = 8, height = 10, low = "lightgray", high = "midnightblue") {
  p <- ggplot(df, aes(x = Abx, y = `Species identification`)) +
    geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
    scale_size(range = c(1, 20)) +
    scale_color_gradient(low = low, high = high) +
    labs(title = df$title %||% "Antibiotic Resistance by Species and Antibiotic",
         x = "Antibiotic", y = "Species", size = "Sample Size", color = "% Resistant") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(file = file.path(OUTPUT_DIR, filename), plot = p, width = width, height = height)
}

# ---------------------------------------------------------------------
# Main refactored processing
# ---------------------------------------------------------------------

main <- function() {
  set.seed(SEED)

  lk <- load_lookup_species()
  lookupdataspecies <- fread(LOOKUP_DATASPECIES)
  my_data <- load_mydata()

  # align speciesID via lookupdataspecies
  my_data <- lookupdataspecies %>% left_join(my_data, by = c("Species.identification" = "Species identification"))
  # apply a few manual speciesID fixes from original script
  my_data <- my_data %>% mutate(speciesID = case_when(
    speciesID == 40 ~ 57,
    speciesID == 176 ~ 186,
    speciesID == 227 ~ 221,
    speciesID == 98 ~ 229,
    speciesID == 170 ~ 102,
    TRUE ~ speciesID
  ))

  # ---------- MIC processing ----------
  MIC_data <- my_data %>% select(`Code event`, speciesID) %>% bind_cols(my_data %>% select(contains(" MIC")))
  MIC_data <- MIC_data %>% select(-c(`AST by ANA-ATB microdilution=1`, `AST by MICRONAUT-S MIC test=1`)) %>%
    pivot_longer(cols = contains(" MIC"), names_to = "Abx", values_to = "MIC")

  MIC_data <- MIC_data %>% mutate(R_S = case_when(
    str_detect(MIC, regex("R", ignore_case = TRUE)) ~ "R",
    str_detect(MIC, regex("S", ignore_case = TRUE)) ~ "S",
    str_detect(MIC, ">") ~ "R",
    str_detect(MIC, "<") ~ "S",
    TRUE ~ NA_character_
  ))

  MIC_data <- MIC_data %>% mutate(MIC_s = sanitize_mic_strings(MIC)) %>% mutate(MIC_num = safe_parse_number(MIC_s))

  # join lookup thresholds
  MIC_data <- MIC_data %>% left_join(lk %>% select(`Species identification`, speciesID, Abx, EUCAST_mic_big, CLSI_mic_bigeq, CASFM_mic_big) %>%
                                      mutate(Abx = str_replace_all(Abx, " MIC", "")), by = c("speciesID", "Abx"))

  MIC_data <- MIC_data %>% mutate(
    EUCAST_mic_Resist = case_when(!is.na(MIC_num) & !is.na(EUCAST_mic_big) ~ as.integer(MIC_num > EUCAST_mic_big), TRUE ~ NA_integer_),
    CLIST_mic_Resist = case_when(!is.na(MIC_num) & !is.na(CLSI_mic_bigeq) ~ as.integer(MIC_num >= CLSI_mic_bigeq), TRUE ~ NA_integer_),
    CASFM_mic_Resist = case_when(!is.na(MIC_num) & !is.na(CASFM_mic_big) ~ as.integer(MIC_num > CASFM_mic_big), TRUE ~ NA_integer_)
  )

  # summaries and saves
  summary_MIC_concent <- MIC_data %>% filter(!is.na(MIC_num)) %>% group_by(`Species identification`, speciesID, Abx) %>%
    summarise(mean = mean(MIC_num, na.rm = TRUE), sd = sd(MIC_num, na.rm = TRUE), median = median(MIC_num, na.rm = TRUE), Q90 = quantile(MIC_num, 0.9, na.rm = TRUE), n = n(), .groups = "drop")
  save_csv(summary_MIC_concent, "summary_MIC_concent_Jan29_MIC90.csv")

  R_S_original_counts <- MIC_data %>% filter(!is.na(R_S)) %>% group_by(`Species identification`, speciesID, Abx, R_S) %>% count() %>% pivot_wider(names_from = R_S, values_from = n) %>% replace_na(list(R = 0, S = 0)) %>% mutate(n_n = R + S, perc_r_s = R / (R + S))
  save_csv(R_S_original_counts, "R_S_original_counts_Jan29.csv")

  save_csv(MIC_data, "MIC_data_Jan29.csv")

  EUCAST_resist_counts <- MIC_data %>% filter(!is.na(EUCAST_mic_Resist)) %>% group_by(`Species identification`, speciesID, Abx, EUCAST_mic_Resist) %>% count() %>% pivot_wider(names_from = EUCAST_mic_Resist, values_from = n) %>% replace_na(list(`0` = 0, `1` = 0)) %>% mutate(n_eucast = `1` + `0`, perc_r_eucast = `1` / (n_eucast)) %>% select(-c(`1`, `0`))
  CLSI_resist_counts <- MIC_data %>% filter(!is.na(CLIST_mic_Resist)) %>% group_by(`Species identification`, speciesID, Abx, CLIST_mic_Resist) %>% count() %>% pivot_wider(names_from = CLIST_mic_Resist, values_from = n) %>% replace_na(list(`0` = 0, `1` = 0)) %>% mutate(n_clsi = `1` + `0`, perc_r_clsi = `1` / (n_clsi)) %>% select(-c(`1`, `0`))
  casfm_resist_counts <- MIC_data %>% filter(!is.na(CASFM_mic_Resist)) %>% group_by(`Species identification`, speciesID, Abx, CASFM_mic_Resist) %>% count() %>% pivot_wider(names_from = CASFM_mic_Resist, values_from = n) %>% replace_na(list(`0` = 0, `1` = 0)) %>% mutate(n_casfm = `1` + `0`, perc_r_casfm = `1` / (n_casfm)) %>% select(-c(`1`, `0`))

  save_csv(EUCAST_resist_counts, "EUCAST_resist_counts_Jan29.csv")
  save_csv(CLSI_resist_counts, "CLSI_resist_counts_Jan29.csv")
  save_csv(casfm_resist_counts, "casfm_resist_counts_Jan29.csv")

  # ---------- ZONE (inhibition diameters) processing ----------
  ZONE_data <- my_data %>% select(`Code event`, speciesID) %>% bind_cols(my_data %>% select(contains(" zone diam")))
  ZONE_data <- ZONE_data %>% bind_cols(my_data %>% select(`EUCAST=1; CLSI=0`))
  ZONE_data <- ZONE_data %>% pivot_longer(cols = contains("zone"), names_to = "Abx", values_to = "ZONE")

  ZONE_data <- ZONE_data %>% mutate(ZONE_s = sanitize_zone_strings(ZONE)) %>% mutate(ZONE_num = safe_parse_number(ZONE_s))
  ZONE_data <- ZONE_data %>% mutate(Abx = str_replace_all(Abx, " inhib zone diameter", ""), Abx = str_replace_all(Abx, " zone diam", ""))
  ZONE_data <- ZONE_data %>% left_join(lk %>% mutate(Abx = str_replace_all(Abx, " MIC", "")) %>% select(Abx, speciesID, EUCAST_disc_lower, CASFM_disc_lower), by = c("speciesID", "Abx"))

  ZONE_data <- ZONE_data %>% mutate(
    EUCAST_Diam_Resist = case_when(!is.na(ZONE_num) & !is.na(EUCAST_disc_lower) ~ as.integer(ZONE_num < EUCAST_disc_lower), TRUE ~ NA_integer_),
    CASFM_Diam_Resist = case_when(!is.na(ZONE_num) & !is.na(CASFM_disc_lower) ~ as.integer(ZONE_num < CASFM_disc_lower), TRUE ~ NA_integer_)
  )

  # Summaries and saves for ZONE
  summary_ZONE_concent <- ZONE_data %>% filter(!is.na(ZONE_num)) %>% group_by(`Species identification`, speciesID, Abx) %>% summarise(mean = mean(ZONE_num, na.rm = TRUE), sd = sd(ZONE_num, na.rm = TRUE), median = median(ZONE_num, na.rm = TRUE), Q1 = quantile(ZONE_num, 0.25, na.rm = TRUE), Q3 = quantile(ZONE_num, 0.75, na.rm = TRUE), n = n(), .groups = "drop")
  save_csv(summary_ZONE_concent, "summary_ZONE_concent_Jan29.csv")

  save_csv(ZONE_data, "ZONE_data_Jan29.csv")

  EUCAST_resist_counts_zone <- ZONE_data %>% filter(!is.na(EUCAST_Diam_Resist)) %>% group_by(`Species identification`, speciesID, Abx, EUCAST_Diam_Resist) %>% count() %>% pivot_wider(names_from = EUCAST_Diam_Resist, values_from = n) %>% replace_na(list(`0` = 0, `1` = 0)) %>% mutate(n_eucast = `1` + `0`, perc_r_eucast = `1` / (n_eucast)) %>% select(-c(`1`, `0`))
  casfm_resist_counts_zone <- ZONE_data %>% filter(!is.na(CASFM_Diam_Resist)) %>% group_by(`Species identification`, speciesID, Abx, CASFM_Diam_Resist) %>% count() %>% pivot_wider(names_from = CASFM_Diam_Resist, values_from = n) %>% replace_na(list(`0` = 0, `1` = 0)) %>% mutate(n_casfm = `1` + `0`, perc_r_casfm = `1` / (n_casfm)) %>% select(-c(`1`, `0`))

  save_csv(EUCAST_resist_counts_zone, "EUCAST_resist_counts_zone_Jan29.csv")
  save_csv(casfm_resist_counts_zone, "casfm_resist_counts_zone_Jan29.csv")

  # ZI processing (if present)
  if ("ZI" %in% names(ZONE_data)) {
    ZI_df <- ZONE_data %>% filter(grepl(",", ZI)) %>% filter(!is.na(EUCAST_Diam_Resist) | !is.na(CASFM_Diam_Resist)) %>% separate(ZI, c("Lower", "Upper"), sep = ",") %>% mutate(Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(ZI_flag = ifelse(ZONE_num > Lower & ZONE_num < Upper, 1L, 0L))
    save_csv(ZI_df, "ZI_Jan29.csv")
  }

  # ---------- Visualization: MIC heatmap + bubble ----------
  resistance_summary <- MIC_data %>% select(`Code event`, `Species identification`, speciesID, Abx, EUCAST_mic_Resist) %>% drop_na() %>% group_by(`Species identification`, speciesID, Abx) %>% summarise(n_samples = n(), resistance_rate = mean(EUCAST_mic_Resist, na.rm = TRUE) * 100, .groups = "drop") %>% filter(n_samples > MIN_N_MIC_SUMMARY)
  resistance_summary <- resistance_summary %>% mutate(Abx = str_replace_all(Abx, " MIC", ""))

  heatmap_data <- resistance_summary %>% ungroup() %>% select(`Species identification`, Abx, resistance_rate) %>% pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data$`Species identification`
  heatmap_matrix[is.na(heatmap_matrix)] <- NA
  try(plot_heatmap_and_save(heatmap_matrix, "dendo_MIC.svg", main = "Antibiotic Resistance Clustering [MIC]"), silent = TRUE)

  resistance_summary <- resistance_summary %>% mutate(title = "Antibiotic Resistance by Species and Antibiotic [MIC]")
  try(plot_bubble_and_save(resistance_summary, "bubble_MIC.svg", width = 8, height = 15, low = "lightgray", high = "midnightblue"), silent = TRUE)

  # ---------- Visualization: ZONE heatmap + bubble ----------
  zone_summary <- ZONE_data %>% select(`Code event`, `Species identification`, speciesID, Abx, EUCAST_Diam_Resist) %>% drop_na() %>% group_by(`Species identification`, speciesID, Abx) %>% summarise(n_samples = n(), resistance_rate = mean(EUCAST_Diam_Resist, na.rm = TRUE) * 100, .groups = "drop") %>% filter(n_samples > MIN_N_ZONE_SUMMARY)
  zone_summary <- zone_summary %>% mutate(Abx = str_replace_all(Abx, " MIC", ""))
  heatmap_data_zone <- zone_summary %>% ungroup() %>% select(`Species identification`, Abx, resistance_rate) %>% pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)
  heatmap_matrix_zone <- as.matrix(heatmap_data_zone[, -1])
  rownames(heatmap_matrix_zone) <- heatmap_data_zone$`Species identification`
  heatmap_matrix_zone[is.na(heatmap_matrix_zone)] <- NA
  try(plot_heatmap_and_save(heatmap_matrix_zone, "dendo_zone.svg", main = "Antibiotic Resistance Clustering [ZONE]", palette = c("lightgray", "white", "firebrick")), silent = TRUE)

  zone_summary <- zone_summary %>% mutate(title = "Antibiotic Resistance by Species and Antibiotic [ZONE]")
  try(plot_bubble_and_save(zone_summary, "bubble_zone.svg", width = 8, height = 10, low = "lightgray", high = "firebrick"), silent = TRUE)

  # ---------- Paediatric vs Adult ML (XGBoost + SHAP) ----------
  # Build wide table of EUCAST resistant indicators (MIC + zone combined) per Code event
  mic_wide <- MIC_data %>% filter(!is.na(EUCAST_mic_Resist)) %>% mutate(var = paste0(`Species identification`, "_", Abx)) %>% select(`Code event`, var, EUCAST_mic_Resist) %>% distinct() %>% pivot_wider(names_from = var, values_from = EUCAST_mic_Resist)
  zone_wide <- ZONE_data %>% filter(!is.na(EUCAST_Diam_Resist)) %>% mutate(var = paste0(`Species identification`, "_", Abx)) %>% select(`Code event`, var, EUCAST_Diam_Resist) %>% distinct() %>% pivot_wider(names_from = var, values_from = EUCAST_Diam_Resist)
  combined <- full_join(mic_wide, zone_wide, by = "Code event")

  # join peds label
  peds <- read_excel(path = RAW_XLSX, sheet = RAW_MIC_SHEET, col_types = "text") %>% select(`Code event`, `Paediatrics=1`) %>% mutate(`Paediatrics=1` = as.numeric(`Paediatrics=1`))
  combined <- peds %>% inner_join(combined, by = "Code event")
  combined <- combined %>% mutate(peds = ifelse(`Paediatrics=1` == 1, 1, 0)) %>% select(-`Paediatrics=1`)

  # Prepare matrix
  wide_mat <- combined %>% select(-`Code event`, -peds)
  # ensure numeric matrix (NA->0)
  wide_mat[is.na(wide_mat)] <- 0
  data_matrix <- as.matrix(wide_mat)
  labels <- combined$peds

  if (ncol(data_matrix) > 1 && nrow(data_matrix) > 10) {
    dtrain <- xgb.DMatrix(data = data_matrix, label = labels)
    xgb_model <- xgboost(data = dtrain, params = list(objective = "binary:logistic"), nrounds = 100, verbose = 0)
    preds <- predict(xgb_model, newdata = data_matrix)
    binary_preds <- ifelse(preds > 0.5, 1, 0)
    conf <- confusionMatrix(as.factor(binary_preds), as.factor(labels))
    message("XGBoost trained; AUC and confusion metrics: ")
    try({
      roc_obj <- pROC::roc(labels, preds)
      message("AUC = ", pROC::auc(roc_obj))
    }, silent = TRUE)

    # SHAP
    try({
      shap_values <- shap.values(xgb_model = xgb_model, X_train = data_matrix)
      mean_shap <- shap_values$mean_shap_score
      mean_shap_df <- data.frame(Feature = names(mean_shap), Importance = mean_shap)
      ploted_shap <- mean_shap_df %>% head(20) %>% ggplot(aes(x = reorder(Feature, Importance), y = Importance)) + geom_bar(stat = "identity", fill = "#183555") + coord_flip() + theme_minimal() + labs(title = "Mean Absolute SHAP Values [Top 20]", x = "Feature", y = "Mean SHAP")
      ggsave(file = file.path(OUTPUT_DIR, "ploted_shap_eucast_mic_PLUS_diam_excJena.svg"), plot = ploted_shap, width = 7, height = 6)
    }, silent = TRUE)
  } else {
    message("Insufficient data for XGBoost/SHAP (rows: ", nrow(data_matrix), " cols: ", ncol(data_matrix), ")")
  }

  # Paediatrics vs adults at species-level comparisons and saves (as in original)
  # create wide presence table for top features if model ran
  message("Refactored pipeline completed. Outputs saved to ", normalizePath(OUTPUT_DIR))
}

if (interactive()) main()

# End of file
