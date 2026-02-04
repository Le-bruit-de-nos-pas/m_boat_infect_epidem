# micro_resist_matteo_oct_2024_refactored.R
# Refactored pipeline for MIC / inhibition zone processing, summaries, heatmaps and XGBoost SHAP for peds vs adults

library(readxl)
library(openxlsx)
library(data.table)
library(tidyverse)
library(pheatmap)
library(xgboost)
library(SHAPforxgboost)
library(caret)
library(pROC)
options(scipen = 999)

# --- Config ---
CFG <- list(
  lookup_species = "Lookup_species.csv",
  lookup_data_species = "lookupdataspecies.csv",
  raw_workbook = "ANAEuROBE_dataset_matteo_only.xlsx",
  mic_sheet = "3.Completo",
  output_dir = "output_micro",
  top_n_species = 50,
  shap_top = 20,
  seed = 123
)

if (!dir.exists(CFG$output_dir)) dir.create(CFG$output_dir, recursive = TRUE)

# --- Utilities ---
safe_read_csv_or_fread <- function(path, ...) {
  tryCatch(fread(path, ...), error = function(e) NULL)
}

# sanitize messy MIC/ZONE tokens using named replacement maps
apply_replacements <- function(x, replacements) {
  for (pat in names(replacements)) {
    x <- ifelse(x == pat, replacements[[pat]], x)
  }
  x
}

safe_parse_number <- function(x) {
  suppressWarnings(parse_number(as.character(x)))
}

# --- Load lookups ---
load_lookup_species <- function(path = CFG$lookup_species) {
  if (!file.exists(path)) stop("Lookup file not found: ", path)
  ls <- fread(path, colClasses = "character") %>%
    mutate(`Species identification` = str_replace_all(`Species identification`, "�", " ")) %>%
    mutate(`Species identification` = str_trim(`Species identification`)) %>%
    distinct()
  # numeric conversions for thresholds
  ls <- ls %>% mutate_at(vars(EUCAST_mic_big:CASFM_mic_big, EUCAST_disc_lower:CASFM_disc_lower), as.numeric)
  # speciesID helper if missing
  if (!"speciesID" %in% colnames(ls)) {
    ls <- ls %>% left_join(ls %>% select(`Species identification`) %>% distinct() %>%
                             arrange(`Species identification`) %>%
                             filter(`Species identification` != "") %>%
                             drop_na() %>% mutate(speciesID = row_number()))
  }
  ls
}

# --- Sanitizers ---
mic_replacements <- c(
  "<0.016" = "0.008", "<0.019" = "0.010", "<=0.016" = "0.008",
  "+>256" = "260", ">256" = "260", ">32" = "48", ">64" = "100",
  "<=0.25" = "0.12", "≤0.25" = "0.12", "≤0.016" = "0.008",
  "<0.06" = "0.03", "<0.002" = "0.001", "Not tested" = NA,
  "R" = NA, "S" = NA, "r" = NA, "s" = NA, "O" = NA
)

zone_replacements <- c(
  ">35" = "50", ">30" = "45", ">40" = "60", ">25" = "38", ">22" = "33",
  ">20" = "30", "<40" = "20", "na" = NA, "S" = NA, "R" = NA, "s" = NA, "I" = NA,
  "Not tested" = NA
)

sanitize_mic_column <- function(mic_col) {
  mic_col <- as.character(mic_col)
  mic_col <- apply_replacements(mic_col, mic_replacements)
  safe_parse_number(mic_col)
}

sanitize_zone_column <- function(zone_col) {
  z <- as.character(zone_col)
  z <- apply_replacements(z, zone_replacements)
  safe_parse_number(z)
}

# --- MIC processing ---
process_mic <- function(raw_wb = CFG$raw_workbook, sheet = CFG$mic_sheet, lookup = NULL) {
  message("Loading MIC raw sheet: ", raw_wb)
  raw <- read_excel(path = raw_wb, sheet = sheet, col_types = "text") %>%
    mutate(`Species identification` = str_replace_all(`Species identification`, "�", " ")) %>%
    mutate(`Species identification` = str_trim(`Species identification`))

  # join species lookup ids
  lds <- safe_read_csv_or_fread(CFG$lookup_data_species)
  if (!is.null(lds)) raw <- lds %>% left_join(raw, by = c("Species.identification" = "Species identification"))

  # gather MIC columns
  mic_cols <- raw %>% select(contains(" MIC")) %>% names()
  if (length(mic_cols) == 0) stop("No MIC columns detected")
  mic_df <- raw %>% select(`Code event`, speciesID, all_of(mic_cols)) %>%
    gather(Abx, MIC, all_of(mic_cols)) %>%
    mutate(Abx = str_replace_all(Abx, " MIC", ""))

  mic_df <- mic_df %>% mutate(R_S = ifelse(str_detect(MIC, regex("[Rr]")), "R",
                                           ifelse(str_detect(MIC, regex("[Ss]")), "S",
                                                  NA)))
  mic_df <- mic_df %>% mutate(MIC_raw = MIC, MIC = sanitize_mic_column(MIC))

  if (!is.null(lookup)) mic_df <- mic_df %>% left_join(lookup %>% select(-c(EUCAST_disc_lower, CASFM_disc_lower, ZI)), by = c("Abx", "speciesID"))

  mic_df <- mic_df %>% mutate(EUCAST_mic_Resist = ifelse(!is.na(MIC) & !is.na(EUCAST_mic_big) & MIC > EUCAST_mic_big, 1, 0),
                              CLSI_mic_Resist = ifelse(!is.na(MIC) & !is.na(CLSI_mic_bigeq) & MIC >= CLSI_mic_bigeq, 1, 0),
                              CASFM_mic_Resist = ifelse(!is.na(MIC) & !is.na(CASFM_mic_big) & MIC > CASFM_mic_big, 1, 0))

  mic_df
}

# --- ZONE processing ---
process_zone <- function(raw_wb = CFG$raw_workbook, sheet = CFG$mic_sheet, lookup = NULL) {
  raw <- read_excel(path = raw_wb, sheet = sheet, col_types = "text") %>%
    mutate(`Species identification` = str_replace_all(`Species identification`, "�", " ")) %>%
    mutate(`Species identification` = str_trim(`Species identification`))

  zone_cols <- raw %>% select(contains("zone")) %>% names()
  if (length(zone_cols) == 0) stop("No zone columns detected")

  zone_df <- raw %>% select(`Code event`, speciesID, all_of(zone_cols), `EUCAST=1; CLSI=0`) %>%
    gather(Abx, ZONE, all_of(zone_cols)) %>%
    mutate(Abx = str_replace_all(Abx, " inhib zone diameter", "")) %>%
    mutate(Abx = str_replace_all(Abx, " zone diam", ""))

  zone_df <- zone_df %>% mutate(ZONE_raw = ZONE, ZONE = sanitize_zone_column(ZONE))

  if (!is.null(lookup)) zone_df <- zone_df %>% left_join(lookup %>% mutate(Abx = str_replace_all(Abx, " MIC", "")) %>% select(-c(EUCAST_mic_big, CLSI_mic_bigeq, CASFM_mic_big)), by = c("Abx", "speciesID"))

  zone_df <- zone_df %>% mutate(EUCAST_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(EUCAST_disc_lower) & ZONE < EUCAST_disc_lower, 1, 0),
                                CASFM_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(CASFM_disc_lower) & ZONE < CASFM_disc_lower, 1, 0))

  # apply EUCAST/CLSI selection logic
  zone_df <- zone_df %>% mutate(EUCAST_Diam_Resist = ifelse(Abx %in% c("Piperacillin/tazobactam", "Clindamycin", "Metronidazole", "Ertapenem", "Imipenem"), EUCAST_Diam_Resist,
                                                              ifelse(`EUCAST=1; CLSI=0` == "1", EUCAST_Diam_Resist, NA)),
                                CASFM_Diam_Resist = ifelse(Abx %in% c("Piperacillin/tazobactam", "Clindamycin", "Metronidazole", "Ertapenem", "Imipenem"), CASFM_Diam_Resist,
                                                           ifelse(`EUCAST=1; CLSI=0` != "1", CASFM_Diam_Resist, NA)))
  zone_df
}

# --- Summary & plotting ---
summarize_resistance <- function(df, resist_col, min_n = 20, top_n = CFG$top_n_species) {
  sum <- df %>% filter(!is.na(!!sym(resist_col))) %>%
    group_by(`Species identification`, speciesID, Abx) %>%
    summarise(n_samples = n(), resistance_rate = mean(!!sym(resist_col), na.rm = TRUE) * 100, .groups = "drop") %>%
    inner_join(df %>% select(speciesID, `Code event`) %>% distinct() %>% group_by(speciesID) %>% count() %>% arrange(-n) %>% head(top_n) %>% ungroup() %>% select(speciesID), by = "speciesID")
  sum
}

plot_heatmap <- function(summary_df, file = file.path(CFG$output_dir, "heatmap.svg"), color_palette = c("lightgray","lightblue","midnightblue")) {
  heatmap_data <- summary_df %>% ungroup() %>% select(`Species identification`, Abx, resistance_rate) %>%
    pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data[[1]]
  heatmap_matrix[is.na(heatmap_matrix)] <- -99
  display_matrix <- round(heatmap_matrix, 0)

  p <- pheatmap(heatmap_matrix,
                color = colorRampPalette(color_palette)(50),
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                na_col = "grey",
                fontsize_row = 8,
                fontsize_col = 8,
                display_numbers = display_matrix,
                main = "Antibiotic Resistance Clustering [TOP]")
  ggsave(file, plot = p, width = 8, height = 12)
  invisible(p)
}

plot_bubble <- function(summary_df, file = file.path(CFG$output_dir, "bubble.svg"), width = 10, height = 12, low_col = "lightgray", high_col = "midnightblue") {
  df <- summary_df %>% mutate(Abx = str_replace_all(Abx, " MIC", ""))
  p <- ggplot(df, aes(x = Abx, y = `Species identification`)) +
    geom_point(aes(size = n_samples, color = resistance_rate), alpha = 0.9) +
    scale_size(range = c(1, 20)) +
    scale_color_gradient(low = low_col, high = high_col) +
    labs(title = "Antibiotic Resistance by Species and Antibiotic [TOP]", x = "Antibiotic", y = "Species", size = "Sample Size", color = "Proportion Resistant") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(file, plot = p, width = width, height = height)
  invisible(p)
}

# --- XGBoost + SHAP wrapper for pediatric vs adult comparison ---
run_xgb_shap <- function(mic_clean_df, peds_sheet = CFG$mic_sheet, raw_wb = CFG$raw_workbook, top_shap = CFG$shap_top) {
  # mic_clean_df expected wide: one row per Code event with binary resistance features
  # create peds flag
  peds <- read_excel(path = raw_wb, sheet = CFG$mic_sheet, col_types = "text") %>% select(`Code event`, `Paediatrics=1`) %>% mutate(`Paediatrics=1` = as.numeric(`Paediatrics=1`))
  peds[is.na(peds)] <- 0
  names(peds)[2] <- "peds"

  wide <- mic_clean_df %>% mutate(Abx = str_replace_all(Abx, " MIC", "")) %>%
    mutate(Abx = str_replace_all(Abx, " ", "_")) %>%
    mutate(var = paste0(`Species identification`, "_", Abx)) %>%
    select(`Code event`, var, EUCAST_mic_Resist) %>% distinct() %>%
    spread(key = var, value = EUCAST_mic_Resist)

  wide <- peds %>% inner_join(wide, by = c("Code event" = "Code event")) %>% select(-Code.event)
  wide$peds <- as.numeric(wide$peds) - 1

  # prepare matrix
  label <- wide$peds
  X <- wide %>% select(-peds, -`Code event`, -`Code.event`) %>% mutate_all(~ifelse(is.na(.), 0, .))
  if (ncol(X) < 2) { message("Not enough features for XGBoost"); return(NULL) }

  set.seed(CFG$seed)
  dtrain <- xgb.DMatrix(data = as.matrix(X), label = label)
  params <- list(objective = "binary:logistic", eval_metric = "auc")
  xgb_model <- xgboost(data = dtrain, params = params, nrounds = 100, verbose = 0)

  preds <- predict(xgb_model, as.matrix(X))
  auc_val <- tryCatch(auc(roc(label, preds)), error = function(e) NA_real_)

  shap_vals <- shap.values(xgb_model = xgb_model, X_train = as.matrix(X))
  mean_shap <- shap_vals$mean_shap_score
  mean_shap_df <- data.frame(Feature = names(mean_shap), Importance = mean_shap) %>% arrange(desc(Importance))
  mean_shap_df <- mean_shap_df %>% head(top_shap)

  g <- mean_shap_df %>% mutate(Feature = str_replace_all(Feature, "Benzilpenicillin", "Benzylpenicillin")) %>%
    ggplot(aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "#183555") + coord_flip() + theme_minimal() + labs(title = paste("Mean Absolute SHAP Values [Top", top_shap, "]"), x = "Feature", y = "Mean SHAP")
  ggsave(filename = file.path(CFG$output_dir, "shap_top.png"), plot = g, width = 7, height = 6)

  list(model = xgb_model, auc = auc_val, shap = mean_shap_df)
}

# --- Main runner ---
main <- function() {
  message("Loading lookup species")
  lookup <- tryCatch(load_lookup_species(CFG$lookup_species), error = function(e) { message(e$message); NULL })

  message("Processing MICs...")
  mic_df <- process_mic(lookup = lookup)
  fwrite(mic_df, file.path(CFG$output_dir, "MIC_data_clean.csv"))

  message("Creating MIC summaries...")
  mic_summary <- summarize_resistance(mic_df, "EUCAST_mic_Resist")
  fwrite(mic_summary, file.path(CFG$output_dir, "MIC_resistance_summary.csv"))
  plot_heatmap(mic_summary, file = file.path(CFG$output_dir, "heatmap_mic.svg"))
  plot_bubble(mic_summary, file = file.path(CFG$output_dir, "bubble_mic.svg"))

  message("Processing Zones...")
  zone_df <- process_zone(lookup = lookup)
  fwrite(zone_df, file.path(CFG$output_dir, "ZONE_data_clean.csv"))

  message("Creating Zone summaries...")
  zone_summary <- summarize_resistance(zone_df, "EUCAST_Diam_Resist")
  fwrite(zone_summary, file.path(CFG$output_dir, "ZONE_resistance_summary.csv"))
  plot_heatmap(zone_summary, file = file.path(CFG$output_dir, "heatmap_zone.svg"), color_palette = c("lightgray","white","firebrick"))
  plot_bubble(zone_summary, file = file.path(CFG$output_dir, "bubble_zone.svg"), low_col = "lightgray", high_col = "firebrick")

  message("Run XGBoost + SHAP (peds vs adults) on MIC features...")
  xgb_res <- tryCatch(run_xgb_shap(mic_df), error = function(e) { message("XGBoost failed: ", e$message); NULL })
  if (!is.null(xgb_res)) {
    fwrite(xgb_res$shap, file.path(CFG$output_dir, "shap_top.csv"))
    write_lines(as.character(xgb_res$auc), file.path(CFG$output_dir, "xgb_auc.txt"))
  }

  message("Refactored pipeline complete. Outputs in: ", CFG$output_dir)
}

if (interactive() || (!interactive() && identical(Sys.getenv("R_SCRIPT_RUNNING"), "1"))) main()
