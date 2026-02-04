# =====================================================================
# r_kbs_csa_boattini_2023_refactored.R.R
# Clean, modular, production-style refactor of R_KBS_CZA_BoattiniM_2023.R
# Features:
# - Central configuration
# - Reusable helper functions for loading & preprocessing
# - Consolidated univariate logistic regression pipeline
# - Random forest wrapper with basic validation and explainability hooks
# - Plotting helpers and metric calculators
# - Saves CSV/RDS outputs for results
# =====================================================================

library(tidyverse)
library(data.table)
library(broom)
library(randomForest)
library(DALEX)
library(factoextra)
library(cluster)
options(scipen = 999)

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
RAW_FILE <- "KBS_CZA_BoattiniM_2023.txt"
OUTPUT_DIR <- "./"
RF_NTREE_DEFAULT <- 500
RF_MTRY_DEFAULT <- NULL
SEED <- 2024

# Columns to coerce to numeric (based on original script)
NUMERIC_COLS <- c(1, 12, 18, 19, 20, 21, 35, 40)

# Columns to drop (common to many script sections)
DROP_COLS <- c("TIGE_susceptible", "MEV_susceptibility",
               "Charlson_big3", "Charlsonbig4", "Charlsonbig5", "Charlsonbig6",
               "PITT_scorebig2","PITT_scorebig3","PITT_scorebig4","PITT_scorebig5",
               "INCREMENT_big_8","INCREMENT_big_9","INCREMENT_big_10","INCREMENT_big_11")

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

load_data <- function(path = RAW_FILE, colClasses = NULL) {
  message("Loading: ", path)
  df <- fread(path, sep = "\t", colClasses = colClasses)
  df
}

coerce_types <- function(df) {
  # Many blocks in the original script repeatedly coerce many columns
  # We'll follow same mapping: all to factor by default, then numeric for key columns
  df[] <- lapply(df, function(x) x)
  # Convert all character columns to factor (safer for modeling downstream)
  df <- df %>% mutate_if(is.character, as.factor)

  # Convert specific numeric index columns if present
  numeric_indices <- intersect(NUMERIC_COLS, seq_along(df))
  for (i in numeric_indices) {
    df[[i]] <- as.numeric(as.character(df[[i]]))
  }
  df
}

patch_missing_indicators <- function(df) {
  # Apply the common NA -> default transformations used in original
  if ("Prolonged_B_lact_infusion" %in% names(df)) df$Prolonged_B_lact_infusion[is.na(df$Prolonged_B_lact_infusion)] <- 0
  if ("KPC_Infection_relapse" %in% names(df)) df$KPC_Infection_relapse[is.na(df$KPC_Infection_relapse)] <- 0
  if ("Dyalisis_30day_preceding" %in% names(df)) df$Dyalisis_30day_preceding[is.na(df$Dyalisis_30day_preceding)] <- 1
  df
}

clean_df <- function(df) {
  df <- df %>% select(-one_of(intersect(DROP_COLS, names(df))))
  df
}

safe_glm_univariate <- function(df, outcome, predictor) {
  # Returns tidy row (predictor coefficient) or NA-row on error
  f <- as.formula(paste(outcome, "~", "`", predictor, "`", sep = ""))
  res <- tryCatch({
    mod <- glm(f, data = df, family = binomial())
    tidy(mod)[2, ] %>% mutate(feature = predictor)
  }, error = function(e) {
    tibble(term = predictor, estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, feature = predictor)
  })
  res
}

run_univariate_pipeline <- function(df, outcome) {
  predictors <- setdiff(names(df), outcome)
  results <- map_dfr(predictors, ~ safe_glm_univariate(df, outcome, .x))
  # Calculate OR & CI where possible
  results <- results %>%
    mutate(OR = exp(estimate),
           Lower = exp(estimate - std.error),
           Upper = exp(estimate + std.error)) %>%
    select(feature, estimate, std.error, statistic, p.value, OR, Lower, Upper)
  results
}

plot_univariate_or <- function(storage, title) {
  storage <- storage %>% drop_na() %>% filter(is.finite(Upper))
  storage %>%
    mutate(p_cat = ifelse(p.value < 0.05, "<0.05", "ns"),
           p_cat = factor(p_cat, levels = c("<0.05","ns"))) %>%
    ggplot(aes(y = fct_reorder(feature, OR), x = OR, xmin = Lower, xmax = Upper, colour = p_cat)) +
    geom_errorbarh(height = .05, size = 1, alpha = 0.8) +
    geom_point(size = 2) +
    labs(title = title, x = "Odds Ratio (log scale)", y = "Predictor") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    theme_minimal() +
    scale_colour_manual(values = c("<0.05" = "lightpink3", "ns" = "lightskyblue4"))
}

run_random_forest <- function(df, outcome, ntree = RF_NTREE_DEFAULT, mtry = RF_MTRY_DEFAULT, seed = SEED) {
  set.seed(seed)
  predictors <- setdiff(names(df), outcome)
  formula <- as.formula(paste(outcome, "~ ."))
  rf <- randomForest(formula, data = df, ntree = ntree, mtry = mtry)
  importance <- as.data.frame(rf$importance) %>% rownames_to_column(var = "feature")
  list(model = rf, importance = importance)
}

rf_validation_metrics <- function(model, train_df, test_df, outcome) {
  preds <- predict(model, test_df, type = "response")
  tab <- table(Actual = test_df[[outcome]], Predicted = preds)
  accuracy <- sum(diag(tab))/sum(tab)
  precision <- ifelse((tab[2,2] + tab[1,2]) == 0, NA, tab[2,2]/(tab[2,2] + tab[1,2]))
  recall <- ifelse((tab[2,2] + tab[2,1]) == 0, NA, tab[2,2]/(tab[2,2] + tab[2,1]))
  f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA, 2 * (precision * recall) / (precision + recall))
  list(confusion = tab, accuracy = accuracy, precision = precision, recall = recall, f1 = f1)
}

# ---------------------------------------------------------------------
# Main analysis flow
# ---------------------------------------------------------------------

main <- function() {
  # Load and basic cleaning
  df <- load_data(RAW_FILE)
  df <- patch_missing_indicators(df)
  df <- coerce_types(df)
  df <- clean_df(df)

  # Outcome variables used in the original: D30_mortality, In_hospital_death
  outcomes <- intersect(c("D30_mortality", "In_hospital_death"), names(df))

  # Run univariate pipelines for each outcome
  for (outcome in outcomes) {
    df_local <- df %>% drop_na(.data[[outcome]])
    # Ensure outcome is factor
    df_local[[outcome]] <- as.factor(df_local[[outcome]])
    res <- run_univariate_pipeline(df_local, outcome)
    out_csv <- file.path(OUTPUT_DIR, paste0("univariate_", outcome, ".csv"))
    write_csv(res, out_csv)
    message("Saved univariate results for ", outcome, " -> ", out_csv)

    # Plot ORs (save when possible)
    p <- plot_univariate_or(res, paste0("Univariate ORs: ", outcome))
    ggsave(filename = file.path(OUTPUT_DIR, paste0("univariate_ORs_", outcome, ".png")), plot = p, width = 8, height = 10)
  }

  # Random forest on In_hospital_death (as in original)
  if ("In_hospital_death" %in% names(df)) {
    df_rf <- df %>% drop_na(In_hospital_death)
    # split small test set using stratified sampling
    set.seed(SEED)
    train_idx <- df_rf %>% group_by(In_hospital_death) %>% sample_frac(0.8) %>% ungroup() %>% pull(row_number())
    # fallback to simple sample if stratified fails
    if (length(train_idx) == 0) train_idx <- sample(nrow(df_rf), floor(0.8 * nrow(df_rf)))
    train_df <- df_rf[train_idx, ]
    test_df <- df_rf[-train_idx, ]

    rf_res <- run_random_forest(train_df, "In_hospital_death")
    saveRDS(rf_res$model, file = file.path(OUTPUT_DIR, "rf_in_hospital_model.rds"))
    write_csv(rf_res$importance, file.path(OUTPUT_DIR, "rf_in_hospital_importance.csv"))

    # Evaluate on test set
    metrics <- rf_validation_metrics(rf_res$model, train_df, test_df, "In_hospital_death")
    saveRDS(metrics, file = file.path(OUTPUT_DIR, "rf_in_hospital_metrics.rds"))
    message("RF trained and metrics saved for In_hospital_death")

    # Provide DALEX explainer (optional; requires DALEX available)
    try({
      explainer <- explain(rf_res$model, data = train_df %>% select(-In_hospital_death), y = train_df$In_hospital_death, label = "rf_in_hospital")
      saveRDS(explainer, file = file.path(OUTPUT_DIR, "rf_in_hospital_explainer.rds"))
    }, silent = TRUE)
  }

  # Random forest for D30_mortality (if present)
  if ("D30_mortality" %in% names(df)) {
    df_rf <- df %>% drop_na(D30_mortality)
    set.seed(SEED)
    train_idx <- df_rf %>% group_by(D30_mortality) %>% sample_frac(0.8) %>% ungroup() %>% pull(row_number())
    if (length(train_idx) == 0) train_idx <- sample(nrow(df_rf), floor(0.8 * nrow(df_rf)))
    train_df <- df_rf[train_idx, ]
    test_df <- df_rf[-train_idx, ]
    rf_res <- run_random_forest(train_df, "D30_mortality")
    saveRDS(rf_res$model, file = file.path(OUTPUT_DIR, "rf_d30_model.rds"))
    write_csv(rf_res$importance, file.path(OUTPUT_DIR, "rf_d30_importance.csv"))
    metrics <- rf_validation_metrics(rf_res$model, train_df, test_df, "D30_mortality")
    saveRDS(metrics, file = file.path(OUTPUT_DIR, "rf_d30_metrics.rds"))
    message("RF trained and metrics saved for D30_mortality")
  }

  message("Refactored pipeline complete. Outputs saved to: ", normalizePath(OUTPUT_DIR))
}

# Run main only when interactive
if (interactive()) main()

# End of file
