# =====================================================================
# kbs_covid_boattini_2023_refactored.R
# =====================================================================

# Libraries
library(tidyverse)
library(data.table)
library(broom)
library(mice)
library(Boruta)
library(randomForest)
library(gbm)
library(glmnet)

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
DATA_RAW <- "KBS_COVID_Boattini23_v2.txt"
IMPUTED_FILE <- "Imputed_data.txt"
OUTPUT_IMPUTED <- "Imputed_data.txt"
REGSUBS_FILE <- "RegSubs.csv"
LASSO_FILE <- "Lasso_Plot.csv"

NUM_MICE <- 50
MAXIT_MICE <- 40
SEED <- 2024

# Numeric columns to coerce (common set used in original script)
NUMERIC_COLS <- c(
  "Age", "Charlson_comorbidity_index", "D_dimer", "LDH", "CPK",
  "NT_proBNP", "Troponin_T", "Ferritin", "Creatinine",
  "Lymphocytes_count", "Procalcitonin", "CRP", "LoS", "ICU_LoS"
)

# Columns to patch values for
PATCH_RULES <- list(
  CN_MV_MR_CPAP = function(x) ifelse(x == "2", "0", x),
  Obesity = function(x) ifelse(x == "", "0", x),
  Pulmonary_embolism = function(x) ifelse(x == "11", "1", x)
)

# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

load_raw <- function(path = DATA_RAW) {
  fread(path, sep = "\t", colClasses = "character")
}

save_imputed <- function(df, path = OUTPUT_IMPUTED) {
  fwrite(df, path, sep = "\t")
}

apply_patch_rules <- function(df) {
  for (col in names(PATCH_RULES)) {
    if (col %in% names(df)) {
      df[[col]] <- PATCH_RULES[[col]](df[[col]])
    }
  }
  df
}

coerce_numeric_cols <- function(df, cols = NUMERIC_COLS) {
  for (col in intersect(cols, names(df))) {
    df[[col]] <- as.numeric(df[[col]])
  }
  df
}

convert_factors <- function(df) {
  df %>% mutate_if(is.character, as.factor)
}

# Run MICE imputation and save result
run_imputation <- function(df, num = NUM_MICE, maxit = MAXIT_MICE, seed = SEED) {
  set.seed(seed)
  # Work on a copy to avoid side effects
  temp <- df
  # Convert obvious numerics to numeric before imputation
  temp <- coerce_numeric_cols(temp)
  imp <- mice(temp, m = num, maxit = maxit, printFlag = FALSE)
  completed <- complete(imp)
  list(imputed = completed, mice_obj = imp)
}

# Summarize numeric means by a grouping variable
summarize_by_group <- function(df, group_var, numeric_vars) {
  df %>% group_by(.data[[group_var]]) %>%
    summarise(across(all_of(numeric_vars), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup()
}

# Run univariate logistic GLMs for each predictor against outcome 'Death'
run_univariate_glms <- function(df, outcome = "Death") {
  predictors <- setdiff(names(df), outcome)
  results <- map_df(predictors, function(pred) {
    # Build formula dynamically
    f <- as.formula(paste(outcome, "~", "`", pred, "`", sep = ""))
    # try-catch to avoid stopping on problematic predictors
    res <- tryCatch({
      model <- glm(f, data = df, family = binomial())
      tidy_res <- tidy(model)
      # take coefficient for predictor (usually second row)
      coef_row <- tidy_res[2, ]
      data.frame(
        feature = pred,
        estimate = coef_row$estimate,
        std.error = coef_row$std.error,
        statistic = coef_row$statistic,
        p.value = coef_row$p.value,
        OR = exp(coef_row$estimate),
        Lower = exp(coef_row$estimate - coef_row$std.error),
        Upper = exp(coef_row$estimate + coef_row$std.error),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(feature = pred, estimate = NA, std.error = NA, statistic = NA, p.value = NA, OR = NA, Lower = NA, Upper = NA)
    })
    res
  })
  results %>% mutate_if(is.numeric, ~ round(., 5))
}

# Run a multivariable logistic regression with an explicit formula or a vector of predictors
run_multivariable_glm <- function(df, predictors, outcome = "Death") {
  f <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  mod <- glm(f, data = df, family = binomial())
  res <- tidy(mod)
  res$OR <- exp(res$estimate)
  res$Lower <- exp(res$estimate - res$std.error)
  res$Upper <- exp(res$estimate + res$std.error)
  res %>% mutate_if(is.numeric, ~ round(., 5))
}

# Run Boruta feature selection (wrap so long prints don't spam)
run_boruta <- function(df, outcome = "Death") {
  set.seed(SEED)
  # Boruta requires factors/numerics - convert characters already done upstream
  Boruta(as.formula(paste(outcome, "~ .")), data = df, doTrace = 0)
}

# Run basic classifiers: Random Forest and GBM (only when packages available)
run_classifiers <- function(df, outcome = "Death") {
  df_factors <- df %>% select(where(is.factor))
  res <- list()
  try({
    rf <- randomForest(as.formula(paste(outcome, "~ .")), data = df_factors, importance = TRUE)
    res$random_forest <- rf
    res$rf_importance <- data.frame(rf$importance) %>% rownames_to_column("feature")
  }, silent = TRUE)

  try({
    gbm_mod <- gbm(as.formula(paste(outcome, "== 1 ~ .")), data = df_factors, n.trees = 5000, distribution = "bernoulli")
    res$gbm <- gbm_mod
  }, silent = TRUE)
  res
}

# LASSO pipeline (glmnet) - expects x matrix and y vector
run_lasso <- function(df, outcome = "Death") {
  x <- model.matrix(as.formula(paste(outcome, "~ .")), df)[, -1]
  y <- as.numeric(df[[outcome]])
  cv <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure = "auc")
  fit <- glmnet(x, y, alpha = 1, lambda = cv$lambda.1se, family = "binomial")
  list(cv = cv, fit = fit)
}

# ---------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------

main <- function() {
  # 1. Load raw data
  raw <- load_raw()

  # 2. Apply patch rules and coerce numerics
  raw <- apply_patch_rules(raw)
  raw <- coerce_numeric_cols(raw)
  raw <- convert_factors(raw)

  # 3. Missingness diagnostics (optional) - user may skip heavy visualizations
  cat("Missing values per column:\n")
  print(data.frame(colSums(is.na(raw))))

  # 4. Run imputation (if not already saved)
  if (!file.exists(IMPUTED_FILE)) {
    imp_res <- run_imputation(raw)
    imputed <- imp_res$imputed
    save_imputed(imputed, IMPUTED_FILE)
  } else {
    imputed <- fread(IMPUTED_FILE, sep = "\t", colClasses = "character")
    imputed <- apply_patch_rules(imputed)
    imputed <- coerce_numeric_cols(imputed)
    imputed <- convert_factors(imputed)
  }

  # 5. Univariate logistic regressions
  univar <- run_univariate_glms(imputed, outcome = "Death")
  write.csv(univar, "univariate_results.csv", row.names = FALSE)
  cat("Univariate regression results saved to univariate_results.csv\n")

  # 6. Identify candidate predictors for multivariable model (simple approach: p<0.1)
  candidates <- univar %>% filter(!is.na(p.value) & p.value < 0.1) %>% pull(feature)
  cat("Candidate predictors for multivariable model (p<0.1):\n")
  print(candidates)

  if (length(candidates) > 0) {
    multiv_res <- run_multivariable_glm(imputed, candidates, outcome = "Death")
    write.csv(multiv_res, "multivariable_results.csv", row.names = FALSE)
    cat("Multivariable regression results saved to multivariable_results.csv\n")
  }

  # 7. Boruta feature selection (optional; can be slow)
  bor <- tryCatch({
    run_boruta(imputed, outcome = "Death")
  }, error = function(e) {
    message("Boruta failed: ", e$message)
    NULL
  })
  if (!is.null(bor)) {
    saveRDS(bor, "boruta_results.rds")
    cat("Boruta results saved to boruta_results.rds\n")
  }

  # 8. Classifiers (RF / GBM)
  classifiers <- run_classifiers(imputed, outcome = "Death")
  if (!is.null(classifiers$rf_importance)) {
    write.csv(classifiers$rf_importance, "rf_importance.csv", row.names = FALSE)
    cat("Random forest importance saved to rf_importance.csv\n")
  }

  # 9. LASSO pipeline
  lasso_res <- tryCatch({
    run_lasso(imputed, outcome = "Death")
  }, error = function(e) {
    message("LASSO failed: ", e$message)
    NULL
  })
  if (!is.null(lasso_res)) {
    saveRDS(lasso_res, "lasso_results.rds")
    cat("LASSO results saved to lasso_results.rds\n")
  }

  # 10. Summary statistics (means / counts)
  numeric_vars <- intersect(NUMERIC_COLS, names(imputed))
  if (length(numeric_vars) > 0) {
    summary_by_death <- summarize_by_group(imputed, "Death", numeric_vars)
    write.csv(summary_by_death, "summary_by_death.csv", row.names = FALSE)
  }

  cat("Refactored pipeline completed. Results saved to working directory.\n")
}

# Run main when script is executed
if (identical(Sys.getenv("RSTUDIO"), "1") || interactive()) {
  main()
}

# End of file
