# multi_carb_short_mat_refactored.R
# Purpose: run univariate tests and simple predictive pipelines (LASSO + Random Forest + LOOCV) for multiple binary outcomes.

library(data.table)
library(tidyverse)
library(readxl)
library(pROC)
library(randomForest)
library(glmnet)
library(caret)
options(scipen = 999)

# --- Configuration ---
CONFIG <- list(
  data_path = "data/db_paulo.xlsx",
  output_dir = "output",
  seed = 123,
  test_frac = 0.3
)

if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

# --- Helpers ---
is_binary <- function(x) {
  ux <- unique(na.omit(x))
  length(ux) == 2
}

detect_vars <- function(df) {
  binary_vars <- df %>% select(where(is_binary)) %>% colnames()
  continuous_vars <- df %>% select(where(~ !is_binary(.))) %>% colnames()
  list(binary = binary_vars, continuous = continuous_vars)
}

run_univariate_tests <- function(df, group_var) {
  vars <- detect_vars(df)
  continuous <- setdiff(vars$continuous, group_var)
  binary <- setdiff(vars$binary, group_var)

  wilcox_results <- map_dfr(continuous, function(var) {
    # Skip variables with all NA or constant
    if (all(is.na(df[[var]])) || length(unique(na.omit(df[[var]]))) < 2) {
      return(tibble(variable = var, p_value = NA_real_, test = "Wilcoxon"))
    }
    safe_test <- tryCatch(wilcox.test(df[[var]] ~ df[[group_var]]), error = function(e) NULL)
    if (is.null(safe_test)) return(tibble(variable = var, p_value = NA_real_, test = "Wilcoxon"))
    tibble(variable = var, p_value = safe_test$p.value, test = "Wilcoxon")
  })

  fisher_results <- map_dfr(binary, function(var) {
    tbl <- table(df[[var]], df[[group_var]])
    if (all(dim(tbl) == c(2, 2))) {
      safe_test <- tryCatch(fisher.test(tbl), error = function(e) NULL)
      p <- if (is.null(safe_test)) NA_real_ else safe_test$p.value
      tibble(variable = var, p_value = p, test = "Fisher")
    } else {
      tibble(variable = var, p_value = NA_real_, test = "Not applicable")
    }
  })

  bind_rows(wilcox_results, fisher_results) %>% arrange(p_value)
}

train_test_split <- function(df, group_var, seed = CONFIG$seed, test_frac = CONFIG$test_frac) {
  set.seed(seed)
  idx <- createDataPartition(df[[group_var]], p = 1 - test_frac, list = FALSE)
  list(train = df[idx, , drop = FALSE], test = df[-idx, , drop = FALSE])
}

fit_lasso <- function(train_df, group_var) {
  x_train <- model.matrix(as.formula(paste(group_var, "~ .")), data = train_df)[, -1, drop = FALSE]
  y_train <- train_df[[group_var]]
  cv <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
  list(model = cv, x_cols = colnames(x_train))
}

predict_lasso_probs <- function(lasso_obj, newdata, group_var) {
  x_new <- model.matrix(as.formula(paste(group_var, "~ .")), data = newdata)[, -1, drop = FALSE]
  probs <- tryCatch(as.vector(predict(lasso_obj$model, newx = x_new, s = lasso_obj$model$lambda.min, type = "response")), error = function(e) rep(NA_real_, nrow(newdata)))
  probs
}

fit_random_forest <- function(train_df, group_var, ntree = 500) {
  randomForest(as.formula(paste(group_var, "~ .")), data = train_df, importance = TRUE, ntree = ntree)
}

loocv_rf_auc <- function(df, group_var, ntree = 500) {
  true_labels <- df[[group_var]]
  probs <- rep(NA_real_, nrow(df))
  for (i in seq_len(nrow(df))) {
    train_df <- df[-i, , drop = FALSE]
    test_df <- df[i, , drop = FALSE]
    rf <- tryCatch(randomForest(as.formula(paste(group_var, "~ .")), data = train_df, ntree = ntree), error = function(e) NULL)
    if (is.null(rf)) { probs[i] <- NA_real_; next }
    p <- tryCatch(predict(rf, test_df, type = "prob")[, 2], error = function(e) NA_real_)
    probs[i] <- p
  }
  auc_val <- tryCatch(auc(roc(as.numeric(true_labels), probs)), error = function(e) NA_real_)
  auc_val
}

plot_feature_importance <- function(rf_model, top_n = 20, title = "Feature Importance") {
  imp <- as.data.frame(importance(rf_model))
  imp$feature <- rownames(imp)
  imp <- imp %>% arrange(desc(MeanDecreaseGini)) %>% head(top_n)
  ggplot(imp, aes(x = reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_col(fill = "#355C7D", alpha = 0.8) +
    coord_flip() +
    labs(title = title, x = "", y = "MeanDecreaseGini") +
    theme_minimal()
}

evaluate_roc <- function(true_vals, probs, title = NULL, outpath = NULL) {
  r <- tryCatch(roc(as.numeric(true_vals), probs), error = function(e) NULL)
  if (is.null(r)) return(list(auc = NA_real_, plot = NULL))
  a <- auc(r)
  p <- ggroc(r) + ggtitle(paste0(title, " (AUC=", round(a, 3), ")")) + theme_minimal()
  if (!is.null(outpath)) ggsave(outpath, p, width = 6, height = 5)
  list(auc = as.numeric(a), plot = p)
}

# --- Main pipeline ---
run_for_outcome <- function(outcome_var, df_raw) {
  message("\n=== Processing outcome: ", outcome_var, " ===\n")
  df <- df_raw
  if (!(outcome_var %in% colnames(df))) stop("Outcome variable not found: ", outcome_var)

  # Ensure outcome is factor/binary
  df[[outcome_var]] <- as.factor(df[[outcome_var]])

  # 1) Univariate tests
  uni <- run_univariate_tests(df, outcome_var)
  write_csv(uni, file.path(CONFIG$output_dir, paste0("univariate_", outcome_var, ".csv")))
  message("Saved univariate results: ", file.path(CONFIG$output_dir, paste0("univariate_", outcome_var, ".csv")))

  # 2) Train/test split
  splits <- train_test_split(df, outcome_var)
  train_df <- splits$train
  test_df <- splits$test

  # Remove non-predictor columns if needed (keep only predictors + outcome)
  # Here we assume all columns except outcome are predictors. Adjust if necessary.

  # 3) LASSO
  lasso_obj <- tryCatch(fit_lasso(train_df, outcome_var), error = function(e) NULL)
  if (!is.null(lasso_obj)) {
    lasso_vars <- rownames(coef(lasso_obj$model, s = lasso_obj$model$lambda.min))[which(coef(lasso_obj$model, s = lasso_obj$model$lambda.min) != 0)]
    lasso_vars <- lasso_vars[lasso_vars != "(Intercept)"]
    write_lines(lasso_vars, file.path(CONFIG$output_dir, paste0("lasso_vars_", outcome_var, ".txt")))
    message("Saved LASSO selected vars: ", file.path(CONFIG$output_dir, paste0("lasso_vars_", outcome_var, ".txt")))
  }

  # 4) Random Forest on training set
  rf_model <- tryCatch(fit_random_forest(train_df, outcome_var), error = function(e) NULL)
  if (!is.null(rf_model)) {
    imp_df <- as.data.frame(importance(rf_model))
    write_csv(imp_df, file.path(CONFIG$output_dir, paste0("rf_importance_", outcome_var, ".csv")))
    p_imp <- plot_feature_importance(rf_model, title = paste("Feature Importance:", outcome_var))
    ggsave(file.path(CONFIG$output_dir, paste0("rf_importance_", outcome_var, ".png")), p_imp, width = 7, height = 6)
    message("Saved RF importance and plot for: ", outcome_var)
  }

  # 5) Evaluate LASSO on test set (ROC)
  if (!is.null(lasso_obj)) {
    probs <- predict_lasso_probs(lasso_obj, test_df, outcome_var)
    roc_res <- evaluate_roc(as.numeric(test_df[[outcome_var]]), probs, title = paste("LASSO ROC -", outcome_var), outpath = file.path(CONFIG$output_dir, paste0("lasso_roc_", outcome_var, ".png")))
    write_lines(as.character(roc_res$auc), file.path(CONFIG$output_dir, paste0("lasso_auc_", outcome_var, ".txt")))
    message("LASSO AUC for ", outcome_var, ": ", roc_res$auc)
  }

  # 6) LOOCV RF AUC (may be slow for big datasets)
  loocv_auc <- tryCatch(loocv_rf_auc(df, outcome_var), error = function(e) NA_real_)
  write_lines(as.character(loocv_auc), file.path(CONFIG$output_dir, paste0("loocv_rf_auc_", outcome_var, ".txt")))
  message("LOOCV RF AUC for ", outcome_var, ": ", loocv_auc)

  invisible(list(univariate = uni, lasso = lasso_obj, rf = rf_model, loocv_auc = loocv_auc))
}

# --- Run for desired outcomes ---
main <- function() {
  message("Loading data: ", CONFIG$data_path)
  df_raw <- read_excel(path = CONFIG$data_path)

  # Define the outcomes to process. Adjust names if your data uses different column names.
  outcomes <- c("group", "In_hospital_mortality", "Mortality_30d")

  results <- map(outcomes, function(outcome) {
    tryCatch(run_for_outcome(outcome, df_raw), error = function(e) { message("Error processing ", outcome, ": ", e$message); NULL })
  })

  names(results) <- outcomes
  saveRDS(results, file = file.path(CONFIG$output_dir, "multi_carb_results.rds"))
  message("All done. Results saved to ", CONFIG$output_dir)
}

# Run when script executed
if (interactive() || (!interactive() && identical(Sys.getenv("R_SCRIPT_RUNNING"), "1"))) {
  main()
}
