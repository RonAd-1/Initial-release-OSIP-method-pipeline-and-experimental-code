
# cardinality_matching.R

# Install
if (!require(causaldata)) install.packages("causaldata")
if (!require(optmatch)) install.packages("optmatch")
if (!require(cobalt)) install.packages("cobalt")
if (!require(sandwich)) install.packages("sandwich")
if (!require(lmtest)) install.packages("lmtest")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")

library(dplyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(cobalt)
library(sandwich)
library(lmtest)
library(optmatch)
library(causaldata)
library(gridExtra)

# Basic usage for NSW data
source("dataset_configs.R")

# ==========================================
# Comparison Metrics Function (Fixed)
# ==========================================

# --- Helper 1: Balance Metrics ---
# comparison_metrics <- function(df, treatment_col, covariates) {
#   
#   # Helper: safe max that never returns -Inf
#   safe_max <- function(x) {
#     x <- x[!is.na(x)]
#     if (length(x) == 0) return(NA_real_)
#     max(x)
#   }
#   
#   # Only keep covariates that actually vary
#   valid_covs <- covariates[
#     vapply(covariates, function(v) {
#       length(unique(df[[v]][!is.na(df[[v]])])) > 1
#     }, logical(1))
#   ]
#   
#   # If nothing varies, return NA metrics immediately
#   if (length(valid_covs) == 0) {
#     return(list(
#       max_smd = NA_real_,
#       max_ks  = NA_real_,
#       tvd     = NA_real_
#     ))
#   }
#   
#   # Compute balance table
#   bal <- suppressWarnings(
#     cobalt::bal.tab(
#       x = df[, valid_covs, drop = FALSE],
#       treat = df[[treatment_col]],
#       estimand = "ATT",
#       stats = c("m", "ks"),
#       un = FALSE
#     )
#   )
#   
#   # Extract adjusted balance statistics safely
#   bal_df <- bal$Balance
#   
#   max_smd <- safe_max(abs(bal_df$Diff.Adj))
#   max_ks  <- safe_max(bal_df$KS.Adj)
#   
#   # Total Variation Distance (only if defined)
#   tvd <- if (!is.null(bal_df$TVD.Adj)) {
#     safe_max(bal_df$TVD.Adj)
#   } else {
#     NA_real_
#   }
#   
#   return(list(
#     max_smd = max_smd,
#     max_ks  = max_ks,
#     tvd     = tvd
#   ))
# }


# --- Helper 2: Statistical Estimation ---
# ==========================================
# 1. FIXED HELPER: calculate_and_format_ate
# ==========================================
# calculate_and_format_ate <- function(matching_res, method_label, data_config) {
#   df <- matching_res$data_matched
#   if (!"weights" %in% names(df)) df$weights <- 1
#   
#   # Only covariates that vary
#   valid_covs <- Filter(
#     function(x) length(unique(df[[x]])) > 1,
#     data_config$NUMERIC_COVARIATES
#   )
#   
#   ate_form <- reformulate(
#     termlabels = c(data_config$TREATMENT_VAR, valid_covs),
#     response   = data_config$OUTCOME_VAR
#   )
#   
#   fit <- tryCatch(
#     lm(ate_form, data = df, weights = weights),
#     error = function(e)
#       lm(
#         reformulate(data_config$TREATMENT_VAR, data_config$OUTCOME_VAR),
#         data = df,
#         weights = weights
#       )
#   )
#   
#   robust_se <- coeftest(fit, vcov = vcovHC(fit, type = "HC2"))
#   row_idx <- grep(data_config$TREATMENT_VAR, rownames(robust_se))[1]
#   
#   est <- robust_se[row_idx, "Estimate"]
#   se  <- robust_se[row_idx, "Std. Error"]
#   
#   # 95% CI
#   ci_low  <- est - 1.96 * se
#   ci_high <- est + 1.96 * se
#   
#   return(data.frame(
#     method    = method_label,
#     n_treated = sum(df[[data_config$TREATMENT_VAR]] == 1),
#     ate       = as.numeric(est),
#     se        = as.numeric(se),
#     ci_low    = as.numeric(ci_low),
#     ci_high   = as.numeric(ci_high),
#     p_value   = as.numeric(robust_se[row_idx, "Pr(>|t|)"]),
#     stringsAsFactors = FALSE
#   ))
# }

# ==========================================
# Updated Cardinality Matching Function
# ==========================================

# run_optimal_1_1_matching <- function(data, data_config) {
#   cat("\n--- Running Optimal (1:1) Matching ---\n")
# 
#   match_formula <- reformulate(termlabels = data_config$NUMERIC_COVARIATES, 
#                                response = data_config$TREATMENT_VAR)
#   dist_matrix <- match_on(match_formula, data = data, method = "mahalanobis")
#   
#   time_start <- Sys.time()
#   m_obj <- pairmatch(dist_matrix, controls = 1, data = data)
#   time_elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
#   
#   # m_obj is an optmatch factor. We convert it to a plain character
#   # to make it "friendly" for your distance calculations later.
#   data$match_id <- as.character(m_obj) 
#   
#   # Filter to only matched units
#   matched_data <- data[!is.na(data$match_id), ]
#   
#   return(list(
#     data_matched = matched_data,
#     time = time_elapsed,
#     m_out = m_obj  # Returning the raw optmatch object just in case
#   ))
# }

# run_optimal_1_1_matching <- function(data, data_config, treatment_col, match_formula, dist_obj = NULL, dist_matrix = NULL) {
#   cat("\n--- Running Optimal (1:1) Matching ---\n")
#   
#   # match_formula <- reformulate(termlabels = data_config$NUMERIC_COVARIATES, 
#   #                              response = treatment_col)
#   # 
#   # # 1. Generate Distance Matrix
#   # # match_on produces an 'optmatch' distance object
#   if (is.null(dist_obj) || is.null(dist_matrix)) {
#     cat("Distance inputs missing. Generating distance object and matrix...\n")
#     dist_obj <- match_on(match_formula, data = data, method = "mahalanobis")
#     dist_matrix <- as.matrix(dist_obj)
#   }
#   
#   # Debug 
#   # browser()
#   
#   # 2. Perform Matching
#   time_start <- Sys.time()
#   m_obj <- pairmatch(dist_obj, controls = 1, data = data)
#   time_elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
#   
#   # 3. Add match_id to data
#   data$match_id <- as.character(m_obj) 
#   matched_data <- data[!is.na(data$match_id), ]
#  
#   match_map_opt <- matched_data %>%
#     group_by(match_id) %>%
#     summarise(
#       # Use which to safely handle indices
#       treated_id = as.character(id[which(get(treatment_col) == 1)]),
#       control_id = as.character(id[which(get(treatment_col) == 0)]),
#       .groups = 'drop'
#     )
#   
#   # Get dimensions for safety
#   rn <- rownames(dist_matrix)
#   cn <- colnames(dist_matrix)
#   
#   # Calculate cost (Mahalanobis distance)
#   match_map_opt$cost <- mapply(function(t, c) {
#     # Case 1: Standard (Treated in rows, Control in cols)
#     if (t %in% rn && c %in% cn) {
#       return(as.numeric(dist_matrix[t, c]))
#     } 
#     # Case 2: Flipped (Control in rows, Treated in cols)
#     else if (c %in% rn && t %in% cn) {
#       return(as.numeric(dist_matrix[c, t]))
#     } 
#     else {
#       # This handles cases where a unit might be missing from the distance matrix
#       return(NA_real_) 
#     }
#   }, match_map_opt$treated_id, match_map_opt$control_id)
#   
#   return(list(
#     data_matched = matched_data,
#     time = time_elapsed,
#     m_out = m_obj,
#     match_map = match_map_opt
#   ))
# }

# ==========================================
# Fixed ATE Estimation Function
# ==========================================

# estimate_ate <- function(matched_data, outcome_var, treatment_col, covariates = NULL) {
#   
#   # FIX 4: Extract outcome as a vector, not data frame
#   outcome_treated <- matched_data[matched_data[[treatment_col]] == 1, outcome_var, drop = TRUE]
#   outcome_control <- matched_data[matched_data[[treatment_col]] == 0, outcome_var, drop = TRUE]
#   
#   # Check if outcome is numeric
#   if (!is.numeric(outcome_treated) || !is.numeric(outcome_control)) {
#     warning(sprintf("Outcome variable '%s' is not numeric. Attempting to convert.", outcome_var))
#     outcome_treated <- as.numeric(outcome_treated)
#     outcome_control <- as.numeric(outcome_control)
#   }
#   
#   # 1. Simple difference in means
#   ate_simple <- mean(outcome_treated, na.rm = TRUE) - mean(outcome_control, na.rm = TRUE)
#   
#   # 2. Regression adjustment
#   if (!is.null(covariates) && length(covariates) > 0) {
#     
#     model_formula <- as.formula(paste(
#       outcome_var, "~", treatment_col, "+", paste(covariates, collapse = " + ")
#     ))
#     
#     model <- lm(model_formula, data = matched_data)
#     ate_adjusted <- coef(model)[treatment_col]
#     
#     vcov_robust <- vcovHC(model, type = "HC2")
#     se_robust <- sqrt(vcov_robust[treatment_col, treatment_col])
#     
#     # Confidence interval
#     ci_lower <- ate_adjusted - 1.96 * se_robust
#     ci_upper <- ate_adjusted + 1.96 * se_robust
#     
#     # P-value
#     t_stat <- ate_adjusted / se_robust
#     p_value <- 2 * pt(abs(t_stat), df = nrow(matched_data) - length(coef(model)), lower.tail = FALSE)
#     
#   } else {
#     # Simple model without covariates
#     model_formula <- as.formula(paste(outcome_var, "~", treatment_col))
#     model <- lm(model_formula, data = matched_data)
#     ate_adjusted <- coef(model)[treatment_col]
#     
#     library(sandwich)
#     vcov_robust <- vcovHC(model, type = "HC2")
#     se_robust <- sqrt(vcov_robust[treatment_col, treatment_col])
#     
#     ci_lower <- ate_adjusted - 1.96 * se_robust
#     ci_upper <- ate_adjusted + 1.96 * se_robust
#     
#     t_stat <- ate_adjusted / se_robust
#     p_value <- 2 * pt(abs(t_stat), df = nrow(matched_data) - 2, lower.tail = FALSE)
#   }
#   
#   results <- list(
#     ate_simple = ate_simple,
#     ate_adjusted = as.numeric(ate_adjusted),  # Ensure numeric
#     se = as.numeric(se_robust),
#     ci_lower = as.numeric(ci_lower),
#     ci_upper = as.numeric(ci_upper),
#     p_value = as.numeric(p_value),
#     n_obs = nrow(matched_data),
#     model = model
#   )
#   
#   return(results)
# }

# ==========================================
# Diagnostic Function (Optional but helpful)
# ==========================================

diagnose_data <- function(data, treatment_col, covariates, outcome_var) {
  
  cat("\n=== Data Diagnostics ===\n")
  
  # Check treatment variable
  cat(sprintf("\nTreatment variable: %s\n", treatment_col))
  cat(sprintf("  Class: %s\n", class(data[[treatment_col]])))
  cat(sprintf("  Unique values: %s\n", paste(unique(data[[treatment_col]]), collapse = ", ")))
  cat(sprintf("  N treated: %d, N control: %d\n", 
              sum(data[[treatment_col]] == 1), 
              sum(data[[treatment_col]] == 0)))
  
  # Check outcome variable
  cat(sprintf("\nOutcome variable: %s\n", outcome_var))
  cat(sprintf("  Class: %s\n", class(data[[outcome_var]])))
  cat(sprintf("  Is numeric: %s\n", is.numeric(data[[outcome_var]])))
  cat(sprintf("  Range: [%.2f, %.2f]\n", 
              min(data[[outcome_var]], na.rm = TRUE),
              max(data[[outcome_var]], na.rm = TRUE)))
  cat(sprintf("  Missing: %d\n", sum(is.na(data[[outcome_var]]))))
  
  # Check covariates
  cat(sprintf("\nCovariates (%d):\n", length(covariates)))
  for (cov in covariates) {
    cat(sprintf("  %s: class=%s, missing=%d\n", 
                cov, 
                class(data[[cov]]),
                sum(is.na(data[[cov]]))))
  }
  
  cat("\n")
}

run_cardinality_matching <- function(data, data_config, treatment_col, match_formula) {
  cat("\n--- Running True Cardinality Matching (IP Solver) ---\n")
  
  # --- STAGE 1: Selection ---
  m.out <- matchit(match_formula,
                   data = data,
                   method = "cardinality", 
                   tols = 0.15,     
                   ratio = 1, 
                   time = 60,
                   solver = "glpk",
                   discard = "none")   
  
  # Get only the units chosen by the solver
  subset_indices <- !is.na(m.out$weights) & m.out$weights > 0
  balanced_subset_data <- data[subset_indices, ]
  
  # Debug
  # browser()
  
  # --- STAGE 2: Pairing ---
  cat("Pairing the balanced subset to resolve distance NAs...\n")
  # This now returns data_matched, m_out, AND match_map
  pairing_res <- run_optimal_1_1_matching(balanced_subset_data, data_config, treatment_col, match_formula)
  
  return(list(
    data_matched = pairing_res$data_matched,
    # This is the map created from the balanced subset
    match_map    = pairing_res$match_map, 
    m_out_stage1 = m.out,        
    m_out_stage2 = pairing_res$m_out  
  ))
}

# run_cardinality_matching <- function(data, data_config) {
#   cat("\n--- Running True Cardinality Matching (IP Solver) ---\n")
#   
#   # --- STAGE 1: Selection ---
#   # Explicitly setting discard = "none" fixes your error
#   m.out <- matchit(reformulate(data_config$NUMERIC_COVARIATES, data_config$TREATMENT_VAR),
#                    data = data,
#                    method = "cardinality", 
#                    tols = 0.1,     
#                    ratio = 1, # For 1:1 matching
#                    solver = "glpk",
#                    discard = "none")   
#   
#   # Get only the units chosen by the solver
#   subset_indices <- !is.na(m.out$weights) & m.out$weights > 0
#   balanced_subset_data <- data[subset_indices, ]
#   
#   # --- STAGE 2: Pairing ---
#   # Use your existing optimal matching logic to create the match_id
#   cat("Pairing the balanced subset to resolve distance NAs...\n")
#   pairing_res <- run_optimal_1_1_matching(balanced_subset_data, data_config)
#   
#   return(list(
#     data_matched = pairing_res$data_matched,
#     m_out_stage1 = m.out,        # The balance-only object
#     m_out_stage2 = pairing_res$m_out  # The PAIRING object (needed for the map)
#   ))
# }
