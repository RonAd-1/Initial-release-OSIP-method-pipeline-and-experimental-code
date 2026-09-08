
# main_comparison_osip.R

library(viridis)
library(Matching)
library(rgenoud)
library(sensitivityfull)
library(sensitivitymv)
library(MASS)
library(PSAgraphics)
library(patchwork)
library(htmlwidgets)
library(webshot2)
library(MatchIt)
library(optmatch)
library(mediation)
library(causaldata)
library(ggplot2)
library(ggridges)
library(tidyr)
library(cobalt)
library(scales)
library(plotly) 
library(dplyr)
library(purrr)
library(sensitivitymult)

source("dataset_configs.R", echo = FALSE)
source("osip_functions.R", echo = FALSE)
source("plots.R", echo = FALSE)
source("complete_comparison.R", echo = FALSE)
source("original_distance_for_intervals.R", echo = FALSE)
source("heuristic_units.R", echo = FALSE)
source("borders_helping_functions.R", echo = FALSE)

# Load your choice from dataset_configs.R file
datasets <- DATASET_CHOICES
params <- PARAMS
delta_values <- params$DELTA_VALUES
methods <- Methods
method_labels <- METHOD_LABELS
matching_non_1_1 <- METHODS_NON_1_TO_1_MATCHING

# dataset_name <- DATASET_CHOICES$RHC
# dataset_name <- DATASET_CHOICES$NSW_MIXTAPE
dataset_name <- DATASET_CHOICES$LINDNER
# dataset_name <- DATASET_CHOICES$JOBS
# dataset_name <- DATASET_CHOICES$IDHP
# dataset_name <- DATASET_CHOICES$NHEFS

# Container for results
delta_sweep_results <- list()

dist_power <- 1

target_dir <- file.path(
  "outputs",
  dataset_name,
  "osip"
)

base_dir <- tryCatch({
  # recursive = TRUE ensures both 'outputs' and the subfolder are created
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  # Add a trailing slash for your sprintf later
  paste0(target_dir, "/") 
  
}, error = function(e) {
  message("Warning: Couldn't create folder structure. Saving in current folder.")
  message(paste("Reason:", e$message))
  "" 
})

debug_glb = params$DEBUG_STATUS[1]

sample_flag = params$SAMPLE[2] # 1 for TRUE, 2 for FALSE 

max_val_for_plot <- params$MAX_VAL_FOR_PLOT

prepared_data <- load_and_prep_data(dataset_name, sample_flag, datasets, params)

# Extract objects for the workspace
dataset_full <- prepared_data$data
data_config  <- prepared_data$data_config
treatment_col <- data_config$TREATMENT_VAR
outcome_var <- data_config$OUTCOME_VAR
id_var <- data_config$ID_VAR

# [IMPROVED] GLOBAL FLIP FOR LINDNER
if (dataset_name == datasets$LINDNER || dataset_name == datasets$JOBS) {
  cat("\n[!] JOBS or LINDNER detected: Creating 'treat_flipped' for matching...\n")
  
  # Ensure the original is numeric first
  dataset_full[[treatment_col]] <- ensure_numeric_col(dataset_full, treatment_col)
  
  # Create the new column (Original 0 becomes 1, Original 1 becomes 0)
  dataset_full$treat_flipped <- 1 - dataset_full[[treatment_col]]
  
  # Update your treatment variable pointer for the analysis functions
  # This tells your matching functions to look at the flipped column instead
  treatment_col <- "treat_flipped"
} 

if (sample_flag) {
  data_subset <- sample_dataset (seed_in, treatment_col, dataset_full, dataset_name, params)
} else {data_subset <- dataset_full}

# 2. Step B: PS on FULL dataset + Plot Full Distribution
message("📊 Calculating PS and plotting Full Population...")
data_subset <- process_propensity_scores(data_subset, dataset_name, data_config, treatment_col, datasets, base_dir, title_suffix = "(Full Population)")

p_sorted <- sort(unique(c(0, 1, data_subset$ps)))

message("📊 Plotting Sampled Distribution...")
plot_ps_distribution(data_subset, dataset_name, treatment_col, base_dir, title_suffix = "(Full Population)", is_trimmed = FALSE)

# ============================================================
# ✅ COMMON SUPPORT TRIMMING (MODERNIZED)
# ============================================================
cat("\n\n=== Common Support Analysis ===\n")

# Use your dynamic treatment column from data_config
t_vec <- data_subset[[treatment_col]]
ps_vec <- data_subset$ps

# Original sample sizes
n_treat_original <- sum(t_vec == 1)
n_control_original <- sum(t_vec == 0)

message(sprintf("Original sample: %d treated, %d controls\n", 
                n_treat_original, n_control_original))

# Calculate PS ranges by group (using the dynamic treatment vector)
ps_range_treated <- range(ps_vec[t_vec == 1], na.rm = TRUE)
ps_range_control <- range(ps_vec[t_vec == 0], na.rm = TRUE)

cat(sprintf("PS range - Treated: [%.4f, %.4f]\n", ps_range_treated[1], ps_range_treated[2]))
cat(sprintf("PS range - Control: [%.4f, %.4f]\n", ps_range_control[1], ps_range_control[2]))

# Define common support region (The Overlap)
ps_min <- max(ps_range_treated[1], ps_range_control[1])
ps_max <- min(ps_range_treated[2], ps_range_control[2])

cat(sprintf("Common support region: [%.4f, %.4f]\n", ps_min, ps_max))

# --- CONTROL PARAMETER ---
# We recommend keeping this TRUE for the osip method to avoid DP "traps"
TRIM_TO_COMMON_SUPPORT <- TRUE 

if (TRIM_TO_COMMON_SUPPORT) {
  # 1. Identify the Overlap
  ps_min <- max(min(ps_vec[t_vec == 1]), min(ps_vec[t_vec == 0]))
  ps_max <- min(max(ps_vec[t_vec == 1]), max(ps_vec[t_vec == 0]))
  
  # 2. Filter the dataframe
  data_subset <- data_subset[data_subset$ps >= ps_min & data_subset$ps <= ps_max, ]
  
  # === ADDED: Print post-trimming unit counts ===
  # Assumes treatment_col is a character string holding the column name (e.g., "abcix")
  n_treated <- sum(data_subset[[treatment_col]] == 1)
  n_control <- sum(data_subset[[treatment_col]] == 0)
  
  message(sprintf("Post-Trimming Counts: Treated = %d, Control = %d (Total N = %d)", 
                  n_treated, n_control, nrow(data_subset)))
  # ==============================================
  
  # 3. CRITICAL: Update p_sorted to the ACTUAL unit boundaries
  # We no longer anchor to 0 and 1 here to keep the DP focused
  p_sorted <- sort(unique(data_subset$ps))
  
  n_p <- length(p_sorted)
  left_border <- p_sorted[1]
  right_border <- p_sorted[n_p]
  
  message(sprintf("DP Space recalibrated to: [%.4f, %.4f]", left_border, right_border))
  
  plot_ps_distribution(data_subset, dataset_name, treatment_col, base_dir, title_suffix = "(Common Support Trimmed)", is_trimmed = TRUE)
  
} else {
  cat("\nNo trimming applied (TRIM_TO_COMMON_SUPPORT = FALSE)\n")
}
cat("==========================================\n\n")
# ============================================================

dataset_post_trimming <- data_subset

treatment_map <- get_treatment_map(data_subset, treatment_col, id_var)

treatment_control_scores <- get_treatment_control_scores(treatment_map, data_subset, id_var, treatment_col)

treatment_scores <- treatment_control_scores$treatment_scores
control_scores <- treatment_control_scores$control_scores

cov_df_standardized <- prepare_standardized_data(data_subset, data_config, treatment_col)

# FIX: Convert to data.frame first (tibbles don't support rownames), then set IDs as rownames
cov_df_standardized <- as.data.frame(cov_df_standardized)
rownames(cov_df_standardized) <- as.character(cov_df_standardized[[data_config$ID_VAR]])

match_formula <- reformulate(data_config$ALL_COVARIATES, treatment_col)
dist_obj <- match_on(match_formula, data = cov_df_standardized, method = "mahalanobis")
dist_matrix <- as.matrix(dist_obj)

# DEBUG: verify alignment
# Are treated unit IDs actually 1-181, or are they something else?
treated_units <- cov_df_standardized[cov_df_standardized[[treatment_col]] == 1, ]
cat("Treated IDs:", head(treated_units[[data_config$ID_VAR]], 10), "\n")
cat("Treated row numbers:", head(as.integer(rownames(treated_units)), 10), "\n")

# Are they the same thing?
cat("IDs match row numbers:", 
    all(as.character(treated_units[[data_config$ID_VAR]]) == rownames(treated_units)), "\n")

# Create a numeric matrix of ALL covariates (handles factors automatically)
# The "-1" removes the intercept column
X_all <- model.matrix(
  as.formula(paste("~", paste(data_config$ALL_COVARIATES, collapse = " + "), "-1")), 
  data = cov_df_standardized
)

S_inv <- solve(cov(X_all) + diag(1e-7, ncol(X_all)))

if (is.null(S_inv)) stop(paste(metric, "matrix inversion failed."))

params$S_inv <- S_inv

# --- SETUP PHASE (Do this once per dataset) ---
# 1. Standard Inputs
X_std <- X_all 
S_inv_std <- S_inv

# 2. Robust Inputs
robust_env <- get_robust_inputs(X_all)
X_rob <- robust_env$X_all
S_inv_rob <- robust_env$S_inv

# 1. Start the timer
start_time = start_timer() 

# Todo: Return later, currently analysis is on single delta value
for (current_delta in delta_values) {

  cat(sprintf("\n\n>>> STARTING ANALYSIS FOR DELTA = %f <<<\n", current_delta))
  
  output_path = sprintf("%smain_%s_output_delta_%g_power_%g.txt", base_dir,
                        dataset_name, current_delta, dist_power)

  # output_path = sprintf("%smain_%s_output.txt", base_dir, dataset_name)

  # Start capturing output to a file
  sink(output_path, type = "output")
  
  delta_path <- sprintf("delta_%g", current_delta)
  
  power <- sprintf("power_%g", dist_power)
  
  base_dir_delta <- file.path(base_dir, delta_path, power)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(base_dir_delta)) {
    dir.create(base_dir_delta, recursive = TRUE)
  }
  
  # 1. Create a fresh copy of global params
  iter_params <- params
  
  if (TRIM_TO_COMMON_SUPPORT) 
    k_bound <- get_max_k(current_delta, left_border, right_border)
  # Use 0 and 1 as the left_border, right_border
  else
    k_bound <- get_max_k(current_delta, 0, 1)
  
  cat(sprintf("[Params Update] DELTA_DP: %.3f | K_DP: %d\n", 
              current_delta, k_bound))
  
  
  # 4. Now proceed to osip's Step 1 with borders
  osip_res_step_1 <- run_osip_step1(cov_df_standardized, iter_params, 
                                    id_var, p_sorted, treatment_col,
                                    treatment_scores, control_scores, 
                                    k_bound, current_delta, 
                                    left_border, right_border, dist_power)
  
  # --- FEASIBILITY CHECK ---
  # A solution is valid ONLY if there is no error AND we have an optimal cost
  is_failed <- !is.null(osip_res_step_1$error) || 
    is.null(osip_res_step_1$cost_dist_optimal) ||
    is.infinite(osip_res_step_1$cost_dist_optimal)
  
  if (is_failed) {
    # 1. Extract the specific message for the user
    error_info <- if (!is.null(osip_res_step_1$error)) {
      if (is.list(osip_res_step_1$error)) osip_res_step_1$error$error else as.character(osip_res_step_1$error)
    }
    # 2. Log the failure
    cat(sprintf("\n[!] SKIPPING Delta = %.2f: %s\n", current_delta, error_info))
    
    # 3. Store a failure record (optional, but helpful for the final summary)
    delta_sweep_results[[as.character(current_delta)]] <- list(delta = current_delta, status = "Infeasible")
    
    # 4. THE JUMP: Skip the rest of this loop iteration
    next 
  }
  
  # Extract DP results
  dp_intervals <- osip_res_step_1$intervals_dist_optimal
  dp_cost <- osip_res_step_1$cost_dist_optimal
  dp_k <- osip_res_step_1$k_dist_optimal
  
  cat("\n--- DEBUG ---\n")
  cat("Delta value:", current_delta, "\n")
  cat("treatment_col:", treatment_col, "\n")
  cat("Class:", class(cov_df_standardized[[treatment_col]]), "\n")
  print(table(cov_df_standardized[[treatment_col]], useNA = "ifany"))
  cat("Rows in cov_df:", nrow(cov_df_standardized), "\n")
  cat("--- END DEBUG ---\n")
  
  
  # Storing the csv file of step 1
  save_partition_csv(
    df        = dp_intervals,
    delta_val = current_delta,
    base_dir  = base_dir_delta,
    step      = 1
  )

  osip_res_strict_step1 <- run_osip_pipeline(
    metric_label = "OSIP-Strict",
    robust_flag = FALSE,
    dp_intervals = dp_intervals,
    cov_df = cov_df_standardized,
    X_working = X_std,
    S_inv_working = S_inv_std,
    dist_matrix = dist_matrix, # Standard Mahalanobis matrix
    iter_params = iter_params,
    data_config = data_config,
    treatment_col = treatment_col,
    outcome_var = outcome_var,
    base_dir = base_dir_delta, 
    delta_dp = current_delta,
    k_bound  = k_bound, 
    left_border = left_border, 
    right_border = right_border,
    max_val_for_plot = max_val_for_plot,
    matching_non_1_1 = matching_non_1_1,
    dataset_name = dataset_name
    # run_step2 = FALSE,
    # osip_res_step1 = NULL
  )
  
  # --- RUN STRICT VERSION ---
  osip_res_strict_step2 <- run_osip_pipeline(
    metric_label = "OSIP-Strict",
    robust_flag = FALSE,
    dp_intervals = dp_intervals,
    cov_df = cov_df_standardized,
    X_working = X_std,
    S_inv_working = S_inv_std,
    dist_matrix = dist_matrix, # Standard Mahalanobis matrix
    iter_params = iter_params,
    data_config = data_config,
    treatment_col = treatment_col,
    outcome_var = outcome_var,
    base_dir = base_dir_delta, 
    delta_dp = current_delta,
    k_bound  = k_bound,
    left_border = left_border, 
    right_border = right_border,
    max_val_for_plot = max_val_for_plot,
    matching_non_1_1 = matching_non_1_1,
    dataset_name = dataset_name,
    # analyze_flag = TRUE,
    # analyze_flag = FALSE,
    run_step2 = TRUE
    # osip_res_step1 = osip_res_strict_step1
  )
  
  # --- RUN ROBUST VERSION ---
  # Note: Generate the robust distance matrix if not already pre-calculated
  dist_matrix_robust <- get_distance_matrix(cov_df_standardized, data_config, treatment_col, robust = TRUE)
  
  osip_res_robust_step1 <- run_osip_pipeline(
    metric_label = "OSIP-Robust",
    robust_flag = TRUE,
    dp_intervals = dp_intervals,
    cov_df = cov_df_standardized,
    X_working = X_rob,
    S_inv_working = S_inv_rob,
    dist_matrix = dist_matrix_robust, # Robust Mahalanobis matrix
    iter_params = iter_params,
    data_config = data_config,
    treatment_col = treatment_col,
    outcome_var = outcome_var,
    base_dir = base_dir_delta, 
    delta_dp = current_delta,
    k_bound  = k_bound,
    left_border = left_border, 
    right_border = right_border,
    max_val_for_plot = max_val_for_plot,
    matching_non_1_1 = matching_non_1_1,
    dataset_name = dataset_name
    # analyze_flag = FALSE,
    # run_step2 = FALSE,
    # osip_res_step1 = NULL
  )
  
  osip_res_robust_step2 <- run_osip_pipeline(
    metric_label = "OSIP-Robust",
    robust_flag = TRUE,
    dp_intervals = dp_intervals,
    cov_df = cov_df_standardized,
    X_working = X_rob,
    S_inv_working = S_inv_rob,
    dist_matrix = dist_matrix_robust,
    iter_params = iter_params,
    data_config = data_config,
    treatment_col = treatment_col,
    outcome_var = outcome_var,
    base_dir = base_dir_delta,
    delta_dp = current_delta,
    k_bound  = k_bound,
    left_border = left_border, 
    right_border = right_border,
    max_val_for_plot = max_val_for_plot,
    matching_non_1_1 = matching_non_1_1,
    dataset_name = dataset_name,
    run_step2 = TRUE
    # analyze_flag = TRUE,
    # osip_res_step1 = osip_res_robust_step1
    # analyze_flag = FALSE
  )
  
  print("\n\n*******************************************\n\n")
  print("Second comparison: complete comparison and sensitivity analysis for OSIP variants\n")
  
  
  inner_folder_name_comparison <- sprintf("files_delta_%.3f_comparison", current_delta)
  full_inner_path_comparison <- file.path(base_dir_delta, inner_folder_name_comparison)
  # --- THE FIX: Create the folder if it doesn't exist ---
  
  if (!dir.exists(full_inner_path_comparison)) {
    dir.create(full_inner_path_comparison, recursive = TRUE)
  }
  
  final_stats_comparison <- compare_methods(
    data = cov_df_standardized,
    matching_non_1_1 = matching_non_1_1,
    data_config = data_config,
    treatment_col = treatment_col,
    
    methods = methods,
    method_labels = method_labels,
  
    osip_step1_strict_res = NULL,
    osip_step1_robust_res = NULL,
    
    # Strict solutions step 2
    # osip_step2_strict_best_balanced = osip_step2_strict_best_balanced,
    osip_step2_strict_res = osip_res_strict_step2, 
    
    # Robust solutions step 2
    # osip_step2_robust_best_balanced = osip_step2_robust_best_balanced,
    osip_step2_robust_res = osip_res_robust_step2,
  )
  
  ate_comparison_plot_comparison <- plot_ate_comparison(
    final_stats_comparison, 
    data_config, 
    current_delta, 
    dataset_name, 
    full_inner_path_comparison
  )
  
  print(ate_comparison_plot_comparison)
  
  # ==========================================
  # RUN THE FULL COMPARISON
  # ==========================================
  
  results_comparison <- complete_comparison(
    results_df     = final_stats_comparison, # The table containing ATE/P-values
    data           = cov_df_standardized,  # The original baseline
    data_config    = data_config,          # Contains ALL_COVARIATES and ID_VAR
    treatment_col  = treatment_col,        # "treat"
    outcome_var    = data_config$OUTCOME_VAR,
    
    # Pass the full result objects, NOT just the matched data
    osip_step1_strict_res = NULL,
    osip_step1_robust_res = NULL,
    
    # Strict solutions step 2
    # osip_step2_strict_best_balanced = osip_step2_strict_best_balanced,
    osip_step2_strict_res = osip_res_strict_step2, 
    
    # Robust solutions step 2
    # osip_step2_robust_best_balanced = osip_step2_robust_best_balanced,
    osip_step2_robust_res = osip_res_robust_step2
  )
  
  print("\n\n*******************************************\n\n")
  print("Results after complete_comparison and sensitivity analysis:\n")
  print(results_comparison)
  
  # 1. Start with the ATE results (The "Anchor" table)
  final_results <- final_stats_comparison %>%
    # 2. Use left_join to keep CEM/Full Matching even if they lack distance data
    left_join(results_comparison$distances, by = "method") %>%
    # 3. Join with Balance Metrics (Max_SMD, Max_KS, etc.)
    left_join(results_comparison$metrics, by = "method")
  
  print(final_results)
  
  # Define the filename based on your current run parameters
  output_name <- sprintf("consolidated_results_%s_delta_%.2f_%s.csv", dataset_name, current_delta, power)
  write.csv(final_results, file.path(base_dir_delta, output_name), row.names = FALSE)
  
  # Pass the specific sub-list of objects to the Love Plot
  love_plot_comparison <- create_love_plot(
    bal_list = results_comparison$bal_objects, 
    delta_dp = current_delta, 
    dataset_name = dataset_name, 
    base_dir = base_dir_delta,
    method_labels = METHOD_LABELS
  )
  
  if (!is.null(love_plot_comparison)) plot(love_plot_comparison)
  
  # Pass the metrics table specifically
  smd_plot <- create_smd_comparison(
    metrics_df   = results_comparison$metrics, 
    delta_dp     = current_delta, 
    dataset_name = dataset_name, 
    base_dir = full_inner_path_comparison,
    method_labels = method_labels
  )
  
  print(smd_plot)
  
  # End of your loop
  cat("\n--- Analysis Complete for delta =", current_delta, "---\n")
 
  # Stop capturing output
  sink() 
}

cat("\n--- OSIP Analysis Complete!!\n")
