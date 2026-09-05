# osip_functions.R

# 1. Start the timer
start_timer <- function() {
  return(Sys.time())
}

# 2. End the timer and print the message
end_timer <- function(start_time, label) {
  end_time <- Sys.time()
  # We calculate seconds first for the message, then minutes if needed
  duration_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
  duration_min <- duration_sec / 60
  
  cat(sprintf("\nStep [%s] took %.2f seconds (%.2f minutes)\n", 
              label, duration_sec, duration_min))
}

#' Create a clean set of default parameters for heuristics
get_default_search_params <- function(params, data_config, treatment_col, X_all,
                                       S_inv, delta_dp, k_bound) {
  
  list(
    # Configuration from PARAMS
    delta_dp = delta_dp,
    K_max    = k_bound,

    # Config from data_config
    id_var        = data_config$ID_VAR,
    treatment_col = treatment_col,
    data_config   = data_config,
    cov_cols      = data_config$ALL_COVARIATES,
    
    # Static Matrices
    X_all = X_all,
    S_inv = S_inv
  )
}

#' Load and Prepare Data
#' @param dataset_name String name of the dataset (e.g., "NSW_MIXTAPE", "LALONDE", "RHC")
#' @param datasets The reference list from your config file
#' @return A list containing 'data' and 'config'
#' 
load_and_prep_data <- function(dataset_name, sample_flag, datasets, params) {
  
  message(sprintf("--- Loading and Prepping: %s ---", dataset_name))
  
  data_full <- NULL
  data_subset <- NULL
  # n_control_sample <- params$N_CONTROL_SAMPLE
  # data_config <- NULL
  
  # ---------------------------------------------------------
  # PATH 1: LALONDE
  # ---------------------------------------------------------
  if (dataset_name == datasets$LALONDE) {
    
    # Debug
    # browser()
    
    # if (!require(MatchIt)) library(MatchIt()
    if (!require(causalsens)) install.packages("Personalized")
    library(causalsens)
    data(lalonde.psid)
    data_full <- lalonde
    data_config <- LALONDE_CONFIG
    
    # ---------------------------------------------------------
    # PATH 2: NSW_MIXTAPE
    # ---------------------------------------------------------
  } else if (dataset_name == datasets$NSW_MIXTAPE) {
    data(nsw_mixtape)
    data_full <- nsw_mixtape
    data_config <- NSW_MIXTAPE_CONFIG
    
    # ---------------------------------------------------------
    # PATH 3: RHC (Requires heavy transformation)
    # ---------------------------------------------------------
    
    
  } else if (dataset_name == datasets$RHC) {
    #dEBUG
    # browser()
    
    rhc_raw <- read.csv("https://hbiostat.org/data/repo/rhc.csv", header = TRUE)
    data_config <- RHC_CONFIG
    
    # 1. Transform categorical strings to binary/numeric
    rhc_clean <- rhc_raw %>%
      mutate(
        # Treatment and Outcome
        !!data_config$TREATMENT_VAR := ifelse(swang1 == "RHC", 1, 0),
        !!data_config$OUTCOME_VAR   := ifelse(dth30 == "Yes", 1, 0),
        # Creating standard Dummies used in configurations
        sex_male   = ifelse(sex == "Male", 1, 0),
        race_black = ifelse(race == "black", 1, 0)
      )
    
    # 2. Select only relevant columns to keep the dataframe lean
    required_cols <- c(
      data_config$TREATMENT_VAR, 
      data_config$OUTCOME_VAR, 
      data_config$NUMERIC_COVARIATES,
      data_config$FACTOR_COVARIATES
    )
    
    data_full <- rhc_clean %>%
      dplyr::select(all_of(required_cols)) %>%
      tidyr::drop_na(all_of(required_cols))
    
    cat(sprintf("RHC Prep: %d rows dropped due to NAs.\n", nrow(rhc_raw) - nrow(data_full)))
    
  } else if (dataset_name == datasets$LINDNER) {
    data(lindner)
    data_config <- LINDNER_CONFIG
    
    # Pre-processing: Literature standard often uses log-costs
    lindner$log_cardbill <- log(lindner$cardbill)
    
    # Update config to use the log version for better sensitivity convergence
    data_config$OUTCOME_VAR <- "log_cardbill"
    
    # Sample (Optional, but Lindner is small enough at N=996 to run full)
    data_full <- lindner
  } else if (dataset_name == datasets$JOBS) {
    
    # Debug
    # browser()
    
    data(jobs)
    
    data_config <- JOBS_CONFIG
    # JOBS II has some missing values; for matching, we typically use complete cases
    # 1. Strip all factor labels and convert everything to numeric at once
    jobs_clean <- as.data.frame(data.matrix(jobs))
    
    # 2. Remove NAs
    jobs_clean <- na.omit(jobs_clean)
    
    # 3. Proceed to your config
    data_full <- jobs_clean
    
    # Convert the factor levels to their numeric equivalent
    # (e.g., if level "1" is High School, it becomes the number 1)
    
    # Debug
    # browser()
    
    # jobs_clean$educ   <- ensure_numeric_col(jobs_clean, "educ")
    # jobs_clean$income <- ensure_numeric_col(jobs_clean, "income")
    
    # Now check the types again to be sure
    cat("\n--- Data Type Check (Post-Fix) ---\n")
    print(sapply(data_full[c("educ", "income")], class))
    
    data_full <- jobs_clean
    cat("JOBS II Dataset Loaded, here are the first few lines:\n")
    
    head(data_full)
    
    # Check which columns are problematic
    cat("\n--- Data Type Check ---\n")
    cols_to_check <- c(data_config$TREATMENT_VAR, data_config$NUMERIC_COVARIATES)
    type_check <- sapply(data_full[cols_to_check], class)
    print(type_check)
    
    # Check if any column is a 'factor' or 'character'
    is_not_numeric <- !sapply(data_full[cols_to_check], is.numeric)
    if(any(is_not_numeric)) {
      cat("Warning: These columns are NOT numeric and will cause NAs:\n")
      print(cols_to_check[is_not_numeric])
    }
    
    # --- 3. Preview Sample Sizes ---
    # n_treated <- sum(data_full$treat == 1)
    # n_control <- sum(data_full$treat == 0)
    # 
    # cat("- Treated units: ", n_treated, "\n")
    # cat("- Control units: ", n_control, "\n")
    # cat("- Total units:   ", nrow(jobs_clean), "\n")
    # print(names(jobs_clean))
  }
  
  else if (dataset_name == datasets$IDHP) {
    #debug
    # browser()
    
    data("ihdp", package = "bartcs")
    
    data_full <- ihdp
    data_config <- IDHP_CONFIG
  }
  
  else if (dataset_name == datasets$NHEFS) {
    #debug
    # browser()
    
    data(nhefs_complete)

    data_full <- nhefs_complete
    data_config <- NHEFS_CONFIG
  }
  
  else {
    stop(sprintf("Dataset '%s' not recognized. Check datasets.", dataset_name))
  }
  
  # ---------------------------------------------------------
  # FINAL STANDARDIZATION (Applies to all)
  # ---------------------------------------------------------
  # treatment_col <- data_config$TREATMENT_VAR
  
  # Ensure an ID column exists based on config
  if (!(data_config$ID_VAR %in% names(data_full))) {
    data_full[[data_config$ID_VAR]] <- 1:nrow(data_full)
  }
  
  # Convert treatment to numeric if it's not already
  # data_full[[treatment_col]] <- ensure_numeric_col(data_full, treatment_col)
  
  #dEBUG
  # browser()
  
  # IMPORTANT: Use double brackets [[...]] for dynamic variable extraction
  # treated_ids <- data_full$id[data_full[[treatment_col]] == 1]
  # control_ids <- data_full$id[data_full[[treatment_col]] == 0]
  
  # sample = TRUE
  # if (sample_flag) {
  #   set.seed(25)
  #   # treat_col = data_config$TREATMENT_VAR
  #   # --- DEFINE LOCAL IDs FOR SAMPLING ---
  #   # We create these here so they are available for the sample() function
  #   treated_ids <- data_full$id[data_full[[treatment_col]] == 1]
  #   control_ids <- data_full$id[data_full[[treatment_col]] == 0]
  #   
  #   n_treated_sample <- params$SAMPLE_SIZE
  #   
  #   #Debug
  #   # browser()
  #   
  #   # 2. Randomly sample n_treated_sample IDs from the treated group
  #   # Check if the desired sample size is feasible
  #   if (n_treated_sample > length(treated_ids)) {
  #     n_treated_sample <- length(treated_ids)
  #     cat(sprintf("INFO: n_treated_sample reduced to %d (max available treated units).\n", n_treated_sample))
  #   }
  #   
  #   sampled_treated_ids <- sample(treated_ids, n_treated_sample, replace = FALSE)
  #   
  #   # Too many contols - sampling them
  #   if (dataset_name == datasets$RHC) {
  #     sampled_control_ids <- sample(control_ids, n_control_sample, replace = FALSE)
  #     control_ids = sampled_control_ids
  #   }
  #   
  #   # 3. Combine sampled treated IDs with all control IDs
  #   subset_ids <- c(sampled_treated_ids, control_ids)
  #   
  #   # 4. Create the final working subset
  #   data_subset <- data_full[data_full$id %in% subset_ids, ]
  #   
  #   # --- ADD THIS LINE TO FIX THE ERROR ---
  #   # data_subset <- as.data.frame(data_subset)
  #   
  #   cat(sprintf("INFO: Original N_Treated: %d. New N_Treated: %d. N_Control: %d.\n", 
  #               length(treated_ids), n_treated_sample, length(control_ids)))    
  #   
  #   # Using the entire dataset, no sampling
  # } else {
  #   data_subset <- data_full
  # }
  
  #debug
  # browser()
  
  return(list(
    # data = data_subset,
    data = data_full,
    data_config = data_config
  ))
}

sample_dataset <- function(seed_in, treat_col, data_full, dataset_name, params) {
    set.seed(seed_in)
    # treat_col = data_config$TREATMENT_VAR
    # --- DEFINE LOCAL IDs FOR SAMPLING ---
    # We create these here so they are available for the sample() function
    treated_ids <- data_full$id[data_full[[treatment_col]] == 1]
    control_ids <- data_full$id[data_full[[treatment_col]] == 0]

    n_treated_sample <- params$SAMPLE_SIZE

    #Debug
    # browser()

    # 2. Randomly sample n_treated_sample IDs from the treated group
    # Check if the desired sample size is feasible
    if (n_treated_sample > length(treated_ids)) {
      n_treated_sample <- length(treated_ids)
      cat(sprintf("INFO: n_treated_sample reduced to %d (max available treated units).\n", n_treated_sample))
    }

    sampled_treated_ids <- sample(treated_ids, n_treated_sample, replace = FALSE)

    # Too many contols - sampling them
    if (dataset_name == datasets$RHC) {
      sampled_control_ids <- sample(control_ids, n_control_sample, replace = FALSE)
      control_ids = sampled_control_ids
    }

    # 3. Combine sampled treated IDs with all control IDs
    subset_ids <- c(sampled_treated_ids, control_ids)

    # 4. Create the final working subset
    data_subset <- data_full[data_full$id %in% subset_ids, ]

    # --- ADD THIS LINE TO FIX THE ERROR ---
    # data_subset <- as.data.frame(data_subset)

    cat(sprintf("INFO: Original N_Treated: %d. New N_Treated: %d. N_Control: %d.\n",
                length(treated_ids), n_treated_sample, length(control_ids)))
    return (data_subset)
}

process_propensity_scores <- function(data_subset, dataset_name, data_config, treatment_col, datasets, base_dir, title_suffix = "") {
  # 1. Determine columns based on dataset_name
  # treat_col <- if(dataset_name == "rhc") "swang1" else "treat"
  # treatment_col <- data_config$TREATMENT_VAR
  # Assuming covariates are known or handled via a different logic
  # For now, using all columns except ID and Treatment
  # TREAT_COL <- data_config$TREATMENT_VAR
  # NUM_COV   <- data_config$NUMERIC_COVARIATES
  # 
  # # 1. Select relevant columns (ID, PS, Treat, and Numerics)
  # cols_to_keep <- c(data_config$ID_VAR, "ps", TREAT_COL, NUM_COV)
  # 
  # ps_formula <- as.formula(paste(treat_col, "~", paste(covs, collapse = " + ")))
  
  # Debug
  # browser()
  
  # FIXED: Use the configured treatment and ALL covariates
  ps_formula <- as.formula(
    paste(treatment_col, "~",
          paste(data_config$ALL_COVARIATES, collapse = " + "))
  )
  
  # Use MatchIt to calculate PS (using logistic regression)
  m_out <- matchit(ps_formula, 
                   data = data_subset, 
                   method = NULL, 
                   distance = "glm")
  
  data_subset$ps <- m_out$distance
  
  # 3. Call the plotter
  # plot_ps_distribution(data_subset, dataset_name, treatment_col, datasets, base_dir, title_suffix)
  
  # return(list(data_subset))
  return(data_subset)
}

#' Propensity Score Processing Wrapper
#' @description Calculates PS, prepares sorted boundaries, and generates interactive Plotly figure.
#' @return A list containing the updated dataframe, the P_sorted vector, and the plotly object.
# process_propensity_scores <- function(data_subset, data_config) {
#   cat("\n=== Propensity Score Calculation & Visualization ===\n")
#   
#   # 1. Calculation using MatchIt
#   # Build formula dynamically using ALL_COVARIATES from config
#   ps_formula <- as.formula(
#     paste(data_config$TREATMENT_VAR, "~",
#           paste(data_config$ALL_COVARIATES, collapse = " + "))
#   )
#   
#   m_out <- matchit(ps_formula, 
#                    data = data_subset, 
#                    method = NULL, 
#                    distance = "glm")
#   
#   data_subset$ps <- m_out$distance
#   cat(sprintf("PS Range: [%.4f, %.4f]\n", min(data_subset$ps), max(data_subset$ps)))
#   
#   # 2. Boundary Preparation (P_sorted)
#   # Unique values + endpoints 0 and 1
#   P_sorted <- sort(unique(c(0, 1, data_subset$ps)))
#   
#   # 3. Preparation for Visualization
#   treat_col <- data_config$TREATMENT_VAR
#   data_subset$treat_label <- ifelse(data_subset[[treat_col]] == 1, "Treated", "Control")
#   data_subset$treat_jitter <- jitter(as.numeric(data_subset[[treat_col]]), amount = 0.05)
#   
#   # 4. Interactive Plotly Figure
#   fig <- plot_ly(
#     data = data_subset,
#     x = ~ps,
#     y = ~treat_jitter,
#     type = "scatter",
#     mode = "markers",
#     color = ~treat_label,
#     colors = c("blue", "red"),
#     hoverinfo = "text",
#     text = ~paste(
#       "Group:", treat_label,
#       "<br>Propensity Score:", round(ps, 4)
#     )
#   ) %>%
#     layout(
#       title = paste("Propensity Score Distribution:", data_config$DATASET_NAME),
#       xaxis = list(title = "Propensity Score", range = c(0, 1.0)),
#       yaxis = list(title = "", showticklabels = FALSE)
#     )
#   
#   # Return all three essential outputs
#   return(list(
#     data = data_subset,
#     P_sorted = P_sorted,
#     plot = fig
#   ))
# }

# Step 1: DP Solver for Initial Partition

# find_initial_partition <- function(treatment_scores, control_scores, k_bound, n_r,
#                                    delta_dp, cov_df, p_sorted, treatment_col, debug = FALSE) {
# 
#   # Inputs:
#   #   treatment_scores: Treated units' propensity scores (PS)
#   #   control_scores: Control units' PS
#   #   K_bound: Max number of allowed intervals
#   #   n_r: Max number of allowed REDUNDANT treated units (SLACK BUDGET)
#   #   delta: Max sub-interval size for counting empty intervals
# 
#   debug_local <- isTRUE(debug_glb)
# 
#   eps <- 1e-12
# 
#   # treatment_scores <- as.numeric(treatment_scores)
#   # control_scores <- as.numeric(control_scores)
# 
#   # Debug
#   # browser()
# 
#   # 🛑 CRITICAL FIX: Remove NAs/NaNs from input data
#   # treatment_scores <- treatment_scores[!is.na(treatment_scores) & !is.nan(treatment_scores)]
#   # control_scores <- control_scores[!is.na(control_scores) & !is.nan(control_scores)]
# 
#   # Debug
#   # browser()
# 
#   n_t <- length(treatment_scores)
#   if (n_t == 0) {
#     return(list(error = "No treated units provided"))
#   }
#   T_sorted <- sort(treatment_scores)
# 
#   # Build P_sorted: candidate endpoints (0, all unique PS from T and C, and 1)
#   # P_sorted <- sort(unique(c(0, treatment_scores, control_scores, 1)))
#   # P_sorted <- P_sorted[!is.na(P_sorted)]
# 
#   n_p <- length(p_sorted)
# 
#   if (debug_local) {
#     cat(sprintf("\n=== DP Setup ===\n"))
#     cat(sprintf("n_t=%d, n_c=%d, n_p=%d, K_bound=%d, n_r=%d, delta_dp=%.6f\n",
#                 n_t, length(control_scores), n_p, k_bound, n_r, delta_dp))
#   }
# 
#   # Debug
#   # browser()
# 
#   # Check 0 at index 1
#   if (abs(p_sorted[1] - 0) > eps) stop("ERROR: p_sorted must start at 0")
# 
#   # -------------------------------------------------------------------------
#   # 🛑 4D DP STATE INITIALIZATION: dp[i+1, v_idx, k+1, j+1]
#   # i: covered_treated_count (0..n_t)
#   # v_idx: p_sorted index (1..n_p)
#   # k: intervals used (0..K_bound)
#   # j: slack used (0..n_r)
#   # -------------------------------------------------------------------------
#   # dp <- array(Inf, dim = c(n_t + 1, n_p, K_bound + 1, n_r + 1))
# 
#   # 1. Re-initialize DP to store LISTS (must use a list array)
#   dp <- array(list(), dim = c(n_t + 1, n_p, k_bound + 1, n_r + 1))
#   backtrack <- array(list(), dim = c(n_t + 1, n_p, k_bound + 1, n_r + 1))
# 
#   # 2. Initialize the Base State with the list structure {cost=0, k_val=0}
#   # Assuming K=0 is stored in the k_val element (k_new=k+1 in the next state)
# 
# 
#   # Initialize all cells in backtrack with sentinel values
#   SENTINEL_DP <- list(cost = Inf, k_val = NA)
# 
#   for (idx in 1:(n_t + 1)) {
#     for (v_idx in 1:n_p) {
#       for (k_idx in 1:(k_bound + 1)) {
#         for (j_idx in 1:(n_r + 1)) {
#           dp[[idx, v_idx, k_idx, j_idx]] <- SENTINEL_DP
#         }
#       }
#     }
#   }
# 
#   # 3. Overwrite the Base State
#   dp[[1, 1, 1, 1]] <- list(cost = 0, k_val = 0)
# 
#   ## --- 2. GUARANTEED BACKTRACK INITIALIZATION ---
# 
#   backtrack <- array(list(), dim = c(n_t + 1, n_p, k_bound + 1, n_r + 1))
# 
#   # Define the sentinel structure for backtrack (e.g., all -1)
#   SENTINEL_BACKTRACK <- list(i_idx = -1, v_idx = -1, k_idx = -1, j_idx = -1)
# 
#   # 2. Initialize ALL cells in backtrack with sentinel values
#   for (idx in 1:(n_t + 1)) {
#     for (v_idx in 1:n_p) {
#       for (k_idx in 1:(k_bound + 1)) {
#         for (j_idx in 1:(n_r + 1)) {
#           backtrack[[idx, v_idx, k_idx, j_idx]] <- SENTINEL_BACKTRACK
#         }
#       }
#     }
#   }
# 
# 
#   # Base state: i=0, v=0 (idx 1), k=0, j=0 -> cost = 0
#   # dp[1, 1, 1, 1] <- 0
#   # backtrack[[1]][[1]][[1]][[1]] <- list(prev_i = NA, prev_v_idx = NA, prev_k = NA, prev_j = NA,
#   #                                       interval_start = 0, interval_end = 0,
#   #                                       n_treated = 0, n_control = 0, is_empty = NA, k_slack_used = 0)
# 
#   T_min_ps <- if (n_t > 0) T_sorted[1] else 1.0
#   T_min_idx <- which(p_sorted >= T_min_ps - eps)[1]
# 
#   # Debug
#   # browser()
# 
#   # Ensure this initialization is present before the loop:
#   # distance_to_T_min <- T_sorted[1]
#   # rho_init <- ceiling(distance_to_T_min / delta_dp)
#   # if (rho_init <= k_bound) {
#   #   # This allows Case A to start directly at the first treated unit
#   #   # by assuming we already paid the 'empty' price to get there.
#   #   dp[[1, T_min_idx, rho_init + 1, 1]] <- list(cost = 0, k_val = rho_init)
#   # }
# 
#   # MAIN DP LOOP
#   for (k in 0:(k_bound - 1)) { # 1. Intervals used so far
#     cat(sprintf("DP Progress: covering k=%d values\n", k))
#     # cat(sprintf("\n[START LOOP K] K=%d (Max K_bound=%d)\n", k, K_bound))
#     # cat(sprintf("[START LOOP I] K=%d, I=%d (Max n_t=%d)\n", k, i, n_t))
#     for (j in 0:n_r) { # 2. "redundant" treatment units (rtu) used so far (GLOBAL), out of the n_r budge
#       for (i in 0:n_t) { # 3. Treated covered so far
# 
#         # cat(sprintf("[START LOOP J] K=%d, I=%d, J=%d (Max n_r=%d)\n", k, i, j, n_r))
#         for (v_idx in 1:n_p) { # 4. Current PS endpoint index
# 
#           # ✅ CORRECT (applying is.finite to the cost value):
#           # Retrieve the cost value from the list structure using double brackets [[...]]
#           prev_cost_value <- dp[[i + 1, v_idx, k + 1, j + 1]]$cost
# 
#           if (!is.finite(prev_cost_value)) next
# 
#           # Debug
#           # browser()
# 
#           v_curr <- p_sorted[v_idx]
#           if (abs(v_curr - 1.0) < eps) next
# 
#           # Try right endpoints
#           for (v_next_idx in (v_idx + 1):n_p) {
#             v_next <- p_sorted[v_next_idx]
#             # cat("v_next is ", v_next, "\n")
#             if (!is.finite(v_next)) next
# 
#             # Check if v_next is the final boundary (1.0)
#             is_at_end <- abs(v_next - 1.0) < eps
# 
#             if (is_at_end) {
#               # Last bin is closed on both sides [L, 1.0] to capture units at 1.0
#               controls_in_interval <- sum(control_scores >= v_curr - eps & control_scores <= v_next + eps)
#             } else {
#               # Standard bins are Left-Closed, Right-Open [L, R)
#               # Using - eps on the right ensures a unit AT v_next belongs to the NEXT bin
#               controls_in_interval <- sum(control_scores >= v_curr - eps & control_scores < v_next - eps)
#             }
# 
#             interval_width <- v_next - v_curr
# 
#             # # --- Check 1: Delta Constraint ---
#             # if (interval_width > delta_dp) {
#             #   cat(sprintf("  [Delta Constraint violated] PS:[%.3f, %.3f]. Delta %.3f > %.3f. Too Wide (Rejecting).\n",
#             #               p_sorted[v_idx], p_sorted[v_next_idx], interval_width, delta_dp))
#             #   next
#             #   # Since P is sorted, if v_next_idx violates delta, we can often break an outer loop, but here, we just skip this transition.
#             #   # Ensure your outer loop over v_next_idx has a break condition if appropriate.
#             #   # For now, we continue to the next i/v_idx combination.
#             # }
# 
#             # cat("interval_width is ", interval_width, "\n")
# 
#             # compute controls in this PS interval once
#             # controls_in_interval <- sum(control_scores > v_curr & control_scores <= v_next)
# 
#             #####################################################################
#             # --- 1. PRE-CALC: DETERMINE n_treated_to_cover for (v_curr, v_next] ---
#             #####################################################################
# 
#             n_treated_to_cover <- 0
#             next_i_candidate <- i
# 
#             # Debug
#             # browser()
# 
#             if (i < n_t) {
# 
#               # Check if the *next* uncovered treated unit (T[i + 1]) lies inside the interval
#               t_next_ps <- T_sorted[i + 1]
# 
#               # This condition is the key to identifying a Case A transition:
#               # T[i+1] must be strictly after the start and at or before the end.
# 
#               # cat("t_next_ps is: ", t_next_ps, "\n")
#               # cat("v_curr is: ", v_curr, "\n")
#               # cat("v_next is: ", v_next, "\n")
#               # cat("v_next - v_curr is: ", v_next - v_curr, "\n")
# 
#               # t_next_ps is:  0.5739357
#               # v_curr is:  0.5766826
#               # v_next is:  0.5767194
#               # v_next - v_curr is:  3.678218e-05
# 
#               # bin_mask <- (cov_df$ps >= vprev - 1e-12) & (cov_df$ps <= vcurr + 1e-12)
#               cond = (t_next_ps >= v_curr - eps) && (t_next_ps <= v_next + eps) && (v_next - v_curr <= delta_dp + eps)
# 
#               # cond = (t_next_ps > v_curr) && (t_next_ps <= v_next + eps) && (v_next - v_curr <= delta_dp + eps)
# 
#               # cat("cond is: ", cond, "\n")
#               # Debug
#               # browser()
# 
#               if (cond) {
# 
#                 # Case A: T[i+1] is covered. Calculate how many contiguous units are covered.
#                 temp_j <- i
# 
# 
#                 # while ((temp_j + 1) <= n_t && is_in_interval_dp_style(T_sorted[temp_j + 1], v_curr, v_next)) {
#                 #   temp_j <- temp_j + 1
#                 # }
# 
#                 while ((temp_j + 1) <= n_t && (T_sorted[temp_j + 1] <= v_next + eps)) {
#                   temp_j <- temp_j + 1
#                 }
#                 n_treated_to_cover <- temp_j - i
#                 next_i_candidate <- temp_j
# 
#               } # end if next t in interval
#             } # end if i < n_t
# 
#             #####################################################################
#             # ------------------ CASE A: non-empty interval ---------------------
#             # We allow covering a contiguous block of newly-covered treated units
#             #####################################################################
# 
#             # cat("DEBUG_TREATED_COUNT: n_treated_to_cover is: ", n_treated_to_cover, "\n")
# 
#             # Check 1: Non-empty treated count AND Delta constraint satisfied (i.e., this interval is valid for Case A)
#             if (n_treated_to_cover > 0 && interval_width <= delta_dp + eps) {
# 
#               # --- Check 2: Fine Balance (Pre-Redundancy) ---
#               if (controls_in_interval >= n_treated_to_cover) {
# 
#                 controls_ps <- control_scores[control_scores > v_curr & control_scores <= v_next]
# 
#                 # Calculate Base Individual Costs (Distance to nearest control in interval)
#                 if (length(controls_ps) == 0) {
#                   individual_costs <- rep(Inf, n_treated_to_cover)
#                 } else {
#                   units_to_cover_ps <- T_sorted[(i + 1):next_i_candidate]
#                   individual_costs <- sapply(units_to_cover_ps, function(t) {
#                     d <- min(abs(controls_ps - t), na.rm = TRUE)
#                     return(d)
#                   })
#                 }
# 
#                 # cat("individual_costs are ", individual_costs, "\n")
# 
#                 # Sort costs descending to find the highest-cost units
#                 sorted_costs <- sort(individual_costs, decreasing = TRUE)
# 
#                 # cat("sorted_costs are ", sorted_costs, "\n")
# 
#                 t_curr <- length(units_to_cover_ps)
# 
#                 # 🛑 Iterate over possible rtu usage in THIS block (k_slack)
#                 max_k_slack <- min(n_r - j, t_curr)
# 
#                 # --- LOOP OVER REDUNDANCY (k_slack/r) ---
#                 for (k_slack in 0:max_k_slack) {
# 
#                   # Check 3: Fine Balance (Post-Redundancy) - Always TRUE here due to pre-check and non-negative k_slack
#                   T_non_redundant <- t_curr - k_slack
#                   # We only need to check T_non_redundant <= controls_in_interval, but since
#                   # controls_in_interval >= t_curr (the pre-check), and k_slack >= 0, this holds.
# 
#                   # Calculate cost for the block
#                   k_adjusted <- min(k_slack, t_curr)
#                   if (k_adjusted >= t_curr) {
#                     current_block_cost <- 0
#                   } else {
#                     current_block_cost <- sum(sorted_costs[(k_adjusted + 1):t_curr], na.rm = TRUE)
#                   }
# 
#                   # Retrieve previous state cost
#                   prev_cost_struct <- dp[[i + 1, v_idx, k + 1, j + 1]] # Structure containing cost, k, j
#                   prev_cost <- prev_cost_struct$cost # Assuming dp stores a list/struct
#                   # cat("DEBUG_COST1: Before tie check: prev_cost is: ", prev_cost, "\n")
# 
# 
#                   if (is.finite(prev_cost)) {
#                     new_cost <- prev_cost + current_block_cost
#                     # cat("DEBUG_COST2: Before tie check: new_cost is: ", new_cost, "current_block_cost is: ", current_block_cost, "\n")
# 
#                     k_new <- k + 1
#                     j_new <- j + k_slack # Update total slack used
# 
#                     # --- Check 4: K_bound & DP Update (Tie-breaker) ---
#                     # cat("DEBUG_VALS: Before tie check: k_new is: ", k_new, "\n")
#                     if (k_new <= k_bound) {
# 
#                       best_cost <- dp[[next_i_candidate + 1, v_next_idx, k_new + 1, j_new + 1]]$cost
#                       # Assuming best_k is stored in the DP structure too, or derived
#                       best_k <- dp[[next_i_candidate + 1, v_next_idx, k_new + 1, j_new + 1]]$k_val
# 
#                       is_better <- FALSE
# 
#                       # 1. Cost check
#                       # if (new_cost < best_cost - eps) {
#                       #   is_better <- TRUE
#                       #   # 2. Tie-breaker on K (prefer smaller K)
#                       # } else if (abs(new_cost - best_cost) < eps) {
#                       #   # If best_k is Inf or NA, assume the path is new/unvisited
#                       #   if (is.na(best_k) || k_new < best_k) {
#                       #     is_better <- TRUE
#                       #   }
#                       #   # NOTE: We DO NOT check j_new for tie-breaking.
#                       # }
# 
#                       # 1. Cost check: strictly better cost is always preferred
#                       if (new_cost < best_cost - eps) {
#                         is_better <- TRUE
#                         # 2. Tie-breaker: If costs are equal, accept if unvisited (NA k_val)
#                       } else if (abs(new_cost - best_cost) < eps) {
#                         # Note: Since the index fixes K_new, we only need to accept the first time
#                         if (is.na(best_k)) { # k_val is NA only if cost=Inf or if it's the base state.
#                           is_better <- TRUE
#                         }
#                       }
# 
#                       # Debug Print Block
#                       # cat(sprintf("\n  [CASE A Check] K'=%d->K=%d, I'=%d->I=%d, J'=%d->J=%d (r=%d)\n",
#                       #             k, k_new, i, next_i_candidate, j, j_new, k_slack))
#                       # cat(sprintf("                 PS:[%.3f, %.3f] | T_non_r=%d, C=%d\n",
#                       #             p_sorted[v_idx], p_sorted[v_next_idx], T_non_redundant, controls_in_interval))
#                       # cat(sprintf("                 Cost=%.4f (Prev=%.4f, Block=%.4f) | Best Cost=%.4f\n",
#                       #             new_cost, prev_cost, current_block_cost, best_cost))
#                       # cat(sprintf("                 Is_Better: %s\n", is_better))
# 
#                       # --- 5. UPDATE DP TABLE ---
#                       if (is_better) {
#                         dp[[next_i_candidate + 1, v_next_idx, k_new + 1, j_new + 1]] <- list(cost = new_cost, k_val = k_new) # Store cost and K
# 
#                         # Backtrack Update (Your existing backtrack list assignment goes here)
#                         backtrack[[next_i_candidate + 1, v_next_idx, k_new + 1, j_new + 1]] <- list(
#                           prev_i = i, prev_v_idx = v_idx, prev_k = k, prev_j = j,
#                           interval_start = p_sorted[v_idx], interval_end = p_sorted[v_next_idx], # Using PS values is clearer
#                           n_treated = t_curr, n_control = controls_in_interval,
#                           is_empty = FALSE, k_slack_used = k_slack
#                         )
#                         # cat("                 *** STATE UPDATED (Case A) ***\n")
#                       }
#                     }
#                   }
#                 } # END k_slack loop
#                 # } else {
#                 #   # Fine Balance (Pre-redundancy) failed
#                 #   cat(sprintf("  [CASE A FAIL] PS:[%.3f, %.3f]. FB fail T=%d > C=%d. Skipping redundancy loop.\n",
#                 #               p_sorted[v_idx], p_sorted[v_next_idx], n_treated_to_cover, controls_in_interval))
#                 # }
#               }
#             } #Ends if (n_treated_to_cover > 0 && delta_dp <= interval_width) (Case A)
# 
#             #####################################################################
#             # ------------------ CASE B: Empty Interval ---------------------
#             #####################################################################
#             else if (n_treated_to_cover == 0) {
# 
#               # Debug
#               # browser()
# 
#               # CRITICAL FIX: For the FINAL completion jump (i == n_t),
#               # we must ENSURE j == n_r (all slack must be used)
#               if (i == n_t) {
#                 # At this point, all treatments are covered.
#                 # Empty intervals are ONLY allowed if we've exhausted the slack budget.
#                 if (j != n_r) {
#                   # Skip this transition: we cannot complete without using all slack
#                   next
#                 }
#               } else if (!large_gaps && i > 0 && i < n_t) {
#                 # Mid-coverage empty intervals are not allowed
#                 next
#               }
# 
#               # Ensure current state is reachable
#               prev_cost_struct <- dp[[i + 1, v_idx, k + 1, j + 1]]
#               prev_cost_dist <- prev_cost_struct$cost # Assuming structure access
#               if (!is.finite(prev_cost_dist)) {
#                 next
#               }
# 
#               # CRITICAL CHECK 1: Already complete (v_curr == 1.0)
#               if (abs(v_curr - 1.0) < eps) {
#                 # The entire space is covered. Exit the v_next_idx loop for this v_idx.
#                 break
#               }
# 
#               # (Your existing NA guards and v_next checks go here)
#               if (is.na(v_next) || !is.finite(v_next)) {
#                 next
#               }
# 
#               # 🛑 NEW CHECK FOR THE INITIAL SEGMENT (0 to T_min)
#               if (i == 0 && v_curr < T_min_ps - eps) {
#                 # Check if we are jumping past an intermediate point before T_min
#                 if (v_next_idx < T_min_idx) {
#                   # This is a jump (v_curr, v_next] that ends before T_min
#                   # We want to skip this jump unless it uses the max allowed delta
#                   # or if v_next is the point T_min_ps (handled by Case A later)
# 
#                   # The logic to enforce the paper's constraint is usually:
#                   # If v_next < T_min, the empty jump must be of length delta.
#                   # But since P_sorted only contains candidate endpoints, we must
#                   # allow the DP to jump as far as possible *up to* the index T_min_idx.
# 
#                   # The simplest way to enforce "efficiency" is to only allow
#                   # the jump to T_min or the largest point < T_min that's a multiple of delta
# 
#                   # Since P_sorted *must* contain T_min (at T_min_idx), the DP
#                   # will naturally allow a large jump to T_min_idx via Case A.
# 
#                   # Here, we only need to prevent empty jumps to non-T-unit points
#                   # between 0 and T_min that waste a K budget slot.
# 
#                   # A simpler, paper-aligned approach is to only allow a jump
#                   # to a point P_sorted[v_next_idx] if P_sorted[v_next_idx]
#                   # is at least T_min, OR if P_sorted[v_next_idx] is the
#                   # smallest multiple of delta greater than v_curr.
# 
#                   # However, since T_min_ps is guaranteed to be covered
#                   # by a Case A transition, the main requirement is just to
#                   # ensure no treated unit is covered (which you have), AND to
#                   # allow the empty jump to reach *up to* T_min.
# 
#                   # A minimal enforcement: if v_next is not T_min and is before T_min,
#                   # we should skip this jump unless it is of length delta.
# 
#                   # Given that Case A will cover T_min, we only need to
#                   # consider empty jumps that skip T_min and land at v_next > T_min.
# 
#                   # *Final simpler fix based on your constraint*: If $v_{curr}=0$, skip
#                   # any jump that ends *before* $T_{min}$ and is not $\delta$-long.
#                   if (abs(v_curr - 0) < eps && v_next < T_min_ps - eps) {
#                     # If we are jumping from 0 to a point before T_min,
#                     # it must be either T_min (Case A) or be exactly delta length
#                     if (prev_cost_dist < delta_dp - eps) next # Skip: wastes K budget slot
#                   }
#                 }
#               }
#               # 🛑 END NEW CHECK
# 
#               # Debug
#               # browser()
# 
#               distance <- v_next - v_curr
#               # Assuming delta is delta_dp globally
#               rho <- ceiling(distance / delta_dp) # Number of delta-sized sub-intervals needed
# 
#               if (is.na(rho) || !is.finite(rho) || rho < 1) next
# 
#               k_new <- k + rho
#               j_new <- j # Slack used does NOT change (k_slack_used = 0)
# 
#               if (k_new <= k_bound) {
# 
#                 # 🛑 FIX: Use double brackets [[...]] to retrieve the specific list element
#                 # This ensures 'best_cost' and 'best_k' are extracted safely from the list structure.
#                 new_cost_dist <- prev_cost_dist # Empty interval cost is 0
#                 best_cost <- dp[[i + 1, v_next_idx, k_new + 1, j_new + 1]]$cost
#                 best_k <- dp[[i + 1, v_next_idx, k_new + 1, j_new + 1]]$k_val
# 
#                 # ========================================
#                 # ADDITIONAL FIX: Empty interval from 0 to first treated unit
#                 # ========================================
# 
#                 # Your code mentions an empty interval from 0 to 0.1090 that doesn't show.
#                 # This is likely because the DP starts at i=0, v_idx=1 (PS=0), k=0, j=0
#                 # and the FIRST transition is a Case A (non-empty) from [0, 0.1090].
# 
#                 # However, the PAPER may require an empty interval if there are no
#                 # treatments in [0, T_min). Check if T_min > 0:
# 
#                 # cat(sprintf("\n=== First Treated Unit ===\n"))
#                 # cat(sprintf("T_min PS: %.6f\n", T_sorted[1]))
#                 # cat(sprintf("If T_min > 0, should there be an initial empty interval?\n"))
# 
#                 # If the paper requires [0, T_min) to be covered by an empty interval,
#                 # you need to add special initialization logic:
# 
#                 if (T_sorted[1] > eps) {
#                   # There's a gap from 0 to the first treated unit
#                   # Initialize the DP to allow an empty interval here
# 
#                   T_min_ps <- T_sorted[1]
#                   T_min_idx <- which(p_sorted >= T_min_ps - eps)[1]
# 
#                   distance_to_T_min <- T_min_ps
#                   rho_initial <- ceiling(distance_to_T_min / delta_dp)
# 
#                   if (rho_initial <= k_bound) {
#                     # Create an empty interval from 0 to T_min
#                     k_new_init <- rho_initial
#                     j_new_init <- 0
# 
#                     dp[[1, T_min_idx, k_new_init + 1, j_new_init + 1]] <- list(
#                       cost = 0,
#                       k_val = k_new_init
#                     )
# 
#                     backtrack[[1, T_min_idx, k_new_init + 1, j_new_init + 1]] <- list(
#                       prev_i = 0,
#                       prev_v_idx = 1,
#                       prev_k = 0,
#                       prev_j = 0,
#                       interval_start = 0,
#                       interval_end = T_min_ps,
#                       n_treated = 0,
#                       n_control = sum(control_scores > 0 & control_scores <= T_min_ps),
#                       is_empty = TRUE,
#                       rho = rho_initial,
#                       k_slack_used = 0
#                     )
# 
#                     # cat(sprintf("✅ Initialized empty interval [0, %.4f] with K=%d\n",
#                     #             T_min_ps, rho_initial))
#                   }
#                 }
# 
# 
#                 # --- Tie-breaking logic (Cost then K) ---
#                 is_better <- FALSE
#                 if (new_cost_dist < best_cost - eps) {
#                   is_better <- TRUE
#                 } else if (abs(new_cost_dist - best_cost) < eps) {
#                   if (is.na(best_k) || k_new < best_k) {
#                     is_better <- TRUE
#                   }
#                 }
# 
#                 # Debug Print: Post-Tie-Breaker
#                 # cat(sprintf("                 Cost=%.4f (Best Cost=%.4f) | Is_Better: %s\n",
#                 #             new_cost_dist, best_cost, is_better))
#                 #
#                 if (is_better) {
#                   # Update DP table
#                   # dp[i + 1, v_next_idx, k_new + 1, j_new + 1] <- list(cost = new_cost_dist, k_val = k_new)
#                   dp[[i + 1, v_next_idx, k_new + 1, j_new + 1]] <- list(cost = new_cost_dist, k_val = k_new)
# 
#                   # Backtrack Update
#                   backtrack[[i + 1, v_next_idx, k_new + 1, j_new + 1]] <- list(
#                     prev_i = i, prev_v_idx = v_idx, prev_k = k, prev_j = j,
#                     interval_start = v_curr, interval_end = v_next,
#                     n_treated = 0, n_control = controls_in_interval,
#                     is_empty = TRUE, rho = rho, k_slack_used = 0
#                   )
#                 }
#                 # Debug
#                 # browser()
#                 break
#               }
#             } # Ends else if (n_treated_to_cover == 0)
#           } # End v_next_idx loop
#         } # End v_idx loop
#       } # END I LOOP (Slack used so far)
#     } # End J loop
#   } # End K loop
# 
#   cat("\n=== Diagnostic: Check States Near PS=1.0 ===\n")
#   v_1_idx <- which(abs(p_sorted - 1.0) < eps)[1]
# 
#   if (!is.na(v_1_idx)) {
#     cat(sprintf("PS=1.0 is at index %d\n", v_1_idx))
# 
#     # Check all possible paths to v_1_idx
#     for (k_try in 1:k_bound) {
#       for (j_try in 0:n_r) {
#         dp_cell <- dp[[n_t + 1, v_1_idx, k_try + 1, j_try + 1]]
#         cost_val <- dp_cell$cost
# 
#         if (!is.infinite(cost_val)) {
#           cat(sprintf("  Reachable at [k=%d, j=%d]: Cost=%.6f\n",
#                       k_try, j_try, cost_val))
#         }
#       }
#     }
#   } else {
#     cat("❌ PS=1.0 not found in p_sorted!\n")
#   }
# 
#   # ========================================
#   # CRITICAL FIX: Prioritize Complete Solutions
#   # ========================================
# 
#   # Replace your "FIND BEST SOLUTION" section with this:
# 
#   # --- Search entire final row (i = n_t) for the optimal endpoint ---
#   cat("\n=== Searching for Optimal Final State (i = n_t) ===\n")
# 
#   # First, identify the index for PS=1.0
#   v_1_idx <- which(abs(p_sorted - 1.0) < 1e-6)[1]
#   if (is.na(v_1_idx)) v_1_idx <- n_p
# 
#   best_cost <- Inf
#   best_k <- NA
#   best_j <- NA
#   best_v_idx <- NA
# 
#   # PRIORITY 1: Find the best COMPLETE solution (reaching PS=1.0)
#   for (k_try in 1:k_bound) {
#     for (j_try in 0:n_r) {
# 
#       dp_cell <- dp[[n_t + 1, v_1_idx, k_try + 1, j_try + 1]]
#       cost_val <- dp_cell$cost
# 
#       if (!is.infinite(cost_val) && !is.na(cost_val)) {
# 
#         # Tie-breaking for complete solutions: 1. Cost, 2. K, 3. J
#         is_better <- FALSE
# 
#         if (cost_val < best_cost - eps) {
#           is_better <- TRUE
#         } else if (abs(cost_val - best_cost) < eps) {
#           if (is.na(best_k) || k_try < best_k) {
#             is_better <- TRUE
#           } else if (k_try == best_k && j_try < best_j) {
#             is_better <- TRUE
#           }
#         }
# 
#         if (is_better) {
#           best_cost <- cost_val
#           best_k <- k_try
#           best_j <- j_try
#           best_v_idx <- v_1_idx  # Complete solution endpoint
#         }
#       }
#     }
#   }
# 
#   if (is.na(best_k)) {
#     return(list(error = "No feasible solution found covering all treated units"))
#   }
# 
#   cat(sprintf("✅ Optimal State Found:\n"))
#   cat(sprintf("   PS Endpoint: %.6f (Index: %d)\n", p_sorted[best_v_idx], best_v_idx))
#   cat(sprintf("   Cost: %.6f, K: %d, J: %d\n", best_cost, best_k, best_j))
# 
#   # ========================================
#   # FIX 2: Update reconstruction call
#   # ========================================
# 
#   # Replace the backtrack extraction section with:
# 
#   final_i_idx <- n_t + 1
#   final_k_idx <- best_k + 1
#   final_j_idx <- best_j + 1
# 
#   cat(sprintf("\n=== Extracting Final Backtrack Record ===\n"))
#   cat(sprintf("State: [i=%d, v_idx=%d, k=%d, j=%d]\n",
#               n_t, best_v_idx, best_k, best_j))
# 
#   final_backtrack_record <- backtrack[[final_i_idx, best_v_idx, final_k_idx, final_j_idx]]
# 
#   if (is.null(final_backtrack_record)) {
#     cat("ERROR: Final backtrack record is NULL.\n")
#     cat("This indicates the DP path was not properly stored.\n")
#     return(list(error = "Backtrack path incomplete"))
#   }
# 
#   cat(sprintf("Source State: i'=%d, v_idx'=%d (PS: %.6f), k'=%d, j'=%d\n",
#               final_backtrack_record$prev_i,
#               final_backtrack_record$prev_v_idx,
#               p_sorted[final_backtrack_record$prev_v_idx],
#               final_backtrack_record$prev_k,
#               final_backtrack_record$prev_j))
# 
#   # ========================================
#   # FIX 3: Handle incomplete paths to 1.0
#   # ========================================
# 
#   # Add this check after finding the optimal state:
# 
#   if (p_sorted[best_v_idx] < 1.0 - 1e-6) {
#     cat(sprintf("\n⚠️ WARNING: Optimal path ends at PS=%.6f (not 1.0)\n",
#                 p_sorted[best_v_idx]))
#     cat("The path does not reach the upper boundary.\n")
#     cat("This may indicate:\n")
#     cat("  1. K_bound is too small to allow the final empty jump\n")
#     cat("  2. The empty interval logic is rejecting valid transitions\n")
#     cat("  3. There's a constraint preventing completion\n")
# 
#     # Calculate what would be needed to complete
#     remaining_distance <- 1.0 - p_sorted[best_v_idx]
#     rho_needed <- ceiling(remaining_distance / delta_dp)
#     k_needed <- best_k + rho_needed
# 
#     cat(sprintf("  To complete: need %d more intervals (total K=%d)\n",
#                 rho_needed, k_needed))
# 
#     if (k_needed > k_bound) {
#       cat(sprintf("  ❌ k_bound=%d is insufficient (need at least %d)\n",
#                   k_bound, k_needed))
#     }
#   }
# 
#   # Modify the reconstruction_result call to use best_v_idx:
# 
#   # Debug
#   # browser()
# 
#   reconstruction_result <- reconstruct_intervals_and_find_redundant(
#     backtrack = backtrack,
#     p_sorted = p_sorted,
#     cov_df = cov_df,
#     n_r = params$N_R,
#     id_var = data_config$ID_VAR,
#     treatment_col = treatment_col,
#     J_bound = params$N_R,
#     dp = dp,
#     final_v_idx = best_v_idx,  # NEW: Pass the actual endpoint
#     final_k = best_k,
#     final_j = best_j
#   )
# 
# 
#   return(list(
#     cost_dist_optimal = best_cost,
#     k_dist_optimal = best_k,
#     n_r_used_optimal = best_j,
#     intervals_dist_optimal = reconstruction_result$intervals,
#     # 🛑 CRITICAL CHANGE: Use the unique ID set 🛑
#     redundant_unit_ids = reconstruction_result$redundant_unit_ids,
#     n_t_covered = n_t
#   ))
# }

get_treatment_map <- function(data_subset, treatment_col, id_var) {
  
  # Create the Treatment Map (Treated Units + IDs + Scores)
  treatment_map <- data_subset %>%
    filter(!!sym(treatment_col) == 1) %>%
    # Use "ps" directly since that is your column name
    dplyr::select(!!sym(id_var), ps) %>% 
    arrange(ps)
  
  return (treatment_map)
  
}

get_treatment_control_scores <- function(treatment_map, data_subset, id_var, treatment_col) {
  # For Treated:
  treatment_ids <- treatment_map[[id_var]]
  treatment_scores <- treatment_map$ps
  # 🆕 Explicitly name the vector
  names(treatment_scores) <- treatment_ids
  
  valid_mask <- !is.na(treatment_scores) & !is.nan(treatment_scores)
  treatment_scores <- treatment_scores[valid_mask]
  
  # For Controls:
  controls_df <- data_subset[data_subset[[treatment_col]] == 0, ]
  control_scores <- controls_df$ps
  # 🆕 Explicitly name the vector
  names(control_scores) <- controls_df[[id_var]]
  
  valid_ctrl <- !is.na(control_scores) & !is.nan(control_scores)
  control_scores <- control_scores[valid_ctrl]
  
  if(length(treatment_scores) == 0) {
    return(list(error = "No treated units provided."))
  }
  
  return(list(treatment_scores = treatment_scores, control_scores = control_scores))
}

# get_treatment_control_scores <- function(treatment_map, data_subset, id_var, treatment_col) {
#   # ... (your existing treatment logic) ...
#   
# # 2. Extract the parallel vectors
#   treatment_ids <- treatment_map[[id_var]]
#   treatment_scores <- treatment_map$ps
# 
#   # 3. Apply your NA/NaN cleaning
#   valid_mask <- !is.na(treatment_scores) & !is.nan(treatment_scores)
#   treatment_ids <- treatment_ids[valid_mask]
#   treatment_scores <- treatment_scores[valid_mask]
# 
#   # Extract controls with IDs
#   controls_df <- data_subset[data_subset[[treatment_col]] == 0, ]
#   control_scores <- controls_df$ps
#   names(control_scores) <- controls_df[[id_var]] # 🆕 Assign names for ID lookup
#   
#   # Cleaning
#   valid_ctrl <- !is.na(control_scores) & !is.nan(control_scores)
#   control_scores <- control_scores[valid_ctrl]
#   
#   cat(sprintf("Debug: Found %d Treated and %d Control units.\n", 
#               length(treatment_scores), length(control_scores)))
#   
#   if(length(treatment_scores) == 0) {
#     return(list(error = "No treated units provided."))
#   }
#   
#   return(list(treatment_scores = treatment_scores, control_scores = control_scores))
# }

# get_treatment_control_scores <- function(treatment_map, data_subset, id_var, treatment_col) {
#   # 2. Extract the parallel vectors
#   treatment_ids <- treatment_map[[id_var]]
#   treatment_scores <- treatment_map$ps
#   
#   # 3. Apply your NA/NaN cleaning
#   valid_mask <- !is.na(treatment_scores) & !is.nan(treatment_scores)
#   treatment_ids <- treatment_ids[valid_mask]
#   treatment_scores <- treatment_scores[valid_mask]
#   
#   control_scores <- data_subset$ps[data_subset[[treatment_col]] == 0]
#   control_scores <- control_scores[!is.na(control_scores) & !is.nan(control_scores)]
#   
#   # DEBUG PRINTS to catch the "No treated units" error
#   cat(sprintf("Debug: Found %d Treated units and %d Control units.\n", 
#               length(treatment_scores), length(control_scores)))
#   
#   if(length(treatment_scores) == 0) {
#     return(list(error = "No treated units provided. Check if TREATMENT_VAR matches data and is numeric 1."))
#   }
#   else {
#     return(list(treatment_scores=treatment_scores, control_scores=control_scores))
#   }
# }

run_osip_step1 <- function(data_subset, params, id_var, p_sorted, 
                            treatment_col, treatment_scores, 
                            control_scores, k_bound, delta_dp, 
                            left_border, right_border, dist_power = 1) {
                            
  cat("\n--- osip Step 1: Solving DP Partition ---\n")
  
  # Debug
  # browser()
  
  # DEBUG: Check if 'ps' and Treatment column exist
  if(!"ps" %in% colnames(data_subset)) stop("CRITICAL: Column 'ps' missing from data_subset")
  
  start_time <- Sys.time()
  
  # Debug
  # browser()
  
  solution <- tryCatch({
    find_initial_partition(
      treatment_scores = treatment_scores, 
      control_scores = control_scores, 
      k_bound  = k_bound, 
      # n_r      = n_r, 
      delta_dp    = delta_dp, 
      treatment_col = treatment_col, 
      cov_df   = data_subset, 
      p_sorted = p_sorted,
      left_border = left_border, 
      right_border = right_border,
      dist_power = dist_power
    )
  }, error = function(e) {
    return(list(error = as.character(e)))
  })
  
  return(solution)
}

#' Step 2: Heuristic Local Search (Refining Covariate Balance)
#' Step 2: Heuristic Local Search Wrapper
run_osip_step2_heuristics <- function(osip_res_step_1, initial_cost, initial_bounds, 
                                       initial_intervals, cov_df, params, 
                                       data_config, p_sorted, 
                                       treatment_col, S_inv, delta_dp, left_border,right_border) {
  
  right_border
  
  # 4. Prepare Search Engine
  opt_funcs <- list(
    "local_search"        = local_search_intervals,
    "enhanced_ls"         = enhanced_local_search,
    "simulated_annealing" = simulated_annealing_intervals
  )
  
  # Debug
  # browser()
  
  initial_k <- length(initial_intervals)
  
  defaults <- get_default_search_params(
    params = params, 
    data_config = data_config, 
    treatment_col = treatment_col, 
    S_inv = S_inv,
    delta_dp = delta_dp
    )
  
  # 5. Run and Process Searches
  final_results <- list()
  
  for (method in params$SEARCHES_TO_USE) {
    cat(sprintf("\nLaunching %s...", method))
    
    # Debug
    # browser()
    
    # Define what changes for this specific run
    overrides <- list(
      initial_intervals = initial_intervals,
      initial_bounds = initial_bounds,
      initial_cost = initial_cost,
      initial_k = initial_k,
      cov_df = cov_df,
      p_sorted = p_sorted,
      fixed_redundant_unit_ids = osip_res_step_1$redundant_unit_ids,
      left_border = left_border,
      right_border = right_border
    )
    
    #Debug
    # browser()
    
    # Execute
    raw_res <- run_search(method, defaults, overrides, opt_funcs)
    
    #Debug
    # browser()
    
    # Process & Store
    final_results[[method]] <- process_and_report_search(
      ls_res = raw_res, 
      search_name = method, 
      initial_cost = initial_cost,
      initial_intervals = initial_intervals,
      initial_k = length(initial_intervals)
    )
  }

  # Debug
  # browser()
  
  # --- Inside run_osip_step2_heuristics ---
  # vapply ensures we get a logical vector and throws an error if a list is returned
  improved_mask <- vapply(final_results, function(res) {
    # Handle potential NULLs or non-logical status values safely
    isTRUE(res$status) 
  }, logical(1))
  
  #Debug
  # browser()
  
  if (!any(improved_mask)) {
    # Return a "Null-Object" list instead of literal NULL
      # 1. Create a basic list 
      fail_res <- list(data_matched = NULL, name = "No Improvement")
      
      # 2. Assign the same attribute you use in the success case
      attr(fail_res, "is_step2") <- FALSE
      
      return(fail_res)
    }
  
  #Debug
  # browser()
  
  # 2. Extract ONLY the results that improved
  improved_results <- final_results[improved_mask]
  
  #Debug
  # browser()
  
  # 3. Find the index of the one with the absolute minimum cost
  best_sub_idx <- which.min(sapply(improved_results, function(res) res$cost))
  
  # 4. Extract the SINGLE best list object
  best_res <- improved_results[[best_sub_idx]]
  best_res_name <- best_res$name
  
  # 5. Assign attributes to the single object for the main program to read
  attr(best_res, "best_method") <- best_res_name
  attr(best_res, "is_step2")     <- TRUE
  
  # 6. RETURN ONLY THE BEST OBJECT
  return(best_res)
}

# visualize_osip_step_results <- function(data_subset, osip_res, title_suffix = "", data_config, delta_dp, dataset_name, orig_cost_step1, base_dir)
#   {
#   
#   # 1. Basic Null Check
#   if (is.null(osip_res)) return(NULL)
#   
#   # Debug 
#   # browser()
#   
#   # Check for Step 2: It will be a named list containing "local_search", "enhanced_ls", etc.
#   # OR it will have the attribute we set above
#   
#   
#   final_intervals <- NULL
#   final_redundant <- NULL
#   method_used     <- ""
#   
#   # Plot step 2 results
#   if (title_suffix == "Optimized") {
#     
#     is_step2 <- !is.null(attr(osip_res, "is_step2")) 
#     
#     if (is_step2) {
#       # Extract costs and find the absolute minimum
#       # This handles the case where Enhanced LS (0.084) beats SA (0.085)
#       costs <- sapply(osip_res, function(x) x$cost)
#       best_idx <- which.min(costs)
#       best_method_name <- names(osip_res)[best_idx]
#       
#       final_intervals <- osip_res[[best_idx]]$intervals
#       method_used <- paste("Step 2 - Best Method:", toupper(best_method_name), 
#                            sprintf("(Cost: %.4f)", costs[best_idx]))
#       } 
#     } else {
#     # Plot step 2 results
#     final_intervals <- osip_res$intervals_dist_optimal
#     method_used <- paste("Step 1 DP", sprintf("(Cost: %.4f)", orig_cost_step1))
#   }
#   
#   # 3. Final Validation & Plot
#   if (is.null(final_intervals)) {
#     message("No intervals found to plot.")
#     return(NULL)
#   }
#   
#   plot_partition_with_bins(
#     data_subset        = data_subset,
#     intervals_df       = final_intervals, 
#     redundant_unit_ids = final_redundant, 
#     title = sprintf("osip Partition for delta = %.2f: %s %s", current_delta, method_used, title_suffix),
#     data_config         = data_config,
#     delta_dp = delta_dp,
#     dataset_name = dataset_name,
#     base_dir = base_dir
#   )
# }

#' Standardize Covariates for Distance-Based Operations
#' @description Calculates pooled within-group SDs and scales numeric covariates.
#' @param data_subset The dataframe containing 'ps' and original covariates.
#' @param data_config The configuration object.
#' @return A list containing 'original' and 'standardized' dataframes.
#' 
prepare_standardized_data <- function(data_subset, data_config, treatment_col) {
  cat("\n--- Preparing data for Step 2 ---\n")
  cat("Standardizing covariates using pooled within-group SDs...\n")
  temp_df <- data_subset
  covariates <- data_config$NUMERIC_COVARIATES
  
  # 1. Select relevant columns (ID, PS, Treat, and Numerics)
  # cols_to_keep <- c(data_config$ID_VAR, "ps", treatment_col, covariates, data_config$OUTCOME_VAR)
  # cov_df_lim <- data_subset[, cols_to_keep, drop = FALSE]
  
  # Ensure Treatment is numeric
  temp_df[[treatment_col]] <- ensure_numeric_col(temp_df, treatment_col)
  
  # --- FIX 1: Force all target covariates to be numeric before checking ---
  for (var in covariates) {
    # as.character then as.numeric handles factors correctly
    temp_df[[var]] <- ensure_numeric_col(temp_df, var)
  }
    
    # as.numeric(as.character(cov_df_original[[treatment_col]]))
  
  # 2. Split data to calculate pooled variance
  control_data <- temp_df[temp_df[[treatment_col]] == 0, covariates, drop = FALSE]
  treated_data <- temp_df[temp_df[[treatment_col]] == 1, covariates, drop = FALSE]
  
  # Safety Check: Numeric only
  if (!all(sapply(control_data, is.numeric)) || !all(sapply(treated_data, is.numeric))) {
    stop("ERROR: Non-numeric columns found in standardization target columns.")
  }
  
  # 3. Calculate Pooled SDs
  # Formula: sqrt((Var_c + Var_t) / 2)
  cov_sds <- sapply(covariates, function(var_name) {
    v_c <- var(control_data[[var_name]], na.rm = TRUE)
    v_t <- var(treated_data[[var_name]], na.rm = TRUE)
    s <- sqrt((v_c + v_t) / 2)
    # Prevent division by zero
    if (is.na(s) || s < 1e-10) return(1) else return(s)
  })
  
  # 4. Create Standardized Version
  cov_df_standardized <- temp_df # Copy structure
  for (var in covariates) {
    cov_df_standardized[[var]] <- temp_df[[var]] / cov_sds[var]
  }
  
  cat("✓ Standardization complete.\n")
  return (cov_df_standardized)
  # return(list(
  #   original = cov_df_original,
  #   standardized = cov_df_standardized,
  #   pooled_sds = cov_sds
  # ))
}
  
#' Execute a search function dynamically
run_search <- function(search_name, default_params, override_params, optimization_functions) {
  selected_func <- optimization_functions[[search_name]]
  
  # Remove NULL placeholders from defaults
  valid_defaults <- default_params[!sapply(default_params, is.null)]
  
  # Merge: overrides take precedence
  run_args <- modifyList(valid_defaults, override_params)
  
  # Execute
  do.call(what = selected_func, args = run_args)
}

process_and_report_search <- function(ls_res, search_name, initial_cost, initial_intervals, initial_k) {
  eps <- 1e-12
  cat(paste0("\n--- **", toupper(search_name), " Report** ---\n"))
  
  if (!is.list(ls_res) || is.null(ls_res$cost)) {
    cat("ERROR: Result structure is invalid.\n")
    return(list(intervals = initial_intervals, k = initial_k, cost = initial_cost))
  }
  
  cat(sprintf("Original DP Cost: %.4f | New Cost: %.4f\n", initial_cost, ls_res$cost))
  
  if (ls_res$has_improved) {
    cat("✓ **Decision:** Using improved intervals from", search_name, "\n")
    print("New intervals by are:\n")
    print(ls_res$intervals)
  } else {
    cat("✗ **Decision:** No improvement; retaining DP solution.\n")
  }
  
  return(list(intervals = ls_res$intervals, k = ls_res$k_actual, cost = ls_res$cost, status = ls_res$has_improved, name = search_name))
}

perform_osip_matching <- function(data_subset, final_intervals,
                                           params, data_config, treatment_col, distance_metric) {

  cat("\n--- Performing Final Bin-Exact Matching ---\n")
  #treatment_col <- data_config$TREATMENT_VAR

  # 1. Create PS Bins

  # 1. Initialize the bin column
  data_subset$ps_bin_final <- NA

  # Extract the list of IDs from the attributes of final_intervals
  # These are the attributes you showed in the screenshot
  t_ids_list <- attr(final_intervals, "treatment_unit_ids")
  c_ids_list <- attr(final_intervals, "control_unit_ids")

  # 2. Map IDs to Bins for both Treated and Control
  for (bin_idx in seq_along(t_ids_list)) {
    bin_label <- paste0("Bin", bin_idx)

    # Assign bin to treated units in this interval
    t_ids <- t_ids_list[[bin_idx]]
    data_subset$ps_bin_final[data_subset$id %in% t_ids] <- bin_label

    # Assign bin to control units in this interval
    c_ids <- c_ids_list[[bin_idx]]
    data_subset$ps_bin_final[data_subset$id %in% c_ids] <- bin_label
  }

  # 3. Clean up orphans (if any units weren't captured by the DP intervals)
  data_subset <- data_subset[!is.na(data_subset$ps_bin_final), ]

  # 2. Prepare Data and Run MatchIt
  data_subset[[treatment_col]] <- factor(data_subset[[treatment_col]], levels = c(0, 1))
  matching_formula <- as.formula(paste(treatment_col, "~", paste(data_config$ALL_COVARIATES, collapse = " + ")))

  m_out <- matchit(
    matching_formula,
    data = data_subset,
    method = "optimal",
    distance = distance_metric,
    exact = "ps_bin_final",
    ratio = 1
  )

  # --- NEW EXTRACTION LOGIC STARTS HERE ---
  # This section creates the 'match_map' that was missing/empty
  mm <- m_out$match.matrix

  # Build map with character conversion to handle LINDER (ints) and Mixtape (chars)
  match_map_osip <- data.frame(
    treated_id = as.character(rownames(mm)),
    control_id = as.character(mm[, 1]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(control_id))

  # Create a lookup for PS values using character keys
  ps_lookup <- data_subset$ps
  names(ps_lookup) <- as.character(data_subset$id)

  # Attach PS and calculate cost (distance)
  match_map_osip <- match_map_osip %>%
    mutate(
      ps_t = ps_lookup[treated_id],
      ps_c = ps_lookup[control_id],
      cost = abs(ps_t - ps_c)
    )
  # --- NEW EXTRACTION LOGIC ENDS HERE ---

  # 4. Extract matched data for ATE calculations
  data_matched <- match.data(m_out)

  return(list(
    m_out = m_out,
    data_matched = data_matched,
    match_map = match_map_osip, # Returning this ensures plots have data
    method_used = "Optimal (Exact Bin-Matching)"
  ))
}

apply_sampling <- function(data_full, dataset_name, data_config, datasets, params, treatment_col) {
  set.seed(25)
  # treat_col <- if(dataset_name == datasets$RHC) "swang1" else "treat"
  # treatment_col <- data_config$TREATMENT_VAR
  
  # Identify IDs
  treated_ids <- data_full$id[data_full[[treatment_col]] == 1]
  control_ids <- data_full$id[data_full[[treatment_col]] == 0]
  
  # Sample Treated
  n_treated_sample <- min(params$SAMPLE_SIZE, length(treated_ids))
  sampled_treated_ids <- sample(treated_ids, n_treated_sample, replace = FALSE)
  
  # Sample Control (RHC only)
  if (dataset_name == datasets$RHC) {
    n_control_sample <- min(params$N_CONTROL_SAMPLE, length(control_ids))
    control_ids <- sample(control_ids, n_control_sample, replace = FALSE)
  }
  
  # Combine and Filter
  subset_ids <- c(sampled_treated_ids, control_ids)
  data_subset <- data_full[data_full$id %in% subset_ids, ]
  
  cat(sprintf("INFO: PS calculated on %d units. Sampled down to %d units for osip.\n", 
              nrow(data_full), nrow(data_subset)))
  
  return(data_subset)
}

ensure_numeric_col <- function(df, col_name) {
  # Extract the column data
  column_data <- df[[col_name]]
  
  # Return as-is if already numeric
  if (is.numeric(column_data)) {
    return(column_data)
  }
  
  # Safely convert binary factors ("0", "1") to numeric (0, 1)
  # using as.character() to avoid the factor level index trap
  return(as.numeric(as.character(column_data)))
}

export_matching_report <- function(map_data, data_matched, target_dir, outfile_name) {
  # 1. Attach the bin info from the matched data to the map
  # We use the treated_id to look up which bin that pair belongs to
  report_df <- map_data %>%
    left_join(data_matched %>% 
                # Convert ID to character to ensure it matches map_data types
                mutate(id_str = as.character(id)) %>%
                dplyr::select(id_str, ps_bin_final), 
              by = c("treated_id" = "id_str"))
  
  full_path <- file.path(target_dir, outfile_name)
  sink(full_path)
  
  cat("====================================================\n")
  cat("          OSIP MATCHING QUALITY REPORT             \n")
  cat("====================================================\n\n")
  
  # --- SUMMARY TABLE ---
  cat("SUMMARY TABLE: COST PER BIN\n")
  cat(sprintf("%-15s | %-8s | %-10s | %-10s\n", "Bin", "Pairs", "Mean Cost", "Max Cost"))
  cat("----------------------------------------------------\n")
  
  summary_stats <- report_df %>%
    group_by(ps_bin_final) %>%
    summarize(n = n(), 
              mean_c = mean(cost, na.rm=TRUE), 
              max_c = max(cost, na.rm=TRUE), 
              .groups = 'drop') %>%
    arrange(ps_bin_final)
  
  for(i in 1:nrow(summary_stats)) {
    cat(sprintf("%-15s | %-8d | %-10.4f | %-10.4f\n", 
                as.character(summary_stats$ps_bin_final[i]), 
                summary_stats$n[i], 
                summary_stats$mean_c[i], 
                summary_stats$max_c[i]))
  }
  cat("\n\n")
  
  # --- DETAILED PER-BIN LOG ---
  bins <- sort(unique(report_df$ps_bin_final))
  for (b in bins) {
    
    bin_data <- report_df %>% filter(ps_bin_final == b)
   
    cat(sprintf("BIN DETAIL: %s\n", b))
    cat("----------------------------------------------------\n")
    
    worst <- bin_data %>% arrange(desc(cost)) %>% head(5)
    cat("Top Max Distances:\n")
    for(j in 1:nrow(worst)) {
      cat(sprintf("  [%d] Cost: %.4f | IDs: %s <-> %s\n", 
                  j, worst$cost[j], worst$treated_id[j], worst$control_id[j]))
    }
    
    cat("\nFull Pair List (Treated -> Control):\n")
    for(k in 1:nrow(bin_data)) {
      cat(sprintf("  %s -> %s (%.3f)\n", 
                  bin_data$treated_id[k], bin_data$control_id[k], bin_data$cost[k]))
    }
    cat("\n\n")
  }
  
  sink()
  message("Report exported to: ", outfile_name)
}

#' Count units in interval using Step 1 DP boundary rules
#' Rule: (v_curr, v_next] - left-open, right-closed
#' This matches the DP's inline logic exactly
# count_units_dp_style <- function(ps_values, v_curr, v_next, eps = 1e-12) {
#   # Exact logic from Step 1 DP
#   mask <- (ps_values > v_curr) & (ps_values <= v_next + eps)
#   return(mask)
# }

# For single values (DP inner loop)
# is_in_interval_dp_style <- function(ps_value, v_curr, v_next, eps = 1e-12) {
#   return(ps_value > v_curr && ps_value <= v_next + eps)
# }

save_partition_csv <- function(df,
                               delta_val,
                               base_dir,
                               step,                 # 1 = DP, 2 = Heuristic
                               heuristic_name = NULL,   # "LS", "ELS", or "SA" (Required if step = 2)
                               step2_flag = NULL       # "strict" or "robust" (Required if step = 2)
                               #cost = NULL
                               ) {

  # 1. Validation check
  if (is.null(df) || nrow(df) == 0) {
    warning("⚠️ Cannot export CSV: Provided df is NULL or empty.")
    return(FALSE)
  }

  # 2. Clone input dataframe
  export_df <- df

  # res_list <- list()

  # 3. Round propensity score bounds to 4 decimal places (if present)
  if ("start_ps" %in% names(export_df)) export_df$start_ps <- round(export_df$start_ps, 4)
  if ("end_ps" %in% names(export_df))   export_df$end_ps   <- round(export_df$end_ps, 4)

  # 4. Remove unwanted columns
  cols_to_remove <- c("is_empty", "k_cost")
  export_df <- export_df[, !(names(export_df) %in% cols_to_remove), drop = FALSE]

  # 5. Direct column renaming
  names(export_df)[names(export_df) == "start_ps"]  <- "start"
  names(export_df)[names(export_df) == "end_ps"]    <- "end"
  names(export_df)[names(export_df) == "n_treated"] <- "treated"
  names(export_df)[names(export_df) == "n_control"] <- "control"

  # # 🆕 ADD COST COLUMN (If cost is provided)
  # if (!is.null(cost)) {
  #   export_df$cost <- round(cost, 4)
  # }

  # 6. Construct dynamic filename based on Step logic
  if (step == 1) {
    # Step 1 is always DP
    filename <- sprintf("step1_dp_partition_delta_%.2f.csv", delta_val)

  } else if (step == 2) {
    # Step 2 validation for Heuristic Name & Mode
    if (is.null(heuristic_name)) {
      warning("⚠️ Step 2 requires 'heuristic_name' to be defibed. Defaulting to 'LS'.")
      heuristic_name <- "local_search"
    }

    step2_mode <- if (step2_flag) "robust" else "strict"

    if (is.null(step2_flag) || !(step2_mode %in% c("strict", "robust"))) {
      warning("⚠️ Step 2 requires 'step2_mode' to be 'strict' or 'robust'. Defaulting to 'strict'.")
      step2_mode <- "strict"
    }

    filename <- sprintf("step2_heuristic_%s_%s_partition_delta_%.2f.csv",
                        heuristic_name, step2_mode, delta_val)
  } else {
    stop("❌ Invalid step provided. 'step' must be 1 or 2.")
  }

  csv_file_path <- file.path(base_dir, filename)
  write.csv(export_df, file = csv_file_path, row.names = FALSE, quote = FALSE)
  cat(sprintf("✅ Partition CSV saved to > %s\n\n", csv_file_path))
  return(TRUE)
}

get_max_k <- function(delta, p0, pk) {
  ceiling_val <- ceiling((pk - p0) / delta)
  k_max <- 2 * ceiling_val
  return(k_max)
}
