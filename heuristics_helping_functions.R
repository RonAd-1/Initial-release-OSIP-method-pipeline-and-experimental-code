# --- heuristics_helping_functions.R ---


source("original_distance_for_intervals.R")

# -----------------------------------------------------------------------------
# Helper: Random Perturbation of Boundaries
# -----------------------------------------------------------------------------
perturb_random <- function(bounds, n, K_max, strength = 5) {
  # bounds: current boundary vector (sorted)
  # n: max index (length of P_sorted)
  # K_max: maximum allowed intervals
  # strength: how aggressive the perturbation (1=mild, 10=aggressive)
  
  move_type <- sample(c("move", "split", "merge"), 1, 
                      prob = c(0.5, 0.3, 0.2))
  
  new_bounds <- bounds
  
  if (move_type == "move" && length(bounds) > 2) {
    # Move a random interior boundary
    b_idx <- sample(2:(length(bounds)-1), 1)
    shift <- sample(-strength:strength, 1)
    new_bounds[b_idx] <- bounds[b_idx] + shift
    
    # Clip to valid range and ensure sorted
    new_bounds[b_idx] <- max(bounds[b_idx-1] + 1, 
                             min(bounds[b_idx+1] - 1, new_bounds[b_idx]))
    
  } else if (move_type == "split" && length(bounds) < K_max + 1) {
    # Add a new boundary in a random interval
    interval_idx <- sample(1:(length(bounds)-1), 1)
    left <- bounds[interval_idx]
    right <- bounds[interval_idx + 1]
    
    if (right - left > 1) {
      # Random position in interval
      new_pos <- sample((left+1):(right-1), 1)
      new_bounds <- sort(c(bounds, new_pos))
    }
    
  } else if (move_type == "merge" && length(bounds) > 2) {
    # Remove a random interior boundary
    b_idx <- sample(2:(length(bounds)-1), 1)
    new_bounds <- bounds[-b_idx]
  }
  
  # Ensure bounds are valid
  new_bounds <- unique(sort(new_bounds))
  new_bounds <- new_bounds[new_bounds >= 0 & new_bounds <= n]
  
  # Must have at least 2 bounds (one interval)
  if (length(new_bounds) < 2) {
    return(bounds)  # Return original if invalid
  }
  
  return(new_bounds)
}

# -----------------------------------------------------------------------------
# Main Simulated Annealing Function (Updated Version)
# -----------------------------------------------------------------------------
simulated_annealing_intervals <- function(
    initial_intervals,
    cov_df,
    p_sorted,
    delta_dp,
    K_max,
    metric,
    
    # 🛑 CRITICAL NEW PARAMETERS (For Cost Function Integration) 🛑
    n_r, # The budget for redundant treated units (slack)
    fixed_redundant_unit_ids, # The fixed set of redundant units (if pre-determined)
    id_var, 
    treatment_col, 
    # data_config,
    initial_bounds,           # 🛑 NEW: Pass bounds directly from DP
    initial_cost,             # 🛑 NEW: Pass cost from DP (was orig_cost_step1)
    initial_k,                # 🛑 NEW: Pass k from DP
    
    left_border,
    right_border,

    # SA Algorithm Parameters
    temp_init = 50,
    cooling_rate = 0.995,
    max_iter = 20000,
    min_temp = 0.01,
    ... # Standard practice to keep '...' for future flexibility
) {
  
  n_ps <- length(p_sorted)
  cov_cols <- data_config$NUMERIC_COVARIATES
  
  
  # Ensure n_r is a required parameter for the cost calculation wrapper
  if (missing(n_r)) {
    stop("Error: 'n_r' (slack budget) must be provided to the Simulated Annealing function.")
  }
  
  # Initial cost
  # --- 1. Initialization (Already Done by DP) ---
  current_bounds <- initial_bounds  # Directly from DP
  current_cost <- initial_cost      # Directly from DP
  current_k <- initial_k            # Directly from DP
  current_intervals <- initial_intervals
  
  best_cost <- current_cost
  best_bounds <- current_bounds
  best_k <- current_k
  best_intervals <- current_intervals
  
  cat(sprintf("\nSimulated Annealing Initialized: Cost = %.6f, K = %d\n", 
              initial_cost, initial_k))
  
  if (!is.finite(initial_cost)) {
    stop(sprintf("ERROR: Initial cost calculation failed (Cost=%s).", initial_cost))
  }

  # --- 2. Setup Centralized Cost Caller ---
  cost_params <- list(
    p_sorted = p_sorted,
    cov_df = cov_df,
    delta_dp = delta_dp,
    K_max = K_max,
    metric = metric,
    cov_cols = cov_cols,
    n_r = n_r,
    fixed_redundant_unit_ids = fixed_redundant_unit_ids,
    id_var = id_var,
    S_inv = S_inv,
    treatment_col = treatment_col,
    data_config = data_config,
    left_border = left_border,
    right_border = right_border

  )

  call_ls_cost <- function(bounds) {
    do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
  }
  
  temp <- temp_init
  accept_count <- 0
  improve_count <- 0
  
  # Main SA loop
  for (iter in 1:max_iter) {
    
    # Generate neighbor
    # Note: Added K_max constraint to perturb_random if it doesn't already exist
    neighbor_bounds <- perturb_random(current_bounds, n_ps, K_max, 
                                      strength = max(1, ceiling(temp/10)))
    
    # 🛑 UPDATE: Use the centralized cost calculation helper 🛑
    neighbor_result <- call_ls_cost(neighbor_bounds)
    neighbor_cost <- neighbor_result$cost
    current_intervals <- neighbor_result$intervals
    
    if (!is.finite(neighbor_cost)) {
      # Infeasible neighbor, skip (Cooling still occurs)
      temp <- temp * cooling_rate
      next
    }
    
    # Acceptance criterion
    delta_E <- neighbor_cost - current_cost # Renamed delta to delta_E to avoid collision with argument 'delta'
    
    if (delta_E < 0) {
      # Always accept improvement
      current_bounds <- neighbor_bounds
      current_cost <- neighbor_cost
      accept_count <- accept_count + 1
      improve_count <- improve_count + 1
      
      # Update best
      if (current_cost < best_cost) {
        best_bounds <- current_bounds
        best_cost <- current_cost
        best_k <- neighbor_result$k_actual
        best_intervals <- current_intervals
        
        cat(sprintf("  ✓ Iter %d: NEW BEST = %.6f (K=%d, temp=%.2f)\n", 
                    iter, best_cost, best_k, temp))
      }
      
    } else {
      # Accept worse solution with probability exp(-delta_E/temp)
      accept_prob <- exp(-delta_E / temp)
      
      if (runif(1) < accept_prob) {
        current_bounds <- neighbor_bounds
        current_cost <- neighbor_cost
        accept_count <- accept_count + 1
      }
    }
    
    # Cool down
    temp <- temp * cooling_rate
    
    # Stop if temperature too low
    if (temp < min_temp) {
      cat(sprintf("Temperature below minimum (%.4f), stopping early at iter %d\n", 
                  min_temp, iter))
      break
    }
    
    # Progress reporting
    if (iter %% 1000 == 0) {
      accept_rate <- accept_count / 1000
      improve_rate <- improve_count / 1000
      cat(sprintf("Iter %d: temp=%.2f, accept=%.1f%%, improve=%.1f%%, best=%.6f\n", 
                  iter, temp, accept_rate*100, improve_rate*100, best_cost))
      accept_count <- 0
      improve_count <- 0
    }
  }
  
  # Build final intervals
  # final_intervals <- build_intervals_from_bounds(best_bounds, p_sorted, cov_df, treatment_col, id_var, delta_dp)
  
  cat(sprintf("\n✓ SA Complete: Final cost = %.6f, K = %d, Iterations = %d\n", 
              best_cost, best_k, iter))
  
  if (best_cost == initial_cost) {
    return(list(
      intervals = initial_intervals,
      cost = initial_cost,
      boundaries = initial_bounds,
      n_intervals = length(initial_bounds) - 1,
      k_actual = initial_k,
      has_improved = FALSE,
      name = "simulated_annealing"
    ))
  }
  
  return(list(
    intervals = best_intervals,
    cost = best_cost,
    boundaries = best_bounds,
    n_intervals = length(best_bounds) - 1,
    k_actual = best_k,
    has_improved = TRUE,
    name = "simulated_annealing"
  ))
}

# original_distance_for_intervals <- function(intervals, cov_df, metric, n_r, 
#                                             fixed_redundant_unit_ids, id_var, 
#                                             treatment_col, data_config) {
#   total_cost <- 0
#   n_t_non_redundant <- 0
#   cov_cols <- data_config$NUMERIC_COVARIATES
#   
#   # 1. Global Pre-check
#   if (nrow(cov_df[cov_df[[treatment_col]] == 0, ]) < 2) {
#     return(list(cost = Inf, n_t_matched = 0))
#   }
#   
#   # 2. Pre-compute S_inv (Mahalanobis Space)
#   if (metric %in% c("mahalanobis", "robust_mahalanobis")) {
#     X_for_S <- as.matrix(cov_df[, cov_cols])
#     S_inv <- solve(cov(X_for_S) + diag(1e-7, length(cov_cols)))
#     if (is.null(S_inv)) stop(paste(metric, "matrix inversion failed."))
#   }
#   
#   # 3. Loop through Intervals
#   for (r in seq_len(nrow(intervals))) {
#     vprev <- intervals$start_ps[r]
#     vcurr <- intervals$end_ps[r]
#     is_empty_flag <- intervals$is_empty[r]
#     
#     # Extract units using DP-style boundary logic
#     bin_mask <- count_units_dp_style(cov_df$ps, vprev, vcurr)
#     bin_data <- cov_df[bin_mask, , drop = FALSE]
#     treated_in_bin <- bin_data[bin_data[[treatment_col]] == 1, , drop = FALSE]
#     control_in_bin <- bin_data[bin_data[[treatment_col]] == 0, , drop = FALSE]
#     
#     # 🆕 VALIDATION: Check counts match what DP stored
#     expected_n_treated <- intervals$n_treated[r]
#     expected_n_control <- intervals$n_control[r]
#     actual_n_treated <- nrow(treated_in_bin)
#     actual_n_control <- nrow(control_in_bin)
#     
#     if (actual_n_treated != expected_n_treated || actual_n_control != expected_n_control) {
#       warning(sprintf(
#         "Bin %d count mismatch! Expected: T=%d, C=%d | Actual: T=%d, C=%d | PS=[%.6f, %.6f], Empty=%s",
#         r, expected_n_treated, expected_n_control, 
#         actual_n_treated, actual_n_control,
#         vprev, vcurr, is_empty_flag
#       ))
#       # Optionally: return Inf or stop execution
#       # return(list(cost = Inf, n_t_matched = n_t_non_redundant))
#     }
#     
#     if (nrow(treated_in_bin) == 0) next
#     
#     # Identify non-redundant units
#     is_joker <- treated_in_bin[[id_var]] %in% fixed_redundant_unit_ids
#     X_t_matching <- as.matrix(treated_in_bin[!is_joker, cov_cols, drop = FALSE])
#     X_c <- as.matrix(control_in_bin[, cov_cols, drop = FALSE])
#     
#     n_t_non_redundant <- n_t_non_redundant + nrow(X_t_matching)
#     if (nrow(X_t_matching) == 0) next
#     
#     # Local Infeasibility Check
#     if (nrow(X_c) == 0) return(list(cost = Inf, n_t_matched = n_t_non_redundant))
#     
#     # --- Vectorized Mahalanobis ---
#     t_diag <- rowSums((X_t_matching %*% S_inv) * X_t_matching)
#     c_diag <- rowSums((X_c %*% S_inv) * X_c)
#     cross_term <- (X_t_matching %*% S_inv) %*% t(X_c)
#     dist_matrix <- outer(t_diag, c_diag, "+") - 2 * cross_term
#     
#     # Row-wise minimum distance
#     min_distances <- apply(dist_matrix, 1, function(x) min(pmax(0, x), na.rm = TRUE))
#     total_cost <- total_cost + sum(min_distances)
#   }
#   
#   # 4. Final Match Count Validation
#   expected_t <- sum(intervals$n_treated) 
#   if (n_t_non_redundant < expected_t) {
#     return(list(cost = Inf, n_t_matched = n_t_non_redundant))
#   }
#   
#   return(list(cost = total_cost, n_t_matched = n_t_non_redundant))
# }

# original_distance_for_intervals <- function(intervals, cov_df, metric, n_r, 
#                                             fixed_redundant_unit_ids, id_var, 
#                                             treatment_col, data_config) {
#   total_cost <- 0
#   n_t_non_redundant <- 0
#   eps_fuzzy <- 1e-10 
#   cov_cols <- data_config$NUMERIC_COVARIATES
#   
#   # 1. Global Pre-check
#   if (nrow(cov_df[cov_df[[treatment_col]] == 0, ]) < 2) {
#     return(list(cost = Inf, n_t_matched = 0))
#   }
#   
#   # 2. Pre-compute S_inv (Mahalanobis Space)
#   if (metric %in% c("mahalanobis", "robust_mahalanobis")) {
#     X_for_S <- as.matrix(cov_df[, cov_cols])
#     S_inv <- solve(cov(X_for_S) + diag(1e-7, length(cov_cols)))
#     if (is.null(S_inv)) stop(paste(metric, "matrix inversion failed."))
#   }
#   
#   # 3. Loop through Intervals
#   for (r in seq_len(nrow(intervals))) {
#     vprev <- intervals$start_ps[r]
#     vcurr <- intervals$end_ps[r]
#     is_empty_flag <- intervals$is_empty[r]
#     
#     # --- The "Boundary Clash" Guard ---
#     if (is_empty_flag) {
#       # Strict window for empty bridges so they don't "steal" boundary units
#       bin_mask <- cov_df$ps > (vprev + eps_fuzzy) & cov_df$ps < (vcurr - eps_fuzzy)
#     } else {
#       # Inclusive window for matching bins to ensure boundary units are caught
#       bin_mask <- cov_df$ps > (vprev - eps_fuzzy) & cov_df$ps <= (vcurr + eps_fuzzy)
#     }
#     
#     bin_data <- cov_df[bin_mask, , drop = FALSE]
#     treated_in_bin <- bin_data[bin_data[[treatment_col]] == 1, , drop = FALSE]
#     control_in_bin <- bin_data[bin_data[[treatment_col]] == 0, , drop = FALSE]
#     
#     if (nrow(treated_in_bin) == 0) next
#     
#     # Identify non-redundant units
#     is_joker <- treated_in_bin[[id_var]] %in% fixed_redundant_unit_ids
#     X_t_matching <- as.matrix(treated_in_bin[!is_joker, cov_cols, drop = FALSE])
#     X_c <- as.matrix(control_in_bin[, cov_cols, drop = FALSE])
#     
#     n_t_non_redundant <- n_t_non_redundant + nrow(X_t_matching)
#     if (nrow(X_t_matching) == 0) next
#     
#     # Local Infeasibility Check
#     if (nrow(X_c) == 0) return(list(cost = Inf, n_t_matched = n_t_non_redundant))
#     
#     # --- Vectorized Mahalanobis (Fixed Syntax) ---
#     t_diag <- rowSums((X_t_matching %*% S_inv) * X_t_matching)
#     c_diag <- rowSums((X_c %*% S_inv) * X_c)
#     
#     # 🛑 CORRECTED: Use standard matrix multiplication %*% and t()
#     cross_term <- (X_t_matching %*% S_inv) %*% t(X_c)
#     
#     dist_matrix <- outer(t_diag, c_diag, "+") - 2 * cross_term
#     
#     # Row-wise minimum distance
#     min_distances <- apply(dist_matrix, 1, function(x) min(pmax(0, x), na.rm = TRUE))
#     total_cost <- total_cost + sum(min_distances)
#   }
#   
#   # 4. Final Match Count Validation
#   expected_t <- sum(intervals$n_treated) 
#   if (n_t_non_redundant < expected_t) {
#     # If we missed units due to boundary logic, flag as Inf
#     return(list(cost = Inf, n_t_matched = n_t_non_redundant))
#   }
#   
#   return(list(cost = total_cost, n_t_matched = n_t_non_redundant))
# }

# original_distance_for_intervals <- function(
#     intervals, 
#     cov_df, 
#     metric, 
#     n_r, # The resource budget for controls/slack (kept for signature completeness)
#     fixed_redundant_unit_ids, 
#     id_var,          
#     treatment_col,   
#     data_config
# ) {
#   # The functions build_intervals_from_bounds, check_feasible, and 
#   # original_distance_for_intervals MUST be available in the environment.
#   
#   # 1. Setup and Pre-calculations
#   # We assume 'ps' is a column in cov_df
#   
#   # 1. Setup and Pre-calculations
#   # TODO: Generalize this, this is specific for nsw_mixtape
#   # cov_names <- setdiff(names(cov_df), c(id_var, "ps", treatment_col, "data_id"))
#   
#   # Debug
#   # browser()
#   
#   total_cost <- 0
#   # 🛑 CRITICAL FIX: Initialize counter for units contributing to the cost 🛑
#   n_t_non_redundant <- 0
#   eps_fuzzy <- 1e-10
#   cov_cols <- data_config$NUMERIC_COVARIATES
#   
#   # At the very start of original_distance_for_intervals:
#   control_rows_global <- cov_df[cov_df[[treatment_col]] == 0, cov_cols, drop = FALSE]
#   
#   if (nrow(control_rows_global) < 2) {
#     # If there aren't at least 2 controls in the WHOLE dataset, 
#     # we can't define a Mahalanobis space.
#     return(list(cost = Inf, n_t_matched = 0))
#   }
#   
#   if (length(cov_cols) == 0) stop("No covariates present in cov_df for computing original distance.")
#   
#   if (metric %in% c("mahalanobis", "robust_mahalanobis")) {
#     # Get control rows for the covariance structure
#     
#     # --- PRE-COMPUTE S_inv ONCE ---{
#       X_for_S <- as.matrix(cov_df[, cov_cols])
#       S_inv <- solve(cov(X_for_S) + diag(1e-7, length(cov_cols)))
#       if (is.null(S_inv)) stop(paste(metric, "matrix inversion failed."))
#     }
#     
#     # Debug
#     # browser()
#     
#   #   if (metric == "robust_mahalanobis") {
#   #     # Use Minimum Covariance Determinant (MCD) for robustness
#   #     # nsamp = "mcd" is the standard high-breakdown estimator
#   #     robust_S <- tryCatch({
#   #       suppressWarnings(MASS::cov.mcd(control_rows)$cov)
#   #     }, error = function(e) {
#   #       cov(control_rows, use = "pairwise.complete.obs")
#   #     })
#   #     S <- robust_S
#   #   } else if (metric == "mahalanobis") {
#   #     # Standard Mahalanobis
#   #     S <- cov(control_rows, use = "pairwise.complete.obs")
#   #   }
#   #   
#   #   # Inversion with regularization for numerical stability
#   #   S_inv <- tryCatch(
#   #     solve(S), 
#   #     error = function(e) solve(S + diag(1e-6, ncol(S)))
#   #   )
#   #   if (is.null(S_inv)) stop(paste(metric, "matrix inversion failed."))
#   # }
#   # 
#   
#   
#   # Precompute Mahalanobis inverse on controls if requested
#   # if (metric == "mahalanobis") {
#   #   control_rows <- cov_df[cov_df[[treatment_col]] == 0, cov_names, drop = FALSE]
#   #   # Use only complete observations for covariance matrix estimation
#   #   S <- cov(control_rows, use = "pairwise.complete.obs")
#   #   # Add a small value to the diagonal for regularization
#   #   S_inv <- tryCatch(solve(S), error = function(e) solve(S + diag(1e-6, ncol(S))))
#   #   if (is.null(S_inv)) stop("Mahalanobis matrix inversion failed.")
#   # }
#   
#   # 2. Iterate through all intervals
#   for (r in seq_len(nrow(intervals))) {
#     vprev <- intervals$start_ps[r]  
#     vcurr <- intervals$end_ps[r]
#     is_empty_flag <- intervals$is_empty[r]
#     
#     # Replace the is.finite check with this:
#     if (!is.numeric(vprev) || !is.numeric(vcurr)) {
#       stop(sprintf("Non-numeric interval bounds detected during distance calc: [%s, %s]", 
#                    as.character(vprev), as.character(vcurr)))
#     }
#     # 🛑 THE "BOUNDARY CLASH" FIX 🛑
#     if (is_empty_flag) {
#       # If the bin is a bridge, use a STRICTOR window so it doesn't "steal" 
#       # the treated unit sitting on the border.
#       bin_mask <- cov_df$ps > (vprev + 1e-9) & cov_df$ps < (vcurr - 1e-9)
#     } else {
#       # If the bin is for matching, use a WIDER window to ensure 
#       # boundary units are captured.
#       bin_mask <- cov_df$ps > (vprev - 1e-9) & cov_df$ps <= (vcurr + 1e-9)
#     }
#     
#     bin_data <- cov_df[bin_mask, , drop = FALSE]
#     
#     # bin_mask <- cov_df$ps > (vprev - eps_fuzzy) & cov_df$ps <= (vcurr + eps_fuzzy)
#     # bin_data <- cov_df[bin_mask, , drop = FALSE]
#     
#     # 1. SKIP EMPTY INTERVALS IMMEDIATELY
#     # If the DP labeled it empty, or there are no treated units, cost is 0
# 
#     treated_in_bin <- bin_data[bin_data[[treatment_col]] == 1, , drop = FALSE]
#     control_in_bin <- bin_data[bin_data[[treatment_col]] == 0, , drop = FALSE]
#     
#     message(sprintf("Processing interval %d: [%.4f, %.4f] - Treated: %d, Controls: %d", 
#                     r, vprev, vcurr, nrow(treated_in_bin), nrow(control_in_bin)))
#     
#     if (nrow(treated_in_bin) == 0) next # Cost remains 0 for this interval
# 
#     # Identify non-redundant units (Units not in the Joker/Slack list)
#     is_joker <- treated_in_bin[[id_var]] %in% fixed_redundant_unit_ids
#     X_t_matching <- as.matrix(treated_in_bin[!is_joker, cov_cols, drop = FALSE])
#     X_c <- as.matrix(control_in_bin[, cov_cols, drop = FALSE])
#     
#     # Update counter for total non-redundant units processed
#     n_t_non_redundant <- n_t_non_redundant + nrow(X_t_matching)
#     
#     # Skip math if there are no units to match in this bin
#     if (nrow(X_t_matching) == 0) next
#     
#     # 🛑 INFEASIBILITY CHECK 🛑
#     # If we have treated units but ZERO controls, the cost is Infinite
#     if (nrow(X_c) == 0) {
#       return(list(cost = Inf, n_t_matched = n_t_non_redundant))
#     }
#     
#     # 🛑 VECTORIZED MAHALANOBIS CALCULATION 🛑
#     # Formula: D^2 = diag(Xt S_inv Xt') + diag(Xc S_inv Xc') - 2(Xt S_inv Xc')
#     
#     # 1. Quadratic terms for rows
#     t_diag <- rowSums((X_t_matching %*% S_inv) * X_t_matching)
#     c_diag <- rowSums((X_c %*% S_inv) * X_c)
#     
#     # 2. Cross term
#     cross_term <- X_t_matching %*% S_inv %t% X_c
#     
#     # 3. Full distance matrix (Treated rows x Control columns)
#     # outer() creates the base matrix, then subtract 2 * cross_term
#     dist_matrix <- outer(t_diag, c_diag, "+") - 2 * cross_term
#     
#     # 4. Find minimum distance for each treated unit
#     # (Using pmax(0, ...) to handle tiny negative numbers from float precision)
#     min_distances <- apply(dist_matrix, 1, function(x) min(pmax(0, x), na.rm = TRUE))
#     
#     # Add to total
#     total_cost <- total_cost + sum(min_distances)
#   }
#   
#   # 3. Final Validation
#   # Ensure we didn't 'lose' any units. Total treated in non-empty bins 
#   # should match the units we actually processed.
#   expected_t <- sum(intervals$n_treated)
#   
#   if (n_t_non_redundant == 0 && expected_t > 0) {
#     return(list(cost = Inf, n_t_matched = 0))
#   }
#   
#   return(list(
#     cost = total_cost,
#     n_t_matched = n_t_non_redundant
#   ))
# }
    
    # Safe bounds handling
    # if (!is.finite(vprev)) vprev <- 0
    # if (!is.finite(vcurr)) vcurr <- 1
    
    # Subset all units in the bin (including ID column)
    
    # bin_data <- cov_df[cov_df$ps > vprev & cov_df$ps <= vcurr, , drop = FALSE]
    
    # Debug 
    # browser()
    # 
    # # Immediate infeasibility check (treated units with no controls)
    # # 🛑 CRITICAL: Check if this interval is the culprit
    # if (nrow(treated_in_bin) > 0 && nrow(control_in_bin) == 0) {
    #   # If this treated unit is NOT redundant, we have no one to match it to.
    #   non_redundant_treated <- setdiff(treated_in_bin[[id_var]], fixed_redundant_unit_ids)
    #   if (length(non_redundant_treated) > 0) {
    #     message(sprintf("Bin [%.4f, %.4f] has no controls for non-redundant treated units.", vprev, vcurr))
    #     return(list(cost = Inf, redundant_unit_ids = fixed_redundant_unit_ids))
    #   }
    # }
    
    # Debug 
  #   browser()
  #   
  #   # 3. Calculate cost for all treated units in bin
  #   if (nrow(treated_in_bin) > 0) {
  #     X_t_all <- as.matrix(treated_in_bin[, cov_names, drop = FALSE])
  #     X_c <- as.matrix(control_in_bin[, cov_names, drop = FALSE])
  #     
  #     for (i in seq_len(nrow(X_t_all))) {
  #       current_id <- treated_in_bin[[id_var]][i] # Get the unique ID for the unit
  #       
  #       # CHECK FOR "JOKER" STATUS (UNIT IS REDUNDANT)
  #       if (current_id %in% fixed_redundant_unit_ids) {
  #         individual_cost <- 0
  #       } else {
  #         # 🛑 UNIT IS NON-REDUNDANT (REGULAR MATCHING) 🛑
  #         n_t_non_redundant <- n_t_non_redundant + 1 # 🛑 INCREMENT COUNTER 🛑
  #         
  #         if (metric == "squared_euclidean") {
  #           diffs <- sweep(X_c, 2, X_t_all[i,], FUN = "-")
  #           d2 <- rowSums(diffs^2)
  #           individual_cost <- min(d2, na.rm = TRUE)
  #         } else if (metric %in% c("mahalanobis", "robust_mahalanobis")) {
  #           # Both use the same quadratic form, the difference is in S_inv
  #           diffs <- sweep(X_c, 2, X_t_all[i,], FUN = "-")
  #           # Matrix multiplication: (diffs %*% S_inv) then multiply by diffs row-wise
  #           qvals <- rowSums((diffs %*% S_inv) * diffs)
  #           individual_cost <- min(qvals, na.rm = TRUE)
  #           
  #           # Debug 
  #           browser()
  #           
  #           message(sprintf("Bin [%.4f, %.4f] has a cost of %.4f.", vprev, vcurr, individual_cost))
  #           
  #           } else {
  #               stop("Error: Undifed distance metric!!")
  #           }
  #       }
  #         
  #         # } else { # mahalanobis
  #         #   # Ensure S_inv exists if mahal is requested
  #         #   if (!exists("S_inv")) stop("Mahalanobis S_inv matrix not computed.")
  #         #   diffs <- sweep(X_c, 2, X_t_all[i,], FUN = "-")
  #         #   qvals <- rowSums((diffs %*% S_inv) * diffs)
  #         #   individual_cost <- min(qvals, na.rm = TRUE)
  #         # }
  #         
  #         # Additional check if a non-redundant unit still found no match (min is Inf)
  #         if (!is.finite(individual_cost)) {
  #           # Debug 
  #           browser()
  #           
  #           return(list(cost = Inf, redundant_unit_ids = fixed_redundant_unit_ids))
  #         }
  #       }
  #     # Debug 
  #     browser()
  #     
  #       total_cost <- total_cost + individual_cost
  #     }
  # }
  
#   cat(sprintf("DEBUG: Total Cost = %.4f, Non-Redundant Units = %d\n", total_cost, n_t_non_redundant))
#   
#   # Debug 
#   browser()
#   
#   if (n_t_non_redundant == 0 && nrow(cov_df[cov_df[[treatment_col]] == 1,]) > 0) {
#     
#     # Debug 
#     browser()
#     
#     # Only return Inf if there WERE treated units but we ignored them all 
#     # (which shouldn't happen with a proper n_r budget)
#     final_cost <- Inf 
#   } else {
#     # Return total_cost (the SUM), matching the DP's objective function
#     
#     # Debug 
#     browser()
#     
#     final_cost <- total_cost
#   }
#   
#   # 4. Final Normalization 🛑 CRITICAL FIX FOR COST=0 BUG 🛑
#   # if (n_t_non_redundant == 0) {
#   #   cat(sprintf("DEBUG: Total Cost = %.4f, Non-Redundant Units = %d\n", total_cost, n_t_non_redundant))
#   #   # If all treated units were successfully identified as redundant, 
#   #   # the solution is trivial and should be penalized, not rewarded with cost=0.
#   #   final_cost <- Inf 
#   # } else {
#   #   # Calculate the AVERAGE cost per non-redundant unit
#   #   final_cost <- total_cost / n_t_non_redundant
#   # }
#   # Debug 
#   browser()
#   
#   # 5. Return the cost and the fixed set of IDs
#   return(list(cost = final_cost, redundant_unit_ids = fixed_redundant_unit_ids))
# }

local_search_intervals <- function(
    initial_intervals,
    initial_bounds,           # 🛑 Pass bounds directly from DP
    initial_cost,             # 🛑 Pass cost from DP
    initial_k,                # 🛑 Pass k from DP
    cov_df,                   # full dataset including ps and treat
    p_sorted,                 # vector of ps sorted
    delta_dp,                 # max width allowed
    K_max,                    # max number of intervals allowed
    metric,                   # distance metric
    S_inv,                    # 🆕 Pre-computed covariance inverse
    
    # Slack parameters
    n_r,                      # The budget for redundant treated units
    fixed_redundant_unit_ids, # The fixed set of redundant units
    id_var,
    treatment_col,
    data_config,

    # LS Algorithm Parameters
    max_iter = 20000,
    max_jump = 10,
    ...
) {
  
  if (missing(n_r)) stop("Error: 'n_r' (slack budget) is required.")
  
  n_ps <- length(p_sorted)  # 🛑 Changed from 'n' to 'n_ps'
  cov_cols <- data_config$NUMERIC_COVARIATES
  
  # Debug
  # browser()
  
  # --- 1. Initialization (Already Done by DP) ---
  best_bounds <- initial_bounds
  best_cost <- initial_cost
  best_k <- initial_k
  best_intervals <- initial_intervals 
  
  cat(sprintf("\nLocal Search Initialized from DP:\n"))
  cat(sprintf("  Cost = %.6f, Boundaries = %d, K_actual = %d\n", 
              best_cost, length(best_bounds) - 1, best_k))
  
  if (!is.finite(best_cost)) {
    stop("ERROR: Initial cost from DP is not finite.")
  }
  
  # Debug
  # browser()
  
  # --- 2. Setup Centralized Cost Caller ---
  # 🆕 Use the new function for heuristic search
  # call_ls_cost <- function(bounds) {
  #   # Call the function and store the list locally first
  #   res_list <- calculate_distance_from_bounds(
  #     bounds = bounds,
  #     cov_df = cov_df,
  #     metric = metric,
  #     fixed_redundant_unit_ids = fixed_redundant_unit_ids,
  #     id_var = id_var,
  #     treatment_col = treatment_col,
  #     S_inv = S_inv,
  #     cov_cols = cov_cols,
  #     delta_dp = delta_dp,
  #     K_max = K_max,
  #     p_sorted = p_sorted,
  #     n_r = n_r
  #   )
  # }
  
  # Debug
  # browser()
  
  cost_params <- list(
    p_sorted = p_sorted,
    cov_df = cov_df,
    delta_dp = delta_dp,
    K_max = K_max,
    metric = metric,
    cov_cols = cov_cols,
    n_r = n_r,
    fixed_redundant_unit_ids = fixed_redundant_unit_ids,
    id_var = id_var,
    treatment_col = treatment_col,
    data_config = data_config,
    S_inv = S_inv
  )
  
  call_ls_cost <- function(bounds) {
    do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
  }
  
  # --- 3. Local Search Loop ---
  iter <- 1
  improved <- TRUE
  b_range <- length(best_bounds) - 1
  
  while(improved && iter <= max_iter) {
    improved <- FALSE
    
    # Debug 
    browser()
    
    # Try moving each interior boundary
    for(b in 2:b_range) {
      
      for(dir in c(-1, 1)) {
        for(jump in 1:max_jump) {
          
          # Debug
          # browser()
          
          cand_bounds <- best_bounds
          cand_bounds[b] <- cand_bounds[b] + dir*jump
          
          # Bounds checks
          if(cand_bounds[b] < 0 || cand_bounds[b] > n_ps) next  # 🛑 Use n_ps
          if(any(diff(cand_bounds) <= 0)) next
          
          
          # Debug
          # browser()
          
          # Evaluate cost
          cand_result <- call_ls_cost(cand_bounds)
          
          # Debug
          # browser()
          
          cand_cost <- cand_result$cost
          cand_k <- cand_result$k_actual
          intervals <- cand_result$intervals
          
          # browser()
          
          # Accept if better
          if(is.finite(cand_cost) && cand_cost < best_cost - 1e-9) {
            # browser()
            
            best_cost <- cand_cost
            best_bounds <- cand_bounds
            best_k <- cand_k
            best_intervals <- intervals 
            
            improved <- TRUE
            cat(sprintf("  ✓ ACCEPT (Iter %d): Cost=%.6f, K=%d\n", 
                        iter, best_cost, best_k))
            break
          }
        } 
        if(improved) break
      } 
      if(improved) break
    } 
    
    iter <- iter + 1
  }
  
  # --- 4. Final Result ---
  # Just build the intervals structure (cost already calculated)
  # Debug
  # browser()
  
  # Debug
  # browser()
  
  # final_intervals <- build_intervals_from_bounds(best_bounds, p_sorted, cov_df, treatment_col, id_var, delta_dp)
  
  cat(sprintf("\n✓ Local Search Complete: Final cost = %.6f, K = %d\n", 
              best_cost, best_k))
  
  # browser()
  
  # No updated has been made
  if (best_cost == initial_cost) {
    return(list(
      intervals = initial_intervals,
      cost = initial_cost,
      boundaries = initial_bounds,
      n_intervals = length(initial_bounds) - 1,
      k_actual = initial_k,
      has_improved = FALSE,
      name = "local_seach"
    ))
  }

  return(list(
    intervals = best_intervals,
    cost = best_cost,
    boundaries = best_bounds,
    n_intervals = length(best_bounds) - 1,
    k_actual = best_k,
    has_improved = TRUE,
    name = "local_seach"
  ))
}

# ---------------------------------------------
# Local Search Function (Updated to Full Parameter Set)
# ---------------------------------------------
# local_search_intervals <- function(
#     initial_intervals,  # data.frame(start_ps, end_ps) from DP
#     cov_df,             # full dataset including ps and treat
#     p_sorted,           # vector of ps sorted
#     delta_dp,              # max width allowed
#     K_max,              # max number of intervals allowed
#     metric,             # distance metric
#     S_inv,              # 🆕 Pre-computed covariance inverse
#     
#     # 🛑 CRITICAL FULL PARAMETER SET (NOW INCLUDED) 🛑
#     n_r,                          # The budget for redundant treated units (slack)
#     fixed_redundant_unit_ids,     # The fixed set of redundant units/Jokers
#     id_var,           # Changed to lowercase
#     treatment_col,   # Changed to lowercase
#     data_config,                   # Full configuration list
#     initial_bounds,           # 🛑 NEW: Pass bounds directly from DP
#     initial_cost,             # 🛑 NEW: Pass cost from DP (was orig_cost_step1)
#     initial_k,                # 🛑 NEW: Pass k from DP
#     
#     # LS Algorithm Parameters
#     max_iter = 20000,
#     max_jump = 10,
#     ...
# ) {
#   
#   # Simple check to ensure required parameters are present
#   if (missing(n_r)) stop("Error: 'n_r' (slack budget) is required.")
#   
#   n <- length(p_sorted)
#   
#   
#   # # --- 1. Initialization ---
#   # initial_bounds <- boundaries_from_intervals(initial_intervals, p_sorted)
#   # initial_bounds <- sort(unique(c(0, initial_bounds, n)))
#   # best_bounds <- initial_bounds
#   
#   # Convert to boundaries (assumes boundaries_from_intervals is available)
#   # initial_bounds <- boundaries_from_intervals(initial_intervals, p_sorted)
#   # initial_bounds <- sort(unique(c(0, initial_bounds, n)))
#   # best_bounds <- initial_bounds
#   
#   # --- 🛑 Setup Centralized Cost Caller (Replacing internal cost_of) 🛑 ---
#   # Store parameters common to all cost calls
#   # cost_params <- list(
#   #   p_sorted = p_sorted,
#   #   cov_df = cov_df, 
#   #   delta_dp = delta_dp, 
#   #   K_max = K_max, 
#   #   metric = metric, 
#   #   n_r = n_r,
#   #   fixed_redundant_unit_ids = fixed_redundant_unit_ids,
#   #   id_var,           # Changed to lowercase
#   #   treatment_col,   # Changed to lowercase
#   #   data_config
#   # )
#   # 
#   # # Helper function to call calculate_ls_cost with fixed parameters
#   # call_ls_cost <- function(bounds) {
#   #   do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
#   # }
#   # ----------------------------------------------------------------------
#   
#   # Initial cost
#   # --- 1. Initialization (Already Done by DP) ---
#   best_bounds <- initial_bounds  # Directly from DP
#   best_cost <- initial_cost      # Directly from DP
#   best_k <- initial_k            # Directly from DP
#   
#   cat(sprintf("\nLocal Search Initialized: Cost = %.6f, Boundaries = %d, K_actual = %d\n", 
#               best_cost, length(best_bounds) - 1, best_k))
#   
#   if (!is.finite(best_cost)) {
#     stop(sprintf("ERROR: Initial cost calculation failed (Cost=%s).", best_cost))
#   }
#   
#   # --- 2. Setup Centralized Cost Caller ---
#   cost_params <- list(
#     p_sorted = p_sorted,
#     cov_df = cov_df, 
#     delta_dp = delta_dp, 
#     K_max = K_max, 
#     metric = metric, 
#     n_r = n_r,
#     fixed_redundant_unit_ids = fixed_redundant_unit_ids,
#     id_var = id_var,
#     treatment_col = treatment_col,
#     data_config = data_config
#   )
#   
#   call_ls_cost <- function(bounds) {
#     do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
#   }
#   
#   # -------------------------------------------------
#   # Local Search Loop
#   # -------------------------------------------------
#   iter <- 1
#   improved <- TRUE
#   
#   while(improved && iter <= max_iter) {
#     improved <- FALSE
#     
#     # Try moving each interior boundary
#     for(b in 2:(length(best_bounds)-1)) {
#       
#       for(dir in c(-1, 1)) {
#         for(jump in 1:max_jump) {
#           
#           cand_bounds <- best_bounds
#           cand_bounds[b] <- cand_bounds[b] + dir*jump
#           
#           # Debug
#           if (!is.numeric(cand_bounds[b]) || !is.numeric(n)) {
#             cat("ERROR: cand_bounds[b] =", cand_bounds[b], "class =", class(cand_bounds[b]), "\n")
#             cat("       n =", n, "class =", class(n), "\n")
#             stop("Type mismatch in bounds check")
#           }
#           
#           # Bounds checks
#           if(cand_bounds[b] < 0 || cand_bounds[b] > n) next
#           if(any(diff(cand_bounds) <= 0)) next # Must remain strictly increasing
#           
#           # Evaluate cost 🛑 UPDATE: Use centralized cost function
#           cand_result <- call_ls_cost(cand_bounds)
#           cand_cost <- cand_result$cost
#           cand_k <- cand_result$k_actual
#           
#           # Accept if better
#           if(is.finite(cand_cost) && cand_cost < best_cost - 1e-9) {
#             best_cost <- cand_cost
#             best_bounds <- cand_bounds
#             best_k <- cand_k
#             improved <- TRUE
#             cat(sprintf("  ✓ ACCEPT (Iter %d): Cost=%.6f, Boundaries=%d, K_actual=%d\n", 
#                         iter, best_cost, length(best_bounds) - 1, best_k))
#             break
#           }
#         } 
#         if(improved) break
#       } 
#       if(improved) break
#     } 
#     
#     iter <- iter + 1
#   }
#   
#   # --- 3. Final Result ---
#   final_intervals <- build_intervals_from_bounds(best_bounds, p_sorted, cov_df, treatment_col)
#   
#   cat(sprintf("\n✓ Local Search Complete: Final cost = %.6f, K = %d\n", 
#               best_cost, best_k))
#   
#   return(list(
#     intervals = final_intervals,
#     cost = best_cost,
#     boundaries = best_bounds,
#     n_intervals = length(best_bounds) - 1,
#     k_actual = best_k
#   ))
# }

# cost_params <- list(
#   p_sorted = p_sorted,
#   cov_df = cov_df_standardized,
#   delta_dp = PARAMS$DELTA_DP,
#   K_max = PARAMS$K_DP,
#   metric = PARAMS$ORIGINAL_DISTANCE_METRIC,
#   cov_cols = cov_cols,
#   n_r = PARAMS$N_R,
#   fixed_redundant_unit_ids = fixed_redundant_unit_ids,
#   id_var = id_var,
#   treatment_col = treatment_col,
#   data_config = data_config,
#   S_inv = S_inv
# )
# # 
# # call_ls_cost <- function(bounds) {
# #   do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
# # }
# call_ls_cost <- function(bounds) {
#       do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
# }
  
# -----------------------------------------------------------------------------
# Enhanced Local Search with Split/Merge
# -----------------------------------------------------------------------------
enhanced_local_search <- function(
    initial_intervals,
    p_sorted, 
    cov_df, 
    delta_dp,
    K_max,
    metric,
    # 🛑 CRITICAL FIX: Add the new parameters 🛑
    fixed_redundant_unit_ids, 
    id_var,           # Changed to lowercase
    treatment_col,   # Changed to lowercase
    data_config,
    n_r, # Need n_r for the cost function (Original Distance)
    # Local Search Algorithm Parameters
    initial_bounds,           # 🛑 NEW: Pass bounds directly from DP
    initial_cost,             # 🛑 NEW: Pass cost from DP (was orig_cost_step1)
    initial_k,                # 🛑 NEW: Pass k from DP
    S_inv,
    # LS Algorithm Parameters
    max_iter = 20000,
    max_jump = 10,
    ...
) {
  
  # n <- length(p_sorted)
  
  # Convert to boundaries (assumes boundaries_from_intervals is available)
  # initial_bounds <- boundaries_from_intervals(initial_intervals, p_sorted)
  # initial_bounds <- sort(unique(c(0, initial_bounds, n)))
  # best_bounds <- initial_bounds
  # 
  # # --- 🛑 REMOVED: Internal cost_of function (now externalized) 🛑 ---
  # 
  # # Initial cost
  # # Assumes calculate_ls_cost is defined and available
  # initial_result <- calculate_ls_cost(
  #   bounds = best_bounds, 
  #   p_sorted = p_sorted,
  #   cov_df = cov_df, 
  #   delta_dp = delta_dp, 
  #   K_max = K_max, 
  #   metric = metric, 
  #   n_r = n_r,
  #   fixed_redundant_unit_ids = fixed_redundant_unit_ids,
  #   id_var,           # Changed to lowercase
  #   treatment_col,   # Changed to lowercase
  #   data_config
  # )
  # 
  # best_cost <- initial_result$cost
  # best_k <- initial_result$k_actual
  # 
  # cat(sprintf("\nEnhanced Local Search Initialized: Cost = %.6f, K = %d\n", 
  #             best_cost, best_k))
  # 
  # if (!is.finite(best_cost)) {
  #   stop("Initial solution is infeasible!")
  # }
  # 
  
  n_ps <- length(p_sorted)
  cov_cols <- data_config$NUMERIC_COVARIATES
  
  # Initial cost
  # --- 1. Initialization (Already Done by DP) ---
  best_bounds <- initial_bounds  # Directly from DP
  best_cost <- initial_cost      # Directly from DP
  best_k <- initial_k            # Directly from DP
  best_intervals <- initial_intervals
  
  cat(sprintf("\nEnhanced LS Initialized: Cost = %.6f, Boundaries = %d, K_actual = %d\n", 
              best_cost, length(best_bounds) - 1, best_k))
  
  if (!is.finite(best_cost)) {
    stop(sprintf("ERROR: Initial cost calculation failed (Cost=%s).", best_cost))
  }
  
  # --- 2. Setup Centralized Cost Caller ---

  cost_params <- list(
    p_sorted = p_sorted,
    cov_df = cov_df,
    delta_dp = delta_dp,
    K_max = K_max,
    metric = metric,
    cov_cols = cov_cols,
    n_r = n_r,
    fixed_redundant_unit_ids = fixed_redundant_unit_ids,
    id_var = id_var,
    treatment_col = treatment_col,
    data_config = data_config,
    S_inv = S_inv
  )
  
  call_ls_cost <- function(bounds) {
    do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
  }
  
  # Main loop
  iter <- 1
  improved <- TRUE
  
  while (improved && iter <= max_iter) {
    improved <- FALSE
    
    # ==============================
    # MOVE: Shift existing boundaries
    # ==============================
    if (length(best_bounds) > 2) {
      for (b in 2:(length(best_bounds)-1)) {
        for (dir in c(-1, 1)) {
          for (jump in 1:max_jump) {
            cand_bounds <- best_bounds
            cand_bounds[b] <- cand_bounds[b] + dir * jump
            
            if (cand_bounds[b] < 0 || cand_bounds[b] > n_ps) next
            if (any(diff(cand_bounds) <= 0)) next
            
            # Call the helper cost function
            cand_result <- calculate_ls_cost(
              bounds = cand_bounds, 
              p_sorted = p_sorted,
              cov_df = cov_df, 
              delta_dp = delta_dp, 
              K_max = K_max, 
              metric = metric, 
              n_r = n_r,
              cov_cols = cov_cols,
              fixed_redundant_unit_ids = fixed_redundant_unit_ids,
              id_var = id_var,           # Changed to lowercase
              treatment_col = treatment_col,   # Changed to lowercase
              data_config = data_config,
              S_inv = S_inv
            )
            
            # Debug
            # browser()
            
            if (is.finite(cand_result$cost) && 
                cand_result$cost < best_cost - 1e-9) {
              best_cost <- cand_result$cost
              best_bounds <- cand_bounds
              best_k <- cand_result$k_actual
              improved <- TRUE
              cat(sprintf("  ✓ MOVE (Iter %d): Cost=%.6f, K=%d\n", 
                          iter, best_cost, best_k))
              break
            }
          }
          if (improved) break
        }
        if (improved) break
      }
    }
    
    # ==============================
    # SPLIT: Add new boundary
    # ==============================
    if (!improved && length(best_bounds) < K_max + 1) {
      for (interval_idx in 1:(length(best_bounds)-1)) {
        left <- best_bounds[interval_idx]
        right <- best_bounds[interval_idx + 1]
        
        if (right - left > 2) {
          # Try midpoint
          mid <- floor((left + right) / 2)
          cand_bounds <- sort(c(best_bounds, mid))
          
          # Call the helper cost function
          cand_result <- calculate_ls_cost(
            bounds = cand_bounds, 
            p_sorted = p_sorted,
            cov_df = cov_df, 
            delta_dp = delta_dp, 
            K_max = K_max, 
            metric = metric, 
            cov_cols= cov_cols,
            n_r = n_r,
            fixed_redundant_unit_ids = fixed_redundant_unit_ids,
            id_var = id_var,           # Changed to lowercase
            treatment_col = treatment_col,   # Changed to lowercase
            S_inv = S_inv
          )
          
          if (is.finite(cand_result$cost) && 
              cand_result$cost < best_cost - 1e-9) {
            best_cost <- cand_result$cost
            best_bounds <- cand_bounds
            best_k <- cand_result$k_actual
            improved <- TRUE
            cat(sprintf("  ✓ SPLIT (Iter %d): Cost=%.6f, K=%d\n", 
                        iter, best_cost, best_k))
            break
          }
        }
      }
    }
    
    # ==============================
    # MERGE: Remove boundary
    # ==============================
    if (!improved && length(best_bounds) > 2) {
      for (b in 2:(length(best_bounds)-1)) {
        cand_bounds <- best_bounds[-b]
        
        # Call the helper cost function
        cand_result <- calculate_ls_cost(
          bounds = cand_bounds, 
          p_sorted = p_sorted,
          cov_df = cov_df, 
          delta_dp = delta_dp, 
          K_max = K_max, 
          metric = metric, 
          n_r = n_r,
          cov_cols = cov_cols,
          fixed_redundant_unit_ids = fixed_redundant_unit_ids,
          id_var = id_var,           # Changed to lowercase
          treatment_col = treatment_col,   # Changed to lowercase
          S_inv = S_inv
        )
        
        if (is.finite(cand_result$cost) && 
            cand_result$cost < best_cost - 1e-9) {
          best_cost <- cand_result$cost
          best_bounds <- cand_bounds
          best_k <- cand_result$k_actual
          best_intervals <- cand_result$intervals
          improved <- TRUE
          cat(sprintf("  ✓ MERGE (Iter %d): Cost=%.6f, K=%d\n", 
                      iter, best_cost, best_k))
          break
        }
      }
    }
    
    iter <- iter + 1
  }
  
  # Build final intervals
  # final_intervals <- build_intervals_from_bounds(best_bounds, p_sorted, cov_df, treatment_col, id_var, delta_dp)
  
  cat(sprintf("\n✓ Enhanced LS Complete: Final cost = %.6f, K = %d\n", 
              best_cost, best_k))
  
  # No updated has been made
  if (best_cost == initial_cost) {
    return(list(
      intervals = initial_intervals,
      cost = initial_cost,
      boundaries = initial_bounds,
      n_intervals = length(initial_bounds) - 1,
      k_actual = initial_k,
      has_improved = FALSE,
      name = "enhanced_local_seach"
    ))
  }
  
  return(list(
    intervals = best_intervals,
    cost = best_cost,
    boundaries = best_bounds,
    n_intervals = length(best_bounds) - 1,
    k_actual = best_k,
    has_improved = TRUE,
    name = "enhanced_local_seach"
  ))
}

# reconstruct_intervals_and_find_redundant <- function(backtrack, p_sorted, cov_df, 
#                                                      n_r, id_var, treatment_col, 
#                                                      J_bound, dp,
#                                                      final_v_idx, final_k, final_j) {
#   
#   n_t <- sum(cov_df[[treatment_col]] == 1)
#   eps_shift <- 1e-10
#   
#   # Use the passed parameters instead of searching again
#   optimal_v_idx_1based <- final_v_idx
#   optimal_k <- final_k
#   optimal_j <- final_j
#   
#   cat(sprintf("\n=== Starting Reconstruction ===\n"))
#   cat(sprintf("Start State: [i=%d, v=%d (PS=%.6f), k=%d, j=%d]\n",
#               n_t, optimal_v_idx_1based, p_sorted[optimal_v_idx_1based],
#               optimal_k, optimal_j))
#   
#   intervals_list <- list()
#   
#   current_i <- n_t
#   current_v_idx <- optimal_v_idx_1based
#   current_k <- optimal_k
#   current_j <- optimal_j
#   
#   step <- 0
#   
#   while (current_i > 0 || current_v_idx > 1 || current_k > 0) {
#     step <- step + 1
#     
#     bp_list <- backtrack[[current_i + 1, current_v_idx, current_k + 1, current_j + 1]]
#     
#     if (is.null(bp_list) || !is.list(bp_list) || length(bp_list) == 0) {
#       cat(sprintf("\n❌ ERROR at Step %d: Backtrack broken\n", step))
#       cat(sprintf("Current State: [i=%d, v=%d, k=%d, j=%d]\n",
#                   current_i, current_v_idx, current_k, current_j))
#       stop("Backtrack path corrupted")
#     }
#     
#     # ---------------------------------------------------------
#     # --- FIX: Boundary Epsilon Shift for EMPTY intervals ---
#     # ---------------------------------------------------------
#     raw_start <- bp_list$interval_start
#     raw_end   <- bp_list$interval_end
#     is_empty_val <- isTRUE(bp_list$is_empty)
#     
#     # 1. DEFINE BOUNDARIES
#     if (is_empty_val) {
#       # Initialize
#       adj_start <- raw_start
#       adj_end   <- raw_end
#       
#       # Shave boundaries for empty bins to avoid double-counting at edges
#       if (raw_start < 1e-9) { 
#         adj_end <- raw_end - 1e-10  
#       } else if (abs(raw_end - 1.0) < 1e-9) {
#         adj_start <- raw_start + 1e-10 
#       } else {
#         adj_end <- raw_end - 1e-10
#       }
#       
#       # 2. RE-CALCULATE COUNTS for Empty Bin [L, R)
#       calc_n_control <- sum(cov_df$ps >= adj_start & cov_df$ps < adj_end & cov_df[[treatment_col]] == 0)
#       
#       new_interval <- list(
#         start_ps = adj_start,
#         end_ps = adj_end,
#         is_empty = TRUE,
#         n_treated = 0,
#         n_control = calc_n_control, 
#         k_cost = bp_list$rho
#       )
#     } else {
#       # 2. RE-CALCULATE COUNTS for Non-Empty Bin [L, R)
#       # Note: We use the same [L, R) logic here to match Step 2
#       calc_n_control <- sum(cov_df$ps >= raw_start & cov_df$ps < raw_end & cov_df[[treatment_col]] == 0)
#       
#       new_interval <- list(
#         start_ps = raw_start,
#         end_ps = raw_end,
#         is_empty = FALSE,
#         n_treated = bp_list$n_treated,
#         n_control = calc_n_control,
#         k_cost = 1
#       )
#     }
#     
#     intervals_list <- append(intervals_list, list(new_interval))
#     
#     cat(sprintf("Step %d: [%.4f, %.4f] %s (k_cost=%d)\n",
#                 step, new_interval$start_ps, new_interval$end_ps,
#                 if (new_interval$is_empty) "EMPTY" else "NON-EMPTY",
#                 new_interval$k_cost))
#     
#     cat("n_treated in interval: ", new_interval$n_treated, "\n")
#     cat("n_control in interval: ", new_interval$n_control, "\n")
#     
#     
#     # Move to previous state
#     current_i <- bp_list$prev_i
#     current_v_idx <- bp_list$prev_v_idx
#     current_k <- bp_list$prev_k
#     current_j <- bp_list$prev_j
#   }
#   
#   cat(sprintf("\n✅ Reconstruction complete: %d intervals traced\n", step))
#   
#   # Verify K cost
#   total_k_cost <- sum(sapply(intervals_list, function(x) x$k_cost))
#   if (abs(total_k_cost - optimal_k) > 1e-9) {
#     warning(sprintf("K cost mismatch: reconstructed=%d, optimal=%d", 
#                     total_k_cost, optimal_k))
#   }
#   
#   # We exit the loop when current_k=0 AND current_i=0 AND current_v_idx=1 (start of PS range)
#   
#   if (sum(sapply(intervals_list, function(x) x$k_cost)) != optimal_k) {
#     # This check is crucial for debugging path failures
#     warning(paste("Reconstructed K cost (", sum(sapply(intervals_list, function(x) x$k_cost)), ") does not match optimal K (", optimal_k, ")"))
#   }
#   
#   # ----------------------------------------------------
#   # --- FIX 1: Final Data Frame Creation with k_cost ---
#   # ----------------------------------------------------
#   
#   final_intervals <- data.frame(
#     start_ps = sapply(rev(intervals_list), function(x) x$start_ps),
#     end_ps = sapply(rev(intervals_list), function(x) x$end_ps),
#     n_treated = sapply(rev(intervals_list), function(x) x$n_treated),
#     n_control = sapply(rev(intervals_list), function(x) x$n_control),
#     is_empty = sapply(rev(intervals_list), function(x) x$is_empty),
#     k_cost = sapply(rev(intervals_list), function(x) x$k_cost)
#   )
# 
#   # --- 2. PREPARE DATA STRUCTURE FOR GLOBAL PS COST COLLECTION ---
#   all_treated_ps_costs <- data.frame(
#     id = character(),
#     ps_cost = numeric()
#   )
#   
#   for (r in seq_len(nrow(final_intervals))) {
#     vprev <- final_intervals$start_ps[[r]]
#     vcurr <- final_intervals$end_ps[[r]]
#     
#     # 🛑 FIX: Use Left-Closed, Right-Open logic [L, R)
#     # This ensures a unit at exactly 0.5 is counted in the bin starting at 0.5
#     if (r < nrow(final_intervals)) {
#       bin_mask <- (cov_df$ps >= vprev) & (cov_df$ps < vcurr)
#     } else {
#       # Last bin remains closed on both ends
#       bin_mask <- (cov_df$ps >= vprev) & (cov_df$ps <= vcurr)
#     }
#     
#     bin_data <- cov_df[bin_mask, ]
#     
#       treated_in_bin <- bin_data[bin_data[[treatment_col]] == 1, ]
#       control_in_bin <- bin_data[bin_data[[treatment_col]] == 0, ]
# 
#       if (nrow(treated_in_bin) == 0) next
# 
#       c_ps_in_bin <- control_in_bin$ps
#       has_controls <- length(c_ps_in_bin) > 0
# 
#       for (i in seq_len(nrow(treated_in_bin))) {
#         t_id <- as.character(treated_in_bin[[id_var]][i]) # Ensure ID is character for rbind consistency
#         t_ps <- treated_in_bin$ps[i]
# 
#         if (has_controls) {
#           ps_cost <- min(abs(t_ps - c_ps_in_bin), na.rm = TRUE)
#         } else {
#           ps_cost <- Inf
#         }
# 
#         all_treated_ps_costs <- rbind(all_treated_ps_costs,
#                                       data.frame(id = t_id, ps_cost = ps_cost))
#       }
#     }
#   
#   # --- 3. ITERATE THROUGH EVERY INTERVAL TO CALCULATE INDIVIDUAL PS COST ---
#   # for (r in seq_len(nrow(final_intervals))) {
#   #   vprev <- final_intervals$start_ps[[r]]
#   #   vcurr <- final_intervals$end_ps[[r]]
#   #   
#   #   # In the all_treated_ps_costs loop:
#   #   # Use epsilon on BOTH sides to be "boundary-blind" for matching
#   #   bin_data <- cov_df[cov_df$ps > (vprev - 1e-11) & cov_df$ps <= (vcurr + 1e-11), ]
#   #   
#   #   treated_in_bin <- bin_data[bin_data[[treatment_col]] == 1, ]
#   #   control_in_bin <- bin_data[bin_data[[treatment_col]] == 0, ]
#   #   
#   #   if (nrow(treated_in_bin) == 0) next
#   #   
#   #   c_ps_in_bin <- control_in_bin$ps
#   #   has_controls <- length(c_ps_in_bin) > 0
#   #   
#   #   for (i in seq_len(nrow(treated_in_bin))) {
#   #     t_id <- as.character(treated_in_bin[[id_var]][i]) # Ensure ID is character for rbind consistency
#   #     t_ps <- treated_in_bin$ps[i]
#   #     
#   #     if (has_controls) {
#   #       ps_cost <- min(abs(t_ps - c_ps_in_bin), na.rm = TRUE)
#   #     } else {
#   #       ps_cost <- Inf
#   #     }
#   #     
#   #     all_treated_ps_costs <- rbind(all_treated_ps_costs,
#   #                                   data.frame(id = t_id, ps_cost = ps_cost))
#   #   }
#   # }
#   
#   # --- 4. IDENTIFY THE n_r WORST-MATCHED UNITS ---
#   fixed_redundant_unit_ids <- character(0)
#   total_initial_cost <- sum(all_treated_ps_costs$ps_cost, na.rm = TRUE)
#   min_cost_recalc <- total_initial_cost
#   
#   if (n_r > 0 && nrow(all_treated_ps_costs) > 0) {
#     # Sort by PS cost in descending order
#     sorted_costs <- all_treated_ps_costs[order(all_treated_ps_costs$ps_cost, decreasing = TRUE), ]
#     
#     n_r_used <- min(n_r, nrow(all_treated_ps_costs))
#     
#     # Select the unique IDs of the top n_r units
#     fixed_redundant_unit_ids <- sorted_costs$id[1:n_r_used]
#     
#     # Calculate the final minimized cost (sum of remaining costs)
#     if (n_r_used < nrow(sorted_costs)) {
#       min_cost_recalc <- sum(sorted_costs$ps_cost[(n_r_used + 1):nrow(sorted_costs)])
#     } else {
#       min_cost_recalc <- 0 # All units were redundant
#     }
#   }
#   
#   return(list(
#     intervals = final_intervals,
#     cost = min_cost_recalc,
#     redundant_unit_ids = fixed_redundant_unit_ids
#   ))
# }

# ============================================================
# UPDATED HELPER: Boundary Index Conversion
# ============================================================
# boundaries_from_intervals <- function(intv, ps_sorted) {
#   b <- c()
#   n_ps <- length(ps_sorted)
#   
#   # Debug
#   # browser()
#   
#   for(i in seq_len(nrow(intv))) {
#     s <- intv$start_ps[i]  # Changed from start_ps to start
#     e <- intv$end_ps[i]    # Changed from end_ps to end
#     
#     # 1. Start boundary: Find index where ps_sorted is closest to s
#     if (s == 0) {
#       idx_s <- 1
#     } else {
#       matches_s <- which(abs(ps_sorted - s) < 1e-9)
#       if (length(matches_s) > 0) {
#         idx_s <- matches_s[1]
#       } else {
#         # Find closest value
#         idx_s <- which.min(abs(ps_sorted - s))
#       }
#     }
#     
#     # Debug
#     # browser()
#     
#     # 2. End boundary: Find index where ps_sorted is closest to e
#     if (e == 1) {
#       idx_e <- n_ps
#     } else {
#       matches_e <- which(abs(ps_sorted - e) < 1e-9)
#       if (length(matches_e) > 0) {
#         idx_e <- matches_e[1]
#       } else {
#         # Find closest value
#         idx_e <- which.min(abs(ps_sorted - e))
#       }
#     }
#     
#     b <- c(b, idx_s, idx_e)
#   }
#   
#   # Debug
#   # browser()
#   
#   # Ensure boundaries are valid, unique, and within range [0, n_ps]
#   b <- sort(unique(b))
#   b <- b[b >= 1 & b <= n_ps]
#   
#   # Debug
#   # browser()
#   
#   return(b)
# }

boundaries_from_intervals <- function(intv, ps_sorted, left_border, right_border) {
  b <- c()
  n_ps <- length(ps_sorted)
  eps <- 1e-9
  
  for(i in seq_len(nrow(intv))) {
    s <- intv$start_ps[i]
    e <- intv$end_ps[i]
    
    # 1. Start boundary: Check against left_border instead of 0
    if (abs(s - left_border) < eps) {
      idx_s <- 1
    } else {
      idx_s <- which.min(abs(ps_sorted - s))
    }
    
    # 2. End boundary: Check against right_border instead of 1
    if (abs(e - right_border) < eps) {
      idx_e <- n_ps
    } else {
      idx_e <- which.min(abs(ps_sorted - e))
    }
    
    b <- c(b, idx_s, idx_e)
  }
  
  # Ensure boundaries are valid and unique
  return(sort(unique(b)))
}

check_feasible <- function(intervals) {
  
  for(i in seq_len(nrow(intervals))) {
  
    n_treated_bin_i <- intervals$n_treated[i]
    n_control_bin_i <- intervals$n_control[i]
  
    if (n_treated_bin_i == 0) {
      next  # Skip to next interval
    }
    
    # Violating the fine-constraints
    if (n_treated_bin_i > n_control_bin_i) {
      return(FALSE)
    }
  }
  return(TRUE)
}

run_sensitivity_analysis <- function(matched_data, outcome_var, treatment_var, match_id_var, gamma_range = seq(1, 3, by = 0.1)) {
  # senm requires matched sets. We need to extract the outcome for treated and controls
  # within the same match_id (subclass)
  
  # Ensure data is sorted for consistent set extraction
  df <- matched_data %>% 
    select(all_of(c(outcome_var, treatment_var, match_id_var))) %>%
    arrange(!!sym(match_id_var), desc(!!sym(treatment_var)))
  
  # Extract vectors for senm
  # y: Outcome, z: Treatment, m: Match ID
  y <- df[[outcome_var]]
  z <- df[[treatment_var]]
  m <- df[[match_id_var]]
  
  results <- data.frame(Gamma = numeric(), P_Value = numeric())
  
  for (g in gamma_range) {
    # inner function to compute upper bound p-value for a given gamma
    # sensitive analysis for M-tests
    res <- senm(y, z, m, gamma = g, inner = 0, trim = 3)
    results <- rbind(results, data.frame(Gamma = g, P_Value = res$pval))
    
    # Optional: Stop if p-value crosses 0.05 to find tipping point faster
    if (res$pval > 0.05) break
  }
  
  # Return the maximum Gamma where p-value was still < 0.05
  tipping_point <- results$Gamma[max(which(results$P_Value < 0.05))]
  return(list(results = results, tipping_point = tipping_point))
}

# source("dataset_configs.R")
# source("levin_functions.R", echo = FALSE)
# 
# # Reconstructing the "better" version (Cost = 157.11)
# best_intervals_old_better <- data.frame(
#   start_ps  = c(0.000, 0.195, 0.308, 0.555, 0.676),
#   end_ps    = c(0.195, 0.308, 0.555, 0.676, 1.000),
#   n_treated = c(0, 10, 154, 21, 0),
#   n_control = c(1, 26, 204, 29, 0),
#   is_empty  = c(FALSE, FALSE, FALSE, FALSE, TRUE),
#   k_cost    = c(1, 1, 1, 1, 2) # Example cost structure from your screenshot
# )
# 
# # Standardized data_subset
# cov_df_standardized = prepare_standardized_data(data_subset, data_config, treatment_col)
# 
# # Define the formula from your config
# match_formula <- reformulate(data_config$NUMERIC_COVARIATES, treatment_col)
# 
# # Generate the standard Mahalanobis distance matrix
# # This is what Optimal, Cardinality (Step B), and Levin will all use.
# 
# dist_obj <- match_on(match_formula, data = cov_df_standardized, method = "mahalanobis")
# 
# dist_matrix <- as.matrix(dist_obj)
# 
# # Pre-compute S_inv ONCE
# X_for_S <- as.matrix(cov_df_standardized[, data_config$NUMERIC_COVARIATES])
# S_inv <- solve(cov(X_for_S) + diag(1e-7, ncol(X_for_S)))
# 
# if (is.null(S_inv)) stop(paste(metric, "matrix inversion failed."))
# 
# # params$S_inv <- S_inv
# 
# # test_params <- PARAMS
# 
# data_config <- NSW_MIXTAPE_CONFIG
# treatment_col <- data_config$TREATMENT_VAR
# 
# cost_params <- list(
#   intervals = best_intervals_old_better,
#   p_sorted = p_sorted,
#   cov_df = cov_df_standardized,
#   delta_dp = PARAMS$DELTA_DP,
#   K_max = PARAMS$K_DP,
#   metric = PARAMS$ORIGINAL_DISTANCE_METRIC,
#   cov_cols = data_config$NUMERIC_COVARIATES,
#   n_r = PARAMS$N_R,
#   fixed_redundant_unit_ids = c(0),
#   id_var = data_config$ID_VAR,
#   S_inv = S_inv,
#   treatment_col = treatment_col,
#   data_config = data_config
# )
# 
# call_ls_cost <- function(bounds) {
#   do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
# }
# 
# old_bounds <- boundaries_from_intervals(best_intervals_old_better, p_sorted)
# 
# # Debug
# browser()
# 
# cand_result <- call_ls_cost(old_bounds)
# 
# # Debug
# browser()
# 
# print("cand_result is: \n")
# print(cand_result)

