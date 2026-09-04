# 
# 
# # original_distance_for_intervals.R


# Simplified version - uses stored IDs, no recalculation needed

# calculate_distance_from_dp_intervals <- function(intervals, cov_df, 
#                                                  # fixed_redundant_unit_ids, 
#                                                  id_var, treatment_col, X_all,
#                                                  S_inv, ...) {
#   
#   total_cost <- 0
#   all_t_indices_captured <- c()
#   n_intervals <- nrow(intervals)
#   
#   # Setup the ID-to-Row Map
#   id_to_row_map <- seq_len(nrow(cov_df))
#   names(id_to_row_map) <- as.character(cov_df[[id_var]])
#   
#   # Get ID lists from attributes
#   t_id_list <- attr(intervals, "treatment_unit_ids")
#   c_id_list <- attr(intervals, "control_unit_ids")
#   
#   # Validate that we have the attributes
#   if (is.null(t_id_list) || is.null(c_id_list)) {
#     stop("intervals must have treatment_unit_ids and control_unit_ids attributes")
#   }
#   
#   for (r in seq_len(n_intervals)) {
#     if (intervals$is_empty[r]) next
#     
#     # ═══════════════════════════════════════════════════════════
#     # Get treatment and control IDs directly from attributes
#     # ═══════════════════════════════════════════════════════════
#     
#     t_ids <- as.character(t_id_list[[r]])
#     c_ids <- as.character(c_id_list[[r]])
#     
#     # Convert IDs to row indices
#     t_indices <- id_to_row_map[t_ids]
#     c_indices <- id_to_row_map[c_ids]
#     
#     # Clean indices (remove NAs)
#     t_indices <- as.integer(t_indices[!is.na(t_indices)])
#     c_indices <- as.integer(c_indices[!is.na(c_indices)])
#     
#     # No need for that when we remove the n_r (slack) vertices
#     # Filter out redundant units (jokers)
#     # t_indices_final <- t_indices[!(as.character(cov_df[[id_var]][t_indices]) %in% fixed_redundant_unit_ids)]
#     # 
#     # all_t_indices_captured <- c(all_t_indices_captured, t_indices_final)
#     
#     # ═══════════════════════════════════════════════════════════
#     # Validate and calculate distances
#     # ═══════════════════════════════════════════════════════════
#     
#     if (length(t_indices) == 0) next
#     
#     if (length(c_indices) == 0) {
#       return(list(cost = Inf, n_t_matched = length(unique(t_indices))))
#     }
#     
#     # Matrix math (unchanged - this part is efficient)
#     X_t <- X_all[t_indices, , drop = FALSE]
#     X_c <- X_all[c_indices, , drop = FALSE]
#     
#     X_t_transformed <- X_t %*% S_inv
#     
#     t_diag <- rowSums(X_t_transformed * X_t)
#     c_diag <- rowSums((X_c %*% S_inv) * X_c)
#     
#     cross_term <- X_t_transformed %*% t(X_c)
#     
#     dist_matrix <- outer(t_diag, c_diag, "+") - 2 * cross_term
#     
#     min_distances <- apply(dist_matrix, 1, function(x) {
#       val <- min(x, na.rm = TRUE)
#       return(sqrt(pmax(0, val)))
#     })
#     
#     total_cost <- total_cost + sum(min_distances)
#   }
#   
#   return(list(cost = total_cost, n_t_matched = length(unique(t_indices))))
# }

calculate_distance_from_bounds_units <- function(bounds, cov_df, metric,
                                           # fixed_redundant_unit_ids, 
                                           id_var, treatment_col, X_all,
                                           S_inv, cov_cols, delta_dp, K_max, 
                                           p_sorted, left_border, right_border,
                                           intervals = NULL) {
  
  if (is.null(intervals)) {
  # Step 1: Rebuild intervals from bounds
    intervals <- build_intervals_from_bounds_units(bounds, p_sorted, cov_df, treatment_col, id_var, delta_dp, left_border)
  } else {
    intervals = intervals
  }
  
  # Step 2: Combined Feasibility Check
  # This catches: 1) Delta-width violations (NULL) 
  #               2) N_treated > N_control violations
  if (is.null(intervals) || !check_feasible(intervals)) {
    return(list(cost = Inf, k_actual = NA))
  }
  
  # Step 3: K calculation
  # TODO: Fix this! Should include empty intervals!!!  
  
  # TODO: Add calculation for the empty intervals (like rho in the dp),

  k_actual <- sum(intervals$k_cost)
  
  if (k_actual > K_max) {
    return(list(cost = Inf, k_actual = k_actual))
  }
  
  cost_result <- calculate_distance_from_dp_intervals_local(
    intervals = intervals,
    cov_df = cov_df,
    # fixed_redundant_unit_ids = fixed_redundant_unit_ids,
    id_var = id_var,
    treatment_col = treatment_col,
    X_working = X_all,
    S_inv_working = S_inv,
    cov_cols = cov_cols,
    # n_r = n_r # Pass n_r through
  )
  
  return(list(cost = cost_result, k_actual = k_actual, intervals = intervals))
}
