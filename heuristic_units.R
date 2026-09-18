# -- heuristic_units.R --

build_intervals_from_unit_indices <- function(
    bounds_idx, 
    treatment_ids, 
    treatment_scores, 
    control_scores, 
    cov_df, 
    treatment_col, 
    id_var, 
    delta_dp
) {
  n_bins <- length(bounds_idx) - 1
  intervals <- data.frame(
    start_ps       = rep(NA_real_, n_bins),
    end_ps         = rep(NA_real_, n_bins),
    n_treated      = rep(0, n_bins),
    n_control      = rep(0, n_bins),
    is_empty       = rep(TRUE, n_bins),
    k_cost         = rep(1, n_bins)
  )
  
  t_ids_list <- vector("list", n_bins)
  c_ids_list <- vector("list", n_bins)
  control_ids <- names(control_scores)
  
  for (i in 1:n_bins) {
    # 1. Identify indices for current bin
    idx_start <- bounds_idx[i] + 1
    idx_end   <- bounds_idx[i+1]
    
    if (idx_start > idx_end) next
    
    # 2. Get scores for current bin
    bin_t_scores <- treatment_scores[idx_start:idx_end]
    v_min <- min(bin_t_scores)
    v_max <- max(bin_t_scores)
    
    # --- NEW: Calculate the neighbors for the boundary math ---
    
    # prev_max is the score of the LAST unit in the PREVIOUS bin
    if (i > 1) {
      prev_max <- treatment_scores[bounds_idx[i]]
    }
    
    # next_min is the score of the FIRST unit in the NEXT bin
    if (i < n_bins) {
      next_min <- treatment_scores[bounds_idx[i+1] + 1]
    }
    
    # 3. Calculate Boundaries with Snap-to-7-decimals
    bin_left  <- if (i > 1) round((prev_max + v_min) / 2, 7) else 0
    bin_right <- if (i < n_bins) round((v_max + next_min) / 2, 7) else 1
    
    # 4. Define the feasible range for controls (PS-space constraints)
    low_lim  <- max(0, v_max - delta_dp, bin_left)
    high_lim <- min(1, v_min + delta_dp, bin_right)
    
    # 5. [L, R) Logic for static controls
    if (i < n_bins) {
      feasible_mask <- control_scores >= low_lim & control_scores < high_lim
    } else {
      feasible_mask <- control_scores >= low_lim & control_scores <= high_lim
    }
    
    bin_c_ids <- control_ids[feasible_mask]
    
    intervals$start_ps[i]  <- bin_left
    intervals$end_ps[i]    <- bin_right
    intervals$n_treated[i] <- length(bin_t_scores)
    intervals$n_control[i] <- length(bin_c_ids)
    intervals$is_empty[i]  <- FALSE
    
    t_ids_list[[i]] <- treatment_ids[idx_start:idx_end]
    c_ids_list[[i]] <- bin_c_ids
  }
  
  attr(intervals, "treatment_unit_ids") <- t_ids_list
  attr(intervals, "control_unit_ids")   <- c_ids_list
  return(intervals)
}

calculate_bin_mahalanobis <- function(t_ids, c_ids, cov_df, cov_cols, S_inv, id_var) {
  if (length(c_ids) == 0) return(Inf)
  
  # 1. Get control centroid (the 'Target' for this bin)
  x_c_mat <- as.matrix(cov_df[cov_df[[id_var]] %in% c_ids, cov_cols])
  centroid_c <- colMeans(x_c_mat)
  
  # 2. Get treated units matrix
  x_t_mat <- as.matrix(cov_df[cov_df[[id_var]] %in% t_ids, cov_cols])
  
  # 3. Calculate distance for each treated unit to that centroid
  # (x_i - centroid) %*% S_inv %*% (x_i - centroid)
  diffs <- sweep(x_t_mat, 2, centroid_c)
  
  # Vectorized Mahalanobis calculation
  dists <- sqrt(rowSums((diffs %*% S_inv) * diffs))
  
  return(sum(dists))
}

build_intervals_from_bounds_units <- function(bounds_indices, p_sorted, cov_df, treatment_col, id_var, delta_dp, left_border) {
  
  eps <- 1e-9
  
  # --- UPDATE 1: ANCHORING ---
  # Map indices to PS values and force the boundaries to be watertight
  bounds_ps <- p_sorted[bounds_indices]
  bounds_ps[1] <- left_border
  bounds_ps[length(bounds_ps)] <- max(p_sorted) 
  
  # --- UPDATE 2: ID-BASED BINNING ---
  # This replaces the need for 'get_units_in_interval'
  # It assigns every row in cov_df to an interval index (1 to n_intervals)
  # 'all.inside = TRUE' ensures units at the absolute max PS are in the last bin
  interval_assignments <- findInterval(cov_df$ps, bounds_ps, all.inside = TRUE)
  
  n_intervals <- length(bounds_ps) - 1
  interval_list <- vector("list", n_intervals)
  t_id_attr <- vector("list", n_intervals)
  c_id_attr <- vector("list", n_intervals)
  
  for (i in 1:n_intervals) {
    vprev <- bounds_ps[i]
    vcurr <- bounds_ps[i+1]
    
    # Select the IDs assigned to this bin
    subset_df <- cov_df[interval_assignments == i, ]
    
    if (nrow(subset_df) == 0) {
      distance <- vcurr - vprev
      rho <- max(1, ceiling(distance / delta_dp))
      
      interval_list[[i]] <- data.frame(
        start_ps = vprev, end_ps = vcurr,
        n_treated = 0, n_control = 0,
        is_empty = TRUE, k_cost = rho,
        stringsAsFactors = FALSE
      )
      t_id_attr[[i]] <- character(0)
      c_id_attr[[i]] <- character(0)
      
    } else {
      # Delta-width Feasibility Check
      if (vcurr - vprev > (delta_dp + eps)) return(NULL)
      
      t_ids <- subset_df[[id_var]][subset_df[[treatment_col]] == 1]
      c_ids <- subset_df[[id_var]][subset_df[[treatment_col]] == 0]
      
      interval_list[[i]] <- data.frame(
        start_ps = vprev, end_ps = vcurr,
        n_treated = length(t_ids), n_control = length(c_ids),
        is_empty = FALSE, k_cost = 1,
        stringsAsFactors = FALSE
      )
      
      t_id_attr[[i]] <- as.character(t_ids)
      c_id_attr[[i]] <- as.character(c_ids)
    }
  }
  
  res <- do.call(rbind, interval_list)
  attr(res, "treatment_unit_ids") <- t_id_attr
  attr(res, "control_unit_ids") <- c_id_attr
  
 # Final Validation
  total_captured_t <- sum(res$n_treated)
  total_actual_t <- sum(cov_df[[treatment_col]] == 1)

  if (total_captured_t != total_actual_t) {

    # 1. Identify which IDs are actually in the data
    actual_ids <- cov_df[[id_var]][cov_df[[treatment_col]] == 1]

    # 2. Identify which IDs were assigned to intervals
    captured_ids <- unlist(t_id_attr)

    # 3. Find the missing and the double-counted
    missing_ids <- setdiff(actual_ids, captured_ids)
    extra_ids   <- captured_ids[duplicated(captured_ids)]

    message("--- Forensic Report ---")
    if(length(missing_ids) > 0) {
      cat("Missing IDs PS scores:\n")
      print(sort(ps_scores[missing_ids]))
    }
    if(length(extra_ids) > 0) {
      cat("Double-counted IDs PS scores:\n")
      print(ps_scores[extra_ids])
    }

    # 4. Check the Bounds vs p_sorted range
    cat(sprintf("\nBounds PS Range: [%.6f, %.6f]\n", min(bounds_ps), max(bounds_ps)))
    cat(sprintf("Data PS Range:   [%.6f, %.6f]\n", min(ps_scores), max(ps_scores)))

    # --- FORENSIC DEBUG START ---
    cat("\n--- PROPENSITY SCORE ALIGNMENT CHECK ---\n")
    cat("Left Border (Parameter): ", sprintf("%.16f", left_border), "\n")
    cat("Min p_sorted (Data):     ", sprintf("%.16f", min(p_sorted)), "\n")
    cat("Difference:              ", left_border - min(p_sorted), "\n")

    cat("\n--- DATASET SIZE ALIGNMENT ---\n")
    cat("Rows in cov_df:          ", nrow(cov_df), "\n")
    cat("Elements in p_sorted:    ", length(p_sorted), "\n")

    # Check for duplicate propensity scores that might be causing index confusion
    cat("Unique p_sorted:         ", length(unique(p_sorted)), "\n")
    # ----------------------------

    cond1_flag <- min(p_sorted) == left_border

    cond2_flag <- length(p_sorted) == nrow(cov_df)

    cat("Cond min(p_sorted) == left_border is: ", cond1_flag, "\n")
    cat("Cond min(p_sorted) == left_border is: ", cond2_flag, "\n")

    warning(sprintf("⚠️ Membership mismatch: %d captured, %d actual", total_captured_t, total_actual_t))
    return(NULL)
  }
  return(res)
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

calculate_distance_from_bounds_units <- function(bounds, cov_df, metric,
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
  
  k_actual <- sum(intervals$k_cost)
  
  if (k_actual > K_max) {
    return(list(cost = Inf, k_actual = k_actual))
  }
  
  cost_result <- calculate_distance_from_dp_intervals_local(
    intervals = intervals,
    cov_df = cov_df,
    id_var = id_var,
    treatment_col = treatment_col,
    X_working = X_all,
    S_inv_working = S_inv,
    cov_cols = cov_cols
  )
  
  return(list(cost = cost_result, k_actual = k_actual, intervals = intervals))
}

run_osip_step2_heuristics_units <- function(osip_res_step1, initial_cost, 
                                             initial_bounds, 
                                             initial_intervals, cov_df, params, 
                                             data_config, p_sorted, 
                                             treatment_col, X_all, S_inv, delta_dp, 
                                             k_bound, distance_metric, left_border, 
                                             right_border, dist_matrix, base_dir, robust_flag, 
                                             matching_non_1_1, dataset_name, analyze_flag) {
  
  id_var <- data_config$ID_VAR
  initial_k <- length(initial_intervals)
  
  cat("distance_metric is: ", distance_metric, "\n")
  cat(sprintf("DEBUG: S_inv trace is %f\n", sum(diag(S_inv))))
  
  # --- 1. PREPARATION ---
  treatment_map_initial <- cov_df %>%
    dplyr::filter(!!sym(treatment_col) == 1) %>%
    dplyr::select(!!sym(id_var), ps) %>%
    dplyr::arrange(ps)
  
  scores_res <- get_treatment_control_scores(treatment_map_initial, cov_df, 
                                             data_config$ID_VAR, treatment_col)
  
  t_scores <- scores_res$treatment_scores
  c_scores <- scores_res$control_scores
  
  # --- 2. Run and Process Searches ---
  opt_funcs <- list(
    "local_search"        = local_search_units,
    "enhanced_ls"         = enhanced_local_search_units,
    "simulated_annealing" = simulated_annealing_units
  )
  
  defaults <- get_default_search_params(params, data_config, treatment_col, X_all,
                                        S_inv, delta_dp, k_bound)
  
  final_results <- list()
  # --- 2. Run and Process Searches ---
  all_trace_reports <- list() # The "Full Picture" collector
  
  osip_res_step_1_row <- extract_step_row_stats(
    m_res            = osip_res_step1,        # Step 1 matched object
    trace_name       = "DP_Step1",
    est_cost         = initial_cost,           # DP initial cost
    data_config      = data_config,
    treatment_col    = treatment_col,
    shifting_point   = osip_res_step1$shifting_point,
    matching_non_1_1 = matching_non_1_1
  )
  
  for (method in params$SEARCHES_TO_USE) {
    cat(sprintf("\nLaunching %s (Unit-Space)...", method))
    
    overrides <- list(
      initial_intervals  = initial_intervals,
      initial_bounds_idx = initial_bounds,
      initial_cost       = initial_cost,
      initial_k          = initial_k,
      p_sorted           = p_sorted,
      metric             = distance_metric,
      treatment_ids      = names(t_scores),
      treatment_scores   = t_scores,
      control_scores     = c_scores,
      cov_df             = cov_df,
      left_border = left_border,
      right_border  = right_border
    )
    
    raw_res <- run_search(method, defaults, overrides, opt_funcs)
    
    # ══════════════════════════════════════════════════════════════════════════
    # 🔍 TRACE ANALYSIS (Supervisor Request)
    # ══════════════════════════════════════════════════════════════════════════
    if (analyze_flag && !is.null(raw_res$trace_log) && length(raw_res$trace_log) > 0) {
      cat(sprintf("\n[Trace] Analyzing %d improvement steps for %s...", length(raw_res$trace_log), method))
      
      # We call the analyzer here. Note: ensure all these variables (cov_df, dist_matrix, etc.) 
      # are accessible in the function scope.
      method_trace_report <- analyze_optimization_trace(
        osip_res_step_1_row = osip_res_step_1_row,
        method_name     = method,
        trace           = raw_res$trace_log,
        cov_df          = cov_df,
        dist_matrix     = dist_matrix, # Ensure this is passed into run_osip_step2
        data_config     = data_config,
        params          = params,
        treatment_col   = treatment_col,
        outcome_var     = data_config$OUTCOME_VAR, # Or your specific outcome variable
        metric_label    = distance_metric,
        delta_dp        = delta_dp,
        base_dir        = base_dir,
        robust_flag     = robust_flag,
        matching_non_1_1 = matching_non_1_1,
        dataset_name = dataset_name
      )
      
      # Store it in our list with a clear name
      all_trace_reports[[method]] <- method_trace_report
    }
    # ══════════════════════════════════════════════════════════════════════════
    
    final_results[[method]] <- process_and_report_search(
      ls_res            = raw_res, 
      search_name       = method, 
      initial_cost      = initial_cost,
      initial_intervals = initial_intervals,
      initial_k         = initial_k
    )
    
    save_partition_csv(
      df             = raw_res$intervals,
      delta_val      = delta_dp,
      base_dir       = base_dir,
      step           = 2,
      heuristic_name = method,
      step2_flag     = robust_flag
    )
  }
  
  # --- 3. Selection Logic ---
  improved_mask <- vapply(final_results, function(res) isTRUE(res$status), logical(1))
  
  if (!any(improved_mask)) {
    fail_res <- list(data_matched = NULL, name = "No Improvement")
    attr(fail_res, "is_step2") <- FALSE
    return(fail_res)
  }
  
  improved_results <- final_results[improved_mask]
  best_sub_idx     <- which.min(sapply(improved_results, function(res) res$cost))
  best_res         <- improved_results[[best_sub_idx]]
  
  # ══════════════════════════════════════════════════════════════════════════
  # 🔑 STANDARDIZATION FIX
  # ══════════════════════════════════════════════════════════════════════════
  
  # Ensure best_res$intervals (which comes from process_and_report_search -> build_intervals_from_bounds)
  # has the required attributes. No manual reconstruction or nudging needed here.
  
  if (is.null(attr(best_res$intervals, "treatment_unit_ids"))) {
    # If the attributes are missing, it means the search/process functions 
    # aren't calling the updated build_intervals_from_bounds yet.
    stop("Search results missing required ID attributes. Check process_and_report_search.")
  }
  
  # We no longer modify viz_data PS scores. We use the actual PS scores 
  # and let get_units_in_interval handle the boundaries.
  # best_res$viz_data <- cov_df 
  
  attr(best_res, "best_method") <- best_res$name
  attr(best_res, "is_step2")    <- TRUE
  
  return(best_res)
}

plot_partition_with_bins <- function(data_subset, intervals_df, 
                                     # redundant_unit_ids, 
                                     title,
                                     data_config, delta_dp, dataset_name,
                                     title_suffix, base_dir, treatment_col, 
                                     show_plot = TRUE, make_interactive = TRUE) {
  
  # 1. Prepare ID mappings from attributes
  t_id_list <- attr(intervals_df, "treatment_unit_ids")
  c_id_list <- attr(intervals_df, "control_unit_ids")
  id_var <- data_config$ID_VAR
  n_bins <- nrow(intervals_df)
  
  # 2. Build plotting data with "Cloud" positioning
  plot_list <- list()
  for (i in 1:n_bins) {
    bin_ids <- c(t_id_list[[i]], c_id_list[[i]])
    if (length(bin_ids) == 0) next
    
    bin_data <- data_subset[data_subset[[id_var]] %in% bin_ids, ]
    
    # ══════════════════════════════════════════════════════════
    # CLUSTER LOGIC:
    # X is the bin index + horizontal jitter to create the "cloud"
    # Y is a fixed level (high for Treated, low for Control)
    # ══════════════════════════════════════════════════════════
    bin_data$bin_idx <- i
    bin_data$x_jittered <- jitter(rep(i, nrow(bin_data)), amount = 0.3)
    bin_data$y_fixed <- ifelse(bin_data[[treatment_col]] == 1, 1.0, 0.2)
    
    bin_data$group_label <- ifelse(bin_data[[treatment_col]] == 1, "Treated", "Control")
    # bin_data$group_label[bin_data[[id_var]] %in% redundant_unit_ids] <- "Redundant Treated"
    
    # Tooltip shows the REAL PS so you don't lose the technical info
    bin_data$tooltip <- sprintf("ID: %s\nPS: %.4f\nBin: %d", 
                                bin_data[[id_var]], bin_data$ps, i)
    
    plot_list[[i]] <- bin_data
  }
  plot_df <- do.call(rbind, plot_list)
  
  # 3. Create Plot
  p <- ggplot(plot_df) +
    # Use the jittered X and fixed Y to create horizontal clusters
    geom_point(aes(x = x_jittered, y = y_fixed, color = group_label, text = tooltip),
               alpha = 0.7, size = 2.5) +
    # Bold solid black vertical lines between bins
    geom_vline(xintercept = seq(0.5, n_bins + 0.5, by = 1), 
               linewidth = 1.0, color = "black", alpha = 0.6) +
    scale_color_manual(values = c("Control" = "blue", "Treated" = "red")) +
    theme_minimal() +
    labs(title = title, x = "Propensity Score Bin (Index)", y = "") +
    # Hide Y axis labels to keep it clean like your example
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(), axis.text.x = element_blank()) +
    coord_cartesian(xlim = c(0.5, n_bins + 0.5), ylim = c(-0.2, 1.5))
  
  # 4. Add Bin Labels and PS Ranges
  for (i in 1:n_bins) {
    ps_range <- sprintf("[%.3f, %.3f]", intervals_df$start_ps[i], intervals_df$end_ps[i])
    label_main <- if(intervals_df$is_empty[i]) sprintf("Bin %d\n(Empty)", i) else 
      sprintf("Bin %d\nT=%d, C=%d", i, intervals_df$n_treated[i], intervals_df$n_control[i])
    
    p <- p + annotate("text", x = i, y = 1.4, label = label_main, color = "darkgreen", fontface = "bold", size = 3) +
      annotate("text", x = i, y = -0.1, label = ps_range, color = "darkgray", size = 2.5)
  }
  
  # Final Export Logic (unchanged)
  filename <- sprintf("bins_plot_%s_step_%s_delta_%.2f.png", 
                      dataset_name, title_suffix, delta_dp)
  
  filepath <- file.path(base_dir, filename)
  
  ggsave(filepath, plot = p, width = 10, height = 6, dpi = 300)
  
  cat(sprintf("✓ Plot saved: %s\n", filepath))
  
  if (show_plot) {
    if (make_interactive && requireNamespace("plotly", quietly = TRUE)) {
      interactive_p <- plotly::ggplotly(p, tooltip = "text")
      print(interactive_p)
      return(interactive_p)
    } else {
      print(p)
      return(p)
    }
  }
  
  return(p)
}


get_treatments_in_interval_dp <- function(treatment_scores, v_left, v_right,
                                          T_sorted, eps = 1e-9) {
 
  #' @param treatment_scores Named numeric vector of treatment PS (names = unit IDs)
  #' @param v_left  Left boundary
  #' @param v_right Right boundary
  #' @param T_sorted Sorted numeric vector of ALL treatment PS values
  #'                 Used to detect when v_left lands exactly on a treatment unit
  #' @param eps Tolerance for floating point comparisons
  
  if (is.null(names(treatment_scores))) {
    stop("treatment_scores must have names (unit IDs)")
  }
  
  is_first_interval    <- abs(v_left - 0.0) < eps
  is_at_treatment_unit <- any(abs(T_sorted - v_left) < eps)
  
  if (is_first_interval || is_at_treatment_unit) {
    # Closed on left: [v_left, v_right]
    mask <- (treatment_scores >= v_left - eps) & (treatment_scores <= v_right + eps)
  } else {
    # Left-open: (v_left, v_right]
    mask <- (treatment_scores > v_left + eps) & (treatment_scores <= v_right + eps)
  }
  
  return(as.character(names(treatment_scores)[mask]))
}


get_controls_in_interval_dp <- function(control_scores, v_left, v_right,
                                        eps = 1e-9, left_border) {
 
  #' @param control_scores Named numeric vector of control PS (names = unit IDs)
  #' @param v_left  Left boundary
  #' @param v_right Right boundary
  #' @param eps Tolerance for floating point comparisons
  
  if (is.null(names(control_scores))) {
    stop("control_scores must have names (unit IDs)")
  }
  
  is_first_interval <- abs(v_left - left_border) < eps
  
  if (is_first_interval) {
    # Closed on left: [0, v_right]
    mask <- (control_scores >= v_left - eps) & (control_scores <= v_right + eps)
  } else {
    # Left-open: (v_left, v_right]
    mask <- (control_scores > v_left + eps) & (control_scores <= v_right + eps)
  }
  
  return(as.character(names(control_scores)[mask]))
}

get_units_in_interval <- function(scores, v_left, v_right, eps, left_border) {
  #' Extract unit IDs that fall within an interval
  #'
  #' @param scores Named numeric vector of propensity scores (names are unit IDs)
  #' @param v_left Left boundary of interval
  #' @param v_right Right boundary of interval
  #' @param eps Tolerance for floating point comparisons
  #'
  #' @return Character vector of unit IDs in the interval
  #'
  #' @details
  #' Intervals use the convention:
  #'   - First interval (v_left ≈ 0): [0, R] (closed on both ends)
  #'   - All other intervals: (L, R] (left-open, right-closed)
  #'
  #' This ensures no overlaps and no gaps between consecutive intervals.
  #
  
  # Validate inputs
  if (is.null(names(scores))) {
    stop("scores must have names (unit IDs)")
  }

  if (v_left > v_right) {
    stop(sprintf("Invalid interval: v_left (%.6f) > v_right (%.6f)", v_left, v_right))
  }

  # Determine if this is the first interval
  is_first_interval <- abs(v_left - left_border) < eps

  if (is_first_interval) {
    # First interval: [0, v_right] - closed on both ends
    mask <- (scores >= v_left - eps) & (scores <= v_right + eps)
  } else {
    # All other intervals: (v_left, v_right] - left-open, right-closed
    mask <- (scores > v_left + eps) & (scores <= v_right + eps)
  }
  
  # Extract IDs
  unit_ids <- names(scores)[mask]
  return(as.character(unit_ids))
}

# Step 1: Calculate costs with power 'p'
# p = 1 is your current version, p = 2 is your squared version
calculate_step1_cost <- function(controls, treated, p = 1) {
  sapply(treated, function(t_val) {
    min_dist <- min(abs(controls - t_val))
    return(min_dist^p)
  })
}

find_initial_partition <- function(treatment_scores, control_scores, k_bound,
                                   delta_dp, cov_df, p_sorted, treatment_col, 
                                   debug = FALSE, left_border, right_border, dist_power) {
  
  debug_local <- isTRUE(debug_glb)
  eps <- 1e-9  # Floating point tolerance
  
  n_t <- length(treatment_scores)
  if (n_t == 0) {
    return(list(error = "No treated units provided"))
  }
  T_sorted <- sort(treatment_scores)
  n_p <- length(p_sorted)
  
  # Debug prints
  if (debug_local) {
    cat(sprintf("\n=== DP Setup ===\n"))
    cat(sprintf("n_t=%d, n_c=%d, n_p=%d, K_bound=%d, delta_dp=%.6f\n",
                n_t, length(control_scores), n_p, k_bound, delta_dp))
  }
  
  # -------------------------------------------------------------------------
  # 3D DP STATE INITIALIZATION 
  # -------------------------------------------------------------------------
  
  dp <- array(list(), dim = c(n_t + 1, n_p, k_bound + 1))
  backtrack <- array(list(), dim = c(n_t + 1, n_p, k_bound + 1))
  
  # Initialize all cells with sentinels
  SENTINEL_DP <- list(cost = Inf, k_val = NA)
  
  for (idx in 1:(n_t + 1)) {
    for (v_idx in 1:n_p) {
      for (k_idx in 1:(k_bound + 1)) {
        dp[[idx, v_idx, k_idx]] <- SENTINEL_DP
      }
    }
  }
  
  # Overwrite the Base State
  dp[[1, 1, 1]] <- list(cost = 0.0, k_val = 0)
  
  # Initialize backtrack with sentinels
  SENTINEL_BACKTRACK <- list(
    prev_i_idx = -1, 
    prev_v_idx = -1, 
    prev_k_idx = -1, 
    interval_start = NA,
    interval_end = NA,
    n_treated = 0,
    n_control = 0,
    is_empty = TRUE,
    treatment_unit_ids = character(0),
    control_unit_ids = character(0)
  )
  
  for (idx in 1:(n_t + 1)) {
    for (v_idx in 1:n_p) {
      for (k_idx in 1:(k_bound + 1)) {
        backtrack[[idx, v_idx, k_idx]] <- SENTINEL_BACKTRACK
      }
    }
  }
  
  # ═══ Initialize treatment ID tracking ═══
  treatment_ids <- names(treatment_scores)
  if (is.null(treatment_ids)) {
    warning("treatment_scores has no names - using indices")
    treatment_ids <- as.character(1:n_t)
    names(treatment_scores) <- treatment_ids
  }
  
  # ═══ Initialize control ID tracking ═══
  if (is.null(names(control_scores))) {
    warning("control_scores has no names - using indices as IDs")
    names(control_scores) <- paste0("C_", seq_along(control_scores))
  }
  
  # MAIN DP LOOP (Iterating over K layers)
  for (k in 0:(k_bound - 1)) {
    cat(sprintf("\n>>> Processing K = %d\n", k))
    cat(sprintf("DP Progress: covering k=%d values\n", k))
    
    for (i in 0:n_t) {

      for (v_idx in 1:n_p) {

        prev_state <- dp[[i + 1, v_idx, k + 1]]

        # IF THIS IS THE VERY START (i=0, v=1, k=0), cost is 0
        if (i == 0 && v_idx == 1 && k == 0) {
          prev_cost_value <- 0
        } else if (is.null(prev_state) || is.infinite(prev_state$cost)) {
          next # Skip unreached or broken states
        } else {
          prev_cost_value <- prev_state$cost
        }

        v_curr <- p_sorted[v_idx]

        if (abs(v_curr - right_border) < eps) next
      
      for (v_next_idx in (v_idx + 1):n_p) {
        
        v_next <- p_sorted[v_next_idx]
        interval_width <- v_next - v_curr
        
        # 🚨 EARLY BREAK: If the interval width exceeds delta_dp, 
        # no further v_next will be valid for this v_curr!
        if (interval_width > delta_dp + eps) {
          break
        }
        
        # Get control IDs in interval
        control_ids_in_interval <- get_units_in_interval(
          scores = control_scores,
          v_left = v_curr,
          v_right = v_next,
          eps = eps,
          left_border = left_border
        )
        controls_in_interval <- length(control_ids_in_interval)
        
        next_unit_ps <- if (i < n_t) T_sorted[i + 1] else Inf
        v_curr_is_on_treatment <- abs(v_curr - next_unit_ps) < eps
        
        units_in_bin <- which(T_sorted >= v_curr - eps & T_sorted <= v_next + eps)
        
        # Only units we haven't covered yet
        new_units <- units_in_bin[units_in_bin > i]
        n_to_capture <- length(new_units)
        
        # Condition Check for Valid Match Bin (Case A)
        is_valid_bin <- (n_to_capture > 0 && 
                           new_units[1] == (i + 1) && 
                           controls_in_interval >= n_to_capture)
        
        if (is_valid_bin) {
          # ═══════════════════════════════════════════════════════
          # VALID NON-EMPTY MATCH BIN (1:1 Match in Bin)
          # ═══════════════════════════════════════════════════════
          
          treatment_ids_in_interval <- get_units_in_interval(
            scores = treatment_scores,
            v_left = v_curr,
            v_right = v_next,
            eps = eps,
            left_border = left_border
          )
          
          n_treated_to_cover <- length(new_units)  
          next_i_candidate <- i + n_treated_to_cover
          
          # Safety check
          if (next_i_candidate > n_t) {
            next
          }
          
          units_to_cover_ps <- T_sorted[(i + 1):next_i_candidate]
          controls_ps_in_interval <- control_scores[control_ids_in_interval]
          
          individual_costs <- calculate_step1_cost(controls_ps_in_interval, units_to_cover_ps, dist_power)
          
          # Sum up squared costs directly
          current_block_cost <- sum(individual_costs^2)
          
          new_cost <- prev_cost_value + current_block_cost
          k_new <- k + 1
          
          # 3D Target Indexing
          ni_idx <- next_i_candidate + 1
          nv_idx <- v_next_idx
          nk_idx <- k_new + 1
          
          if (k_new <= k_bound) {
            
            # Check existing target state in 3D DP
            target_state <- dp[[ni_idx, nv_idx, nk_idx]]
            existing_record <- backtrack[[ni_idx, nv_idx, nk_idx]]
            existing_prev_i <- if(is.null(existing_record)) -1 else existing_record$prev_i_idx
            
            current_stored_cost <- if(is.null(target_state)) Inf else target_state$cost
            
            is_better_cost   <- new_cost < current_stored_cost - eps
            is_equal_cost    <- abs(new_cost - current_stored_cost) < eps
            is_more_progress <- i > existing_prev_i 
            
            if (is.infinite(current_stored_cost) || is_better_cost || (is_equal_cost && is_more_progress)) {
              
              captured_ids <- treatment_ids[(i + 1):next_i_candidate]
              
              # Update 3D DP table
              dp[[ni_idx, nv_idx, nk_idx]] <- list(
                cost = new_cost,
                k_val = k_new
              )
              
              # Update 3D Backtrack Table
              backtrack[[ni_idx, nv_idx, nk_idx]] <- list(
                prev_i_idx = i,
                prev_v_idx = v_idx,
                prev_k_idx = k,
                interval_start = v_curr,
                interval_end   = v_next,
                n_treated      = (next_i_candidate - i),
                n_control      = controls_in_interval,
                is_empty       = FALSE,
                treatment_unit_ids = as.character(captured_ids),
                control_unit_ids   = control_ids_in_interval
              )
            }
          }
        } # End is_valid_bin
      } # End v_next_idx loop
      } 
    } 
  }    
      
  # ═══════════════════════════════════════════════════════════════════
  # 3D DP TERMINAL REACHABILITY & POST-MORTEM INSPECTION
  # ═══════════════════════════════════════════════════════════════════
  
  n_t <- length(treatment_scores)
  n_p <- length(p_sorted)
  v_term_idx <- n_p 
  i_term_idx <- n_t + 1
  
  cat("\n--- Checking 3D DP Terminal Reachability ---\n")
  
  # Check terminal slice across all K values
  terminal_slice_costs <- c()
  for (k_idx in 1:(k_bound + 1)) {
    cell <- dp[[i_term_idx, v_term_idx, k_idx]]
    if (!is.null(cell) && is.finite(cell$cost)) {
      terminal_slice_costs <- c(terminal_slice_costs, cell$cost)
    }
  }
  
  n_valid_paths <- length(terminal_slice_costs)
  cat(sprintf("Valid paths reaching Target (i=%d, PS=right_border): %d\n", n_t, n_valid_paths))
  
  # Handle Stalled DP
  if (n_valid_paths == 0) {
    max_i_reached <- 0
    for (i_idx in 1:(n_t + 1)) {
      reached_this_i <- FALSE
      for (v_idx in 1:n_p) {
        for (k_idx in 1:(k_bound + 1)) {
          if (is.finite(dp[[i_idx, v_idx, k_idx]]$cost)) {
            reached_this_i <- TRUE; break
          }
        }
        if (reached_this_i) break
      }
      if (reached_this_i) max_i_reached <- i_idx - 1
    }
    
    cat(sprintf("⚠️ DP stalled. Furthest treatment unit reached: %d / %d\n", max_i_reached, n_t))
    return(NULL) 
  }
  
  # ═══════════════════════════════════════════════════════════════════
  # FIND OPTIMAL TERMINAL STATE (Minimizing Cost, then Minimizing K)
  # ═══════════════════════════════════════════════════════════════════
  
  best_final_cost  <- Inf
  best_final_k     <- Inf
  best_final_state <- NULL
  
  for (k_val in 0:k_bound) {
    state <- dp[[i_term_idx, v_term_idx, k_val + 1]]
    
    if (!is.null(state) && is.finite(state$cost)) {
      curr_cost <- state$cost
      
      # 1. Strictly better cost
      is_better_cost <- (curr_cost < best_final_cost - eps)
      
      # 2. Equal cost (within tolerance), but fewer intervals (lower k)
      is_equal_cost_fewer_k <- (abs(curr_cost - best_final_cost) < eps) && (k_val < best_final_k)
      
      if (is_better_cost || is_equal_cost_fewer_k) {
        best_final_cost  <- curr_cost
        best_final_k     <- k_val
        best_final_state <- list(
          final_cost  = curr_cost,
          final_i_idx = i_term_idx,
          final_v_idx = v_term_idx,
          final_k_idx = k_val + 1
        )
      }
    }
  }
  
  if (is.null(best_final_state)) {
    cat(sprintf("DEBUG: Checking dp[[%d, %d, ...]]\n", i_term_idx, v_term_idx))
    stop("❌ Search loop failed to capture the terminal state.")
  }
  
  # ═══════════════════════════════════════════════════════════════════
  # POST-MORTEM TERMINAL SUMMARY TABLE
  # ═══════════════════════════════════════════════════════════════════
  
  cat("\n┌─────────────────────────────────────────────────────────────┐\n")
  cat("│         POST-MORTEM TERMINAL SLICE INSPECTION (v_end)       │\n")
  cat("├──────┬─────────────────┬────────────────────────────────────┤\n")
  cat("│  k   │      Cost       │ Selection Status                   │\n")
  cat("├──────┼─────────────────┼────────────────────────────────────┤\n")
  
  for (k_val in 0:k_bound) {
    state <- dp[[i_term_idx, v_term_idx, k_val + 1]]
    if (!is.null(state) && is.finite(state$cost)) {
      is_selected <- (k_val == (best_final_state$final_k_idx - 1))
      status_str  <- if (is_selected) "<-- CHOSEN (Min Cost, Min K)" else "Valid path"
      cat(sprintf("│ %4d │ %15.6f │ %-34s │\n", k_val, state$cost, status_str))
    } else {
      cat(sprintf("│ %4d │ %15s │ %-34s │\n", k_val, "Inf / Unreached", ""))
    }
  }
  cat("└──────┴─────────────────┴────────────────────────────────────┘\n\n")
  
  cat(sprintf("✅ Success! Optimal Final Cost: %.6f at State [k_val=%d, k_idx=%d]\n", 
              best_final_state$final_cost, 
              best_final_state$final_k_idx - 1,
              best_final_state$final_k_idx))
  
  
  # ═══════════════════════════════════════════════════════════════════
  # EXTRACT BACKTRACK AND RECONSTRUCT
  # ═══════════════════════════════════════════════════════════════════
  
  best_i_idx <- i_term_idx
  best_v_idx <- best_final_state$final_v_idx
  best_k_idx <- best_final_state$final_k_idx
  
  best_k <- best_k_idx - 1  # Convert index to value
  
  cat(sprintf("\n=== Extracting Final Backtrack Record ===\n"))
  cat(sprintf("State: [i=%d, v_idx=%d, k=%d]\n",
              n_t, best_v_idx, best_k))
  
  final_backtrack_record <- backtrack[[best_i_idx, best_v_idx, best_k_idx]]
  
  if (is.null(final_backtrack_record)) {
    cat("ERROR: Final backtrack record is NULL.\n")
    return(list(error = "Backtrack path incomplete"))
  }
  
  # Reconstruct intervals (using the updated 3D signature)
  reconstruction_result <- reconstruct_intervals_and_find_redundant(
    backtrack     = backtrack,
    p_sorted      = p_sorted,
    cov_df        = cov_df,
    id_var        = data_config$ID_VAR,
    treatment_col = treatment_col,
    n_t           = n_t,
    final_v_idx   = best_v_idx,
    final_k_idx   = best_k_idx,
    eps           = eps
  )
  
  return(list(
    cost_dist_optimal      = best_final_state$final_cost,
    k_dist_optimal         = best_k,
    intervals_dist_optimal = reconstruction_result$intervals,
    n_t_covered            = n_t
  ))
}

reconstruct_intervals_and_find_redundant <- function(backtrack, p_sorted, cov_df,
                                                     id_var, treatment_col,
                                                     n_t, final_v_idx,
                                                     final_k_idx, eps) {

  n_p <- length(p_sorted)
  k_bound <- dim(backtrack)[3] - 1

  # Current indices (1-based for R)
  current_i_idx <- n_t + 1
  current_v_idx <- final_v_idx
  current_k_idx <- final_k_idx

  cat(sprintf("\n=== Starting Reconstruction ===\n"))

  intervals_list <- list()
  step <- 0

  # Main backtracking loop
  while (current_i_idx >= 1 && current_v_idx >= 1) {
    step <- step + 1

    # --- 2. DATA ACCESS (Using 3D array index) ---
    bp_list <- backtrack[[current_i_idx, current_v_idx, current_k_idx]]

    if (is.null(bp_list)) {
      cat(sprintf("  Stopping: Found NULL at [i=%d, v=%d, k=%d]\n",
                  current_i_idx, current_v_idx, current_k_idx))
      break
    }

    # Check if this is the origin state
    if (is.null(bp_list$prev_i_idx) || bp_list$prev_i_idx == -1) {
      cat(sprintf("  Reached Origin sentinel at Step %d. Done.\n", step))
      break
    }

    # --- 3. ID EXTRACTION & TYPE CHECK ---
    # Extract treatment IDs
    treatment_ids <- bp_list$treatment_unit_ids
    if (is.null(treatment_ids)) {
      treatment_ids <- character(0)
    } else {
      treatment_ids <- as.character(unlist(treatment_ids))
    }

    # Extract control IDs
    control_ids <- bp_list$control_unit_ids
    if (is.null(control_ids)) {
      control_ids <- character(0)
    } else {
      control_ids <- as.character(unlist(control_ids))
    }

    # --- 4. PRINTING ---
    s_ps  <- as.numeric(bp_list$interval_start)
    e_ps  <- as.numeric(bp_list$interval_end)
    k_c   <- as.integer(if(length(treatment_ids) == 0) bp_list$rho else 1)
    label <- if(length(treatment_ids) == 0) "EMPTY" else "NON-EMPTY"

    cat(sprintf("Step %d: [%.4f, %.4f] %s (k_cost=%d)\n",
                step, s_ps, e_ps, label, k_c))

    # Verify control count matches attribute expectation
    n_c_stored <- as.integer(bp_list$n_control %||% 0)

    if (length(control_ids) != n_c_stored) {
      cat(sprintf("  ⚠️ Warning: Control count mismatch - stored: %d, IDs: %d\n",
                  n_c_stored, length(control_ids)))
    }

    # ═══ SAVE THE INTERVAL DATA ═══
    new_interval <- list(
      start_ps      = s_ps,
      end_ps        = e_ps,
      is_empty      = (length(treatment_ids) == 0),
      n_control     = n_c_stored,
      n_treated     = length(treatment_ids),
      k_cost        = k_c,
      treatment_ids = treatment_ids,
      control_ids   = control_ids
    )

    intervals_list[[step]] <- new_interval

    # --- 5. MOVE TO PREVIOUS STATE ---
    next_i_logical <- bp_list$prev_i_idx
    next_v_idx     <- bp_list$prev_v_idx
    next_k_logical <- bp_list$prev_k_idx

    if (next_i_logical == -1) {
      cat(sprintf("  Success: Reached the origin (Sentinel -1) at Step %d.\n", step))
      break
    }

    # Convert logical values to R indices (1-based index translation)
    current_i_idx <- next_i_logical + 1
    current_v_idx <- next_v_idx
    current_k_idx <- next_k_logical + 1
  }

  # Reverse the list to get chronological order (PS 0 -> 1)
  final_intervals_list <- rev(intervals_list)

  # Check if it's empty to avoid downstream crashes
  if (length(final_intervals_list) == 0) {
    return(list(intervals = data.frame(), redundant_ids = NULL))
  }

  # Stack elements into a clean dataframe
  intervals_df <- do.call(rbind, lapply(final_intervals_list, function(x) {
    data.frame(
      start_ps  = as.numeric(x$start_ps),
      end_ps    = as.numeric(x$end_ps),
      n_treated = as.integer(x$n_treated),
      n_control = as.integer(x$n_control),
      is_empty  = as.logical(x$is_empty),
      k_cost    = as.integer(x$k_cost),
      stringsAsFactors = FALSE
    )
  }))

  # Attach ID tracking arrays as attributes to the df
  attr(intervals_df, "treatment_unit_ids") <- lapply(final_intervals_list, `[[`, "treatment_ids")
  attr(intervals_df, "control_unit_ids")   <- lapply(final_intervals_list, `[[`, "control_ids")

  # ═══ ASSIGNMENT VERIFICATION ═══
  cat("\n=== Verification ===\n")
  total_treatments <- sum(sapply(attr(intervals_df, "treatment_unit_ids"), length))
  total_controls   <- sum(sapply(attr(intervals_df, "control_unit_ids"), length))

  cat(sprintf("Treatment units assigned: %d (expected: %d)\n", total_treatments, n_t))
  cat(sprintf("Control units assigned: %d\n", total_controls))

  # Duplicate unit validation
  all_treatment_ids <- unlist(attr(intervals_df, "treatment_unit_ids"))
  all_control_ids   <- unlist(attr(intervals_df, "control_unit_ids"))

  dup_treatments <- all_treatment_ids[duplicated(all_treatment_ids)]
  dup_controls   <- all_control_ids[duplicated(all_control_ids)]

  if (length(dup_treatments) > 0) {
    cat(sprintf("⚠️ WARNING: %d treatment units assigned to multiple bins!\n",
                length(dup_treatments)))
  }

  if (length(dup_controls) > 0) {
    cat(sprintf("⚠️ WARNING: %d control units assigned to multiple bins!\n",
                length(dup_controls)))
  }

  if (total_treatments == n_t && length(dup_treatments) == 0 && length(dup_controls) == 0) {
    cat("✅ All units accounted for with no duplicates!\n")
  }

  return(list(intervals = intervals_df))
}

calculate_ls_cost <- function(
    bounds, 
    p_sorted,
    cov_df, 
    delta_dp, 
    K_max, 
    metric, 
    cov_cols,
    id_var, 
    treatment_col, 
    data_config,
    X_all,
    S_inv, 
    intervals = NULL,
    left_border,
    right_border
) {
  
  # 5. Final Cost Calculation
  cost_result <- calculate_distance_from_bounds_units(
    bounds = bounds, 
    cov_df = cov_df, 
    metric = metric,
    id_var = id_var,
    treatment_col = treatment_col,
    X_all = X_all,                   
    S_inv = S_inv,           
    cov_cols = cov_cols,    
    delta_dp = delta_dp,     
    K_max = K_max,
    p_sorted = p_sorted,
    intervals = intervals,
    left_border = left_border,
    right_border = right_border
  )
  
  return(list(cost = cost_result$cost, k_actual = cost_result$k_actual, intervals = cost_result$intervals))
}

local_search_units <- function(
    initial_intervals, 
    initial_bounds_idx,
    treatment_ids,
    treatment_scores,
    control_scores,
    initial_cost, 
    initial_k, 
    cov_df, 
    metric, 
    K_max, 
    delta_dp, 
    data_config, 
    X_all,
    S_inv, 
    id_var, 
    treatment_col,
    p_sorted,
    left_border,
    right_border,
    max_iter = 20000, 
    max_jump = 10,
    enable_tie_breaking = TRUE,
    ...
) {
  
  # 1. Standardize Constants
  COST_EPS <- 1e-9
  
  # trace keep all improved solutions
  trace <- list()
  trace_len <- 0
  
  # Initialize state
  best_bounds_idx <- initial_bounds_idx
  best_cost <- initial_cost
  best_k <- initial_k
  best_intervals <- initial_intervals
  curr_ids <- treatment_ids
  curr_scores <- treatment_scores
  cov_cols <- data_config$ALL_COVARIATES
  
  n_ps <- length(p_sorted)
  
  cat(sprintf("\n=== Local Search (Units) Initialized ===\n"))
  cat(sprintf("Initial Cost: %.6f, K: %d\n", best_cost, best_k))
  
  # Helper for cost calculation
  cost_params <- list(
    p_sorted = p_sorted, cov_df = cov_df, delta_dp = delta_dp,
    K_max = K_max, metric = metric, cov_cols = cov_cols,
    id_var = id_var, treatment_col = treatment_col,
    data_config = data_config, X_all = X_all, S_inv = S_inv,
    left_border = left_border, right_border = right_border
  )
  
  call_ls_cost <- function(bounds) {
    do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
  }
  
  # ═══════════════════════════════════════════════════════════════════
  # MAIN LOCAL SEARCH LOOP
  # ═══════════════════════════════════════════════════════════════════
  iter <- 1
  improved <- TRUE
  best_bounds_len <- length(best_bounds_idx)
  
  while (improved && iter <= max_iter) {
    improved <- FALSE
    
    for (b in 2:(best_bounds_len - 1)) {
      
      # LAYER 1: Boundary Shift
      boundary_shifted <- FALSE
      for (dir in c(-1, 1)) {
        for (jump in 1:max_jump) {
          cand_bounds <- best_bounds_idx
          cand_bounds[b] <- cand_bounds[b] + dir * jump
          
          if (cand_bounds[b] < 0 || cand_bounds[b] > n_ps) next
          if (any(diff(cand_bounds) <= 0)) next
          
          cand_result <- call_ls_cost(cand_bounds)
          
          if (is.finite(cand_result$cost) && cand_result$cost < best_cost - COST_EPS) {
            best_cost <- cand_result$cost
            best_bounds_idx <- cand_bounds
            best_k <- cand_result$k_actual
            best_intervals <- cand_result$intervals
            
            boundary_shifted <- TRUE
            improved <- TRUE
            cat(sprintf("  ✓ Iter %d: Boundary %d shifted, cost %.6f\n", iter, b, best_cost))
            
            # Log the "improvement event"
            trace_len <- length(trace) + 1
            trace[[trace_len]] <- list(
              index = trace_len,
              cost = best_cost,
              indices = best_bounds_idx,
              intervals = best_intervals,
              k = best_k
            )
            break
          }
        }
        if (boundary_shifted) break
      }
      
      if (improved) break # Restart sweep from first interior boundary
    }
    iter <- iter + 1
  }
  
  # ═══════════════════════════════════════════════════════════════════
  # FINAL RETURN
  # ═══════════════════════════════════════════════════════════════════
  cat(sprintf("\n✓ Local Search (Units) Complete. Final Cost: %.4f\n", best_cost))
  
  return(list(
    intervals = best_intervals, # Contains updated unit IDs in attributes
    cost = best_cost,
    boundaries = best_bounds_idx,
    treatment_ids = curr_ids,
    treatment_scores = curr_scores,
    k_actual = best_k,
    has_improved = (best_cost < initial_cost - COST_EPS),
    name = "local_search_units",
    trace_log = trace
  ))
}

enhanced_local_search_units <- function(
    initial_intervals, 
    initial_bounds_idx, 
    initial_cost, 
    initial_k,
    treatment_ids, 
    treatment_scores, 
    control_scores,      
    cov_df, 
    delta_dp, 
    K_max, 
    metric, 
    X_all,
    S_inv, 
    id_var, 
    treatment_col, 
    data_config,
    p_sorted,
    left_border,
    right_border,
    max_iter = 20000, 
    max_jump = 10, 
    enable_tie_breaking = TRUE,
    ...                  
) {
  
  # 1. Standardize Constants
  COST_EPS <- 1e-9
  n_ps <- length(p_sorted)
  cov_cols <- data_config$ALL_COVARIATES
  
  # trace keep all improved solutions
  trace <- list()
  trace_len <- 0
  
  # Initialize state
  best_bounds_idx <- initial_bounds_idx
  best_cost <- initial_cost
  best_k <- initial_k
  best_intervals <- initial_intervals
  curr_ids <- treatment_ids
  curr_scores <- treatment_scores
  
  cat(sprintf("\n=== Enhanced Local Search (Shift/Split/Merge) Initialized ===\n"))
  cat(sprintf("Initial Cost: %.6f, K: %d, Bins: %d\n", 
              best_cost, best_k, length(best_bounds_idx)-1))
  
  # Cost Calculation Wrapper
  cost_params <- list(
    p_sorted = p_sorted, cov_df = cov_df, delta_dp = delta_dp,
    K_max = K_max, metric = metric, cov_cols = cov_cols,
    # n_r = n_r, fixed_redundant_unit_ids = fixed_redundant_unit_ids,
    id_var = id_var, treatment_col = treatment_col,
    data_config = data_config, X_all = X_all, S_inv = S_inv, 
    left_border = left_border, right_border = right_border
  )
  
  call_ls_cost <- function(bounds) {
    do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
  }
  
  # ═══════════════════════════════════════════════════════════════════
  # MAIN SEARCH LOOP
  # ═══════════════════════════════════════════════════════════════════
  iter <- 1
  improved <- TRUE
  
  while (improved && iter <= max_iter) {
    improved <- FALSE
    
    # --- MOVE 1: SHIFT Existing Boundaries ---
    if (length(best_bounds_idx) > 2) {
      for (b in 2:(length(best_bounds_idx) - 1)) {
        current_val <- best_bounds_idx[b]
        candidates <- unique(pmax(1, pmin(n_ps - 1, current_val + c(-max_jump:max_jump))))
        
        for (cand_val in candidates) {
          if (cand_val == current_val) next
          cand_bounds <- best_bounds_idx
          cand_bounds[b] <- cand_val
          if (any(diff(cand_bounds) <= 0)) next
          
          res <- call_ls_cost(cand_bounds)
          if (is.finite(res$cost) && res$cost < best_cost - COST_EPS) {
            best_cost <- res$cost
            best_bounds_idx <- cand_bounds
            best_k <- res$k_actual
            best_intervals <- res$intervals
            improved <- TRUE
            cat(sprintf("  ✓ Iter %d: Shift boundary %d, cost %.6f\n", iter, b, best_cost))
            
            # Log the "improvement event"
            trace_len <- length(trace) + 1
            trace[[trace_len]] <- list(
              index = trace_len,
              cost = best_cost,
              indices = best_bounds_idx,
              intervals = best_intervals,
              k = best_k
            )
            break
          }
        }
        if (improved) break
      }
    }
    
    # --- MOVE 2: SPLIT (Add Boundary) ---
    if (!improved && (length(best_bounds_idx) - 1) < K_max) {
      for (i in 1:(length(best_bounds_idx) - 1)) {
        left <- best_bounds_idx[i]; right <- best_bounds_idx[i+1]
        if (right - left > 2) {
          mid <- floor((left + right) / 2)
          cand_bounds <- sort(c(best_bounds_idx, mid))
          
          res <- call_ls_cost(cand_bounds)
          if (is.finite(res$cost) && res$cost < best_cost - COST_EPS) {
            best_cost <- res$cost
            best_bounds_idx <- cand_bounds
            best_k <- res$k_actual
            best_intervals <- res$intervals
            improved <- TRUE
            cat(sprintf("  ✓ Iter %d: Split bin %d, cost %.6f\n", iter, i, best_cost))
            
            # Log the "improvement event"
            trace_len <- length(trace) + 1
            trace[[trace_len]] <- list(
              index = trace_len,
              cost = best_cost,
              indices = best_bounds_idx,
              intervals = best_intervals,
              k = best_k
            )
            break
          }
        }
      }
    }
    
    # --- MOVE 3: MERGE (Remove Boundary) ---
    if (!improved && (length(best_bounds_idx) - 1) > 1) {
      for (b in 2:(length(best_bounds_idx) - 1)) {
        cand_bounds <- best_bounds_idx[-b]
        
        res <- call_ls_cost(cand_bounds)
        if (is.finite(res$cost) && res$cost < best_cost - COST_EPS) {
          best_cost <- res$cost
          best_bounds_idx <- cand_bounds
          best_k <- res$k_actual
          best_intervals <- res$intervals
          improved <- TRUE
          cat(sprintf("  ✓ Iter %d: Merge bins at boundary %d, cost %.6f\n", iter, b, best_cost))
          
          # Log the "improvement event"
          trace_len <- length(trace) + 1
          trace[[trace_len]] <- list(
            index = trace_len,
            cost = best_cost,
            indices = best_bounds_idx,
            intervals = best_intervals,
            k = best_k
          )
          break
        }
      }
    }
    
    iter <- iter + 1
  }
  
  cat(sprintf("\n✓ Enhanced Local Search Complete. Final Cost: %.4f\n", best_cost))
  
  return(list(
    intervals = best_intervals, # Metadata contains final control/treatment ID lists
    cost = best_cost, 
    boundaries = best_bounds_idx, 
    treatment_ids = curr_ids,
    treatment_scores = curr_scores,
    k_actual = best_k, 
    has_improved = (best_cost < initial_cost - COST_EPS), 
    name = "enhanced_ls_units",
    trace_log = trace
  ))
}

simulated_annealing_units <- function(
    initial_intervals, 
    initial_bounds_idx, 
    initial_cost, 
    initial_k,
    treatment_ids, 
    treatment_scores, 
    control_scores,      
    cov_df, 
    delta_dp, 
    K_max, 
    metric, 
    X_all,
    S_inv, 
    id_var, 
    treatment_col, 
    data_config,
    p_sorted,
    left_border,
    right_border,
    temp_init = 150, 
    temp_scale_factor = 10, 
    cooling_rate = 0.995, 
    max_iter = 20000, 
    min_temp = 0.01, 
    enable_tie_breaking = TRUE,
    ...                                     
) {
  
  # Standardize Constants
  COST_EPS <- 1e-9
  n_ps <- length(p_sorted)
  cov_cols <- data_config$ALL_COVARIATES
  
  
  # trace keep all improved solutions
  trace <- list()
  trace_len <- 0
  
  # Current State
  current_bounds <- initial_bounds_idx
  current_cost <- initial_cost
  current_intervals <- initial_intervals
  curr_ids <- treatment_ids
  curr_scores <- treatment_scores
  
  # Best State
  best_cost <- current_cost
  best_bounds <- current_bounds
  best_intervals <- current_intervals
  best_ids <- curr_ids
  best_scores <- curr_scores
  
  cost_params <- list(
    p_sorted = p_sorted, cov_df = cov_df, delta_dp = delta_dp,
    K_max = K_max, metric = metric, cov_cols = cov_cols,
    id_var = id_var, treatment_col = treatment_col,
    data_config = data_config, X_all = X_all, S_inv = S_inv, 
    left_border = left_border, right_border = right_border
  )
  
  call_ls_cost <- function(bounds) {
    do.call(calculate_ls_cost, c(list(bounds = bounds), cost_params))
  }
  
  temp <- temp_init
  cat(sprintf("\n=== Simulated Annealing Initialized ===\n"))
  
  for (iter in 1:max_iter) {
    temp_buffer <- ceiling(temp/temp_scale_factor)
    
    # 1. PERTURB
    # Randomly shift, add, or remove boundaries based on temp
    neighbor_bounds <- perturb_random_units(
      current_bounds, 
      n_ps, 
      K_max, 
      strength = max(1, temp_buffer)
    )
    
    # 2. EVALUATE
    res <- call_ls_cost(neighbor_bounds)
    neighbor_cost <- res$cost
    
    if (is.finite(neighbor_cost)) {
      delta_E <- neighbor_cost - current_cost
      
      # 3. METROPOLIS CRITERION
      if (delta_E < 0 || runif(1) < exp(-delta_E / temp)) {
        current_bounds <- neighbor_bounds
        current_cost <- neighbor_cost
        current_intervals <- res$intervals
        
        # Final check against global best
        if (current_cost < (best_cost - COST_EPS)) {
          best_cost <- current_cost
          best_bounds <- current_bounds
          best_intervals <- current_intervals
          best_ids <- curr_ids
          best_scores <- curr_scores
          cat(sprintf("  ✓ Iter %d: New Global Best! Cost %.6f (Temp: %.2f)\n",
                      iter, best_cost, temp))

          # Log the "improvement event"
          trace_len <- length(trace) + 1
          trace[[trace_len]] <- list(
            index = trace_len,
            cost = best_cost,
            indices = best_bounds,
            intervals = best_intervals,
            k = length(best_intervals),
            ids = best_ids,
            scores = best_scores
            )
        }
      }
    }
    
    # COOLING
    temp <- temp * cooling_rate
    if (temp < min_temp) break
  }
  
  cat(sprintf("\n✓ Simulated Annealing Complete. Final Cost: %.4f\n", best_cost))
  
  return(list(
    intervals = best_intervals, 
    cost = best_cost, 
    boundaries = best_bounds, 
    treatment_ids = best_ids,
    treatment_scores = best_scores,
    has_improved = (best_cost < initial_cost - COST_EPS), 
    name = "sa_units",
    trace_log = trace
  ))
}

perturb_random_units <- function(bounds_idx, n_units, K_max, strength = 5) {
  # bounds_idx: current boundary vector (sorted integers from 0 to n_units)
  # n_units: total count of treated individuals
  # K_max: maximum allowed intervals
  # strength: magnitude of the index shift
  
  # Probabilities for move types
  move_type <- sample(c("move", "split", "merge"), 1, 
                      prob = c(0.5, 0.3, 0.2))
  
  new_bounds <- bounds_idx
  
  if (move_type == "move" && length(bounds_idx) > 2) {
    # 1. MOVE: Shift an interior index boundary
    # Pick a boundary that isn't the fixed start (0) or fixed end (n_units)
    b_idx <- sample(2:(length(bounds_idx)-1), 1)
    shift <- sample(-strength:strength, 1)
    
    # Propose new index
    proposed_val <- bounds_idx[b_idx] + shift
    
    # Clip: Ensure it stays between the neighboring boundaries 
    # and within the 0 to n_units range
    new_bounds[b_idx] <- max(bounds_idx[b_idx-1] + 1, 
                             min(bounds_idx[b_idx+1] - 1, proposed_val))
    
  } else if (move_type == "split" && length(bounds_idx) < K_max + 1) {
    # 2. SPLIT: Insert a new boundary index
    interval_idx <- sample(1:(length(bounds_idx)-1), 1)
    left <- bounds_idx[interval_idx]
    right <- bounds_idx[interval_idx + 1]
    
    # We need at least one integer 'slot' available between left and right
    if (right - left > 1) {
      new_pos <- sample((left+1):(right-1), 1)
      new_bounds <- sort(c(bounds_idx, new_pos))
    }
    
  } else if (move_type == "merge" && length(bounds_idx) > 2) {
    # 3. MERGE: Remove an interior boundary index
    b_idx <- sample(2:(length(bounds_idx)-1), 1)
    new_bounds <- bounds_idx[-b_idx]
  }
  
  # Final Safety Checks
  new_bounds <- unique(sort(new_bounds))
  new_bounds <- new_bounds[new_bounds >= 0 & new_bounds <= n_units]
  
  # If perturbation resulted in an invalid state, revert to original
  if (length(new_bounds) < 2) {
    return(bounds_idx)
  }
  
  return(new_bounds)
}

calc_mahalanobis_to_group <- function(t_id, c_ids, cov_df, X_all, S_inv, id_var) {
  if (length(c_ids) == 0) return(Inf)
  
  # 1. Map IDs to row positions in the original matrix
  # This is much faster and ensures we get the dummy-coded numeric versions
  t_idx <- which(as.character(cov_df[[id_var]]) == as.character(t_id))
  c_idxs <- which(as.character(cov_df[[id_var]]) %in% as.character(c_ids))
  
  if (length(t_idx) == 0 || length(c_idxs) == 0) return(Inf)
  
  # 2. Extract numeric vectors from X_all
  # x_t is a vector (1 x K)
  x_t <- X_all[t_idx, ]
  
  # x_c_mat is a matrix (N_controls x K)
  x_c_mat <- X_all[c_idxs, , drop = FALSE]
  
  # 3. Calculate Centroid and Difference
  # colMeans on a numeric matrix is safe and fast
  centroid_c <- colMeans(x_c_mat)
  diff <- x_t - centroid_c
  
  # 4. Compute Mahalanobis distance to centroid
  # This formula applies seamlessly to both standard and robust ranked matrices
  dist <- sqrt(as.numeric(t(diff) %*% S_inv %*% diff))
  
  return(dist)
}

get_controls_for_bin_heuristic <- function(bin_idx, bounds_idx, t_scores, control_scores, delta_dp) {
  # 1. Determine the PS territory of this bin based on the boundaries
  # If it's the first bin, lower bound is the minimum possible PS (0 or min(t_scores))
  # If it's the last bin, upper bound is the maximum possible PS (1 or max(t_scores))
  
  lower_bound_ps <- if(bin_idx == 1) 0 else t_scores[bounds_idx[bin_idx]]
  upper_bound_ps <- if(bin_idx == length(bounds_idx)) 1 else t_scores[bounds_idx[bin_idx + 1]]
  
  # 2. Filter controls based on this territory
  # We look for controls that are within delta of the bin's PS range.
  # A control is 'relevant' to this bin if it could potentially match 
  # a treated unit anywhere in this bin's PS territory.
  
  feasible_mask <- control_scores >= (lower_bound_ps - delta_dp) & 
    control_scores <= (upper_bound_ps + delta_dp)
  
  # 3. Refinement: To ensure we get the 4 vs 4 split you expect visually:
  # We can further restrict to controls that actually fall INSIDE or near the bin range
  # so that the centroids for Bin 1 and Bin 2 are distinct.
  
  return(names(control_scores)[which(feasible_mask)])
}

calculate_and_format_ate <- function(matching_res, matching_non_1_1, method_label, data_config, treatment_col, gamma_val) {
  
  df <- matching_res$data_matched
  outcome_col <- data_config$OUTCOME_VAR
  covariates <- data_config$NUMERIC_COVARIATES
  
  # 1. Handle Weights Safely
  if (!"weights" %in% names(df)) {
    df$weights <- 1
  }
 
  # 2. Regression on matched data
  ate_formula <- as.formula(
    paste(outcome_col, "~", treatment_col, "+", 
          paste(covariates, collapse = " + "))
  )
  
  fit <- lm(ate_formula, data = df, weights = weights)
  summary_fit <- summary(fit)
  
  # 3. Extract the Treatment Coefficient (ATE)
  all_coeffs <- summary_fit$coefficients
  treat_row_name <- grep(treatment_col, rownames(all_coeffs), value = TRUE)[1]
  
  ate_val <- all_coeffs[treat_row_name, "Estimate"]
  se_val  <- all_coeffs[treat_row_name, "Std. Error"]
  p_val   <- all_coeffs[treat_row_name, "Pr(>|t|)"]
  
  # 4. Calculate 95% Confidence Intervals
  ci_low  <- ate_val - (1.96 * se_val)
  ci_high <- ate_val + (1.96 * se_val)
  
  # 5. Sample Sizes
  n_treated <- sum(df[[treatment_col]] == 1)
  n_control <- sum(df[[treatment_col]] == 0)
  
  # 6. Formatter Helper
  format_pval <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.0001) return("< 0.0001")
    return(sprintf("%.4f", p))
  }
  
  p_val_formated = format_pval(p_val) 
  
  # Logic for your results table
  if (method_label %in% matching_non_1_1) {
    # Stratification/Weighting Methods (CEM, Full Matching)
    gamma_shift <- NA
  } else {
    # 1:1 Methods (Optimal, Osip, etc.)
    gamma_shift <- if (!is.na(gamma_val)) round(gamma_val, 4) else 1.0
  }
  
  return(data.frame(
    method      = method_label,
    n_treated   = n_treated,
    n_control   = n_control,
    ate         = round(ate_val, 4),
    se          = round(se_val, 4),
    p_value     = p_val_formated,
    # Use the safe_gamma we just calculated
    gamma_shift = gamma_shift, 
    ci_low      = round(ci_low, 4),
    ci_high     = round(ci_high, 4),
    stringsAsFactors = FALSE
  ))
}

extract_step_row_stats <- function(m_res, trace_name, est_cost, 
                                   data_config, treatment_col, 
                                   shifting_point, matching_non_1_1) {
  if (is.null(m_res)) return(NULL)
  
  # 1. Row/Effect Stats (ATE, SE, p-value, t-stat, etc.)
  row_stats <- calculate_method_row_stats(
    df               = m_res$data_matched, 
    label            = trace_name, 
    matching_non_1_1 = matching_non_1_1,
    data_config      = data_config,
    treatment_col    = treatment_col,
    gamma_val        = shifting_point
  )
  
  # 2. Covariate Balance Stats (SMDs, Max-SMD, Mean-SMD)
  bal_stats <- generate_balance_stats(
    matched_df    = m_res$data_matched, 
    covariates    = data_config$ALL_COVARIATES, 
    treatment_col = treatment_col
  )
  
  # 3. Retention Count
  n_treated <- sum(m_res$data_matched[[treatment_col]] == 1)
  
  # Combine base stats
  full_trace_row <- cbind(
    row_stats, 
    bal_stats, 
    N_Matched = n_treated,
    est_cost  = est_cost
  )
  
  # 4. Distance Stats (Mahalanobis, KS, etc.)
  if (!is.null(m_res$match_map)) {
    dist_stats <- get_map_stats(m_res$match_map, trace_name)
    dist_stats$method <- NULL # Remove duplicate label column
    full_trace_row <- cbind(full_trace_row, dist_stats)
  }
  
  # save_solution_full(step, m_res, step_path, dist_type_name, index, method_name)
  
  return(as.data.frame(full_trace_row, stringsAsFactors = FALSE))
}

analyze_optimization_trace <- function(osip_res_step_1_row, method_name, trace, 
                                       cov_df, dist_matrix, matching_non_1_1,
                                       data_config, params, treatment_col, 
                                       outcome_var, metric_label, delta_dp, 
                                       base_dir, robust_flag, dataset_name) {
  
  message(sprintf("📈 Analyzing Optimization Trace: %d improvement steps found.", length(trace)))
  
  # res_list <- list()
  dist_type      <- if (robust_flag) "robust_mahalanobis" else "mahalanobis"
  dist_type_name <- if (robust_flag) "robust" else "strict"
  
  # 1. Define the directory path (e.g., .../optimization_traces_delta_0.1)
  
  # This is clean, safe, and handles the slashes for you
  trace_dir_path <- file.path(
    base_dir, 
    method_name, 
    "optimization_traces"
  )
  
  # 🚨 Just make sure to create the directory before saving!
  if (!dir.exists(trace_dir_path)) {
    dir.create(trace_dir_path, recursive = TRUE)
  }
  
  # Step 1-3: Loop through each improvement in the trace
  trace_rows <- lapply(seq_along(trace), function(i) {
    step <- trace[[i]]
    
    step_label = sprintf("Trace_Step_%d", i)
    full_trace_row <- NULL
    
    # This is clean, safe, and handles the slashes for you
    trace_dir_path_solutions <- file.path(
      trace_dir_path,
      "solutions"
    )
    
    # 🚨 Just make sure to create the directory before saving!
    if (!dir.exists(trace_dir_path_solutions)) {
      dir.create(trace_dir_path_solutions, recursive = TRUE)
    }
    
    # 3. Construct the filename 
    # Example: trace_Strict_Trace_Step_1.csv
    csv_path <- file.path(trace_dir_path_solutions, sprintf("solution_step_%d_%s.csv", step$index, dist_type_name))
    # csv_path    <- file.path(trace_dir_path_solutions, csv_filename)
    
    # 1. Start with the intervals (the 8 rows from your image)
    save_df <- step$intervals
    
    # 2. Extract and flatten the Unit IDs from the attributes
    # This is required for Step 3 of your comparison function
    save_df$treatment_ids <- sapply(attr(step$intervals, "treatment_unit_ids"), paste, collapse = ";")
    save_df$control_ids   <- sapply(attr(step$intervals, "control_unit_ids"), paste, collapse = ";")
    
    # 3. Add global metadata to every row
    save_df$step_index <- step$index
    save_df$total_cost <- step$cost
    
    # 4. Save to CSV
    write.csv(save_df, 
              file = csv_path, 
              row.names = FALSE)
    
    # This is clean, safe, and handles the slashes for you
    step_path <- file.path(
      trace_dir_path, 
      step_label,
      dist_type_name
    )
    
    # 🚨 Just make sure to create the directory before saving!
    if (!dir.exists(step_path)) {
      dir.create(step_path, recursive = TRUE)
    }
    
    # --- PART 1: RECONSTRUCT MATCHING ---
    # We call your existing matcher but suppress the plots/saving to keep it fast
    m_res <- create_osip_matching(
      final_intervals = step$intervals,
      cov_df          = cov_df,
      dist_matrix     = dist_matrix,
      params          = params,
      data_config     = data_config,
      treatment_col   = treatment_col,
      dist_type       = dist_type,
      dist_type_name  = dist_type_name,
      metric_label    = metric_label,
      current_delta   = delta_dp,
      base_dir        = base_dir, 
      outcome_var     = outcome_var,
      step_label      = step_label
    )
   
    shifting_point = m_res$sens$exact_threshold

    # Determine the label based on heuristic type
    is_strict <- tolower(dist_type_name) == "strict"
    base_method_label <- if(is_strict) method_labels[methods$osip_step2_strict] else method_labels[methods$osip_step2_robust]
    current_trace_name <- paste(base_method_label, step_label, sep = "_")
    
    full_trace_row <- extract_step_row_stats(
      m_res            = m_res,        
      trace_name       = current_trace_name,
      est_cost         = step$cost,           
      data_config      = data_config,
      treatment_col    = treatment_col,
      shifting_point   = shifting_point,
      matching_non_1_1 = matching_non_1_1
    )
    return(full_trace_row)
    })
  
  # Bind all heuristic trace steps into a single dataframe
  heuristic_df <- do.call(rbind, trace_rows)
  
  # Prepend the Step 1 (DP) baseline row as Row 1
  comparison_df <- rbind(osip_res_step_1_row, heuristic_df)
  
  # Calculate Efficiency Metrics
  # Efficiency is inversely proportional to the variance (SE^2)
  display_df <- comparison_df %>%
    mutate(
      t_stat = ate / se
    )
  
  # =========================================================================
  # PRINT SUMMARY
  # =========================================================================
  cat("\n\n═══════════════════════════════════════════════════════════════\n")
  cat("                 FINAL COMPARISON TABLE\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  
  # display_df <- comparison_df %>%
  #   # Select columns for console viewing
  #   dplyr::select(any_of(c("method", "n_treated", "n_control", "ate", "se", "p_value", 
  #                          "ci_low", "ci_high", "gamma_shift", "t_stat", "Max_SMD", 
  #                          "N_Imbalanced_01", "N_Imbalanced_02"))) 
  
  print(display_df)
  # =========================================================================
  # 💾 SAVE TRACE REPORT TO CSV
  # =========================================================================
  
  # 3. Construct the filename 
  # Example: trace_Strict_Trace_Step_1.csv
  csv_filename <- paste0(method_name, "_", dist_type_name, ".csv")
  full_path    <- file.path(trace_dir_path, csv_filename)
  
  # 4. Write the CSV
  # We save the full comparison_df to keep high precision for later math, 
  # but you can use display_df if you want the rounded version.
  write.csv(display_df, full_path, row.names = FALSE)
  
  inner_folder_name_comparison <- sprintf("files_delta_%.2f_comparison", delta_dp)
  full_inner_path_comparison <- file.path(trace_dir_path, inner_folder_name_comparison)
  
  message(sprintf("💾 Trace report saved to: %s", csv_filename))
  
  return(display_df)
}

extract_treatment_assignments <- function(solution) {
  # 1. Access ONLY the treatment attribute from the intervals dataframe
  # Based on image_5d7fec.png, this is 'treatment_unit_ids'
  id_list <- attr(solution$intervals, "treatment_unit_ids")
  
  if (is.null(id_list)) return(numeric(0))
  
  # Initialize a list to collect the ID-to-Bin mappings
  collection <- list()
  
  for (i in seq_along(id_list)) {
    # Extract the character vector for the i-th bin (e.g., Bin 1, 2, etc.)
    # Use [[i]] to reach inside the 'List of 8'
    current_ids <- id_list[[i]]
    
    # Check for empty bins or chr(0) as shown in image_5d6944.png
    if (length(current_ids) > 0 && !identical(current_ids, "")) {
      
      # Ensure ids are treated as a flat character vector for naming
      clean_ids <- as.character(unlist(current_ids))
      
      # Value = Bin Index (i), Name = Unit ID String
      collection[[i]] <- setNames(rep(i, length(clean_ids)), clean_ids)
    }
  }
  
  # Combine into one master dictionary
  return(do.call(c, collection))
}

compare_solutions <- function(solution_a, solution_b, name_a = "Solution A", name_b = "Solution B") {
  cat(sprintf("\n═══ COMPARISON: %s vs %s ═══\n", name_a, name_b))
  
  # --- Step 1: Compare Number of Bins ---
  n_bins_a <- length(solution_a$intervals$n_treated)
  n_bins_b <- length(solution_b$intervals$n_treated)
  
  cat(sprintf("\n[Step 1: Bin Count]\n"))
  cat(sprintf("  %s: %d bins\n", name_a, n_bins_a))
  cat(sprintf("  %s: %d bins\n", name_b, n_bins_b))
  if(n_bins_a != n_bins_b) cat("  ⚠ Note: Bin counts differ.\n")
  
  # --- Step 2: Compare Bin Boundaries ---
  cat(sprintf("\n[Step 2: Boundary Alignment]\n"))
  # We compare up to the shorter bin list to avoid out-of-bounds errors
  max_iter <- min(n_bins_a, n_bins_b)
  
  for (i in 1:max_iter) {
    bin_a <- solution_a$intervals[i, ]
    bin_b <- solution_b$intervals[i, ]
    
    l_match <- identical(bin_a$start_ps, bin_b$start_ps)
    r_match <- identical(bin_a$end_ps, bin_b$end_ps)
    
    cat(sprintf("  Bin %d:\n", i))
    cat(sprintf("    Left  (%s vs %s): %s\n", bin_a$start_ps, bin_b$start_ps, ifelse(l_match, "Identical", "DIFFERENT")))
    cat(sprintf("    Right (%s vs %s): %s\n", bin_a$end_ps, bin_b$end_ps, ifelse(r_match, "Identical", "DIFFERENT")))
  }
  
  # --- Step 3: Compare Unit Assignments (Differences Only) ---
  cat(sprintf("\n[Step 3: Unit Assignment Shifts]\n"))
  
  # Inside your compare_solutions function:
  ids_a <- extract_treatment_assignments(solution_a)
  ids_b <- extract_treatment_assignments(solution_b)
  
  # Step 3: Compare assignments
  common_ids <- intersect(names(ids_a), names(ids_b))
  diff_ids <- common_ids[ids_a[common_ids] != ids_b[common_ids]]
  
  if (length(diff_ids) > 0) {
    cat(sprintf("  Detected %d unit shifts:\n", length(diff_ids)))
    for (id in diff_ids) {
      cat(sprintf("    Unit %s: %s -> Bin %d, %s -> Bin %d\n", 
                  id, name_a, ids_a[id], name_b, ids_b[id]))
    }
  } else {
    cat("  ✓ All units are assigned to the same bins in both solutions.\n")
  } 
}

load_solution_from_csv <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # 1. Read the flat CSV
  raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # 2. Extract global metadata from the first row
  step_index <- raw_data$step_index[1]
  total_cost <- raw_data$total_cost[1]
  
  # 3. Reconstruct the Unit ID lists from the strings
  # We split by ";" and handle empty strings/NAs safely
  t_ids_list <- lapply(raw_data$treatment_ids, function(x) {
    if (is.na(x) || x == "") return(character(0))
    unlist(strsplit(as.character(x), ";"))
  })
  
  c_ids_list <- lapply(raw_data$control_ids, function(x) {
    if (is.na(x) || x == "") return(character(0))
    unlist(strsplit(as.character(x), ";"))
  })
  
  # 4. Clean the dataframe to match the original 'intervals' structure
  # Remove the flattened columns and metadata columns
  clean_intervals <- raw_data[, !(names(raw_data) %in% 
                                    c("treatment_ids", "control_ids", "step_index", "total_cost"))]
  
  # 5. Restore the attributes as seen in image_94e1ef.png
  attr(clean_intervals, "treatment_unit_ids") <- t_ids_list
  attr(clean_intervals, "control_unit_ids")   <- c_ids_list
  
  # 6. Reassemble the 'step' list object
  step_object <- list(
    index     = step_index,
    cost      = total_cost,
    intervals = clean_intervals
  )
  
  return(step_object)
}

save_solution_full <- function(step, m_res, trace_dir, dist_name, index, method_name) {
  
  # 1. Ensure the directory exists
  if (!dir.exists(trace_dir)) {
    dir.create(trace_dir, recursive = TRUE)
  }
  
  # 2. Define the dynamic file path
  file_name <- sprintf("solution_full_step_%d_%s.rds", index, tolower(dist_name))
  file_path <- file.path(trace_dir, file_name)
  
  # 3. Create a dynamic label for the comparison table (e.g., "Strict_Step_5")
  dynamic_label <- sprintf("OSIP_%s_Step_%d", method_name, index)
  
  # 4. Wrap the step so it matches the pipeline's return structure
  # This ensures data_matched and shifting_point are where they need to be
  wrapped_sol <- list(
    label          = dynamic_label,
    data_matched   = m_res$data_matched,
    match_map      = m_res$match_map,
    shifting_point = m_res$sens$exact_threshold,
    bal_object     = m_res$bal_obj,
    m_out          = m_res$m_out,
    best_method    = method_name,
    best_intervals = step$intervals
  )
  
  # 5. Save the wrapped object
  saveRDS(wrapped_sol, file = file_path)
  
  return(file_path)
}

load_solution_full <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("Full solution RDS not found at: ", file_path)
  }
  
  # Read the RDS back into R
  step_object <- readRDS(file_path)
  
  return(step_object)
}

wrap_step_for_comparison <- function(step, m_res, metric_label, best_method_name) {
  # This builds the exact structure your run_osip_pipeline returns
  wrapped_obj <- list(
    label          = metric_label,
    data_matched   = m_res$data_matched,
    match_map      = m_res$match_map,
    shifting_point = m_res$sens$exact_threshold, # Mapping your shifting_point
    bal_object     = m_res$bal_obj,
    m_out          = m_res$m_out,
    best_method    = best_method_name,
    best_intervals = step$intervals
  )
  
  return(wrapped_obj)
}

get_best_balanced_row <- function(trace_df) {
  # 1. Filter for rows where BOTH imbalance metrics are zero
  # Using column names from your MIXTAPE spreadsheet
  balanced_rows <- trace_df[trace_df$N_Imbalanced_01 == 0 & 
                              trace_df$N_Imbalanced_02 == 0, ]
  
  # 2. Safety check: if no rows are perfectly balanced, return NULL or a message
  if (nrow(balanced_rows) == 0) {
    message("No perfectly balanced steps (0/0) found in this trace.")
    return(NULL)
  }
  
  # 3. From the balanced candidates, pick the one with the MINIMUM cost
  # which corresponds to your Step 8 (Row 9) in the image
  best_row <- balanced_rows[which.min(balanced_rows$est_cost), ]
  
  return(best_row)
}

extract_step_info <- function(row) {
  method_string <- as.character(row$method)
  
  # Extract the step number (digits after the last "Step_")
  step_num <- as.numeric(gsub(".*_Step_", "", method_string))
  
  # Extract the heuristic name
  # If it's stored in a column named 'heuristic_source' from our tournament function
  # or if we need to pull it from the row attribute/index:
  
  heuristic_name <- ifelse(!is.null(row$source_heuristic), 
                           row$source_heuristic, 
                           "unknown")
  
  return(list(
    heuristic = heuristic_name,
    step      = step_num
  ))
}

get_global_best_balanced_solution <- function(base_dir_delta, dist_type_name) {
  # 1. Define the heuristics to check
  heuristics <- c("local_search", "enhanced_ls", "simulated_annealing")
  heuristic_winners <- list()
  
  for (h_name in heuristics) {
    run_list_path <- file.path(base_dir_delta, h_name, "optimization_traces")
    file_name <- sprintf("%s_%s.csv", h_name, tolower(dist_type_name))
    trace_path <- file.path(run_list_path, file_name)
    
    if (file.exists(trace_path)) {
      trace_df <- read.csv(trace_path, stringsAsFactors = FALSE)
      n_rows <- nrow(trace_df)
      
      if (n_rows < 1) {
        message(sprintf("[%s] Trace file is empty.", h_name))
        next
      }
      
      # 2. Check the LAST row first
      last_row <- trace_df[n_rows, ]
      is_last_balanced <- (last_row$N_Imbalanced_01 == 0 && last_row$N_Imbalanced_02 == 0)
      
      if (is_last_balanced) {
        # Final row is already balanced -> handled by final summary, ignore here to prevent duplication
        message(sprintf("[%s] Final solution is already balanced. Skipping to avoid duplicate counting.", h_name))
        next
      }
    
      # 3. If last row isn't balanced, look for an intermediate balanced row (rows 1 to n-1)
      if (n_rows > 1) {
        intermediate_df <- trace_df[1:(n_rows - 1), ]
        best_intermediate_row <- get_best_balanced_row(intermediate_df)
        
        if (is.null(best_intermediate_row)) {
          message(sprintf("[%s] No balanced (0/0) intermediate state found prior to termination.", h_name))
        } else {
          message(sprintf("[%s] Found candidate intermediate balanced solution (Cost: %s).", 
                          h_name, round(best_intermediate_row$est_cost, 4)))
          best_intermediate_row$source_heuristic <- h_name
          heuristic_winners[[h_name]] <- best_intermediate_row
        }
      } else {
        message(sprintf("[%s] Single-row trace with unbalanced solution.", h_name))
      }
      
    } else {
      message(sprintf("Skipping: Trace for %s not found.", h_name))
    }
  }
  
  # 4. Tournament: Pick the candidate intermediate solution with minimum cost
  if (length(heuristic_winners) == 0) {
    message("No extra intermediate balanced solutions found across any heuristics.")
    return(NULL)
  }
  
  tournament_df <- do.call(rbind, heuristic_winners)
  global_winner <- tournament_df[which.min(tournament_df$est_cost), ]
  
  message(sprintf("--> Global intermediate winner selected from: %s", global_winner$source_heuristic))
  return(global_winner)
}

get_best_rds_path <- function(base_dir, dist_type_name, winner_row) {
  # 1. Extract metadata using the extraction logic we discussed
  info <- extract_step_info(winner_row)
  
  winner_heuristic <- info$heuristic
  winner_step      <- info$step
  
  # 2. Construct the folder and file names based on your naming convention
  step_folder <- sprintf("Trace_Step_%d", winner_step)
  file_name   <- sprintf("solution_full_step_%d_%s.rds", winner_step, tolower(dist_type_name))
  
  # 3. Build the full path
  # Structure: base_dir / heuristic / optimization_traces / Trace_Step_X / dist_type / file.rds
  rds_full_path <- file.path(
    base_dir,
    winner_heuristic,
    "optimization_traces",
    step_folder,
    tolower(dist_type_name),
    file_name
  )
  
  # 4. Safety Check
  if (!file.exists(rds_full_path)) {
    warning(sprintf("RDS object not found at expected path: %s", rds_full_path))
    return(NULL)
  }
  
  return(rds_full_path)
}
