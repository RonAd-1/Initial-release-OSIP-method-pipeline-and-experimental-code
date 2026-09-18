# borders_helping_functions.R

calculate_distance_from_dp_intervals_local <- function(intervals, cov_df, 
                                                       id_var, X_working, S_inv_working, ...) {
  total_cost <- 0
  n_intervals <- nrow(intervals)
  id_to_row_map <- seq_len(nrow(cov_df))
  names(id_to_row_map) <- as.character(cov_df[[id_var]])
  
  t_id_list <- attr(intervals, "treatment_unit_ids")
  c_id_list <- attr(intervals, "control_unit_ids")
  
  for (r in seq_len(n_intervals)) {
    if (intervals$is_empty[r]) next
    t_indices <- id_to_row_map[as.character(t_id_list[[r]])]
    c_indices <- id_to_row_map[as.character(c_id_list[[r]])]
    t_indices <- as.integer(t_indices[!is.na(t_indices)])
    c_indices <- as.integer(c_indices[!is.na(c_indices)])
    
    if (length(t_indices) == 0 || length(c_indices) == 0) next
    
    X_t <- X_working[t_indices, , drop = FALSE]
    X_c <- X_working[c_indices, , drop = FALSE]
    X_t_transformed <- X_t %*% S_inv_working
    t_diag <- rowSums(X_t_transformed * X_t)
    c_diag <- rowSums((X_c %*% S_inv_working) * X_c)
    cross_term <- X_t_transformed %*% t(X_c)
    dist_sq <- outer(t_diag, c_diag, "+") - 2 * cross_term
    
    min_distances <- apply(dist_sq, 1, function(x) sqrt(pmax(0, min(x, na.rm = TRUE))))
    total_cost <- total_cost + sum(min_distances)
  }
  return(total_cost)
}

create_osip_matching <- function(final_intervals,
                                 cov_df,
                                 dist_matrix,
                                 params,
                                 data_config,
                                 treatment_col,
                                 dist_type,
                                 dist_type_name,
                                 metric_label,
                                 current_delta,
                                 base_dir,
                                 outcome_var,
                                 step_label) {
                                 
  # ------------------------------------------------------------------
  # 1. Perform bin-exact matching
  # ------------------------------------------------------------------
  matching_res <- perform_osip_matching(
    data_subset      = cov_df,
    final_intervals  = final_intervals,
    params           = params,
    data_config      = data_config,
    treatment_col    = treatment_col,
    distance_metric  = dist_type
  )
  
  matching_res$data_matched <- match.data(matching_res$m_out) %>%
    rename(match_id = subclass) %>%
    mutate(match_id = as.character(match_id))
  
  # ------------------------------------------------------------------
  # 2. Build and validate match map
  #    Overwrite 'cost' with true Mahalanobis from pre-computed matrix
  # ------------------------------------------------------------------
  map_osip <- matching_res$match_map %>%
    mutate(
      cost = mapply(
        function(t, c) dist_matrix[t, c],
        as.character(treated_id),
        as.character(control_id)
      )
    )
  
  # ------------------------------------------------------------------
  # 3. Export matching report
  # ------------------------------------------------------------------
  matching_outfile <- sprintf(
    "matching_output_%s_%s_delta_%.2f.txt",
    tolower(dist_type_name), step_label, current_delta
  )
  
  export_matching_report(
    map_osip, 
    matching_res$data_matched, 
    base_dir, 
    matching_outfile
  )
  
  # ------------------------------------------------------------------
  # 4. Rosenbaum sensitivity analysis
  # ------------------------------------------------------------------
  map_sens <- calculate_rosenbaum_gamma(
    match_map   = map_osip,
    data_subset = cov_df,
    outcome_var = outcome_var
  )
  
  title <- paste0("osip ", dist_type_name, " [", step_label, "]")
  
  cat("\n========================================\n")
  cat(" ROSENBAUM SENSITIVITY:", title, "\n")
  cat("========================================\n")
  print(map_sens$results_table)
  cat("----------------------------------------\n")
  cat("The result is robust up to Gamma =", map_sens$threshold, "\n")
  cat("========================================\n\n")
  cat("Gamma shifting point by model =", map_sens$exact_threshold, "\n")
  cat("========================================\n\n")
  
  # ------------------------------------------------------------------
  # 5. Rosenbaum plot
  # ------------------------------------------------------------------
  label_slug    <- gsub(" ", "_", tolower(metric_label))
  plot_title    <- sprintf(
    "Sensitivity Analysis: %s [%s] (Delta %.2f)",
    metric_label, step_label, current_delta
  )
  sens_plot_path <- file.path(
    base_dir,
    sprintf("sensitivity_%s_%s_delta_%.2f.png",
            label_slug, step_label, current_delta)
  )
  
  plot_rosenbaum(
    results   = map_sens$results_table,
    alpha     = 0.05,
    title     = plot_title,
    save_path = sens_plot_path
  )
  
  # ------------------------------------------------------------------
  # 6. Distance distribution plot
  # ------------------------------------------------------------------
  dist_plot <- plot_match_barplot_hist(
    distances_vector = map_osip$cost,
    method_name      = sprintf("%s_%s", label_slug, step_label),
    base_dir         = base_dir
  )
  print(dist_plot)
  
  # ------------------------------------------------------------------
  # 7. Balance object
  # ------------------------------------------------------------------
  bal_obj <- cobalt::bal.tab(
    matching_res$m_out,
    un         = TRUE,
    stats      = c("m", "ks"),
    thresholds = c(m = .1)
  )
  
  # ------------------------------------------------------------------
  # 8. Return all outputs as a named list
  # ------------------------------------------------------------------
  return(list(
    data_matched = matching_res$data_matched,
    match_map    = map_osip,
    # total_cost   = total_actual_cost,
    m_out        = matching_res$m_out,
    sens         = map_sens,
    bal_obj      = bal_obj,
    dist_plot    = dist_plot
  ))
}

run_osip_pipeline <- function(metric_label, robust_flag, dp_intervals, cov_df,
                               X_working, S_inv_working, dist_matrix,
                               iter_params, data_config,
                               treatment_col, outcome_var, base_dir, delta_dp,
                               k_bound, left_border, right_border, max_val_for_plot, 
                               matching_non_1_1, dataset_name, 
                              analyze_flag = FALSE, run_step2 = FALSE, osip_res_step1 = NULL) {
                               
  cat(sprintf("\n--- Starting Pipeline for: %s ---\n", metric_label))
  
  baseline <- calculate_distance_from_dp_intervals_local(
    intervals = dp_intervals, cov_df = cov_df,
    id_var = data_config$ID_VAR, X_working = X_working, S_inv_working = S_inv_working
  )
  
  if (!is.finite(baseline)) {
    message(sprintf("[!] %s has infinite cost! Skipping...", metric_label))
    return(NULL)
  }
  
  cat(sprintf("[✓] %s Baseline Cost: %.4f\n", metric_label, baseline))
  
  dist_type      <- if (robust_flag) "robust_mahalanobis" else "mahalanobis"
  dist_type_name <- if (robust_flag) "Robust" else "Strict"
  
  visualize_osip_step_results(
    data_subset        = cov_df,
    best_intervals     = dp_intervals,
    title_suffix       = paste0("Baseline_", dist_type_name),
    method_used        = "Step 1 DP",
    data_config        = data_config,
    delta_dp           = delta_dp,
    dataset_name       = dataset_name,
    cost               = baseline,
    base_dir           = base_dir,
    treatment_col      = treatment_col
  )
  
  # ------------------------------------------------------------------
  # Resolve best intervals: Step 1 only, or Step 2 heuristics
  # ------------------------------------------------------------------
  best_intervals   <- dp_intervals   # default — overwritten below if Step 2 improves
  best_method_name <- "dp"
  step_label       <- "Step1_Baseline"
  
  if (run_step2) {
    start_time <- start_timer()
    dp_bounds  <- boundaries_from_intervals(dp_intervals, p_sorted, left_border, right_border)
    
    osip_res_step_2 <- run_osip_step2_heuristics_units(
      osip_res_step1    = osip_res_step1,
      initial_cost      = baseline,
      initial_bounds    = dp_bounds,
      initial_intervals = dp_intervals,
      cov_df            = cov_df,
      params            = params,
      data_config       = data_config,
      p_sorted          = p_sorted,
      treatment_col     = treatment_col,
      X_all             = X_working,
      S_inv             = S_inv_working,
      delta_dp          = delta_dp,
      k_bound           = k_bound,
      distance_metric   = dist_type,
      left_border       = left_border,
      right_border      = right_border,
      dist_matrix       = dist_matrix, 
      base_dir          = base_dir,
      robust_flag       = robust_flag,
      matching_non_1_1  = matching_non_1_1,
      dataset_name      = dataset_name,
      analyze_flag      = analyze_flag
    )
    end_timer(start_time, "osip_res_step_2")
    
    if (isTRUE(attr(osip_res_step_2, "is_step2"))) {
      cat(sprintf("\n[Final Decision] Step 2 Improved via: %s", attr(osip_res_step_2, "best_method")))
      best_method_name <- osip_res_step_2$name
      best_intervals   <- osip_res_step_2$intervals
      best_cost        <- osip_res_step_2$cost
      
      visualize_osip_step_results(
        data_subset        = cov_df,
        best_intervals     = best_intervals,
        title_suffix       = sprintf("Optimized_%s", dist_type_name),
        method_used        = dist_type_name,
        data_config        = data_config,
        delta_dp           = delta_dp,
        dataset_name       = dataset_name,
        cost               = best_cost,
        base_dir           = base_dir,
        treatment_col      = treatment_col
      )
    } else {
      cat("\n[Final Decision] No improvement in Step 2. Using Step 1 DP results.")
    }
    
    step_label <- "Step2_Optimized"
  }
  
  # ------------------------------------------------------------------
  # Matching + diagnostics on whichever intervals won
  # ------------------------------------------------------------------
 
  osip_matching <- create_osip_matching(
    final_intervals     = best_intervals,
    cov_df              = cov_df,
    dist_matrix         = dist_matrix,
    params              = params,
    data_config         = data_config,
    treatment_col       = treatment_col,
    dist_type           = dist_type,
    dist_type_name      = dist_type_name,
    metric_label        = metric_label,
    current_delta       = delta_dp,
    base_dir            = base_dir,
    outcome_var         = outcome_var,
    step_label          = step_label
  )
  
  return(list(
    label          = metric_label,
    data_matched   = osip_matching$data_matched,
    match_map      = osip_matching$match_map,
    shifting_point = osip_matching$sens$exact_threshold,
    bal_object     = osip_matching$bal_obj,
    m_out          = osip_matching$m_out,
    best_method    = best_method_name,
    best_intervals = best_intervals
  ))
}

execute_and_standardize <- function(matching_func, 
                                    label, 
                                    data, 
                                    config, 
                                    treatment_col,
                                    ...) {
  
  # 1. Start with the standard required parameters
  standard_args <- list(
    data = data,
    config = config,
    treatment_col = treatment_col,
    label = label
  )
  
  # 1. INTERCEPT: Separate plotting params from matching params
  all_dots <- list(...)
  
  # Define which names belong to the Wrapper/Plotting logic
  plot_param_names <- c("max_val_for_plot", "base_dir", "current_delta")
  
  # Extract plotting params
  plot_args <- all_dots[names(all_dots) %in% plot_param_names]
  
  # Extract matching params (everything else)
  matching_args <- all_dots[!(names(all_dots) %in% plot_param_names)]
  
  # 2. Build the clean argument list for the 'run' function
  standard_args <- list(data = data, config = config, 
                        treatment_col = treatment_col, label = label)
  
  # Merge only the relevant matching arguments
  func_args <- modifyList(standard_args, matching_args)
  
  cat(paste0("\n--- Executing: ", label, " ---\n"))
  
  # 1. Run the core matching logic
  res <- do.call(matching_func, func_args)
  
  # DEBUG: check cost NAs immediately out of the method function
  if (!is.null(res$match_map) && "cost" %in% names(res$match_map)) {
    cat("DEBUG - cost NAs straight out of method:", sum(is.na(res$match_map$cost)), "/", nrow(res$match_map), "\n")
    cat("DEBUG - cost class:", class(res$match_map$cost), "\n")
    cat("DEBUG - first 5 costs:", head(res$match_map$cost, 5), "\n")
  } else {
    cat("DEBUG - no cost column in res$match_map yet\n")
  }
 
  if (is.null(res) || nrow(res$data_matched) == 0) {
    cat(paste0("Warning: ", label, " failed to find matches.\n"))
    return(NULL)
  }
  
  # 2. Extract or Build the 1:1 Match Map
  # We look for match_map, then m_out$match.matrix, then generate one if possible
  # 2. Extract or Build the 1:1 Match Map (KEEP YOUR EXISTING LOGIC)
  match_map <- NULL
  if (!is.null(res$match_map)) {
    match_map <- res$match_map
  } else if (!is.null(res$m_out$match.matrix)) {
    match_map <- data.frame(
      treated_id = as.character(rownames(res$m_out$match.matrix)),
      control_id = as.character(res$m_out$match.matrix[, 1]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(control_id))
  } else if (!is.null(res$m_out_stage2$match.matrix)) {
    match_map <- data.frame(
      treated_id = as.character(rownames(res$m_out_stage2$match.matrix)),
      control_id = as.character(res$m_out_stage2$match.matrix[, 1]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(control_id))
  }
  
  # 3. Calculate Sensitivity (shifting_point)
  shifting_point <- NA
  if (!is.null(match_map)) {
    
    # This is debugging analysis
    
    cat("\n=== Match Map Debug [", label, "] ===\n")
    cat("Total pairs in match_map:", nrow(match_map), "\n")
    cat("NA control_ids:", sum(is.na(match_map$control_id)), "\n")
    cat("NA treated_ids:", sum(is.na(match_map$treated_id)), "\n")
    
    # If cost/distance column already exists at this stage
    if ("cost" %in% names(match_map)) {
      cat("NA costs:  ", sum(is.na(match_map$cost)), "\n")
      cat("Inf costs: ", sum(is.infinite(match_map$cost)), "\n")
      cat("Cost summary:\n")
      print(summary(match_map$cost))
      
      # Show the actual offending rows
      bad_rows <- match_map[!is.finite(match_map$cost), ]
      if (nrow(bad_rows) > 0) {
        cat("Offending pairs (non-finite cost):\n")
        print(bad_rows)
        
        # CHECK 1: Are these IDs present in the distance matrix at all?
        if (exists("dist_matrix")) {
          treated_in_matrix <- bad_rows$treated_id %in% rownames(dist_matrix)
          control_in_matrix <- bad_rows$control_id %in% colnames(dist_matrix)
          cat("Bad pairs - treated IDs found in dist_matrix:", sum(treated_in_matrix), "/", nrow(bad_rows), "\n")
          cat("Bad pairs - control IDs found in dist_matrix:", sum(control_in_matrix), "/", nrow(bad_rows), "\n")
        } else {
          cat("(dist_matrix not found in environment - skipping matrix lookup)\n")
        }
        
        # CHECK 2: Are NA-cost pairs clustered in specific match_id ranges?
        if ("match_id" %in% names(match_map)) {
          cat("match_id prefix distribution (finite vs NA cost):\n")
          match_map %>%
            mutate(
              cost_status = ifelse(is.finite(cost), "finite", "NA"),
              match_id_prefix = substr(match_id, 1, 3)
            ) %>%
            count(cost_status, match_id_prefix) %>%
            print()
        }
      }
    
  
      # --- NEW: Matched sample quality summary ---
      cat("\n--- Matched Sample Quality ---\n")
      cat("Pairs with measurable cost (in dist_matrix):", sum(!is.na(match_map$cost)), "\n")
      cat("Pairs used but unmeasurable (NA cost):      ", sum(is.na(match_map$cost)), "\n")
      cat("Effective high-confidence match rate:       ", 
          round(sum(!is.na(match_map$cost)) / nrow(match_map) * 100, 1), "%\n")
      
      # CHECK 3: Are the NA-cost pairs actually used downstream?
      # If your analysis data has a 'weights' or 'subclass' column, check if bad pairs appear there
      if (exists("matched_data") && is.data.frame(matched_data)) {
        bad_treated_in_analysis <- bad_rows$treated_id %in% matched_data$id
        cat("NA-cost treated units found in matched_data:", sum(bad_treated_in_analysis), "/", nrow(bad_rows), "\n")
        if (sum(bad_treated_in_analysis) > 0) {
          cat("WARNING: unmeasurable pairs are being used in downstream analysis!\n")
        }
      } else {
        cat("(matched_data not found - cannot verify downstream usage)\n")
      }
    }
    
    # --- NEW: Generate and Save Distance Distribution Plot ---
    label_slug <- gsub(" ", "_", tolower(label))
    
    dist_plot <- plot_match_barplot_hist(
      distances_vector = match_map$cost,
      method_name      = label,
      base_dir         = base_dir
    )
    
    print(dist_plot)
    
    sens_res <- calculate_rosenbaum_gamma(
      match_map = match_map, 
      data_subset = data, 
      outcome_var = config$OUTCOME_VAR
    )
    shifting_point <- sens_res$exact_threshold
    
    # Update the internal result object for consistency
    res$match_map <- match_map
    res$shifting_point <- shifting_point
  }
  else {
      cat("WARNING: match_map is NULL for", label, "\n")
  }
  
  return(res)
}

run_quintile <- function(data, config, treatment_col, label, n_subclasses = 5, ...) {
  df_internal <- as.data.frame(data)
  
  # 1. Subclassification Stage
  m_sub <- matchit(as.formula(paste(treatment_col, "~ ps")), 
                   data = df_internal, 
                   method = "subclass", 
                   subclass = n_subclasses)
  
  df_internal$subclass_id <- m_sub$subclass
  df_for_optimal <- df_internal[!is.na(df_internal$subclass_id), ]
  
  # 2. Optimal Matching within Subclasses
  match_formula <- reformulate(termlabels = config$ALL_COVARIATES, response = treatment_col)
  m_quint <- matchit(match_formula, 
                     data = df_for_optimal,
                     method = "optimal",
                     distance = "mahalanobis",
                     exact = ~ subclass_id)
  
  # 3. Create the Match Map
  mm <- m_quint$match.matrix
  match_map_quintile <- data.frame(
    treated_id = rownames(mm),
    control_id = mm[, 1],
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(control_id))
  
  # Generate distance matrix for the units in this specific matching run
  dist_matrix <- as.matrix(match_on(match_formula, data = df_for_optimal, method = "mahalanobis"))
  
  # Use high-speed matrix indexing to pull costs
  idx <- as.matrix(match_map_quintile[, c("treated_id", "control_id")])
  
  # Ensure indices exist in the matrix before assignment
  valid_idx <- idx[,1] %in% rownames(dist_matrix) & idx[,2] %in% colnames(dist_matrix)
  match_map_quintile$cost <- NA_real_
  match_map_quintile$cost[valid_idx] <- dist_matrix[idx[valid_idx, , drop=FALSE]]
  
  # 4. Finalize matched data
  matched_data <- df_for_optimal[df_for_optimal[[config$ID_VAR]] %in% 
                                   c(match_map_quintile$treated_id, match_map_quintile$control_id), ]
  
  return(list(
    data_matched = matched_data,
    match_map    = match_map_quintile, 
    m_out        = m_quint,
    method       = label,
    n_original   = nrow(df_internal),
    n_matched    = nrow(matched_data)
  ))
}

run_refined_quintile <- function(data, config, treatment_col, label, quintile_res = NULL, ...) {
  # 1. Input Handling
  if (is.null(quintile_res)) {
    stop("Refined Quintile requires the initial quintile_res object.")
  }
  
  df_internal <- quintile_res$data_matched

  # Only pull subclasses for the units that exist in df_internal
  if ("subclass_id" %in% colnames(df_internal)) {
    # Keep only units that have a valid quintile assignment
    df_internal <- df_internal[!is.na(df_internal$subclass_id), ]
    base_strata <- as.factor(df_internal$subclass_id)
  } else {
    df_internal <- df_internal[!is.na(df_internal$subclass), ]
    base_strata <- as.factor(df_internal$subclass) 
  }
  
  # Optimization: Refine strata
  z <- as.numeric(as.character(df_internal[[treatment_col]]))
  
  # --- THE ROBUST FIX: Use model.matrix to ensure X is numeric ---
  # This converts factors into dummy (0/1) variables automatically
  X_formula <- as.formula(paste("~", paste(config$ALL_COVARIATES, collapse = " + "), "-1"))
  X <- model.matrix(X_formula, data = df_internal)
  
  refined_output <- refine(
    X = X, z = z, strata = base_strata, 
    options = list(criterion = "combo", wMax = 5, minsplit = 2)
  )
  
  print("Printing refined_output details:\n")
  
  print_refined_summary <- function(ref_out) {
    cat("\n================================================\n")
    cat("   REFINED STRATIFICATION SUMMARY\n")
    cat("================================================\n")
    
    # 1. Original Strata (Quintile) Counts
    cat("\n[1] ORIGINAL COUNTS (5 QUINTILES):\n")
    # Pulling base_strata directly from the object
    original_counts <- table(Treatment = ref_out$z, Original_Quintile = ref_out$base_strata)
    print(original_counts)
    
    # 2. Refined Strata (10 Bins) Counts
    cat("\n[2] REFINED COUNTS (10 BINS):\n")
    refined_counts <- table(Treatment = ref_out$z, Refined_Bin = ref_out$refined_strata)
    print(refined_counts)
    
    # 3. Covariate Balance (using the S3 summary method)
    cat("\n[3] BALANCE COMPARISON (SMDs):\n")
    # If ref_out$summary is NULL, we use the default summary method
    sum_obj <- summary(ref_out)
    if(!is.null(sum_obj)) {
      print(sum_obj)
    } else {
      cat("Detailed SMD summary not available in object.\n")
    }
    
    # 4. Optimization Metrics
    cat("\n[4] OPTIMIZATION METRICS:\n")
    cat("LP Objective Value: ", ref_out$details$valueLP, "\n")
    cat("IP Objective Value: ", ref_out$details$valueIP, "\n")
    cat("Number of Fractional Units: ", ref_out$details$n_fracs, "\n")
    cat("================================================\n")
  }
  
  # Run the updated print
  print_refined_summary(refined_output)
  # print(refined_output$details)
  
  df_internal$refined_subclass_id <- refined_output$refined_strata
  
  # 3. Matching: Optimal Mahalanobis within new strata
  match_formula <- reformulate(termlabels = config$ALL_COVARIATES, response = treatment_col)
  m_refined_match <- matchit(
    match_formula, 
    data = df_internal,
    method = "optimal",
    distance = "mahalanobis",
    exact = ~ refined_subclass_id
  )
  
  # 4. Create the Match Map
  mm <- m_refined_match$match.matrix
  match_map_refined <- data.frame(
    treated_id = rownames(mm),
    control_id = mm[, 1],
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(control_id))
  
  # 5. NEW: Optimized Cost Calculation (Matrix Indexing)
  # Instead of re-calculating the matrix, we use match_on which is standardized
  dist_matrix <- as.matrix(match_on(match_formula, data = df_internal, method = "mahalanobis"))
  
  # Fast lookup using matrix coordinates [row, col]
  idx <- as.matrix(match_map_refined[, c("treated_id", "control_id")])
  
  match_map_refined$cost <- dist_matrix[idx]
  
  # 6. Extract Matched Subset
  # Important: Use the original ID_VAR for the subsetting to avoid row-name confusion
  matched_data <- df_internal[df_internal[[config$ID_VAR]] %in% 
                                c(match_map_refined$treated_id, match_map_refined$control_id), ]
  
  # 7. Return Standardized List
  return(list(
    data_matched   = matched_data,
    match_map      = match_map_refined,
    m_out          = m_refined_match,
    method         = label,
    shifting_point = NA,
    n_original     = nrow(df_internal),
    n_matched      = nrow(matched_data)
  ))
}

run_optimal_1_1 <- function(data, 
                            config, 
                            treatment_col, 
                            label, 
                            dist_obj = NULL, 
                            dist_matrix = NULL, 
                            ...) {
  
  # 1. Distance Input Handling
  if (is.null(dist_obj) || is.null(dist_matrix)) {
    match_formula <- reformulate(termlabels = config$ALL_COVARIATES, response = treatment_col)
    dist_obj <- match_on(match_formula, data = data, method = "mahalanobis")
    dist_matrix <- as.matrix(dist_obj)
  }
  
  # 2. Perform Optimal Matching
  time_start <- Sys.time()
  m_obj <- pairmatch(dist_obj, controls = 1, data = data)
  time_elapsed <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  
  # 3. Process Matches
  data$match_id <- as.character(m_obj) 
  matched_data <- data[!is.na(data$match_id), ]
  
  # 4. Create Match Map
  match_map_opt <- matched_data %>%
    group_by(match_id) %>%
    summarise(
      treated_id = as.character(.data[[config$ID_VAR]][which(.data[[treatment_col]] == 1)]),
      control_id = as.character(.data[[config$ID_VAR]][which(.data[[treatment_col]] == 0)]),
      .groups = 'drop'
    )
  
  # 5. Calculate Costs
  rn <- as.character(rownames(dist_matrix))
  cn <- as.character(colnames(dist_matrix))
  
  match_map_opt$cost <- mapply(function(t, c) {
    if (t %in% rn && c %in% cn) return(as.numeric(dist_matrix[t, c]))
    if (c %in% rn && t %in% cn) return(as.numeric(dist_matrix[c, t]))
    return(NA_real_) 
  }, match_map_opt$treated_id, match_map_opt$control_id)
  
  # 6. Return Standardized List
  return(list(
    data_matched = matched_data,
    time         = time_elapsed,
    m_out        = m_obj, 
    match_map    = match_map_opt,
    method       = label        # Using the explicit parameter
  ))
}

run_cem <- function(data, col_config, treatment_col, formula_cem) {

  # 1. Attempt the match using tryCatch to handle "No units matched" errors
  m.out <- tryCatch({
    matchit(formula_cem, 
            data = data, 
            method = "cem", 
            k2k = FALSE)
  }, error = function(e) {
    # If matchit fails (e.g., no units matched), return the error object or NULL
    return(NULL)
  })
  
  # 2. Check if m.out is NULL or if zero weights were assigned
  if (is.null(m.out) || sum(m.out$weights > 0) == 0) {
    cat("CEM failed: No units were matched. Skipping CEM for this run.\n")
    return(NULL) # Return NULL so the calling script can skip this method
  }
  
  # 3. If we have matches, proceed as normal
  matched_data <- match.data(m.out, data = data)
  
  # Standardize naming for your comparison function
  if ("subclass" %in% names(matched_data)) {
    matched_data <- matched_data %>% rename(match_id = subclass)
  }
  
  return(list(
    data_matched = matched_data,
    m_out = m.out
  ))
}

run_cardinality <- function(data, config, treatment_col, label = "Cardinality (1:1)", ...) {
  cat("\n--- Running Cardinality (IP Solver) ---\n")
  
  # 1. Define formula based on explicit treatment_col
  match_formula <- reformulate(termlabels = config$ALL_COVARIATES, response = treatment_col)
  
  # 2. STAGE 1: Selection (Balance Optimization)
  m_out_card <- matchit(match_formula,
                        data = data,
                        method = "cardinality", 
                        tols = 0.15,     
                        ratio = 1, 
                        time = 60,
                        solver = "glpk",
                        discard = "none")   
  
  # Get only the units chosen by the solver
  subset_indices <- !is.na(m_out_card$weights) & m_out_card$weights > 0
  balanced_subset_data <- data[subset_indices, ]
  
  if (nrow(balanced_subset_data) == 0) {
    warning("Cardinality solver failed to find a balanced subset.")
    return(NULL)
  }
  
  # 3. STAGE 2: Pairing
  cat("Pairing the balanced subset...\n")
  # Pass the label and ... to ensure run_optimal_1_1 has everything it needs
  pairing_res <- run_optimal_1_1(balanced_subset_data, config, treatment_col, label = label, ...)
  
  # 4. Standardized Return
  return(list(
    data_matched = pairing_res$data_matched,
    match_map    = pairing_res$match_map, 
    m_out        = pairing_res$m_out,     
    m_out_stage1 = m_out_card,                       
    method       = label # Use the explicit parameter
  ))
}

run_genetic <- function(data, 
                        config, 
                        treatment_col, 
                        label = "Genetic (1:1)", 
                        S_inv = NULL,     # Pass your global S_inv
                        X_all = NULL,     # Pass your global X_all
                        pop.size = 100, 
                        max.generations = 10, 
                        seed = 123, 
                        ...) {

  cat(paste0("\n--- Running ", label, " (Standardized Distances) ---\n"))
  
  start_time <- Sys.time()
  set.seed(seed)
  
  # 1. Internal Numeric Matrix for GenMatch
  # Note: GenMatch needs its own matrix to evolve weights, 
  # but we use X_all later for the standardized cost.
  X_internal <- model.matrix(
    as.formula(paste("~", paste(config$ALL_COVARIATES, collapse = " + "), "-1")), 
    data = data
  )
  
  # 2. Genetic Search
  genout <- Matching::GenMatch(
    Tr = data[[treatment_col]],
    X = X_internal,
    BalanceMatrix = X_internal,
    pop.size = pop.size,
    max.generations = max.generations,
    print.level = 0,
    ... 
  )
  
  # 3. Perform 1:1 Matching
  matchout <- Matching::Match(
    Tr = data[[treatment_col]],
    X = X_internal,
    Weight.matrix = genout,
    M = 1,
    replace = FALSE
  )
  
  # --- 4. THE WORKING FIX: Standardize Distance Calculation ---
  
  # Get the raw integer positions from the match result
  treated_pos <- as.numeric(matchout$index.treated)
  control_pos <- as.numeric(matchout$index.control)
  
  # Map positions to Actual IDs using the global ID vector
  all_ids <- as.character(data[[config$ID_VAR]])
  treated_ids_clean <- all_ids[treated_pos]
  control_ids_clean <- all_ids[control_pos]
  
  # Use the positions directly to slice the GLOBAL numeric model matrix X_all
  # This guarantees ncol(diffs) == nrow(S_inv)
  if (!is.null(X_all) && !is.null(S_inv)) {
    diffs <- X_all[treated_pos, , drop = FALSE] - X_all[control_pos, , drop = FALSE]
    standard_costs <- sqrt(rowSums((diffs %*% S_inv) * diffs))
  } else {
    # Fallback to local calculation if global objects aren't passed
    warning("Global S_inv or X_all missing. Using local Mahalanobis.")
    S_local_inv <- solve(cov(X_internal) + diag(1e-7, ncol(X_internal)))
    diffs <- X_internal[treated_pos, , drop = FALSE] - X_internal[control_pos, , drop = FALSE]
    standard_costs <- sqrt(rowSums((diffs %*% S_local_inv) * diffs))
  }
  
  # 5. Build Result Map
  match_map_gen <- data.frame(
    treated_id = treated_ids_clean,
    control_id = control_ids_clean,
    cost       = as.numeric(standard_costs),
    stringsAsFactors = FALSE
  )
  
  # 6. Extract Matched Subset
  matched_data <- data[data[[config$ID_VAR]] %in% c(treated_ids_clean, control_ids_clean), ]
  matched_data$weights <- 1.0
  
  end_time <- Sys.time()
  
  return(list(
    data_matched = matched_data,
    match_map    = match_map_gen,
    m_out        = matchout, 
    m_out_gen    = genout,   
    time         = as.numeric(difftime(end_time, start_time, units = "secs")),
    method       = label
  ))
}

run_quintile_matching <- function(data, 
                                  config = col_config, 
                                  n_subclasses = 5, 
                                  method_label = "Quintile (1:1)") {
  
  df_internal <- as.data.frame(data)
  treat_var <- config$TREATMENT_VAR
  
  # Step A: Subclassification
  m_sub <- matchit(as.formula(paste(treat_var, "~ ps")), 
                   data = df_internal, 
                   method = "subclass", 
                   subclass = n_subclasses)
  
  df_internal$subclass_id <- m_sub$subclass
  df_for_optimal <- df_internal[!is.na(df_internal$subclass_id), ]
  
  if (length(unique(df_for_optimal[[treat_var]])) < 2) {
    warning("Insufficient common support.")
    return(NULL)
  }
  
  match_formula <- reformulate(termlabels = config$ALL_COVARIATES, response = treat_var)
  
  # Step B: Optimal Matching within subclasses
  # This produces a MatchIt object with an underlying 'optmatch' structure
  m_quint <- matchit(match_formula, 
                     data = df_for_optimal,
                     method = "optimal",
                     distance = "mahalanobis",
                     exact = ~ subclass_id)
  
  # Map weights back
  df_internal$weights <- 0 
  df_internal[names(m_quint$weights), "weights"] <- m_quint$weights
  df_matched <- df_internal[df_internal$weights > 0, ]
  
  # NEW: Create the match_map immediately using the internal match.matrix
  # This follows the pattern of your other 'pairing_res' objects
  matches <- m_quint$match.matrix
  match_map <- data.frame(
    treated_id = rownames(matches),
    control_id = matches[, 1],
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(control_id))
  
  # Return standardized list
  return(list(
    data_matched = df_matched,
    match_map    = match_map,        # Explicit map
    m_out_stage1 = m_sub,            # The subclassification stage
    m_out_stage2 = m_quint,          # The optimal matching stage
    method       = method_label,
    n_original   = nrow(df_internal),
    n_matched    = nrow(df_matched)
  ))
}

run_full_match <- function(
    cov_df,
    treatment_col,
    ps_col,
    covariates,
    seed = 123
) {
  
  start_time <- Sys.time()
  set.seed(seed)
  
  # Ensure MatchIt and optmatch are loaded
  if (!requireNamespace("MatchIt", quietly = TRUE)) stop("Please install 'MatchIt'")
  if (!requireNamespace("optmatch", quietly = TRUE)) stop("Please install 'optmatch'")
  
  # Formula for matching
  formula_str <- as.formula(paste(treatment_col, "~", 
                                  paste(c(covariates, ps_col), collapse = " + ")))
  
  # Step 1: Run Full Matching
  # We use the Mahalanobis distance within the Full Match
  m.out <- MatchIt::matchit(
    formula = formula_str,
    data = cov_df,
    method = "full",
    distance = "mahalanobis"
  )
  
  # Step 2: Extract matched data
  matched_data <- MatchIt::match.data(m.out)
  
  # Step 3: Build the match_map
  # In Full Match, we map treated units to their corresponding subclass
  # For Rosenbaum Sensitivity, we need to treat subclasses as the 'pairs'
  # Note: Standard Rosenbaum works best on pairs, but can be adapted for sets.
  
  match_map <- data.frame(
    id = rownames(matched_data),
    subclass = matched_data$subclass,
    treated = matched_data[[treatment_col]],
    weights = matched_data$weights
  )
  
  end_time <- Sys.time()
  
  return(list(
    data_matched = matched_data,
    match_map    = match_map,
    m_out        = m.out,
    time         = as.numeric(difftime(end_time, start_time, units = "secs"))
  ))
}

get_imbalanced_covariates <- function(bal_obj, threshold = 0.1) {
  if (is.null(bal_obj)) return(character(0))
  
  # Extract the balance stats table
  stats <- as.data.frame(bal_obj$Balance)
  
  # Get the relevant SMD (Adjusted if matched, Unadjusted if baseline)
  smd_vector <- abs(dplyr::coalesce(stats[["Diff.Adj"]], stats[["Diff.Un"]]))
  
  # Find indices where SMD exceeds threshold
  imbalanced_idx <- which(smd_vector > threshold)
  
  # Get the rownames (covariate names)
  imbalanced_vars <- rownames(stats)[imbalanced_idx]
  
  # Clean names (remove distance/prop.score and strip dummy suffixes like you did in love_plot)
  imbalanced_vars <- imbalanced_vars[!imbalanced_vars %in% c("distance", "prop.score")]
  imbalanced_vars <- gsub("[_\\.]\\d.*", "", imbalanced_vars)
  
  return(unique(imbalanced_vars))
}

get_robust_inputs <- function(X_all) {
  X_robust <- apply(X_all, 2, rank)
  S_inv_robust <- solve(cov(X_robust) + diag(1e-7, ncol(X_robust)))
  return(list(X_all = X_robust, S_inv = S_inv_robust))
}

get_distance_matrix <- function(cov_df, data_config, treatment_col, robust = FALSE) {
  match_formula <- reformulate(data_config$ALL_COVARIATES, treatment_col)
  method_type <- if (robust) "rank_mahalanobis" else "mahalanobis"
  dist_obj <- match_on(match_formula, data = cov_df, method = method_type)
  return(as.matrix(dist_obj))
}

# Supporting Function: Plotting the Distance Distribution

analyze_covariate_distribution <- function(data, cov_name, treatment_var = "treat", 
                                           bins = 30, method_label = "OSIP", step_label = "") {
  
  # 1. Check variable type
  is_categorical <- is.factor(data[[cov_name]]) || 
    is.character(data[[cov_name]]) || 
    length(unique(na.omit(data[[cov_name]]))) <= 2
  
  # 2. Prepare data for the mirror effect
  plot_df <- data %>%
    dplyr::select(all_of(c(cov_name, treatment_var))) %>%
    dplyr::mutate(
      Group = ifelse(!!sym(treatment_var) == 1, "Treated", "Control"),
      val = !!sym(cov_name)
    )
  
  title_text <- paste(method_label, step_label, "-", cov_name)
  
  if (is_categorical) {
    # Categorical: Stick to dodge bars for clarity
    p <- ggplot(plot_df, aes(x = as.factor(val), fill = Group)) +
      geom_bar(position = "dodge", aes(y = after_stat(prop), group = Group)) +
      scale_y_continuous(labels = scales::percent) +
      labs(title = title_text, x = cov_name, y = "Percentage", fill = "Group") +
      theme_minimal() +
      scale_fill_manual(values = c("Control" = "#7FB3D5", "Treated" = "#E67E22"))
    
  } else {
    # --- MIRROR HISTOGRAM (Back-to-Back) ---
    p <- ggplot(plot_df) +
      # Treated group grows UP
      geom_histogram(data = subset(plot_df, Group == "Treated"),
                     aes(x = val, y = after_stat(density), fill = "Treated"),
                     bins = bins, alpha = 0.9) +
      # Control group grows DOWN (note the minus sign on density)
      geom_histogram(data = subset(plot_df, Group == "Control"),
                     aes(x = val, y = -after_stat(density), fill = "Control"),
                     bins = bins, alpha = 0.9) +
      # Ensure Y-axis shows absolute values (removes the minus signs from labels)
      scale_y_continuous(labels = abs) + 
      labs(title = title_text,
           subtitle = "Mirror Distribution: Treated (Up) vs Control (Down)",
           x = cov_name, y = "Density", fill = "Group") +
      theme_minimal() +
      scale_fill_manual(values = c("Control" = "#7FB3D5", "Treated" = "#E67E22")) +
      # Add a clear baseline at zero
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
  }
  
  return(p)
}
