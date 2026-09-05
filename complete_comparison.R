
# -- complete_comparison.R --

source("dataset_configs.R", echo = FALSE)


get_map_stats <- function(map_obj, method_name) {
  
  # If there is no map or the map is empty, return NULL
  # rbind will ignore this NULL result entirely
  if (is.null(map_obj) || nrow(map_obj) == 0) {
    return(NULL)
  }
  
  return(data.frame(
    method       = method_name,
    Mean_Mahal   = mean(map_obj$cost, na.rm = TRUE),
    Median_Mahal = median(map_obj$cost, na.rm = TRUE),
    Max_Mahal    = max(map_obj$cost, na.rm = TRUE),
    SD_Mahal     = sd(map_obj$cost, na.rm = TRUE)
  ))
}

calculate_method_row_stats <- function(df, label, matching_non_1_1, data_config, treatment_col, gamma_val = NA) {
  
  # Use NAMED arguments to avoid positional confusion
  row <- calculate_and_format_ate(
    matching_res  = list(data_matched = df), 
    matching_non_1_1 = matching_non_1_1, 
    method_label  = label, 
    data_config   = data_config, 
    treatment_col = treatment_col, # Ensure this variable is available in scope
    gamma_val     = gamma_val      # Now R knows exactly what this is
  )
  
  # Ensure sample sizes in the table are correct
  row$n_treated <- sum(df[[treatment_col]] == 1)
  row$n_control <- sum(df[[treatment_col]] == 0)
  
  return(row)
}

# --- UTILITY 1: Extract statistics from a bal.tab object ---
get_balance_metrics <- function(bal_obj) {
  if (is.null(bal_obj)) return(NULL)
  stats <- as.data.frame(bal_obj$Balance)
  
  # Filter out distance/PS rows
  stats <- stats[!rownames(stats) %in% c("distance", "prop.score"), ]
  
  # Handle Adjusted vs Unadjusted columns
  smd_val <- coalesce(stats[["Diff.Adj"]], stats[["Diff.Un"]])
  ks_val  <- coalesce(stats[["KS.Adj"]], stats[["KS.Un"]])
  
  smd_clean <- abs(smd_val[!is.na(smd_val)])
  ks_clean  <- ks_val[!is.na(ks_val)]
  
  return(data.frame(
    Max_SMD = if(length(smd_clean) > 0) max(smd_clean) else NA,
    Mean_SMD = if(length(smd_clean) > 0) mean(smd_clean) else NA,
    Max_KS = if(length(ks_clean) > 0) max(ks_clean) else NA,
    N_Imbalanced_01 = sum(smd_clean > 0.1),
    N_Imbalanced_02 = sum(smd_clean > 0.2)
  ))
}

generate_balance_stats <- function(matched_df, covariates, treatment_col) {
  # 1. Identify active covariates (those with variation)
  has_variation <- sapply(covariates, function(cn) {
    val_count <- length(unique(na.omit(matched_df[[cn]])))
    return(val_count > 1)
  })
  active_covariates <- covariates[has_variation]
  
  # Extract weights safely (NULL if not found)
  weights_vec <- if ("weights" %in% names(matched_df)) matched_df$weights else NULL
  
  # 2. Run bal.tab with weights parameter
  bal_obj <- bal.tab(
    x = matched_df[, active_covariates, drop = FALSE], 
    treat = matched_df[[treatment_col]],
    weights = weights_vec,                # <--- CRITICAL FIX
    stats = c("m", "ks"), 
    s.d.denom = "pooled"
  )
  
  # 3. Get the 1-row dataframe of metrics
  return(get_balance_metrics(bal_obj))
}

create_bal_object <- function(res_obj, covariates, treatment_col) {
  if (is.null(res_obj)) return(NULL)
  
  matched_df <- if(is.data.frame(res_obj)) res_obj else res_obj$data_matched
  
  has_variation <- sapply(covariates, function(cn) {
    length(unique(na.omit(matched_df[[cn]]))) > 1
  })
  active_covs <- covariates[has_variation]
  
  # Pass weights if available in matched_df
  weights_vec <- if ("weights" %in% names(matched_df)) matched_df$weights else NULL
  
  return(
    bal.tab(
      x = matched_df[, active_covs, drop = FALSE],
      treat = matched_df[[treatment_col]],
      weights = weights_vec,                # <--- EXPLICIT WEIGHTS ADDED HERE
      stats = c("m", "ks"), 
      s.d.denom = "pooled"
    )
  )
}

# ==========================================
# 3. THE MAIN WRAPPER: compare_methods
# ==========================================
compare_methods <- function(data, 
                            matching_non_1_1,
                            data_config, 
                            treatment_col, 
                            method_labels,
                            methods,
                            
                            osip_step1_strict_res = NULL, 
                            osip_step1_robust_res = NULL,
                            
                            # Strict solutions step 2
                            osip_step2_strict_best_balanced = NULL,
                            osip_step2_strict_res = NULL, 
                            
                            # Robust solutions step 2
                            osip_step2_robust_best_balanced = NULL,
                            osip_step2_robust_res = NULL,
                           
                            cem_res = NULL,
                            fm_res = NULL,
                            quintile_res = NULL,
                            refined_res = NULL,
                            optimal_1_1_res = NULL,
                            card_res = NULL,
                            gm_res = NULL,
                            unmatched_data = NULL) {
  
  # Check for required columns
  if (!"ps" %in% names(data)) {
    stop("Input data must contain a 'ps' (propensity score) column.")
  }
  
  get_row <- function(df, label, gamma_val = NA) {
    # Call the external global function with all required parameters
    return (calculate_method_row_stats(
      df                = df, 
      label             = label, 
      matching_non_1_1  = matching_non_1_1, 
      data_config       = data_config, 
      treatment_col     = treatment_col, 
      gamma_val         = gamma_val
    ))
  }
  
  # =========================================================================
  # PROCESS ALL METHODS
  # =========================================================================
  res_list <- list()
  
  # --- UNMATCHED BASELINE ---
  if (!is.null(unmatched_data)) {
  
    cat("\nProcessing: Unmatched baseline...")
    
    # If it's a list (which you are passing), extract the dataframe.
    # If it's already a dataframe, use it as is.
    df_to_process <- if(is.list(unmatched_data)) unmatched_data$data_matched else unmatched_data
    
    # Ensure PS and weights are present
    if (!"ps" %in% names(df_to_process)) {
      df_to_process$ps <- data$ps[match(rownames(df_to_process), rownames(data))]
    }
    if (!"weights" %in% names(df_to_process)) {
      df_to_process$weights <- 1
    }
    
    r_unmatched <- get_row(df_to_process, method_labels[methods$unmatched])
    res_list <- append(res_list, list(r_unmatched))
  }
  
  # --- Optimal MATCHING (1:1) ---
  if (!is.null(optimal_1_1_res)) {
    cat("\nProcessing: Optimal (1:1)...")
    
    # Restore PS column if missing
    if (!"ps" %in% names(optimal_1_1_res$data_matched)) {
      optimal_1_1_res$data_matched$ps <- data$ps[match(rownames(optimal_1_1_res$data_matched), 
                                                       rownames(data))]
    }
    
    r_card <- get_row(optimal_1_1_res$data_matched, method_labels[methods$optimal_1_1], optimal_1_1_res$shifting_point)
    res_list <- append(res_list, list(r_card))
  }
  
  # --- CARDINALITY MATCHING (1:1) ---
  if (!is.null(card_res)) {
    cat("\nProcessing: Cardinality (1:1)...")
    
    # Restore PS column if missing
    if (!"ps" %in% names(card_res$data_matched)) {
      card_res$data_matched$ps <- data$ps[match(rownames(card_res$data_matched), 
                                                rownames(data))]
    }
    
    r_card <- get_row(card_res$data_matched, method_labels[methods$cardinality_1_1], card_res$shifting_point)
    res_list <- append(res_list, list(r_card))
  }
  
  # --- OSIP STEP 1 STRICT---
  if (!is.null(osip_step1_strict_res)) {
    cat("\nProcessing: Osip Strict Step 1: Baseline...")
    r_osip_step1_strict <- get_row(osip_step1_strict_res$data_matched, method_labels[methods$osip_step1_strict], osip_step1_strict_res$shifting_point)
    res_list <- append(res_list, list(r_osip_step1_strict))
  }
  
  # --- OSIP STEP 2 STRICT---
  
  # BEST OVERALL
  if (!is.null(osip_step2_strict_res)) {
    cat("\nProcessing: Osip Step 2 Strict: Optimized...")
    r_osip_step2_strict <- get_row(osip_step2_strict_res$data_matched, method_labels[methods$osip_step2_strict], osip_step2_strict_res$shifting_point)
    res_list <- append(res_list, list(r_osip_step2_strict))
  }
  
  # BEST BALANCED
  if (!is.null(osip_step2_strict_best_balanced)) {
    cat("\nProcessing: Osip Step 2 Strict Balanced: Optimized...")
    r_osip_step2_strict_balanced <- get_row(osip_step2_strict_best_balanced$data_matched, method_labels[methods$osip_step2_strict_balanced], osip_step2_strict_best_balanced$shifting_point)
    res_list <- append(res_list, list(r_osip_step2_strict_balanced))
  }
  
  # --- OSIP STEP 1 Robust---
  if (!is.null(osip_step1_robust_res)) {
    cat("\nProcessing: Osip Robust Step 1: Baseline...")
    r_osip_step1_robust <- get_row(osip_step1_robust_res$data_matched, method_labels[methods$osip_step1_robust], osip_step1_robust_res$shifting_point)
    res_list <- append(res_list, list(r_osip_step1_robust))
  }
  
  # --- OSIP STEP 2 Robust---
  
  # BEST OVERALL
  if (!is.null(osip_step2_robust_res)) {
    cat("\nProcessing: Osip Step 2 Robust: Optimized...")
    r_osip_step2_robust <- get_row(osip_step2_robust_res$data_matched, method_labels[methods$osip_step2_robust], osip_step2_robust_res$shifting_point)
    res_list <- append(res_list, list(r_osip_step2_robust))
  }
  
  # BEST BALANCED
  if (!is.null(osip_step2_robust_best_balanced)) {
    cat("\nProcessing: Osip Step 2 Robust Balanced: Optimized...")
    r_osip_step2_robust_balanced <- get_row(osip_step2_robust_best_balanced$data_matched, method_labels[methods$osip_step2_robust_balanced], osip_step2_robust_best_balanced$shifting_point)
    res_list <- append(res_list, list(r_osip_step2_robust_balanced))
  }
  
  # --- CEM ---
  if (!is.null(cem_res)) {
    cat("\nProcessing: CEM...")
    
    # CRITICAL CHECK: Verify CEM actually matched
    n_original <- nrow(data)
    n_cem <- nrow(cem_res$data_matched)
    
    if (n_cem == n_original) {
      warning("⚠️  CEM kept ALL units - matching may have failed!")
      cat(sprintf("    Original: %d units, CEM: %d units (no change)\n",
                  n_original, n_cem))
    }
   
    r_cem <- get_row(cem_res$data_matched, method_labels[methods$cem])
    res_list <- append(res_list, list(r_cem))
  }
  
  # --- Full Match ---
  if (!is.null(fm_res)) {
    cat("\nProcessing: Full Match...")
    
    # Data extraction
    df_fm <- fm_res$data_matched
    n_original <- nrow(cov_df_standardized)
    n_fm <- nrow(df_fm)
    
    # Use your get_row helper (ensure it handles weights if that's part of your metric)
    r_fm <- get_row(fm_res$data_matched, method_labels[methods$full_matching])
    res_list <- append(res_list, list(r_fm))
  }
  
  # --- GENETIC ---
  if (!is.null(gm_res)) {
    cat("\nProcessing: GENETIC...")
    r_gm <- get_row(gm_res$data_matched, method_labels[methods$genetic_1_1], gm_res$shifting_point)
    res_list <- append(res_list, list(r_gm))
  }
  
  # --- QUINTILE ---
  if (!is.null(quintile_res)) {
    cat("\nProcessing: QUINTILE...")
    r_qn <- get_row(quintile_res$data_matched, method_labels[methods$quintile], quintile_res$shifting_point)
    res_list <- append(res_list, list(r_qn))
  }
  
  # --- QUINTILE ---
  if (!is.null(refined_res)) {
    cat("\nProcessing: REFINED QUINTILE...")
    r_rqn <- get_row(refined_res$data_matched, method_labels[methods$refined_quintile], refined_res$shifting_point)
    res_list <- append(res_list, list(r_rqn))
  }
  
  # =========================================================================
  # COMBINE RESULTS
  # =========================================================================
  
  comparison_df <- do.call(rbind, res_list)
  
  # Removed the efficiency metrics everywhere
  
  # Calculate Efficiency Metrics
  # Efficiency is inversely proportional to the variance (SE^2)
  display_df <- comparison_df %>%
    mutate(
      t_stat = round(ate / se, 4)
    )
  
  # =========================================================================
  # PRINT SUMMARY
  # =========================================================================
  cat("\n\n═══════════════════════════════════════════════════════════════\n")
  cat("                 FINAL COMPARISON TABLE\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  
  
  # display_df <- comparison_df %>%
  #   # 1. Select your columns in the desired order
  #   dplyr::select(method, n_treated, n_control, ate, se, p_value, ci_low, ci_high, 
  #                 gamma_shift, t_stat)  
  
   
      
      # 3. Format rel_eff to scientific (keeping it 1.94E-06)
      # rel_eff_formatted = ifelse(is.na(rel_eff), NA, format(rel_eff, scientific = TRUE, digits = 4)),
      
      # 4. Round the rest of the numeric metrics to 4 digits
      # across(any_of(c("ate", "se", "p_value", "ci_low", "ci_high", 
      #                 "gamma_shift", "t_stat")), 
            # ~round(., 4))
    # ) %>%
    # 5. Bring back the formatted relative efficiency and clean up
    # mutate(rel_eff = rel_eff_formatted) %>%
    # dplyr::select(-rel_eff_formatted)
  
  print(display_df)
  return(display_df)
}

complete_comparison <- function(
    results_df,         # The summary table with ATE/P-values
    data,               # The original baseline data
    data_config,             # Contains ALL_COVARIATES, ID_VAR, etc.
    treatment_col,      # e.g., "treat"
    outcome_var,        # e.g., "outcome"
    unmatched_res   = NULL, # Now safe to skip
    optimal_res     = NULL, 
    cardinality_res = NULL, 
    genetic_res     = NULL, 
    cem_res         = NULL, 
    fm_res          = NULL, 
    quintile_res    = NULL, 
    refined_res     = NULL,
    
    osip_step1_strict_res = NULL,
    osip_step1_robust_res = NULL,
    
    # Strict solutions step 2
    osip_step2_strict_best_balanced = NULL,
    osip_step2_strict_res = NULL, 
    
    # Robust solutions step 2
    osip_step2_robust_best_balanced = NULL,
    osip_step2_robust_res = NULL
) {
  
  covariates <- data_config$ALL_COVARIATES
  
  # --- HELPER 2: Safely add metrics to the summary table ---
  # add_metric <- function(method_key, res_obj) {
  #   
  #   if (is.null(res_obj)) return(NULL)
  #   
  #   # Inside your complete_comparison.R or add_metric function
  #   matched_df <- res_obj$data_matched
  #   
  #   # 1. Identify which covariates actually have more than one level in this specific match
  #   # This prevents the "contrasts can be applied only to factors with 2 or more levels" error
  #   has_variation <- sapply(covariates, function(cn) {
  #     # Drop NAs then check unique values
  #     val_count <- length(unique(na.omit(matched_df[[cn]])))
  #     return(val_count > 1)
  #   })
  #   
  #   # 2. Subset your covariates list for the balance check
  #   active_covariates <- covariates[has_variation]
  #   
  #   # 3. Log if any variables were dropped from the balance table
  #   dropped_covs <- covariates[!has_variation]
  #   if (length(dropped_covs) > 0) {
  #     cat(paste0("\n[Note] Dropping constant covariates from balance check (", 
  #                res_obj$method, "): ", paste(dropped_covs, collapse = ", "), "\n"))
  #   }
  #   
  #   # 4. Now run bal.tab safely
  #   bal_obj <- bal.tab(
  #     x = matched_df[, active_covariates, drop = FALSE], 
  #     treat = matched_df[[treatment_col]],
  #     stats = c("m", "ks"), 
  #     s.d.denom = "pooled"
  #   )
  #   
  #   stats <- get_balance_metrics(bal_obj)
  #   return(data.frame(method = method_labels[[method_key]], stats))
  # }
  
  add_metric <- function(method_key, res_obj) {
    if (is.null(res_obj)) return(NULL)
    
    matched_df <- res_obj$data_matched
    
    has_variation <- sapply(covariates, function(cn) {
      val_count <- length(unique(na.omit(matched_df[[cn]])))
      return(val_count > 1)
    })
    
    active_covariates <- covariates[has_variation]
    
    dropped_covs <- covariates[!has_variation]
    if (length(dropped_covs) > 0) {
      cat(paste0("\n[Note] Dropping constant covariates from balance check (", 
                 res_obj$method, "): ", paste(dropped_covs, collapse = ", "), "\n"))
    }
    
    # Pass weights if available in matched_df
    weights_vec <- if ("weights" %in% names(matched_df)) matched_df$weights else NULL
    
    bal_obj <- bal.tab(
      x = matched_df[, active_covariates, drop = FALSE], 
      treat = matched_df[[treatment_col]],
      weights = weights_vec,                  # <--- EXPLICIT WEIGHTS ADDED HERE
      stats = c("m", "ks"), 
      s.d.denom = "pooled"
    )
    
    stats <- get_balance_metrics(bal_obj)
    return(data.frame(method = method_labels[[method_key]], stats))
  }
  
  # --- STEP 1: Build the summary tables ---
  balance_summary <- rbind(
    add_metric(methods$unmatched,         unmatched_res),
    add_metric(methods$optimal_1_1,       optimal_res),
    add_metric(methods$cem,               cem_res),
    add_metric(methods$genetic_1_1,       genetic_res),
    add_metric(methods$quintile,          quintile_res),
    add_metric(methods$refined_quintile,  refined_res),
    add_metric(methods$full_matching,     fm_res),
    add_metric(methods$cardinality_1_1,   cardinality_res),
    add_metric(methods$osip_step1_strict, osip_step1_strict_res),
    add_metric(methods$osip_step2_strict, osip_step2_strict_res),
    add_metric(methods$osip_step2_strict_balanced, osip_step2_strict_best_balanced),
    add_metric(methods$osip_step1_robust, osip_step1_robust_res),
    add_metric(methods$osip_step2_robust, osip_step2_robust_res),
    add_metric(methods$osip_step2_robust_balanced, osip_step2_robust_best_balanced)
  )

  # --- STEP 2: Build the distance metrics table ---
  # Only for methods that return a 1:1 match_map
  distance_metrics <- rbind(
    get_map_stats(optimal_res$match_map,     methods$optimal_1_1),
    get_map_stats(cardinality_res$match_map, methods$cardinality_1_1),
    get_map_stats(genetic_res$match_map,     methods$genetic_1_1),
    get_map_stats(quintile_res$match_map,    methods$quintile),
    get_map_stats(refined_res$match_map,     methods$refined_quintile),
    get_map_stats(osip_step1_strict_res$match_map,  methods$osip_step1_strict),
    get_map_stats(osip_step2_strict_res$match_map,  methods$osip_step2_strict), 
    get_map_stats(osip_step2_strict_best_balanced$match_map,  methods$osip_step2_strict_balanced),
    get_map_stats(osip_step1_robust_res$match_map,  methods$osip_step1_robust),
    get_map_stats(osip_step2_robust_res$match_map,  methods$osip_step2_robust),
    get_map_stats(osip_step2_robust_best_balanced$match_map,  methods$osip_step2_robust_balanced)
  )
  
  # Map internal keys to Display Labels
  distance_metrics$method <- method_labels[distance_metrics$method]
  
  # --- STEP 3: Sample Retention ---
  # We count treated units in the matched data (works for k:k and 1:1)
  get_retention <- function(key, res_obj) {
    if (is.null(res_obj)) return(NULL)
    n_treated <- sum(res_obj$data_matched[[treatment_col]] == 1)
    return(data.frame(method = method_labels[[key]], N_Matched = n_treated))
  }
  
  retention <- rbind(
    get_retention(methods$optimal_1_1,       optimal_res),
    get_retention(methods$cem,               cem_res),
    get_retention(methods$genetic_1_1,       genetic_res),
    get_retention(methods$full_matching,     fm_res),
    get_retention(methods$cardinality_1_1,   cardinality_res),
    get_retention(methods$quintile,          quintile_res),
    get_retention(methods$refined_quintile,  refined_res),
    get_retention(methods$osip_step1_strict, osip_step1_strict_res),
    get_retention(methods$osip_step2_strict, osip_step2_strict_res),
    get_retention(methods$osip_step2_strict_balanced, osip_step2_strict_best_balanced),
    get_retention(methods$osip_step1_robust, osip_step1_robust_res),
    get_retention(methods$osip_step2_robust, osip_step2_robust_res),
    get_retention(methods$osip_step2_robust_balanced, osip_step2_robust_best_balanced)
  )
  
  # Display results
  cat("\n╔═══════════════════════════════════════════════════════════╗\n")
  cat("║         COMPREHENSIVE METHOD COMPARISON (FINAL)           ║\n")
  cat("╚═══════════════════════════════════════════════════════════╝\n")
  cat("\n═══ COVARIATE BALANCE ═══\n"); print(balance_summary)
  cat("\n═══ MAHALANOBIS DISTANCE (Only 1:1 methods) ═══\n"); print(distance_metrics)
  cat("\n═══ SAMPLE RETENTION (Treated Units) ═══\n"); print(retention)
  
  return(list(
    metrics = balance_summary, 
    distances = distance_metrics,
    bal_objects = list(
      unmatched         = create_bal_object(unmatched_res, covariates, treatment_col),
      optimal_1_1       = create_bal_object(optimal_res, covariates, treatment_col),
      cem               = create_bal_object(cem_res, covariates, treatment_col),
      genetic_1_1       = create_bal_object(genetic_res, covariates, treatment_col),
      full_matching     = create_bal_object(fm_res, covariates, treatment_col),
      cardinality_1_1   = create_bal_object(cardinality_res, covariates, treatment_col),
      quintile          = create_bal_object(quintile_res, covariates, treatment_col),
      refined_quintile  = create_bal_object(refined_res, covariates, treatment_col),
      osip_step1_strict = create_bal_object(osip_step1_strict_res, covariates, treatment_col),
      osip_step2_strict = create_bal_object(osip_step2_strict_res, covariates, treatment_col),
      osip_step2_strict_balanced = create_bal_object(osip_step2_strict_best_balanced, covariates, treatment_col),
      osip_step1_robust = create_bal_object(osip_step1_robust_res, covariates, treatment_col),
      osip_step2_robust = create_bal_object(osip_step2_robust_res, covariates, treatment_col),
      osip_step2_robust_balanced = create_bal_object(osip_step2_robust_best_balanced, covariates, treatment_col)
    )
  ))
}

create_smd_comparison <- function(metrics_df, delta_dp, dataset_name, base_dir, method_labels = NULL) {
  
  # 1. Robust Column Handling (Case Insensitive)
  col_names_lower <- tolower(colnames(metrics_df))
  method_col_idx <- which(col_names_lower == "method")
  
  if (length(method_col_idx) > 0) {
    metrics_df$method <- as.character(metrics_df[[method_col_idx[1]]])
  } else {
    metrics_df$method <- as.character(rownames(metrics_df))
  }
  
  # 2. Crash-Proof Label Mapping
  if (!is.null(method_labels)) {
    metrics_df$method <- vapply(metrics_df$method, function(x) {
      # Check if x exists as a valid key/name in method_labels
      if (!is.null(names(method_labels)) && x %in% names(method_labels)) {
        val <- method_labels[[x]]
        if (length(val) > 0 && !is.na(val)) return(as.character(val))
      }
      return(as.character(x)) # Fallback to original string
    }, FUN.VALUE = character(1))
  }
  
  # 3. Clean Types (Guarantees no functions passed to ggplot)
  metrics_df$Max_SMD <- as.numeric(metrics_df$Max_SMD)
  metrics_df$method  <- as.character(metrics_df$method)
  
  delta_label <- sprintf("%.2f", delta_dp)
  smd_title <- sprintf("Maximum Absolute SMD by Method (Delta = %s)", delta_label)
  
  # 4. Build ggplot
  p_smd <- ggplot(metrics_df, aes(x = reorder(method, Max_SMD), y = Max_SMD, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "darkred", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", Max_SMD)), vjust = -0.5, fontface = "bold") +
    labs(
      title    = smd_title,
      subtitle = "Thresholds at 0.1 (Target) and 0.2 (Acceptable)",
      y        = "Maximum Absolute SMD",
      x        = "Matching Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title      = element_text(face = "bold", size = 14),
      axis.text.x     = element_text(angle = 45, hjust = 1) 
    )
  
  # 5. Save Output
  file_name <- sprintf("smd_bar_%s_delta_%s.png", dataset_name, delta_label)
  full_storing_path <- file.path(base_dir, file_name)
  
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }
  
  ggsave(full_storing_path, p_smd, width = 10, height = 8, dpi = 300)
  
  return(p_smd)
}

create_love_plot <- function(bal_list, delta_dp, dataset_name, base_dir, method_labels) {
  
  love_data <- data.frame()
  
  for (method_key in names(bal_list)) {
    bal_obj <- bal_list[[method_key]]
    if (is.null(bal_obj)) next
    
    stats <- as.data.frame(bal_obj$Balance)
    smd_vector <- dplyr::coalesce(stats[["Diff.Adj"]], stats[["Diff.Un"]])
    
    # 1. Check if method_labels provided an explicit override for this key
    if (method_key %in% names(method_labels)) {
      method_display_name <- method_labels[[method_key]]
    } else if (grepl("osip", method_key, ignore.case = TRUE)) {
      
      # Extract mode: "strict" or "robust"
      dist_tag <- if (grepl("robust", method_key, ignore.case = TRUE)) "robust" else "strict"
      
      # Extract step info: "step1" or "step2"
      step_str <- if (grepl("step1", method_key, ignore.case = TRUE)) "Step 1" else "Step 2"
      
      # Extract balanced suffix
      is_bal <- grepl("balanced", method_key, ignore.case = TRUE)
      bal_suffix <- if (is_bal) ", Bal" else ""
      
      # Build unique display name, e.g., "OSIP(robust, 0.15)" vs "OSIP(robust, 0.15, Bal)"
      # Or: sprintf("OSIP(%s, %s, %.2f%s)", step_str, dist_tag, delta_dp, bal_suffix)
      method_display_name <- sprintf("OSIP(%s, %.2f%s)", dist_tag, delta_dp, bal_suffix)
      
    } else {
      # Fallback for standard baseline methods (e.g., PSM, Mahalanobis, Unmatched)
      method_display_name <- method_key
    }
    
    # Build df
    df <- data.frame(
      Method    = method_display_name,
      Covariate = rownames(stats),
      SMD       = round(as.numeric(smd_vector), 4),
      stringsAsFactors = FALSE
    )
    
    love_data <- rbind(love_data, df)
  }
  
  # Pivot long format to wide format (Now guaranteed unique Method names!)
  wide_balance_df <- love_data %>%
    pivot_wider(
      names_from = Method,
      values_from = SMD
    )
  
  # Export to CSV
  file_name <- sprintf("love_plot_table_delta_%.2f.csv", delta_dp)
  full_path <- file.path(base_dir, file_name)
  write.csv(wide_balance_df, file = full_path, row.names = FALSE)
  
  # Clean Covariate names for ggplot
  love_data$Covariate <- gsub("[_\\.]\\d.*", "", love_data$Covariate)
  love_data <- love_data[!love_data$Covariate %in% c("distance", "prop.score"), ]
  
  all_methods <- unique(love_data$Method)
  love_data$Is_OSIP <- grepl("OSIP", love_data$Method)
  
  # --- 1. Dynamic Palette & Color Mapping ---
  standard_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(all_methods))
  custom_colors <- setNames(standard_palette, all_methods)
  
  # Pattern-based color assignment matching generated display names
  for (m in all_methods) {
    if (grepl("strict.*Bal", m, ignore.case = TRUE)) {
      custom_colors[m] <- "#0044BB" # Blue (Strict Balanced)
    } else if (grepl("strict", m, ignore.case = TRUE)) {
      custom_colors[m] <- "#000000" # Black (Strict)
    } else if (grepl("robust.*Bal", m, ignore.case = TRUE)) {
      custom_colors[m] <- "#009E73" # Bluish Green (Robust Balanced)
    } else if (grepl("robust", m, ignore.case = TRUE)) {
      custom_colors[m] <- "#FF0000" # Red (Robust)
    } else if (grepl("Unmatched", m, ignore.case = TRUE)) {
      custom_colors[m] <- "#999999" # Gray
    }
  }
  
  # --- 2. Dynamic Shape Mapping ---
  shape_values <- setNames(rep(16, length(all_methods)), all_methods)
  for (m in all_methods) {
    if (grepl("strict.*Bal", m, ignore.case = TRUE)) {
      shape_values[m] <- 18 # Diamond
    } else if (grepl("strict", m, ignore.case = TRUE)) {
      shape_values[m] <- 17 # Triangle
    } else if (grepl("robust.*Bal", m, ignore.case = TRUE)) {
      shape_values[m] <- 8  # Star
    } else if (grepl("robust", m, ignore.case = TRUE)) {
      shape_values[m] <- 15 # Square
    } else if (grepl("Unmatched", m, ignore.case = TRUE)) {
      shape_values[m] <- 1  # Hollow Circle
    }
  }
  
  p_love <- ggplot(love_data, aes(x = abs(SMD), y = Covariate, color = Method, shape = Method)) +
    geom_vline(xintercept = c(0, 0.1), linetype = c("solid", "dashed"), color = "gray") +
    geom_point(aes(size = Is_OSIP), alpha = 0.8, position = position_dodge(width = 0.5)) +
    theme_minimal() +
    labs(title = paste("Covariate Balance -", dataset_name),
         subtitle = paste("Delta =", delta_dp),
         x = "Absolute Standardized Mean Difference (ASMD)",
         y = "",
         color = "Method", 
         shape = "Method") +
    scale_color_manual(values = custom_colors, breaks = all_methods) +
    scale_shape_manual(values = shape_values, breaks = all_methods) +
    scale_size_manual(values = c("TRUE" = 4.5, "FALSE" = 2.5), guide = "none") +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 5)),
      shape = "legend" 
    )
  
  file_name <- sprintf("love_plot_%s_delta_%.2f.png", dataset_name, delta_dp)
  ggsave(file.path(base_dir, file_name), p_love, width = 12, height = 8, dpi = 300)
  
  return(p_love)
}

calculate_rosenbaum_gamma <- function(match_map, data_subset, outcome_var, gamma_range = seq(1, 5, by = 0.5)) {
  
  # 1. Prepare the outcome vector (ensuring alignment with the map IDs)
  outcomes <- data_subset[[outcome_var]]
  names(outcomes) <- as.character(data_subset$id)
  
  # 2. Create the Pair-Difference Matrix (Required for Rosenbaum)
  # Rows = Matched Sets, Col 1 = Treated Outcome, Col 2 = Control Outcome
  y_matrix <- matrix(NA, nrow = nrow(match_map), ncol = 2)
  
  y_matrix[, 1] <- outcomes[as.character(match_map$treated_id)]
  y_matrix[, 2] <- outcomes[as.character(match_map$control_id)]
  
  # Remove any rows with NAs (unmatched units)
  y_matrix <- y_matrix[complete.cases(y_matrix), ]
  
  # Check the average difference. If it's negative, flip the matrix 
  # so we are testing the magnitude of the effect correctly.
  avg_diff <- mean(y_matrix[, 1] - y_matrix[, 2])
  if (avg_diff < 0) {
    y_matrix <- -y_matrix
  }
  # ----------------
  
  # 3. Iterate through Gamma values
  sens_results <- data.frame(Gamma = gamma_range, P_Value_Bound = NA)
  
  for(i in 1:nrow(sens_results)) {
    g <- sens_results$Gamma[i]
    # senmv calculates the upper bound p-value for a given Gamma
    res <- senmv(y_matrix, gamma = g, method = "h")
    # Explicitly grab only the p-value to avoid the length warning
    sens_results$P_Value_Bound[i] <- res$pval
  }
  
  # Recommended Logic Flow
  # --- Logic Flow ---
  sig_indices <- which(sens_results$P_Value_Bound <= 0.05)
  
  if (length(sig_indices) == 0) {
    # Scenario A: Not significant at Gamma = 1.0
    threshold_gamma <- 1.0
    exact_threshold <- 1.0  # Initialize here so it exists for the return list
  } else {
    # Scenario B: Significant! Now find the exact interpolated point
    exact_val <- estimate_gamma_threshold(sens_results)
    
    # Fallback logic
    exact_threshold <- if (!is.na(exact_val)) exact_val else max(sens_results$Gamma)
    threshold_gamma <- exact_threshold # You can use the same value or the discrete one
  }
  
  # Now both variables are guaranteed to exist
  return(list(
    results_table = sens_results,
    threshold = threshold_gamma,
    exact_threshold = exact_threshold
  ))
}

estimate_gamma_threshold <- function(results, alpha = 0.05) {
  
  # Find first p-value above alpha
  idx <- which(results$P_Value_Bound > alpha)[1]
  
  if (is.na(idx) || idx == 1) {
    return(NA)
  }
  
  g1 <- results$Gamma[idx - 1]
  g2 <- results$Gamma[idx]
  p1 <- results$P_Value_Bound[idx - 1]
  p2 <- results$P_Value_Bound[idx]
  
  # Linear interpolation
  g_star <- g1 + (alpha - p1) * (g2 - g1) / (p2 - p1)
  
  return(g_star)
}

get_match_map <- function(m_input, dist_mat) {
  
  # CASE 1: Standard MatchIt object
  if (inherits(m_input, "matchit")) {
    matches <- m_input$match.matrix
    map <- data.frame(
      treated_id = as.character(rownames(matches)),
      control_id = as.character(matches[, 1]),
      stringsAsFactors = FALSE
    ) %>% dplyr::filter(!is.na(control_id))
    
  } else if (inherits(m_input, "optmatch")) {
    # CASE 2: The 'optmatch' factor vector
    map <- data.frame(
      id = as.character(names(m_input)),
      group = as.character(m_input),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(!is.na(group)) %>%
      mutate(is_treated = id %in% rownames(dist_mat)) %>%
      group_by(group) %>%
      summarize(
        treated_id = id[is_treated == TRUE][1],
        control_id = id[is_treated == FALSE][1],
        .groups = 'drop'
      ) %>%
      dplyr::filter(!is.na(treated_id) & !is.na(control_id))
    
  } else {
    # CASE 3: Already a dataframe (Levin)
    map <- as.data.frame(m_input)
    map$treated_id <- as.character(map$treated_id)
    map$control_id <- as.character(map$control_id)
  }
  
  # --- UPDATED SYMMETRIC LOOKUP STARTS HERE ---
  # We use a character-safe, direction-agnostic check to handle 'flipped' IDs
  map$cost <- mapply(function(t, c) {
    t_chr <- as.character(t)
    c_chr <- as.character(c)
    
    # 1. Try Standard Orientation: Treated in rows, Control in columns
    if (t_chr %in% rownames(dist_mat) && c_chr %in% colnames(dist_mat)) {
      return(dist_mat[t_chr, c_chr])
    } 
    
    # 2. Try Flipped Orientation: Control in rows, Treated in columns
    # This solves the Jobs dataset issue where IDs like "1" are columns in the matrix
    # but appear in the treated_id column of the quintile map.
    else if (c_chr %in% rownames(dist_mat) && t_chr %in% colnames(dist_mat)) {
      return(dist_mat[c_chr, t_chr])
    } 
    
    # 3. Fallback
    else { 
      return(NA) 
    }
  }, map$treated_id, map$control_id)
  # --- UPDATED SYMMETRIC LOOKUP ENDS HERE ---
  
  return(as.data.frame(map))
}

print_efficiency_rankings <- function(efficiency_vector) {
  cat("\n--- Efficiency Rankings (Most Efficient to Least) ---\n")
  
  # 1. Create a data frame for easier manipulation
  eff_df <- data.frame(
    Method = names(efficiency_vector),
    Efficiency = as.numeric(efficiency_vector)
  )
  
  # 2. Rank: Higher is better (usually Efficiency = 1/Variance)
  # If your Efficiency metric is "Higher = Better", sort descending
  eff_df <- eff_df[order(-eff_df$Efficiency), ]
  
  # 3. Add a "Ratio to Best" column to make differences intuitive
  # This shows how much efficiency you lose compared to the top performer
  best_val <- eff_df$Efficiency[1]
  eff_df$Ratio_to_Best <- eff_df$Efficiency / best_val
  
  # 4. Format for printing
  eff_df$Efficiency_Fixed <- format(eff_df$Efficiency, scientific = FALSE, digits = 8)
  eff_df$Rank <- 1:nrow(eff_df)
  
  # 5. Final Display
  for(i in 1:nrow(eff_df)) {
    cat(sprintf("%d. %-15s | Value: %s | Score: %.2f%%\n", 
                eff_df$Rank[i], 
                eff_df$Method[i], 
                eff_df$Efficiency_Fixed[i],
                eff_df$Ratio_to_Best[i] * 100))
  }
}
