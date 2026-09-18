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
 
  # ---------------------------------------------------------
  # PATH 1: LALONDE
  # ---------------------------------------------------------
  if (dataset_name == datasets$LALONDE) {
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
  # PATH 3: RHC 
  # ---------------------------------------------------------
    
    
  } else if (dataset_name == datasets$RHC) {
    
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
    
  # ---------------------------------------------------------
  # PATH 4: LINDNER 
  # ---------------------------------------------------------
    
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
  }
  
  # ---------------------------------------------------------
  # PATH 5: IHDP 
  # ---------------------------------------------------------
  
  else if (dataset_name == datasets$IHDP) {
 
    data("ihdp", package = "bartcs")
    
    data_full <- ihdp
    data_config <- IDHP_CONFIG
  }
  
  # ---------------------------------------------------------
  # PATH 6: NHEFS 
  # ---------------------------------------------------------
  
  else if (dataset_name == datasets$NHEFS) {
  
    data(nhefs_complete)

    data_full <- nhefs_complete
    data_config <- NHEFS_CONFIG
  }
  
  else {
    stop(sprintf("Dataset '%s' not recognized. Check datasets.", dataset_name))
  }
  
  # Ensure an ID column exists based on config
  if (!(data_config$ID_VAR %in% names(data_full))) {
    data_full[[data_config$ID_VAR]] <- 1:nrow(data_full)
  }

  return(list(
    data = data_full,
    data_config = data_config
  ))
}

sample_dataset <- function(seed_in, treat_col, data_full, dataset_name, params) {
    set.seed(seed_in)
    # We create these here so they are available for the sample() function
    treated_ids <- data_full$id[data_full[[treatment_col]] == 1]
    control_ids <- data_full$id[data_full[[treatment_col]] == 0]

    n_treated_sample <- params$SAMPLE_SIZE

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
  
  # Use the configured treatment and ALL covariates
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
  
  return(data_subset)
}

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
  # Explicitly name the vector
  names(treatment_scores) <- treatment_ids
  
  valid_mask <- !is.na(treatment_scores) & !is.nan(treatment_scores)
  treatment_scores <- treatment_scores[valid_mask]
  
  # For Controls:
  controls_df <- data_subset[data_subset[[treatment_col]] == 0, ]
  control_scores <- controls_df$ps
  # Explicitly name the vector
  names(control_scores) <- controls_df[[id_var]]
  
  valid_ctrl <- !is.na(control_scores) & !is.nan(control_scores)
  control_scores <- control_scores[valid_ctrl]
  
  if(length(treatment_scores) == 0) {
    return(list(error = "No treated units provided."))
  }
  
  return(list(treatment_scores = treatment_scores, control_scores = control_scores))
}

run_osip_step1 <- function(data_subset, params, id_var, p_sorted, 
                            treatment_col, treatment_scores, 
                            control_scores, k_bound, delta_dp, 
                            left_border, right_border, dist_power = 1) {
                            
  cat("\n--- osip Step 1: Solving DP Partition ---\n")
  
  # DEBUG: Check if 'ps' and Treatment column exist
  if(!"ps" %in% colnames(data_subset)) stop("CRITICAL: Column 'ps' missing from data_subset")
  
  start_time <- Sys.time()
  
  solution <- tryCatch({
    find_initial_partition(
      treatment_scores = treatment_scores, 
      control_scores = control_scores, 
      k_bound  = k_bound, 
      delta_dp  = delta_dp, 
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
    
    # Execute
    raw_res <- run_search(method, defaults, overrides, opt_funcs)
    
    # Process & Store
    final_results[[method]] <- process_and_report_search(
      ls_res = raw_res, 
      search_name = method, 
      initial_cost = initial_cost,
      initial_intervals = initial_intervals,
      initial_k = length(initial_intervals)
    )
  }

  # --- Inside run_osip_step2_heuristics ---
  # vapply ensures we get a logical vector and throws an error if a list is returned
  improved_mask <- vapply(final_results, function(res) {
    # Handle potential NULLs or non-logical status values safely
    isTRUE(res$status) 
  }, logical(1))
  
  if (!any(improved_mask)) {
    # Return a "Null-Object" list instead of literal NULL
      # 1. Create a basic list 
      fail_res <- list(data_matched = NULL, name = "No Improvement")
      
      # 2. Assign the same attribute you use in the success case
      attr(fail_res, "is_step2") <- FALSE
      
      return(fail_res)
    }
  
  # 2. Extract ONLY the results that improved
  improved_results <- final_results[improved_mask]
  
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

prepare_standardized_data <- function(data_subset, data_config, treatment_col) {
  cat("\n--- Preparing data for Step 2 ---\n")
  cat("Standardizing covariates using pooled within-group SDs...\n")
  temp_df <- data_subset
  covariates <- data_config$NUMERIC_COVARIATES
  
  # Ensure Treatment is numeric
  temp_df[[treatment_col]] <- ensure_numeric_col(temp_df, treatment_col)
  
  # --- FIX 1: Force all target covariates to be numeric before checking ---
  for (var in covariates) {
    # as.character then as.numeric handles factors correctly
    temp_df[[var]] <- ensure_numeric_col(temp_df, var)
  }
    
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
  
  cat("Standardization complete!.\n")
  return (cov_df_standardized)
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

  # This section creates the 'match_map' (mm) 
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
  

  # Extract matched data for ATE calculations
  data_matched <- match.data(m_out)

  return(list(
    m_out = m_out,
    data_matched = data_matched,
    match_map = match_map_osip, # Returning this ensures plots have data
    method_used = "Optimal (Exact Bin-Matching)"
  ))
}

apply_sampling <- function(data_full, dataset_name, data_config, datasets, 
                           params, treatment_col) {
  
  # In general, we should put this 25 as a parameter. But, as we aren't using 
  # sampling we keep it like that
  set.seed(25)
 
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
  
  # Convert binary factors ("0", "1") to numeric (0, 1)
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
    warning("Cannot export CSV: Provided df is NULL or empty.")
    return(FALSE)
  }

  # 2. Clone input dataframe
  export_df <- df

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
      warning("Step 2 requires 'step2_mode' to be 'strict' or 'robust'. Defaulting to 'strict'.")
      step2_mode <- "strict"
    }

    filename <- sprintf("step2_heuristic_%s_%s_partition_delta_%.2f.csv",
                        heuristic_name, step2_mode, delta_val)
  } else {
    stop("nvalid step provided. 'step' must be 1 or 2.")
  }

  csv_file_path <- file.path(base_dir, filename)
  write.csv(export_df, file = csv_file_path, row.names = FALSE, quote = FALSE)
  cat(sprintf("Partition CSV saved to > %s\n\n", csv_file_path))
  return(TRUE)
}

get_max_k <- function(delta, p0, pk) {
  ceiling_val <- ceiling((pk - p0) / delta)
  k_max <- 2 * ceiling_val
  return(k_max)
}

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
