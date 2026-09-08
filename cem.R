cem_ps <- function(data, 
                   config = col_config,
                   custom_endpoints = NULL,
                   method_label = "CEM") {
  
  # 1. Force to base data.frame to stop tibble subsetting issues
  df_internal <- as.data.frame(data)
  
  # 2. Extract strictly what is needed for matching
  match_cols <- c(config$TREATMENT_VAR, "ps")
  df_for_cem <- df_internal[, match_cols, drop = FALSE]
  
  # 3. Setup the cutpoints list
  cp_list <- NULL  # Default: let CEM choose bins automatically
  
  if (!is.null(custom_endpoints)) {
    # Ensure PS boundaries are unique and sorted
    ends <- sort(unique(as.numeric(custom_endpoints)))
    
    # Validate endpoints
    if (length(ends) < 2) {
      warning("custom_endpoints must have at least 2 values. Using automatic binning.")
    } else if (min(ends) > min(df_for_cem$ps) || max(ends) < max(df_for_cem$ps)) {
      warning(sprintf("Endpoints [%.3f, %.3f] don't cover PS range [%.3f, %.3f]",
                      min(ends), max(ends), min(df_for_cem$ps), max(df_for_cem$ps)))
    } else {
      cp_list <- list(ps = ends)
    }
  }
  
  # 4. Run CEM on the minimized data frame
  cat(sprintf("Running CEM matching (%s)...\n", method_label))
  cat(sprintf("  Input: %d units (T=%d, C=%d)\n",
              nrow(df_for_cem),
              sum(df_for_cem[[config$TREATMENT_VAR]] == 1),
              sum(df_for_cem[[config$TREATMENT_VAR]] == 0)))
  
  if (!is.null(cp_list)) {
    cat(sprintf("  Using %d custom bins: [%.3f, %.3f, ..., %.3f]\n",
                length(cp_list$ps) - 1,
                cp_list$ps[1],
                cp_list$ps[2],
                cp_list$ps[length(cp_list$ps)]))
  } else {
    cat("  Using automatic binning (Scott's rule)\n")
  }
  
  # Always drop unmatched units for proper CEM matching
  # The cem package has conflicting behavior between 'drop' and 'keep.all'
  # Best practice: ONLY use 'drop' parameter
  cem_out <- cem::cem(
    treatment = config$TREATMENT_VAR,
    data = df_for_cem,
    cutpoints = cp_list,
    drop = "both"  # Always drop unmatched from both treatment and control
    # Do NOT specify keep.all - let it use default (FALSE)
  )
  
  # 5. Diagnostic output
  n_matched <- sum(cem_out$matched)
  n_dropped <- nrow(df_for_cem) - n_matched
  
  cat(sprintf("  Matched: %d units (%.1f%% kept)\n",
              n_matched, 100 * n_matched / nrow(df_for_cem)))
  cat(sprintf("  Dropped: %d units\n", n_dropped))
  
  # Check weights
  weights_nonzero <- cem_out$w[cem_out$w > 0]
  if (length(weights_nonzero) > 0) {
    cat(sprintf("  Weights: min=%.3f, max=%.3f, mean=%.3f, unique=%d\n",
                min(weights_nonzero),
                max(weights_nonzero),
                mean(weights_nonzero),
                length(unique(weights_nonzero))))
  }
  
  # Check balance metrics from CEM
  if (!is.null(cem_out$imbalance)) {
    cat(sprintf("  L1 distance: %.4f\n", cem_out$imbalance$L1$L1))
  }
  
  # 6. Attach weights back to the FULL original data frame
  df_internal$weights <- cem_out$w
  
  # 7. Filter to matched units only
  df_matched <- df_internal[df_internal$weights > 0, ]
  
  # CRITICAL CHECK: Verify we actually changed something
  if (nrow(df_matched) == nrow(df_internal) && 
      all(df_matched$weights == 1)) {
    warning("⚠️  CEM did not drop any units or assign differential weights!")
    warning("    This suggests matching failed or all units were trivially matched.")
    warning("    Result will be identical to 'Unmatched' baseline.")
  }
  
  cat(sprintf("  Final matched data: %d units\n\n", nrow(df_matched)))
  
  # 8. Return the list formatted for calculate_and_format_ate
  return(list(
    data_matched = df_matched,
    cem_object = cem_out,
    method = method_label,
    n_original = nrow(df_internal),
    n_matched = nrow(df_matched),
    n_dropped = nrow(df_internal) - nrow(df_matched)
  ))
}


# ============================================================================
# HELPER: Create CEM with automatic binning (for comparison)
# ============================================================================
cem_ps_auto <- function(data, 
                        config = col_config,
                        n_bins = NULL,  # Specify number of bins (K from DP)
                        method_label = "CEM (Auto Bins)") {
  
  df_internal <- as.data.frame(data)
  
  # If n_bins specified, create equal-width bins
  if (!is.null(n_bins)) {
    ps_range <- range(df_internal$ps, na.rm = TRUE)
    
    # Create n_bins equal-width intervals
    auto_endpoints <- seq(from = ps_range[1], 
                          to = ps_range[2], 
                          length.out = n_bins + 1)
    
    # Ensure endpoints cover [0, 1] if PS is in that range
    if (ps_range[1] >= 0 && ps_range[2] <= 1) {
      auto_endpoints[1] <- 0
      auto_endpoints[length(auto_endpoints)] <- 1
    }
    
    method_label <- sprintf("CEM (%d Equal Bins)", n_bins)
    
    cat(sprintf("Creating %d equal-width bins from PS range [%.3f, %.3f]\n",
                n_bins, ps_range[1], ps_range[2]))
    
    return(cem_ps(
      data = data,
      config = config,
      custom_endpoints = auto_endpoints,
      method_label = method_label
    ))
  } else {
    # No bins specified - use CEM's automatic binning (Scott's rule)
    return(cem_ps(
      data = data,
      config = config,
      custom_endpoints = NULL,
      method_label = method_label
    ))
  }
}

run_cem <- function(data, col_config, treatment_col, formula_cem) {
  cat("\n--- Running CEM (k2k) ---\n")
  
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
    cat("⚠️ CEM failed: No units were matched. Skipping CEM for this run.\n")
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
