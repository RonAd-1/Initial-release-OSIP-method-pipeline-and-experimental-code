# -- full_matching.R --

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