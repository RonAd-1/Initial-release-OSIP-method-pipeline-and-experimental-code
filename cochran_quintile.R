# cochran_quintile.R

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

# run_quintile_matching <- function(data, 
#                                   config = col_config, 
#                                   n_subclasses = 5, 
#                                   method_label = "Quintile (1:1)") {
#   
#   # 1. Force to base data.frame
#   df_internal <- as.data.frame(data)
#   treat_var <- config$TREATMENT_VAR
#   
#   # 2. Step A: Subclassification to define the bins
#   # This identifies which units fall into which quintile based on PS
#   m_sub <- matchit(as.formula(paste(treat_var, "~ ps")), 
#                    data = df_internal, 
#                    method = "subclass", 
#                    subclass = n_subclasses)
#   
#   # Assign subclass IDs (contains NAs for units outside common support)
#   df_internal$subclass_id <- m_sub$subclass
#   
#   # 3. FIX: Handle the 'Symbol' error by filtering NAs
#   # We only pass units with a valid subclass to the optimal matcher
#   # This prevents the VECTOR_ELT() error in the C++ backend
#   df_for_optimal <- df_internal[!is.na(df_internal$subclass_id), ]
#   
#   # Check if we still have treated and control units after filtering
#   if (length(unique(df_for_optimal[[treat_var]])) < 2) {
#     warning("Insufficient common support: One group has no units in any subclass.")
#     return(NULL)
#   }
#   
#   # 4. Build the formula for Mahalanobis distance
#   match_formula <- reformulate(termlabels = config$ALL_COVARIATES, 
#                                response = treat_var)
#   
#   cat(sprintf("Running %s (1:1 Optimal within %d strata)...\n", 
#               method_label, n_subclasses))
#   
#   # 5. Step B: Perform 1:1 Optimal Matching within those subclasses
#   # exact = ~ subclass_id ensures matches only happen within the same quintile
#   m_quint <- matchit(match_formula, 
#                      data = df_for_optimal,
#                      method = "optimal",
#                      distance = "mahalanobis",
#                      exact = ~ subclass_id)
#   
#   # 6. Extract Results and map weights back to the ORIGINAL data frame
#   # This ensures we don't "lose" the dropped units in the final return object
#   df_internal$weights <- 0 # Initialize all to 0
#   
#   # Match weights back to df_internal by row names (or your id_var)
#   # m_quint$weights only contains values for df_for_optimal
#   df_internal[names(m_quint$weights), "weights"] <- m_quint$weights
#   
#   df_matched <- df_internal[df_internal$weights > 0, ]
#   
#   n_total <- nrow(df_internal)
#   n_matched <- nrow(df_matched)
#   
#   cat(sprintf("  Matched: %d units (T=%d, C=%d)\n",
#               n_matched,
#               sum(df_matched[[treat_var]] == 1),
#               sum(df_matched[[treat_var]] == 0)))
#   
#   return(list(
#     data_matched = df_matched,
#     matchit_object = m_quint,
#     method = method_label,
#     n_original = n_total,
#     n_matched = n_matched,
#     n_dropped = n_total - n_matched
#   ))
# }