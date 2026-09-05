# -- comparison_methods.R --


# source("dataset_configs.R", echo = FALSE)
# 
# calculate_rosenbaum_gamma <- function(match_map, data_subset, outcome_var, gamma_range = seq(1, 5, by = 0.5)) {
#   
#   # 1. Prepare the outcome vector (ensuring alignment with the map IDs)
#   outcomes <- data_subset[[outcome_var]]
#   names(outcomes) <- as.character(data_subset$id)
#   
#   # 2. Create the Pair-Difference Matrix (Required for Rosenbaum)
#   # Rows = Matched Sets, Col 1 = Treated Outcome, Col 2 = Control Outcome
#   y_matrix <- matrix(NA, nrow = nrow(match_map), ncol = 2)
#   
#   y_matrix[, 1] <- outcomes[as.character(match_map$treated_id)]
#   y_matrix[, 2] <- outcomes[as.character(match_map$control_id)]
#   
#   # Remove any rows with NAs (unmatched units)
#   y_matrix <- y_matrix[complete.cases(y_matrix), ]
#   
#   # Check the average difference. If it's negative, flip the matrix 
#   # so we are testing the magnitude of the effect correctly.
#   avg_diff <- mean(y_matrix[, 1] - y_matrix[, 2])
#   if (avg_diff < 0) {
#     y_matrix <- -y_matrix
#   }
#   # ----------------
#   
#   # 3. Iterate through Gamma values
#   sens_results <- data.frame(Gamma = gamma_range, P_Value_Bound = NA)
#   
#   for(i in 1:nrow(sens_results)) {
#     g <- sens_results$Gamma[i]
#     # senmv calculates the upper bound p-value for a given Gamma
#     res <- senmv(y_matrix, gamma = g, method = "h")
#     # Explicitly grab only the p-value to avoid the length warning
#     sens_results$P_Value_Bound[i] <- res$pval
#   }
#   
#   # Recommended Logic Flow
#   # --- Logic Flow ---
#   sig_indices <- which(sens_results$P_Value_Bound <= 0.05)
#   
#   if (length(sig_indices) == 0) {
#     # Scenario A: Not significant at Gamma = 1.0
#     threshold_gamma <- 1.0
#     exact_threshold <- 1.0  # Initialize here so it exists for the return list
#   } else {
#     # Scenario B: Significant! Now find the exact interpolated point
#     exact_val <- estimate_gamma_threshold(sens_results)
#     
#     # Fallback logic
#     exact_threshold <- if (!is.na(exact_val)) exact_val else max(sens_results$Gamma)
#     threshold_gamma <- exact_threshold # You can use the same value or the discrete one
#   }
#   
#   # Now both variables are guaranteed to exist
#   return(list(
#     results_table = sens_results,
#     threshold = threshold_gamma,
#     exact_threshold = exact_threshold
#   ))
# }
# 
# estimate_gamma_threshold <- function(results, alpha = 0.05) {
#   
#   # Find first p-value above alpha
#   idx <- which(results$P_Value_Bound > alpha)[1]
#   
#   if (is.na(idx) || idx == 1) {
#     return(NA)
#   }
#   
#   g1 <- results$Gamma[idx - 1]
#   g2 <- results$Gamma[idx]
#   p1 <- results$P_Value_Bound[idx - 1]
#   p2 <- results$P_Value_Bound[idx]
#   
#   # Linear interpolation
#   g_star <- g1 + (alpha - p1) * (g2 - g1) / (p2 - p1)
#   
#   return(g_star)
# }
# 
# get_match_map <- function(m_input, dist_mat) {
#   
#   # CASE 1: Standard MatchIt object
#   if (inherits(m_input, "matchit")) {
#     matches <- m_input$match.matrix
#     map <- data.frame(
#       treated_id = as.character(rownames(matches)),
#       control_id = as.character(matches[, 1]),
#       stringsAsFactors = FALSE
#     ) %>% dplyr::filter(!is.na(control_id))
#     
#   } else if (inherits(m_input, "optmatch")) {
#     # CASE 2: The 'optmatch' factor vector
#     map <- data.frame(
#       id = as.character(names(m_input)),
#       group = as.character(m_input),
#       stringsAsFactors = FALSE
#     ) %>%
#       dplyr::filter(!is.na(group)) %>%
#       mutate(is_treated = id %in% rownames(dist_mat)) %>%
#       group_by(group) %>%
#       summarize(
#         treated_id = id[is_treated == TRUE][1],
#         control_id = id[is_treated == FALSE][1],
#         .groups = 'drop'
#       ) %>%
#       dplyr::filter(!is.na(treated_id) & !is.na(control_id))
#     
#   } else {
#     # CASE 3: Already a dataframe (Levin)
#     map <- as.data.frame(m_input)
#     map$treated_id <- as.character(map$treated_id)
#     map$control_id <- as.character(map$control_id)
#   }
#   
#   # --- UPDATED SYMMETRIC LOOKUP STARTS HERE ---
#   # We use a character-safe, direction-agnostic check to handle 'flipped' IDs
#   map$cost <- mapply(function(t, c) {
#     t_chr <- as.character(t)
#     c_chr <- as.character(c)
#     
#     # 1. Try Standard Orientation: Treated in rows, Control in columns
#     if (t_chr %in% rownames(dist_mat) && c_chr %in% colnames(dist_mat)) {
#       return(dist_mat[t_chr, c_chr])
#     } 
#     
#     # 2. Try Flipped Orientation: Control in rows, Treated in columns
#     # This solves the Jobs dataset issue where IDs like "1" are columns in the matrix
#     # but appear in the treated_id column of the quintile map.
#     else if (c_chr %in% rownames(dist_mat) && t_chr %in% colnames(dist_mat)) {
#       return(dist_mat[c_chr, t_chr])
#     } 
#     
#     # 3. Fallback
#     else { 
#       return(NA) 
#     }
#   }, map$treated_id, map$control_id)
#   # --- UPDATED SYMMETRIC LOOKUP ENDS HERE ---
#   
#   return(as.data.frame(map))
# }
# 
# print_efficiency_rankings <- function(efficiency_vector) {
#   cat("\n--- Efficiency Rankings (Most Efficient to Least) ---\n")
#   
#   # 1. Create a data frame for easier manipulation
#   eff_df <- data.frame(
#     Method = names(efficiency_vector),
#     Efficiency = as.numeric(efficiency_vector)
#   )
#   
#   # 2. Rank: Higher is better (usually Efficiency = 1/Variance)
#   # If your Efficiency metric is "Higher = Better", sort descending
#   eff_df <- eff_df[order(-eff_df$Efficiency), ]
#   
#   # 3. Add a "Ratio to Best" column to make differences intuitive
#   # This shows how much efficiency you lose compared to the top performer
#   best_val <- eff_df$Efficiency[1]
#   eff_df$Ratio_to_Best <- eff_df$Efficiency / best_val
#   
#   # 4. Format for printing
#   eff_df$Efficiency_Fixed <- format(eff_df$Efficiency, scientific = FALSE, digits = 8)
#   eff_df$Rank <- 1:nrow(eff_df)
#   
#   # 5. Final Display
#   for(i in 1:nrow(eff_df)) {
#     cat(sprintf("%d. %-15s | Value: %s | Score: %.2f%%\n", 
#                 eff_df$Rank[i], 
#                 eff_df$Method[i], 
#                 eff_df$Efficiency_Fixed[i],
#                 eff_df$Ratio_to_Best[i] * 100))
#   }
# }