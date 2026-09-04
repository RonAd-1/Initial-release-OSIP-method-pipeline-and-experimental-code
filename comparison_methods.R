# -- comparison_methods.R --


source("dataset_configs.R", echo = FALSE)

# # ==========================================
# # 3. THE MAIN WRAPPER: compare_methods
# # ==========================================
# compare_methods <- function(data, 
#                             data_config, 
#                             treatment_col, 
#                             levin_res = NULL, 
#                             cem_res = NULL,
#                             fm_res = NULL,
#                             # cem_res_auto = NULL, 
#                             optimal_1_1_res = NULL,
#                             card_res = NULL,
#                             gm_res = NULL,
#                             unmatched_data = NULL) {
#   
#   # Check for required columns
#   if (!"ps" %in% names(data)) {
#     stop("Input data must contain a 'ps' (propensity score) column.")
#   }
#   
#   # =========================================================================
#   # HELPER: Process one method and return metrics row
#   # =========================================================================
#   get_row <- function(df, label) {
#     
#     # Debug
#     # browser()
#     
#     # Now calculate ATE - it will always be negative and clinically correct!
#     row <- calculate_and_format_ate(list(data_matched = df), label, data_config, treatment_col)
#     
#     # Ensure sample sizes in the table are correct
#     row$n_treated <- sum(df[[treatment_col]] == 1)
#     row$n_control <- sum(df[[treatment_col]] == 0)
#     
#     return(row)
#   }
#   
#   # =========================================================================
#   # PROCESS ALL METHODS
#   # =========================================================================
#   res_list <- list()
#   
#   # --- UNMATCHED BASELINE ---
#   if (!is.null(unmatched_data)) {
#     
#     # Debug
#     # browser()
#     
#     cat("\nProcessing: Unmatched baseline...")
#     
#     # If it's a list (which you are passing), extract the dataframe.
#     # If it's already a dataframe, use it as is.
#     df_to_process <- if(is.list(unmatched_data)) unmatched_data$data_matched else unmatched_data
#     
#     # Ensure PS and weights are present
#     if (!"ps" %in% names(df_to_process)) {
#       df_to_process$ps <- data$ps[match(rownames(df_to_process), rownames(data))]
#     }
#     if (!"weights" %in% names(df_to_process)) {
#       df_to_process$weights <- 1
#     }
#     r_unmatched <- get_row(df_to_process, METHOD_LABELS[Methods$unmatched])
#     res_list <- append(res_list, list(r_unmatched))
#   }
#   
#   # Debug
#   
#   # browser()
#   
#   # --- Optimal MATCHING (1:1) ---
#   if (!is.null(optimal_1_1_res)) {
#     cat("\nProcessing: Optimal (1:1)...")
#     
#     # Restore PS column if missing
#     if (!"ps" %in% names(optimal_1_1_res$data_matched)) {
#       optimal_1_1_res$data_matched$ps <- data$ps[match(rownames(optimal_1_1_res$data_matched), 
#                                                        rownames(data))]
#     }
#     
#     r_card <- get_row(optimal_1_1_res$data_matched, METHOD_LABELS[Methods$optimal_1_1])
#     res_list <- append(res_list, list(r_card))
#   }
#   
#   # --- CARDINALITY MATCHING (1:1) ---
#   if (!is.null(card_res)) {
#     cat("\nProcessing: Cardinality (1:1)...")
#     
#     # Restore PS column if missing
#     if (!"ps" %in% names(card_res$data_matched)) {
#       card_res$data_matched$ps <- data$ps[match(rownames(card_res$data_matched), 
#                                                 rownames(data))]
#     }
#     
#     r_card <- get_row(card_res$data_matched, METHOD_LABELS[Methods$cardinality_1_1])
#     res_list <- append(res_list, list(r_card))
#   }
#   
#   #Debug
#   #browser()
#   
#   # --- LEVIN OPTIMIZED ---
#   if (!is.null(levin_res)) {
#     cat("\nProcessing: Levin Optimized...")
#     r_levin <- get_row(levin_res$data_matched, METHOD_LABELS[Methods$levin])
#     res_list <- append(res_list, list(r_levin))
#   }
#   
#   # --- CEM ---
#   if (!is.null(cem_res)) {
#     cat("\nProcessing: CEM...")
#     
#     # CRITICAL CHECK: Verify CEM actually matched
#     n_original <- nrow(data)
#     n_cem <- nrow(cem_res$data_matched)
#     
#     if (n_cem == n_original) {
#       warning("⚠️  CEM kept ALL units - matching may have failed!")
#       cat(sprintf("    Original: %d units, CEM: %d units (no change)\n",
#                   n_original, n_cem))
#     }
#     
#     r_cem <- get_row(cem_res$data_matched, METHOD_LABELS[Methods$cem])
#     res_list <- append(res_list, list(r_cem))
#   }
#   
#   # --- Full Match ---
#   if (!is.null(fm_res)) {
#     cat("\nProcessing: Full Match...")
#     
#     # Data extraction
#     df_fm <- fm_res$data_matched
#     n_original <- nrow(cov_df_standardized)
#     n_fm <- nrow(df_fm)
#     
#     # Use your get_row helper (ensure it handles weights if that's part of your metric)
#     r_fm <- get_row(fm_res$data_matched, METHOD_LABELS[Methods$full_matching])
#     res_list <- append(res_list, list(r_fm))
#   }
#   
#   # --- GENETIC ---
#   if (!is.null(gm_res)) {
#     cat("\nProcessing: GENETIC...")
#     r_gm <- get_row(gm_res$data_matched, METHOD_LABELS[Methods$genetic_1_1])
#     res_list <- append(res_list, list(r_gm))
#   }
#   
#   # =========================================================================
#   # COMBINE RESULTS
#   # =========================================================================
#   
#   #Debug
#   # browser()
#   
#   comparison_df <- do.call(rbind, res_list)
#   
#   # =========================================================================
#   # PRINT SUMMARY
#   # =========================================================================
#   cat("\n\n═══════════════════════════════════════════════════════════════\n")
#   cat("                 FINAL COMPARISON TABLE\n")
#   cat("═══════════════════════════════════════════════════════════════\n\n")
#   
#   print(comparison_df %>% 
#           dplyr::select(method, n_treated, n_control, ate, p_value) %>%
#           mutate(across(where(is.numeric), ~round(., 4))))
#   
#   return(comparison_df)
# }

# complete_comparison <- function(results_df, data, treatment_col, covariates, covariates_for_balance,
#                                 outcome_var, unmatched_with_weights, cem_matched, optimal_1_1_matched, 
#                                 cardinality_1_1_matched, genetic_1_1_matched, fm_matched, 
#                                 levin_matched, global_worst_levin, S_inv) {
#   # cem_matched, cem_auto_matched) {
#   
#   cat("\n╔═══════════════════════════════════════════════════════════╗\n")
#   cat("║        COMPREHENSIVE METHOD COMPARISON                   ║\n")
#   cat("╚═══════════════════════════════════════════════════════════╝\n")
#   
#   # ==========================================
#   # 1. COVARIATE BALANCE
#   # ==========================================
#   
#   # Debug
#   # browser()
#   
#   # Instead of balance_formula, pass columns directly
#   bal_unmatched <- bal.tab(
#     x = unmatched_with_weights[, covariates_for_balance, drop = FALSE],   # The data frame of covariates_for_balance
#     treat = unmatched_with_weights[[treatment_col]], # The treatment vector
#     stats = c("m", "ks"), # <--- Add this to get KS statistics
#     un = FALSE, 
#     s.d.denom = "pooled"
#   )
#   
#   bal_optimal_1_1_matched <- bal.tab(
#     x = optimal_1_1_matched[, covariates_for_balance, drop = FALSE],   # The data frame of covariates_for_balance
#     treat = optimal_1_1_matched[[treatment_col]], # The treatment vector
#     stats = c("m", "ks"), # <--- Add this to get KS statistics
#     un = FALSE, 
#     s.d.denom = "pooled"
#   )
#   
#   bal_cardinality_1_1_matched <- bal.tab(
#     x = cardinality_1_1_matched[, covariates_for_balance, drop = FALSE],   # The data frame of covariates_for_balance
#     treat = cardinality_1_1_matched[[treatment_col]], # The treatment vector
#     stats = c("m", "ks"), # <--- Add this to get KS statistics
#     un = FALSE, 
#     s.d.denom = "pooled"
#   )
#   
#   bal_levin <- bal.tab(
#     x = levin_matched[, covariates_for_balance, drop = FALSE],   # The data frame of covariates_for_balance
#     treat = levin_matched[[treatment_col]], # The treatment vector
#     stats = c("m", "ks"), # <--- Add this to get KS statistics
#     un = FALSE, 
#     s.d.denom = "pooled"
#   )
#   
#   bal_cem <- bal.tab(
#     x = cem_matched[, covariates_for_balance, drop = FALSE],   # The data frame of covariates_for_balance
#     treat = cem_matched[[treatment_col]], # The treatment vector
#     stats = c("m", "ks"), # <--- Add this to get KS statistics
#     un = FALSE,
#     s.d.denom = "pooled"
#   )
#   
#   bal_genetic_1_1_matched <- bal.tab(
#     x = genetic_1_1_matched[, covariates_for_balance, drop = FALSE],   # The data frame of covariates_for_balance
#     treat = genetic_1_1_matched[[treatment_col]], # The treatment vector
#     stats = c("m", "ks"), # <--- Add this to get KS statistics
#     un = FALSE,
#     s.d.denom = "pooled"
#   )
#   
#   bal_fm_matched <- bal.tab(
#     x = fm_matched[, covariates_for_balance, drop = FALSE],   # The data frame of covariates_for_balance
#     treat = fm_matched[[treatment_col]], # The treatment vector
#     stats = c("m", "ks"), # <--- Add this to get KS statistics
#     un = FALSE,
#     s.d.denom = "pooled"
#   )
#   
#   # bal_cem_auto <- bal.tab(
#   #   x = cem_auto_matched[, covariates, drop = FALSE],   # The data frame of covariates_for_balance
#   #   treat = cem_auto_matched[[treatment_col]], # The treatment vector
#   #   stats = c("m", "ks"), # <--- Add this to get KS statistics
#   #   un = FALSE, 
#   #   s.d.denom = "pooled"
#   # )
#   
#   # Extract metrics
#   get_balance_metrics <- function(bal_obj) {
#     stats <- bal_obj$Balance
#     
#     # This takes Diff.Adj if available/not-NA, otherwise falls back to Diff.Un
#     smd_val <- coalesce(stats[["Diff.Adj"]], stats[["Diff.Un"]])
#     ks_val  <- coalesce(stats[["KS.Adj"]], stats[["KS.Un"]])
#     
#     # Remove NAs to avoid math errors
#     smd_clean <- abs(smd_val[!is.na(smd_val)])
#     ks_clean  <- ks_val[!is.na(ks_val)]
#     
#     return(data.frame(
#       Max_SMD = if(length(smd_clean) > 0) max(smd_clean) else NA,
#       Mean_SMD = if(length(smd_clean) > 0) mean(smd_clean) else NA,
#       Max_KS = if(length(ks_clean) > 0) max(ks_clean) else NA,
#       N_Imbalanced_01 = sum(smd_clean > 0.1),
#       N_Imbalanced_02 = sum(smd_clean > 0.2)
#     ))
#   }
#   
#   balance_metrics <- rbind(
#     data.frame(Method = METHOD_LABELS[Methods$unmatched],       get_balance_metrics(bal_unmatched)),
#     data.frame(Method = METHOD_LABELS[Methods$optimal_1_1],     get_balance_metrics(bal_optimal_1_1_matched)),
#     data.frame(Method = METHOD_LABELS[Methods$cem],     get_balance_metrics(bal_cem)),
#     data.frame(Method = METHOD_LABELS[Methods$genetic_1_1],     get_balance_metrics(bal_genetic_1_1_matched)),
#     data.frame(Method = METHOD_LABELS[Methods$full_matching],     get_balance_metrics(bal_fm_matched)),
#     data.frame(Method = METHOD_LABELS[Methods$cardinality_1_1], get_balance_metrics(bal_cardinality_1_1_matched)),
#     data.frame(Method = METHOD_LABELS[Methods$levin],           get_balance_metrics(bal_levin))
#   )
#   
#   cat("\n═══ COVARIATE BALANCE ═══\n")
#   print(balance_metrics)
#   
#   # ==========================================
#   # 2. MATCHING DISTANCE QUALITY
#   # ==========================================
#   
#   #  calc_pair_distances <- function(matched_data, treatment_col, covariates, s_inv) {
#     
#     # 1. Initialize empty vectors at the top
#     distances_mahal <- numeric(0)
#     distances_ps <- numeric(0)
#     
#     # 2. Safety: If no match_id, return NAs immediately
#     if (!"match_id" %in% colnames(matched_data)) {
#       return(data.frame(Mean_Mahal=NA, Median_Mahal=NA, Max_Mahal=NA, SD_Mahal=NA, Mean_PS=NA, Max_PS=NA))
#     }
#     
#     match_ids <- unique(matched_data$match_id)
#     match_ids <- match_ids[!is.na(match_ids)]
#     
#     # 3. Covariance matrix for Mahalanobis (calculated once)
#     # Use only the matched data covariates
#     # cov_mat <- cov(matched_data[, covariates, drop = FALSE])
#     # Safety: If matrix is singular, we can't solve it
#     # inv_cov <- tryCatch(solve(cov_mat), error = function(e) return(NULL))
#     
#     for (mid in match_ids) {
#       pair <- matched_data[matched_data$match_id == mid, ]
#       
#       if (nrow(pair) == 2 && !is.null(s_inv)) {
#         t_idx <- which(pair[[treatment_col]] == 1)
#         c_idx <- which(pair[[treatment_col]] == 0)
#         
#         
#         # TODO: Make the logic similar to LEvin matching, use the matrix values.
#         
#         # Ensure we have one of each
#         if(length(t_idx) == 1 && length(c_idx) == 1) {
#           t_covs <- as.numeric(pair[t_idx, covariates])
#           c_covs <- as.numeric(pair[c_idx, covariates])
#           
#           diff <- t_covs - c_covs
#           dist_m <- sqrt(t(diff) %*% s_inv %*% diff)
#           distances_mahal <- c(distances_mahal, as.numeric(dist_m))
#           
#           if ("ps" %in% colnames(pair)) {
#             distances_ps <- c(distances_ps, abs(pair$ps[t_idx] - pair$ps[c_idx]))
#           }
#         }
#       }
#     }
#     
#     # 4. Return results (with check to avoid Mean of empty vector error)
#     has_dist <- length(distances_mahal) > 0
#     
#     return(data.frame(
#       Mean_Mahal   = if(has_dist) mean(distances_mahal, na.rm = TRUE) else NA,
#       Median_Mahal = if(has_dist) median(distances_mahal, na.rm = TRUE) else NA,
#       Max_Mahal    = if(has_dist) max(distances_mahal, na.rm = TRUE) else NA,
#       SD_Mahal     = if(has_dist && length(distances_mahal) > 1) sd(distances_mahal, na.rm = TRUE) else NA,
#       Mean_PS      = if(length(distances_ps) > 0) mean(distances_ps, na.rm = TRUE) else NA,
#       Max_PS       = if(length(distances_ps) > 0) max(distances_ps, na.rm = TRUE) else NA
#     ))
#   }
#   
#   # Debug
#   # browser()
# 
# calc_pair_distances_from_map <- function(map_data, data_matched) {
#   
#   # 1. Safety check: Ensure we have data
#   if (is.null(map_data) || nrow(map_data) == 0) {
#     return(data.frame(Mean_Mahal=NA, Median_Mahal=NA, Max_Mahal=NA, 
#                       SD_Mahal=NA, Mean_PS=NA, Max_PS=NA))
#   }
#   
#   # 2. We already have the Mahalanobis costs in the 'cost' column!
#   distances_mahal <- map_data$cost
#   
#   # 3. Handle Propensity Score distances
#   # We need to join with data_matched to get the 'ps' values for both IDs
#   # if they aren't already in the map
#   distances_ps <- numeric(0)
#   if ("ps_t" %in% colnames(map_data) && "ps_c" %in% colnames(map_data)) {
#     distances_ps <- abs(map_data$ps_t - map_data$ps_c)
#   } else {
#     # Fallback: Join to get PS if not in map
#     ps_lookup <- data_matched %>% 
#       mutate(id_str = as.character(id)) %>% 
#       dplyr::select(id_str, ps)
#     
#     map_with_ps <- map_data %>%
#       left_join(ps_lookup, by = c("treated_id" = "id_str")) %>%
#       rename(ps_t = ps) %>%
#       left_join(ps_lookup, by = c("control_id" = "id_str")) %>%
#       rename(ps_c = ps)
#     
#     distances_ps <- abs(map_with_ps$ps_t - map_with_ps$ps_c)
#   }
#   
#   # 4. Return results using the map's validated costs
#   has_dist <- length(distances_mahal) > 0
#   
#   return(data.frame(
#     Mean_Mahal   = if(has_dist) mean(distances_mahal, na.rm = TRUE) else NA,
#     Median_Mahal = if(has_dist) median(distances_mahal, na.rm = TRUE) else NA,
#     Max_Mahal    = if(has_dist) max(distances_mahal, na.rm = TRUE) else NA,
#     SD_Mahal     = if(has_dist && length(distances_mahal) > 1) sd(distances_mahal, na.rm = TRUE) else NA,
#     Mean_PS      = if(length(distances_ps) > 0) mean(distances_ps, na.rm = TRUE) else NA,
#     Max_PS       = if(length(distances_ps) > 0) max(distances_ps, na.rm = TRUE) else NA
#   ))
# }
#   
#   distance_metrics <- rbind(
#     # Unmatched
#     data.frame(Method = METHOD_LABELS[Methods$unmatched], Mean_Mahal = NA, Median_Mahal = NA,
#                Max_Mahal = NA, SD_Mahal = NA, Mean_PS = NA, Max_PS = NA),
#     
#     
#     # TODO: Do I really need CEM for distances? What does it mean? 
#     # CEM (Numeric) - Note: If you run CEM, replace na_dist_row with calc_pair_distances
#     data.frame(Method = METHOD_LABELS[Methods$cem], calc_pair_distances(cem_matched ,treatment_col, covariates, S_inv)),
#     
#     
#     # TODO: Same question for full matching? 
#     # CEM (Numeric) - Note: If you run CEM, replace na_dist_row with calc_pair_distances
#     data.frame(Method = METHOD_LABELS[Methods$full_matching], calc_pair_distances(fm_matched ,treatment_col, covariates, S_inv)),
#     
#     # Optimal (1:1)
#     data.frame(Method = METHOD_LABELS[Methods$optimal_1_1], 
#                calc_pair_distances(optimal_1_1_matched, treatment_col, covariates, S_inv)),
#     
#     # Cardinality (1:1)
#     data.frame(Method = METHOD_LABELS[Methods$cardinality_1_1], 
#                calc_pair_distances(cardinality_1_1_matched, treatment_col, covariates, S_inv)),
#     
#     # Genetic (1:1)
#     data.frame(Method = METHOD_LABELS[Methods$genetic_1_1], 
#                calc_pair_distances(genetic_1_1_matched, treatment_col, covariates, S_inv)),
#     
#     # Levin Optimized
#     data.frame(Method = METHOD_LABELS[Methods$levin], 
#                calc_pair_distances(levin_matched, treatment_col, covariates, S_inv))
#   )
#   
#   cat("\n═══ MATCHING DISTANCE QUALITY ═══\n")
#   print(distance_metrics)
#   
#   # ==========================================
#   # 3. SAMPLE SIZE & OVERLAP
#   # ==========================================
#   
#   get_sample_metrics <- function(matched_data, treatment_col) {
#     n_t <- sum(matched_data[[treatment_col]] == 1)
#     n_c <- sum(matched_data[[treatment_col]] == 0)
#     
#     # if (!is.null(dataset_name) && dataset_name == "LINDNER") {
#     #   # Clinical Control was flipped to 1, Clinical Treated was flipped to 0
#     #   n_c <- sum(matched_data[[treatment_col]] == 1)
#     #   n_t <- sum(matched_data[[treatment_col]] == 0)
#     # } else {
#     #   # Standard behavior for other datasets
#     #   n_t <- sum(matched_data[[treatment_col]] == 1)
#     #   n_c <- sum(matched_data[[treatment_col]] == 0)
#     # }
#     
#     if ("ps" %in% colnames(matched_data)) {
#       ps_t <- matched_data[matched_data[[treatment_col]] == 1, "ps"]
#       ps_c <- matched_data[matched_data[[treatment_col]] == 0, "ps"]
#       
#       ps_overlap <- min(max(ps_t), max(ps_c)) - max(min(ps_t), min(ps_c))
#       ps_range_ratio <- (max(ps_t) - min(ps_t)) / (max(ps_c) - min(ps_c))
#     } else {
#       ps_overlap <- NA
#       ps_range_ratio <- NA
#     }
#     
#     return(data.frame(
#       N_Control = n_c,
#       N_Treatment = n_c,
#       PS_Overlap = ps_overlap,
#       PS_Range_Ratio = ps_range_ratio
#     ))
#   }
#   
#   sample_metrics <- rbind(
#     data.frame(Method = METHOD_LABELS[Methods$unmatched], 
#                get_sample_metrics(unmatched_with_weights, treatment_col)),
#     
#     data.frame(Method = METHOD_LABELS[Methods$cem], 
#                get_sample_metrics(cem_matched, treatment_col)),
#     
#     data.frame(Method = METHOD_LABELS[Methods$full_matching], 
#                get_sample_metrics(fm_matched, treatment_col)),
#     
#     data.frame(Method = METHOD_LABELS[Methods$optimal_1_1], 
#                get_sample_metrics(optimal_1_1_matched, treatment_col)),
#     
#     data.frame(Method = METHOD_LABELS[Methods$cardinality_1_1], 
#                get_sample_metrics(cardinality_1_1_matched, treatment_col)),
#     
#     # Genetic (1:1)
#     data.frame(Method = METHOD_LABELS[Methods$genetic_1_1], 
#                get_sample_metrics(genetic_1_1_matched, treatment_col)),
#     
#     data.frame(Method = METHOD_LABELS[Methods$levin], 
#                get_sample_metrics(levin_matched, treatment_col))
#   )
#   
#   cat("\n═══ SAMPLE SIZE & OVERLAP ═══\n")
#   print(sample_metrics)
#   
#   # ==========================================
#   # 4. COMBINED RESULTS TABLE
#   # ==========================================
#   
#   # Define the order using the Map values
#   method_order <- c(
#     METHOD_LABELS[Methods$unmatched],
#     METHOD_LABELS[Methods$cem],
#     METHOD_LABELS[Methods$full_matching],
#     METHOD_LABELS[Methods$optimal_1_1],
#     METHOD_LABELS[Methods$cardinality_1_1],
#     METHOD_LABELS[Methods$genetic_1_1],
#     METHOD_LABELS[Methods$levin]
#   )
#   
#   results_df <- results_df[
#     match(method_order, results_df$method),
#   ]
#   
#   rel_eff <- 1 / (results_df$se^2)
#   eff_rank <- rank(-rel_eff, ties.method = "min")
#   best_eff <- max(rel_eff, na.rm = TRUE)
#   eff_score <- (rel_eff / best_eff) * 100
#   
#   combined <- data.frame(
#     Method = method_order,
#     N_Treated = results_df$n_treated,
#     N_Control = sample_metrics$N_Control,
#     
#     # Balance
#     Max_SMD = balance_metrics$Max_SMD,
#     Mean_SMD = balance_metrics$Mean_SMD,
#     Max_KS = balance_metrics$Max_KS,
#     N_Imbal_01 = balance_metrics$N_Imbalanced_01,
#     
#     # Distance
#     Mean_Mahal = distance_metrics$Mean_Mahal,
#     Max_Mahal = distance_metrics$Max_Mahal,
#     Mean_PS_Dist = distance_metrics$Mean_PS,
#     
#     # ATE 
#     ATE = results_df$ate,
#     SE = results_df$se,
#     T_Stat = results_df$ate / results_df$se,
#     P_Value = results_df$p_value,
#     CI_Lower = results_df$ci_low,
#     CI_Upper = results_df$ci_high,
#     
#     # Efficiency Columns
#     Relative_Efficiency = rel_eff,
#     Eff_Rank = eff_rank,
#     Eff_Score_Pct = eff_score
#   )
#   
#   cat("\n═══ COMBINED RESULTS TABLE ═══\n")
#   print(combined)
#   
#   return(list(
#     balance = balance_metrics,
#     distances = distance_metrics,
#     samples = sample_metrics,
#     combined = combined,
#     balance_objects = list(
#       unmatched = bal_unmatched,
#       cem = bal_cem,
#       full_matching = bal_fm_matched,
#       optimal_1_1 = bal_optimal_1_1_matched,
#       cardinality_1_1 = bal_cardinality_1_1_matched,
#       genetic_1_1 = bal_genetic_1_1_matched,
#       levin = bal_levin
#       # cem_auto = bal_cem_auto
#     )
#   ))
# }

# run_delta_sensitivity_analysis <- function(data_subset, params, data_config, p_sorted, delta_values) {
#   
#   all_results <- list()
#   
#   # Loop through Delta values (e.g., 0.1, 0.15, 0.2...)
#   for (current_delta in delta_values) {
#     cat(sprintf("\n\n>>> STARTING ANALYSIS FOR DELTA = %.2f <<<\n", current_delta))
#     
#     # 1. Dynamically calculate K upper bound
#     # Formula: K_max = floor(1 / Delta) + 2
#     k_max <- floor(1 / current_delta) + 2
#     cat(sprintf("Setting K upper bound to: %d\n", k_max))
#     
#     # 2. Update parameters for this specific run
#     current_params <- params
#     current_params$DELTA <- current_delta
#     current_params$K_MAX <- k_max
#     
#     # 3. Launch Step 1: DP Baseline
#     # This finds the optimal partition for the current Delta constraint
#     step1_res <- run_levin_step1(data_subset, params, data_config, p_sorted)
#     
#     # 4. Launch Step 2: Heuristic Refinement
#     # This optimizes the original covariate distance (Mahalanobis)
#     step2_res <- run_levin_step2_heuristics(step1_res, data, current_params, data_config)
#     
#     # 5. Extract the "Best" result from Step 2
#     # Using the attr-based logic discussed previously
#     best_method <- attr(step2_res, "best_method")
#     final_match <- step2_res[[best_method]]
#     
#     # Store results indexed by Delta
#     all_results[[as.character(current_delta)]] <- list(
#       delta = current_delta,
#       k_used = final_match$k,
#       cost = final_match$cost,
#       intervals = final_match$intervals
#     )
#   }
#   
#   return(all_results)
# }

calculate_rosenbaum_gamma <- function(match_map, data_subset, outcome_var, gamma_range = seq(1, 5, by = 0.5)) {
  
  # Debug
  # browser()
  
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
  
  # Debug
  # browser()
  
  # --- THE FIX ---
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

# find_gamma_root <- function(y_matrix, alpha = 0.05,
#                             lower = 1, upper = 2) {
#   
#   f <- function(g) {
#     senmv(y_matrix, gamma = g, method = "h")$pval - alpha
#   }
#   
#   uniroot(f, lower = lower, upper = upper)$root
# }

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

# get_match_map <- function(m_input, dist_mat) {
# 
#   # CASE 1: Standard MatchIt object
#   if (inherits(m_input, "matchit")) {
#     matches <- m_input$match.matrix
#     map <- data.frame(
#       treated_id = rownames(matches),
#       control_id = matches[, 1],
#       stringsAsFactors = FALSE
#     ) %>% filter(!is.na(control_id))
# 
#   } else if (inherits(m_input, "optmatch")) {
#     # CASE 2: The 'optmatch' factor vector from your image
#     # We find IDs that share the same factor level
#     map <- data.frame(
#       id = names(m_input),
#       group = as.character(m_input),
#       stringsAsFactors = FALSE
#     ) %>%
#       filter(!is.na(group)) %>%
#       # We need to distinguish who is Treated and who is Control
#       # We look up their treatment status in the dist_mat rownames
#       mutate(is_treated = id %in% rownames(dist_mat)) %>%
#       group_by(group) %>%
#       summarize(
#         treated_id = id[is_treated == TRUE][1],
#         control_id = id[is_treated == FALSE][1],
#         .groups = 'drop'
#       ) %>%
#       filter(!is.na(treated_id) & !is.na(control_id))
# 
#   } else {
#     # CASE 3: Already a dataframe (Levin)
#     map <- as.data.frame(m_input)
#   }
# 
#   # Final Step: Attach Costs from the Matrix
#   map$cost <- mapply(function(t, c) {
#     if (t %in% rownames(dist_mat) && c %in% colnames(dist_mat)) {
#       return(dist_mat[as.character(t), as.character(c)])
#     } else { return(NA) }
#   }, map$treated_id, map$control_id)
# 
#   return(as.data.frame(map))
# }

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