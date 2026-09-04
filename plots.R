# --- plots.R ---

visualize_osip_step_results <- function(data_subset, cost, best_intervals, 
                                         title_suffix, method_used, 
                                         data_config, delta_dp, dataset_name, 
                                         base_dir, treatment_col) 
{
 
  bins_plot_title <- sprintf("osip Partition for delta = %.2f: %s %s (est_cost = %.2f)", delta_dp, method_used, title_suffix, cost)
  
  plot_partition_with_bins(
    data_subset = data_subset,
    intervals_df = best_intervals, 
    title = bins_plot_title,
    data_config = data_config,
    delta_dp = delta_dp,
    dataset_name = dataset_name,
    title_suffix = title_suffix, 
    base_dir = base_dir,
    treatment_col = treatment_col
  )
}

# ==========================================
# CREATE COMPARISON VISUALIZATIONS
# ==========================================

# create_smd_comparison <- function(results, delta_dp, dataset_name, base_dir) {
#   # 1. SMD Reduction Bar Chart
#   smd_comparison <- results$balance
#   
#   # 1. Define the parameters for this run
#   delta_label <- sprintf("%.2f", delta_dp)
#   
#   smd_title <- sprintf("Maximum Absolute SMD by Method (Delta = %s)", delta_label)
#   
#   p_smd <- ggplot(smd_comparison, aes(x = Method, y = Max_SMD, fill = Method)) +
#     geom_bar(stat = "identity") +
#     geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
#     geom_text(aes(label = sprintf("%.3f", Max_SMD)), vjust = -0.5) +
#     labs(title = smd_title,
#          subtitle = "Lower is better (dashed line at 0.1)",
#          y = "Max |SMD|",
#          x = "") +
#     theme_minimal() +
#     theme(legend.position = "none")
#   
#   storing_path <- sprintf("%ssmd_plot_data_%s_delta_%.2f.png", base_dir, dataset_name, delta_dp)
#   ggsave(storing_path, p_smd, width = 10, height = 8, dpi = 300)
#   return(p_smd)
# }

# ==========================================
# CREATE LOVE PLOT (Visual Balance Comparison)
# ==========================================

# create_love_plot <- function(bal_list, delta_dp, dataset_name, base_dir) {
#   
#   # 1. Define the titles
#   delta_label <- sprintf("%.2f", delta_dp)
#   love_plot_title <- sprintf("Covariate Balance Comparison (Delta = %s)", delta_label)
#   
#   love_data <- data.frame()
#   
#   # Loop through the list of bal.tab objects
#   for (method_key in names(bal_list)) {
#     bal_obj <- bal_list[[method_key]]
#     
#     # Extract the balance table as a data frame
#     stats_df <- as.data.frame(bal_obj$Balance)
#     
#     # 2. Pick the correct SMD column (Adj if possible, else Un)
#     # We use explicit column checking instead of coalesce for safety
#     if ("Diff.Adj" %in% colnames(stats_df)) {
#       smd_vals <- stats_df$Diff.Adj
#     } else {
#       smd_vals <- stats_df$Diff.Un
#     }
#     
#     # 3. Build the data frame for this method
#     df <- data.frame(
#       Method = method_key,
#       Covariate = rownames(stats_df),
#       SMD = as.numeric(smd_vals),
#       stringsAsFactors = FALSE
#     )
#     
#     love_data <- rbind(love_data, df)
#   }
#   
#   # --- CLEANING ---
#   # Remove the 'distance' rows so they don't clutter the plot
#   love_data <- love_data[!love_data$Covariate %in% c("distance", "prop.score"), ]
#   
#   # Strip suffixes (e.g., .0, .1) from factor variables
#   love_data$Covariate <- gsub("[_\\.]\\d.*", "", love_data$Covariate)
#   
#   # Filter out rows with NA SMD (in case a method failed)
#   love_data <- love_data[!is.na(love_data$SMD), ]
#   
#   # Map keys to pretty labels (e.g., "levin" -> "Levin Optim")
#   love_data$Method <- factor(love_data$Method, 
#                              levels = names(METHOD_LABELS), 
#                              labels = METHOD_LABELS)
#   
#   manual_shapes <- c(16, 17, 15, 18, 19, 8, 4)
#   
#   # Plotting
#   p <- ggplot(love_data, aes(x = abs(SMD), y = reorder(Covariate, abs(SMD)), 
#                              color = Method, shape = Method)) +
#     geom_point(size = 3, alpha = 0.8, position = position_dodge(width = 0.5)) +
#     scale_shape_manual(values = c(16, 17, 15, 18, 19, 8, 4)) +
#     geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
#     geom_vline(xintercept = 0.2, linetype = "dotted", color = "darkred") +
#     labs(title = love_plot_title,
#          subtitle = "Dashed lines at |SMD| = 0.1 and 0.2",
#          x = "Absolute Standardized Mean Difference",
#          y = "Covariates") +
#     theme_minimal() +
#     theme(legend.position = "bottom")
#   
#   # Saving
#   storing_path <- file.path(base_dir, sprintf("love_plot_%s_delta_%.2f.png", dataset_name, delta_dp))
#   ggsave(storing_path, p, width = 10, height = 8, dpi = 300)
#   
#   return(p)
# }

# create_love_plot <- function(balance_results, delta_dp, dataset_name, base_dir) {
  
  # Extract SMD for each covariate from each method
#   bal_list <- balance_results$balance_objects
#   
#   # 1. Define the parameters for this run
#   delta_label <- sprintf("%.2f", delta_dp)
#   
#   love_plot_title <- sprintf("Covariate Balance Comparison (Delta = %s)", delta_label)
#   
#   love_data <- data.frame()
#   
#   for (method_name in names(bal_list)) {
#     bal_obj <- bal_list[[method_name]]
#     
#     # 1. Identify the available columns
#     # We use coalesce to pick Adj if it's there, otherwise Un
#     # If Diff.Adj is a column of NAs, coalesce will fill it with Diff.Un
#     smd_vector <- dplyr::coalesce(bal_obj$Balance$Diff.Adj, bal_obj$Balance$Diff.Un)
#     
#     # 2. Build the data frame
#     df <- data.frame(
#       Method = method_name,
#       Covariate = rownames(bal_obj$Balance),
#       SMD = smd_vector
#     )
#     
#     love_data <- rbind(love_data, df)
#   }
#   
#   # 1. Strip the dummy suffixes (e.g., everything after the underscore or dot)
#   # This regex looks for an underscore or dot followed by numbers and replaces it with nothing
#   love_data$Covariate <- gsub("[_\\.]\\d.*", "", love_data$Covariate)
#   
#   # 2. (Optional) Force specific clean names if simple stripping isn't enough
#   # # This ensures "marr" becomes "Married", etc.
#   # clean_cov_map <- c(
#   #   "marr"     = "Married",
#   #   "nodegree" = "No Degree",
#   #   "hisp"     = "Hispanic",
#   #   "black"    = "Black",
#   #   "re74"     = "Real Earnings '74",
#   #   "re75"     = "Real Earnings '75",
#   #   "stent"    = "Stent",
#   #   "acutemi"  = "Acute MI"
#   # )
#   
#   # Apply the mapping only to names that exist in your map
#   # love_data$Covariate <- ifelse(love_data$Covariate %in% names(clean_cov_map), 
#   #                               clean_cov_map[love_data$Covariate], 
#   #                               love_data$Covariate)
#   
#   #Debug
#   # browser()
#   
#   # 3. Now run your Method label cleaning
#   love_data$Method <- factor(love_data$Method, 
#                              levels = names(METHOD_LABELS), 
#                              labels = METHOD_LABELS)
#   
#   #debug
#   # browser()
#   
#   manual_shapes <- c(16, 17, 15, 18, 19, 8, 4)
#   
#   # Create plot
#   p <- ggplot(love_data, aes(x = abs(SMD), y = Covariate, 
#                              color = Method, shape = Method)) +
#     # Add position_dodge here to shift points vertically
#     geom_point(size = 3, alpha = 0.7, 
#                position = position_dodge(width = 0.5)) +
#     
#     # --- ADD THIS LINE HERE ---
#     scale_shape_manual(values = manual_shapes) +
#     # --------------------------
#   
#     geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", linewidth = 0.5) +
#     geom_vline(xintercept = 0.2, linetype = "dashed", color = "darkred", linewidth = 0.5) +
#     labs(title = love_plot_title,
#          subtitle = "Dashed lines at |SMD| = 0.1 and 0.2",
#          x = "Absolute Standardized Mean Difference",
#          y = "") +
#     theme_minimal() +
#     theme(legend.position = "bottom",
#           plot.title = element_text(face = "bold", size = 14))
#   
#   storing_path <- sprintf("%slove_plot_data_%s_delta_%.2f.png", base_dir, dataset_name, delta_dp)
#   ggsave(storing_path, p, width = 10, height = 8, dpi = 300)
#   return(p)
# }

plot_ps_distribution <- function(df, dataset_name, treatment_col, datasets, base_dir, title_suffix = "") {
  # Map dataset name to the correct column
  # Using the global datasets list logic
  # treatment_col <- data_config$TREATMENT_VAR
  
  # Debug
  # browser()
  
  df$treat_label <- ifelse(df[[treatment_col]] == 1, "Treated", "Control")
  df$treat_jitter <- jitter(as.numeric(df[[treatment_col]]), amount = 0.05)
  
  fig <- plot_ly(
    data = df,
    x = ~ps,
    y = ~treat_jitter,
    type = "scatter",
    mode = "markers",
    color = ~treat_label,
    colors = c("blue", "red"),
    hoverinfo = "text",
    text = ~paste(
      "Group:", treat_label,
      "<br>Propensity Score:", round(ps, 4)
    )
  ) %>%
    layout(
      title = paste("PS Distribution:", dataset_name, title_suffix),
      xaxis = list(title = "Propensity Score", range = c(0, 1.0)),
      yaxis = list(title = "", showticklabels = FALSE)
    )
  
  # 1. Update your path to end in .jpg
  storing_path <- sprintf("%sps_distribution_%s.jpg", base_dir, dataset_name)
  
  # 1. Create the plot using ggplot
  # We use 'text' as a dummy aesthetic for the hover info
  p <- ggplot(df, aes(x = ps, y = treat_jitter, color = treat_label, 
                      text = paste("Group:", treat_label, 
                                   "\nPropensity Score:", round(ps, 4)))) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("blue", "red")) +
    labs(
      title = paste("PS Distribution:", dataset_name, title_suffix),
      x = "Propensity Score",
      y = "",
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    xlim(0, 1.0)
  
  # 2. Save it directly and easily
  storing_path <- sprintf("%sps_distribution_%s.png", base_dir, dataset_name)
  ggsave(storing_path, plot = p, width = 10, height = 8, dpi = 300)
  
  # Display the plot directly
  print(fig)
}

# plot_partition_with_bins <- function(data_subset, intervals_df, redundant_unit_ids, 
#                                      title, title_suffix, data_config, delta_dp, dataset_name, base_dir, treatment_col,                                    
#                                      show_plot = TRUE, make_interactive = TRUE) {
#   
#   # treatment_col = data_config$TREATMENT_VAR
#   
#   # Debug 
#   # browser()
#   
#   stage_label = title_suffix
#   # 1. Define the parameters for this run
#   delta_label <- sprintf("%.2f", delta_dp)
#   
#   if (!is.data.frame(intervals_df)) {
#     stop("intervals_df must be a data frame from DP reconstruction.")
#   }
#   
#   if (nrow(intervals_df) == 0) {
#     message("intervals_df is empty. Returning NULL to skip partition plotting.")
#     return(NULL)
#   }
#   
#   # --- 1. Prepare Intervals ---
#   intervals_df$bin_index <- 1:nrow(intervals_df) 
#   intervals_df$ps_display_center <- 1:nrow(intervals_df)
#   intervals_df$ps_display_start <- intervals_df$ps_display_center - 0.5
#   intervals_df$ps_display_end <- intervals_df$ps_display_center + 0.5
#   
#   # --- 2. Prepare Data ---
#   data_subset$bin_index <- NA
#   n_intervals <- nrow(intervals_df)
#   
#   # DEBUG: Print input data summary
#   cat(sprintf("\n=== Plotting Debug ===\n"))
#   cat(sprintf("Total units: %d (T=%d, C=%d)\n", 
#               nrow(data_subset),
#               sum(data_subset[[treatment_col]] == 1),
#               sum(data_subset[[treatment_col]] == 0)))
#   cat(sprintf("Intervals: %d\n", nrow(intervals_df)))
#   
#   # Create breaks from intervals
#   all_boundaries <- c(intervals_df$start_ps, intervals_df$end_ps)
#   breaks <- sort(unique(all_boundaries))
#   
#   # --- 2. Prepare Data (Fixed for Epsilon Synchronization) ---
#   # data_subset$bin_index <- NA
#   # eps_plot <- 1e-9 # Consistent with your distance function
#   
#   # Debug
#   # browser()
#   
#   # n_intervals <- nrow(intervals_df)
#   # eps <- 1e-12
#   
#   # Identify neighbors for boundary logic
#   has_empty_before <- c(FALSE, intervals_df$is_empty[-n_intervals])
#   
#   # Debug
#   # browser()
#   
#   # In plot_partition_with_bins:
#   
#   # In plot_partition_with_bins:
#   
#   cat("\n=== Assigning bins by ID ===\n")
#   
#   for (i in 1:nrow(intervals_df)) {
#     # Skip empty intervals
#     if (intervals_df$is_empty[i]) {
#       cat(sprintf("Bin %d: EMPTY - skipping\n", i))
#       next
#     }
#     
#     # 🔑 SAFE extraction for list columns or regular columns
#     ids_value <- tryCatch({
#       # Try list column first
#       if (is.list(intervals_df$ids)) {
#         intervals_df$ids[[i]]  # Use [[ for list columns
#       } else {
#         intervals_df$ids[i]     # Use [ for regular columns
#       }
#     }, error = function(e) {
#       NA
#     })
#     
#     # Debug
#     cat(sprintf("Bin %d: ids_value type = %s, length = %d\n", 
#                 i, class(ids_value), length(ids_value)))
#     
#     # Check if we have valid IDs
#     if (length(ids_value) == 0 || (length(ids_value) == 1 && is.na(ids_value))) {
#       cat(sprintf("  ⚠️ Bin %d: No IDs - using boundary fallback\n", i))
#       
#       # Fallback to boundary logic
#       vprev <- intervals_df$start_ps[i]
#       vcurr <- intervals_df$end_ps[i]
#       eps <- 1e-12
#       
#       # Use same logic as calculate_distance_from_dp_intervals
#       has_empty_before <- i > 1 && intervals_df$is_empty[i-1]
#       
#       if (has_empty_before) {
#         bin_mask <- (data_subset$ps >= vprev - eps) & (data_subset$ps <= vcurr + eps)
#       } else {
#         bin_mask <- (data_subset$ps > vprev + eps) & (data_subset$ps <= vcurr + eps)
#       }
#       
#       data_subset$bin_index[bin_mask] <- intervals_df$bin_index[i]
#       cat(sprintf("  Assigned %d units via boundaries\n", sum(bin_mask, na.rm = TRUE)))
#       
#     } else {
#       # We have IDs - parse them
#       
#       # If it's a character string, split it
#       if (is.character(ids_value)) {
#         unit_ids <- as.numeric(strsplit(ids_value, ",")[[1]])
#       } else {
#         # Already numeric vector
#         unit_ids <- as.numeric(ids_value)
#       }
#       
#       # Remove NAs
#       unit_ids <- unit_ids[!is.na(unit_ids)]
#       
#       cat(sprintf("  Parsed %d IDs\n", length(unit_ids)))
#       
#       if (length(unit_ids) > 0) {
#         matches <- data_subset[[data_config$ID_VAR]] %in% unit_ids
#         data_subset$bin_index[matches] <- intervals_df$bin_index[i]
#         cat(sprintf("  ✓ Assigned %d units by ID\n", sum(matches)))
#       }
#     }
#   }
#   
#   # Validation
#   cat("\n=== Bin Assignment Validation ===\n")
#   for (i in 1:nrow(intervals_df)) {
#     assigned_t <- sum(data_subset$bin_index == i & 
#                         data_subset[[treatment_col]] == 1, na.rm = TRUE)
#     assigned_c <- sum(data_subset$bin_index == i & 
#                         data_subset[[treatment_col]] == 0, na.rm = TRUE)
#     expected_t <- intervals_df$n_treated[i]
#     expected_c <- intervals_df$n_control[i]
#     
#     match_status <- ifelse(assigned_t == expected_t && assigned_c == expected_c, 
#                            "✓", "❌")
#     
#     cat(sprintf("Bin %d: Expected T=%d C=%d, Got T=%d C=%d %s\n",
#                 i, expected_t, expected_c, assigned_t, assigned_c, match_status))
#   }
#   
#   
#   
#   
#   # 🛑 CHANGE START: Update to Left-Closed, Right-Open [L, R) logic
#   # for (i in 1:n_intervals) {
#   #   vprev <- intervals_df$start_ps[i]
#   #   vcurr <- intervals_df$end_ps[i]
#   #   
#   #   if (i < n_intervals) {
#   #     # Standard bins: [vprev, vcurr) 
#   #     # Include left boundary exactly, exclude right boundary exactly
#   #     bin_mask <- (data_subset$ps >= vprev) & (data_subset$ps < vcurr)
#   #   } else {
#   #     # Final bin: [vprev, vcurr]
#   #     # Include both to ensure PS = 1.0 is captured
#   #     bin_mask <- (data_subset$ps >= vprev) & (data_subset$ps <= vcurr)
#   #   }
#   #   
#   #   data_subset$bin_index[bin_mask] <- intervals_df$bin_index[i]
#   # }
#   
#   # for (i in 1:nrow(intervals_df)) {
#   #   vprev <- intervals_df$start_ps[i]
#   #   vcurr <- intervals_df$end_ps[i]
#   #   
#   #   
#   #   # TODO: I removed the epsilones, but that should be tested again
#   #   if (has_empty_before[i]) {
#   #     bin_mask <- (data_subset$ps >= vprev - eps) & (data_subset$ps <= vcurr + eps)
#   #   } else {
#   #     bin_mask <- (data_subset$ps > vprev + eps) & (data_subset$ps <= vcurr + eps)
#   #   }
#   #   
#   #   data_subset$bin_index[bin_mask] <- intervals_df$bin_index[i]
#   # }
# 
#   # 🛑 CRITICAL: Store original order and key columns BEFORE merge
#   original_order <- 1:nrow(data_subset)
#   data_subset$original_order <- original_order
#   data_subset$original_id <- data_subset[[data_config$ID_VAR]]
#   data_subset$original_treat <- data_subset[[data_config$TREATMENT_VAR]]
#   data_subset$original_ps <- data_subset$ps
#   
#   # Merge with display positions
#   data_subset <- merge(
#     data_subset, 
#     intervals_df[, c("bin_index", "ps_display_center")], 
#     by = "bin_index",
#     all.x = TRUE
#   )
#   
#   # 🛑 CRITICAL: Restore original order to prevent misalignment
#   data_subset <- data_subset[order(data_subset$original_order), ]
#   
#   # Verify no data loss
#   cat(sprintf("After merge: %d rows (should be same as before)\n", nrow(data_subset)))
#   
#   # Define y-level mapping to separate the three groups vertically
#   # Control -> 0.0 (Bottom), Treated -> 0.8 (Middle), Redundant Treated -> 1.1 (Top)
#   y_level_mapping <- c("Control" = 0.0, "Treated" = 0.8, "Redundant Treated" = 1.1)
#   
#   # Create plot data using the preserved original columns
#   plot_data <- data.frame(
#     id = data_subset$original_id,  # Use preserved ID
#     ps_display = jitter(data_subset$ps_display_center, amount = 0.2), 
#     ps_original = data_subset$original_ps,  # Use preserved PS
#     treat = data_subset$original_treat,  # Use preserved treatment
#     is_redundant = data_subset$original_treat == 1 & 
#       data_subset$original_id %in% redundant_unit_ids,
#     color_group = factor(
#       ifelse(data_subset$original_treat == 0,
#              "Control",
#              ifelse(data_subset$original_id %in% redundant_unit_ids,
#                     "Redundant Treated", 
#                     "Treated")),
#       levels = c("Control", "Treated", "Redundant Treated")
#     ),
#     stringsAsFactors = FALSE
#   )
#   
#   # 🛑 FIX: Calculate the Y-axis center based on the color group and apply jitter
#   plot_data$y_center_value <- y_level_mapping[as.character(plot_data$color_group)]
#   plot_data$treat_jitter <- jitter(plot_data$y_center_value, amount = 0.05)
#   
#   # Create plot with explicit aesthetic mapping
#   p <- ggplot(plot_data) +
#     geom_point(aes(x = ps_display, 
#                    y = treat_jitter, 
#                    color = color_group,
#                    text = paste("PS:", round(ps_original, 4), "<br>Group:", color_group)),
#                alpha = 0.8, 
#                size = 2.5) +
#     scale_color_manual(
#       name = "Group",
#       values = c("Control" = "blue", "Treated" = "red", "Redundant Treated" = "purple"),
#       drop = FALSE
#     ) +
#     labs(title = title, 
#          x = "Propensity Score Bin (Original PS Range)", 
#          y = "", 
#          color = "Group") +
#     theme_minimal() +
#     theme(axis.text.y = element_blank(), 
#           axis.ticks.y = element_blank(),
#           plot.title = element_text(face = "bold", size = 14))
#   
#   cat("ggplot created successfully\n")
#   
#   # Print the ggplot to see if colors appear correctly BEFORE plotly conversion
#   # cat("\n=== Printing static ggplot (before plotly) ===\n")
#   # print(p)
#   
#   # --- 4. Add Vertical Lines and Labels for ALL bins (including empty ones) ---
#   boundaries_display <- intervals_df$ps_display_end[-nrow(intervals_df)]
#   
#   # Add left boundary
#   p <- p + geom_vline(xintercept = 0.5, linetype = "dashed", 
#                       color = "black", linewidth = 0.7, alpha = 0.7)
#   
#   # Add remaining boundaries
#   for (boundary in boundaries_display) {
#     p <- p + geom_vline(xintercept = boundary, linetype = "dashed",
#                         color = "black", linewidth = 0.7, alpha = 0.7)
#   }
#   
#   # 🛑 FIX: Add labels for ALL bins in intervals_df, not just those with data
#   for (i in 1:nrow(intervals_df)) {
#     bin_center <- intervals_df$ps_display_center[i]
#     ps_range_label <- sprintf("[%.3f, %.3f]", intervals_df$start_ps[i], intervals_df$end_ps[i])
#     
#     # Debug
#     # browser()
#     
#     # Main bin label
#     if (intervals_df$is_empty[i]) {
#       label_text <- sprintf("Bin %d\n(Empty)", intervals_df$bin_index[i])
#     } else {
#       label_text <- sprintf("Bin %d\nT=%d, C=%d",
#                             intervals_df$bin_index[i],
#                             intervals_df$n_treated[i],
#                             intervals_df$n_control[i])
#     }
#     
#     # High position label
#     p <- p + annotate("text", x = bin_center, y = 1.25, # Adjusted y-position for label
#                       label = label_text, size = 3,
#                       color = "darkgreen", fontface = "bold")
#     
#     # Low position label (PS range)
#     p <- p + annotate("text", x = bin_center, y = -0.15,
#                       label = ps_range_label, size = 2.5, color = "darkgray")
#   }
#   
#   # --- 5. Final Aesthetics ---
#   # Adjusted y-limits to accommodate the new 1.1 y-center and labels
#   p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#     coord_cartesian(xlim = c(0.5, nrow(intervals_df) + 0.5), ylim = c(-0.2, 1.3))
#   
#   storing_path <- sprintf("%sbins_plot_%s_step_%s_delta_%s.png", base_dir, dataset_name, stage_label, delta_label)
#   ggsave(storing_path, p, width = 10, height = 8, dpi = 300)
#   
#   if (show_plot) {
#     if (make_interactive && requireNamespace("plotly", quietly = TRUE)) {
#       interactive_p <- plotly::ggplotly(p, tooltip = "text")
#       print(interactive_p)
#       return(interactive_p)
#     } else {
#       print(p)
#       return(p)
#     }
#   }
#   
#   return(p)
# }

plot_ate_comparison <- function(results_df, data_config, delta_dp, dataset_name, base_dir) {
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    message("No ATE data available to plot.")
    return(NULL)
  }
  
  # debug
  # browser()
  
  delta_label <- sprintf("%.2f", delta_dp)
  
  ate_title <- sprintf("Comparison of Average Treatment Effects (ATE) (Delta = %s)", delta_label)
  
  ate_plot <- ggplot(results_df, aes(x = reorder(method, ate), y = ate, fill = method)) +
    # Bar representing the ATE point estimate
    
    geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.7) +
    
    # Error bars representing the 95% Confidence Interval
    # CHANGED: 'linewidth' instead of 'size'
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.8) +
    
    # Reference line at zero (No Effect)
    # CHANGED: 'linewidth' instead of 'size'
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    
    labs(
      title = ate_title,
      subtitle = paste("Outcome Variable:", data_config$OUTCOME_VAR, "| 95% Confidence Intervals"),
      x = "Matching Method",
      y = "Estimated Treatment Effect"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank()
    )
  
  
  storing_path <- file.path(base_dir, sprintf("ate_plot_%s_delta_%.2f.png", dataset_name, delta_dp))
  ggsave(storing_path, ate_plot, width = 10, height = 8, dpi = 300)
  
  return(ate_plot)
}

plot_matching_map <- function(match_map, cov_df, x_var = "ps", bin_col = NULL, title = "Match Analysis") {
  
  # s <- sum(is.na(cov_df$ps))
  # 
  # s1 <- sum(map_quintile$treated_id %in% cov_df$id)
  
  # Debug
  # browser()
  
  # 1. Protection & Filtering
  match_map <- match_map %>% filter(!is.na(cost))
  if (nrow(match_map) == 0) {
    message("Warning: No matches to plot for ", title)
    return(NULL)
  }
  
  # 2. Join Logic with Type Safety
  ref_data <- cov_df %>%
    mutate(id_str = as.character(id)) %>%
    dplyr::select(id_str, x_val = !!sym(x_var))
  
  plot_data <- match_map %>%
    mutate(treated_id = as.character(treated_id),
           control_id = as.character(control_id)) %>%
    left_join(ref_data, by = c("treated_id" = "id_str")) %>%
    rename(x_treated = x_val) %>%
    left_join(ref_data, by = c("control_id" = "id_str")) %>%
    rename(x_control = x_val)
  
  # --- DEBUG CHECK: If this triggers, your IDs are not matching ---
  if(all(is.na(plot_data$x_treated)) | all(is.na(plot_data$x_control))) {
    stop("Error: No matching IDs found between map and data. Check if 'id' columns are identical.")
  }
  
  # 3. Plotting
  p <- ggplot(plot_data) +
    # Use linewidth and pmax to prevent Inf/Negative errors
    geom_segment(aes(x = x_treated, xend = x_control, y = 1, yend = 0, 
                     color = cost, 
                     linewidth = 1 / (pmax(cost, 0, na.rm = TRUE) + 0.01)), 
                 alpha = 0.6) +
    geom_point(aes(x = x_treated, y = 1), color = "red", size = 2) +
    geom_point(aes(x = x_control, y = 0), color = "blue", size = 2) +
    scale_color_viridis_c(option = "plasma") +
    scale_linewidth_continuous(range = c(0.2, 2.5), guide = "none") +
    theme_minimal() +
    labs(title = title, 
         subtitle = paste("Matches:", nrow(plot_data)),
         x = paste("Propensity Score (", x_var, ")", sep=""), 
         y = "", 
         color = "Mahalanobis Cost",
         linewidth = "Match Quality") + # Fixed the label warning here
    scale_y_continuous(breaks = c(0, 1), labels = c("Control", "Treated"))
  
  if (!is.null(bin_col) && bin_col %in% colnames(plot_data)) {
    p <- p + geom_text(aes(x = (x_treated + x_control)/2, y = 0.5, label = !!sym(bin_col)), 
                       size = 3, vjust = -1, check_overlap = TRUE)
  }
  
  p <- p + coord_cartesian(clip = "off")
  
  message("Plot generated successfully for: ", title)
  return(p)
}

plot_match_distances <- function(map_data, data_matched, target_dir, outfile_name, plot_title, delta_dp = NULL) {
  
  full_path <- file.path(target_dir, outfile_name)
  
  # 1. PREPARE DATA
  is_osip <- "ps_bin_final" %in% colnames(data_matched)
  
  # Ensure IDs are characters to prevent join type-mismatch errors
  dist_data <- map_data %>%
    mutate(treated_id = as.character(treated_id)) %>%
    left_join(data_matched %>% 
                mutate(id_str = as.character(id)), 
              by = c("treated_id" = "id_str"))
  
  # Create a unified 'plot_group' column
  if (is_osip) {
    dist_data <- dist_data %>%
      filter(!is.na(ps_bin_final)) %>%
      mutate(plot_group = factor(ps_bin_final))
    y_axis_label <- "Osip Bin"
  } else {
    dist_data <- dist_data %>%
      mutate(plot_group = factor("Global Match"))
    y_axis_label <- "Matching Method"
  }
  
  # --- SAFETY CHECK ---
  if (nrow(dist_data) == 0) {
    warning("Skipping plot - No rows found in dist_data after join!")
    return(NULL)
  }
  
  # 2. PLOT LOGIC
  bin_counts <- table(dist_data$plot_group)
  
  if (any(bin_counts < 3)) {
    mode_label <- "(Boxplot Mode)"
    p_match <- ggplot(dist_data, aes(x = cost, y = plot_group)) +
      geom_boxplot(aes(fill = plot_group), alpha = 0.3, outlier.shape = NA) + 
      # Fixed point layer:
      geom_point(aes(color = cost), 
                 position = position_jitter(height = 0.1, width = 0), 
                 alpha = 0.8, size = 2) +
      scale_color_viridis_c(name = "Mahal. Cost", option = "plasma") +
      guides(fill = "none")
  } else {
    mode_label <- ""
    p_match <- ggplot(dist_data, aes(x = cost, y = plot_group, fill = after_stat(x))) +
      ggridges::geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
      scale_fill_viridis_c(name = "Mahal. Cost", option = "plasma") +
      theme_minimal()
  }
  
  # 3. METADATA
  display_subtitle <- paste("Max Mahalanobis:", round(max(dist_data$cost, na.rm=TRUE), 4))
  if (!is.null(delta_dp)) {
    display_subtitle <- paste(display_subtitle, "| Delta:", delta_dp)
  }
  
  # 4. FINAL LABELS
  p_match <- p_match +
    labs(
      title = paste(plot_title, mode_label),
      subtitle = display_subtitle,
      x = "Mahalanobis Distance (Cost)", 
      y = y_axis_label
    )
  
  # 5. SAVE
  ggsave(full_path, p_match, width = 10, height = 8, dpi = 300)
  message("Report exported to: ", full_path)
  
  return(p_match)
}

# plot_match_distances <- function(map_data, data_matched, target_dir, outfile_name, plot_title, delta_dp = NULL) {
#   
#   full_path <- file.path(target_dir, outfile_name)
#   
#   # 1. PREPARE DATA: Handle the Binning vs. Global logic
#   # Determine if we are in 'Levin Mode' or 'Global Mode'
#   is_levin <- "ps_bin_final" %in% colnames(data_matched)
#   
#   dist_data <- map_data %>%
#     mutate(treated_id = as.character(treated_id)) %>%
#     left_join(data_matched %>% 
#                 mutate(id_str = as.character(id)), 
#               by = c("treated_id" = "id_str"))
#   
#   # Create a unified 'plot_group' column
#   if (is_levin) {
#     dist_data <- dist_data %>%
#       filter(!is.na(ps_bin_final)) %>%
#       mutate(plot_group = as.factor(ps_bin_final),
#              y_label = "Levin Bin")
#   } else {
#     # For Genetic, Optimal, Cardinality: Everything goes into one 'Global' group
#     dist_data <- dist_data %>%
#       mutate(plot_group = factor("Global Match"),
#              y_label = "Matching Method")
#   }
#   
#   # --- SAFETY CHECK ---
#   if (nrow(dist_data) == 0) {
#     warning(paste("Skipping plot - No rows found in dist_data!"))
#     return(NULL)
#   }
#   
#   # 2. PLOT LOGIC
#   # Check counts for the 'Heat-Boxplot' vs 'Ridge' decision
#   bin_counts <- table(dist_data$plot_group)
#   
#   # For Global models, usually one group has many points, so it triggers Ridge.
#   # For Levin, sparse bins might trigger Boxplot.
#   if (any(bin_counts < 3)) {
#     mode_label <- "(Boxplot Mode)"
#     p_match <- ggplot(dist_data, aes(x = cost, y = plot_group)) +
#       geom_boxplot(aes(fill = plot_group), alpha = 0.3, outlier.shape = NA) + 
#       geom_jitter(aes(color = cost), height = 0.1, alpha = 0.8, size = 2) +
#       scale_color_viridis_c(name = "Mahal. Cost", option = "plasma") +
#       guides(fill = "none")
#   } else {
#     mode_label <- ""
#     p_match <- ggplot(dist_data, aes(x = cost, y = plot_group, fill = after_stat(x))) +
#       ggridges::geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
#       scale_fill_viridis_c(name = "Mahal. Cost", option = "plasma")
#   }
#   
#   # 3. DYNAMIC METADATA (Levin Delta vs Global)
#   display_subtitle <- paste("Max Mahalanobis:", round(max(dist_data$cost, na.rm=TRUE), 4))
#   if (!is.null(delta_dp)) {
#     display_subtitle <- paste(display_subtitle, "| Delta:", delta_dp)
#   }
#   
#   # 4. FINAL THEME & LABELS
#   p_match <- p_match +
#     theme_minimal() +
#     labs(
#       title = paste(plot_title, mode_label),
#       subtitle = display_subtitle,
#       x = "Mahalanobis Distance (Cost)", 
#       y = ifelse(is_levin, "Levin Bin", "Global Match")
#     )
#   
#   # 5. SAVE
#   ggsave(full_path, p_match, width = 10, height = 8, dpi = 300)
#   message("Report exported to: ", full_path)
#   
#   return(p_match)
# }


plot_rosenbaum <- function(results, alpha = 0.05, title = "Rosenbaum Sensitivity Analysis", save_path = NULL) {
  
  # 1. Ensure the folder exists
  if (!is.null(save_path)) {
    plot_dir <- dirname(save_path)
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  }
  
  # 2. Create the ggplot object
  p <- ggplot(results, aes(x = Gamma, y = P_Value_Bound)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "blue", size = 3) +
    # Significance line
    geom_hline(yintercept = alpha, linetype = "dashed", color = "red", linewidth = 1) +
    # Label for alpha
    annotate("text", x = min(results$Gamma), y = alpha, label = paste("α =", alpha), 
             vjust = -1, color = "red", fontface = "bold") +
    # Dynamic styling
    scale_y_continuous(limits = c(0, max(min(max(results$P_Value_Bound) * 1.2, 1), 0.1))) +
    labs(title = title,
         x = expression(Gamma),
         y = "Upper-bound p-value") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 3. Print the plot to the console/viewer
  print(p)
  
  # 4. Save if a path is provided
  if (!is.null(save_path)) {
    ggsave(
      filename = save_path, 
      plot = p, 
      width = 8, 
      height = 6, 
      dpi = 300 # Higher DPI for JASA-quality figures
    )
    message("✓ Sensitivity plot saved to: ", save_path)
  }
}


# plot_rosenbaum <- function(results, alpha = 0.05, title = "Rosenbaum Sensitivity Analysis", save_path = NULL) {
#   
#   # 1. Ensure the folder exists
#   plot_dir <- dirname(save_path)
#   if (!dir.exists(plot_dir)) {
#     dir.create(plot_dir, recursive = TRUE)
#   }
#   
#   # Generate the plot
#   plot(results$Gamma, results$P_Value_Bound,
#        type = "b",
#        pch = 19,
#        col = "blue",
#        ylim = c(0, max(min(max(results$P_Value_Bound) * 1.2, 1), 0.1)), # Dynamic Y limit
#        xlab = expression(Gamma),
#        ylab = "Upper-bound p-value",
#        main = title)
#   
#   # Add significance threshold line
#   abline(h = alpha, col = "red", lty = 2, lwd = 1.5)
#   
#   # Label the alpha line
#   text(x = min(results$Gamma), y = alpha, 
#        labels = paste("α =", alpha), 
#        pos = 3, col = "red", cex = 0.8)
#   
#   # If saving, close the device
#   if (!is.null(save_path)) {
#     
#     # 2. If it's a ggplot object named 'p'
#     ggsave(
#       filename = save_path, 
#       plot = p, 
#       width = 8,   # ggsave uses inches by default
#       height = 6, 
#       dpi = 120    # This matches your 'res = 120'
#     )
# 
#     message("✓ Sensitivity plot saved to: ", save_path)
#   }
# }

# plot_rosenbaum <- function(results, alpha = 0.05) {
#   
#   plot(results$Gamma, results$P_Value_Bound,
#        type = "b",
#        pch = 19,
#        ylim = c(0, 1),
#        xlab = expression(Gamma),
#        ylab = "Upper-bound p-value",
#        main = "Rosenbaum Sensitivity Analysis")
#   
#   abline(h = alpha, col = "red", lty = 2)
#   text(max(results$Gamma), alpha,
#        labels = "α = 0.05",
#        pos = 3, col = "red")
# }


# Run the diagram
# dist_diag <- plot_match_distances(map_levin, matching_res_levin$data_matched)
# print(dist_diag)

# Generate it
# dist_diag <- plot_match_distances(map_levin, data_subset)
# print(dist_diag)