# --- plots.R ---

plot_ps_distribution <- function(df, dataset_name, treatment_col, datasets, base_dir, title_suffix = "") {
  # Map dataset name to the correct column
  # Using the global datasets list logic
  # treatment_col <- data_config$TREATMENT_VAR
  
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

plot_ps_distribution_trimmed <- function(df, dataset_name, treatment_col, datasets, base_dir, title_suffix = "", is_trimmed = FALSE) {
  
  df$treat_label <- ifelse(df[[treatment_col]] == 1, "Treated", "Control")
  df$treat_jitter <- jitter(as.numeric(df[[treatment_col]]), amount = 0.05)
  
  # Determine Plot Boundaries
  # If trimmed, we zoom into the actual data range plus a small buffer
  if (is_trimmed) {
    x_min <- min(df$ps, na.rm = TRUE) - 0.02
    x_max <- max(df$ps, na.rm = TRUE) + 0.02
    file_tag <- "trimmed"
  } else {
    x_min <- 0
    x_max <- 1.0
    file_tag <- "full"
  }
  
  # 1. GGPLOT Version (for Saving)
  p <- ggplot(df, aes(x = ps, y = treat_jitter, color = treat_label)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Control" = "blue", "Treated" = "red")) +
    labs(
      title = paste("PS Distribution:", dataset_name, title_suffix),
      x = "Propensity Score",
      y = "",
      color = "Group"
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(x_min, x_max)) + # Use coord_cartesian to zoom without dropping data
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  # Save with a unique name so you don't overwrite the original
  storing_path <- sprintf("%sps_distribution_%s_%s.png", base_dir, dataset_name, file_tag)
  ggsave(storing_path, plot = p, width = 10, height = 6, dpi = 300)
  
  # 2. PLOTLY Version (for Interactive Display)
  fig <- plot_ly(data = df, x = ~ps, y = ~treat_jitter, type = "scatter", mode = "markers",
                 color = ~treat_label, colors = c("blue", "red")) %>%
    layout(
      title = paste("PS Distribution:", dataset_name, title_suffix),
      xaxis = list(title = "Propensity Score", range = c(x_min, x_max)),
      yaxis = list(title = "", showticklabels = FALSE)
    )
  
  print(fig)
}



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

plot_ate_comparison <- function(results_df, data_config, delta_dp, dataset_name, base_dir) {
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    message("No ATE data available to plot.")
    return(NULL)
  }
  
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
      dpi = 300 
    )
    message("✓ Sensitivity plot saved to: ", save_path)
  }
}
