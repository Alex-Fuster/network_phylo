
### Function that takes as input a list of results for each simulation run, and computes the mean diversification rate curves (renders a plot and dataframes with the data)


compute_div_curves_from_list_res <- function(list_simulations) {
  
  
  # Initialize an empty list to store dataframes for each simulation
  df_list <- list()
  
  # Loop through each simulation in list_results_fac
  for (sim_index in seq_along(list_simulations)) {
    res <- list_simulations[[sim_index]]
    
    # count number of timesteps where there are spp
    n_steps <- length(res$network_list)
    
    pres <- res$presence_matrix
    
    # Number of speciation and extinction events
    t0 <- pres[1:(n_steps-1), ]
    t1 <- pres[2:n_steps, ] 
    spec_mat <- pres[1:(n_steps-1), ] * 0
    ext_mat <- pres[1:(n_steps-1), ] * 0
    spec_mat[t1 - t0 == 1] <- 1 
    ext_mat[t1 - t0 == -1] <- 1
    spec <- apply(spec_mat, 1, sum)
    ext <- apply(ext_mat, 1, sum)
    
    S <- apply(pres, 1, sum)[2:n_steps]
    
    time <- 1:(n_steps-1)
    
    # Create dataframe for the simulation
    df_divrates <- data.frame(spec, ext, S, time)
    df_divrates$simulation <- sim_index  # Add simulation index as a column
    
    # Reshape data using melt function from reshape2 package
    df_plot_divrates <- melt(df_divrates, id = c("simulation", "time", "S"))
    
    # Add to the list
    df_list[[sim_index]] <- df_plot_divrates
  }
  
  # Combine all dataframes into a single dataframe
  df_combined <- do.call(rbind, df_list)
  
  
  # summarize mean values between simulations
  
  summary_data <- df_combined %>%
    group_by(S, variable) %>%
    summarise(mean_value = mean(value/S),
              sd_value = sd(value/S),
              se_value = sd_value / sqrt(n()),  # Standard error of the mean
              ci_low = mean_value - 1.96 * se_value,  # Lower bound of 95% CI
              ci_high = mean_value + 1.96 * se_value)  # Upper bound of 95% CI
  
  # Plotting
  plot <- ggplot(summary_data, aes(x = S, y = mean_value, color = variable, group = variable)) +
    geom_smooth(method = "loess") +
    #  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = variable), alpha = 0.3) +
    scale_color_manual(values = c("blue4", "red2")) +
    scale_fill_manual(values = c("blue4", "red2"), guide = FALSE) +
    theme_classic() +
    labs(x = "Species richness", y = "Rate") +
    my.theme
  
  list_result <- list("plot" = plot, 
                      "df_res_sims" = df_combined, 
                      "df_mean_sims" = summary_data)
  
  return(list_result)
  
}
