


track_events <- function(ancestry_list) {
  species_events <- data.frame(spp = character(), timestep = integer(), event = character(), ancestor = character(), stringsAsFactors = FALSE)
  all_spp <- unique(unlist(lapply(ancestry_list, function(x) x$spp)))
  
  for (sp in all_spp) {
    found_speciation <- FALSE
    found_extinction <- FALSE
    
    for (t in 1:length(ancestry_list)) {
      ancestry_table <- ancestry_list[[t]]
      
      if (sp %in% ancestry_table$spp) {
        sp_row <- ancestry_table[ancestry_table$spp == sp, ]
        
        if (!found_speciation) {
          species_events <- rbind(species_events, data.frame(spp = sp, timestep = t, event = "speciation", ancestor = sp_row$ancestor, stringsAsFactors = FALSE))
          found_speciation <- TRUE
        }
        
        if (!found_extinction && sp_row$`A/E` == "E") {
          species_events <- rbind(species_events, data.frame(spp = sp, timestep = t, event = "extinction", ancestor = NA, stringsAsFactors = FALSE))
          found_extinction <- TRUE
        }
      }
    }
  }
  

  
  # Calculate longevity
  longevity <- data.frame(spp = unique(species_events$spp), longevity = integer(length(unique(species_events$spp))), stringsAsFactors = FALSE)
  
  for (sp in unique(species_events$spp)) {
    speciation_time <- species_events[species_events$spp == sp & species_events$event == "speciation", "timestep"]
    extinction_time <- species_events[species_events$spp == sp & species_events$event == "extinction", "timestep"]
    
    if (length(speciation_time) > 0 && length(extinction_time) > 0) {
      longevity[longevity$spp == sp, "longevity"] <- extinction_time - speciation_time
    } else {
      # Handle cases where extinction time is missing
      longevity[longevity$spp == sp, "longevity"] <- NA
    }
  }
  

  
  # Calculate number of speciations
  speciation_counts <- data.frame(spp = unique(species_events$spp), num_speciations = integer(length(unique(species_events$spp))), stringsAsFactors = FALSE)
  
  for (sp in unique(species_events$spp)) {
    num_speciations <- sum(species_events$ancestor == sp & species_events$event == "speciation")
    speciation_counts[speciation_counts$spp == sp, "num_speciations"] <- num_speciations
  }
  

  
  species_events <- merge(species_events, longevity, by = "spp", all.x = TRUE)
  species_events <- merge(species_events, speciation_counts, by = "spp", all.x = TRUE)
  
  return(species_events)
}

events <- track_events(list_anc_dist_letters)
print(events)

######################################################

# Example function to compute degree centrality using igraph
compute_degree_centrality_igraph <- function(interaction_matrix) {
  # Create graph object from interaction matrix
  g <- graph_from_adjacency_matrix(as.matrix(interaction_matrix), mode = "directed")
  
  # Compute degree centrality measures
  degree_centrality <- degree(g, mode = "all")  # Total degree (in + out)
  in_degree <- degree(g, mode = "in")  # In-degree
  out_degree <- degree(g, mode = "out")  # Out-degree
  
  # Extract species names
  species <- V(g)$name
  
  # Create dataframe for degree centrality
  degree_centrality_df <- data.frame(
    spp = species,
    in_degree = in_degree,
    out_degree = out_degree,
    total_degree = degree_centrality,
    stringsAsFactors = FALSE
  )
  
  return(degree_centrality_df)
}

# Example: Apply the function across timesteps
degree_centrality_list <- lapply(list_net_present_spp.letters, compute_degree_centrality_igraph)


#######################################################

# Function to compute relative degree centrality
compute_relative_degree_centrality <- function(degree_centrality_df) {
  max_in_degree <- max(degree_centrality_df$in_degree)
  max_out_degree <- max(degree_centrality_df$out_degree)
  max_total_degree <- max(degree_centrality_df$total_degree)
  
  
  degree_centrality_df$in_relative_degree <- degree_centrality_df$in_degree / max_in_degree
  degree_centrality_df$out_relative_degree <- degree_centrality_df$out_degree / max_out_degree
  degree_centrality_df$relative_degree <- degree_centrality_df$total_degree / max_total_degree
  
  return(degree_centrality_df)
}

# Apply across the list of degree centrality dataframes
relative_degree_centrality_list <- lapply(degree_centrality_list, compute_relative_degree_centrality)


########################################################

# Initialize an empty list to store species dataframes
species_data <- list()

# Iterate over timesteps and handle NA values
species_data <- lapply(seq_along(relative_degree_centrality_list), function(i) {
  timestep <- i + 1  # Assuming timestep index starts from 0 or 1
  
  # Extract dataframe for current timestep
  current_relative_degree <- relative_degree_centrality_list[[i]]
  
  # Add timestep information
  current_relative_degree$timestep <- timestep
  
  # Remove NA values
  current_relative_degree <- na.omit(current_relative_degree)
  
  return(current_relative_degree)
})

# Combine into a single dataframe
relative_degree_df <- do.call(rbind, species_data)




# Compute mean degrees per species
mean_degrees <- relative_degree_df %>%
  group_by(spp) %>%
  summarize(
    mean_in_rel_degree = mean(in_relative_degree ),
    mean_out_rel_degree = mean(out_relative_degree ),
    mean_total_rel_degree = mean(relative_degree),
    mean_in_degree = mean(in_degree),
    mean_out_degree = mean(out_degree ),
    mean_total_degree = mean(total_degree),
    
    .groups = 'drop'  # Use 'drop' to avoid grouped data in subsequent operations
  )



events <- merge(events, mean_degrees, by = "spp", all.x = TRUE)




# Define the timestep segments
max_timestep <- max(events$timestep)
timestep_segments <- seq(0, max_timestep, by = 80)
timestep_segments <- c(timestep_segments, max_timestep)  # Ensure to include the last segment

# Create a factor for timestep segments
events$timestep_segment <- cut(events$timestep, breaks = timestep_segments, labels = FALSE)

# Plotting using ggplot2
ggplot(events, aes(x = mean_total_rel_degree, y = longevity)) +
  geom_point() +
  geom_smooth(method = "gam", se = FALSE, color = "blue") +  # Add linear regression line
  facet_wrap(~ timestep_segment, scales = "free") +
  labs(x = "Mean Relative Degree Centrality", y = "Longevity") +
  theme_minimal()

ggplot(events, aes(x = mean_in_rel_degree , y = longevity)) +
  geom_point() +
  geom_smooth(method = "gam", se = FALSE, color = "blue") +  # Add linear regression line
  facet_wrap(~ timestep_segment, scales = "free") +
  labs(x = "Mean Relative Degree Centrality", y = "Longevity") +
  theme_minimal()

ggplot(events, aes(x = mean_out_rel_degree  , y = longevity)) +
  geom_point() +
  geom_smooth(method = "gam", se = FALSE, color = "blue") +  # Add linear regression line
  facet_wrap(~ timestep_segment, scales = "free") +
  labs(x = "Mean Relative Degree Centrality", y = "Longevity") +
  theme_minimal()


##### total (no segments)

ggplot(events, aes(x = mean_out_rel_degree  , y = longevity)) +
  geom_point() +
  geom_smooth(method = "gam", se = FALSE, color = "blue") +  # Add linear regression line
  labs(x = "Mean Relative Degree Centrality", y = "Longevity") +
  theme_minimal()

ggplot(events, aes(x = mean_in_rel_degree  , y = longevity)) +
  geom_point() +
  geom_smooth(method = "gam", se = FALSE, color = "blue") +  # Add linear regression line
  labs(x = "Mean Relative Degree Centrality", y = "Longevity") +
  theme_minimal()



# Assuming your data frame is named events with columns mean_out_rel_degree and longevity
# Adjust these values based on your actual data distribution
low_threshold <- 0.25
high_threshold <- 0.75

# Categorize mean_out_rel_degree
events$centrality_category <- cut(events$mean_in_rel_degree,
                                  breaks = c(-Inf, low_threshold, high_threshold, Inf),
                                  labels = c("low", "medium", "high"))

# Create a boxplot
ggplot(events, aes(x = centrality_category, y = longevity)) +
  geom_boxplot() +
  labs(x = "Mean Relative Degree Centrality", y = "Longevity") +
  theme_minimal()