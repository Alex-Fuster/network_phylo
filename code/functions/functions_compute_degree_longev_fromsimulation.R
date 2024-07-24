############# Functions to compute degree centrality from a simulation output

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




compute_degree_centrality_igraph <- function(interaction_matrix) {
  # Create graph object from interaction matrix
  g <- graph_from_adjacency_matrix(as.matrix(interaction_matrix), mode = "directed")
  
  # Compute degree centrality measures
  degree_centrality <- igraph::degree(g, mode = "all")  # Total degree (in + out)
  in_degree <- igraph::degree(g, mode = "in")  # In-degree
  out_degree <- igraph::degree(g, mode = "out")  # Out-degree
  
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