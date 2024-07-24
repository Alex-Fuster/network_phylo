##########################################################

compute_df_centrality_longevity <- function(results_simulation, int, Smax, sim) {
  
  presence_matrix <- results_simulation$presence_matrix
  
  # number of timesteps
  n_steps <- length(results_simulation$network_list)
  
  #### Crop presence matrix to n_steps
  presence_matrix <- presence_matrix[1:n_steps,]
  
  #### Identify timesteps where phylogenetic distances cant be calculated
  non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
  # until what timestep need to discard:
  final.discarded_timestep <- non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)]
  
  #### homogenize elements to start from valid timesteps
  # ancestry table
  list_anc_dist <- results_simulation$list_anc_dist[(final.discarded_timestep+1):length(results_simulation$list_anc_dist)]
  # network list
  network_list <- results_simulation$network_list[(final.discarded_timestep+1):length(results_simulation$list_anc_dist)]
  
  
  ## ------------------- Check
  if (length(which(is.null(results_simulation$network_list))) > 0) {
    print("PROBLEM - null network somewhere")
  } 
  if(length(list_anc_dist) != length(network_list)){
    print("PROBLEM - length list_anc_dist != length(network_list)")
  }
  ## -------------------
  
  
  
  if(int == "foodweb"){
    #### Eliminate basal species
    Sbasals <- length(results_simulation$basal)
    network_list <- lapply(network_list, eliminate_basals, nbasals = Sbasals)
  }
  
  #### Convert spp names from numbers to letters
  ## ancestry-distances table
  list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
  ## Network list 
  list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
  list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
  
  
  
  ## ------------------- Check
  if(length(list_anc_dist_letters) != length(list_networks_sppnames_letters)){
    print("PROBLEM - list_anc_dist_letters != list_networks_sppnames_letters")
  } 
  ## -------------------
  
  
  ## Loop for obtaining phylogenetic distances:
  list_newick <- list()
  list_trees <- list()
  list_newick_tails <- list()
  list_dist.phylo <- list()
  
  for (i in 1:length(list_anc_dist_letters)) {
    list_newick[[i]] <- ToPhylo(list_anc_dist_letters[[i]])
    list_newick_tails[[i]] <- paste(list_newick[[i]], "root")
    list_trees[[i]] <- read.tree(text = sub("A root",";",list_newick_tails[[i]]))
    list_dist.phylo[[i]] <- cophenetic.phylo(list_trees[[i]])
  }
  
  
  ## ------------------- Check
  if (length(which(unlist(lapply(list_dist.phylo, is.null)) == TRUE)) != 0) {
    print("PROBLEM - nulls in list_dist.phylo")
  }
  if (length(list_dist.phylo) != length(list_networks_sppnames_letters)){
    print("PROBLEM - length(list_dist.phylo) != length(list_networks_sppnames_letters)")
  }
  ## -------------------
  
  
  #### convert spp names to letters in presence matrix
  colnames(presence_matrix) <- seq(1:Smax)
  colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
  
  
  
  ## ------------------- Check
  
  if(length(which(colnames(presence_matrix) != colnames(list_networks_sppnames_letters[[1]]))) != 0){
    print("PROBLEM - spp names in presence_matrix dont correspond to spp names in list_networks_sppnames_letterst")
  }
  ## -------------------
  
  
  
  # Discard same timesteps (rows) than the discarted phylogenetic distance matrices
  presence_matrix <- presence_matrix[(final.discarded_timestep+1):length(results_simulation$list_anc_dist),]
  
  
  
  ## ------------------- Check
  
  if(nrow(presence_matrix) != length(list_networks_sppnames_letters) | 
     nrow(presence_matrix) != length(list_dist.phylo)){
    print("PROBLEM - presence matrix, list phylo dist and list interaction networks dont have the same n_steps")
  }
  ## -------------------
  
  
  #### Retain only present species in phylogenetic distance matrices
  list_dist.phylo_pres <- list()
  
  for (i in 1:length(list_dist.phylo)) {
    list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
  }
  
  
  ## ------------------- Check
  #Check that phylogenetic distance matrices retain present species:
  vec_error <- c()
  for (i in 1:nrow(presence_matrix)) {
    vec_tf <-  names(presence_matrix[i, which(presence_matrix[i,] == 1)]) == colnames(list_dist.phylo_pres[[i]])
    if(length(which(vec_tf == FALSE)) > 0){
      vec_error[i] <- "error"
    } else if (length(which(vec_tf == FALSE)) == 0 ){
      vec_error[i] <- "g"
    }
  }
  if(length(which(vec_error == "error")) > 0){
    print("PROBLEM - list_dist.phylo_pres dont retain present species")
  }
  ## -------------------
  
  
  #### Retain only present species in network matrices
  list_net_present_spp.letters <- list()
  
  for (i in 1:length(list_networks_sppnames_letters)) {
    list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
  }
  
  
  
  
  #########################################################
  #### compute species degree centrality and longevity ####
  
  # note timesteps for extinction and speciation, and longevity
  events <- track_events(list_anc_dist_letters)
  
  # compute in, out and total relative degree per species at each timestep
  degree_centrality_list <- lapply(list_net_present_spp.letters, compute_degree_centrality_igraph)
  
  # compute relative degree
  relative_degree_centrality_list <- lapply(degree_centrality_list, compute_relative_degree_centrality)
  
  
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
  
  events$simulation <- rep(sim, times = nrow(events))
  
  
  return(events)
  
}
