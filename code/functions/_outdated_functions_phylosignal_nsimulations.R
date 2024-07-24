################### Compute phylogenetic signal for a list of simulation results








### Check what matrices of interactions have all 0 and discard them


check_matrix.with.values <- function(matrix) {
  
  result <- any(matrix != 0)
  
  return(result)
  
}




eliminate_basals <- function(matrix, nbasals) {
  
  mat <- matrix[(nbasals+1):nrow(matrix),]
  
  return(mat)
}




############## FACILIATION AND COMPETITION #################



compute_cor_phylosig_time_comp.fac <- function(list_sim) {
  
  
  res = data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("timesteps",
                                                                "phylo cor mean",
                                                                "phylo sign mean",
                                                                "phylo cor pred",
                                                                "phylo sign pred",
                                                                "phylo cor prey",
                                                                "phylo sign prey",
                                                                "sim")))) #makes an empty dataframe
  
  df_richness = data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("timesteps",
                                                                              "nspp",
                                                                              "sim"))))
  
  
  
  cor_signal_time <- numeric(length(list_sim))
  
  
  
  
  
  for (sim in 1:length(list_sim)) {
    
    
    print(paste("simulation", sim, "of", length(list_sim)))
    
    
    
    path <- list_sim[sim]
    
    list_res <- readRDS(path)
    
    
    # count number of timesteps where there were spp
    
    list_simulation1 <- list_res
    
    n_steps <- length(list_simulation1$network_list)
    
    
    presence_matrix <- list_simulation1$presence_matrix
    
    presence_matrix <- presence_matrix[1:n_steps,]
    
    
    
    
    #### Identify timesteps where phylogenetic distances cant be calculated (those with less than 3 spp)
    
    
    non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
    
    # until what timestep need to discard:
    
    final.discarded_timestep <- non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)]
    
    
    
    
    
    
    
   
    
    
    
    
    
    
    #### homogenize elements to start from valid timesteps
    
    
    # ancestry-distances table
    
    list_anc_dist <- list_simulation1$list_anc_dist[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
    
    # Network list
    
    network_list <- list_simulation1$network_list[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
    
   
    
    
    
    
    ## Convert spp names from numbers to letters
    
    
    ## ancestry-distances table
    
    list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
    
    
    ## Network list 
    
    
    list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
    
    
    #### convert numbers to letters
    
    list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
    
    
    
    
    
    
    
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
    
    
    
    ######### Crop phylogenetic and interaction distance 
    
    
    # Set the same spp names for the presence_matrix than for the interacion matrices
    
    
    colnames(presence_matrix) <- seq(1:Smax)
    
    colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
    
    
    # Discard same timesteps (rows) than the discarted phylogenetic distance matrices
    
    presence_matrix <- presence_matrix[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist),]
    
    
    
    
    
    ## crop the phylogenetic distance matrix with present spp
    
    list_dist.phylo_pres <- list()
    
    for (i in 1:length(list_dist.phylo)) {
      
      
      list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    
    ## Retain only present species in network matrices
    
    
    
    ## crop the interaction matrix with present spp
    
    list_net_present_spp.letters <- list()
    
    for (i in 1:length(list_networks_sppnames_letters)) {
      
      list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    
    
    
    
    
    
    ## Compute interaction distances (NMI)
    
    
    ## DISTANCES IN DONNORS (columns)
    
    
    #list_interact_distances_pred <- lapply(list_net_present_spp.letters, FUN = compute_nmi_cols)
    list_interact_distances_pred <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_pred)
    
    
    # thosa that are all 0 will have NaN - I need to convert them into 0
    
    
    ## DISTANCES AS RECEPTORS (rows)
    
    #list_interact_distances_prey <- lapply(list_net_present_spp.letters, FUN = compute_nmi_rows)
    list_interact_distances_prey <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_prey)
    
    
    
    # Set matrix diagonals to 0
    
    list_interact_distances_pred <- lapply(list_interact_distances_pred, FUN = diag_to0)
    
    list_interact_distances_prey <- lapply(list_interact_distances_prey, FUN = diag_to0)
    
    
    
    
    
    # set all Na to 0 (spp that compared vectors with all 0)
    
    
    list_interact_distances_pred_corrected <- list()
    
    for (i in 1:length(list_interact_distances_pred)) {
      
      list_interact_distances_pred_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_pred[[i]],
                                                                            marg = 2) # 1 (rows), 2 (col), or c(1,2)
      
    }
    
    
    list_interact_distances_prey_corrected <- list()
    
    for (i in 1:length(list_interact_distances_prey)) {
      
      list_interact_distances_prey_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_prey[[i]],
                                                                            marg = 1) # 1 (rows), 2 (col), or c(1,2)
      
    }
    
    
    
    
    
    
    # compute mean distances
    
    list_interact_distances_mean_corrected <- list()
    
    
    for (i in 1:length(list_interact_distances_pred_corrected)) {
      
      #pair_mat <- list(list_interact_distances_pred_corrected[[i]], list_interact_distances_prey_corrected[[i]])
      
      #list_interact_distances_mean_corrected[[i]] <- compute_mean_two_mat_from_list(list = pair_mat)
      
      list_interact_distances_mean_corrected[[i]] <- (list_interact_distances_pred_corrected[[i]] + list_interact_distances_prey_corrected[[i]]) / 2
      
    }
    
    
    
    
    
    
    # Principal Coordinate Analyses
    
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! - do not discard them. 
    
#    vec_timesteps_all0 <- which(lapply(list_interact_distances_pred_corrected,check_matrix.with.values) != TRUE)
    
    
#    if(length(vec_timesteps_all0) > 0){
      
#      n_discard_timesteps_all0 <- vec_timesteps_all0[length(vec_timesteps_all0)]
      
 #     list_interact_distances_pred_corrected <- list_interact_distances_pred_corrected[-(1:n_discard_timesteps_all0)]
      
 #     list_interact_distances_prey_corrected <- list_interact_distances_prey_corrected[-(1:n_discard_timesteps_all0)]
      
  #    list_interact_distances_mean_corrected <- list_interact_distances_mean_corrected[-(1:n_discard_timesteps_all0)]
      
   #   list_dist.phylo_pres <- list_dist.phylo_pres[-(1:n_discard_timesteps_all0)]
      
      

      
    }
    
  #  list_dist_dist.phylo_pres <- list()
  #  list_dist_interact_distances_pred_corrected <- list()
  #  list_dist_interact_distances_prey_corrected <- list()
  #  list_dist_interact_distances_mean_corrected <- list()
   
    
    ## Convert matrix to distance objects
    
    for (i in 1:length(list_dist.phylo_pres)) {
      
      list_dist_dist.phylo_pres[[i]] <- as.matrix(list_dist.phylo_pres[[i]])
      list_dist_interact_distances_pred_corrected[[i]] <- as.matrix(list_interact_distances_pred_corrected[[i]])
      list_dist_interact_distances_prey_corrected[[i]] <- as.matrix(list_interact_distances_prey_corrected[[i]])
      list_dist_interact_distances_mean_corrected[[i]] <- as.matrix(list_interact_distances_mean_corrected[[i]])
      
      
    
      
    }
    
      #  list_dist_dist.phylo_pres <- list_dist.phylo_pres
      # list_dist_interact_distances_pred_corrected <- list_interact_distances_pred_corrected
      # list_dist_interact_distances_prey_corrected <- list_interact_distances_prey_corrected
      # list_dist_interact_distances_mean_corrected <- list_interact_distances_mean_corrected
    
    
    list_pco.phy <- list()
    list_pco.int_pred <- list()
    list_pco.int_prey <- list()
    list_pco.int_mean <- list()
    
   
    
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      
      pco_phy <- cmdscale(list_dist_dist.phylo_pres[[i]], k = 2)
      
      list_pco.phy[[i]] <- pco_phy
      
      pco_pred <- cmdscale(list_dist_interact_distances_pred_corrected[[i]], k = 2)
      
      list_pco.int_pred[[i]] <- pco_pred
      
      
      pco_prey <- cmdscale(list_dist_interact_distances_prey_corrected[[i]], k = 2)
      
      list_pco.int_prey[[i]] <- pco_prey
      
      
      pco_mean <- cmdscale(list_dist_interact_distances_mean_corrected[[i]], k = 2)
      
      list_pco.int_mean[[i]] <- pco_mean
      
      
      
      
    }
    
    
    
    
    
    
    protest_pred <- list()
    procrustes_pred <- list()
    
    protest_prey <- list()
    procrustes_prey <- list()
    
    protest_mean <- list()
    procrustes_mean <- list()
    
   
    
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_pred[[i]] <- protest(list_pco.int_pred[[i]], list_pco.phy[[i]])
      
      procrustes_pred[[i]] <- procrustes(list_pco.int_pred[[i]],list_pco.phy[[i]])
      
      
    }
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_prey[[i]] <- protest(list_pco.int_prey[[i]], list_pco.phy[[i]])
      
      procrustes_prey[[i]] <- procrustes(list_pco.int_prey[[i]],list_pco.phy[[i]])
      
      
    }
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_mean[[i]] <- protest(list_pco.int_mean[[i]], list_pco.phy[[i]])
      
      procrustes_mean[[i]] <- procrustes(list_pco.int_mean[[i]],list_pco.phy[[i]])
      
      
    }
    
    
    
    
    ## Create dataframe results
    
    
    
    
    
    
    protest_pval_pred <- c()
    protest_corr_pred <- c()
    protest_t_pred <- c()
    
    protest_pval_prey <- c()
    protest_corr_prey <- c()
    protest_t_prey <- c()
    
    protest_pval_mean <- c()
    protest_corr_mean <- c()
    protest_t_mean <- c()
    
    
  
    
    for (i in 1:length(protest_mean)) {
      
      protest_pval_pred[i] <- protest_pred[[i]]$signif
      
      protest_corr_pred[i] <-protest_pred[[i]]$t0
      
      protest_t_pred[i] <-mean(protest_pred[[i]]$t)
      
      
      protest_pval_prey[i] <- protest_prey[[i]]$signif
      
      protest_corr_prey[i] <-protest_prey[[i]]$t0
      
      protest_t_prey[i] <-mean(protest_prey[[i]]$t)
      
      
      protest_pval_mean[i] <- protest_mean[[i]]$signif
      
      protest_corr_mean[i] <-protest_mean[[i]]$t0
      
      protest_t_mean[i] <-mean(protest_mean[[i]]$t)
      
      
    
      
      
    }
    
    
    if(length(vec_timesteps_all0) > 0){
      
      n_max_step_discarded <- n_discard_timesteps_all0 + final.discarded_timestep
      
      timesteps <- (n_max_step_discarded+1):n_steps
      
    }else{
      
      timesteps <- (final.discarded_timestep+1):n_steps
    }
    
    
    
    # Save results
    
    df_signal_time <- data.frame(timesteps,
                                 protest_pval_pred,protest_corr_pred,protest_t_pred,
                                 protest_pval_prey,protest_corr_prey,protest_t_prey,
                                 protest_pval_mean,protest_corr_mean,protest_t_mean)
    
    sign_pred <- with(df_signal_time, ifelse(protest_pval_pred < 0.051, 'sign', 'non.sign'))
    sign_prey <- with(df_signal_time, ifelse(protest_pval_prey < 0.051, 'sign', 'non.sign'))
    sign_mean <- with(df_signal_time, ifelse(protest_pval_mean < 0.051, 'sign', 'non.sign'))
   
    
    res.add <- data.frame("timesteps" = timesteps,
                          "phylo cor mean" = protest_corr_mean,
                          "phylo sign mean" = sign_mean,
                          "phylo cor pred" = protest_corr_pred,
                          "phylo sign pred" = sign_pred,
                          "phylo cor prey" = protest_corr_prey,
                          "phylo sign prey" = sign_prey,
                          "sim" = rep(sim, times = length(timesteps))
                          
    )
    
    
    res<-rbind(res,res.add)
    
    
    
    # dataframe spp richness with time
    
    n_spp_time <- c()
    
    presence_matrix_cropped <- list_simulation1$presence_matrix[timesteps,]
    
    for (i in 1:nrow(presence_matrix_cropped)) {
      
      n_spp_time[i] <- length(which(presence_matrix_cropped[i,] == 1))
      
    }
    
    
    df_richness.add <- data.frame("timesteps" = timesteps,
                                  "nspp" = n_spp_time,
                                  "sim" = rep(sim, times = length(timesteps)))
    
    df_richness<-rbind(df_richness,df_richness.add)
    
    list_res <- list(res, df_richness)
    
    return(list_res)
    
  }
  











############## FOODWEB #################



compute_cor_phylosig_time_fw <- function(list_sim) {
  
  
  res = data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("timesteps",
                                                                "phylo cor mean",
                                                                "phylo sign mean",
                                                                "phylo cor pred",
                                                                "phylo sign pred",
                                                                "phylo cor prey",
                                                                "phylo sign prey",
                                                                "sim")))) #makes an empty dataframe
  
  df_richness = data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("timesteps",
                                                                              "nspp",
                                                                              "sim"))))
  
  
  
  cor_signal_time <- numeric(length(list_sim))
  
  
  
  
  
  for (sim in 1:length(list_sim)) {
    
    
    print(paste("simulation", sim, "of", length(list_sim)))
    
    
    
    path <- list_sim[sim]
    
    list_res <- readRDS(path)
    
    
    # count number of timesteps where there were spp
    
    list_simulation1 <- list_res
    
    n_steps <- length(list_simulation1$network_list)
    
    
    presence_matrix <- list_simulation1$presence_matrix
    
    presence_matrix <- presence_matrix[1:n_steps,]
    
    
    
    
    #### Identify timesteps where phylogenetic distances cant be calculated (those with less than 3 spp)
    
    
    non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
    
    # until what timestep need to discard:
    
    final.discarded_timestep <- non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)]
    
    
    
   
    
    
    #### homogenize elements to start from valid timesteps
    
    
    # ancestry-distances table
    
    list_anc_dist <- list_simulation1$list_anc_dist[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
    
    # Network list
    
    network_list <- list_simulation1$network_list[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
   
    
    
    
    
    
    
    
    
    
    #### Eliminate basal species
    
    Sbasals <- length(list_simulation1$basal)
    
    network_list <- lapply(network_list, eliminate_basals, nbasals = Sbasals)
    
    
    ## Convert spp names from numbers to letters
    
    
    ## ancestry-distances table
    
    list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
    
    
    ## Network list 
    
    
    list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
    
    
    #### convert numbers to letters
    
    list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
    
    
    
    
    
    
    
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
    
    
    
    ######### Crop phylogenetic and interaction distance 
    
    
    # Set the same spp names for the presence_matrix than for the interacion matrices
    
    
    colnames(presence_matrix) <- seq(1:Smax)
    
    colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
    
    
    # Discard same timesteps (rows) than the discarted phylogenetic distance matrices
    
    presence_matrix <- presence_matrix[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist),]
    
    
    
    
    
    ## crop the phylogenetic distance matrix with present spp
    
    list_dist.phylo_pres <- list()
    
    for (i in 1:length(list_dist.phylo)) {
      
      
      list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    
    ## Retain only present species in network matrices
    
    
    
    ## crop the interaction matrix with present spp
    
    list_net_present_spp.letters <- list()
    
    for (i in 1:length(list_networks_sppnames_letters)) {
      
      list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    
    
    
 
    
    
    ## Compute interaction distances (NMI)
    
    
    ## DISTANCES IN DONNORS (columns)
    
    
    #list_interact_distances_pred <- lapply(list_net_present_spp.letters, FUN = compute_nmi_cols)
    list_interact_distances_pred <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_pred)
    
    
    # thosa that are all 0 will have NaN - I need to convert them into 0
    
    
    ## DISTANCES AS RECEPTORS (rows)
    
    #list_interact_distances_prey <- lapply(list_net_present_spp.letters, FUN = compute_nmi_rows)
    list_interact_distances_prey <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_prey)
    
    
    
    # Set matrix diagonals to 0
    
    list_interact_distances_pred <- lapply(list_interact_distances_pred, FUN = diag_to0)
    
    list_interact_distances_prey <- lapply(list_interact_distances_prey, FUN = diag_to0)
    
    

    
    
    # set all Na to 0 (spp that compared vectors with all 0)
    
    
    list_interact_distances_pred_corrected <- list()
    
    for (i in 1:length(list_interact_distances_pred)) {
      
      list_interact_distances_pred_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_pred[[i]],
                                                                            marg = 2) # 1 (rows), 2 (col), or c(1,2)
      
    }
    
    
    list_interact_distances_prey_corrected <- list()
    
    for (i in 1:length(list_interact_distances_prey)) {
      
      list_interact_distances_prey_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_prey[[i]],
                                                                            marg = 1) # 1 (rows), 2 (col), or c(1,2)
      
    }
    
    
    
    
    
    
    # compute mean distances
    
    list_interact_distances_mean_corrected <- list()
    
    
    for (i in 1:length(list_interact_distances_pred_corrected)) {
      
      #pair_mat <- list(list_interact_distances_pred_corrected[[i]], list_interact_distances_prey_corrected[[i]])
      
      #list_interact_distances_mean_corrected[[i]] <- compute_mean_two_mat_from_list(list = pair_mat)
      
      list_interact_distances_mean_corrected[[i]] <- (list_interact_distances_pred_corrected[[i]] + list_interact_distances_prey_corrected[[i]]) / 2
      
    }
    
    
    
    
    
    
    
    # compute Principal Coordinate Analyses
    
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    vec_timesteps_all0 <- which(lapply(list_interact_distances_pred_corrected,check_matrix.with.values) != TRUE)
    
    
    if(length(vec_timesteps_all0) > 0){
      
      n_discard_timesteps_all0 <- vec_timesteps_all0[length(vec_timesteps_all0)]
      
      list_interact_distances_pred_corrected <- list_interact_distances_pred_corrected[-(1:n_discard_timesteps_all0)]
      
      list_interact_distances_prey_corrected <- list_interact_distances_prey_corrected[-(1:n_discard_timesteps_all0)]
      
      list_interact_distances_mean_corrected <- list_interact_distances_mean_corrected[-(1:n_discard_timesteps_all0)]
      
      list_dist.phylo_pres <- list_dist.phylo_pres[-(1:n_discard_timesteps_all0)]
      
      
    
      
    }
    
    list_dist_dist.phylo_pres <- list()
    list_dist_interact_distances_pred_corrected <- list()
    list_dist_interact_distances_prey_corrected <- list()
    list_dist_interact_distances_mean_corrected <- list()
    
    
   
    
    ## Convert matrix to distance objects
    
    for (i in 1:length(list_dist.phylo_pres)) {
      
      list_dist_dist.phylo_pres[[i]] <- as.dist(list_dist.phylo_pres[[i]])
      list_dist_interact_distances_pred_corrected[[i]] <- as.dist(list_interact_distances_pred_corrected[[i]])
      list_dist_interact_distances_prey_corrected[[i]] <- as.dist(list_interact_distances_prey_corrected[[i]])
      list_dist_interact_distances_mean_corrected[[i]] <- as.dist(list_interact_distances_mean_corrected[[i]])

    }
    
    
    
    
    list_pco.phy <- list()
    list_pco.int_pred <- list()
    list_pco.int_prey <- list()
    list_pco.int_mean <- list()
    

    
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      
      pco_phy <- dudi.pco(list_dist_dist.phylo_pres[[i]], scannf = FALSE, full = TRUE) 
      
      list_pco.phy[[i]] <- pco_phy$li
      
      pco_pred <- dudi.pco(list_dist_interact_distances_pred_corrected[[i]], scannf = FALSE, full = TRUE)
      
      list_pco.int_pred[[i]] <- pco_pred$li
      
      pco_prey <- dudi.pco(list_dist_interact_distances_prey_corrected[[i]], scannf = FALSE, full = TRUE)
      
      list_pco.int_prey[[i]]  <- pco_prey$li
      
      pco_mean <- dudi.pco(list_dist_interact_distances_mean_corrected[[i]], scannf = FALSE, full = TRUE)
      
      list_pco.int_mean[[i]]  <- pco_mean$li
      

      
      
    }
    
    
    
    
    
    
    protest_pred <- list()
    procrustes_pred <- list()
    
    protest_prey <- list()
    procrustes_prey <- list()
    
    protest_mean <- list()
    procrustes_mean <- list()
    

    
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_pred[[i]] <- protest(list_pco.int_pred[[i]], list_pco.phy[[i]])
      
      procrustes_pred[[i]] <- procrustes(list_pco.int_pred[[i]],list_pco.phy[[i]])
      
      
    }
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_prey[[i]] <- protest(list_pco.int_prey[[i]], list_pco.phy[[i]])
      
      procrustes_prey[[i]] <- procrustes(list_pco.int_prey[[i]],list_pco.phy[[i]])
      
      
    }
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_mean[[i]] <- protest(list_pco.int_mean[[i]], list_pco.phy[[i]])
      
      procrustes_mean[[i]] <- procrustes(list_pco.int_mean[[i]],list_pco.phy[[i]])
      
      
    }
    

    
    
    
    
    ## Create dataframe results
    
    
    
    
    
    
    protest_pval_pred <- c()
    protest_corr_pred <- c()
    protest_t_pred <- c()
    
    protest_pval_prey <- c()
    protest_corr_prey <- c()
    protest_t_prey <- c()
    
    protest_pval_mean <- c()
    protest_corr_mean <- c()
    protest_t_mean <- c()
    
    

    
    
    
    for (i in 1:length(protest_mean)) {
      
      protest_pval_pred[i] <- protest_pred[[i]]$signif
      
      protest_corr_pred[i] <-protest_pred[[i]]$t0
      
      protest_t_pred[i] <-mean(protest_pred[[i]]$t)
      
      
      protest_pval_prey[i] <- protest_prey[[i]]$signif
      
      protest_corr_prey[i] <-protest_prey[[i]]$t0
      
      protest_t_prey[i] <-mean(protest_prey[[i]]$t)
      
      
      protest_pval_mean[i] <- protest_mean[[i]]$signif
      
      protest_corr_mean[i] <-protest_mean[[i]]$t0
      
      protest_t_mean[i] <-mean(protest_mean[[i]]$t)
      
      

      
      
    }
    
    
    
    
    timesteps <- (final.discarded_timestep+1):n_steps
    
    
    
    # Save results
    
    df_signal_time <- data.frame(timesteps,
                                 protest_pval_pred,protest_corr_pred,protest_t_pred,
                                 protest_pval_prey,protest_corr_prey,protest_t_prey,
                                 protest_pval_mean,protest_corr_mean,protest_t_mean)
    
    sign_pred <- with(df_signal_time, ifelse(protest_pval_pred < 0.051, 'sign', 'non.sign'))
    sign_prey <- with(df_signal_time, ifelse(protest_pval_prey < 0.051, 'sign', 'non.sign'))
    sign_mean <- with(df_signal_time, ifelse(protest_pval_mean < 0.051, 'sign', 'non.sign'))
  
    
    
    res.add <- data.frame("timesteps" = timesteps,
                          "phylo cor mean" = protest_corr_mean,
                          "phylo sign mean" = sign_mean,
                          "phylo cor pred" = protest_corr_pred,
                          "phylo sign pred" = sign_pred,
                          "phylo cor prey" = protest_corr_prey,
                          "phylo sign prey" = sign_prey,
                          "sim" = rep(sim, times = length(timesteps))
                          
    )
    
    
    res<-rbind(res,res.add)
    
    
    
    
    # dataframe spp richness with time
    
    n_spp_time <- c()
    
    presence_matrix_cropped <- list_simulation1$presence_matrix[timesteps,]
    
    for (i in 1:nrow(presence_matrix_cropped)) {
      
      n_spp_time[i] <- length(which(presence_matrix_cropped[i,] == 1))
      
    }
    
    
    df_richness.add <- data.frame("timesteps" = timesteps,
                                  "nspp" = n_spp_time,
                                  "sim" = rep(sim, times = length(timesteps)))
    
    df_richness<-rbind(df_richness,df_richness.add)
    
    
  }
  
  
  
  
  list_res <- list(res, df_richness)
  
  return(list_res)
  
  
}