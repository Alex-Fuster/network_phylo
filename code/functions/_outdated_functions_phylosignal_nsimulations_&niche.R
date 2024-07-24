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
  
  
  res = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("timesteps",
                                                                "phylo cor niche",
                                                                "phylo sign niche",
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
    
    
    
    
    
    
    
    ### obtain list of traits per timestep
    
    niche_vec <- list_simulation1$traits[,1]
    
    # names niche vector
    
    names(niche_vec) <- seq(1:1000)
    names(niche_vec) <- chartr("0123456789", "ABCDEFGHIJ", names(niche_vec))
    
    
    # names presence matrix
    colnames(presence_matrix) <- seq(1:1000)
    
    colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
    
    list_niche_vec <- list()
    
    for (i in 1:nrow(presence_matrix)) {
      
      list_niche_vec[[i]] <- niche_vec[names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    
    
    
    
    #### homogenize elements to start from valid timesteps
    
    
    # ancestry-distances table
    
    list_anc_dist <- list_simulation1$list_anc_dist[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
    
    # Network list
    
    network_list <- list_simulation1$network_list[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
    
    # vector niche
    
    list_niche_vec <- list_niche_vec[(final.discarded_timestep+1):length(list_niche_vec)]
    
    
    
    
    
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
    
    
    
    
    
    ## Retain only present species in network matrices
    
    
    
    # Set the same spp names for the presence_matrix than for the interacion matrices
    
    
    colnames(presence_matrix) <- seq(1:1000)
    
    colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
    
    
    # Discard same timesteps (rows) than the discarted phylogenetic distance matrices
    
    presence_matrix <- presence_matrix[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist),]
    
    
    
    
    ## crop the interaction matrix with present spp
    
    list_net_present_spp.letters <- list()
    
    for (i in 1:length(list_networks_sppnames_letters)) {
      
      list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    
    
    ## compute niche distances
    
    
    list_niche_dist <- list()
    
    for (i in 1:length(list_niche_vec)) {
      
      list_niche_dist[[i]] <- abs(outer(list_niche_vec[[i]],list_niche_vec[[i]],'-'))
      
    } 
    
    
    
    
  
    ## For phylogenetic distance matrices, retain only those species present
    
    
    list_dist.phylo_pres <- list()
    
    for (i in 1:length(list_dist.phylo)) {
      
      spp_present <- colnames(list_net_present_spp.letters[[i]])
      list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][spp_present ,spp_present] 
      
      
    }
    
    
    
    
    
  
   
    
    # make sure order of names are correct
    
    for (i in length(list_dist.phylo_pres)) {
      
      
      colnames(list_niche_dist[[i]])[order(match(colnames(list_dist.phylo_pres[[i]]),colnames(list_niche_dist[[i]])))]
    }
    
    # compute MDS
    
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    vec_timesteps_all0 <- which(lapply(list_niche_dist,check_matrix.with.values) != TRUE)
    
    
    if(length(vec_timesteps_all0) > 0){
      
      n_discard_timesteps_all0 <- vec_timesteps_all0[length(vec_timesteps_all0)]
      
      
      list_dist.phylo_pres <- list_dist.phylo_pres[-(1:n_discard_timesteps_all0)]
      
      list_niche_dist <- list_niche_dist[-(1:n_discard_timesteps_all0)]
      
    }
    
    list_dist_dist.phylo_pres <- list()
    
    list_dist_niche_dist <- list()
    
    ## Convert matrix to distance objects
    
    for (i in 1:length(list_dist.phylo_pres)) {
      
      list_dist_dist.phylo_pres[[i]] <- as.dist(list_dist.phylo_pres[[i]])
      
      list_dist_niche_dist[[i]] <- as.dist(list_niche_dist[[i]])
      
    }
    
    
    
    
    list_pco.phy <- list()
    
    list_pco.niche <- list()
    
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      
      pco_phy <- dudi.pco(list_dist_dist.phylo_pres[[i]], scannf = FALSE, full = TRUE) 
      
      list_pco.phy[[i]] <- pco_phy$li
      
      
      
      pco_niche <- dudi.pco(list_dist_niche_dist[[i]], scannf = FALSE, full = TRUE) 
      
      list_pco.niche[[i]] <- pco_niche$li
      
      
    }
    
    
    
    
    
    protest_niche <- list()
    procrustes_niche <- list()
    
  
    
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_niche[[i]] <- protest(list_pco.niche[[i]], list_pco.phy[[i]])
      
      procrustes_niche[[i]] <- procrustes(list_pco.niche[[i]],list_pco.phy[[i]])
      
      
    }
    
    
    
    
    ## Create dataframe results
    
    
    

    
    
    protest_niche_pval <- c()
    protest_niche_corr <- c()
    protest_niche_t <- c()
    
    for (i in 1:length(protest_niche)) {
      
    
      
      protest_niche_pval[i] <- protest_niche[[i]]$signif
      
      protest_niche_corr[i] <-protest_niche[[i]]$t0
      
      protest_niche_t[i] <-mean(protest_niche[[i]]$t)
      
      
    }
    
    
    if(length(vec_timesteps_all0) > 0){
      
      n_max_step_discarded <- n_discard_timesteps_all0 + final.discarded_timestep
      
      timesteps <- (n_max_step_discarded+1):n_steps
      
    }else{
      
      timesteps <- (final.discarded_timestep+1):n_steps
    }
    
    
    
    # Save results
    
    df_signal_time <- data.frame(timesteps,
                                 protest_niche_pval,protest_niche_corr,protest_niche_t)
    
    sign_niche <- with(df_signal_time, ifelse(protest_niche_pval < 0.051, 'sign', 'non.sign'))
    
    
    res.add <- data.frame("timesteps" = timesteps,
                          "phylo cor niche" = protest_niche_corr,
                          "phylo sign niche" = sign_niche,
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










############## FOODWEB #################



compute_cor_phylosig_time_fw <- function(list_sim) {
  
  
  
  res = data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("timesteps",
                                                               "phylo cor niche",
                                                               "phylo sign niche",
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
    
    
    
    ### obtain list of traits per timestep
    
    niche_vec <- list_simulation1$traits[,1]
    
    # names niche vector
    
    names(niche_vec) <- seq(1:1000)
    names(niche_vec) <- chartr("0123456789", "ABCDEFGHIJ", names(niche_vec))
    
    
    # names presence matrix
    colnames(presence_matrix) <- seq(1:1000)
    
    colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
    
    list_niche_vec <- list()
    
    for (i in 1:nrow(presence_matrix)) {
      
      list_niche_vec[[i]] <- niche_vec[names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    #### homogenize elements to start from valid timesteps
    
    
    # ancestry-distances table
    
    list_anc_dist <- list_simulation1$list_anc_dist[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
    
    # Network list
    
    network_list <- list_simulation1$network_list[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist)]
    
    # vector niche
    
    list_niche_vec <- list_niche_vec[(final.discarded_timestep+1):length(list_niche_vec)]
    
    
    
    
    
    
    
    
    
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
    
    
    
    
    
    ## Retain only present species in network matrices
    
    
    
    # Set the same spp names for the presence_matrix than for the interacion matrices
    
    
    colnames(presence_matrix) <- seq(1:1000)
    
    colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
    
    
    # Discard same timesteps (rows) than the discarted phylogenetic distance matrices
    
    presence_matrix <- presence_matrix[(final.discarded_timestep+1):length(list_simulation1$list_anc_dist),]
    
    
    
    
    ## crop the interaction matrix with present spp
    
    list_net_present_spp.letters <- list()
    
    for (i in 1:length(list_networks_sppnames_letters)) {
      
      list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
      
    }
    
    
    
    
    
    ## compute niche distances
    
    list_niche_dist <- list()
    
    for (i in 1:length(list_niche_vec)) {
      
      list_niche_dist[[i]] <- abs(outer(list_niche_vec[[i]],list_niche_vec[[i]],'-'))
      
    } 
    
    
    
    
    
    ## For phylogenetic distance matrices, retain only those species present
    
    
    
    list_dist.phylo_pres <- list()
    
    for (i in 1:length(list_dist.phylo)) {
      
      spp_present <- colnames(list_net_present_spp.letters[[i]])
      list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][spp_present ,spp_present] 
      
      
    }
    
    
    

    # make sure order of names are correct
    
    for (i in length(list_dist.phylo_pres)) {
      
      
      colnames(list_niche_dist[[i]])[order(match(colnames(list_dist.phylo_pres[[i]]),colnames(list_niche_dist[[i]])))]
    }
    
    
    
    
    
    
    # compute MDS
    
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    # identify what matrices have all interaction distances = 0 and discard them
    
    vec_timesteps_all0 <- which(lapply(list_niche_dist,check_matrix.with.values) != TRUE)
    
    
    if(length(vec_timesteps_all0) > 0){
      
      n_discard_timesteps_all0 <- vec_timesteps_all0[length(vec_timesteps_all0)]
      
      
      list_dist.phylo_pres <- list_dist.phylo_pres[-(1:n_discard_timesteps_all0)]
      
      list_niche_dist <- list_niche_dist[-(1:n_discard_timesteps_all0)]
      
    }
    
    list_dist_dist.phylo_pres <- list()
    
    list_dist_niche_dist <- list()
    
    ## Convert matrix to distance objects
    
    for (i in 1:length(list_dist.phylo_pres)) {
      
      list_dist_dist.phylo_pres[[i]] <- as.dist(list_dist.phylo_pres[[i]])
      
      list_dist_niche_dist[[i]] <- as.dist(list_niche_dist[[i]])
      
    }
    
    
    
    
    list_pco.phy <- list()
    
    list_pco.niche <- list()
    
    
    for (i in 1:length(list_dist_dist.phylo_pres)) {
      
      pco_phy <- dudi.pco(list_dist_dist.phylo_pres[[i]], scannf = FALSE, full = TRUE) 
      
      list_pco.phy[[i]] <- pco_phy$li
      
      
      
      pco_niche <- dudi.pco(list_dist_niche_dist[[i]], scannf = FALSE, full = TRUE) 
      
      list_pco.niche[[i]] <- pco_niche$li
      
      
    }
    
    
    
    
    protest_niche <- list()
    procrustes_niche <- list()
    
    
    
    
    for (i in 1:length(list_pco.phy)) {
      
      protest_niche[[i]] <- protest(list_pco.niche[[i]], list_pco.phy[[i]])
      
      procrustes_niche[[i]] <- procrustes(list_pco.niche[[i]],list_pco.phy[[i]])
      
      
    }
    
    
    
    
    ## Create dataframe results
    
    
    
    
    
    
    protest_niche_pval <- c()
    protest_niche_corr <- c()
    protest_niche_t <- c()
    
    for (i in 1:length(protest_niche)) {
      
      
      
      protest_niche_pval[i] <- protest_niche[[i]]$signif
      
      protest_niche_corr[i] <-protest_niche[[i]]$t0
      
      protest_niche_t[i] <-mean(protest_niche[[i]]$t)
      
      
    }
    
    
    if(length(vec_timesteps_all0) > 0){
      
      n_max_step_discarded <- n_discard_timesteps_all0 + final.discarded_timestep
      
      timesteps <- (n_max_step_discarded+1):n_steps
      
    }else{
      
      timesteps <- (final.discarded_timestep+1):n_steps
    }
    
    
    
    
    # Save results
    
    df_signal_time <- data.frame(timesteps,
                                 protest_niche_pval,protest_niche_corr,protest_niche_t)
    
    sign_niche <- with(df_signal_time, ifelse(protest_niche_pval < 0.051, 'sign', 'non.sign'))
    
    
    res.add <- data.frame("timesteps" = timesteps,
                          "phylo cor niche" = protest_niche_corr,
                          "phylo sign niche" = sign_niche,
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


