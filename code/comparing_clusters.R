
library(igraph)     
library(cluster)    
library(fpc)        


k_clusters = 5


# Example data for timestep 20 (interaction matrix and phylogenetic correlation matrix)
interaction_matrix <- list_net_present_spp.letters[[20]]
phylo_corr_matrix <- list_phylo.corr_cropped[[20]]



# --------------------- Interaction-based clustering ------------------------


###### Method 1

# k-means clustering on the interaction matrix
set.seed(123) 
k_clusters_interaction <- kmeans(interaction_matrix, centers = k_clusters)  # Assuming 3 groups
interaction_clusters <- k_clusters_interaction$cluster

###### Method 2

# interaction_matrix <- list_net_present_spp.letters[[20]]
# g <- graph.adjacency(interaction_matrix, mode = "undirected", diag = FALSE)
# 
# # Leading eigenvector clustering
# clusters_egein <- cluster_leading_eigen(g)
# interaction_clusters <- membership(clusters_egein)


# --------------------- Phylogeny-based clustering ------------------------

# Hierarchical clustering on the phylogenetic correlation matrix
phylo_dist <- as.dist(1 - phylo_corr_matrix)  # Convert correlation to distance
hc_phylo <- hclust(phylo_dist, method = "ward.D2")

# Cut the phylogenetic tree into k clusters
phylo_clusters <- cutree(hc_phylo, k = k_clusters)


# --------------------- Compare cluster assignments ------------------------

# Adjusted Rand Index (ARI)
ari_value <- cluster.stats(d = phylo_dist, phylo_clusters, interaction_clusters)$corrected.rand

# Normalized Mutual Information (NMI)
nmi_value <- compare(as.integer(interaction_clusters), as.integer(phylo_clusters), method = "nmi")

# Print results
cat("Adjusted Rand Index (ARI):", ari_value, "\n")
cat("Normalized Mutual Information (NMI):", nmi_value, "\n")

# ---------------------------------------------

# Dendrogram for phylogenetic clustering
plot(hc_phylo, labels = colnames(phylo_corr_matrix), main = "Phylogenetic Clustering")
rect.hclust(hc_phylo, k = k_clusters, border = "red")




















######################


compute_cluster_metrics <- function(results_simulation, Smax, nbasals, num_clusters = 3) {
  
  
  presence_matrix <- results_simulation$presence_matrix
  
  # Number of timesteps
  n_steps <- length(results_simulation$network_list)
  
  # Identify timesteps where phylogenetic distances can't be calculated
  non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
  final.discarded_timestep <- (non.valid_timesteps_phylo_distance[length(non.valid_timesteps_phylo_distance)]) + 1
  
  # Homogenize elements to start from valid timesteps
  list_anc_dist <- results_simulation$list_anc_dist[(final.discarded_timestep + 1):n_steps]
  network_list <- results_simulation$network_list[(final.discarded_timestep + 1):n_steps]
  
  # Convert species names from numbers to letters
  list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
  list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
  list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
  
  colnames(presence_matrix) <- seq(1:Smax)
  colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
  
  presence_matrix <- presence_matrix[(final.discarded_timestep + 1):n_steps, ]
  
  # Initialize lists to store clustering results
  ari_values <- numeric(length(list_networks_sppnames_letters))
  nmi_values <- numeric(length(list_networks_sppnames_letters))
  
  for (i in 1:length(list_anc_dist_letters)) {
    
    # Process phylogenetic data
    newick <- ToPhylo2(list_anc_dist_letters[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root",";", newick_tail))
    
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    phylo.vcv <- vcv(tree)
    phylo.corr <- cov2cor(phylo.vcv)
    
    alive_species <- names(which(presence_matrix[i, ] == 1))
    phylo_corr_cropped <- phylo.corr[alive_species, alive_species]
    
    # Clustering for phylogeny using hierarchical clustering
    phylo_dist <- as.dist(1 - phylo_corr_cropped)  # Convert correlation to distance
    phylo_clusters <- cutree(hclust(phylo_dist, method = "ward.D2"), k = num_clusters)
    
    # Clustering for interactions using k-means
    interaction_matrix <- list_networks_sppnames_letters[[i]][alive_species, alive_species]
    kmeans_result <- kmeans(interaction_matrix, centers = num_clusters, nstart = 25)
    interaction_clusters <- kmeans_result$cluster
    
    # Compute ARI and NMI
    ari_values[i] <- cluster.stats(d = phylo_dist, phylo_clusters, interaction_clusters)$corrected.rand
    nmi_values[i] <- compare(interaction_clusters, phylo_clusters, method = "nmi")
  }
  
  return(data.frame(
    timestep = 1:length(ari_values),
    ARI = ari_values,
    NMI = nmi_values
  ))
}


cluster_results <- compute_cluster_metrics(results_simulation = res_sim, 
                                           Smax = pars$Smax, 
                                           nbasals = pars$Sbasal, 
                                           num_clusters = 3)

# View the results
head(cluster_results)
