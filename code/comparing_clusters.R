
library(igraph)     
library(cluster)    
library(fpc)        
library(tidyr)


k_clusters = 10


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











#####################

# testing with 1 timestep

# Example data for timestep 20 (interaction matrix and phylogenetic correlation matrix)
interaction_matrix <- list_net_present_spp.letters[[1]]
phylo_corr_matrix <- list_phylo.corr_cropped[[1]]


# Step 1: Clustering interaction matrix using SBM (Stochastic Block Model)
sbm_fit <- sbm::estimateSimpleSBM(interaction_matrix, model = "bernoulli")
interaction_clusters <- sbm_fit$memberships

# Step 2: Clustering phylogenetic correlation matrix using k-means with optimal clusters
phylo_dist <- as.dist(1 - phylo_corr_matrix)  # Convert correlation to distance

# Use silhouette analysis to determine optimal number of clusters
sil_widths <- numeric()
for (k in 2:5) {
  clust <- kmeans(phylo_corr_matrix, centers = k, nstart = 25)
  sil_widths[k] <- mean(silhouette(clust$cluster, phylo_dist)[, 3])
}
optimal_k <- which.max(sil_widths)

# Perform k-means clustering with optimal number of clusters
phylo_clusters <- kmeans(phylo_corr_matrix, centers = optimal_k, nstart = 25)$cluster

# Step 3: Compare clusters using ARI and NMI
ari_value <- cluster.stats(phylo_dist, phylo_clusters, interaction_clusters)$corrected.rand
nmi_value <- compare(interaction_clusters, phylo_clusters, method = "nmi")

# Step 4: Print results
print(paste("Adjusted Rand Index (ARI):", ari_value))
print(paste("Normalized Mutual Information (NMI):", nmi_value))

# Step 5: Visualize the matrices
heatmap(phylo_corr_matrix, main = "Phylogenetic Correlation Matrix", col = heat.colors(10))
heatmap(interaction_matrix, main = "Interaction Matrix", col = blues9)


######################


compute_optimal_cluster_metrics <- function(results_simulation, Smax, nbasals) {
  
  presence_matrix <- results_simulation$presence_matrix
  
  # Number of timesteps
  n_steps <- length(results_simulation$network_list)
  
  # Identify valid timesteps
  non.valid_timesteps_phylo_distance <- which(rowSums(presence_matrix) < 3)
  final.discarded_timestep <- ifelse(length(non.valid_timesteps_phylo_distance) > 0, 
                                     max(non.valid_timesteps_phylo_distance) + 1, 1)
  
  list_anc_dist <- results_simulation$list_anc_dist[(final.discarded_timestep + 1):n_steps]
  network_list <- results_simulation$network_list[(final.discarded_timestep + 1):n_steps]
  
  # Convert species names to letters
  list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
  list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
  list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
  
  colnames(presence_matrix) <- seq(1:Smax)
  colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
  
  presence_matrix <- presence_matrix[(final.discarded_timestep + 1):n_steps, ]
  
  # Store results
  ari_values <- numeric(length(list_networks_sppnames_letters))
  nmi_values <- numeric(length(list_networks_sppnames_letters))
  
  for (i in seq_along(list_anc_dist_letters)) {
    
    cat("Processing timestep:", i, "\n")
    
    # Process phylogenetic correlation matrix
    newick <- ToPhylo2(list_anc_dist_letters[[i]])
    newick_tail <- paste(newick, "root")
    tree <- read.tree(text = sub("A root",";", newick_tail))
    
    tree$edge.length <- sapply(tree$edge.length, function(x) ifelse(x == 0, 1e-5, x))
    phylo.vcv <- vcv(tree)
    phylo.corr <- cov2cor(phylo.vcv)
    
    alive_species <- names(which(presence_matrix[i, ] == 1))
    phylo_corr_cropped <- phylo.corr[alive_species, alive_species]
    
    # Clustering interaction matrix using SBM
    interaction_matrix <- list_networks_sppnames_letters[[i]][alive_species, alive_species]
    
    diag(interaction_matrix) <- 0  
    
    if (sum(interaction_matrix) == 0) {
      cat("Sparse or empty interaction matrix at timestep:", i, "- Assigning all species to one cluster.\n")
      interaction_clusters <- rep(1, nrow(interaction_matrix))
    } else {
      sbm_fit <- sbm::estimateSimpleSBM(interaction_matrix, model = "bernoulli")
      interaction_clusters <- sbm_fit$memberships
    }
    
    
    
    # Clustering phylogeny using k-means with automatic cluster number selection (gap statistic)
    phylo_dist <- as.dist(1 - phylo_corr_cropped)
    num_species <- nrow(phylo_corr_cropped)
    
    if (num_species > 2) {
      max_k <- min(4, num_species - 1)
      
      set.seed(123)
      gap_stat <- tryCatch({
        result <- clusGap(phylo_corr_cropped, FUN = kmeans, K.max = max_k, B = 50)
        if (is.null(result) || !is.list(result) || !("Tab" %in% names(result))) NULL else result
      }, error = function(e) {
        cat("Gap statistic failed at timestep:", i, "- Assigning all species to a single cluster.\n")
        NULL
      })
      
      if (!is.null(gap_stat) && is.matrix(gap_stat$Tab)) {
        optimal_k <- which.max(gap_stat$Tab[, "gap"])
        phylo_clusters <- kmeans(phylo_corr_cropped, centers = max(2, optimal_k), nstart = 25)$cluster
      } else {
        phylo_clusters <- rep(1, num_species)
      }
      
    } else {
      phylo_clusters <- rep(1, num_species)
    }
    
    # Compare clusters using ARI and NMI
    ari_values[i] <- cluster.stats(phylo_dist, phylo_clusters, interaction_clusters)$corrected.rand
    nmi_values[i] <- compare(interaction_clusters, phylo_clusters, method = "nmi")
  }
  
  
  return(data.frame(
    timestep = 1:length(ari_values),
    ARI = ari_values,
    NMI = nmi_values
  ))
}

# Run the function with appropriate inputs
cluster_results <- compute_optimal_cluster_metrics(results_simulation = res_sim, 
                                                   Smax = pars$Smax, 
                                                   nbasals = pars$Sbasal)


ggplot(cluster_results, aes(x = timestep)) +
  geom_point(aes(y = ARI), color = 'blue', alpha = 0.7, size = 3) +
  geom_line(aes(y = ARI), color = 'blue') +
  geom_point(aes(y = NMI), color = 'red', alpha = 0.7, size = 3) +
  geom_line(aes(y = NMI), color = 'red') +
  labs(x = "Time Step", y = "Clustering Metrics", 
       title = "ARI and NMI Over Time") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), name = "Metric Value") +
  scale_x_continuous(name = "Time Step")+
  theme(legend.position = "right")


