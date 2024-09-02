

timestep <- 98

tree_raw <- ToPhylo(list_anc_dist_letters[[timestep]])
tree_raw <- paste(tree_raw, "root")
tree_raw <- read.tree(text = sub("A root",";",tree_raw))

tree_corrected<-tree_raw
tree_corrected$edge.length<-sapply(tree_corrected$edge.length,function(x) ifelse(x==0,1e-5,x))
phylo.vcv<-vcv(tree_corrected)
phylo.corr<-cov2cor(phylo.vcv)

phylo.corr_cropped <- phylo.corr[names(which(presence_matrix[timestep,] == 1)), 
                                        names(which(presence_matrix[timestep,] == 1))]


list_phylo.corr_cropped[[timestep]]

#match_phylogeny_interaction <- match(rownames(data.inter), rownames(phylo.corr))

svd_eigen.phy<-eigen(list_phylo.corr_cropped[[timestep]], symmetric = T)$vec



# Perform SVD
matrix <- list_net_present_spp.letters[[timestep]]

svd_result <- svd(matrix)

# Extract left (U), middle (D), and right (V) matrices
U <- svd_result$u       # Left matrix (n_prey x kept axes)
D <- diag(svd_result$d) # Middle matrix (kept axes x kept axes)
V <- svd_result$v       # Right matrix (n_predators x kept axes)

# Decide on the number of axes to keep
kept_axes <- ncol(matrix)  # For example, keep the first 3 axes

# Select the kept axes from U, D, and V
U_kept <- U[, 1:kept_axes]   # n_prey x kept axes
D_kept <- D[1:kept_axes, 1:kept_axes] # kept axes x kept axes
V_kept <- V[, 1:kept_axes]   # n_predators x kept axes

# Transpose the right matrix V to match the desired output
V_kept_transposed <- t(V_kept)  # Transpose to get (kept axes x n_predators)


protest(svd_eigen.phy, V_kept_transposed)

diff_matrix <- svd_eigen.phy - V_kept_transposed
print(max(abs(diff_matrix)))  # Check the maximum absolute difference


#### try different dim

kept_axes <- 30
svd_eigen.phy_kept <- svd_eigen.phy[, 1:kept_axes]
V_kept_transposed_kept <- V_kept_transposed[, 1:kept_axes]

protest_result <- protest(svd_eigen.phy_kept, V_kept_transposed_kept)
print(protest_result)

#######################################################

list_network[[15]]

# Define a function to classify species roles and plot the network
plot_network_at_timestep <- function(matrix, timestep) {
  # Create an igraph object from the adjacency matrix
  g <- graph_from_adjacency_matrix(as.matrix(matrix), mode = "directed")
  
  # Identify species roles
  in_degree <- igraph::degree(g, mode = "in")   # Number of incoming edges (eaten by others)
  out_degree <- igraph::degree(g, mode = "out") # Number of outgoing edges (eat others)
  
  # Assign roles
  V(g)$role <- ifelse(in_degree > 0 & out_degree == 0, "Herbivore", 
                      ifelse(in_degree > 0 & out_degree > 0, "Intermediate", 
                             "Top Predator"))
  
  # Set colors for the roles
  V(g)$color <- ifelse(V(g)$role == "Herbivore", "green",
                       ifelse(V(g)$role == "Intermediate", "orange", "red"))
  
  # Plot the network
  plot <- ggplot(igraph::as_data_frame(g, what = "edges"), aes(x = from, y = to)) +
    geom_point(data = data.frame(x = V(g)$name, role = V(g)$role), aes(x = x, y = x, color = role)) +
    geom_segment(aes(x = from, xend = to, y = from, yend = to), arrow = arrow(length = unit(0.2, "cm"))) +
    scale_color_manual(values = c("green", "orange", "red")) +
    labs(title = paste("Network at Timestep", timestep), color = "Species Role") +
    theme_minimal()
  
  return(plot)
}

# Timesteps to plot
timesteps_to_plot <- c(5, 10, 20, 30, 50, 60, 80, 90, 98)

# List to store the plots
plot_list <- list()

# Loop through the timesteps and create the plots
for (timestep in timesteps_to_plot) {
  plot_list[[as.character(timestep)]] <- plot_network_at_timestep(list_network[[timestep]], timestep)
}

# Arrange the plots using ggarrange
final_plot <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)

# Display the final plot
print(final_plot)



#################################

# Define a function to classify species roles and plot the network using igraph
plot_network_at_timestep <- function(matrix, timestep) {
  # Create an igraph object from the adjacency matrix
  g <- graph_from_adjacency_matrix(as.matrix(matrix), mode = "directed")
  
  # Identify species roles
  in_degree <- igraph::degree(g, mode = "in")   # Number of incoming edges (eaten by others)
  out_degree <- igraph::degree(g, mode = "out") # Number of outgoing edges (eat others)
  
  # Assign roles
  V(g)$role <- ifelse(in_degree > 0 & out_degree == 0, "Herbivore", 
                      ifelse(in_degree > 0 & out_degree > 0, "Intermediate", 
                             "Top Predator"))
  
  # Set colors for the roles
  colFG <- c("green", "orange", "red") # Colors for Herbivore, Intermediate, Top Predator
  V(g)$color <- colFG[as.numeric(factor(V(g)$role, levels = c("Herbivore", "Intermediate", "Top Predator")))]
  
  # Plot the food web with custom colors
  plot(g, vertex.color = V(g)$color, 
       vertex.label = NA, # Remove labels for clarity
       edge.width = 0.3, edge.arrow.size = 0.3, 
       main = paste("Network at Timestep", timestep))
}

# Timesteps to plot
timesteps_to_plot <- c(5, 10, 20, 30, 50, 60, 80, 90, 98)

# Create plots and store them in a list
plot_list <- lapply(timesteps_to_plot, function(timestep) {
  plot_network_at_timestep(list_network[[timestep]], timestep)
})


# Use ggarrange to arrange the plots in a single panel
ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)


##############################

library(ggraph)

# Define a function to classify species roles and create a network plot using ggraph
plot_network_at_timestep_ggraph <- function(matrix, timestep) {
  # Create an igraph object from the adjacency matrix
  g <- graph_from_adjacency_matrix(as.matrix(matrix), mode = "directed")
  
  # Identify species roles
  in_degree <- igraph::degree(g, mode = "in")   # Number of incoming edges (eaten by others)
  out_degree <- igraph::degree(g, mode = "out") # Number of outgoing edges (eat others)
  
  # Assign roles
  V(g)$role <- ifelse(in_degree > 0 & out_degree == 0, "Herbivore", 
                      ifelse(in_degree > 0 & out_degree > 0, "Intermediate", 
                             "Top Predator"))
  
  # Create a ggraph plot
  plot <- ggraph(g, layout = "fr") + # Use 'fr' (Fruchterman-Reingold) layout for better spacing
    geom_edge_link(aes(edge_width = 0.08), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) + 
    geom_node_point(aes(color = role), size = 5) + # Plot nodes with color based on role
    scale_color_manual(values = c("Herbivore" = "green", "Intermediate" = "orange", "Top Predator" = "red")) +
    labs(title = paste("Network at Timestep", timestep), color = "Species Role") +
    theme_void() + # Use a minimal theme
    theme(legend.position = "bottom") # Place legend at the bottom
  
  return(plot)
}

# Timesteps to plot
timesteps_to_plot <- c(5, 10, 20, 30, 50, 60, 80, 90, 98)

# Create plots for the specified timesteps using ggraph
plot_list <- lapply(timesteps_to_plot, function(timestep) {
  plot_network_at_timestep_ggraph(list_network[[timestep]], timestep)
})

# Arrange the plots using ggarrange
final_plot <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)

# Display the final plot
print(final_plot)

################################

# Define a function to classify species roles and create a stratified network plot using ggraph
plot_network_at_timestep_stratified <- function(matrix, timestep) {
  # Create an igraph object from the adjacency matrix
  g <- graph_from_adjacency_matrix(as.matrix(matrix), mode = "directed")
  
  # Identify species roles
  in_degree <- igraph::degree(g, mode = "in")   # Number of incoming edges (eaten by others)
  out_degree <- igraph::degree(g, mode = "out") # Number of outgoing edges (eat others)
  
  # Assign roles
  V(g)$role <- ifelse(in_degree > 0 & out_degree == 0, "Herbivore", 
                      ifelse(in_degree > 0 & out_degree > 0, "Intermediate", 
                             "Top Predator"))
  
  # Define y-coordinates for stratification based on species role
  V(g)$y <- ifelse(V(g)$role == "Herbivore", 1, 
                   ifelse(V(g)$role == "Intermediate", 2, 3))
  
  # Generate a custom layout with stratified y-coordinates and evenly spaced x-coordinates
  layout_df <- data.frame(x = 1:vcount(g), y = V(g)$y)
  
  # Create a ggraph plot with a custom layout
  plot <- ggraph(g, layout = layout_df) + 
    geom_edge_link(edge_width = 0.01, arrow = arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) + 
    geom_node_point(aes(color = role), size = 3) + # Plot nodes with color based on role
    scale_color_manual(values = c("Herbivore" = "green", "Intermediate" = "orange", "Top Predator" = "red")) +
    labs(title = paste("Network at Timestep", timestep), color = "Species Role") +
    theme_void() + # Use a minimal theme
    theme(legend.position = "bottom") # Place legend at the bottom
  
  return(plot)
}

# Timesteps to plot
timesteps_to_plot <- c(5, 10, 20, 30, 50, 60, 80, 90, 98)

# Create plots for the specified timesteps using ggraph with stratified layout
plot_list <- lapply(timesteps_to_plot, function(timestep) {
  plot_network_at_timestep_stratified(list_network[[timestep]], timestep)
})

# Arrange the plots using ggarrange
final_plot <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)

# Display the final plot
print(final_plot)