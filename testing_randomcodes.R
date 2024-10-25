
t <- 5

present_species_t <- list_anc_dist_letters[[t]][which(list_anc_dist_letters[[t]]$"A/E" == "A"), "spp"]

list_networks_sppnames_letters[[t]][present_species_t, present_species_t]

which(presence_matrix[t,]==1)


###################

newick <- ToPhylo(list_anc_dist_letters[[t]])
newick_tail <- paste(newick, "root")
tree <- read.tree(text = sub("A root",";",newick_tail))
#list_dist.phylo[[i]] <- cophenetic.phylo(tree)

tree$edge.length<-sapply(tree$edge.length,function(x) ifelse(x==0,1e-5,x))
phylo.vcv<-vcv(tree)
phylo.corr<-cov2cor(phylo.vcv)
phylo.corr[names(which(presence_matrix[i,] == 1)), 
                                           names(which(presence_matrix[t,] == 1))]


####################

phy<-ToPhylo(list_anc_dist_letters[[5]])
tree<-read.tree(text = gsub("A",";",phy))
plot(tree)





####################

ToPhylo <- function(data) {
  # Make a copy of the data
  data.2 <- data
  data.2$repr <- data$spp
  
  # Identify unique species (tips) and ancestors (nodes)
  sisters <- unique(data$spp)
  mothers <- unique(data$ancestor)
  tips <- setdiff(sisters, mothers)
  root <- setdiff(mothers, sisters)
  
  # Initialize the nodes to focus on
  foc.nodes <- unique(data[which(data$spp %in% tips), "ancestor"])
  n <- length(foc.nodes)
  
  # Correct assignment of tips
  data.2$repr[data.2$spp %in% tips] <- as.character(data.2$spp[data.2$spp %in% tips])
  
  # Build the tree iteratively
  while (n > 0) {
    # Find next set of focus nodes
    foc.nodes2 <- unique(data.2[which(data.2$spp %in% foc.nodes), "ancestor"])
    
    # Process each focal node
    for (i in 1:n) {
      daughters <- data.2[which(data.2$ancestor == foc.nodes[i]), "repr"]
      daughters.dist <- data.2[which(data.2$ancestor == foc.nodes[i]), "distance"]
      
      # Combine daughter species into a Newick format string
      data.2$repr[data.2$spp == foc.nodes[i]] <- paste0(sister.group(daughters, daughters.dist), foc.nodes[i])
    }
    
    # Update the tips and focal nodes for the next iteration
    tips <- foc.nodes
    foc.nodes <- foc.nodes2
    n <- length(foc.nodes)
  }
  
  # Construct the final Newick format string for the tree
  daughters <- data.2[which(data.2$ancestor == root), "repr"]
  daughters.dist <- data.2[which(data.2$ancestor == root), "distance"]
  newick_string <- paste0(sister.group(daughters, daughters.dist), root)
  
  return(newick_string)
}


######################

ToPhylo <- function(data) {
  # Copy data and prepare for processing
  data.2 <- data
  data.2$repr <- data$spp
  
  # Identify alive species (tips) and ancestors (nodes)
  alive_species <- data[data$"A/E" == "A", "spp"]  # Only keep species marked as "A" (alive)
  sisters <- unique(data$spp)
  mothers <- unique(data$ancestor)
  
  # Identify tips (alive species not listed as an ancestor)
  tips <- setdiff(alive_species, mothers)
  
  # Identify the root (species that are ancestors but not in species)
  root <- setdiff(mothers, sisters)
  
  # Initialize the nodes to focus on
  foc.nodes <- unique(data.2[data.2$spp %in% tips, "ancestor"])
  n <- length(foc.nodes)
  
  # Correct assignment of tips
  data.2$repr[data.2$spp %in% tips] <- as.character(data.2$spp[data.2$spp %in% tips])
  
  # Build the tree iteratively
  while (n > 0) {
    # Find next set of focus nodes
    foc.nodes2 <- unique(data.2[which(data.2$spp %in% foc.nodes), "ancestor"])
    
    # Process each focal node
    for (i in 1:n) {
      daughters <- data.2[which(data.2$ancestor == foc.nodes[i]), "repr"]
      daughters.dist <- data.2[which(data.2$ancestor == foc.nodes[i]), "distance"]
      
      # Combine daughter species into a Newick format string
      data.2$repr[data.2$spp == foc.nodes[i]] <- paste0(sister.group(daughters, daughters.dist), foc.nodes[i])
    }
    
    # Update the tips and focal nodes for the next iteration
    tips <- foc.nodes
    foc.nodes <- foc.nodes2
    n <- length(foc.nodes)
  }
  
  # Construct the final Newick format string for the tree
  daughters <- data.2[which(data.2$ancestor == root), "repr"]
  daughters.dist <- data.2[which(data.2$ancestor == root), "distance"]
  newick_string <- paste0(sister.group(daughters, daughters.dist), root)
  
  return(newick_string)
}



phy<-ToPhylo(list_anc_dist_letters[[10]])
phy <- gsub("A", ";", phy) 
tree <- read.tree(text = phy)
plot(tree)




##########################


# Required Libraries
library(ggplot2)

# Parameters
set.seed(123)  # Set seed for reproducibility
n_samples <- 100  # Number of samples to generate
beta_n <- 2  # Shape parameter for beta distribution
pars <- list(beta_n = beta_n, sd = 0.1, r = 0.5)  # Example parameters

# Function to generate mutant traits
rand_traits_mut <- function(traits_anc, pars) {
  with(as.list(c(traits_anc, pars)), {
    
    # Compute the alpha parameter for the beta distribution
    a <- beta_n * n / (1 - n)
    
    # Determine the direction of divergence (either side of the ancestor trait)
    if (runif(1) < 0.5) {  # 50% chance for each direction
      # Generate mutant trait greater than ancestor trait
      n_m <- rbeta(1, shape1 = a, shape2 = beta_n)
      n_m <- n + abs(n_m - n)  # Ensure n_m is on the upper side of n
      n_m <- min(n_m, 1)  # Make sure n_m does not exceed 1
    } else {
      # Generate mutant trait less than ancestor trait
      n_m <- rbeta(1, shape1 = a, shape2 = beta_n)
      n_m <- n - abs(n_m - n)  # Ensure n_m is on the lower side of n
      n_m <- max(n_m, 0)  # Make sure n_m is not below 0
    }
    
    # Calculate the optimum (o_m) for the mutant
    o_m <- n_m / 2
    
    # Combine the new traits
    traits_mut <- c(n = n_m, r = r, o = o_m)
    
    traits_mut
  })
}

# Generate Ancestor Traits and Mutant Traits
ancestor_traits <- runif(n_samples, 0.1, 0.9)  # Random ancestor traits between 0.1 and 0.9
mutant_traits <- sapply(ancestor_traits, function(n) rand_traits_mut(list(n = n), pars)[1])

# Create DataFrame for Plotting
df <- data.frame(
  Ancestor = ancestor_traits,
  Mutant = mutant_traits
)

# Plot Ancestor vs. Mutant Traits
ggplot(df, aes(x = Ancestor, y = Mutant)) +
  geom_point(color = 'blue', alpha = 0.6) +  # Scatter plot of ancestor vs. mutant traits
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Reference line y = x
  labs(
    title = "Divergent Selection of Mutant Traits",
    x = "Ancestor Trait (n)",
    y = "Mutant Trait (n_m)"
  ) +
  theme_minimal() +
  coord_fixed(ratio = 1)  # Equal scaling for x and y axes


##################

singular_values <- c(3.2, 2.8, 1.5, 0.9, 0.5, 0.3, 0.2)
eigenvalues_phy <- c(4.5, 3.2, 2.5, 1.8, 1.2, 0.7, 0.4)


cum_var_pred <- cumsum(singular_values^2) / sum(singular_values^2)
cum_var_phy <- cumsum(eigenvalues_phy) / sum(eigenvalues_phy)

# Create a dataframe for plotting
df_variance <- data.frame(
  Dimensions = 1:length(cum_var_pred),
  CumulativeVarianceInteraction = cum_var_pred,
  CumulativeVariancePhylogeny = cum_var_phy
)


ggplot(df_variance, aes(x = Dimensions)) +
  geom_line(aes(y = CumulativeVarianceInteraction, color = "Interaction Matrix")) +
  geom_line(aes(y = CumulativeVariancePhylogeny, color = "Phylogenetic Matrix")) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
  labs(x = "Number of Dimensions", y = "Cumulative Explained Variance",
       title = "Cumulative Explained Variance by Number of Dimensions") +
  scale_color_manual(values = c("Interaction Matrix" = "blue", "Phylogenetic Matrix" = "green")) +
  theme_minimal() +
  theme(legend.title = element_blank())



####################

# TESTING METHODS TO ANALYZE THE CORRESPONDENCE OF THE NETWORK WITH THE PHYLOGENY

* adonis from vegan (or a similar function) can make an analysis similar to a distance-based manova, i.e. trying to assess how much a particular grouping of species explains a pattern of distances. Trying to explain Jaccard distances based on interactions using groups obtained from the phylogeny might be something to try. Or, because adonis is flexible, I think we might also want to explain the distance matrix of interaction by the distance matrix obtained from the phylogeny

1) Obtain distance matrix by interactions (by NMI and Jaccard)
2) groups obtained from the phylogeny (?)
3) distance matrix from phyogeny


```{r}
library(vegan)
```

1) Obtain distance matrix by interactions (by NMI and Jaccard)



###############################

# Define the probability of interaction function (example: Gaussian function)
prob_interaction <- function(distance, sigma = 0.04) {
  exp(- (distance^2) / (2 * sigma^2))
}

# Generate a range of distances
distance_values <- seq(0, 1, length.out = 100)

# Compute interaction probabilities for each distance
sigma_values <- c(0.01, 0.02, 0.04, 0.1)  # Example sigma values for comparison
probabilities <- sapply(sigma_values, function(sigma) sapply(distance_values, prob_interaction, sigma = sigma))

# Convert to data frame for plotting
prob_df <- data.frame(distance = rep(distance_values, times = length(sigma_values)),
                      interaction_prob = as.vector(probabilities),
                      sigma = factor(rep(sigma_values, each = length(distance_values))))

# Plot using ggplot2
library(ggplot2)
ggplot(prob_df, aes(x = distance, y = interaction_prob, color = sigma)) +
  geom_line() +
  labs(x = "Distance Between Predator and Prey Optima",
       y = "Probability of Interaction",
       color = "Sigma") +
  ggtitle("Probability of Interaction vs. Distance") +
  theme_minimal()





# Function to compute the full interaction matrix for all species
get_L_mat <- function(basal, pars, traits_mat) {
  with(as.list(pars), {
    # Initialize the interaction matrix
    L <- matrix(0, nr = Smax + Sbasal, nc = Smax)
    
    # Lower and upper boundaries for niches based on predator's optimum and range
    low <- traits_mat$o - traits_mat$r
    low_mat <- matrix(low, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)
    high <- traits_mat$o + traits_mat$r
    high_mat <- matrix(high, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)
    
    # Matrix of niche positions for all species
    n_mat <- matrix(traits_mat$n, nr = Smax, nc = Smax, byrow = FALSE)
    
    # Add basal species to the prey matrix
    n_basal <- matrix(basal, nr = Sbasal, nc = Smax, byrow = FALSE)
    n_mat <- rbind(n_basal, n_mat)
    
    # Define the probability function (e.g., Gaussian probability)
    prob_interaction <- function(distance, sigma = 0.04) {
      exp(- (distance^2) / (2 * sigma^2))
    }
    
    # Calculate distances and interaction probabilities
    for (i in 1:(Smax + Sbasal)) {
      for (j in 1:Smax) {
        # Check for NA values before comparing
        if (!is.na(n_mat[i, j]) && !is.na(low_mat[i, j]) && !is.na(high_mat[i, j])) {
          # First check if prey falls within the predator's range
          if (n_mat[i, j] > low_mat[i, j] && n_mat[i, j] < high_mat[i, j]) {
            
            # Compute the distance from the predator's optimal niche to prey's niche
            distance <- abs(n_mat[i, j] - traits_mat$o[j])
            
            # Compute the initial probability of interaction based on distance
            interaction_prob <- prob_interaction(distance)
            
            # Apply a penalty if prey's niche value is higher than the predator's
            if (n_mat[i, j] > traits_mat$n[j]) {
              interaction_prob <- interaction_prob * 0.2
            }
            
            # Assign interaction based on probability
            L[i, j] <- rbinom(1, 1, interaction_prob)  # Binomial draw: 1 interaction with probability 'interaction_prob'
          }
        }
      }
    }
    
    # Set diagonal to 0 (no self-interaction)
    if (Smax > 1) diag(L[(Sbasal + 1):(Sbasal + Smax), ]) <- 0
    L
  })
}