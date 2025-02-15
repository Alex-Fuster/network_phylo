
library(ggplot2)
# Function to test the topology generated by the niche model
test_niche_topology <- function(pars, num_species) {
  
  # Generate random traits for the species
  traits_mat <- as.data.frame(matrix(NA, nrow = num_species, ncol = 3))
  colnames(traits_mat) <- c("n", "r", "o")
  for (i in 1:num_species) {
    traits_mat[i, ] <- rand_traits_anc(pars)  # Generate random traits
  }
  
  # Generate basal species traits
  basal <- runif(pars$Sbasal, 0, 0.005)
  
  # Compute the interaction network
  L <- get_L_mat(basal, pars, traits_mat)
  
  # Identify trophic levels
  herbivores <- which(colSums(L[1:pars$Sbasal, ]) > 0 & colSums(L[(pars$Sbasal + 1):nrow(L), ]) == 0)
  predators <- which(colSums(L[(pars$Sbasal + 1):nrow(L), ]) > 0 & colSums(L[1:pars$Sbasal, ]) == 0)
  top_predators <- which(colSums(L[(pars$Sbasal + 1):nrow(L), ]) > 0 & rowSums(L) == 0)
  
  # Return the number of top predators
  length(top_predators)
}

# Parameters
set.seed(42)
av_r_values <- seq(0.05, 0.4, length.out = 8)  # Vary av_r values from 0.05 to 0.4
num_species <- 60
num_simulations <- 100

# Create a data frame to store results for all simulations
results_df <- data.frame(av_r = numeric(), top_predators = numeric())

# Loop over av_r values and simulate
for (av_r in av_r_values) {
  pars <- list(Sbasal = 5, Smax = 60, av_r = av_r)
  
  # Run 100 simulations for each av_r value and store number of top predators
  predator_counts <- sapply(1:num_simulations, function(x) test_niche_topology(pars, num_species))
  
  # Combine results into a data frame
  results_df <- rbind(results_df, data.frame(av_r = av_r, top_predators = predator_counts))
}

# Plot the results using boxplot
ggplot(results_df, aes(x = factor(av_r), y = top_predators)) +
  geom_boxplot() +
  xlab("av_r (Niche Range)") +
  ylab("Number of Top Predators") +
  ggtitle("Distribution of Top Predators vs av_r") +
  theme_minimal()



###################

pars <- list(Sbasal = 5, Smax = 60, av_r = 0.1)
test_result <- test_niche_topology(pars, num_species = 60)
print(test_result)


####################


library(ggplot2)
library(dplyr)

# Function to generate the predator-prey plot from a specific timestep's traits data
plot_predator_prey_relationship <- function(traits_df, basal_species_count) {
  # Extract body size, range, and optimum for each species from traits_df
  traits_mat <- traits_df %>% 
    rename(n = n, r = r, o = o) %>%
    mutate(r = 0.12,            # Fixed range
           o = n / 2)           # Optimum as half of niche position
  
  # Generate basal species (prey only)
  basal <- traits_mat[1:basal_species_count, ]
  traits_mat <- rbind(basal, traits_mat[-(1:basal_species_count), ])
  
  # Convert basal niche positions to a vector
  basal_vec <- as.numeric(basal$n)
  
  # Compute the interaction matrix using get_L_mat
  L <- get_L_mat(basal_vec, pars, traits_mat)
  
  # Extract predator-prey pairs based on L matrix
  interactions <- expand.grid(predator = 1:nrow(L), prey = 1:ncol(L)) %>%
    filter(L[cbind(predator, prey)] == 1) %>%
    mutate(
      predator_size = traits_mat$n[predator],
      prey_size = traits_mat$n[prey],
      predator_range = traits_mat$r[predator],
      predator_optimum = traits_mat$o[predator]
    )
  
  # Generate 1:1 line for reference
  one_to_one <- data.frame(x = seq(0, 1, length.out = 100), y = seq(0, 1, length.out = 100))
  
  # Plot the predator-prey body size relationship
  ggplot(interactions, aes(x = predator_size, y = prey_size)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # 1:1 line
    geom_point(aes(color = "Predator-Prey Pair"), alpha = 0.6) +
    geom_errorbar(aes(ymin = predator_optimum - predator_range / 2, 
                      ymax = predator_optimum + predator_range / 2), width = 0.02, color = "blue") +
    geom_smooth(method = "lm", color = "black", se = FALSE) +  # Linear fit (centroid)
    stat_quantile(quantiles = c(0.05, 0.95), linetype = "dotted", color = "red") +  # Range boundaries
    labs(x = "Predator Body Size", y = "Prey Body Size",
         title = "Predator-Prey Body Size Relationship at Given Timestep") +
    scale_color_manual(values = c("Predator-Prey Pair" = "darkblue")) +
    theme_minimal() +
    theme(legend.position = "none")
}


plot_predator_prey_relationship(res_sim$traits_df, basal_species_count = pars$Sbasal)
