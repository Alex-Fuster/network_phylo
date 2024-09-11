
########################################
# Function to generate a new set of traits for ancestors
# Each species is characterized by a set of 3 traits: n, o and r
rand_traits_anc = function(pars) {
  with(as.list(pars),{
    n = runif(1, 0, 1)
    r = av_r
    o = runif(1,0,n)
    traits = c(n = n, r = r, o = o)
    traits
  })
}

########################################
rand_traits_mut <- function(traits_anc, pars, direction = "random") {
  with(as.list(c(traits_anc, pars)), {
    
    # Compute the alpha parameter for the beta distribution
    a <- beta_n * n / (1 - n)
    
    # Determine the direction of divergence
    if (direction == "random") {
      direction <- ifelse(runif(1) < 0.5, "greater", "lesser")
    }
    
    if (direction == "greater") {
      # Generate mutant trait greater than ancestor trait
      n_m <- rbeta(1, shape1 = a, shape2 = beta_n)
      n_m <- n + abs(n_m - n)  # Ensure n_m is on the upper side of n
      n_m <- min(n_m, 1)  # Make sure n_m does not exceed 1
    } else if (direction == "lesser") {
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

########################################
# Function to compute the interaction network from a set of traits
# get_L_mat = function(basal, pars, traits_mat) {
#   with(as.list(pars),{
#     L = matrix(0, nr = Smax+Sbasal, nc = Smax)
#     
#     # Lower boundary
#     low = traits_mat$o - traits_mat$r
#     low_mat = matrix(low, nr = Smax+Sbasal, nc = Smax, byrow = TRUE)
#     
#     # Upper boundary
#     high = traits_mat$o + traits_mat$r
#     high_mat = matrix(high, nr = Smax+Sbasal, nc = Smax, byrow = TRUE)	
#     S = nrow(traits_mat)
#     
#     # Matrix of niche positions
#     n_mat = matrix(traits_mat$n, nr = Smax, nc = Smax, byrow = FALSE)
#     
#     # Add the basal species
#     n_basal = matrix(basal, nr = Sbasal, nc = Smax, byrow = FALSE)
#     n_mat = rbind(n_basal, n_mat)
#     
#     # Test interactions
#     L[n_mat > low_mat & n_mat < high_mat] = 1
#     if(Smax > 1) diag(L[(Sbasal+1):(Sbasal+Smax),]) = 0
#     L
#   })
# }

########################################
# get_L_mat <- function(basal, pars, traits_mat) {
#   with(as.list(pars), {
#     # Initialize the interaction matrix
#     L <- matrix(0, nr = Smax + Sbasal, nc = Smax)
#     
#     # Lower and upper boundaries for niches
#     low <- traits_mat$o - traits_mat$r
#     low_mat <- matrix(low, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)
#     high <- traits_mat$o + traits_mat$r
#     high_mat <- matrix(high, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)  
#     S <- nrow(traits_mat)
#     
#     # Matrix of niche positions
#     n_mat <- matrix(traits_mat$n, nr = Smax, nc = Smax, byrow = FALSE)
#     
#     # Add the basal species
#     n_basal <- matrix(basal, nr = Sbasal, nc = Smax, byrow = FALSE)
#     n_mat <- rbind(n_basal, n_mat)
#     
#     # Define the probability function (e.g., Gaussian probability)
#     prob_interaction <- function(distance, sigma = 0.1) {
#       exp(- (distance^2) / (2 * sigma^2))
#     }
#     
#     # Calculate distances and interaction probabilities
#     for (i in 1:(Smax + Sbasal)) {
#       for (j in 1:Smax) {
#         if (n_mat[i, j] > low_mat[i, j] && n_mat[i, j] < high_mat[i, j]) {
#           # Compute the distance from the optimal niche
#           distance <- abs(n_mat[i, j] - traits_mat$o[j])
#           # Compute probability of interaction
#           interaction_prob <- prob_interaction(distance)
#           # Assign interaction based on probability
#           L[i, j] <- rbinom(1, 1, interaction_prob)  # Binomial draw: 1 interaction with probability 'interaction_prob'
#         }
#       }
#     }
#     
#     # Set diagonal to 0 (no self-interaction)
#     if (Smax > 1) diag(L[(Sbasal + 1):(Sbasal + Smax), ]) <- 0
#     L
#   })
# }


get_L_mat <- function(basal, pars, traits_mat) {
  with(as.list(pars), {
    # Initialize the interaction matrix
    L <- matrix(0, nr = Smax + Sbasal, nc = Smax)
    
    # Lower and upper boundaries for niches
    low <- traits_mat$o - traits_mat$r
    low_mat <- matrix(low, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)
    high <- traits_mat$o + traits_mat$r
    high_mat <- matrix(high, nr = Smax + Sbasal, nc = Smax, byrow = TRUE)  
    S <- nrow(traits_mat)
    
    # Matrix of niche positions
    n_mat <- matrix(traits_mat$n, nr = Smax, nc = Smax, byrow = FALSE)
    
    # Add the basal species
    n_basal <- matrix(basal, nr = Sbasal, nc = Smax, byrow = FALSE)
    n_mat <- rbind(n_basal, n_mat)
    
    # Define the probability function (e.g., Gaussian probability)
    prob_interaction <- function(distance, sigma = 0.1) {
      exp(- (distance^2) / (2 * sigma^2))
    }
    
    # Calculate distances and interaction probabilities
    for (i in 1:(Smax + Sbasal)) {
      for (j in 1:Smax) {
        # Check for NA values before comparing
        if (!is.na(n_mat[i, j]) && !is.na(low_mat[i, j]) && !is.na(high_mat[i, j])) {
          if (n_mat[i, j] > low_mat[i, j] && n_mat[i, j] < high_mat[i, j]) {
            # Compute the distance from the optimal niche
            distance <- abs(n_mat[i, j] - traits_mat$o[j])
            # Compute probability of interaction
            interaction_prob <- prob_interaction(distance)
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



########################################
# Function to compute the interactions of a given species
get_L_vec = function(basal, pars, traits_mat, traits_mut) {
  with(as.list(pars),{
    L_vec = numeric(Smax+Sbasal)
    
    # Lower boundary
    low = traits_mut["o"] - traits_mut["r"]
    
    # Upper boundary
    high = traits_mut["o"] + traits_mut["r"]
    
    # Vector of niche positions
    n_vec = c(basal, traits_mat$n)
    
    # Test interactions
    L_vec[n_vec > as.numeric(low) & n_vec < as.numeric(high)] = 1
    L_vec
  })
}



########################################

sim_model = function(seed, pars, nsteps) {
  
  with(pars, {
    
    set.seed(seed)
    
    # Draw the traits of the producers
    basal = runif(pars$Sbasal, 0, 0.2)
    
    # Draw the first species traits
    traits_mat = matrix(nr = Smax, nc = 3)
    traits_mat[1,] = rand_traits_anc(pars)
    traits_mat = as.data.frame(traits_mat)
    names(traits_mat) = c("n", "r", "o")
    
    # Set the presence/absence matrix
    pres = matrix(0, nr = nsteps, nc = Smax)
    pres[1, 1] = 1
    
    # Set the ancestry object
    anc = matrix(0, nr = Smax, nc = 3)
    
    # Set the extinctions object
    extinct = matrix(NA,nr = Smax, nc = 2)
    
    #  [TABLE DIST_ANC] Set the distance ancestry matrix -----------------------------------
    # (new_spp, ancestor, tip, distance)
    dist_anc = as.data.frame(matrix(NA, nr = Smax, nc = 4))
    colnames(dist_anc) = c("spp", "ancestor", "A/E", "distance")
    dist_anc[1,] <- c(1, 0, "A", 0)
    dist_anc$distance <- as.numeric(dist_anc$distance)
    dist_anc$distance[is.na(dist_anc$distance)] <- 0
    
    #  [TABLE DIST_ANC] record the dist_anc at each timestep
    
    list_dist_anc <- list()
    
    # Record the matrices 
    L_list = list()
    L_cropped_list = list()
    
    
    # Species count
    S = 1
    
    
    ##########################
    # MAIN LOOP
    for(step in 2:nsteps) {
      ActualS = sum(pres[step - 1,])
      cat("Step = ", step - 1, " Total S = ", S, " Actual S = ", ActualS, '\n')
      
      
      
      
      # Test for successful speciation probability
      for(i in 1:Smax) {
        
        if(S >= Smax) break
        
        # Speciation occurs if the species is present
        if(pres[step - 1, i] == 1) {
          
          
          # [TABLE DIST_ANC] Define distances of Alive spp -------------------
          
          dist_anc[i,"A/E"] <- "A"
          dist_anc[i,"distance"] <- as.numeric(dist_anc[i,"distance"])+1
          
          
          # Species is maintained
          pres[step, i] = 1
          
          # Test if there is mutation
          test_number = runif(1, 0, 1)
          speciation_prob = u_max / (1 + exp(d * (ActualS - I_max)))
          
          if(test_number < speciation_prob) {
            
            # Ancestor traits
            traits_anc <- traits_mat[i,]
            
            # Pick new parameters
            # traits_mut = rand_traits_mut(traits_mat[i, ], pars) 
            
            
            # ------------- 1st mutant
            
            # Generate the first mutant with a random direction
            first_mutant_direction <- ifelse(runif(1) < 0.5, "greater", "lesser")
            traits_mut1 <- rand_traits_mut(traits_anc, pars, direction = first_mutant_direction)
            
            
            # Recompute the interactions 
            I = get_L_vec(basal, pars, traits_mat, traits_mut1)
            
            # Compute the number of interactions among present species
            sum_I = sum(I * c(rep(1, Sbasal), pres[step, ]))	
            
            # # Compute the probability of establishment	
            # if(int == 0) {
            #   estab_prob_sel = u_0neg + u_1neg * exp(-a_uneg * sum_I)
            #   estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)
            # }
            # 
            # if(int == 1 | int == 2) {
            #   estab_prob_sel = u_0pos + u_1pos * exp(-a_upos * sum_I)
            #   estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)
            # }
            
            # Test if there is speciation
            if(runif(1) < 1) {
              S = S + 1
              traits_mat[S, ] = traits_mut1 
              pres[step, S] = 1
              anc[S,] = c(step, i, S) #Time step, ancestor, new_species)
              
              # [TABLE DIST_ANC] add new spp to table dist_anc ------------
              
              dist_anc[S, ] <- c(S, i, "A", 1)
              
              if(S >= Smax) break
              
              pres[step,i] = 0 #the ancestor disapear because of trait deplacement
              
            }
            
            
            
            
            # ------------- 2nd mutant
            
            # Determine the opposite direction for the second mutant
            second_mutant_direction <- ifelse(first_mutant_direction == "greater", "lesser", "greater")
            traits_mut2 <- rand_traits_mut(traits_anc, pars, direction = second_mutant_direction)
            
            # Recompute the interactions 
            I = get_L_vec(basal, pars, traits_mat, traits_mut2)
            
            # Compute the number of interactions among present species
            sum_I = sum(I * c(rep(1, Sbasal), pres[step, ]))	
            
            # # Compute the probability of establishment	
            # if(int == 0) {
            #   estab_prob_sel = u_0neg + u_1neg * exp(-a_uneg * sum_I)
            #   estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)
            # }
            # 
            # if(int == 1 | int == 2) {
            #   estab_prob_sel = u_0pos + u_1pos * exp(-a_upos * sum_I)
            #   estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)
            # }
            
            # Test if there is speciation
            if(runif(1) < 1) {
              
              print("2nd spp established")
              
              S = S + 1
              traits_mat[S, ] = traits_mut2
              pres[step, S] = 1
              anc[S,] = c(step, i, S) #Time step, ancestor, new_species)
              
              # [TABLE DIST_ANC] add new spp to table dist_anc ------------
              
              dist_anc[S, ] <- c(S, i, "A", 1)
              
              if(S >= Smax) break
              
              pres[step,i] = 0 #the ancestor disapear because of trait deplacement
              
            }
            
            
          }
        }
      }
      
      if(S >= Smax) break
      # Test for extinction
      
      # Compute the interaction matrix among present species
      pres_vec = pres[step, ]
      cooc = matrix(pres_vec, nr = Smax, nc = Smax, byrow = TRUE) * 
        matrix(pres_vec, nr = Smax, nc = Smax, byrow = FALSE)
      L = get_L_mat(basal, pars, traits_mat)
      L[c((Sbasal + 1):(Sbasal + Smax)), ] = L[c((Sbasal + 1):(Sbasal + Smax)), ] * cooc
      L_list[[step]] = L
      
      
      
      
      # Identify present species using the presence matrix
      present_species <- which(pres[step, ] == 1)  # Indices of species that are present
      
      if(length(present_species) > 1){
        
        # Crop the interaction matrix L to only include present species
       # L_cropped <- L[present_species, present_species]
        L_cropped <- L[c(1:Sbasal, Sbasal + present_species), present_species]
        
        L_cropped_list[[step]] <- L_cropped
        
        #  # Perform SVD on the cropped matrix
          svd_result <- svd(L_cropped)
          
         # svd_result$u %*% diag(svd_result$d) %*% t(V)
          
        #  # Compute the similarity matrix from the right singular vectors
          similarity_matrix <- svd_result$v %*% diag(svd_result$d) %*% t(svd_result$v)
        #  
        # set the diagonal to 0
          diag(similarity_matrix) <- 0
        #  # Calculate average similarity for each species
          avg_similarity <- rowMeans(similarity_matrix)*(ncol(similarity_matrix)/(ncol(similarity_matrix)-1))
          
          # Identify species with no interactions (rows with all zeros)
          no_interaction_species <- rowSums(L_cropped) == 0
          
          # Set average similarity values smaller than e^-5 to 0
          avg_similarity[avg_similarity < exp(-5)] <- 0
          
        #  
        # Calculate extinction probabilities based on avg_similarity
         # ext_prob_sel <- e_0neg + e_1neg * exp(-a_eneg * avg_similarity)
          
          
          ############## compute extinction probability
           in_I = colSums(L_cropped)
           out_I = rowSums(L_cropped)
           ext_prob_sel = e_0neg + e_1neg * exp(-a_eneg * out_I) + 
             e_0pos + e_1pos * exp(-a_epos * in_I)
           
         #  print(paste("ext_prob_sel formula", ext_prob_sel))
           
           # Adjust extinction probability with competition
           competition_factor = 1 + competition_coefficient * avg_similarity
           
         #  print(paste("competition_factor", competition_factor))
           
           ext_prob_sel = ext_prob_sel * competition_factor
        #   print(paste("in_I", in_I))
        #   print(paste("out_I", out_I))
           
        #  print(paste("ext_prob_sel", ext_prob_sel))
        #  
         #  print(L_cropped)
         # print(avg_similarity)
         #  print(ext_prob_sel)
        #  
          # Initialize extinction probabilities vector for all species
          ext_prob_sel_full <- numeric(Smax)
        #  
          # Assign calculated probabilities to the correct positions for present species
          ext_prob_sel_full[present_species] <- ext_prob_sel
          
          # Identify herbivores and predators
          # Herbivores: Species consuming only basal species
          herbivores <- rowSums(L[(Sbasal + 1):(Sbasal + Smax), 1:Sbasal]) > 0 & 
            rowSums(L[(Sbasal + 1):(Sbasal + Smax), (Sbasal + 1):Smax]) == 0
          
          # Predators: Species consuming other non-basal species
          predators <- rowSums(L[(Sbasal + 1):(Sbasal + Smax), (Sbasal + 1):Smax]) > 0
          
          # Apply the extinction rule for predators with no prey
          predators_with_no_prey <- predators & (rowSums(L[(Sbasal + 1):(Sbasal + Smax), (Sbasal + 1):Smax]) == 0)
          
          # Set extinction probability to 1 for predators with no prey
          ext_prob_sel_full[predators_with_no_prey] <- 1
          
        # Neutral + selection
          
          ext_prob = SN * (ext_prob_neutral) + (1 - SN) * (ext_prob_sel_full)
        

        
        # Perform extinctions
        
        
        present_spe <- grep(1,pres[step,])
        test_extprob <- rep(0,Smax)
        random_number <- runif(length(present_spe),0,1)
        test_extprob[present_spe] <- random_number
        
        
        
        pres[step, pres[step-1,] & test_extprob < ext_prob] = 0
        
      }
      
      for(i in 1:Smax) {
        if(pres[step,i] != pres[step-1,i] & pres[step,i] == 0){
          
          # [TABLE DIST_ANC] Define a new E spp -------------------
          
          dist_anc[i,"A/E"] <- "E"
          
          
          
        }
      }
      
      # Validation step (optional)
      dist_anc1 <- dist_anc
      colnames(dist_anc1)[colnames(dist_anc1) == "A/E"] <- "Status"
      extinct_in_pres <- which(pres[step,] == 0 & pres[step-1,] == 1)
      extinct_in_dist_anc <- which(dist_anc1$Status == "E")
      if(!all(extinct_in_pres %in% extinct_in_dist_anc)) {
        warning("Mismatch between pres and dist_anc: Check species extinction assignments.")
      }
      
      
      # [TABLE DIST_ANC] save table dist_anc at timestep step
      
      list_dist_anc[[step]] <- na.omit(dist_anc)
      
      
      if(step == 30 & sum(pres[step,]) <= 3){
        
        cat("Simulation stopped due to low species presence.\n")
        break
      }
      
      
    } # End of main loop
    
    list(pres = pres, 
         traits = traits_mat, 
         anc = anc, 
         L_list = L_list, 
         L_cropped_list = L_cropped_list,
         basal = basal,
         dist_anc = dist_anc,
         list_dist_anc = list_dist_anc)
    
  })
}





# competition_coefficient = 0.025
# Sbasal = 25 # number of basal species
# Smax = 1000 # Maximal number of species allowed
# #Bspe = 4 # Basal species impact on the speciation probality
# #Bext = 4 # Basal species impact on extinction probability
# av_r = 0.185 # 0.2 range of the niche
# sd = 0.5*av_r + 0.0001 # Standard deviation of the normal distribution used to calculate the niche optimum trait of a new species
# #sd = 0.5*av_r + 0.0001 # Standard deviation of the normal distribution used to calculate the niche optimum trait of a new species
# 
# # PROBABILITY OF MUTATION
# 
# u_max = 0.23#0.15 #0.15 # mutation probability (0.075, )
# d = 0.5 # Decrease speed of the establishment probability
# I_max = 60 # Maximal number of interactioning species
# beta_n =  1 # parameter of the beta distribution 
# 
# # STRENGTH OF NEUTRAL-DRIVEN EVOLUTION
# 
# SN = 0 # strength for neutral-driven evolution
# 
# # PROBABILITY OF ESTABLISHMENT
# 
# # strength of selection-driven selection is 1 - SN
# estab_prob_neutral = rep(0, Smax) # neutral probability of establishment
# 
# # Facilitation & Foodweb
# 
# u_0pos = 1  
# u_1pos = -1 
# a_upos = 0.45 
# 
# # PROBABILITY OF EXTINCTION
# 
# ext_prob_neutral = 0 # neutral probability of extinction
# 
# # Competition
# 
# e_0neg = 0.1 #0.15 # Asymptotic extinction probability with infinite negative interactions
# a_eneg = 0.025 # Shape of the exponential decay of the negative extinction - interaction relationship
# e_1neg = -e_0neg  # Extinction probability with absence of interactions
# 
# # Facilitation & Foodweb
# 
# e_0pos = 0.075 
# e_1pos = 5.19 
# a_epos = 1.2 
# 
# 
# #########################################
# 
# # Logistic function parameters
# k <- 10  # Steepness of the logistic curve
# midpoint <- 0.5  # Midpoint for the logistic curve