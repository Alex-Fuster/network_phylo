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
# Function to generate a new set of traits from mutants
rand_traits_mut = function(traits_anc, pars) {
  
  with(as.list(c(traits_anc, pars)),{
  
    if(int == 2) {
     # a = beta_n*n/(1-n)
      n_m = EnvStats::rnormTrunc(1, mean = n, sd = pars$sd, min = 0, max = 1)
      #n_m = EnvStats::rnormTrunc(1, mean = n, sd = pars$sd, min = 0, max = 1)
      #			o_m = runif(1,0,n_m)
      o_m = n_m/2
      traits_mut = c(n = n_m, r = r, o = o_m)
    }
    traits_mut 
  })
}

########################################
# Function to compute the interaction network from a set of traits
get_L_mat = function(basal, pars, traits_mat) {
  with(as.list(pars),{
    L = matrix(0, nr = Smax+Sbasal, nc = Smax)
    
    # Lower boundary
    low = traits_mat$o - traits_mat$r
    low_mat = matrix(low, nr = Smax+Sbasal, nc = Smax, byrow = TRUE)
    
    # Upper boundary
    high = traits_mat$o + traits_mat$r
    high_mat = matrix(high, nr = Smax+Sbasal, nc = Smax, byrow = TRUE)	
    S = nrow(traits_mat)
    
    # Matrix of niche positions
    n_mat = matrix(traits_mat$n, nr = Smax, nc = Smax, byrow = FALSE)
    
    # Add the basal species
    n_basal = matrix(basal, nr = Sbasal, nc = Smax, byrow = FALSE)
    n_mat = rbind(n_basal, n_mat)
    
    # Test interactions
    L[n_mat > low_mat & n_mat < high_mat] = 1
    if(Smax > 1) diag(L[(Sbasal+1):(Sbasal+Smax),]) = 0
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
    anc = matrix(NA,nr = Smax, nc = 3)
    
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
    
    
    # Record speciation and extinction events per species
    
    speciation_mat = matrix(0, nr = nsteps, nc = Smax)
    extinction_mat = matrix(0, nr = nsteps, nc = Smax)
    
    
    # Species count
    S = 1
    
    Stotact <- matrix(NA, nsteps, 3)
    Stotact <- as.data.frame(Stotact)
    names(Stotact) = c("Step","TotalS","ActualS")
    
    ##########################
    # MAIN LOOP
    for(step in 2:nsteps) {
      ActualS = sum(pres[step - 1,])
      cat("Step = ", step - 1, " Total S = ", S, " Actual S = ", ActualS, '\n')
      Stotact[step-1,] <- c(step-1, S, ActualS)
      
      
      
      # Test for successful speciation probability
      for(i in 1:Smax) {
        
        if(S >= Smax) break
        
        # Speciation occurs if the species is present
        if(pres[step-1, i] == 1) {
          
          
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
            
            # Pick new traits for the first mutant
            traits_mut = rand_traits_mut(traits_mat[i,], pars)
            
            # Simulation of divergente selection
            ## Define which "side" of the ancestor the second new species will be
            if(traits_mut[1]<traits_anc[1]){
              limit_inf <- as.numeric(traits_anc[1])
              limit_sup <- 1
            }else{
              limit_inf <- 0
              limit_sup <- as.numeric(traits_anc[1])
            }
            
            # niche avlue is different from the other new species, but the optimum and the interaction range are the same
            
            n_tmp <- EnvStats::rnormTrunc(1, mean = as.numeric(traits_anc[1]), sd = sd, min = limit_inf, max = limit_sup)
            o_tmp <- traits_mut[2]
            r_tmp <- traits_mut[3]
            
            # Asign new traits to the ancestor to fake the new second species
            traits_mat[i,] <- c(n_tmp, o_tmp, r_tmp)
            
            # Recompute the interactions 
            I = get_L_vec(basal, pars, traits_mat, traits_mut)
            
            # Compute the number of interactions among present species
            sum_I = sum(I * c(rep(1, Sbasal), pres[step, ]))	
            
            # Compute the probability of establishment	
            if(int == 0) {
              estab_prob_sel = u_0neg + u_1neg * exp(-a_uneg * sum_I)
              estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)
            }
            
            if(int == 1 | int == 2) {
              estab_prob_sel = u_0pos + u_1pos * exp(-a_upos * sum_I)
              estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)

            }
            
            # Test if there is establishment of the two new species
            
            ## We start to test for one of them
            
            if(runif(1) < estab_prob) {
              S <- S + 1
              traits_mat[S,] = traits_mut
              pres[step,S] = 1
              anc[S,] = c(step, i, S) #Time step, ancestor, new_species
              
              # [TABLE DIST_ANC] add new spp to table dist_anc ------------
              
              dist_anc[S, ] <- c(S, i, "A", 1)
              
              
              if(S >= Smax) break
              
              speciation_mat[step,i] = 1 # record speciation for ancestor
              
              pres[step,i] = 0 #the ancestor disapear because of trait deplacement
            }
            
            
            # Reasign original traits to the ancestor
            traits_mat[i,] <- traits_anc
            
            # Test for the second one
            
            # Pick new traits for the second mutuant
            traits_mut <- c(n_tmp, o_tmp, r_tmp)
            
            # Recompute the interactions 
            I = get_L_vec(basal, pars, traits_mat, traits_mut)
            
            # Compute the number of interactions among present species
            sum_I = sum(I * c(rep(1, Sbasal), pres[step, ]))
            
            # Compute the probability of establishment	
            if(int == 0) {
              estab_prob_sel = u_0neg + u_1neg * exp(-a_uneg * sum_I)
              estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)
            }
            
            if(int == 1 | int == 2) {
              estab_prob_sel = u_0pos + u_1pos * exp(-a_upos * sum_I)
              estab_prob = SN * (estab_prob_neutral[i]) + (1 - SN) * (estab_prob_sel)

            }
            
            # Test if there is establishment of the second species
            
            if(runif(1) < estab_prob) {
              S <- S + 1
              traits_mat[S,] = traits_mut
              pres[step,S] = 1
              anc[S,] = c(step, i, S) #Time step, ancestor, new_species
              
              # [TABLE DIST_ANC] add new spp to table dist_anc ------------
              
              dist_anc[S, ] <- c(S, i, "A", 1)
              
              
              if(S >= Smax) break
              
              speciation_mat[step,i] = 1 # record speciation for ancestor
              
              pres[step,i] = 0 #the ancestor disapear because of trait deplacement
            }
          }
        }
      }
      
      if(S >= Smax) break
      
      
      # Compute the interaction matrix among present species
      pres_vec = pres[step, ]
      cooc = matrix(pres_vec, nr = Smax, nc = Smax, byrow = TRUE) * 
        matrix(pres_vec, nr = Smax, nc = Smax, byrow = FALSE)
      L = get_L_mat(basal, pars, traits_mat)
      L[c((Sbasal + 1):(Sbasal + Smax)), ] = L[c((Sbasal + 1):(Sbasal + Smax)), ] * cooc
      L_list[[step]] = L
      
      
      
      # Test for extinction
      if(int == 0) {
        in_I = colSums(L)
        ext_prob_sel = e_0neg + e_1neg * (1 - exp(-a_eneg * in_I)) 
        ext_prob = SN * (ext_prob_neutral) + (1 - SN) * (ext_prob_sel)
      }
      
      if(int == 1) {
        in_I = colSums(L)
        ext_prob_sel = e_0pos + e_1pos * exp(-a_epos * in_I) 
        ext_prob = SN * (ext_prob_neutral) + (1 - SN) * (ext_prob_sel)
      }
      
      if(int == 2) {
        in_I = colSums(L)
        out_I = rowSums(L)[(Sbasal + 1):(Sbasal + Smax)]
        ext_prob_sel = e_0neg + e_1neg * exp(-a_eneg * out_I) + 
          e_0pos + e_1pos * exp(-a_epos * in_I)
        
        
        # Compute competition index using matrix operations
        shared_victims = L[(Sbasal + 1):(Sbasal + Smax), ] %*% t(L[(Sbasal + 1):(Sbasal + Smax), ])
        diag(shared_victims) = 0  # Remove self-competition by setting diagonal to 0
        
        # Compute the number of victims for each species
        num_victims = rowSums(L[(Sbasal + 1):(Sbasal + Smax), ])
        
        # Avoid division by zero
        num_victims[num_victims == 0] = 1
        
        # Compute the proportion of shared victims
        prop_shared_victims = shared_victims / num_victims
        
        # Compute the competition index as the row sums of the proportion of shared victims
        competition_index = rowSums(prop_shared_victims)
        
        # Adjust extinction probability with competition
        competition_factor = 1 + competition_coefficient * competition_index
        ext_prob_sel = ext_prob_sel * competition_factor
        
        ext_prob = SN * (ext_prob_neutral) + (1 - SN) * (ext_prob_sel)
 
      }
      
      
      # Perform extinctions
      
      
      present_spe <- grep(1,pres[step,])
      test_extprob <- rep(0,Smax)
      random_number <- runif(length(present_spe),0,1)
      test_extprob[present_spe] <- random_number
      
      pres[step, pres[step-1,] & test_extprob < ext_prob] = 0
      
      for(i in 1:Smax) {
        if(pres[step,i] != pres[step-1,i] & pres[step,i] == 0){
          extinct[i,] = c(step, i)
          
          # [TABLE DIST_ANC] Define a new E spp -------------------
          
          dist_anc[i,"A/E"] <- "E"
          
          extinction_mat[step,i] = 1 # record extinction
          
        }
      }
      
      
      # [TABLE DIST_ANC] save table dist_anc at timestep step
      
      list_dist_anc[[step]] <- na.omit(dist_anc)
      
      
      # check sim
      
      if(step == 30 & sum(pres[step,]) <= 5){
        
        cat("Simulation stopped due to low species presence.\n")
        break
      }
      
      
    } # End of main loop
    
    list(pres = pres, 
         traits = traits_mat, 
         anc = anc, 
         extinct = extinct, 
         dynamic = Stotact, 
         L_list = L_list, 
         basal = basal,
         dist_anc = dist_anc,
         list_dist_anc = list_dist_anc,
         speciation_matrix = speciation_mat,
         extinction_matrix = extinction_mat)
    
  })
}