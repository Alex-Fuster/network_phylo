########################################


# Function to append a list into an existing list
list_append <- function(lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}


# Function to generate a new set of traits for ancestors
# Each species is characterized by a set of 3 traits: n (niche), o (optimum) and r (niche's range)
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
    
    if(int == 0) {
      #n_m = -1
      #while((n_m < 0) || (n_m > 1)){n_m = n + rnorm(1, mean = 0, sd = pars$sd)}
      n_m = EnvStats::rnormTrunc(1, mean = n, sd = pars$sd, min = 0, max = 1)
      r_m = -1
      while((r_m < 0) || (r_m > 1)){r_m = r + rnorm(1, mean = 0, sd = 0.02)}
      traits_mut = c(n = n_m, r = r_m, o = n)
    }
    
    if(int == 1) {
      #n_m = -1
      #while((n_m < 0) || (n_m > 1)){n_m = n + rnorm(1, mean = 0, sd = pars$sd)}
      n_m = EnvStats::rnormTrunc(1, mean = n, sd = pars$sd, min = 0, max = 1)
      r_m = -1
      while((r_m < 0) || (r_m > 1)){r_m = r + rnorm(1, mean = 0, sd = 0.02)}
      traits_mut = c(n = n_m, r = r_m, o = n)
    }
    traits_mut
  })
}

########################################
# Function to compute the interaction network from a set of traits
get_L_mat = function(pars, traits_mat) {
  with(as.list(pars),{
    L = matrix(0, nr = pars$Smax, nc = pars$Smax)
    
    # Lower boundary
    low = traits_mat$n - traits_mat$r
    low_mat = matrix(low, nr = pars$Smax, nc = pars$Smax, byrow = TRUE)
    
    # Upper boundary
    high = traits_mat$n + traits_mat$r
    high_mat = matrix(high, nr = pars$Smax, nc = pars$Smax, byrow = TRUE)
    S = nrow(traits_mat)
    
    # Matrix of niche positions
    n_mat = matrix(traits_mat$n, nr = pars$Smax, nc = pars$Smax, byrow = FALSE)
    
    ## Add the basal species
    #n_basal = matrix(basal, nr = pars$Sbasal, nc = pars$Smax, byrow = FALSE)
    #n_mat = rbind(n_basal, n_mat)
    
    # Test interactions
    L[n_mat > low_mat & n_mat < high_mat] = 1
    if(pars$Smax > 1) diag(L) = 0
    L
  })
}

########################################
# Function to compute the interactions of a given species
get_L_vec = function(pars, traits_mat, traits_mut) {
	with(as.list(pars),{
		L_vec = numeric(pars$Smax)

		# Lower boundary
		low = as.numeric(traits_mut["n"]) - as.numeric(traits_mut["r"])

		# Upper boundary
		high = as.numeric(traits_mut["n"]) + as.numeric(traits_mut["r"])

		# Vector of niche positions
		n_vec = c(traits_mat$n)

	  	# Test interactions
	  	L_vec[n_vec > as.numeric(low) & n_vec < as.numeric(high)] = 1
	  	L_vec
  })
}

########################################
sim_model_bif_fc = function(seed, pars, nsteps) {

	with(pars,{

	set.seed(seed)

	# Draw the first species traits
	traits_mat = matrix(nr = Smax, nc = 3)
	traits_mat[1,] = rand_traits_anc(pars)
	traits_mat = as.data.frame(traits_mat)
	names(traits_mat) = c("n","r","o")

	# Set the presence/absence matrix
	pres = matrix(0, nr = nsteps, nc = Smax)
	pres[1,1] = 1

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
		ActualS = sum(pres[step-1,])
		#cat("Step = ", step-1," Total S = ", S," Actual S = ", ActualS,'\n')
		Stotact[step-1,] <- c(step-1, S, ActualS)
		
	
		
		

		# Test for successful origination probability
		for(esp in 1:Smax) {

			if(S >= Smax) break

			# Speciation occurs if the species is present
			if(pres[step-1,esp] == 1) {
			  
			  
			  # [TABLE DIST_ANC] Define distances of Alive spp -------------------
			 
			  dist_anc[esp,"A/E"] <- "A"
			  dist_anc[esp,"distance"] <- as.numeric(dist_anc[esp,"distance"])+1
			  

				# Species is maintained
				pres[step, esp] = 1

				# Test if there is speciation
				test_number = runif(1, 0, 1)
				speciation_prob = u_max/(1 + exp(d * (ActualS - I_max)))
				if(test_number < speciation_prob) {

					# Ancestor traits
					traits_anc <- traits_mat[esp,]

					# Pick new traits for the first mutant
					traits_mut = rand_traits_mut(traits_mat[esp,], pars)

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
					traits_mat[esp,] <- c(n_tmp, o_tmp, r_tmp)

					# Recompute the interactions
					I = get_L_vec(pars, traits_mat, traits_mut)

					# Compute the number of interactions among present species
					sum_I = sum(I*pres[step,]) + Bspe

					# Compute the probability of successful establishment #---------------------- P(establishment)
					if(int == 0){
					#	estab_prob = u_0neg + u_1neg*exp(-a_uneg * sum_I)
					  estab_prob_sel = u_0neg + u_1neg * exp(-a_uneg * sum_I)
					  estab_prob = SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
					}

					if(int == 1){
					#	estab_prob = (u_0 + u_1*exp(-a_u*sum_I))#/(1 + exp(pars$d * (sum_I - pars$I_max)))
						##print(paste("estab prob =",spec_prob, " sum of interactions =", sum_I))
						estab_prob_sel = u_0 + u_1 * exp(-a_u * sum_I)
						estab_prob = SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
					}

					# Test if there is establishment of the two new species

					## We start to test for one of them
					if(runif(1) < estab_prob) {
						S <- S + 1
						traits_mat[S,] = traits_mut
						pres[step,S] = 1
						anc[S,] = c(step, esp, S) #Time step, ancestor, new_species
						
						# [TABLE DIST_ANC] add new spp to table dist_anc ------------
						
						dist_anc[S, ] <- c(S, esp, "A", 1)
						

						if(S >= Smax) break
						
						speciation_mat[step,esp] = 1 # record speciation for ancestor

						pres[step,esp] = 0 #the ancestor disapear because of trait deplacement
					}

					# Reasign original traits to the ancestor
					traits_mat[esp,] <- traits_anc

					# Test for the second one

					# Pick new traits for the second mutuant
					traits_mut <- c(n_tmp, o_tmp, r_tmp)

					# Recompute the interactions
					I = get_L_vec(pars, traits_mat, traits_mut)

					# Compute the number of interactions among present species
					sum_I = sum(I*pres[step,]) + Bspe

					# Compute the probability of successful establishment #---------------------- P(establishment)
					if(int == 0){
					  #estab_prob = u_0neg + u_1neg*exp(-a_uneg * sum_I)
					  estab_prob_sel = u_0neg + u_1neg * exp(-a_uneg * sum_I)
					  estab_prob = SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
					}

					if(int == 1){
					  #estab_prob = (u_0 + u_1*exp(-a_u*sum_I))#/(1 + exp(pars$d * (sum_I - pars$I_max)))
						##print(paste("estab prob =",spec_prob, " sum of interactions =", sum_I))
					  estab_prob_sel = u_0 + u_1 * exp(-a_u * sum_I)
					  estab_prob = SN * estab_prob_neutral + (1 - SN) * estab_prob_sel
					}

					# Test if there is establishment of the second	species
					if(runif(1) < estab_prob) {
						S <- S + 1
						traits_mat[S,] = traits_mut
						pres[step,S] = 1
						anc[S,] = c(step, esp, S) #Time step, ancestor, new_species
						
						# [TABLE DIST_ANC] add new spp to table dist_anc ------------
						
						dist_anc[S, ] <- c(S, esp, "A", 1)

						if(S >= Smax) break

						speciation_mat[step,esp] = 2 # record speciation for ancestor
						pres[step,esp] = 0
					}
				}
			}
		}

		if(S >= Smax) break
		if(ActualS == 0) break

		# Compute the interaction matrix among present species
		pres_vec = pres[step,]
		cooc = matrix(pres_vec, nr=Smax, nc=Smax, byrow=TRUE)*
		matrix(pres_vec, nr=Smax, nc=Smax, byrow=FALSE)
		L = get_L_mat(pars, traits_mat)
		L = L*cooc
		L_list[[step]] = L                            #---------------------- List of networks/time-step

		# Test for extinction                         #---------------------- P(extinction)
		if(int == 0) {
			in_I = colSums(L) + Bext*pres_vec
		#	ext_prob = e_0neg * (1 - exp(-a_eneg*in_I))
			ext_prob_sel = e_0neg * (1 - exp(-a_eneg*in_I)) 
			ext_prob = SN * (ext_prob_neutral) + (1 - SN) * (ext_prob_sel)
			
		}

	 	if(int == 1) {
			in_I = colSums(L) + Bext*pres_vec
		#	ext_prob = e_0pos + e_1pos*exp(-a_epos*in_I)
			##print(paste("extinction probability =", ext_prob[1]))
			ext_prob_sel = e_0pos + e_1pos * exp(-a_epos * in_I) 
			ext_prob = SN * (ext_prob_neutral) + (1 - SN) * (ext_prob_sel)
		}

		
		# Perform extinctions
		
		#test_extprob = runif(S,0,1)
		#test_extprob = c(test_extprob, rep(0,pars$Smax - S))
		##print(paste("random number for the extinction prob = ", test_extprob[1]))
		present_spe <- grep(1,pres[step,])
		test_extprob <- rep(0,Smax)
		random_number <- runif(length(present_spe),0,1)
		test_extprob[present_spe] <- random_number
		#cat("S actuel = ", ActualS," longueur random_number = ", length(random_number),'\n')
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
		
	} # End of main loop

	list(pres = pres, 
	     traits = traits_mat, 
	     anc = anc, 
	     extinct = extinct, 
	     dynamic = Stotact, 
	     L = L, 
	     L_list = L_list,
	     dist_anc = dist_anc,
	     list_dist_anc = list_dist_anc,
	     speciation_matrix = speciation_mat,
	     extinction_matrix = extinction_mat)
	})
}


