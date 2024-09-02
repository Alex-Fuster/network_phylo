

tree_raw <- ToPhylo(list_anc_dist_letters[[50]])
tree_raw <- paste(tree_raw, "root")
tree_raw <- read.tree(text = sub("A root",";",tree_raw))

tree_corrected<-tree_raw
tree_corrected$edge.length<-sapply(tree_corrected$edge.length,function(x) ifelse(x==0,1e-5,x))
phylo.vcv<-vcv(tree_corrected)
phylo.corr<-cov2cor(phylo.vcv)

phylo.corr_cropped <- phylo.corr[names(which(presence_matrix[50,] == 1)), 
                                        names(which(presence_matrix[50,] == 1))]


list_phylo.corr_cropped[[50]]

#match_phylogeny_interaction <- match(rownames(data.inter), rownames(phylo.corr))

svd_eigen.phy<-eigen(list_phylo.corr_cropped[[50]], symmetric = T)$vec



# Perform SVD
matrix <- list_net_present_spp.letters[[50]]

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
