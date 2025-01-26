#setwd("C:\\Massol\\Encadrement\\2021- Thèse Fuster")


#data<-data.frame("spp"=c("A","B","C","D","E","F"),"ancestor"=c("C","C","G","F","F","G"),"A/E"=c(1,1,0,1,1,0),"distance"=c(2,2,3,4,4,1))

#write.tree(phy, file = "newick1.txt")

sister.group<-function(sis.names,sis.dist){
	n<-length(sis.dist)
	distances<-paste0(rep(":",n),sis.dist)
	sis<-paste0(sis.names,distances)
	res<-paste0(sis,collapse=",")
	res<-paste0("(",res,")")
	res
}

ToPhylo<-function(data){
	data.2<-data
	data.2$repr<-data$spp
	sisters<-levels(as.factor(data$spp))
	mothers<-levels(as.factor(data$ancestor))
	tips<-setdiff(sisters,mothers)
	root<-setdiff(mothers,sisters)
	foc.nodes<-unique(data[which(data$spp%in%tips),"ancestor"])
	n<-length(foc.nodes)
	data.2$repr[data.2$spp%in%tips]<-data.2$repr[data.2$spp%in%tips]
	while(n>1){
		foc.nodes2<-unique(data.2[which(data.2$spp%in%foc.nodes),"ancestor"])
		for(i in 1:n){
			daughters<-data.2[which(data.2$ancestor==foc.nodes[i]),"repr"]
			#print(daughters)
			daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[i]),"distance"]
			data.2$repr[data.2$spp==foc.nodes[i]]<-paste0(sister.group(daughters,daughters.dist),foc.nodes[i])
		}
		tips<-foc.nodes
		foc.nodes<-foc.nodes2
		n<-length(foc.nodes)
	}
	daughters<-data.2[which(data.2$ancestor==foc.nodes[1]),"repr"]
			#print(daughters)
	daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[1]),"distance"]
	paste0(sister.group(daughters,daughters.dist),root)
}

# phy<-ToPhylo(data)
# phy
# 
# 
# library(ape)
# tree<-read.tree(text = gsub("G",";",phy))
# plot(tree)
# write.tree(tree, file = "newick1.txt")


############################

ToPhylo2 <- function(data){
  data.2 <- data
  data.2$repr <- data$spp
  sisters <- levels(as.factor(data$spp))
  mothers <- levels(as.factor(data$ancestor))
  tips <- setdiff(sisters, mothers)
  root <- setdiff(mothers, sisters)
  
  # the root might be ancestor=0
  if(length(root) == 0) root <- 0
  foc.nodes <- unique(data[which(data$spp %in% tips), "ancestor"])
  n <- length(foc.nodes)
  data.2$repr[data.2$spp %in% tips] <- data.2$repr[data.2$spp %in% tips]
  
  while(n > 1){
    foc.nodes2 <- unique(data.2[which(data.2$spp %in% foc.nodes), "ancestor"])
    for(i in 1:n){
      daughters <- data.2[which(data.2$ancestor == foc.nodes[i]), "repr"]
      daughters.dist <- data.2[which(data.2$ancestor == foc.nodes[i]), "distance"]
      
      # This block handles the case where the ancestor is still extant (coexists with mutants)
      if(foc.nodes[i] %in% data.2$spp){
        daughters <- c(daughters, foc.nodes[i])
        daughters.dist <- c(daughters.dist, 0)  # Set ancestor's distance as 0
      }
      
      data.2$repr[data.2$spp == foc.nodes[i]] <- paste0(sister.group(daughters, daughters.dist), foc.nodes[i])
    }
    tips <- foc.nodes
    foc.nodes <- foc.nodes2
    n <- length(foc.nodes)
  }
  
  daughters <- data.2[which(data.2$ancestor == foc.nodes[1]), "repr"]
  daughters.dist <- data.2[which(data.2$ancestor == foc.nodes[1]), "distance"]
  
  paste0(sister.group(daughters, daughters.dist), root, ";")  # Ensuring the tree ends with a semicolon
}



# data <- res_sim$list_anc_dist[[10]]
# phy <- ToPhylo2(data)
# print(phy)
# 
# # Read the tree and plot
# library(ape)
# tree <- read.tree(text = phy)
# plot(tree)
# 
# phylo_distances <- cophenetic.phylo(tree)
# 
# # View the matrix of distances
# phylo_distances
