
install.packages("https://cran.r-project.org/src/contrib/Archive/geomorph/geomorph_3.3.1.tar.gz", repos = NULL, type = "source", dependencies = TRUE)

install.packages("RRPP")
install.packages("devtools")
library(devtools)
install_github("mlcollyer/RRPP")
library(geomorph)
data(plethspecies) 
Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment    

#Test for phylogenetic signal in shape
PS.shape <- physignal(A = Y.gpa$coords, phy = plethspecies$phy)
summary(PS.shape)
plot(PS.shape)
plot(PS.shape$PACA, phylo = TRUE)
PS.shape$K.by.p # Phylogenetic signal profile

#Test for phylogenetic signal in size
PS.size <- physignal(A = Y.gpa$Csize, phy = plethspecies$phy)
summary(PS.size)
plot(PS.size)

install.packages("installr")
library(installr)
updateR()


