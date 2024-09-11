setwd("C:\\Massol\\Encadrement\\2021- Thèse Fuster")


data<-data.frame("spp"=c("A","B","C","D","E","F"),"ancestor"=c("C","C","G","F","F","G"),"A/E"=c(1,1,0,1,1,0),"distance"=c(2,2,3,4,4,1))

write.tree(phy, file = "newick1.txt")

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

phy<-ToPhylo(data)
phy


library(ape)
tree<-read.tree(text = gsub("G",";",phy))
plot(tree)
write.tree(tree, file = "newick1.txt")
