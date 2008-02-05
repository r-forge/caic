"clade.members" <-
function(x, phy, tip.labels=FALSE){
    
    # NEW2OLD: CONVERTED...
	
	# returns a vector of the tips that descend from an identified node
	if(class(phy) != "phylo") stop("Phylogeny required")
	
	NallNodes <- max(phy$edge)
	Ntips <- max(phy$edge) - phy$Nnode
	
	if(!(x %in% 1:NallNodes)) stop("Node not in range for phylogeny")
	
	# find the children of the node, append them to the vector of nodes (x)
	# and remove the parent, until all the nodes in the vector are tips...
	while(length(intN <- x[x>Ntips]) > 0){
		
		minIntN <- min(intN)
		childOfMinIntN <- with(phy, edge[,2][which(edge[,1]==minIntN)])
		
		x <- c(x[x != minIntN], childOfMinIntN)
	}
	
	x <- unique(x)
	
	if(tip.labels) {
		return(phy$tip.label[x])
	} else {
		return(x)
	}
}

