"all.clades" <-
function(phy, tips=FALSE, tip.labels=FALSE){
    
    # OLD2NEW CONVERTED
    
	# returns a list of vectors showing the tips
	# subtending from each node in the tree
	if(class(phy) != "phylo") stop("Phylogeny required")

	nodes <- 1:max(phy$edge)
	
	if(!tips) nodes <- nodes[nodes > length(nodes) - phy$Nnode]
	
	clade.list <- mapply(clade.members, nodes, MoreArgs=list(phy=phy, tip.labels=tip.labels))
	names(clade.list) <- nodes
	
	return(clade.list)
}

