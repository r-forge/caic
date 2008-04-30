`node2tip` <-
function(phy){
	
	# OLD2NEW: CONVERTED
	
	# returns the maximum number of nodes between each internal node and the tips of the tree
	# including the internal node itself and the tip (i.e. is 2 for the parent of dichotomous tip)
	
	clm <- clade.matrix(phy)$clade.matrix
	# want the maximum number of times a tip descending from an internal node has appeared in
	# less nested nodes
	
	clm <- clm[order(rowSums(clm)),] # order by number of descendants
	nd <- apply(apply(clm, 2, cumsum) * clm, 1, max)
	names(nd) <- rownames(clm)
	nd <- nd[order(as.numeric(names(nd)))]
	return(nd)
}

