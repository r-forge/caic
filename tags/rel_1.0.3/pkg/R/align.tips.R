align.tips <- function(phy, at=NULL){
	
    
	vcvPhy <- vcv.phylo(phy)
	
	if(is.null(at)) {
	    at <- max(diag(vcvPhy))
	} else {
	    if(at < max(vcvPhy[lower.tri(vcvPhy)])){
	        stop("User specified 'at' is closer to the root than at least one internal node.")
	    }
	}
	
	shortfall <- at - diag(vcvPhy)
	
	tipEdgeInd <- match(1:length(shortfall), phy$edge[,2])
	
	phy$edge.length[tipEdgeInd] <- phy$edge.length[tipEdgeInd] + shortfall
	
	return(phy)
	
}