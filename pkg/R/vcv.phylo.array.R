vcv.phylo.array <- function(phy){
    
    ## turns a phylogeny into a 3d array similar to a VCV matrix
    ## but keeping each beanch length separate. This is useful for 
    ## handling branch length transformations in functions where
    ## VCVs are used to handle phylogenetic structure
    
    ## rewritten to use new ape, via the clade matrix structure
    
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")

    cm <- clade.matrix(phy)
    cmM <- cm$clade.matrix
    cmE <- cm$edge.length
    nTip <- dim(cmM)[2]
    max.node.depth <- max(colSums(cmM))
    
    edge.array <- array(NA, dim=c(nTip, nTip, max.node.depth), dimnames=list(phy$tip.label, phy$tip.label, NULL))
    
     # must be a way of 'applying' or 'outering' this next bit off the clade matrix
    for(i in 1:nTip){
        for(j in 1:nTip){
        	
        	# get the shared edge lengths and insert into the array
        	shared.edge.lengths <- cmE[as.logical(cmM[,i] * cmM[,j])]
        	edge.array[i,j,seq(along=shared.edge.lengths)] <- shared.edge.lengths
        }
    }
    

    class(edge.array) <- "phylo.array"
    return(edge.array)
    
}