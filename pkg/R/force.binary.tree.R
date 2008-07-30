"force.binary.tree" <-
function(phylo, arb.branch =0, node.suffix="M.MassPoly"){

## Written by David Orme
## Takes a polytomous tree and arbitrarily resolves the polytomies
## into a sequence of dichotomies with edge length equal to arb.branch
## I can't see any reason for this to be anything other than zero but
## using no-zero values was helpful in debugging the function.

## - rewritten completely to handle the new ape phylo structure
## - amended to handle node labels
if(class(phylo) != "phylo") stop(substitute(phylo), " is not of class 'phylo'.")

if(is.binary.tree(phylo)){
    
    cat(substitute(phylo), " is already binary.")
    invisible(phylo)
    
} else {
    
    ## get vector of polytomous node ids in the order of root to tip...
    polytomies <- sort(as.numeric(names(which(table(phylo$edge[,1]) > 2))))
    
    ## loop through the polytomies:
    while(length(polytomies) > 0) {
        
        current.polytomy <- polytomies[1]
        ## get the indices of the edges forming the polytomy
        poly.edges <- which(phylo$edge[,1] == current.polytomy)
        nPol <- length(poly.edges)
        
        ## insert a gap into the edge numbering that takes account of the number of nodes
        phylo$edge <- phylo$edge + ifelse(phylo$edge > current.polytomy, nPol - 2, 0)
        
        ## add in the internal node numbers of the new parents leading to each tip 
        phylo$edge[poly.edges,1] <- phylo$edge[poly.edges,1] + c(0:(nPol-2), nPol-2)
        
        ## add the internal links between the new parents
        new.parents <- matrix(current.polytomy + c(0:(nPol-3), 1:(nPol-2)), ncol=2)
        phylo$edge <- rbind(phylo$edge, new.parents)
        
        ## handle node labels
        if(! is.null(phylo$node.label)){
            
            ## duplicate existing (or default N1234 style) into new labels with suffixes
            nTip <- length(phylo$tip.label)
            currNodeLabInd <- current.polytomy - nTip
            currNodeLab <- phylo$node.label[currNodeLabInd]
            if(length(currNodeLab) == 0) currNodeLab <- paste("N", current.polytomy, sep="")
            newNodeLab <- paste(currNodeLab, node.suffix, 1:(nPol-1), sep="")
            
            ## insert into the node label vector, allowing for a root polytomy or final node polytomy
            nodeID <- seq(along=phylo$node.label) + nTip 
            nodeGrp <- cut(nodeID, breaks=c(nTip, current.polytomy-0.5, current.polytomy+0.5,
                           nTip+ phylo$Nnode +2), labels=c("Below","Node","Above"))
            splitNode <- split(phylo$node.label,nodeGrp)
            splitNode[[2]] <- newNodeLab
            phylo$node.label <- unlist(splitNode)
            names(phylo$node.label) <- NULL # remove names - simply internal coding 
        }
        ## and add the new branch lengths and Nnodes
        phylo$edge.length <- c(phylo$edge.length, rep(arb.branch, nPol-2))
        phylo$Nnode <- phylo$Nnode + nPol-2
        
        ## drop the one that has been handled and make sure the rest match their new names...
        polytomies <- polytomies[-1] + nPol - 2
        
    }

    ## the resulting order causes plot.phylo to hang, so to be defensive
    phylo <- reorder(phylo, "cladewise")
    attr(phylo, "force.binary") <- TRUE
    return(phylo)
}   
}

