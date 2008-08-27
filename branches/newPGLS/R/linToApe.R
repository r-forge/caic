linToApe <- function(linObj){
    
    lineages <- linObj$lineages
    clade <- linObj$clade
    
    if(dim(lineages)[1] > 1){
        # extract information (excluding root node)
        linNoR <- lineages[-1,]
            
        # now need to coerce the numbering into ape style
        parentMap <- data.frame(linPar=unique(linNoR$parent.id), apePar=with(clade, (nTip+1):nLin))
        childMap <- data.frame(linPar=with(linNoR, id[tip]), apePar=with(clade, 1:nTip))
        nodeMap <- rbind(parentMap, childMap)
        linNoR$pnode <- nodeMap$apePar[match(linNoR$parent.id, nodeMap$linPar)]
        linNoR$node <- nodeMap$apePar[match(linNoR$id, nodeMap$linPar)]
        
        edge <- as.matrix(linNoR[, c("pnode", "node")])
        dimnames(edge) <- NULL
        
        # lets be honest... the caic.code column is really only there as a cheap
        # route to a cladewise sorting of the edge matrix!
        ord <- order(linNoR$caic.code)
        edge <- edge[ord,]
        edge.length <- linNoR$lin.age[ord]
        
        phy <- list(edge=edge, edge.length=edge.length, tip.label=1:clade$nTip, 
                    Nnode=with(clade, nLin-nTip), root.edge=lineages$lin.age[1])              
        class(phy) <- "phylo"
        
        lastRules <- linObj$epochRules[[length(linObj$epochRules)]]
        
        if(! is.null(lastRules$ct.set)){
            phy$ct.data <- subset(linNoR, select=c("node", names(lastRules$ct.set$ct.start)))[ord,]
        }
        
        if(! is.null(lastRules$dt.rates)){
            phy$dt.data <- subset(linNoR, select=c("node", names(lastRules$dt.rates)))[ord,]
        }

    } else {
        # have an unspeciated root so put in slightly differently as a single tip
        phy <- list(edge=matrix(c(1,2), ncol=2), edge.length=lineages$lin.age[1], tip.label=1)
        class(phy) <- "phylo"
        if(! is.null(linObj$ct.set)){
            phy$ct.data <- subset(lineages, select=c(node="id", names(linObj$ct.set$ct.start)))
            phy$ct.set <- linObj$ct.set
        }
        if(! is.null(linObj$dt.rates)){
            phy$dt.data <- subset(lineages, select=c("node", names(linObj$dt.rates)))
            phy$dt.rates <- linObj$dt.rates
        }
   }
     
    phy$epochRules <- linObj$epochRules
    
    return(phy)
}