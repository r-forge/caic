`caic.table` <-
function(caicObj, validNodes=TRUE, nodalValues=FALSE, ultrametric.tol=0.0001, CAIC.codes=FALSE){
    # simple code to create a table of the contrasts from the caic object
    
        nodeNum <- matrix(as.numeric(names(caicObj$contrast.data$contrVar)),
                          ncol=1, dimnames=list(NULL, "nodeNumber"))
        contr <- with(caicObj$contrast.data$contr, cbind(response, explanatory))
        # colnames(contr) <- paste("C_", colnames(contr), sep="")
        if(nodalValues){
            nv <- with(caicObj$contrast.data$nodalVals, cbind(response, explanatory))
            colnames(nv) <- paste("nodal.", colnames(nv), sep="")
            tab <- as.data.frame(cbind(nodeNum, contr, nv))
        } else {
            tab <- as.data.frame(cbind(nodeNum, contr))
        }
        
        tab$contrVar <- caicObj$contrast.data$contrVar
        tab$validNodes <- caicObj$contrast.data$validNodes
        tab$nChild <- caicObj$contrast.data$nChild
        tab$nodeDepth <- caicObj$contrast.data$nodeDepth
        if(is.ultrametric(caicObj$phy, tol=ultrametric.tol)) tab$nodeAge <- branching.times(caicObj$phy) else tab$nodeAge <- NA
        stRes <- rstudent(caicObj$mod)
        tab$studResid <- NA
        tab$studResid[match(as.numeric(names(stRes)), tab$nodeNumber)] <- stRes
        if(validNodes) tab <- subset(tab, validNodes, select=-validNodes)
       
       
       if(CAIC.codes){
           Cphy <- caic.label(caicObj$phy)
           tab$CAIC.code <- Cphy$edge.caic.code[match(tab$nodeNumber, names(Cphy$edge.caic.code))]
       }
       
       return(tab)

}

