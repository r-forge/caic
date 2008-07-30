`caic.table` <-
function(caicObj, validNodes=TRUE, nodalValues=FALSE, ultrametric.tol=0.0001){
    # simple code to create a table of the contrasts from the caic object
    
        nodeNum <- matrix(as.numeric(names(caicObj$contrast.data$contrVar)),
                          ncol=1, dimnames=list(NULL, "nodeNumber"))
        contr <- with(caicObj$contrast.data$contr, cbind(response, explanatory))
        # colnames(contr) <- paste("CM.Mass", colnames(contr), sep="")
        if(nodalValues){
            nv <- with(caicObj$contrast.data$nodalVals, cbind(response, explanatory))
            colnames(nv) <- paste("NVM.Mass", colnames(nv), sep="")
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
        tab$studResid[match(as.numeric(names(stRes)), tab$nodeNumber)] <- stRes
        if(validNodes) tab <- subset(tab, validNodes, select=-validNodes)
       
        return(tab)

}

