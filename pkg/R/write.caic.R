`write.caic` <-
function(phy, filebase, equal.brlen=NULL, charset=LETTERS){
    
    # OLD2NEW status: CONVERTED...
    
    phy <- caic.label(phy, action="insert", charset)
    
    # sort out branch lengths...
    if(! is.null(equal.brlen)){
        brlen <- rep(equal.brlen, length(phy$edge.length) + 1)
    } else {
        if(is.null(phy$edge.length)){
            brlen <- rep(2, length(phy$edge.length) + 1)
            warning("Phylogeny has no branch lengths and no equal branch length provided: using equal branches of length 2.")
        } else {
            brlen <- c(if(is.null(phy$root.edge)) 0 else phy$root.edge, phy$edge.length)
        }
    }
    
     # create a data.frame of CAIC codes, branch lengths and tip.labels including the root edge
    treeStruct <- data.frame( tip.label = I(c(NA, with(phy, tip.label[match(edge[,2], seq(along=tip.label))]))),
                              code = I(c("", with(phy, edge.caic.code[match(edge[,2], names(edge.caic.code))]))),
                              brlen = brlen,
                              node.depth = -9)
    
    # write a branch length file - the branch length associated with each code
    blen <- treeStruct[order(treeStruct[,2]), 2:4]
    
    write.table(blen, file=paste(filebase, ".blen", sep=""),
                sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    # write the .phyl file - the codes and names of the tips
    phyl <- subset(treeStruct, ! is.na(treeStruct$tip.label))
    phyl <- phyl[order(phyl[,2]), 2:1]

    phyl <- with(phyl, paste(code, tip.label, sep="\n"))
    cat(phyl, file=paste(filebase, ".phyl", sep=""), sep="\n")
    
    invisible(NULL)
}

