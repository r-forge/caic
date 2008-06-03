`caic.label` <-
function(phy, charset=NULL, action="insert", style="CAIC"){
    
    # OLD2NEW STATUS: CONVERTED...

    if(! inherits(phy, "phylo")) 
         stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
        
    match.arg(action, c("insert", "replace", "append"))
    match.arg(style, c("RLE", "CAIC"))
    
    # OLD2NEW: contrGp <- rev(split(phy$edge[,2], as.numeric(phy$edge[,1]))) # edge is already numeric and root to tip order is not reversed
    # OLD2NEW: caicLab <- character(length(unique(as.character(phy$edge)))) # max(phy$edge) is the number of tips and nodes...
    # OLD2NEW: names(caicLab) <- c((-1:min(as.numeric(phy$edge))), (1:max(as.numeric(phy$edge)))) # names are 1:max(phy$edge)

    contrGp <- split(phy$edge[,2], f=phy$edge[,1]) # handily, split retains numeric order not alphabetic...
    caicLab <- character(max(phy$edge)) 
    names(caicLab) <- 1:max(phy$edge)


    if(is.null(charset)) charset <- LETTERS

    # loop the nodes
    for(nd in seq(along=contrGp)){

        parent <- names(contrGp)[nd]
        children <- contrGp[[nd]]
        if(length(children) > length(charset)) stop("Insufficient characters to code polytomies")
        caicLab[children] <- paste(caicLab[parent], charset[1:length(children)], sep="")

    }

    if(style=="RLE"){
        caicLab <- strsplit(caicLab, split="")
        caicLab <- sapply(caicLab, function(X) with(rle(X), paste(ifelse(lengths > 1, lengths, ""), values, sep="", collapse="")))
    }

    # OLD2NEW: intBool <- as.numeric(names(caicLab)) < 0 # internal nodes now from max(phy$edge)-phy$Nnode +1 to  max(phy$edge)
    intBool <- with(phy, 1:max(edge) > (max(edge) - Nnode))

    switch(action, 
        "replace" = { phy$tip.label <- caicLab[! intBool]
                      phy$node.label <- caicLab[intBool]},
        "append"  = { if(is.null(phy$node.label)) phy$node.label <- with(phy, (max(edge)-Nnode +1):max(edge))
                      phy$tip.label <- paste(phy$tip.label, caicLab[! intBool])
                      phy$node.label <- paste(phy$node.label, caicLab[intBool])},
        "insert"  =   phy$edge.caic.code <- caicLab[match(phy$edge[,2], names(caicLab))])

    return(phy)
}

