"clade.stats" <-
function(dat, phy, fun, ..., tips=FALSE){

	fun <- match.fun(fun)

	if(!is.data.frame(dat)) stop("Requires data block to be in a data frame")
	
	if(class(phy) != "phylo") stop("Phylogeny required")

	name.check <- match(phy$tip.label, rownames(dat))
	
	if(any(is.na(name.check))){
		
		warning("Some tips in phylogeny not present in data frame: aborting and returning missing tips.")
		return(phy$tip.label[is.na(name.check)])
	
	} else {
	
		clades <- clade.members.list(phy, tips=tips, tip.labels=TRUE)
		
		
		# get the rownames 
		dat.names <- rownames(dat)
		# ditch non-numeric columns
		dat <- dat[,unlist(lapply(dat, is.numeric)), drop=FALSE]

		# ugly line that computes 'fun' for each numeric column in each subset of tips defined by internal nodes
		cstats <- t(sapply(clades, function(x){apply(dat[match(x, dat.names),,drop=FALSE], 2, fun, ...)}, simplify=TRUE))
		
		# cope with the way that sapply deals with column vectors
		if(dim(dat)[2] == 1)  cstats <- t(cstats) 
        
        cstats <- as.data.frame(cstats)
        rownames(cstats) <- names(clades)
        colnames(cstats) <- colnames(dat)
		return(cstats)
	}
	
}

