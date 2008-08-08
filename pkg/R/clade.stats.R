"clade.stats" <-
function(dataf, phyl, fun, ..., tips=FALSE){

	fun <- match.fun(fun)

	if(!is.data.frame(dataf)) stop("Requires data block to be in a dataframe")
	
	if(class(phyl) != "phylo") stop("Phylogeny required")

	name.check <- match(phyl$tip.label, rownames(dataf))
	
	if(any(is.na(name.check))){
		
		warning("Some tips in phylogeny	not present in data frame: aborting and returning missing tips.")
		return(phyl$tip.label[is.na(name.check)])
	
	} else {
	
		clades <- clade.members.list(phyl, tips=tips, tip.labels=TRUE)
		
		
		# get the rownames 
		dataf.names <- rownames(dataf)
		# ditch non-numeric columns
		dataf <- dataf[,unlist(lapply(dataf, is.numeric)), drop=FALSE]

		# ugly line that computes 'fun' for each numeric column in each subset of tips defined by internal nodes
		cstats <- t(sapply(clades, function(x){apply(dataf[match(x, dataf.names),,drop=FALSE], 2, fun, ...)}, simplify=TRUE))
		return(cstats)
	}
	
}

