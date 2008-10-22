"pd.calc" <-
function(cm,  tip.subset=NULL, method="TBL", root.edge=FALSE){

	# check we have a valid method
	method <- match.arg(method, c("TBL", "MST", "UEH","SBL", "TIP"))
	
	# check we have a clade matrix and, if not, get one
	if(class(cm) != "clade.matrix"){
		if(class(cm) == "phylo"){
			warning("Converting phylo object to clade.matrix object")
			cm <- clade.matrix(cm)
		} else { 
			stop("pd.calc requires a phylogeny")	
		}
	}

    # if requested, drop the root edge
    nSpp <- dim(cm$clade.matrix)[2]
    if(! root.edge) cm$edge.length[nSpp + 1] <- 0
    	
	# subset the tips if requested
	if(!is.null(tip.subset)){
		# could be names - in which case they must match tip labels
		# could be numbers - in which case they must be in range
		switch(mode(tip.subset),
			"character" = {
				tip.subset <- match(tip.subset, cm$tip.label)
				if(any(is.na(tip.subset))){
					stop("Unmatched names in tip.subset")}},
			"numeric" = {
				if(any(tip.subset %in% 1:dim(cm$clade.matrix)[2] == FALSE)){
					stop("numeric tip.subset contains outside the range 1 to number of tips")}},
			stop("tip.subset must be either a vector of names or numbers"))
	} else {
		tip.subset <- 1:dim(cm$clade.matrix)[2]
	}
	
	#choose method
	switch(method,
		"TBL" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 0 
		},
		"MST" = {
			edge.in.matrix <- cm$clade.matrix[, tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 0 & edge.in.matrix < length(tip.subset)
		},
		"TIP" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 0 
			edge.in.matrix[(dim(cm$clade.matrix)[2]+1):length(edge.in.matrix)] <- FALSE
		},
		"UEH" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix == 1 
		},
        "SBL" = {
			edge.in.matrix <- cm$clade.matrix[,tip.subset]
			if(is.array(edge.in.matrix)){
				edge.in.matrix <- rowSums(edge.in.matrix)}
			edge.in.matrix <- edge.in.matrix > 1 
		})
	
	pd <- sum(cm$edge.len[edge.in.matrix])
	RET <- list(pd=pd, method=method)
		
	return(RET)
	
}

