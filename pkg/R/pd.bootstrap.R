"pd.bootstrap" <-
function(cm, ntips, reps=1000, method="TBL", tip.weights=NULL){
    
	# check we have a valid method
	method <- match.arg(method, c("TBL", "MST", "UEH", "SBL", "TIP"))

	# check we have a clade matrix and, if not, get one
	if(class(cm) != "clade.matrix"){
		if(class(cm) == "phylo"){
			warning("Converting phylo object to clade.matrix object")
			cm <- clade.matrix(cm)
		} else { 
			stop("pd.calc requires a phylogeny")	
		}
	}
	
	# check for sensible sample
	total.nb.tips <- dim(cm$clade.matrix)[2]
	if(!(ntips %in% 1:(total.nb.tips - 1))){
		stop("'sample' must be a positive integer lower than the number of tips")}
	
	#set up the store
	pd.store <- numeric(reps)
	tips <- 1:total.nb.tips
	
	# if there are weights make sure they go in the right place...
	if( ! is.null(tip.weights)){
	        
        # if the vector is named then match the order to tips
        if(! is.null(names(tip.weights))) {
            wght.match <- match(cm$tip.label, names(tip.weights))
            
            if(any(is.na(wght.match))){ # this is not elegant but can't work out how to stop and return
                warning("The returned tip labels have no matching named element in tip.weights")
                return(cm$tip.label[is.na(wght.match)])
                }
            
            tip.weights <- tip.weights[wght.match]
            
        } else {
             stop("'weights' must be a vector of weights, named to match the tip labels")
        }
    }
    
	# get the pd values
    for(rep in seq(along=pd.store)){
        which.tips <- sample(tips, ntips, prob=tip.weights)
        pd.store[rep] <- pd.calc(cm, tip.subset=which.tips, method=method)$pd
	}
	
	return(list(pd.distrib=pd.store, method=method))
}

