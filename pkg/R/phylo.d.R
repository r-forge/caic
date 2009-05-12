phylo.d <- function(ds, phy, binvar, tipnames, permut=1000, verbose=TRUE) {
	
	# data checking
	if (!is.data.frame(ds)) (stop("'", deparse(substitute(ds)), "' must be an object of class 'data.frame'."))
	if (!inherits(phy, "phylo")) (stop("'", deparse(substitute(phy)), "' is not of class 'phylo'."))
	binvar <- deparse(substitute(binvar))
    bininds <- match(binvar, names(ds))
    if (is.na(bininds)) (stop("'", binvar, "' is not a column name in ds."))
    tipnames <- deparse(substitute(tipnames))
    tipsinds <- match(tipnames, names(ds))
    if (is.na(tipsinds)) (stop("'", tipnames, "' is not a column name in ds."))
    if (!is.numeric(permut)) (stop("'", permut, "' is not numeric.")) 
    phyName <- deparse(substitute(phy))
    dsName <- deparse(substitute(ds))

	d.data <- function (ds, phy, bininds, tipsinds, permut=1000, verbose=FALSE) {
		
	    # make up 1000 random data huffles
		if (verbose) (print("Computing phylogenetically random permutations..."))
		kseq <- seq(length(ds)+1, by=1, length.out=permut)
		for (k in kseq) { 
			ds[,k] <- sample(ds[,bininds])
			names(ds)[k] <- paste("Var", k, sep="")
		}
		names(ds)[tipsinds] <- "tipnames"
		explvarsran <- paste(paste("Var", kseq, sep=""), collapse="+")
		
		# make up 1000 simulated characters, brownian, lambda=0.999
		if (verbose) (print("Computing phylogenetically clumped permutations..."))
	    lamtree <- lambdaTree(phy, 1)
	    simmat <- matrix(rep(0,permut^2), nrow=permut, ncol=permut) # covariance between traits is zero
	    for (i in 1:permut) (simmat[i,i] <- 1) # trait variance is 1
	    simchar <- sim.char(lamtree, simmat)
	    simds <- as.data.frame(simchar[1:length(ds[,1]),1:permut,1])
	    simds$tipnames <- rownames(simds)
	    simds$resp <- sample(simds[,1])
	    propone <- length(which(ds[,bininds]==1)) / length(ds[,bininds]) # proportion of ones in the response variable
	    
	    # convert the simulated traits into binary traits using the proportion of ones in the response variable as threshold
	    for (k in 1:permut) {
		    tipvals <- simds[,k]
		    names(tipvals) <- rownames(simds)
		    tipvals <- sort(tipvals)
			ones <- names(tipvals)[(length(tipvals)-round(length(tipvals)*propone, 0)+1):length(tipvals)]
			simds[which(is.element(rownames(simds), ones)),k] <- round(1,0)
			simds[which(!is.element(rownames(simds), ones)),k] <- round(0,0)
		}
		explvarssim <- paste(paste("V", c(1:permut), sep=""), collapse="+")
		
		return(list(ds.ran=ds, formula.ran=paste(names(ds)[bininds], "~", explvarsran, sep=""), ds.sim=simds, formula.sim=paste("resp~", explvarssim, sep="")))
		
	}

	
	d.fit <- function(dds, phy, verbose=FALSE) {
		
		if (verbose) (print("Getting nodal value estimates..."))
		caicran <- crunch(as.formula(dds[[2]]), dds[[1]], phy, names.col=tipnames, polytomy.brlen=0)
		caicsim <- crunch(as.formula(dds[[4]]), dds[[3]], phy, names.col=tipnames, polytomy.brlen=0)
		return(list(caicran=caicran, caicsim=caicsim, ds.ran=dds$ds.ran, ds.sim=dds$ds.sim))
		
	}
	
	sumofdiffs <- function (caicobj, ds, tipnames, resp=FALSE) {
		
		# data checking
		if (!is.data.frame(ds)) (stop("'", deparse(substitute(ds)), "' must be an object of class 'data.frame'."))
		if (!inherits(caicobj, "caic")) (stop("'", deparse(substitute(caicobj)), "' is not of class 'caic'."))
	    tipnames <- deparse(substitute(tipnames))
	    tipsinds <- match(tipnames, names(ds))
	    if (is.na(tipsinds)) (stop("'", tipnames, "' is not a column name in ds."))
	    
		edges <- as.data.frame(caicobj$phy$edge)
		edges$bl <- caicobj$phy$edge.length
		responsev <- attr(caicobj$contrast.data$contr$response, "dimnames")[[2]]
		variables <- attr(caicobj$contrast.data$contr$explanatory, "dimnames")[[2]]
		
		# response variable
		if (resp) {
			nvs <- caicobj$contrast.data$nodalVals$response[,1]
			edges$startnv <- nvs[match(edges[,1], names(nvs))]
			tips <- ds[match(caicobj$phy$tip.label, ds[,tipsinds]),which(names(ds)==responsev)]
			names(tips) <- c(1:length(tips))
			allnvs <- c(nvs, tips)
			allnvs <- allnvs[order(as.numeric(names(allnvs)))]
			edges$endnv <- allnvs[match(edges[,2], names(allnvs))]
			edges$change <- abs(edges$startnv-edges$endnv)
			rtvr <- sum(edges$change)
			names(rtvr) <- responsev
			obsnvs <- nvs
		} else {
			rtvr <- NA
			obsnvs <- NA
		}
		
		# explanatory variables
		rtv <- numeric(0)
		for (i in seq(along=variables)) {
			nvs <- caicobj$contrast.data$nodalVals$explanatory[,which(attr(caicobj$contrast.data$nodalVals$explanatory, "dimnames")[[2]]==variables[i])]
			edges$startnv <- nvs[match(edges[,1], names(nvs))]
			tips <- ds[,which(names(ds)==variables[i])][match(caicobj$phy$tip.label, ds[,tipsinds])]
			names(tips) <- c(1:length(tips))
			allnvs <- c(nvs, tips)
			allnvs <- allnvs[order(as.numeric(names(allnvs)))]
			edges$endnv <- allnvs[match(edges[,2], names(allnvs))]
			edges$change <- abs(edges$startnv-edges$endnv)
			rtv <- append(rtv, sum(edges$change))
			names(rtv)[i] <- variables[i]
			if (i==1) (nvsds <- data.frame(t(nvs))) else (nvsds[i,] <- nvs)
		}
		
		return(list(sumoc.resp=rtvr, sumoc.expl=rtv, obs.nvs=obsnvs, perm.nvs=nvsds))
	}
	
	
	d.compute <- function(dfit, verbose=FALSE) {
		
		if (verbose) (print("Computing output..."))
		ransocc <- sumofdiffs(dfit$caicran, dfit$ds.ran, tipnames, resp=TRUE)
		simsocc <- sumofdiffs(dfit$caicsim, dfit$ds.sim, tipnames)
		soccratio <- (as.numeric(ransocc$sumoc.resp)-mean(simsocc$sumoc.expl)) / (mean(ransocc$sumoc.expl)-mean(simsocc$sumoc.expl))
		soccpval1 <- length(which(ransocc$sumoc.expl<as.numeric(ransocc$sumoc.resp))) / length(ransocc$sumoc.expl)
		soccpval0 <- length(which(simsocc$sumoc.expl>as.numeric(ransocc$sumoc.resp))) / length(simsocc$sumoc.expl)	
		ret <- list(DEstimate=soccratio, Pval1=soccpval1, Pval0=soccpval0, Parameters=list(Observed=as.numeric(ransocc$sumoc.resp), MeanRandom=mean(ransocc$sumoc.expl), MeanBrownian=mean(simsocc$sumoc.expl)), Permutations=list(random=ransocc$sumoc.expl, brownian=simsocc$sumoc.expl), NodalVals=list(observed=ransocc$obs.nvs, random=ransocc$perm.nvs, brownian=simsocc$perm.nvs))
		# if (verbose) (cat(" D estimate - ", ret$DEstimate, "\n P value for D<1 - ", ret$Pval1, "\n P value for D>0 - ", ret$Pval0, "\n"))
		return(ret)
	}
	

	dds <- d.data(ds=ds, phy=phy, bininds=bininds, tipsinds=tipsinds, permut=permut, verbose=verbose)
	dfit <- d.fit(dds=dds, phy=phy, verbose=verbose)
	dvals <- d.compute(dfit, verbose=verbose)
	dvals$binvar <- binvar
	dvals$phyName <- phyName
	dvals$dsName <- dsName
	dvals$nPermut <- permut
	
	class(dvals) <- 'phylo.d'
	return(dvals)
	
}

print.phylo.d <- function(x, ...){
    summary(x)
}

summary.phylo.d <- function(object, ...){
    cat('\nCalculation of D statistic for the phylogenetic structure of a binary variable\n')
    cat('\n  Data : ', object$dsName)
    cat('\n  Binary variable : ', object$binvar)
    cat('\n  Phylogeny : ', object$phyName)
    cat('\n  Number of permutations : ', object$nPermut)
    
    cat("\n\nEstimated D : ", object$DEstimate)
    cat("\nProbability of E(D) resulting from no (random) phylogenetic structure : ", object$Pval1)
    cat("\nProbability of E(D) resulting from Brownian phylogenetic structure    : ", object$Pval0)
    cat("\n\n")
}