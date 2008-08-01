## CDLO - minor spacing and bracket alignment edits to allow function folding
## CDLO - comment out headers
##
## library(ape)
## library(mvtnorm) # need to download from CRAN
## library(MASS)
## 
## THISVERSION = 3.3
## 
## # ----------------------------------------------------------------------
## # This generates random data
## #
## rand.dat <- function(lambda, Vmat) {
## 	imat <- matrix(0, nrow = dim(Vmat)[1], ncol = dim(Vmat)[2] )
## 	diag(imat) <- 1	
## 	Vmat <- lambda * Vmat + (1 - lambda) * imat
## 	dat <-  rmvnorm(1, matrix( 0, nrow = dim(Vmat)[1] ), Vmat) 
## 	dat <- as.matrix(t(dat))
## 	return(dat)
## 	}
## 
## # ----------------------------------------------------------------------
## # This first section of code fits lambda to data on single traits
## #


## prune.single has been made redundant by expanding the existing order.D and order.V and adding code to lik.lambda

## prune.single <- function(x, data, V) {
## 
## 	dat <- data.frame(V1 = x, row.names=rownames(dat))
## 
## 	nms <- row.names(dat)
## 	if(length(nms) == 0) stop("Need to supply row names for the data")
## 	idx <- sort(nms, index.return = TRUE)$ix
## 	sort.nms <- sort(nms, index.return = TRUE)$x
## 	dat <- dat[idx,]
## 	Vnms <- row.names(V)
## 	if(length(Vnms) == 0) stop("Need to supply row names for the Variance matrix")
## 	idx <- sort(Vnms, index.return = TRUE)$ix
## 	V <- V[idx, idx]
## 	
## 	nms <- sort(nms)
## 	Vn <- sort(row.names(V))
## 	idx <- which(nms != Vn)
## 	
## 	
## 	if(length(idx) > 0) stop("Error, taxon names do not match")
## 	
## 	complete <- complete.cases(dat)
## 	idx <- which(complete == TRUE)
## 	V <- V[idx, idx]
## 	x <- dat[idx]
## 	sort.prune.names <- sort.nms[idx]
## 	dat <- data.frame(V1 = x, row.names = sort.prune.names)
## 
## 	return(list(dat = dat , V = V))
## }

lik.lambda <- function(x, V, data=NULL, lambda) {
	
	# rewritten CDLO 30/07/08
	# data checking and sorting
	if(is.null(data)){
	    x <- order.D(x)
	} else {
	    if(! is.character(x)) stop("If not a named vector, 'x' must be a character string identifying a column 'data'")

	    if(! x %in% names(data)) stop("Variable '", x, "' not found in dataset.")
	    x <- subset(data, select=x, drop=TRUE)
	    names(x) <- row.names(data)
	    x <- order.D(x)
	}
    
    # handle the possibility of a VCV array
    if(inherits(V, 'vcv.array')) V <- apply(V, c(1,2), FUN=sum, na.rm=TRUE) # reduces the array to a standard VCV

	V <- order.V(V)
	
    if(! all(names(x) == rownames(V))) stop("Taxon names do not match between the sorted data and the sorted V matrix.")

	# ditch missing data
	missing <- which(is.na(x))
	if(length(missing) > 0 ){
    	x <- x[-missing]
    	V <- V[-missing, -missing]
	}
	
	lamTrans <- function(V, lambda) {
		dV <- diag(V) 
		V <- V * lambda
		diag(V) <- dV
		return(V)
		}
	
	est.mean <- function(y, V) {
		iV <- solve(V, tol = .Machine$double.eps)
		xdum <- matrix(1, nrow = length(y))
		xVix <- crossprod(xdum, iV %*% xdum)
		xViy <- crossprod(xdum, iV %*% y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy # This is a bad thing!
		return(mu[1])
		}

	est.var <- function(x, V) {
		mu <- est.mean(x, V)
		iV <- solve(V, tol = .Machine$double.eps)
		e <- x - mu
		s2 <- crossprod(e, iV%*%e)
		n <- length(x) 
		return(s2 / (n - 1) )
		}

	mv.lik <- function(x, V) {
		mu <- est.mean(x, V)
		s2 <- est.var(x, V)
		n <- length(x)
		logDetV <- determinant(V, logarithm = TRUE)$modulus[1]
		ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n-1) / 2.0
		return( list(ll = ll, mu = mu, s2 = s2) )
		}

	V <- lamTrans(V, lambda)
	return(mv.lik(x, V))
}


lik.kappa <- function(x, V.arr, data=NULL, kappa) {
	
	# rewritten CDLO 30/07/08
	# data checking and sorting
	if(is.null(data)){
	    x <- order.D(x)
	} else {
	    if(! is.character(x)) stop("If not a named vector, 'x' must be a character string identifying a column 'data'")

	    if(! x %in% names(data)) stop("Variable '", x, "' not found in dataset.")
	    x <- subset(data, select=x, drop=TRUE)
	    names(x) <- row.names(data)
	    x <- order.D(x)
	}

    # needs an array form of the VCV
    if(! inherits(V.arr, 'vcv.array')) stop("Kappa transformations require an array form of the VCV matrix storing indvidual branch lengths")
	V.arr <- order.V(V.arr)
	
    if(! all(names(x) == rownames(V.arr))) stop("Taxon names do not match between the sorted data and the sorted V matrix.")

	# ditch missing data
	missing <- which(is.na(x))
	if(length(missing) > 0 ){
    	x <- x[-missing]
    	V <- V[-missing, -missing]
	}
	

	kappaTrans <- function(V.arr, kappa) {
	    
	    ## In  a move I find surprising, NA^0 = 1, so need to catch that case
	    if(kappa == 0) V.arr <- (V.arr > 0) else V.arr <-  V.arr ^ kappa # V.arr becomes logical in the case kappa = 0 but sum handles that in the next lin
		V <- apply(V.arr, c(1,2), FUN=sum, na.rm=TRUE) # reduces the array to a standard VCV
		return(V)
	}
		
	est.mean <- function(y, V) {
		iV <- solve(V, tol = .Machine$double.eps)
		xdum <- matrix(1, nrow = length(y))
		xVix <- crossprod(xdum, iV %*% xdum)
		xViy <- crossprod(xdum, iV %*% y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy # This is a bad thing!
		return(mu[1])
	}

	est.var <- function(x, V) {
		mu <- est.mean(x, V)
		iV <- solve(V, tol = .Machine$double.eps)
		e <- x - mu
		s2 <- crossprod(e, iV%*%e)
		n <- length(x) 
		return(s2 / (n - 1) )
	}

	mv.lik <- function(x, V) {
		mu <- est.mean(x, V)
		s2 <- est.var(x, V)
		n <- length(x)
		logDetV <- determinant(V, logarithm = TRUE)$modulus[1]
		ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n-1) / 2.0
		return( list(ll = ll, mu = mu, s2 = s2) )
	}

	V <- kappaTrans(V.arr, kappa) 
	return(mv.lik(x, V))
}

# renamed to avoid the profile() generic for model coefficient
profileLik.lambda <- function(x, V, data=NULL, at=seq(0, 1, by = 0.01), plot=TRUE, ...) {
    
    if(any(at < 0)  | any(is.na(at))) stop("At least one lambda values is either NA or less than zero. ")
        
    ll <- sapply(at, function(l) return(lik.lambda(x, V, data, l)$ll ) )

    if(plot) plot(at, ll, type = "l", col =  "red", ...)

    invisible(data.frame(at=at, ll=ll))
}

profileLik.kappa <- function(x, V.arr, data=NULL, at=seq(0, 1, by = 0.01), plot=TRUE, ...) {

        if(any(at < 0) | any(is.na(at))) stop("At least one kappa values is either NA or less than zero. ")

		ll <- sapply(at, function(l) return(lik.kappa(x, V.arr, data, l )$ll ) )
		
		plot(at, ll, type = "l", col =  "red", ...)

        invisible(data.frame(at=at, ll=ll))

}

profile.lambda <- function(x, dat, V){
    .Deprecated("profileLik.lambda")
    profileLik.lambda(x, dat, V)
}

profile.kappa <- function(x, V){
    .Deprecated("profileLik.kappa")
    profileLik.kappa(x, V)
}

# renamed to avoid clashes with the max() generic 
maxLik.lambda <- function(x, V, data=NULL) {
		ll <- function(lambda) return(-1* lik.lambda(x, V, data, lambda)$ll )
		optimize(ll, c(0,1)) 
}

maxLik.kappa <- function(x, V.arr, data=NULL) {
		ll <- function(kappa) return(-1* lik.kappa(x, V.arr, data, kappa)$ll )
		optimize(ll, c(0,3)) 
}		

max.lik.lambda <- function(x, V, dat){
    .Deprecated("maxLik.lambda")
    maxLik.lambda(x, V, dat)
}

max.lik.kappa <- function(x, V){
    .Deprecated("maxLik.kappa")
    maxLik.kappa(x,V)
}


lam.test.single <- function(x, data, V, pretty=TRUE) {

	prune.dat <- prune.single(x, dat, V)
	
	Vmat <- prune.dat$V
	data <- prune.dat$dat

	ml.lam <- maxLik.lambda(data[,1], Vmat, data)
	lam0 <- lik.lambda(data[,1], Vmat, data, 0)
	lam1 <- lik.lambda(data[,1], Vmat, data, 1)
	
	lrt0 <- 2 * (-ml.lam$objective - lam0$ll)
	lrt1 <- 2 * (-ml.lam$objective - lam1$ll)
	
	p0 <- 1 - pchisq(lrt0, 1)
	p1 <- 1 - pchisq(lrt1, 1)
	
	cat("____________________________\n")
	cat("\n")
	cat("ML estimate of lambda: ", ml.lam$minimum, "\n")
	cat("Maximised log-likelihood: ", -ml.lam$objective, "\n")
	cat("\n")
	cat("Likelihood at lamdba = 0: ", lam0$ll, "\n")
	cat("Test of Lambda = 0: chisq = ", lrt0, " P = ", p0, "\n")
	cat("\n")
	cat("Likelihood at lamdba = 1: ", lam1$ll, "\n")
	cat("Test of Lambda = 1: chisq = ", lrt1, " P = ", p1, "\n")
	cat("____________________________\n")
}


# This fits a GLM, correcting for phylogeny, with the option
# of setting the value of lambda, the index of phylogenetic 
# dependence (default is 1.0)

## TODO - rename phylomat to V throughout to maintain consistenct with profile functions...

pglm <- function(formula, data, phylomat, lambda = 1.0, ...) {

	prune <- function(dat, Vmat) {
	# ------ Delete from here if you want to take the risk! ---------------------	
	# Makes sure data and matrix are in the same order
		nms <- row.names(dat)
		if(length(nms) == 0) stop("Need to supply row names for the data")
		idx <- sort(nms, index.return = TRUE)$ix
		dat <- dat[idx,]
		Vnms <- row.names(Vmat)
		if(length(Vnms) == 0) stop("Need to supply row names for the Variance matrix")
		idx <- sort(Vnms, index.return = TRUE)$ix
		Vmat <- Vmat[idx, idx]
		
		nms <- sort(nms)
		Vn <- sort(row.names(Vmat))
		idx <- which(nms != Vn)

		
		if(length(idx) > 0) stop("Error, taxon names do not match")
    # ------ Delete to here if you want to take the risk! ---------------------
		
		complete <- complete.cases(dat)
		idx <- which(complete == TRUE)
		Vmat <- Vmat[idx, idx]
		dat <- dat[idx,]

		return(list(dat = dat , Vmat = Vmat))
	}
	
	Dfun <- function(Cmat) {
		iCmat <- solve(Cmat,  tol = .Machine$double.eps)
		svdCmat <- La.svd(iCmat)
		D <- svdCmat$u %*% diag(sqrt( svdCmat$d )) %*% t(svdCmat$v)
		return( t(D) )
	}
			
	lamTrans <- function(Vmat, lambda) {
		V1 <- Vmat
		diag(V1) <- 0
		V2 <- diag( diag(Vmat), ncol = length(Vmat[1,]), nrow = length(Vmat[,1]))
		Vmat <- V1 * lambda + V2
		return(Vmat)
	}
	
	
	## CDLO - commented out. This function is currently unused within pglm() and contains
	##        the mystery variable nx, which is not declared anywhere
	
	## resVar <- function(y, Vmat, p) {
	## 	iV <- solve(Vmat, tol = .Machine$double.eps)
	## 	e <- y - p
	## 	s2 <- crossprod(e, iV %*% e)
	## 	if( s2 < 0) { cat("ERROR -- negative variance\n")
	## 		s2 <- s2 * -1}
	## 	n <- length(y) 
	## 	return( s2 / (n- nx) )
	## }
		
	# Estimates the GLS parameters for given data
	get.coeffs <- function(Y, V, X) {
		iV <- solve(V, tol = .Machine$double.eps)
		xVix <- crossprod(X, iV %*% X)
		xViy <- crossprod(X, iV %*% Y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy 	#This is  a bad thing to do!!!!
		return(mu)
	}

	# Estimates the variance of a given trait (accounting for phylogeny)
	est.var <- function(y, V, x, mu ) {
		iV <- solve(V, tol = .Machine$double.eps)
		e <- y - x %*% mu
		s2 <- crossprod(e, iV %*% e)
		n <- length(y) 
		k <- length(x[1,])
		return( s2 / (n- k) )
	}
	
	# Full ML estimation for given x and V
	log.likelihood <- function(y, x, V) {
		mu <- get.coeffs(y, V, x)
		s2 <- est.var(y, V, x, mu)
		n <- length(x[,1])
		logDetV <- determinant(Vmat, logarithm = TRUE)$modulus[1]
		ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - 1)/2.0
		ypred <- x%*%mu	
		return( list(ll = ll, mu = mu, s2 = s2) )
	}
		
		
	null.var <- function(y, V) {
		X <- matrix(1, nrow = length(y))
		mu <- get.coeffs(y, V, X)
		return(est.var(y, V, X, mu))
	}

	
	Vmat <- as.matrix(phylomat)
	Vmat <- lamTrans(phylomat, lambda)
	
	prune.dat <- prune(data, Vmat)
	Vmat <- prune.dat$Vmat
	data <- prune.dat$dat
	nm <- names(data)
	
	n <- length(data[,1])
	
	# Get the design matrix
	m <- model.frame(formula, data)
	y <- m[,1]
	x <- model.matrix(formula, m)
	k <- length(x[1,])
	
	namey <- names(m)[1]
	
	ll <- log.likelihood(y, x, Vmat)
	
	log.lik <- ll$ll	

	aic <- -2 * log.lik + 2 * k
	aicc <- -2 * log.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
	
	coeffs <- ll$mu
	coeffs <- data.frame(t(coeffs))
	names(coeffs) <- colnames(x)
	varNames = names(m)

	
	pred <- x %*% ll$mu 
	
	res <- y - pred
	D <- Dfun(Vmat)
	pres <- D %*% res
	
	fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
	
	logDetV <- determinant(Vmat, logarithm = TRUE)$modulus[1]
 	
	logLikY <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(ll$s2) - logDetV / 2.0  - (n - 1 )/ 2.0
	
	RMS <- ll$s2
	RSSQ <- ll$s2 * (n - k)
	NMS <- RMS
	NSSQ <- RSSQ
	
	if(k > 0) {
		NMS <- null.var(y, Vmat)
		NSSQ <- NMS * (n - 1)
		}

	# Bits for parameter errors	
	errMat <- t(x)%*% solve(Vmat) %*% x  
	errMat <- solve(errMat) * RMS[1] 
	sterr <- diag(errMat)
	sterr <- sqrt(sterr)
	
	
	ret <- list(model = fm, formula = formula, logLikY = logLikY, RMS = RMS, NMS = NMS,
	            NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, aicc = aicc, n = n, k = k,
	            sterr = sterr, vcv = errMat, fitted = pred, residuals = res, phyres = pres,
	            x = x, data = data,  varNames = varNames, y = y, V = Vmat, lambda = lambda,
	            L0 = NULL, L1 = NULL, LamOptimised = FALSE, namey = namey)

	class(ret) <- "pglm"
	return(ret)
	
}

# This returns the coefficients from the model
## CDLO - argument name changed for consistency with S3 generic
coef.pglm <- function(object, ...) {
    ret <- object$model$coef 
	return(ret) 
}

# This returns the residuals from the model
## CDLO - argument name changed for consistency with S3 generic
residuals.pglm <- function(object, phylo = FALSE, ...) {
    ret <- NULL
	if(phylo == FALSE){ret <- object$res} else {ret <- object$phyres}
	return(ret)
}

# This returns the fitted values
## CDLO - argument name changed for consistency with S3 generic
fitted.pglm <- function(object, ...){
    ret <- object$fitted
    return(ret)
}

# This predicts for given x
## CDLO - argument name changed for consistency with S3 generic
## CDLO - argument name of x changed to discriminate from generic to plot and print

predict.pglm <- function(object, pred.x, ...) {
    mu <- as.matrix(coef(object) )
    ret <- cbind(1,  pred.x) %*% t(mu)
    return(ret)
}

# This returns the AIC
## CDLO - argument name changed for consistency with S3 generic
##      - TODO because aic is calculated automatically, k cannot be changed
##        in line with the generic method

AIC.pglm <- function(object, k, ...) {
    ret <- object$aic
    return(ret[1])
}

# This returns the AICc
## CDLO - argument name changed for consistency with S3 generic
AICc.pglm <- function(object) {
    ret <- object$aicc
    return(ret[1])
}

# This returns the value of lambda at which the pglm was evaluated
## CDLO - argument name changed for consistency with S3 generic
lambda.pglm <- function(object) {
    ret <- object$lambda
    return(ret[1])
}

# Very rough function for displaying a pglm output
## CDLO - argument name changed for consistency with S3 generic

summary.pglm <- function(object,...) {
		
		testLambda <- function(pobj) {
			
			lrt0 <- 2 * (pobj$logLikY - pobj$L0)
			lrt1 <- 2 * (pobj$logLikY - pobj$L1)
			
			p0 <- 1 - pchisq(lrt0, 1)
			p1 <- 1 - pchisq(lrt1, 1)
			
			cat("     Test of Lambda = 0: chisq = ", lrt0, " P = ", p0, "\n")
			cat("     Test of Lambda = 1: chisq = ", lrt1, " P = ", p1, "\n")
			}
			
	
	cat("\n\n--------------------------------------------------------\n")
	cat("Summary of Generalised Least Squares, correcting for \n")
	cat("Phylogeny:\n\n")
	cat("Number of parameters = ", object$k,"\n")
	cat("Number of data points = ", object$n,"\n\n")
	cat("Lambda statistic = ", object$lambda, "\n")
	if(object$LamOptimised == TRUE) { testLambda(object)}
	cat("Maximised log-likelihood = ", object$logLikY,"\n\n")
	cat("Model AIC = ", object$aic, "\n")
	cat("Model AICc = ", object$aicc, "\n\n")

	cat("Null Mean Square = ", object$NSSQ,"\n")
	cat("Residual Mean Square = ", object$RSSQ, "\n\n")
	cat("Raw R^2 = ", (object$NSSQ - object$RSSQ) / object$NSSQ, "\n")
	cat("Adjusted R^2 = ", (object$NMS - object$RMS) / object$NMS, "\n")
	
	Fstat <- ((object$NSSQ - object$RSSQ) / object$RMS) / (object$k - 1)
	
	cat("F statistic = ",  ((object$NSSQ - object$RSSQ) / object$RMS) / (object$k - 1), " ")
	cat("P model = ", pf(Fstat, object$k - 1, object$n - object$k,  ncp=0, lower.tail = FALSE, log.p = FALSE), "\n\n")
	cat("Summary of coefficients:\n\n")
	coeffs <- coef(object)
	errs <- object$sterr
	cat("Term\t\tEstimate\t\tStd Err\t\tT-value\tP\n")
	storet<-c()
	for(i in 1:length(coeffs) ) {
		est <- coeffs[1,i]
		nm <- names(coeffs)[i]
		se <- errs[i]
		Tstat <- est / se
		storet<-c(storet,Tstat)	
		Pval <- 2 * ( 1 - pt( abs(Tstat), object$n - object$k) )
		cat(nm,"\t")
		cat(est, "\t", se, "\t", Tstat, "\t", Pval, "\n")
		}
	
	cat("\n\n--------------------------------------------------------\n")
}

# Simple print function
## CDLO - argument name changed for consistency with S3 generic
print.pglm <- function(x, ...) {
	return(summary(x))
}


# Diagnostic plots
## CDLO - argument name changed for consistency with S3 generic
plot.pglm <- function(x, ...) {
	layout(matrix(c(1,2,3,4), 2, 2, byrow = FALSE))
	res <- residuals(x, phylo = TRUE)
	res <- res / sqrt(var(res))[1]
	truehist(res, xlab = "Residual value (corrected for phylogeny)")
	qqnorm(res)
	abline(0, 1)
	plot(fitted(x), res, xlab = "Fitted value", ylab = "Residual value (corrected for phylogeny)"  )
	plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value")
}

# ----------------------------------------------------------------------
# This optimises the value of lambda for a given dataset
# and phylogeny. The return is the model at the ML value of
# lambda together with the ML estimate of lambda and the maximised
# log-likelihood. Also returned are the logLikelihoods at Lambda = 0
# and Lambda = 1.0

## CDLO - phylomat argument renamed for consistenct with pglm

pglmEstLambda <- function(formula, data, phylomat, plotit = FALSE,  ...) {
	
	ll.fun <- function(lam) {
		pg <- pglm(formula, data, phylomat, lam)
		ll <- pg$logLikY
		return( ll )
	}
	
	oL <- optimize( ll.fun, interval = c(0,1), maximum = TRUE )
	L1 <- ll.fun(1)
	L0 <- ll.fun(0)
	fm <- pglm(formula, data, phylomat, oL$maximum)
	
	if(plotit == TRUE) {
		lambda <- seq(0, 1, by = 0.01)
		logLikelihood <- sapply(lambda, ll.fun)
		plot(lambda, logLikelihood, type = "l")
		maxlikelambda<-which(logLikelihood==max(logLikelihood))[1]
		sigdiff<-logLikelihood[maxlikelambda]-(.5*3.8414587)
		if(lambda[maxlikelambda]<0.01){lowerlamb<-0}
		if(lambda[maxlikelambda]>0.99){upperlamb<-1}
		if(lambda[maxlikelambda]>=0.01){
		    lowerlamb <- lambda[which(abs(logLikelihood[1:maxlikelambda] - sigdiff)
		                 == min(abs(logLikelihood[1:maxlikelambda]-sigdiff)))]
		}
		if(lambda[maxlikelambda]<=0.99){
		    upperlamb <- lambda[maxlikelambda+which(abs(logLikelihood[maxlikelambda:length(lambda)] - sigdiff)
		                 == min(abs(logLikelihood[maxlikelambda:length(lambda)]-sigdiff)))]
		}
	}
		
	
	fm$L1 <- L1
	fm$L0 <- L0
	fm$LamOptimised <- TRUE
	ret <- fm
	return( ret )
}


# ----------------------------------------------------------------------
# Performs an ANOVA on the pglm Object using Sequential Sums of Squares
#
## CDLO - argument name changed for consistency with S3 generic
anova.pglm <- function(object, ...) {
    
	V <- object$V
	dat <- object$data
	x <- data.frame( object$x )
	nm.x <- attr( terms(object$formula), "term.labels")
	k <- length(nm.x)
	nm.x <- nm.x[1:k]
	SSQ <- matrix(0, nrow = k + 1)
	remDF<- matrix(0, nrow = k + 1)
	SSQ[1] <- object$NSSQ
	remDF[1] <- object$n - 1
	for( i in 1: k ) {
			xv <- paste(nm.x[1:i], sep = "")
			fmla <- as.formula(paste( object$namey, " ~ ", paste(xv, collapse = "+") ))
			plm <- pglm(fmla, dat, V)
			SSQ[i + 1] <- plm$RSSQ
			remDF[i + 1] <- remDF[1] - plm$k + 1
	}
	errorDF <- object$n - object$k
	errorSSQ <- object$RSSQ
	errorMS <- object$RSSQ / errorDF
	termDF <- remDF[1:k] - remDF[2:(k+1)]
	termSSQ <- SSQ[1:k] - SSQ[2:(k+1)]
	MS <- termSSQ / termDF
	F <- MS / errorMS
	pF <-  1 - pf(F, termDF, errorDF)
	cat("Source\tDF\t\tSSQ\t\tMS\t\tF\n")
			for(i in 1:k) {
			cat(nm.x[i], "\t\t", termDF[i], "\t", termSSQ[i], "\t",
			MS[i], "\t", F[i])
			cat("   { P = ", pF[i], " }\n")
	}
	cat("Error\t", errorDF, "\t", errorSSQ, "\t", errorMS, "\n")
}


# ----------------------------------------------------------------------
# Functions for getting the phylo matrix and data in alphabetical order
#
order.V <- function(V) { # CDLO argument changed from Vmat to V for consistency
    
    if(! is.array(V)) stop("V is not a matrix or array")
    nmV <- dimnames(V)
    if(is.null(nmV)) stop("V does not have dimnames.")
    
    # now copes with differently ordered rows and columns
    if(! setequal(nmV[[1]], nmV[[2]])) stop("The row and column names in V do not match")
	ordR <- order(nmV[[1]])
	ordC <- order(nmV[[2]])
	
	nDim <- length(nmV)
	switch(as.character(nDim), "2" = {V <- V[ordR, ordC]},
	             "3" = {V <- V[ordR, ordC, ]},
	             stop("Unexpectedly large number of dimensions in V"))
	return(V)
	
}

order.D <- function(data) {
    
    if(is.vector(data)){
        
        if(! is.numeric(data)) stop("x must be a numeric vector")
	    if(is.null(names(data))) stop("Need to supply names for the data vector")
	    
        data <- data[order(names(data))]
        
    } else if(is.data.frame(data)) {
        
	    if(is.null(row.names(data))) stop("Need to supply row names for the data frame")

    	data <- data[order(row.names(data)),]
	} else {
	    stop("Data not a named vector or data frame.")
	}

	return(data)
}


# ------------------------------------------------------------------------



