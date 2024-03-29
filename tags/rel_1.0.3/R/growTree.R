growTree <- function(b=1,d=0,halt=20, grain=0.1, linObj=NULL,
                     ct.start=NULL, ct.change=NULL, ct.var=NULL, dt.rates=NULL,
                     inheritance=NULL, trace=FALSE, output.phylo=TRUE, 
                     neg.rates="abort", inf.rates="abort", stall.time=10, extend.proportion=0){

    # CHANGES 0.2 to 0.3 - use a data.frame rather than a list for lineages. 
    # TODO - dt.rates can be modified based on lineage properties, but there is
    # currently no mechanism for modifiying continuous traits based on lineage properties.
    
    # 0.3.2 - attempt to introduce inheritance rules...
    #       - cut clade properties out to their own list again?
    #       - multiple speciation/extinction rules ?
    
	# rexp can't handle rates of zero or Inf without a warning so handle that in producing waiting times
    waitTime <- function(rates) {
        zero <- rates == 0
        inf <- rates == Inf
        rates[zero | inf] <- 1
        wait <- rexp(length(rates), rates)
        wait[zero] <- Inf
        wait[inf] <- 0
        return(wait)
    }
	
	# The simulation has the following clade level properties:
	# - clade.age, nLin, nTip, nExtantTip, nExtinctTip
	# and the following lineage properties:
	# - id, parent.id, lin.age, birth.time, death.time, extinct, tip

    # CREATE A LINEAGE OBJECT OR USE THE PROVIDED ONE
    if(is.null(linObj)){
        
	    # IF NO EXISTING LINEAGES ARE PROVIDED, THEN INTIALIZE A NEW ONE      
    	lineages <- data.frame(parent.id=0, id=1, lin.age=0, birth.time=0, 
    	                       death.time=as.numeric(NA), extinct=FALSE, tip=TRUE, 
    	                       caic.code="", stringsAsFactors=FALSE) 
    	clade <- list(clade.age=0, nLin=1, nTip=1, nExtantTip=1, nExtinctTip=0)
        
        # CONTINUOUS TRAIT SETUP
        if(! is.null(ct.start)){
            
            require(MASS)
            ctFlag <- TRUE
        
            # check we have the right sort of info to pass to 
            if(! is.numeric(ct.start)) stop("'ct.start' must be a numeric vector providing the starting values for continuous traits.")
            if(! is.numeric(ct.change)) stop("'ct.change' must be a numeric vector providing the mean change per unit time for continuous traits.")
            if(! is.numeric(ct.var)) stop("'ct.var' must be a numeric vector  or matrix providing the (co-)variances for continuous traits.")
        
            if( length(ct.start) != length(ct.change)) stop("ct.start and ct.change must be the same length")
        
            if( is.matrix(ct.var)){
                if(dim(ct.var)[1] != dim(ct.var)[2]) stop("ct.var must be a square matrix")
                if(dim(ct.var)[1] != length(ct.start)) stop("The dimensions of ct.var must match the length of ct.start")
            } else {
                if( length(ct.start) != length(ct.var)) stop("ct.start and ct.var must be the same length")
                ct.varMat <- array(0, dim=rep(length(ct.var),2))
                diag(ct.varMat) <- ct.var
                ct.var <- ct.varMat
            }
        
            if(is.null(names(ct.start))) names(ct.start) <- paste("ct", seq(along=ct.start), sep="")
            
            # add the traits to the lineages
            lineages <- cbind(lineages, t(ct.start))
        
        } else {
            ctFlag <- FALSE
        }
        
        # DISCRETE TRAIT SETUP
        if(! is.null(dt.rates)){
            # check that the transition matrix matches the number of levels
            if( ! is.list(dt.rates) || ! all(sapply(dt.rates, is.matrix))) {
                stop("dt.rates must be a list of matrices giving the rates of change between states")}
            rateDims <- sapply(dt.rates, function(X) dim(X))
            if( ! all( rateDims[1,] - rateDims[2,] == 0)) {
                stop("dt.rates must be a list of _square_ matrices giving the rates of change between states")}
            if(any(unlist(dt.rates)) < 0) stop("Negative values in dt.rates matrix")
            dtFlag <- TRUE
            
            # name the traits
            if(is.null(names(dt.rates))) names(dt.rates) <- paste("dt", seq(along=dt.rates), sep="")
            # name the states
            dt.rates <- lapply(dt.rates,  function(X) {
                                                if(is.null(dimnames(X))){
                                                    dimnmX <- paste("st", seq(to=dim(X)[1]),sep="")
                                                    dimnames(X) <- list(dimnmX, dimnmX)
                                                }
                                                return(X)
                                                })
                                                
            # add the traits to the lineages
            # - get a data.frame of one row with first state in each matrix 
            #   as a factor with levels set by the trait states...
            dt.start <- lapply(dt.rates, function(X) factor(dimnames(X)[[1]][1], levels = dimnames(X)[[1]]))
            lineages <- cbind(lineages, dt.start)

        } else {
            dtFlag <- FALSE
        }
        
    } else {
        # TODO - some checking to make sure that possible traits are catered for
        # or assume that they are with the ct slot on the lineages object
        if(! inherits(linObj, "growTree")) stop("Lineage object is not of class 'growTree'.")
        clade <- linObj$clade
        lineages <- linObj$lineage
        if(is.null(linObj$ct.set)) ctFlag <- FALSE else ctFlag <- TRUE
        if(is.null(linObj$dt.set)) dtFlag <- FALSE else dtFlag <- TRUE
    }
    
    
	# CHECKING RATE INPUTS
	# - can be a numeric constant, an expression, a list of expressions
	#   or a named object containing one of those things...
	
	# - expressions can only use the available properties of lineages and clades
    validVar <- c(names(lineages), names(clade))
	validExpr <- function(X){ if(any(is.na(match(all.vars(X), validVar)))) FALSE else TRUE }
	
	# convert them all into a list of 'rules'
	# SPECIATION RULE(S)
	switch(mode(b),
    	"numeric" = if(length(b) > 1 | b < 0) stop("Speciation rate 'b' is numeric but not a positive scalar") else b <- list(as.expression(b)),
    	"expression"= if( ! validExpr(b)) stop("Speciation expression 'b' contains unknown variables.") else b <- list(b),
    	"list" =   {if(! all(sapply(b, mode) == "expression")) stop("The list 'b' must be a list of expressions")
    				if(! all(sapply(b, validExpr))) stop("One or more expressions in the list 'b' contain unknown variables")},
    				stop("Speciation rate 'b' not in a recognized format."))

    # assign names to the list
    if(is.null(names(b))) names(b) <- paste("b", seq(along=b), sep="")
	
	# EXTINCTION RULE(S)
	switch(mode(d),
    	"numeric" = if(length(d) > 1 | d < 0) stop("Extinction rate 'd' is numeric but not a positive scalar") else d <- list(as.expression(d)),
    	"expression"= if( ! validExpr(d)) stop("Extinction expression 'd' contains unknown variables.") else d <- list(d),
    	"list" =   {if(! all(sapply(d, mode) == "expression")) stop("The list 'd' must be a list of expressions")
    				if(! all(sapply(d, validExpr))) stop("One or more expressions in the list 'd' contain unknown variables")},
    				stop("Extinction rate 'd' not in a recognized format."))

    # assign names to the list
    if(is.null(names(d))) names(d) <- paste("d", seq(along=d), sep="")
 	
	# HALT RULE(S)
	# default is to convert a single numeric value into a number of extant tups
	switch(mode(halt),	 
	    "numeric"=if(length(halt) > 1 | halt <= 1) {
	                         stop("If numeric, 'halt' must be a single number greater than 1 giving the number of extant tips to create.")
                    } else { halt <- list(as.expression(substitute(nExtantTip >= XXX, list(XXX=halt))))},
     	"expression"= if( ! validExpr(halt)) stop("Stopping expression 'halt' contains unknown variables.") else halt <- list(halt),
    	"list" =   {if(! all(sapply(halt, mode) == "expression")) stop("The list 'halt' must be a list of expressions")
    				if(! all(sapply(halt, validExpr))) stop("One or more expressions in the list 'halt' contain unknown variables")},
    				stop("Stopping expression 'halt' not in a recognized format."))

    # assign names to the list
    if(is.null(names(halt))) names(halt) <- paste("halt", seq(along=halt), sep="")

	# check for behaviour on encountering negative or infinite rates
	neg.rates <- match.arg(neg.rates, c("abort", "warn", "quiet"))
	inf.rates <- match.arg(inf.rates, c("abort", "warn", "quiet"))
	
	
	# now start simulation - keep the simulation running whilst all of the halt
	# expressions are FALSE and while something is alive and , approximately, 
	# whilst anything is actually happening
	status <- "complete"; lastRealEvent <- 0
	
	# A function to get a rate for each lineage from either b or d
	# - is also used to extend the simulation after the last speciation
	ratesCheck <- function(rateExp, lineages, inf.rates, neg.rates){

        # - force multiplication by unity to extend constants across the clade
	    unitConst <- rep(1,length(lineages$id))
        rates <- lapply(rateExp, function(X) eval(X, env=c(lineages, clade)) * unitConst)
	    
        # negative rates?
        if(any(unlist(rates) < 0)){
                    switch(neg.rates,
                    abort=stop("Negative rates produced"),
                    warn=warning("Negative rates produced - setting to zero"))
        }
        
        rates <- lapply(rates, function(X) ifelse(X < 0, 0, X))
        
        # infinite rates
	    if(any(is.infinite(unlist(rates)))){
            switch(inf.rates,
                    abort=stop("Infinite rates produced"),
                    warn=warning("Infinite rates produced"))
        }
        
        # ensure extinct stay dead
	    rates <- lapply(rates, function(X, ext) ifelse(ext, 0, X), ext=lineages$extinct)
	    
	    return(rates)
	}
	
    while( ! any(haltStatus <-sapply(halt, eval, env=c(lineages, clade)))){
        
        # evaluate the birth and death rate expressions
        bRates <- ratesCheck(b, lineages, inf.rates, neg.rates)
        dRates <- ratesCheck(d, lineages, inf.rates, neg.rates)
        
        allRates <- c(unlist(bRates), unlist(dRates))
        
        # check to see if the stall criteria are met...
        if(all(allRates == 0)){
            if(grain==Inf) {
                warning("All rates are zero and grain is set to infinity giving no finite waiting times: exiting stalled simulation")
                status <- "stalled"
                break
            }
            if((clade$clade.age - lastRealEvent) > stall.time ){
                status <- "stalled"
                warning("Rates are all zero and stall.time is exceeded: exiting stalled simulation")
                break # so end the simulation
            }
        }
        
       # make sure something is alive...
       if(all(lineages$extinct)) {
            status <- "extinct"
            warning("All lineages extinct: exiting simuation")
            break # everything is dead so end the simulation
       }
              

       # turn rates into waiting times...
       bWait <- lapply(bRates, waitTime)
       dWait <- lapply(dRates, waitTime)
       
       # look at discrete trait changes
       # TODO - rethink this along the lapply lines - not sure it can be easily done...
       if(dtFlag){
           # for each trait, take the relevant column in lineages and
           # use it to sample columns from the appropriate rate matrix
           dtWait <- numeric(length(dt.rates))
           names(dtWait) <- names(dt.rates)
           dtWhich <- rbind(lineage=dtWait, state=dtWait)
           
           for(dt in names(dt.rates)){
               
               currDtRate <- dt.rates[[dt]][,lineages[,dt], drop=FALSE]
               currDtRate[,lineages$extinct] <- 0
               currDtWait <- as.matrix(apply(currDtRate,1, waitTime)) # because apply drops dimensions on the root
               dtWait[dt]  <- min(currDtWait)
               dtWhich[,dt] <- which(currDtWait == dtWait[dt], arr.ind=TRUE)[1,] 
               # TODO - think about that [1,] - removes ties but these are always likely to be between Inf so this is reasonable
               
           }
           
           # firstDT<- names(which.min(dtWait))
           # dtWait <- dtWait[firstDT]
           # dtWhich <- dtWhich[,firstDT]
           
       } else {
           dtWait <- Inf
           dtWhich <- 0
      }
       
       # ... look for the winning event...
       firstB_ID <- sapply(bWait, which.min)
       firstD_ID <- sapply(dWait, which.min)
       firstB_Time <- sapply(bWait, min)
       firstD_Time <- sapply(dWait, min)
       
       competWait <- c(firstB_Time, firstD_Time, dtWait, grain) # order meaningful here - breaks ties in this order
       competID <- c(firstB_ID, firstD_ID, dtWhich, 0) # if grain wins the race then no row will be selected...
       competType <- c(rep("Spec", length(bWait)), rep("Ext", length(dWait)), rep("Discrete", length(dtWait)), "Grain")
       competName <- c(names(b), names(d), names(dtWait), "Grain"  )
       
       winner <- which.min(competWait)
       winnerWait <- competWait[winner]
       winnerID <- competID[winner]
       winnerType <- competType[winner]
       winnerName <- competName[winner]

       # ... allow time to pass for the clade and extant lineages...
       clade$clade.age <- clade$clade.age + winnerWait
       lineages$lin.age[! lineages$extinct] <- lineages$lin.age[! lineages$extinct] + winnerWait

       # trait changes?
       if(ctFlag){
           # need some code in here
           delta <- as.matrix(mvrnorm(n=dim(lineages)[1], mu=ct.change, Sigma=ct.var) * winnerWait)
           delta[lineages$extinct,] <- 0 # the dead don't evolve
           lineages[,names(ct.start)] <- lineages[,names(ct.start)] + delta
       }
       
       #... and maybe something happens...
       if(winnerType == "Spec"){ # a birth!

           # copy parent into the daughters (incidentally inheriting any traits...)
           daughterID <- clade$nLin + (1:2)
           parent <- lineages[winnerID,]
           daughterInfo <- rbind(parent,parent)
           
           # inheritance rules go in here
           if(! is.null(inheritance)){
               
               for(x in seq(along=inheritance)){
                   daughterInfo[, names(inheritance)[x]] <- eval(inheritance[[x]], env=daughterInfo)
               }
           }
           
           daughterInfo$id <- daughterID
           daughterInfo$parent.id <- parent$id
           daughterInfo$birth.time <- clade$clade.age
           daughterInfo$lin.age <- 0
           daughterInfo$caic.code <- paste(daughterInfo$caic.code, c("A","B"), sep="")
           rownames(daughterInfo) <- daughterID
           lineages <- rbind(lineages, daughterInfo)
           
           # kill the parent
           lineages$extinct[winnerID] <- TRUE
           lineages$tip[winnerID] <- FALSE
           lineages$death.time[winnerID] <- clade$clade.age
           
           # update the clade
           clade$nLin <- clade$nLin+2
           clade$nTip <- clade$nTip+1
           clade$nExtantTip <- clade$nExtantTip+1
           # record that something happened here
           lastRealEvent <- clade$clade.age
        }
        
        if(winnerType == "Ext"){ #an extinction
            # kill the lineage
           lineages$extinct[winnerID] <- TRUE
           lineages$death.time[winnerID] <- clade$clade.age
           # update the clade
           clade$nExtantTip <- clade$nExtantTip-1
           clade$nExtinctTip <- clade$nExtinctTip+1
           # record that something happened here
           lastRealEvent <- clade$clade.age
        }
        
        if(winnerType == "Discrete"){ # a discrete trait changes

            lin <- dtWhich["lineage", winnerName]
            state <- dtWhich["state", winnerName]
            lineages[lin, winnerName] <- dimnames(dt.rates[[winnerName]])[[1]][state]
            lastRealEvent <- clade$clade.age
        }
    
    }

    if(trace) cat("Simulation halted at T=", clade$clade.age, "by halt rule '", names(halt)[haltStatus], "'\n")
    
    # create a lineage structure
	RET <- list(lineages=lineages, clade=clade, rules=list(b=b, d=d, halt=halt, inheritance=inheritance), status=status)
	if(ctFlag) RET <- c(RET, list(ct.set=list(ct.start=ct.start, ct.change=ct.change, ct.var=ct.var)))
	if(dtFlag) RET <- c(RET, list(dt.rates=dt.rates))
    class(RET) <- "growTree"
    
    # if the user wants to allow clade growth to continue for a further proportion (p)
    # of the waiting time to the next clade. Generate a waiting time (W) and then run again with
    # speciation set to zero, allowing the tree to grow to the current clade age + W*p, 
    # with the same extinction and character change rules

    if(extend.proportion > 0){
        

        
        # generate a waiting time
        bRates <- ratesCheck(b, lineages, inf.rates, neg.rates)
        bWait <- lapply(bRates, waitTime)
        timeToStop <- clade$clade.age + min(unlist(bWait)) * extend.proportion
        haltExpr <- as.expression(substitute(clade.age >= XXX, list(XXX= timeToStop)))
        
        RET <- growTree(b=0, d=d, halt= haltExpr, grain=grain, linObj=RET,
                     ct.start=ct.start, ct.change=ct.change, ct.var=ct.var, dt.rates=dt.rates,
                     inheritance=NULL, trace=FALSE, output.phylo=FALSE, 
                     neg.rates=neg.rates, inf.rates=inf.rates, stall.time=stall.time, extend.proportion=0)

        
        
    }
    
    
	if(output.phylo) RET <- linToApe(RET)

	return(RET)
}


