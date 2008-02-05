`caic` <-
function(formula, data, phy, names.col, contr.type = "Crunch",
                  stand.contr = TRUE, ref.var=NULL, node.depth=NULL,
                  crunch.brlen=NULL, equal.branch.length=FALSE)
{


        # OLD2NEW STATUS: Converted to use new ape phylo structure
        # All geographic based contrast calculation code has been removed in order to hasten a first release of the package
        # and the caic wrapper function to contrCalc has been split to explicit caic and macrocaic wrappers 
        
        # Program Flow:
        #   1) setup - check arguments, 
        #   2) get the union set of tips and data rows and reduce the phylogeny and data to this union
        #   3) use model functions to get design and response matrices, including all NA data
        #   4) feed the model matrices into a function to calculate nodal values and contrasts
        #   5) feed the returned contrast versions of the design and response matrices into lm.fit
        
        # TODO - return node age/height
        # TODO - allow caic to be used as a contrast calculator
        # TODO - farm out common data setup in caic and macrocaic to a single function.
        
        # CHECKS AND SETUP

        # record call and dataset names
        cl <- match.call()
        phyName <- deparse(substitute(phy))
        dataName <- deparse(substitute(data))

        # check inputs are what they should be TODO FILL THIS OUT
        # if(! is.formula(formula)) stop("'formula' must be an object of class 'formula'.")
        if(! is.data.frame(data)) stop("'data' must be an object of class 'data.frame'.")
        

        # match arguments to acceptable values
        match.arg(contr.type, c("Crunch", "Brunch")) 
        
        # if not imposed, set intermediate branch length to use at polytomies - CAIC used 1 but this is not fixed
        if(is.null(crunch.brlen)) crunch.brlen <- 1
        
        # check node.depth
        if(! is.null(node.depth)){
            if(node.depth%%1 != 0 || node.depth < 1) stop("node.depth must be a positive integer greater than 1.")
        }
                
        # check the phylogeny is a rooted phylogeny and set branch lengths...
        if(! inherits(phy, "phylo")) 
            stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
        if(! is.rooted(phy))
            stop("'", deparse(substitute(phy)), "' is not rooted.")
        
        if(as.logical(equal.branch.length)) {# doesn't get evaluated if FALSE or zero
            phy$edge.length <- rep(2, length(phy$edge.length))
        }
        
        # identify the name column and make sure it is of mode character
        names.col <- deparse(substitute(names.col))
        namesInd  <- match(names.col, names(data))
        
        if(is.na(namesInd)) {
            stop("Names column '",  names.col, "' not found in data frame '", dataName, "'")
        }
        data[,namesInd] <- as.character(data[,namesInd])
        

    # DATA MATCHING AND REDUCTION
        # store original dataset size
        origTips <- with(phy, max(edge) - Nnode)
        origData <- dim(data)[1]

        # find the intersection between tip.labels and names in data frame
        in.both <- intersect(data[,namesInd], phy$tip.label)
        if(length(in.both) < 2) stop("Fewer than two tips are common to the dataset and phylogeny")

        # Label the internal nodes by their node number in the original tree to provide a backreference
        phy$node.label <- with(phy, ((max(edge)-Nnode) +1):max(edge)) 

        # i >> ditch rows with no tip
        row.in.tree <- match(data[,namesInd], in.both)
        data <- subset(data, !is.na(row.in.tree))
         
        # ii >> ditch tips which have no rows.
        tip.in.data <-  match(phy$tip.label, in.both)
        to.drop <- phy$tip.label[is.na(tip.in.data)]
        #  get subset of phylogeny to be used
        if(length(to.drop) > 0) analysisPhy <- drop.tip(phy, to.drop) else analysisPhy <- phy
        
        # useful info...
        root <- with(analysisPhy, (max(edge) - Nnode) + 1)
        

         # get the data into the same order as the tips
        tip.order <- match(analysisPhy$tip.label, data[,namesInd])
        if(any(is.na(tip.order))) stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
        data <- data[tip.order,]
        # Label the data frame rows by tip number to allow the tree to be traversed
        rownames(data) <- 1:dim(data)[1]
        
        # Size of conjunction of tree and dataset
        unionData <- dim(data)[1] 

    # GET THE BASIC MODEL MATRICES
        
        # drop any intercept from the formula
        formula <- update(formula, . ~ . - 1) # no effect if the interecept is already omitted         
        # Get the model frame including missing data
        # and check the number of complete cases in the model frame
        initMf <- model.frame(formula, data, na.action="na.pass")
        initMfComplete <- complete.cases(initMf)
        # TODO - think whether this check is always sufficient...
        if(sum(initMfComplete) < 2 ) stop("Fewer than two taxa contain complete data for this analysis")
        
    # CALCULATE MODEL 
    # GET THE MODEL MATRIX and Model Response
        
        # now we have the union of the phylogeny and data
        # get the model frame, matrix and response
        # these show the values at the tips for each term 
        mf <- model.frame(formula, data, na.action=na.pass)
        

        # HANDLE CATEGORICAL VARIABLES:
        # find the factors
        varClass <- attributes(attributes(mf)$terms)$dataClasses
        termFactors <- attributes(attributes(mf)$terms)$factors
    
        factorCols <- names(varClass)[varClass %in% c("ordered","factor")]

        # currently, want to exclude the possibility of interactions in a factor
        # in Brunch because I don't have a clue how that should be handled
        # and also need a vector showing which terms are categorical and numeric
        # in order to allow correct standardization of contrasts 
        if(contr.type=="Brunch"){
            if(any(varClass %in% c("ordered","factor") & rowSums(termFactors) > 1)){
                stop("Interactions using categorical variables not supported in Brunch analyses")}
            termClass <- apply(termFactors,2,function(X) unique(varClass[as.logical(X)]))
         } else {
             # there could be no terms with a model of the form ~ 1.
             if(length(termFactors) == 0){
                 termClass <- numeric(0)
             } else {
                 termClass <- rep("numeric", dim(termFactors)[2]) # all variables treated as numeric...
             }             
         }
     
        for(fact in factorCols){
            # - check whether all factors are ordered or binary
            currFact <- with(mf, get(fact))
            lev <-  levels(currFact)
            ord <- is.ordered(currFact)
            if(length(lev) > 2 & ! ord) stop("Unordered non-binary factor included in model formula.")
             # - In both crunch and brunch, all factors become numeric but are
            #   then treated differently by the contrast calculation engine. 
            # - modify the model frame object to make the factors numeric
            # - quote the names of the variables to assign to
             eval(parse(text=paste("mf$'", fact, "'<- as.numeric(currFact)", sep=""))) 
             attr(mf, "dataClasses") <- rep("numeric", dim(termFactors)[2])
        }
            
        # get the design matrix
        md <- model.matrix(formula, mf) 

        # sort out the reference variable
        if(! deparse(substitute(ref.var)) == "NULL"){ # this is a bit ugly but can't think of the canonical method right now...
            ref.var <- deparse(substitute(ref.var))
            if(is.na(match(ref.var, colnames(md)))) stop("Reference variable not found in design matrix")
        } else {
            ref.var <- colnames(md)[1]
        }

        # MODEL RESPONSE
        mr <- model.response(mf)
        # turn into a column matrix
        mr <- as.matrix(mr)
        colnames(mr) <- as.character(formula[2])
        # now that we have the model response for CAIC style contrasts we can substitute the reference variable
        # for empty models (i.e. models specified as X~1)
        if(is.empty.model(formula)) ref.var <- colnames(mr)
        # add to the design matrix
        md <- cbind(mr, md)

    # NOW SETUP TO GET CONTRASTS AND NODAL VALUES
        # We know the tip values, the analysis tree         
        contr <- contrCalc(md, analysisPhy, ref.var, contr.type, crunch.brlen)

    # GET RESPONSE MATRIX
        # first column of contrasts is response
        mrC <- contr$contr[,1,drop=FALSE]
        mdC <- contr$contr[,-1,drop=FALSE]

    # standardize the contrasts if required 
    # (never standardize categorical variables in Brunch)
    
        if(stand.contr){
            notCateg <- ! termClass %in% c("factor","ordered")
            mdC[,notCateg] <- mdC[,notCateg, drop=FALSE]/sqrt(contr$var.contr)
            mrC <- mrC/sqrt(contr$var.contr)
        }
    
    # FEED THE RETURNED DESIGN AND RESPONSE MATRICES INTO THE MODELLING FUNCTIONS
        # assemble the data into a finished contrast object

        ContrObj <- list()
        ContrObj$contr$response <- mrC[,,drop=FALSE]
        ContrObj$contr$explanatory <- mdC
        ContrObj$nodalVals$response <- contr$nodVal[as.numeric(rownames(contr$nodVal)) >= root,1,drop=FALSE]
        ContrObj$nodalVals$explanatory <- contr$nodVal[as.numeric(rownames(contr$nodVal)) >= root,-1,drop=FALSE]            
        ContrObj$contrVar <- contr$var.contr
        ContrObj$nChild <- contr$nChild
        
        
        # get the node depth using the tree for which we have complete data
        # and then match those nodes up against the analysis tree
        tipsWithNAdata <- data[,names.col][! complete.cases(mf)]
        if(length(tipsWithNAdata) > 0){
        	compPhy <- drop.tip(analysisPhy, tipsWithNAdata) 
    	}	else { compPhy <- analysisPhy }
    	

		nd <- node2tip(compPhy) # slow because it uses clade matrix
		names(nd) <- with( compPhy, c(tip.label, node.label))
        nd <- nd[match(analysisPhy$node.label, names(nd))]
        ContrObj$nodeDepth <- nd
        
        # gather the row ids of NA nodes to drop from the model
        validNodes <- with(ContrObj$contr, complete.cases(explanatory) & complete.cases(response))
       
        # enforce any node depth requirement
        if(! is.null(node.depth)){
            validNodes[ContrObj$nodeDepth > node.depth] <- FALSE
        }
                 
        # save for the user
        ContrObj$validNodes <- validNodes
        
        # feed the contr.model.response and contr.model.matrix
        # into lm.fit to get the model and then set up the lm object
        # - need to use lm.fit here rather than calling the model on 
        #   data=contrData because any functions in the formula are now
        #   set in the column names - don't want lm to try and reinterpret
        #   them in parsing the formula.
        # - the problem then becomes how to get the model to refer to the dataset
        mod <- with(ContrObj$contr, lm.fit(explanatory[validNodes,,drop=FALSE], response[validNodes,,drop=FALSE]))
       class(mod) <-  "lm"
       
       # assemble the output
       # return fitted model and contrasts
       RET <- list(contrast.data=ContrObj, phy=analysisPhy, mod=mod)
       class(RET) <- c("caic")
       
        # convert the ContrObj into a data frame...
       contrData <- caic.table(RET, valid=TRUE)
       RET$mod$call <- substitute(lm(FORM, data=contrData, subset=validNodes), list(FORM=formula))
       RET$mod$terms <- attr(mf, "terms")
       
       # put the model.frame in to the lm object so that predict, etc. calls work
       RET$mod$model <- contrData
       attr(RET$mod$model, "terms") <- attr(mf, "terms")
       
       # add some attributes
       attr(RET, "origTips") <- origTips
       attr(RET, "origData") <- origData
       attr(RET, "unionData") <- unionData
       attr(RET, "phyName") <- phyName 
       attr(RET, "dataName") <- dataName
       attr(RET, "contr.method") <- "caic"
       attr(RET, "contr.type") <- contr.type
       attr(RET, "stand.contr") <- stand.contr
       
       return(RET)

}