\name{growTree}
\alias{growTree}
\alias{linToApe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tree simulation with traits.  }
\description{
This function is still in development but is intended to provide a very general environment in which to simulate trees. The basic philosophy is that the user provides a series of expressions that define speciation rates, extinction rates and trait evolution. These expressions can make use of information about the internal state of the tree, allowing for very flexible definitions of rules for tree growth.}

\usage{
growTree(b = 1, d = 0, halt = 20, grain = 0.1, linObj = NULL, ct.start = NULL, ct.change = NULL, ct.var = NULL, dt.rates = NULL, inheritance = NULL, trace = FALSE, output.phylo = TRUE, neg.rates = "abort", inf.rates = "abort", stall.time = 10, extend.proportion=0)
linToApe(linObj)
}

\arguments{
  \item{b}{A speciation rate. This can be a numeric constant, as in the default, which specifies a single speciation rate for the simulation. Alternatively, this can be an expression, or a list of expressions which define speciation rate in terms of the properties of the tree. See details for discussion of those properties.}
  \item{d}{An extinction rate, described as above.}
  \item{halt}{A rule use to halt the simulation. The default is the number of tips in the simulation, specified as a single integer, but this can also be an expression or list of expressions on the properties of the tree. The simulation is halted when any of these expressions becomes true.}
  \item{grain}{Where rates depend on time or trait values, it becomes necessary to allow time to pass discretely in order to re-evaluate waiting times under the changing values. This sets the amount of time that is allowed to pass before re-evaluation. If rates do not depend on such changing parameters, it is sensible to set this to infinity - this will ensure that the flow of the simulation is not slowed by checking.}
  \item{linObj}{This can be used to supply an existing simulation object, which will then continue to grow under the provided rules. This allows the user to simulate trees with different sets of rules operating in different epochs. The function \code{linToApe} will convert such an object to a 'phylo' object, retaining additional trait data as extra components of the 'phylo' object list. }
  \item{ct.start}{A numeric vector specifying the starting values for continuous traits. If unnamed these will be sequentially named as 'ct1', 'ct2' etc. Knowing the names of traits becomes important if their values influence tree growth rules.}
  \item{ct.change}{A numeric vector describing the mean change per unit time in continuous trait values.}
  \item{ct.var}{A square matrix describing the covariance amongst the continuous traits.}
  \item{dt.rates}{A list of matrices describing the rate of transition between discrete character traits. Each matrix defines a trait and, as with ct.start, the list names are used to identify the traits in the simulation and default to 'dt1', 'dt2', etc. The dimnames of the matrix are used to identify the states of the trait and default to 'st1', 'st2', etc. The matrix need not be symmetrical: the rates are defined from the states in the columns to the states in the rows, hence the diagonal should probably be zero.}
  \item{inheritance}{A list of rules that are applied after a speciation and can be used to modify trait values for the descendent lineages. The names of the list specify which traits are to be modified and, for each trait specified, should return a vector of length two which replaces the existing values. }
  \item{trace}{A logical flag, indicating whether or not to provide detailed trace information on the progress of the simulation. NOT FULLY IMPLEMENTED}
  \item{output.phylo}{A logical flag indicating whether to convert the internal lineages object to a 'phylo' object.}
  \item{neg.rates}{One of 'abort', 'warn' or 'quiet', defining the behaviour when a rate calculation produces a negative number. With 'warn'  and 'quiet', negative rates are set to zero and the simulation continues.}
  \item{inf.rates}{One of 'abort', 'warn' or 'quiet', as for\code{neg.rates}. With 'warn' or 'quiet', infinite rates are left in place, resulting in events happening instantly. This may, in some cases, be desirable!}
  \item{stall.time}{If the all rates within the simulation are zero then only this length of time is allowed to pass before the simulation exits with a 'stalled' status. If \code{grain} is infinite, then the simulation stalls immediately when all rates are zero.}
  \item{extend.proportion}{This option allows the simulation to continue running for a given proportion of the time to the next speciation. This makes sense when growing a clade to a given number of extant taxa; with the default setting of zero, the resulting tree ends at a bifurcation with zero branch lengths and this option allows the tree to grow (and taxa to go extinct and traits to evolve).}
}
\details{
The main idea behind this function (which is still heavily in development) is to provide a flexible framework for simulating tree growth and trait evolution. The user provides expressions for the main arguments (\code{b}, \code{d} and \code{halt}) which act as rules defining speciation and extinction and the ending of the simulation. These can be simple constants, but can also make use of the properties of the environment of the evolving tree. This includes both lineage specific properties (as described in the \code{lineages} section of the returned value) or properties of the clade as a whole (as described in the \code{clade} section of the returned value). For example, a extinction rate might increase with lineage age (\code{d=0.01*lin.age}) or a speciation rate might decrease according to a density dependent process (\code{b=1 - (nExtantTip/500)}). Halt expressions will typically use clade properties (\code{halt=clade.age >= 5} or \code{nTips >= 50}) but could use lineage properties, for example stopping when a trait value hits a certain value (\code{halt=any(ct1 >= 10)}). It is currently not an explicit error but it is very unwise to use '==' in a halt function because it will allow the simulation to run away if the actual value steps over the test value.

Discrete traits are defined using a matrices of rates for transitions between states for each trait. At present, these are fixed for the duration of a simulation epoch and cannot be set as expressions of tree variables.

Continuous trait evolution currently employs a simple Brownian model, given a starting value and variance per unit time. The traits can have defined co-variance (the simulation uses mvrnorm at present) and can also have a defined mean change, allowing for a directional walk in the trait values. At present, it is not possible for the trait variance to vary according to the internal state of the tree; continuous characters retain the same variance and covariance for the whole of the simulation epoch.

Whilst none of the \code{halt} rules are TRUE, then the function evaluates the birth, death and discrete trait rates and converts these to waiting times using random variates from a exponential distribution with the calculated rates. These competing waiting times are compared both to each other and the \code{grain} of the simulation, the shortest waiting time is found and the relevant event is then triggered. The winning event is identified in the character vector \code{winnerName} in order to allow inheritance rules to differentiate events.
}
\value{
Depending on the value of \code{output.phylo}, either an object of class 'phylo' or an object of class 'growTree' with the following structure:
  \item{lineages}{A data frame with a row for each lineage in the tree. Each row identifies the \code{parent.id} and \code{id} of the row along with the total age of the lineage (\code{lin.age}) and the times at which the lineage was born (\code{birth.time}). If the species went extinct (or speciated) then the \code{death.time} is recorded and \code{extinct} is set to TRUE. Speciating lineages have  \code{tip} set to FALSE. Each row also records the \code{caic.code} of the lineage - this is used as a sorting code for conversion to a 'phylo' object and is a kludge. If traits are defined in the simulation then the values or states are recorded in this table. These are the current values for extant tips and the values at extinction for extinct tips and internal nodes.}
  \item{clade}{A list containing:\code{clade.age}, the total age of the simulation; \code{nLin}, the total number of lineages; \code{nTip}, the total number of tips, differentated into \code{nExtantTip} and \code{nExtinctTip}.}
  \item{rules}{A list reporting the birth (\code{b}), death (\code{d}) and stopping (\code{halt}) rules and any inheritance rules.}
  \item{ct.set}{If continuous characters were simulated, a list of the \code{ct.start}, \code{ct.change} and \code{ct.var} details provided.}
  \item{dt.rates}{If discrete characters were simulated, a list containing the \code{dt.rates} details provided.}
}
\references{ ~put references to the literature/web site here ~ }
\author{David Orme, drawing heavily on discussions with Paul-Michael Agapow.}
\examples{

# Note that some of these example could be slow, particularly if the clades speciates early.
# Alternatively, nothing may happen, returning the a tree with a single edge!

# A basic 200 tip tree, output as a 'phylo' object
    tree <- growTree(halt=200, grain=Inf)
    plot(tree)

# A basic tree of age 4 time units, output as a 'phylo' object
    tree <- growTree(halt=expression(clade.age >= 4), grain=Inf)
    plot(tree)

# A simulation in which both speciation and extinction rates increase with lineage age
# The grain is set here to allow the rates to vary with lineage age
    bUser <- expression(lin.age * 0.1)
    dUser <- expression(lin.age * 0.01)
    haltUser <- expression(clade.age >= 20)
    tree <- growTree(b=bUser, d=dUser, halt=haltUser, grain=0.01)
    plot(tree)

# Discrete character evolution between three states
    dtR <- list(matrix(c(0,0.4,0,0.2,0,0.4,0,0.2,0), ncol=3, byrow=TRUE))
    tree <- growTree(halt=100, grain=Inf, dt.rates=dtR)
    plot(tree, edge.col=unclass(tree$dt.data$dt1))

# Continuous character evolution
    tree <- growTree(halt=100, ct.start=5, ct.change=0, ct.var=0.1)

    # colour code on a tree
    cols <- as.character(cut(tree$ct.data$ct1,20, labels=rainbow(20, start=0, end=4/6)))
    plot(tree, edge.col=cols)

# Both continuous and discrete together
    tree <- growTree(halt=100, grain=Inf, ct.start=5, ct.change=0, ct.var=0.1, dt.rates=dtR)

    # colour code continuous trait and code discrete trait as edge width
    cols <- as.character(cut(tree$ct.data$ct1,20, labels=rainbow(20, start=0, end=4/6)))
    plot(tree, edge.col=cols, edge.width= 1.5 * unclass(tree$dt.data$dt1))

# Multiple expressions and named states
    dtR <- list(breedSyst=matrix(c(0,0.2,0.1,0), ncol=2, byrow=TRUE, 
                dimnames=list(c("ovip", "vivip"), c("ovip", "vivip"))))

    tree <- growTree(b=list(expression(ifelse(breedSyst == "ovip", 0.5, 0.25))),
                     d=list(expression(ifelse(breedSyst == "ovip", 0.01, 0.05))),
                     halt=list(expression(nExtantTip >= 200)), dt.rates=dtR)

    # colour code tree
    plot(tree, edge.col= unclass(tree$dt.data$breedSyst), show.tip.label=FALSE)

# Epochs... 
    tree <- growTree(halt=100, output.phylo=FALSE, grain=Inf)

    # currently there is no function for a mass extinction 
    # but here we lose extant tips with an 80% probability
    tipRow <- which(tree$lineages$tip)
    dieOff <- tipRow[runif(length(tipRow)) >= 0.2]
    tree$lineages$extinct[dieOff] <- TRUE
    tree$lineages$death.time[dieOff] <- tree$clade$clade.age

    # let the survivors repopulate
    tree <- growTree(linObj=tree, halt=200, grain=Inf)
    plot(tree)
    
# inheritance -  two regions but speciation occurs in only one region and at different rates

    # probabilities of dispersing between regions or to no regions 
    regions <- list(region = matrix(c(0,0.05,0.05,0.0001,0.1,0,0,0.01,0.1,0,0,0.01,0,0,0,0), 
                        ncol=4, dimnames=list(c("AB", "A", "B", "0"), c("AB", "A", "B", "0"))))

    # probability of speciating in each region
    spec <- list(regA_spec = expression(ifelse(region == "AB" | region == "A", 1, 0)),
                 regB_spec = expression(ifelse(region == "AB" | region == "B", 0.5, 0)))

    # inheritance rules - see which region speciated and one of the daughters inherits only that region
    inherit <- list(region = expression(if(winnerName=="regB_spec" && region[1] == "AB") c("AB","A") else region),
                    region = expression(if(winnerName=="regA_spec" && region[1] == "AB") c("AB","B") else region))

    # extinction rule - need to enforce extinction when the species is in neither region.
    extinct <- expression(ifelse(region == "0", Inf, 0))

    # grow the tree - we are deliberately setting infinite rates 
    tree <- growTree(b=spec, d=extinct, halt=300, inherit=inherit, dt.rates=regions, inf="quiet")

    # colour code tree
    plot(tree, edge.col= unclass(tree$dt.data$region), show.tip.label=FALSE)
    

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
