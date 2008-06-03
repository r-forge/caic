\name{piclm}
\alias{piclm}
\alias{piclmCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparative analysis by independent contrasts using the 'piclm' algorithm.}
\description{
Calculate a linear model using phylogenetic independent contrasts using a formula interface.
}
\usage{
piclm(formula, data, phy, names.col, stand.contr = TRUE, ref.var = NULL, node.depth = NULL,  equal.branch.length = FALSE, factor.action="abort")
}


\arguments{
  \item{formula}{A model formula.}
  \item{data}{A data frame containing the variables specified in the model formula. The data frame can contain missing data for explanatory variables.}
  \item{phy}{An object of class 'phylo' giving the phylogenetic structure in the data.}
  \item{names.col}{A name specifying the column in 'data' that matches rows to tips in 'phy'.}
  \item{stand.contr}{A logical flag indicating whether or not to standardize contrasts}.
  \item{ref.var}{A reference variable present in the model that is used to specify the direction of calculation of contrasts. If null, this is assumed to be the first explanatory variable.}
  \item{node.depth}{A positive integer greater than 1 used to restrict the model to contrasts with a node depth less than or equal to the specified depth. Tips have a depth of 1.}
  \item{equal.branch.length}{If set to 'TRUE' then all branch lengths are set to 2.}
 \item{factor.action}{One of "abort", "warn" or "allow", describing whether to stop if the formula contains a factor ("abort"), or continue after converting the factor to a numeric variable, either with ("warn") or without ("allow") a warning.}
}

\details{

Not functional or complete - will use the method implemented in the 'pendek' function piclm()...



}
\value{
  A list of class 'caic' containing the following:
	\item{contrast.data}{ A list of the following:
		\item{contr}{A list containing matrices of the contrasts in the response variable (contr\$response) and explanatory variables (contr\$explanatory).}
		\item{nodalVals}{A list containing matrices of the nodal values in the response variable (contr\$response) and explanatory variables (contr\$explanatory).}
		\item{contrVar}{A numeric vector of the calculated variances for each contrast.}
		\item{nChild}{A vector showing the number of nodes descending from each internal node}
		\item{nodeDepth}{A vector showing the maximum number of nodes between each internal node and the tips of the phylogeny (including both the node in question and the tip and hence always >=2) }
		\item{validNodes}{A logical vector showing which internal nodes on the tree have valid contrasts, given the available data and any user constraints.}

  }
  \item{phy}{A 'phylo' object describing the phylogeny used to calculate contrasts. This is the phylogeny in the original call reduced to those tips that are also named in the provided data but possibly including tips with missing data.}
  \item{lm}{An 'lm' object containing a regression model through the origin for the calculated contrasts}
}
\references{
Felsenstein, J.  (1985).  Phylogenies and the comparative method.  Am. Nat.  125, 1-15
Pagel, M. D. (1992). A method for the analysis of comparative data.  J. theor. Biol. 156, 431-442.
Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent contrasts (CAIC): an Apple Macintosh application for analysing comparative data.  Computer Appl. Biosciences 11, 247-251.
}

\author{David Orme; Andy Purvis}

\seealso{\code{\link{caic.table}} for extracting a data table, \code{\link{summary.caic}} for a data and model summary and \code{\link{caic.diagnostics}} for model checking.}
\examples{

## no examples until code working - do not use

## data(SmallTree)
## piclmMod <- piclm(othervar ~ weight, phy=diTree, data=SmallTreeDat, names.col=tip)
## 
## summary(piclmMod)
## 
## # plot the contrasts
## caicTab <- caic.table(piclmMod)
## plot(othervar ~ weight, caicTab)
## 
## # for the actual model diagnostics
## par(mfrow=c(3,1))
## caic.diagnostics(piclmMod)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
