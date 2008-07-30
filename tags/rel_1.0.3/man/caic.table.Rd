\name{caic.table}
\alias{caic.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract a data frame from a CAIC object.}
\description{
A CAIC object contains a relatively complex list containing the contrasts, phylogeny and linear model. This function extracts the contrast values to a data frame, optionally including information on all nodes in the phylogeny, not just those for which valid contrasts can be calculated. There is also an option to include nodal values for the internal nodes in the data frame. The table also includes studentized residuals from the linear model and calculated node ages, variances, number of descendent lineages and node depths.
}
\usage{
caic.table(caicObj, validNodes = TRUE, nodalValues=FALSE, ultrametric.tol=0.0001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{caicObj}{A CAIC object produced by either \code{caic()} or \code{macrocaic()}.}
  \item{validNodes}{ }
  \item{nodalValues}{ ~~Describe \code{validNodes} here~~ }
  \item{ultrametric.tol}{Passed to is.ultrametric to determine the tolerance for detecting whether the analysis tree is ultrametric.}

}
\details{
The function is also called internally by caic and macrocaic to provide a model frame for the 'lm' object. 

}
\value{
A data frame containing TODO
}
\author{David Orme}
\seealso{ \code{\link{crunch}},\code{\link{brunch}},\code{\link{macrocaic}}}
\examples{
data(SmallTree)
crunchMod <- crunch(othervar ~ weight, phy=diTree, data=SmallTreeDat, names.col=tip)
caic.table(crunchMod)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{utilities}% M.MassM.MassONLY ONEM.MassM.Mass keyword per line
