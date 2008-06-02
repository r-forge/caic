\name{align.tips}
\alias{align.tips}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Align the tips on a phylogeny}
\description{
Creates an ultrametric phylogeny by forcing all the tips to align at a given time subsequent to the last speciation.
}
\usage{
align.tips(phy, at = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{A object of class 'phylo'.}
  \item{at}{The time at which to align the tips.}
}
\details{
If the value of \code{at} predates the last speciation in the phylogeny then the function will return an error
}
\value{
  The function returns \code{phy}, with terminal branch lengths altered so that the tips align at the time provided.
}
\author{David Orme}
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\examples{
phy <- rtree(20)
par(mfrow=c(2,1))
plotPhy <- plot(phy)
phyAlign <- align.tips(phy, at=plotPhy$x.lim[2])
plot(phyAlign)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
