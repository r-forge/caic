\name{node2tip}
\alias{node2tip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate node depth in a phylogeny }
\description{
This function calculates node depth as the maximum number of nodes between each node in a phylogeny and the tips. 
}
\usage{
node2tip(phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{A 'phylo' object.}
}
\details{
The tips are counted and so all tips have a node to tip distance of 1.
}
\value{
A vector of the node to tip counts, named using the node numbers from the phylogeny edge matrix.
}

\author{David Orme}
\examples{
data(SmallTree)
node2tip(diTree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
