\name{force.binary.tree}
\alias{force.binary.tree}
\title{Ensure a phylogeny contains no polytomies}
\description{
Many phylogenies contain polytomies, but many functions require bifurcating trees. This 
function replaces each polytomy by an arbitrary series of bifurcations separated by edges 
of length \code{arb.branch} (zero by default).}
\usage{
force.binary.tree(phylo, arb.branch = 0, node.suffix="_Poly")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{An object of class \code{"phylo"} }
  \item{arb.branch}{A number to which new branches are set. Zero is the only sensible value for 
  most purposes.}
  \item{node.suffix}{The node labels of the phylogeny are amended to indicate which are resolved polytomies - this suffix is used to identify originally polytomous nodes and the shared polychotomy from which they have been split.}
}
\details{
}

\value{
An object of class \code{"phylo"} with no polytomies.
}

\author{David Orme}


\examples{
tdi <- rcoal(20)
tmult <- di2multi(tdi, quantile(tdi$edge.length, 0.1))

par(mfrow=c(1,2))

plot(tmult, show.node.label=TRUE)

# use a non-zero branch length to show resolution of polytomies 
tforce <- force.binary.tree(tmult, arb.branch=0.1)

## TODO - this crashes out with that stupid ape plot bug when running examples.
## plot(tforce, show.node.label=TRUE)


}
\keyword{ manip}
