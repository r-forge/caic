\name{force.binary.tree}
\alias{force.binary.tree}
\title{Ensure a phylogeny contains no polytomies}
\description{
Many phylogenies contain polytomies, but many functions require bifurcating trees. This 
function replaces each polytomy by an arbitrary series of bifurcations separated by edges 
of length \code{arb.branch} (zero by default).}
\usage{
force.binary.tree(phylo, arb.branch = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{An object of class \code{"phylo"} }
  \item{arb.branch}{A number to which new branches are set. Zero is the only sensible value for 
  most purposes.}
}
\details{
If polytomies are hard, the result can be used straightforwardly in \code{pic}. If the 
polytomies are soft, use \code{caic} with the polytomous phylogeny. The function \code{map.nodes} 
maps nodes between the polytomous and bifurcating trees.}

\value{
An object of class \code{"phylo"} with no polytomies.
}

\author{David Orme}

\section{Version }{1.0}
\section{Date }{13 January 2006}

\seealso{\code{\link{map.nodes}}}

\examples{
tdi <- rcoal(20)
tmult <- di2multi(tdi, quantile(tdi$edge.length, 0.1))


plot(tmult)

par(mfrow=c(1,2))
plot(td, cex=2, show.nod=TRUE)
plot.new()
text(rep(c(0.25,0.75), each=16), rep(seq(0.9,0.1, length=16), 2), as.numeric(td$edge), cex=3)


}
\keyword{ manip}
