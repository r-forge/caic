\name{vcv.phylo.array}
\alias{vcv.phylo.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a 3D variance-covariance matrix for a phylogeny}
\description{
The function turns a phylogeny into a 3d array similar to a VCV matrix but retaining the individual branch lengths using a third dimension to the matrix. This is useful for handling some branch length transformations, such as kappa, in functions where VCVs are used to handle phylogenetic structure, and has a lower overhead than repeatedly calling \code{vcv.phylo} on a phylogeny that has has the edge lengths transformed.}
\usage{
vcv.phylo.array(phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{An object of class 'phylo'.}
}
\value{
A 3 dimensional array of class 'vcv.phylo.array' with dimensions of the number of taxa in the phylogeny for the rows and columns and the maximum number of branches on root to tip path as the depth. The rows and columns are named using the tip labels from the phylogeny and the depth is not named. A pair of row and column coordinates will extract a vector of all branch lengths shared by a pair of taxa, padded with NA in order to accomodate the longest such vector. The apply function can be easily used to collapse the array down to a standard VCV matrix, as in the example.
}
\author{David Orme}
\seealso{ \code{\link{vcv.phylo}}, \code{\link{lik.kappa}},}
\examples{
tree <- rcoal(8)
tree.VCV <- vcv.phylo(tree)
tree.VCVA <- vcv.phylo.array(tree)

# reconstruct a simple VCV array
tree.VCVA.reduced <- apply(tree.VCVA, c(1,2), sum, na.rm=TRUE)

all(tree.VCVA.reduced == tree.VCV)

# a kappa transformation of 0.5
apply(tree.VCVA ^ 0.5, c(1,2), sum, na.rm=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ utilities }% __ONLY ONE__ keyword per line
