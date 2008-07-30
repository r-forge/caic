\name{clade.members}
\alias{clade.members}
\alias{all.clades}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify tips descended from a node}
\description{
Obtains a vector of the tips subtending from either one node or all nodes in a phylogeny.}
\usage{
clade.members(x, phy, tip.labels = FALSE)
all.clades(phy, tips = FALSE, tip.labels = FALSE)
}

\arguments{
  \item{x}{A integer in the range minus number of internal of nodes to positive number of tips. }
  \item{phy}{An object of class 'phylo'. }
  \item{tips}{ A logical indicating whether to include external node membership in the list.}
  \item{tip.labels}{A logical flag indicating whether to return the node numbers of the tips or their tip labels.}
}
\details{
The function \code{all.clades} runs \code{clade.members} over each node in the phylogeny, possibly including the external nodes as indicated by the tips argument, and returns a list of vectors showing the members of the clade defined by each node.
}

\value{
A numeric vector of external node (i.e. tip) numbers or a character vector of tip labels for a single internal node or, for \code{all.clades}, a list of such vector for all nodes in the phylogeny.
}
\author{David Orme}


\seealso{} %% \code{\link{clade.stats}} }
\examples{
data(SmallTree)
clade.members(22, diTree, tip.labels=TRUE)
all.clades(diTree, tip.labels=FALSE)
}
\keyword{ manip }% at least one, from doc/KEYWORDS
\keyword{ utilities}% M.MassM.MassONLY ONEM.MassM.Mass keyword per line
