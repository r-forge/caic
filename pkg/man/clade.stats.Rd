\name{clade.stats}
\alias{clade.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize data at nodes in a phylogeny.}
\description{
Calculates a summary statistic at each node of a phylogeny, given data for the tips. Data is provided for each tip in a phylogeny and the function supplied by the user is applied to subsets of that data for each node with the phylogeny.
}
\usage{
clade.stats(dataf, phyl, fun, ..., tips = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{A data frame of numeric variables with rownames showing the tips to which each row of data belongs. Non-numeric variables are dropped.}
  \item{phyl}{ An object of class 'phylo'.}
  \item{fun}{ A function to be applied to subsets of numeric data at each node in the tree.}
  \item{\dots}{ Further arguments to be passed to 'fun'. }
  \item{tips}{A logical indicating whether to return the statistics for external nodes.}
}
\details{
The function runs \code{all.clades} to get a list of tips present at each node and then applies the provided function to subsets of the data frame defined by the tips descending from each node. The row names of the data frame provided are matched against the tip labels of the phylogeny. If any tips are not represented in the data frame then the function aborts and returns a character vector of the missing tip names.
}
\value{
A data frame with a row for each internal node and possibly each external node as indicated by the 'tips' argument. The values in each row are the results of applying 'fun' to subsets of the tips for each numeric variable in the data frame.
}
\author{David Orme}

\seealso{\code{\link{all.clades}}}
\examples{
data(bird.orders)
ntip <- length(bird.orders$tip.label)
dat <- data.frame(a=rnorm(ntip), b=runif(ntip),
                  row.names=bird.orders$tip.label)
clade.stats(dat, bird.orders, mean)
}

\keyword{utilities}% at least one, from doc/KEYWORDS
