\name{pd.bootstrap}
\alias{pd.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a bootstrapped phylogenetic diversity distribution.}
\description{
Takes a phylogeny and repeatedly calculates the phylogenetic diversity (PD) of a randomly selected set of tips. The size of the sample of tips and the number of replicates are specified by the user, as is the method used to calculate PD. The probabilities of selecting each tip can be controlled by providing weights.}
\usage{
pd.bootstrap(cm, ntips, reps = 1000, method = "TBL", tip.weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cm}{A 'clade.matrix' object (or a 'phylo' object that will be converted).}
  \item{ntips}{A single integer giving the number of tips to be selected.}
  \item{reps}{The number of replicate values to calculate.}
  \item{method}{The method used to calculate PD (see \code{\link{pd.calc}}).}
  \item{tip.weights}{A numeric vector containing weights for all the tips in the phylogeny. Each element must be named in order to match weights to the tips.}
}
\details{
If used, weights are passed to \code{\link{sample}} to select a subset of tips. Note that, since this is sampling without replacement, the weights "are applied sequentially,  that is the probability of choosing the next item is proportional to the probabilities amongst the remaining items". The named 'tip.weights' must contain matches to each tip, otherwise the function will exit and return a character vector of the missing tips. A zero weight will prevent a tip from being selected and the number of non-zero weights must be equal to or greater than the number of tips to be sampled.}
\value{
A list of:
  \item{pd.distrib }{A numeric vector of length 'reps' giving the PD distribution for a sample of size 'ntips' from the phylogeny.}
  \item{method }{The method passed to \code{\link{pd.calc}}.}
  ...
}
\author{David Orme}

\seealso{\code{\link{pd.calc}}}
\examples{
tre <- read.tree(text="((((A:1,B:1):1.5,C:2.5):0.5,(D:0.6,E:0.6):2.4):0.5,((F:1.9,G:1.9):0.8,(H:1.6,I:1.6):1.1):0.8):0.2;")
clmat <- clade.matrix(tre)
pd.boot <- pd.bootstrap(clmat, 5)
hist(pd.boot$pd.distrib)

}
\keyword{ utilities }% at least one, from doc/KEYWORDS
