\name{phylo.d}
\alias{phylo.d}
\alias{print.phylo.d}
\alias{summary.phylo.d}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates the phylogenetic D statistic}
\description{
Calculates the D value, a measure of phylogenetic signal in a binary trait, and tests the estimated D value for significant departure from one (phylogenetically random) and zero (phylogenetically as clumped as expected under a Brownian evolution threshold model).
}
\usage{
phylo.d(ds, phy, binvar, tipnames, permut = 1000, verbose = TRUE)
\method{print}{phylo.d}(x, ...)
\method{summary}{phylo.d}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ds}{An object of class 'data.frame'.}
  \item{phy}{An object of class 'phylo'.}
  \item{binvar}{The name of the column in \code{ds} holding the binary variable of interest.}
  \item{tipnames}{The name of the column in \code{ds} holding the names to be used to match rows in the dataset (\code{ds}) to the tipnames in the phylogeny (\code{phy})}
  \item{permut}{Number of permutations to be used in the randomisation test.}
  \item{verbose}{If TRUE, the function reports progress and other output likely to be important to users}
  \item{x}{An object of class 'phylo.d'}
  \item{object}{An object of class 'phylo.d'}
}

\details{
D is one if there is no phylogenetic signal in the trait, and zero if the phylogenetic signal is the same as that expected under a Brownian threshold model (see also reference for details).
Calculation of D is based on the observed sum of sister-clade differences for the binary trait of interest, which is scaled using two expectations: a phylogenetically random distribution of zeros and ones at the tips of the phylogeny, and a phylogenetically clumped distribution of ones that is simulated under a Brownian model. The scaling and significance testing uses 1000 permutations each (default) for the random and Brownian expectations, which are based on the phylogeny and binary data of interest, so that D itself becomes independent of phylogeny size and proportion of ones in the binary trait (prevalence).
For the random expectation, trait values are randomly shuffled along the tips of the phylogeny. For the Brownian expectation, a continuous trait is evolved under Brownian evolution in each permutation, and then the sum of sister-clade differences is computed for a binary trait ranked on the basis of that continuous trait (values higher than a threshold are ones, and the threshold is defined so that the simulated binary trait has the same prevalence as the observed one).
}
\value{
Returns an object of class 'phylo.d', which is a list of the following:
\item{DEstimate} {The estimated D value}
\item{Pval1} {A p value, giving the result of testing whether D is significantly different from one}
\item{Pval0} {A p value, giving the result of testing whether D is significantly different from zero}
\item{Parameters} {A list of the Observed, MeanRandom and MeanBrownian sums of sister-clade differences}
\item{Permutations} {A list with elements random and brownian, containing the sums of sister-clade differences from random permutations and simulations of Brownian evolution under a threshold model}
\item{NodalVals} {A list with elements observed, random and brownian, containing the nodal values estimated for the observed trait (as a numeric with names that correspond to node numbers in the phylo object) and those estimated in the permutations (as dataframes with a column for each node, named e.g. X17 for node 17)}
\item{binvar}{The binary variable used}
\item{phyName}{The name of the phylogeny object used}
\item{dsName}{The name of the dataframe used}
\item{nPermut}{The number of permutations used}
}
}
\references{ S. A. Fritz & A. Purvis (2009) In prep. }
\author{Susanne Fritz (SFritz@bio.ku.dk)}
\examples{
data(shorebird)

# add threatened status (threat data from IUCN Red List 2008)
# and a phylogenetically clumped dummy example to the data frame
shorebird.data$threat <- rep(c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0), times=c(22,2,2,2,1,1,4,1,8,1,3,1,1,2,20))
shorebird.data$clustered <- rep(c(1,0,1,0,1,0,1,0,1,0,1,0), times=c(4,13,2,12,1,2,1,7,2,2,2,23))

# plot the data
par(mfrow=c(1,2))
tipOrd <- match(shorebird.data$Species, shorebird.tree$tip.label)
plot(shorebird.tree, show.tip.label=FALSE) 
tiplabels(tip=tipOrd, pch=21, bg=ifelse(shorebird.data$threat==0, 'white','black'))
plot(shorebird.tree, show.tip.label=FALSE) 
tiplabels(tip=tipOrd, pch=21, bg=ifelse(shorebird.data$clustered==0, 'white','black'))

# an example where D is not significantly different from a brownian threshold model of evolution (0)
ex1 <- phylo.d(shorebird.data, shorebird.tree, threat, Species)
# an example where D is not significantly different from random (1)
ex2 <- phylo.d(shorebird.data, shorebird.tree, clustered, Species)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ htest }% __ONLY ONE__ keyword per line
