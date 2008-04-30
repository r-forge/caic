\name{fuscoData}
\alias{fuscoBirdData}
\alias{fuscoBirdTree}

\title{Example dataset for Fusco imbalance calculations}
\description{
This dataset contains the phylogeny of bird families and species richness originally included with the FUSCO imbalance calculation programs.
}
\usage{
data(fuscoData)
}
\format{
The dataset provides a phylogeny (\code{BigTree}) and a data frame (\code{BigTreeDat}). The phylogeny is a 200 tip tree containing some polytomies and the dataframe provides four continuous variables (yv, xc1, xc2, xc3) and one binary categorical variable (xf) for each of the tips in the tree. The dataset contains some NA values.

One reason this dataset  has been preserved as a dataset is that it also provides a reference dataset back to analyses performed in the original CAIC and MacroCAIC programs for Mac Classic. 
}


\seealso{ caic, macrocaic }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}


