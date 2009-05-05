\name{SmallTree}
\alias{SmallTree}
\alias{diTree}
\alias{polyTree}
\alias{SmallTreeDat}
\alias{SmallTreeDatNA}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Example dataset for CAIC analyses}
\description{
SmallTree is an artificially created dataset used to illustrate methods in the CAIC package.
}
\usage{
data(SmallTree)
}
\format{
SmallTree contains two phylogenies stored in the ape package 'phylo' format. The first (\code{diTree}) is a 15 tip, fully bifurcating phylogeny. The second (\code{polyTree}) is based on \code{diTree} but has had three nodes collapsed to give three way polytomies. The are also two accompanying datasets (\code{SmallTreeDat} and \code{SmallTreeDatNA}), which contain data for each of the tips in the tree. They contain a column of tip names (tip), two continuous variables (weight and othervar), a column of the number of species in the tip group (nSpp) and two categorical variables (catX2 and catX3). \code{SmallTreeDatNA} differs only in having two lines for which no data is available.

One reason these comparatively trivial toy examples have been preserved as a dataset is that they are also used to provide a reference dataset back to analyses performed in the original CAIC and MacroCAIC programs for Mac Classic. 
}


\seealso{ caic, macrocaic }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}

