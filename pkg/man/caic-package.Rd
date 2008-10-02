\name{caic-package}
\alias{caic-package}
\alias{CAIC}

%\alias{caicPack}
\docType{package}
\title{
A package to carry out comparative analysis using independent contrasts.
}
\description{The \pkg{CAIC} package implements methods CAIC and MacroCAIC for comparative analysis of independent contrasts and to derive correlates of species richness.
}
\details{
\tabular{ll}{
Package: \tab CAIC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-11-26\cr
License: \tab GPL\cr
}
A package to carry out comparative analysis using independent contrasts.
}
\author{
David Orme

Maintainer: David Orme <d.orme@imperial.ac.uk>
}
\references{
Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent contrasts (CAIC): an Apple Macintosh application for analysing comparative data.  Computer Appl. Biosciences 11, 247---251\\
Agapow, P.-M. and Isaac, N. J. B. (2002) MacroCAIC: correlates of species richness. Diversity \& Distributions, 8, 41---43\\
}

\keyword{ package }
\seealso{
\code{\link[ape:ape-package]{ape}}
}
\examples{
data(SmallTree)
crunch(othervar ~ weight, phy=diTree, data=SmallTreeDat, names.col=tip)
}
