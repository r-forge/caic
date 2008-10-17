\name{benchTestInputs}

\alias{benchTestInputs}
\alias{BenchData}
\alias{BENCH}
\alias{BENCHPoly}

\alias{benchBrunchOutputs}
\alias{CAIC.BrDi813}
\alias{CAIC.BrDi913}
\alias{CAIC.BrDi1057}
\alias{CAIC.BrDi1157}
\alias{CAIC.BrPl813}
\alias{CAIC.BrPl913}
\alias{CAIC.BrPl1057}
\alias{CAIC.BrPl1157}

\alias{benchCrunchOutputs}
\alias{CAIC.CrDi213}
\alias{CAIC.CrDi657}
\alias{CAIC.CrPl213}
\alias{CAIC.CrPl413}
\alias{CAIC.CrPl657}

\alias{benchFuscoOutputs}
\alias{FuscoDiSpp}  
\alias{FuscoDiTax}  
\alias{FuscoPolySpp}
\alias{FuscoPolyTax}

\alias{benchMacroCAICOutputs}
\alias{MacroCAIC.DiSpp23}
\alias{MacroCAIC.DiSpp67}
\alias{MacroCAIC.DiTax23}
\alias{MacroCAIC.DiTax67}
\alias{MacroCAIC.PolySpp23}
\alias{MacroCAIC.PolySpp67}
\alias{MacroCAIC.PolyTax23}
\alias{MacroCAIC.PolyTax67}

\alias{benchMesaOutputs}
\alias{MeSA.I}

\title{Benchmark datasets for the CAIC package}
\description{
These datasets provide the input data used in benchmarking the implementations in this package against other implementations. They also provide data files of the outputs of those other implementations. Further details are provided in the benchmark vignette and more files, including logs and original input files to other programs, are stored in the package repository but are not included in binary packages.
}
\usage{
data(benchTestInputs)
data(benchBrunchOutputs)
data(benchCrunchOutputs)
data(benchFuscoOutputs)
data(benchMacroCAICOutputs)
data(benchMesaOutputs)
data(benchTestInputs)
}
\format{
The data file \code{benchTestInputs} contains two \code{phylo} objects: a 200 tip dichotomous, ultrametric tree (\code{BENCH}) and a version with the introduction of a few distal polytomies (\code{BENCHPoly}). It also contains a data frame (\code{BenchData}) of character data for the tips of the tree.

The data file \code{benchBrunchOutputs} contains eight data frames (\code{CAIC.BrDi813}, \code{CAIC.BrDi913}, \code{CAIC.BrDi1057}, \code{CAIC.BrDi1157}, \code{CAIC.BrPl813}, \code{CAIC.BrPl913}, \code{CAIC.BrPl1057} and \code{CAIC.BrPl1157}). These contain the outputs of CAIC running with the Brunch algorithm on the \code{BENCH} tree (\code{CAIC.BrDixxxx}) and  \code{BENCHPoly} tree (\code{CAIC.BrPlxxxx}). The numbers indicate the columns in \code{BenchData} used in the analysis. Similarly, the data file \code{benchCrunchOutputs} contains 5 data frames (\code{CAIC.CrDi213}, \code{CAIC.CrDi657}, \code{CAIC.CrPl213}, \code{CAIC.CrPl413} and \code{CAIC.CrPl657}) of the outputs of CAIC running with the Crunch algorithm.


The data file \code{benchMacroCAICOutputs} contains eight data frames  (\code{MacroCAIC.DiSpp23}, \code{MacroCAIC.DiSpp67}, \code{MacroCAIC.DiTax23}, \code{MacroCAIC.DiTax67}, \code{MacroCAIC.PolySpp23}, \code{MacroCAIC.PolySpp67}, \code{MacroCAIC.PolyTax23} and \code{MacroCAIC.PolyTax67}) of the outputs of the MacroCAIC program, again running on the \code{BENCH} tree (\code{MacroCAIC.Di\dots}) and  \code{BENCHPoly} tree (\code{MacroCAIC.Poly\dots}). The frames differ in whether the phylogeny was treated as a species-level tree (e.g. \code{MacroCAIC.DiSpp23}) or as a higher taxon tree with known numbers of species at each tip (e.g. \code{MacroCAIC.DiTax23}).


The data file \code{benchFuscoOutputs} contains the output of four sets of phylogenetic imbalance calculations (\code{FuscoDiSpp}, \code{FuscoDiTax}, \code{FuscoPolySpp} and \code{FuscoPolyTax}) from the routine `IMB\_CALC' in the  FUSCO program, again using both the dichotomous and polytomous trees and, as with the MacroCAIC outputs, treating the tips as species or higher taxa of known richness. The dataset \code{benchMesaOutputs} contains a single data frame (\code{MeSA.I}) of imbalance calculation using the modified FUSCO methods implemented in the program MeSA.
}


\seealso{crunch, brunch, macrocaic, fusco.test}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}

