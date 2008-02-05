\name{write.caic}
\alias{write.caic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Writes a phylogeny to a CAIC format file}
\description{
The programs CAIC (Purvis and Rambaut, 1995) and MacroCAIC (Agapow and Isaac, 2002) require input phylogenies to be stored in CAIC format, consisting of two files. This function exports phylogenies in 'phylo' format into CAIC format files.
}
\usage{
write.caic(phy, filebase, equal.brlen, charset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{A phylo object}
  \item{filebase}{A text string used as the file base name for the .phyl and .blen files}
  \item{equal.brlen}{If provided, all branches in the branch length file are set to this value. By default, the branch lengths in the phylo object are used or, if the phylogeny has no branch lengths, equal branch lengths of 2 are used.}
  \item{charset}{Can be used to extend the number of coding characters available for extremely polytomous trees. Defaults to LETTERS.}
}
\value{
Returns NULL.
}
\details{
Both files are simple text files with the following formats. The first file, with the suffix '.phyl', contains an alternating list of the CAIC codes (see \code{\link{caic.label}}) of the tips of the phylogeny followed by a label for the tip.  The second file, with the suffix '.blen', contains a tab delimited list of all the branches in the phylogeny with each row giving the CAIC code of the branch, followed by length of the branch and finally the height of the node. The height is not meaningful for non-ultrametric trees but is not currently calculated for ultrametric trees either. In both files, the lines are presented in alphabetical order of the CAIC codes.
}
\references{
Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent contrasts (CAIC): an Apple Macintosh application for analysing comparative data.  Computer Appl. Biosciences 11, 247-251.
Agapow, P.-M. and Isaac, N. J. B. (2002) MacroCAIC: correlates of species richness. Diversity & Distributions, 8, 41-43
}
\author{David Orme}
\seealso{\code{\link{caic.label}}}
\examples{
data(SmallTree)
write.caic(diTree, filebase="diTree")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
