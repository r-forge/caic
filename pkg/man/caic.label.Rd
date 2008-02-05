\name{caic.label}
\alias{caic.label}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create CAIC style node labels for a phylogeny}
\description{
This function returns a phylogeny object to  which CAIC style node labels have been added. The function can replace existing tip and internal node labels, append the CAIC codes to existing labels or insert a new item into the phylo structure that contains the codes in the order they appear in the phylo object edge matrix.
}
\usage{
caic.label(phy, charset = NULL, action = "insert", style = "CAIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{A 'phylo' object.}
  \item{charset}{A vector of characters to be used to construct CAIC codes. By default, LETTERS is used and a error will be generated if there are insufficient characters to code a particular polytomy.}
  \item{action}{One of 'replace', 'append' or 'insert'.}
  \item{style}{Either 'CAIC', which is the default, or 'RLE' for run-length encoded codes.}
}
\details{
CAIC codes are a way of describing a tree structure using simple text codes for each node, developed by Purvis and Rambaut (1995) for describing phylogenies for use in the program CAIC. The root node is labeled as an empty character string, and then daughter nodes are labelled with unique characters. All subsequent descendent nodes inherit their parents code and a new unique letter. With deeply nested trees, these codes can be very long and so the option 'RLE' is provided to reduce the length of codes for plotting onto trees: for example, the code AABAAAC would be converted to 2AB3AC. If \code{action} is 'replace', then the phylogeny tip and node labels are replaced with the CAIC codes. If \code{action} is 'append', then the codes are appended onto the end of existing labels, with the internal node numbers used if no node labels exist. Finally, if \code{action} is 'insert', then a character vector called 'edge.caic.code', ordered to match the rows in the phylogeny edge matrix, is inserted into the 'phylo' object list.
}
\value{
An 'phylo' object with either modified 'tip.label' and 'node.label' components or an extra character vector component called 'edge.caic.code'.
}

\references{Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent contrasts (CAIC): an Apple Macintosh application for analysing comparative data.  Computer Appl. Biosciences 11, 247-251.
}
\author{David Orme}
\seealso{ write.caic }
\examples{
data(SmallTree)
par(mfrow=c(1,2))
diCAIC <- caic.label(diTree, action="append")
plot(diCAIC, show.node=TRUE)
polyCAIC <- caic.label(polyTree, action="replace", style="RLE")
plot(polyCAIC, show.node=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

