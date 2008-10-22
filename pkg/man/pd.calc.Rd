\name{pd.calc}
\alias{pd.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates phylogenetic diversity measurements.}
\description{
Takes a phylogeny  and returns the phylogenetic diversity (PD) associated with a set of tips. Four different measures of PD are implemented as described in the details below.}
\usage{
pd.calc(cm, tip.subset = NULL, method = "TBL", root.edge=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cm}{A object of class 'clade matrix'. Alternatively an object of class 'phylo', which will be converted to a clade.matrix.}
  \item{tip.subset}{An optional vector identifying the subset of tips to use for PD calculations. If no tip.subset is provided the method is applied to the complete phylogeny [Hmm.. this might be undesirable]. Can either be a character vector, in which case the elements are matched against tip labels, or a vector of positive integers in the range 1 to the number of tips, in which case the tips with those numbers are used.}
  \item{method}{ One of 'TBL', 'MST', 'UEH', 'SBL', defaulting to 'TBL'. See details.}
  \item{root.edge}{Logical indicating whether to include the root edge length in calculations, defaulting to FALSE.}
}
\details{
There are five implemented methods. The default is Total Branch Length (TBL), defined as the sum of all the edge lengths in the subtree given by the tip subset. This measure can be partitioned in two ways: Shared Branch Length (SBL) calculates the sum of all edges in the subtree that are shared by more than one tip; Unique Evolutionary History (UEH) calculates the sum of the edge lengths that give rise to only one tip in the subtree. Note that TBL is the sum of SBL and UEH for a particular subtree and that UEH effectively sums the edges that form tips in the \strong{subtree}. The TIP option calculates the sum of the edges leading to the subset tips in the \strong{complete tree}.  Minimum Spanning Tree (MST) calculate the length of the edges for the smallest tree that links the subset tips, excluding any edges below that node. These options are illustrated in the examples below: the red edges are those summed by each method for the tips labelled in bold. 
}
\value{
  A list giving:
  \item{pd}{The phylogenetic diversity of the tips.}
  \item{method }{The method used to calculate that PD.}
}
\references{There are some references! Insert them...}
\author{David Orme, Gavin Thomas}

\seealso{ \code{\link{pd.bootstrap}}}
\examples{
tre <- read.tree(text="((((A:1,B:1):1.5,C:2.5):0.5,(D:0.6,E:0.6):2.4):0.5,((F:1.9,G:1.9):0.8,(H:1.6,I:1.6):1.1):0.8):0.2;")
clmat <- clade.matrix(tre)
tips <- c("A","C","D","E","G","H")
pd.calc(clmat, tip.subset=tips)
pd.calc(clmat, tip.subset=c(1,3,4,5,7,8))
pd.calc(clmat, tip.subset=tips, root.edge=TRUE)

# Methods Illustrations. This is a little slow because 
# plot.phylo needs to be called repeatedly to 
# draw the correct horizontal lines only.


# demo function to colour correct branches
pd.plot <- function(tree, edges, tips){
    tip.font <- rep(1, length=length(tree$tip.label))
    tip.font[tips] <- 2
    plot(tree, label.offset=0.1, cex=1.2, font=tip.font, x.lim=3.8)
    for(edge in edges){
        col <- rep(NA, 16)
        col[edge] <- "red"
        par(new=TRUE)
        plot(tree, edge.col=col, label.offset=0.1,
             cex=1.2, font=tip.font, x.lim=3.8)
    }
}

par(mfrow=c(2,3), mar=rep(1,4))

# TBL
pd.plot(tre, c(1,2,3,4,6,7,8,9,10,11,13,14,15), c(1,3,4,5,7,8))
text(0,9, "TBL", cex=1.5, adj=c(0,0))

# UEH
pd.plot(tre, c(3,4,6,8,9,11,14,13,15), c(1,3,4,5,7,8))
text(0,9, "UEH", cex=1.5, adj=c(0,0))

# SBL
pd.plot(tre, c(1,2,7,10), c(1,3,4,5,7,8))
text(0,9, "SBL", cex=1.5, adj=c(0,0))

# TIP
pd.plot(tre, c(4,6,8,9,13,15), c(1,3,4,5,7,8))
text(0,9, "TIP", cex=1.5, adj=c(0,0))

# MST
pd.plot(tre, c(3,4,6), c(1,3))
text(0,9, "MST", cex=1.5, adj=c(0,0))

}
\keyword{utilities}% at least one, from doc/KEYWORDS
