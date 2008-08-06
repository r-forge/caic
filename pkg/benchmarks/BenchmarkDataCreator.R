# code to create a benchmark dataset for use in testing code against
# existing implementations

library(CAIC)

# pure birth model with three continuous traits, one of which depends on the other two
# a binary trait and an ordered three level (can't go from A to C directly)

contTraits <- c(5,9,12)
names(contTraits) <- c("contResp","contExp1","contExp2")
contVar <- matrix(c(1,0.5,0.25,0.5,1,0,0.25,0,1), ncol=3)
discTraits <- list(binFact=matrix(c(0,0.05,0.05,0), ncol=2, dimnames=list(c("A","B"), c("A","B"))),
                   triFact=matrix(c(0,0.05,0,0.05,0,0.05,0,0.05,0), ncol=3, dimnames=list(c("A","B","C"), c("A","B","C"))))

BENCH <- growTree(b=1, halt=200,
                  ct.start=contTraits, ct.change=c(1,2,2), ct.var=contVar,
                  dt.rates=discTraits, extend.proportion=1, grain=Inf)


# make seven internal polytomies...
collapseNodes <- with( BENCH, sort(edge.length[edge[,2] > Nnode +1]))[7]
BENCHPoly <- align.tips(di2multi(BENCH, collapseNodes))

write.tree(BENCH, file="BenchTreeDi.tre")
write.tree(BENCHPoly, file="BenchTreePoly.tre")

write.caic(BENCH, filebase="BenchTreeDi")
write.caic(BENCHPoly, filebase="BenchTreePoly")

# basic data types 
# - continuous response
# - two continuous explanatory
# - integer response for species richness contrasts using a broken stick
# - binary factor
# - three-level factor

BenchData <- with(BENCH, cbind(ct.data, dt.data[,-1])) # drop the node labels from dt.data
BenchData <- BenchData[BenchData$node <= 200,] # drop the internal node data

# - variable that has no variance at two polytomies of the polytomous tree
#   simply by making all descendents of a polytomy share the same value as the start fr that trait
polyTab <- data.frame(brtime=branching.times(BENCHPoly), nDesc=table(BENCHPoly$edge[,1]))
polyTab <- subset(polyTab, nDesc.Freq > 2)
polyNoVar <- as.numeric(as.character(polyTab[order(polyTab$brtime), 2][1:2]))

clMat <- clade.matrix(BENCHPoly)
polyDesc <- clMat$clade.matrix[polyNoVar,]
polyDesc <- apply(polyDesc, 2, sum)
BenchData$polyNoVar <- ifelse(polyDesc, 9, BenchData$contExp1) 

# CAIC needs alphabetic order
CAICData <- BenchData[order(as.character(BenchData$node)),]

# also needs factors as numerics ..
CAICData$binFact <- unclass(CAICData$binFact)
CAICData$triFact <- unclass(CAICData$triFact)

# ...and as the last columns in the file.
CAICData <- CAICData[,c(1:4,7,5,6)]

write.table(CAICData, file="BenchCAIC.dat", quote=FALSE, sep="\t", row.names=FALSE)
write.table(BenchData, file="BenchData.txt", quote=FALSE, sep="\t", row.names=FALSE)

save(BenchData, BENCH, BENCHPoly, file="CAIC_Benchmark.Rda")