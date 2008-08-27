`print.caic` <-
function(x, ...){

    cat("Phylogenetic Independent Contrasts analysis using ",  attr(x, "contr.method"), ".\n\n", sep="")

    cat("Phylogeny: ", attr(x, "phyName"), " (",  attr(x, "origTips")  ," tips)\n", sep="")
    cat("Data: ",  attr(x, "dataName"), " (",  attr(x, "origData")  ," rows)\n", sep="")
    cat("Number of matching tips and rows: ",  attr(x, "unionData"),"\n", sep="")
    cat("Number of valid contrasts: ", sum(x$contrast.data$validNodes), "\n", sep="")

    print(summary(x))

}

