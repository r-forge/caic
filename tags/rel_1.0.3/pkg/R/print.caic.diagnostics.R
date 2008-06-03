print.caic.diagnostics <- function(x, ...){
    
    for(vars in seq(along=x)){
        cat(names(x)[vars], ":\n")
        
        datf <- array(NA, dim=c(length(x[[vars]]), 4), dimnames=list(names(x[[vars]]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
        
        for(tests in seq(along=x[[vars]])){
            datf[tests, ] <- coef(summary(x[[vars]][[tests]]))[2,]
        }
        
        print(datf)
        
    }
}