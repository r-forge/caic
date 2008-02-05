`predict.caic` <-
function(object, ...){
    
    # need to force the model to get predictions using the contrast table rather than the original data table...
    # don't completely hijack the newdata argument...
    
    dots <- list(...)
    newdataProv <- pmatch(names(dots), "newdata")
    if(all(is.na(newdataProv))) nD <- caic.table(object) else nD <- dots[[newdataProv]]
    predict(object$mod, newdata=nD)
    
    
}

