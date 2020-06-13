################################################################################
## Turning matrices and ISMs into edge lists
################################################################################

setGeneric("edgelist", function(x, y) { stop("Not implemented.") })

setMethod("edgelist", c(x = "InfinitySparseMatrix"), function(x, y) {
    el  <- data.frame(i = factor(x@rownames[x@rows], levels=y),
                      j = factor(x@colnames[x@cols], levels=y),
                      dist = x@.Data)
    ccs  <- complete.cases(el) # to remove rows involving i/j not in y
    el[ccs,]
})

setMethod("edgelist", c(x = "matrix"), function(x, y) {
    return(edgelist(as.InfinitySparseMatrix(x), y))
})

setMethod("edgelist", c(x = "data.frame"), function(x, y){
    stopifnot(ncol(x)==3, setequal(colnames(x), c('control', 'treated', 'distance')),
              is.numeric(x$distance))
    el <- data.frame(i = factor(x$treated, levels=y),
                     j= factor(x$control, levels=y),
                     dist=x$distance)
    ccs  <- complete.cases(el) # to remove rows involving i/j not in y
    el[ccs,]    
    })
