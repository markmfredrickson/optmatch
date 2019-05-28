################################################################################
## Turning matrices and ISMs into edge lists
################################################################################

setGeneric("edgelist", function(x) { stop("Not implemented.") })

setMethod("edgelist", c(x = "InfinitySparseMatrix"), function(x) {

    return(data.frame(i = x@rownames[x@rows],  j = x@colnames[x@cols], dist = x@.Data))
})

setMethod("edgelist", c(x = "matrix"), function(x) {
    return(edgelist(as.InfinitySparseMatrix(x)))
})

setMethod("edgelist", c(x = "data.frame"), function(x){
    stopifnot(ncol(x)==3, setequal(colnames(x), c('control', 'treated', 'distance')),
              is.numeric(x$distance))
    data.frame(i = as.character(x$treated), j=as.character(x$control), dist=x$distance)
    })
