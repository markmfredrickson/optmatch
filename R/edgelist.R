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


