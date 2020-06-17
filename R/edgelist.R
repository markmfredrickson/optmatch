################################################################################
## Turning matrices and ISMs into edge lists
################################################################################
setOldClass("data.frame")
setClass("EdgeList", contains="data.frame",
         prototype=prototype(data.frame(i=factor(), j=factor(),
                                        dist=double(0)
                                        )
                             )
         )
setValidity("EdgeList", function(object) {
    errors <- character(0)
    if (!all(colnames(object)[1:3]==c("i", "j", "dist")))
        errors  <- c(errors,
                     'Cols 1-3 should be c("i", "j", "dist")')
    if (!is.factor(object[['i']]) || !is.factor(object[['j']]) ||
        identical(levels(object[['i']]), levels(object[['']]))
        )
        errors  <- c(errors,
                     'Cols "i", "j" must be factors w/ equal level sets')
    if (!is.numeric(object[['dist']]))
        errors  <- c(errors,
                     'Col "dist" must be numeric')
    if (length(errors)==0) TRUE else errors
})
#' @export
t.EdgeList  <- function(x) new("EdgeList", data.frame(dist=x[['dist']], i=x[['j']], j=x[['i']]))

setGeneric("edgelist", function(x, y=NULL) { stop("Not implemented.") })

setMethod("edgelist", c(x = "InfinitySparseMatrix"), function(x, y=unlist(dimnames(x))) {
    elist  <- data.frame(i = factor(x@rownames[x@rows], levels=y),
                      j = factor(x@colnames[x@cols], levels=y),
                      dist = x@.Data)
    ccs  <- complete.cases(elist) # to remove rows involving i/j not in y
    new("EdgeList", elist[ccs,])
})

setMethod("edgelist", c(x = "matrix"), function(x, y=unlist(dimnames(x))) {
    return(edgelist(as.InfinitySparseMatrix(x), y))
})

setMethod("edgelist", c(x = "EdgeList"), function(x, y=levels(x[['i']])) {
    if (isTRUE(all.equal(levels(x[['i']]),
                         y, check.attributes=FALSE)
               )
        )
        return(x)

    elist  <- data.frame(i = factor(x[['i']], levels=y),
                      j = factor(x[['j']], levels=y),
                      dist = x[['dist']])
    ccs  <- complete.cases(elist) # to remove rows involving i/j not in y
    new("EdgeList", elist[ccs,])
    })
setMethod("edgelist", c(x = "data.frame"), function(x, y=levels(x[['i']])){
    stopifnot(ncol(x)==3, setequal(colnames(x), c('i', 'j', 'dist')),
              is.numeric(x$dist))
    elist <- data.frame(i = factor(x[['i']], levels=y),
                     j= factor(x[['j']], levels=y),
                     dist=x[['dist']])
    ccs  <- complete.cases(elist) # to remove rows involving i/j not in y
    new("EdgeList", elist[ccs,])
    })
