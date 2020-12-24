################################################################################
## Turning matrices and ISMs into edge lists
################################################################################
setOldClass("tbl_df")
setClass("EdgeList", contains="tbl_df",
         prototype=prototype(tibble::tibble(i=factor(), j=factor(),
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
#' @importFrom tibble tibble
#' @export
t.EdgeList  <- function(x) new("EdgeList", tibble(dist=x[['dist']], i=x[['j']], j=x[['i']]))
#' @export
dim.EdgeList  <- base::dim.data.frame
setGeneric("edgelist", function(x, y=NULL) { stop("Not implemented.") })

setMethod("edgelist", c(x = "InfinitySparseMatrix"), function(x, y=unlist(dimnames(x))) {
    elist  <- tibble::tibble(i = factor(x@rownames[x@rows], levels=y),
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

    elist  <- tibble::tibble(i = factor(x[['i']], levels=y),
                      j = factor(x[['j']], levels=y),
                      dist = x[['dist']])
    ccs  <- complete.cases(elist) # to remove rows involving i/j not in y
    new("EdgeList", elist[ccs,])
    })
setMethod("edgelist", c(x = "data.frame"), function(x, y=levels(x[['i']])){
    stopifnot(ncol(x)==3, setequal(colnames(x), c('i', 'j', 'dist')),
              is.numeric(x$dist))
    elist <- tibble::tibble(i = factor(x[['i']], levels=y),
                     j= factor(x[['j']], levels=y),
                     dist=x[['dist']])
    ccs  <- complete.cases(elist) # to remove rows involving i/j not in y
    new("EdgeList", elist[ccs,])
    })


setGeneric("is_matchable",
           function(node_names, distspec, which=c("rows", "cols", "either")[3]) {
               stop("Not implemented")
           }
           )
setMethod("is_matchable", c(distspec="EdgeList"),
          function(node_names, distspec, which){
              stopifnot(is.character(node_names), is.character(which),
                        length(which)==1, nchar(which)>0,
                        sum(which==c("rows", "cols", "either"))==1
                        )
              ans  <- logical(length(node_names))
              if (which=="rows" | which=="either")
                  ans  <- ans | (node_names %in% distspec$'i')
              if (which=="cols" | which=="either")
                  ans  <- ans | (node_names %in% distspec$'j')
              ans
          })

filter.EdgeList <- function(.data, ...) {
    .data  <- asS3(.data)
    filter(.data, ...)
    }
