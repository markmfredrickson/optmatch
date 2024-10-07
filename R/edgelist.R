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
##'
##' If y is a named character vector, then the names should
##' correspond to whatever in x would otherwise (i.e. if y were NULL)
##' translate to the levels set of the nodes-representing columns, while
##' the values themselves give the new levels.
##' @title Create EdgeList object
##' @param x object to convert to edgelist
##' @param y named character vector giving levels for nodes-representing columns, or NULL
##' @return EdgeList
##' @author Ben Hansen
##' @keywords internal
##' @importFrom stats complete.cases
setGeneric("edgelist", function(x, y=NULL) { stop("Not implemented.") })
setMethod("edgelist", c(x = "InfinitySparseMatrix"), function(x, y=NULL) {
    if (is.null(y))
        y  <- unlist(dimnames(x), use.names=FALSE)
    if (is.null(names(y)))
        names(y)  <- unname(y)
    elist  <-
        tibble::tibble(i = factor(x@rownames[x@rows], levels=names(y), labels=unname(y)),
                       j = factor(x@colnames[x@cols], levels=names(y), labels=unname(y)),
                       dist = x@.Data
                       )
    ccs  <- complete.cases(elist) # to remove rows involving i/j not in y
    new("EdgeList", elist[ccs,])
})

setMethod("edgelist", c(x = "matrix"), function(x, y=NULL) {
    return(edgelist(as.InfinitySparseMatrix(x), y))
})

setMethod("edgelist", c(x = "EdgeList"), function(x, y=NULL) {
    if (is.null(y) || identical(levels(x[['i']]), y) )
        return(x)

    if (is.null(names(y)))
        names(y)  <- unname(y)
    elist  <-
        tibble::tibble(i = factor(x[['i']], levels=names(y), labels=unname(y)),
                       j = factor(x[['j']], levels=names(y), labels=unname(y)),
                       dist = x[['dist']]
                       )
    ccs  <- stats::complete.cases(elist) # to remove rows involving i/j not in y
    new("EdgeList", elist[ccs,])
    })
setMethod("edgelist", c(x = "tbl_df"), function(x, y=NULL) {
    stopifnot(ncol(x)==3,
              setequal(colnames(x), c('i', 'j', 'dist')),
              is.factor(x$i) || !is.null(y),
              is.numeric(x$dist)
              )

    ccs <- stats::complete.cases(x)# to remove rows involving i/j not in y

    if (identical(colnames(x), c('i', 'j', 'dist')) &&
        is.factor(x$i) && is.factor(x$j) && is.numeric(x$dist) &&
        (is.null(y) || identical(levels(x[['i']]), y)) &&
        identical(levels(x$i), levels(x$j)) &
        all(ccs)
        )
        return(new("EdgeList", x))
    if (is.null(y))
        y  <- levels(x[['i']])
    if (is.null(names(y)))
        names(y)  <- unname(y)
    elist <-
        tibble::tibble(i = factor(x[['i']], levels=names(y), labels=unname(y)),
                       j= factor(x[['j']], levels=names(y), labels=unname(y)),
                       dist=x[['dist']]
        )
    new("EdgeList", elist[ccs,])
    })
setMethod("edgelist", c(x = "data.frame"), function(x, y=NULL){
    stopifnot(ncol(x)==3,
              setequal(colnames(x), c('i', 'j', 'dist')),
              is.factor(x$i) || !is.null(y),
              is.numeric(x$dist)
              )
    if (is.null(y))
        y  <- levels(x[['i']])
    if (is.null(names(y)))
        names(y)  <- unname(y)
    elist <-
        tibble::tibble(i = factor(x[['i']], levels=names(y), labels=unname(y)),
                       j= factor(x[['j']], levels=names(y), labels=unname(y)),
                       dist=x[['dist']]
                       )
    ccs  <- stats::complete.cases(elist) # to remove rows involving i/j not in y
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

#' @keywords internal
filter.EdgeList <- function(.data, ...) {
    .data  <- asS3(.data)
    filter(.data, ...)
    }

#' @export
as.data.frame.EdgeList <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x@.Data, col.names = x@names, ...)
}
