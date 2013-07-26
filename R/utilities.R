################################################################################
# Utility functions used throughout the package
################################################################################

### toZ: turn various objects into a logical vector indicating treatment
### accepts a variety of inputs and keeps names/rownames if present

setGeneric("toZ", function(x) {

  if (any(is.na(x))) {
    stop("NAs not allowed in treatment indicator.")
  }

  if (is.data.frame(x) | is.matrix(x)) {
    if (dim(x)[2] > 1) {
      stop("Treatment indicators must be vectors or single column
      matrices/data.frames")
    }

    nms <- rownames(x)
    x <- x[,1]
    names(x) <- nms
  }

  if (length(unique(x)) != 2) {
    stop(paste("Treatment indicator must have exactly 2 levels not",
      length(unique(x))))
  }

  nms <- names(x)
  tmp <- standardGeneric("toZ")
  names(tmp) <- nms
  return(tmp)
})

# a noop
setMethod("toZ", "logical", identity)

# we already know it has two levels, so just call as.logical
setMethod("toZ", "numeric", function(x) as.logical(x))

setMethod("toZ", "character", function(x) toZ(as.factor(x)))

setMethod("toZ", "factor", function(x) {
  toZ(as.numeric(x) - 1)
})

#' (Internal) Remove the call before digesting a distance so things like
#' omitting caliper and calling caliper=NULL give the same digest
#'
#' @param dist Distance object to hash. Must be one of
#' \code{InfinitySparseMatrix}, \code{BlockedInfinitySparseMatrix},
#' \code{DenseMatrix}, \code{matrix}, or \code{distmatch.dlist}.
#' @return Hash on the distance object with a null \code{call}
#' @author Josh Errickson
#' @import digest
dist_digest <- function(dist) {
  if (class(dist)[1] %in% c("InfinitySparseMatrix", "BlockedInfinitySparseMatrix", "optmatch.dlist", "DenseMatrix", "matrix")) {
    csave <- attr(dist, "call")
    attr(dist, "call") <- NULL
    out <- digest(dist)
    attr(dist, "call") <- csave
    return(out)
  }
  stop("Must pass distance object")
}
