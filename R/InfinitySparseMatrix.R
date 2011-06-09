################################################################################
### InfinitySparseMatrix: A sparsem matrix class where infinity is the default
###                       instead of zero.
################################################################################

# instantiate with:
# a <- new("InfinitySparseMatrix", c(1,2,3), cols = c(1,2, 2), rows = c(1,2,4))
# To get a matrix like:
# 1   2
# Inf 3
setClass("InfinitySparseMatrix", 
  representation(cols = "numeric", rows = "numeric"),
  contains = "numeric")

### Basic Matrix-like Operations and Conversions ###
dim.InfinitySparseMatrix <- function(x) { c(length(x@rows) - 1, max(x@cols)) }

setMethod("as.matrix", "InfinitySparseMatrix", function(x) {
  dims <- dim(x) ; nrow <- dims[1] ; ncol <- dims[2]
  v <- rep(Inf, nrow * ncol)
  
  # we already have the column ids of the data, we need the row ids
  rowids <- rep(1:nrow, diff(x@rows))

  # combine them to get the ID of the matrix (in row major order)
  ids <- (rowids - 1) * ncol + x@cols
  v[ids] <- x

  return(matrix(v, nrow = nrow, ncol = ncol, byrow = T))
})

# setMethod("as", "matrix", "InfinitySparseMatrix", function(x) {
#  return(1)  
# })

as.InfinitySparseMatrix <- function(x) { 
  dims <- dim(x) ; nrow <- dims[1] ; ncol <- dims[2]  
  vxt <- as.vector(t(x)) # make it into row major order
  finite <- is.finite(vxt)
  
  colids <- rep(1:ncol, nrow)[finite]
  rowids <- rep(1:nrow, each = ncol)[finite]
  
  inf.by.row <- apply(is.finite(x), 1, sum)
  rows <- cumsum(c(1, inf.by.row))

  return(new("InfinitySparseMatrix", vxt[finite], cols = colids, rows = rows))
}
