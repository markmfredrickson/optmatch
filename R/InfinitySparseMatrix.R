################################################################################
### InfinitySparseMatrix: A sparsem matrix class where infinity is the default
###                       instead of zero.
################################################################################

# instantiate with:
# a <- new("InfinitySparseMatrix", c(1,2,3), cols = c(1,2, 2), rows = c(1,3,4))
# To get a matrix like:
# 1   2
# Inf 3
setClass("InfinitySparseMatrix", 
  representation(cols = "numeric", rows = "numeric", dimension = "numeric"),
  contains = "numeric")

### Basic Matrix-like Operations and Conversions ###
dim.InfinitySparseMatrix <- function(x) { 
  if (length(x@dimension) == 0) { return(c(length(x@rows) - 1, max(x@cols))) }
  else { return(x@dimension) }
}

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

setAs("InfinitySparseMatrix", "matrix", function(from) { as.matrix(from) })

# setMethod("as", "matrix", "InfinitySparseMatrix", function(x) {
#  return(1)  
# })

setAs("matrix", "InfinitySparseMatrix", function(from) { 
  dims <- dim(from) ; nrow <- dims[1] ; ncol <- dims[2]  
  vfromt <- as.vector(t(from)) # make it into row major order
  finite <- is.finite(vfromt)
  
  colids <- rep(1:ncol, nrow)[finite]
  rowids <- rep(1:nrow, each = ncol)[finite]
  
  inf.by.row <- apply(is.finite(from), 1, sum)
  rows <- cumsum(c(1, inf.by.row))

  return(new("InfinitySparseMatrix", vfromt[finite], cols = colids, rows = rows, dimension = dims))
})

as.InfinitySparseMatrix <- function(x) { as(x, "InfinitySparseMatrix") }

