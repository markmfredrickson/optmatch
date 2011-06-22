################################################################################
### InfinitySparseMatrix: A sparsem matrix class where infinity is the default
###                       instead of zero.
################################################################################

# instantiate with:
# a <- new("InfinitySparseMatrix", c(1,2,3), cols = c(1,2, 2), rows = c(1,1,2))
# data is the data, cols are column IDs, rows are row IDs
# To get a matrix like:
# 1   2
# Inf 3
setClass("InfinitySparseMatrix", 
  representation(cols = "numeric", rows = "numeric", dimension = "numeric",
    colnames = "character", rownames = "character"),
  contains = "numeric")

# using a maker function for now, probably should be an initialize function
makeInfinitySparseMatrix <- function(data, colids, rowids, colnames = NULL, rownames = NULL, dimension = NULL) {
  if (!all.equal(length(data), length(colids), length(rowids))) {
    stop("Data and column/row ids must be vectors of the same length")  
  }

  if(is.null(dimension)) {
    dimension <- c(max(rowids), max(colids))  
  }  

  if(is.null(rownames)) {
    rownames <- paste("T", 1:dimension[1], sep = "")  
  }

  if(is.null(colnames)) {
    colnames <- paste("C", 1:dimension[2], sep = "")  
  }


  return(new("InfinitySparseMatrix", data, cols = colids, rows = rowids, colnames = colnames, rownames =
    rownames, dimension = dimension))
}

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

  return(matrix(v, nrow = nrow, ncol = ncol, byrow = T, dimnames = list(x@rownames, x@colnames)))
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
  
  x <- new("InfinitySparseMatrix", vfromt[finite], cols = colids, rows = rows, dimension = dims)
  if (!is.null(rownames(from))) {
    x@rownames <- rownames(from)
  }
  if (!is.null(colnames(from))) {
    x@colnames <- colnames(from)
  }
  return(x)
})

as.InfinitySparseMatrix <- function(x) { as(x, "InfinitySparseMatrix") }

# dimnames implementation

setMethod("colnames<-", "InfinitySparseMatrix", function(x, value) {
    x@colnames <- value
    return(x)
})

setMethod("colnames", "InfinitySparseMatrix", function(x) {
  return(x@colnames)  
})

setMethod("rownames<-", "InfinitySparseMatrix", function(x, value) {
    x@rownames <- value
    return(x)
})

setMethod("rownames", "InfinitySparseMatrix", function(x) {
  return(x@rownames)  
})

