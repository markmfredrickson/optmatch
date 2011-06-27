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
setClassUnion("OptionalCharacter", c("character", "NULL"))
setClass("InfinitySparseMatrix", 
  representation(cols = "numeric", rows = "numeric", dimension = "numeric",
    colnames = "OptionalCharacter", rownames = "OptionalCharacter"),
  contains = "numeric")

# using a maker function for now, probably should be an initialize function
makeInfinitySparseMatrix <- function(data, cols, rows, colnames = NULL, rownames = NULL, dimension = NULL) {
  if (!all.equal(length(data), length(cols), length(rows))) {
    stop("Data and column/row ids must be vectors of the same length")  
  }

  if(is.null(dimension)) {
    dimension <- c(max(rows), max(cols))  
  }  

  # if(is.null(rownames)) {
  #   rownames <- paste("T", 1:dimension[1], sep = "")  
  # }

  # if(is.null(colnames)) {
  #   colnames <- paste("C", 1:dimension[2], sep = "")  
  # }


  return(new("InfinitySparseMatrix", data, cols = cols, rows = rows, colnames = colnames, rownames =
    rownames, dimension = dimension))
}

### Basic Matrix-like Operations and Conversions ###
dim.InfinitySparseMatrix <- function(x) { 
  return(x@dimension)
}

setMethod("as.matrix", "InfinitySparseMatrix", function(x) {
  dims <- dim(x) ; nrow <- dims[1] ; ncol <- dims[2]
  v <- matrix(Inf, nrow = nrow, ncol = ncol, dimnames = list(x@rownames, x@colnames))

  # There might be a better, vectorized way to do this, but this is correct, if slow
  n <- length(x)
  for (i in 1:n) {
    v[x@rows[i], x@cols[i]] <- x[i]  
  }

  return(v)
})

setAs("InfinitySparseMatrix", "matrix", function(from) { as.matrix(from) })

# setMethod("as", "matrix", "InfinitySparseMatrix", function(x) {
#  return(1)    
# })

setAs("matrix", "InfinitySparseMatrix", function(from) { 
  dims <- dim(from) ; nrow <- dims[1] ; ncol <- dims[2]  
  finite <- is.finite(from)
  
  colids <- rep(1:ncol, each = nrow)[finite]
  rowids <- rep(1:nrow, ncol)[finite]
  
  x <- makeInfinitySparseMatrix(as.vector(from[finite]),
                                cols = colids,
                                rows = rowids,
                                rownames = rownames(from),
                                colnames = colnames(from),
                                dimension = dims)
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

################################################################################
# Infinity Sparse Matrix: Math Ops
################################################################################
setMethod("Arith", signature(e1 = "InfinitySparseMatrix", e2 = "InfinitySparseMatrix"), 
function(e1, e2) {
  if (!identical(dim(e1), dim(e2))) {
    stop(paste("non-conformable matrices"))  
  }  

  pairs1 <- mapply(c, e1@rows, e1@cols, SIMPLIFY = F)
  pairs2 <- mapply(c, e2@rows, e2@cols, SIMPLIFY = F)
  
  # Note: This might be expensive. There may be a way to do this in one pass if the data
  # were pre-sorted in into row/column order
  idx1 <- which(pairs1 %in% pairs2)
  idx2 <- which(pairs2 %in% pairs1)

  data1 <- e1@.Data[idx1]
  data2 <- e2@.Data[idx2]

  res <- callGeneric(data1, data2)

  return(makeInfinitySparseMatrix(res, cols = e1@cols[idx1], rows = e1@rows[idx1], 
    dimension = e1@dimension, colnames = e1@colnames, rownames = e1@rownames))
  
})

setMethod("Arith", signature(e1 = "InfinitySparseMatrix", e2 = "matrix"), 
function(e1, e2) {
  callGeneric(e1, as.InfinitySparseMatrix(e2))
})

setMethod("Arith", signature(e1 = "matrix", e2 = "InfinitySparseMatrix"), 
function(e1, e2) {
  callGeneric(as.InfinitySparseMatrix(e1), e2)
})



################################################################################
# Manipulating matrices: subset, cbind, rbind, etc
################################################################################

subset.InfinitySparseMatrix <- function(x, subset, ...) {
  y <- x # just in case x is passed by ref, I don't think it is, but safety first
  ss <- subset & !is.na(subset) # from default subset implementation
  y@.Data <- y[ss]
  y@cols <- y@cols[ss]
  y@rows <- y@rows[ss]

  # should we purge cols and rows? Otherwise, we have the same size matrix, just 
  # with fewer finit entries. 
  # this might be better thought of as a replace method

  return(y)
}




