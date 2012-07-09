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
    dimension <- c(max(c(rows, length(rownames))), max(c(cols, length(colnames))))  
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
  v <- matrix(Inf, nrow = nrow, ncol = ncol, dimnames = list(treated = x@rownames, control = x@colnames))

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

  # copy the original data order, if any from the matrix
  attr(x, "order") <- attr(from, "order")
  return(x)
})

as.InfinitySparseMatrix <- function(x) { as(x, "InfinitySparseMatrix") }

# dimnames implementation

setMethod("dimnames", "InfinitySparseMatrix", function(x) {
  if (is.null(x@rownames) & is.null(x@colnames)) {
    return(NULL) 
  }
  list(treated = x@rownames, control = x@colnames)
})

setMethod("dimnames<-", "InfinitySparseMatrix", function(x, value) {
  if (length(value) != 2) {
    # message copied from matrix method
    stop(paste("length of 'dimnames' [", length(value), "] not equal to dims [2]", sep = ""))
  }
  x@rownames <- value[[1]]
  x@colnames <- value[[2]]
  x
})
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

  tmp <- e1
  tmp@.Data <- res
  tmp@cols <- e1@cols[idx1]
  tmp@rows <- e1@rows[idx1]
  
  # if either has an order, use the first
  if (!is.null(attr(e1, "order"))) {
    attr(tmp, "order") <- attr(e1, "order") 
  } else {
    attr(tmp, "order") <- attr(e2, "order") 
  }

  return(tmp)
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

# we match the subset.matrix semmantics: subset = rows, select = columns, both logical
subset.InfinitySparseMatrix <- function(x, subset, select, ...) {
  
  xdim <- dim(x)
  
  if (missing(subset)) {
    subset <- rep(TRUE, xdim[1])  
  }

  if (missing(select)) {
    select <- rep(TRUE, xdim[2])  
  }

  if (!is.logical(subset) | !is.logical(select)) {
    stop("Subset and select arguments must be logical")  
  }

  if (length(subset) != xdim[1] | length(select) != xdim[2]) {
    stop("Subset and select must be same length as rows and columns, respectively.")  
  }

  # get the indexes of the selected rows and columns
  selectedRows <- which(subset)
  selectedCols <- which(select)

  # combine the two indexes
  idx <- (x@rows %in% selectedRows) & (x@cols %in% selectedCols)
  newRowIdx <- cumsum(subset)
  newColIdx <- cumsum(select)
  return(makeInfinitySparseMatrix(x[idx],
                                  newColIdx[x@cols[idx]],
                                  newRowIdx[x@rows[idx]],
                                  colnames = x@colnames[select],
                                  rownames = x@rownames[subset]))
}

# a slightly different method of subsetting. Return a matrix of the same size, but with
# entries not in the index listed as zero
discardOthers <- function(x, index) {
  y <- x # just in case x is passed by ref, I don't think it is, but safety first
  ss <- index & !is.na(index) # from default subset implementation
  y@.Data <- y[ss]
  y@cols <- y@cols[ss]
  y@rows <- y@rows[ss]
  
  return(y)
}


cbind.InfinitySparseMatrix <- function(x, y, ...) {
  y <- as.InfinitySparseMatrix(y) # this is a noop if y is an ISM

  if(!is.null(x@rownames) & !(is.null(y@rownames))) {
    if (!all(x@rownames %in% y@rownames) | !all(y@rownames %in% x@rownames)) {
      stop("Matrices must have matching rownames.")  
    }
  }

  if (!is.null(x@colnames)) {
    xcols <- length(x@colnames)

    if (any(y@colnames %in% x@colnames)) {
      warning("Matrices share column names. This is probably a bad idea.")  
    }

  } else {
    xcols <- max(x@cols)
  }

  ymatch <- match(y@rownames, x@rownames)
  yorder <- ymatch[y@rows]

  z <- makeInfinitySparseMatrix(
    c(x@.Data, y@.Data),
    rows = c(x@rows, yorder),
    cols = c(x@cols, y@cols + xcols),
    rownames = c(x@rownames),
    colnames = make.unique(c(x@colnames, y@colnames))
  )
  return(z)

}

# I don't like the duplication between cbind and rbind, but
# I don't know that calling t(cbind(t(x), t(y))) is the correct solution
# (at least the error messages will be wrong)

rbind.InfinitySparseMatrix <- function(x, y, ...) {
  y <- as.InfinitySparseMatrix(y) # this is a noop if y is an ISM

  if(!is.null(x@colnames) & !(is.null(y@colnames))) {
    if (!all(x@colnames %in% y@colnames) | !all(y@colnames %in% x@colnames)) {
      stop("Matrices must have matching rownames")  
    }
  }

  if (!is.null(x@rownames)) {
    xrows <- length(x@rownames)

    if (any(y@rownames %in% x@rownames)) {
      warning("Matrices share row names. This is probably a bad idea.")  
    }

  } else {
    xrows <- max(x@rows)
  }

  ymatch <- match(y@colnames, x@colnames)
  yorder <- ymatch[y@cols]

  z <- makeInfinitySparseMatrix(
    c(x@.Data, y@.Data),
    rows = c(x@rows, y@rows + xrows),
    cols = c(x@cols, yorder),
    rownames = make.unique(c(x@rownames, y@rownames)),
    colnames = c(x@colnames)
  )
  return(z)

}

t.InfinitySparseMatrix <- function(x) {
  makeInfinitySparseMatrix(x@.Data, cols = x@rows, rows = x@cols, 
                           colnames = x@rownames, rownames = x@colnames)
}

setMethod("show", "InfinitySparseMatrix", function(object) { show(as.matrix(object)) })

################################################################################
# Blocked Infinity Sparse Matrix
# Just like an ISM, but keeps track of which group a unit is in
################################################################################

setClass("BlockedInfinitySparseMatrix", 
  representation(groups = "factor"),
  contains = "InfinitySparseMatrix")

# in both of the next two methods I use callGeneric(as(...), ...)
# I would have prefered callNextMethod, but I kept getting errors,
# so I manually made the call to the parent class.
setMethod("Arith", signature(e1 = "BlockedInfinitySparseMatrix", 
                             e2 = "BlockedInfinitySparseMatrix"), 
function(e1, e2) {
  tmp <- callGeneric(as(e1, "InfinitySparseMatrix"), as(e2, "InfinitySparseMatrix"))
  tmp <- as(tmp, "BlockedInfinitySparseMatrix")
  if (length(e1@groups) > length(e2@groups)) {
    tmp@groups <- e1@groups  
  } else {
    tmp@groups <- e2@groups  
  }

  return(tmp)
})

# the case where BISM is first is covered above
setMethod("Arith", signature(e1 = "InfinitySparseMatrix", 
                             e2 = "BlockedInfinitySparseMatrix"), 
function(e1, e2) {
  tmp <- callGeneric(e1, as(e2, "InfinitySparseMatrix"))
  tmp <- as(tmp, "BlockedInfinitySparseMatrix")
  tmp@groups <- e2@groups
  return(tmp)
})

# BISMs need to maintain grouping info when getting flipped
t.BlockedInfinitySparseMatrix <- function(x) {
  tmp <- t.InfinitySparseMatrix(x)
  tmp <- as(tmp, "BlockedInfinitySparseMatrix")
  tmp@groups <- x@groups
  return(tmp)
}

### Cbind/rbind
cbind.BlockedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the blocked representation to a regular ISM and call the usual cbind method
  cbind(as.InfinitySparseMatrix(x), y, ...)
}

rbind.BlockedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the blocked representation to a regular ISM and call the usual cbind method
  rbind(as.InfinitySparseMatrix(x), y, ...)
}

# Splits out the blocked matrix into its consitutent parts
setMethod("show", "BlockedInfinitySparseMatrix", function(object) { show(findSubproblems(object)) })


