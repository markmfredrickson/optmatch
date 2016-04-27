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

# for the OptionalCall data type
#' @include DenseMatrix.R
NA

# NB: this use of setOldClass also appears in mdist.R; I tried using @include
# in an roxygen block, but couldn't get it to properly order the files. Having
# two calls to setOldClass does not appear to have any harm.
setOldClass(c("optmatch.dlist", "list"))

setClassUnion("OptionalCharacter", c("character", "NULL"))

#' Objects for sparse matching problems.
#'
#' \code{InfinitySparseMatrix} is a special class of distance specifications. Finite entries
#' indicate possible matches, while infinite entries indicated non-allowed
#' matches. This data type can be more space efficient for sparse matching
#' problems.
#'
#' Usually, users will create distance specification using \code{\link{match_on}}, \code{\link{caliper}}, or
#' \code{\link{exactMatch}}.
#'
#' @author Mark M. Fredrickson
#' @seealso \code{\link{match_on}}, \code{\link{caliper}}, \code{\link{exactMatch}}, \code{\link{fullmatch}},  \code{\link{pairmatch}}
setClass("InfinitySparseMatrix",
  representation(cols = "integer",
                 rows = "integer",
                 dimension = "integer",
                 colnames = "OptionalCharacter",
                 rownames = "OptionalCharacter",
                 call = "OptionalCall"),
  contains = "vector")

# using a maker function for now, probably should be an initialize function
makeInfinitySparseMatrix <- function(data, cols, rows, colnames = NULL,
                                     rownames = NULL, dimension = NULL,
                                     call = NULL) {
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


  return(new("InfinitySparseMatrix",
             data,
             cols = as.integer(cols),
             rows = as.integer(rows),
             colnames = colnames, rownames = rownames,
             dimension = as.integer(dimension),
             call = call))
}

### Basic Matrix-like Operations and Conversions ###
#' @export
dim.InfinitySparseMatrix <- function(x) {
  return(x@dimension)
}

#' @export
as.matrix.InfinitySparseMatrix <- function(x, ...) {
  dims <- dim(x) ; nrow <- dims[1] ; ncol <- dims[2]
  v <- matrix(Inf, nrow = nrow, ncol = ncol, dimnames = list(treated = x@rownames, control = x@colnames))

  # There might be a better, vectorized way to do this, but this is correct, if slow
  n <- length(x)
  for (i in 1:n) {
    v[x@rows[i], x@cols[i]] <- x[i]
  }

  return(v)
}

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

##' Convert an object to InfinitySparseMatrix
##'
##' @param x An object which can be coerced into InfinitySparseMatrix, typically a matrix.
##' @return An InfinitySparseMatrix
##' @export
as.InfinitySparseMatrix <- function(x) { as(x, "InfinitySparseMatrix") }

# dimnames implementation
#' Get and set dimnames for InfinitySparseMatrix objects
#'
#' InfinitySparseMatrix objects represent sparse matching problems
#' with treated units as rows of a matrix and controls units as
#' the columns of the matrix. The names of the units can be retrieved
#' and set using these methods.
#'
#' @param x An InfinitySparseMatrix object.
#' @param value A list with two entries: the treated names and control names, respectively.
#' @return A list with treated and control names.
#' @rdname dimnames-InfinitySparseMatrix
#' @export
setMethod("dimnames", "InfinitySparseMatrix", function(x) {
  if (is.null(x@rownames) & is.null(x@colnames)) {
    return(NULL)
  }
  list(treated = x@rownames, control = x@colnames)
})

#' @rdname dimnames-InfinitySparseMatrix
#' @export
setMethod("dimnames<-", c("InfinitySparseMatrix", "list"), function(x, value) {
  if (length(value) != 2) {
    # message copied from matrix method
    stop(paste("length of 'dimnames' [", length(value), "] not equal to dims [2]", sep = ""))
  }
  x@rownames <- value[[1]]
  x@colnames <- value[[2]]
  x
})

#' @rdname dimnames-InfinitySparseMatrix
#' @export
setMethod("dimnames<-", c("InfinitySparseMatrix", "NULL"), function(x, value) {
  x@rownames <- NULL
  x@colnames <- NULL
  x
})

#' @importFrom Rcpp sourceCpp
NULL

################################################################################
# Infinity Sparse Matrix: Binary Ops
################################################################################
ismOpHandler <- function(binOp, e1, e2) {
  if (!identical(dim(e1), dim(e2))) {
    stop(paste("non-conformable matrices"))
  }

  # re-order e2 by the row and column names in e1
  if (!is.null(e1@rownames) && !is.null(e2@rownames)) {

    if (!all(e1@rownames == e2@rownames)) {
      # re-order e2 by e1 if they are not exactly same
      e1e2order <- match(e2@rownames, e1@rownames)
      e2@rows <- e1e2order[e2@rows]
      e2@rownames <- e1@rownames
    }
  }

  if (!is.null(e1@colnames) && !is.null(e2@colnames)) {

    if (!all(e1@colnames == e2@colnames)) {
      # re-order e2 by e1 if they are not exactly same
      e1e2order <- match(e2@colnames, e1@colnames)
      e2@cols <- e1e2order[e2@cols]
      e2@colnames <- e1@colnames
    }
  }

  return(
    ismOps(binOp, e1, e2)
  )
}

#' Element-wise addition
#'
#' \code{e1 + e2} returns the element-wise sum of
#'   two InfinitySparseMatrix objects.
#'   If either element is inf then
#'   the resulting element will be inf.
#'
#' @param e1 an InfinitySparseMatrix object
#' @param e2 an InfinitySparseMatrix object
#' @return an InfinitySparseMatrix object representing
#'   the element-wise sum of the two ISM objects
#' @docType methods
#' @rdname ismBinaryOps
#' @export
setMethod("+", signature(e1 = "InfinitySparseMatrix", e2 = "InfinitySparseMatrix"),
  function(e1, e2) ismOpHandler('+', e1, e2))

#' Elementwise subtraction
#'
#' \code{e1 - e2} returns the element-wise subtraction of
#'   two InfinitySparseMatrix objects.
#'   If either element is inf then
#'   the resulting element will be inf.
#'
#' @docType methods
#' @rdname ismBinaryOps
#' @export
setMethod("-", signature(e1 = "InfinitySparseMatrix", e2 = "InfinitySparseMatrix"),
  function(e1, e2) ismOpHandler('-', e1, e2))

#' Elementwise multiplication
#'
#' \code{e1 * e2} returns the element-wise multiplication of
#'   two InfinitySparseMatrix objects.
#'   If either element is inf then
#'   the resulting element will be inf.
#'
#' @docType methods
#' @rdname ismBinaryOps
#' @export
setMethod("*", signature(e1 = "InfinitySparseMatrix", e2 = "InfinitySparseMatrix"),
  function(e1, e2) ismOpHandler('*', e1, e2))

#' Elementwise division
#'
#' \code{e1 / e2} returns the element-wise division of
#'   two InfinitySparseMatrix objects.
#'   If either element is inf then
#'   the resulting element will be inf.
#'
#' @docType methods
#' @rdname ismBinaryOps
#' @export
setMethod("/", signature(e1 = "InfinitySparseMatrix", e2 = "InfinitySparseMatrix"),
  function(e1, e2) ismOpHandler('/', e1, e2))

setMethod("Ops", signature(e1 = "InfinitySparseMatrix", e2 = "matrix"),
function(e1, e2) {
  callGeneric(e1, as.InfinitySparseMatrix(e2))
})

setMethod("Ops", signature(e1 = "matrix", e2 = "InfinitySparseMatrix"),
function(e1, e2) {
  callGeneric(as.InfinitySparseMatrix(e1), e2)
})

setMethod("Ops", signature(e1 = "optmatch.dlist", e2 = "InfinitySparseMatrix"),
function(e1, e2) {
  callGeneric(as.matrix(e1), e2)
})

setMethod("Ops", signature(e1 = "InfinitySparseMatrix", e2 = "optmatch.dlist"),
function(e1, e2) {
  callGeneric(e1, as.matrix(e2))
})

################################################################################
# Operations over vectors and ISMs
################################################################################

binaryOpHelperVector <- function(i, v) {
  vn <- length(v)
  vpos <- ((i@cols - 1) * i@dimension[1] + i@rows) %% vn
  vpos[vpos == 0] <- vn
  return(v[vpos])
}

setMethod("Ops", signature(e1 = "InfinitySparseMatrix", e2 = "vector"),
function(e1, e2) {
  newv <- binaryOpHelperVector(e1, e2)
  e1@.Data <- callGeneric(e1@.Data, newv)
  return(e1)
})

setMethod("Ops", signature(e1 = "vector", e2 = "InfinitySparseMatrix"),
function(e1, e2) {
  newv <- binaryOpHelperVector(e2, e1)
  e2@.Data <- callGeneric(newv, e2@.Data)
  return(e2)
})

################################################################################
# Manipulating matrices: subset, cbind, rbind, etc
################################################################################

##' This matches the syntax and semantics of
##' subset for matrices.
##'
##' @title Subsetting for InfinitySparseMatrices
##' @param x InfinitySparseMatrix to be subset or bound.
##' @param subset Logical expression indicating rows to keep.
##' @param select Logical expression indicating columns to keep.
##' @param ... Other arguments are ignored.
##' @return An InfinitySparseMatrix with only the selected elements.
##' @author Mark Fredrickson
##' @export
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

  subset.data <- subsetInfSparseMatrix( as.integer(subset),
                                       select,
                                       x)
  return(makeInfinitySparseMatrix(subset.data[, 3],
                                  subset.data[, 2],
                                  subset.data[, 1],
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


##' This matches the syntax and semantics of
##' cbind and rbind for matrices.
##'
##' @title Combine InfinitySparseMatrices or
##'   BlockedInfinitySparseMatrices by row or column
##' @param x An InfinitySparseMatrix or BlockedInfinitySparseMatrix,
##'   agreeing with \code{y} in the appropriate dimension.
##' @param y An InfinitySparseMatrix or BlockedInfinitySparseMatrix,
##'   agreeing with \code{x} in the appropriate dimension.
##' @param ... Other arguments ignored.
##' @return A combined InfinitySparseMatrix or BlockedInfinitySparseMatrix
##' @author Mark Fredrickson
##' @export
##' @rdname cbindrbind
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

##' @export
##' @rdname cbindrbind
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

#' @export
t.InfinitySparseMatrix <- function(x) {
  makeInfinitySparseMatrix(x@.Data, cols = x@rows, rows = x@cols,
                           colnames = x@rownames, rownames = x@colnames)
}

##' Displays an InfinitySparseMatrix
##'
##' Specifically, displays an InfinitySparseMatrix by converting it to
##' a matrix first.
##' @param object An InfinitySparseMatrix to print.
##' @return NULL
##' @export
setMethod("show", "InfinitySparseMatrix", function(object) { show(as.matrix(object)) })

##' Sort the internal structure of an InfinitySparseMatrix.
##'
##' Internally, an InfinitySparseMatrix (Blocked or non) comprises of
##' vectors of values, row positions, and column positions. The
##' ordering of these vectors is not enforced. This function sorts the
##' internal structure, leaving the external structure unchanged
##' (e.g. `as.matrix(ism)` and `as.matrix(sort(ism))` will look
##' identical despite sorting.)
##'
##' By default, the InfinitySparseMatrix is row-dominant, meaning the
##' row positions are sorted first, then column positions are sorted
##' within each row. Use argument `byCol` to change this.
##' @param x An InfinitySparseMatrix or BlockedInfinitySparseMatrix.
##' @param decreasing Logical. Should the sort be increasing or
##'   decreasing?
##' @param ... Additional arguments ignored.
##' @param byCol Logical. Defaults to FALSE, so the returned ISM is
##'   row-dominant. TRUE returns a column-dominant ISM.
##' @return An object of the same class as `x` which is sorted
##'   according to `byCol`.
##' @rdname sort.ism
##' @export
sort.InfinitySparseMatrix <- function(x,
                                      decreasing=FALSE,
                                      ...,
                                      byCol = FALSE) {
  byCol <- as.logical(byCol)
  if (is.na(byCol)) {
    stop("byCol must be TRUE or FALSE.")
  }

  sorter <- order(attr(x, ifelse(byCol, "cols", "rows")),
                  attr(x, ifelse(byCol, "rows", "cols")),
                  decreasing=decreasing)

  makeInfinitySparseMatrix(x[sorter],
                           cols = attr(x, "cols")[sorter],
                           rows = attr(x, "rows")[sorter],
                           rownames = attr(x, "rownames"),
                           colnames = attr(x, "colnames"),
                           dimension = attr(x, "dimension"),
                           call = attr(x, "call"))

}


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
setMethod("Ops", signature(e1 = "BlockedInfinitySparseMatrix",
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
setMethod("Ops", signature(e1 = "InfinitySparseMatrix",
                             e2 = "BlockedInfinitySparseMatrix"),
function(e1, e2) {
  tmp <- callGeneric(e1, as(e2, "InfinitySparseMatrix"))
  tmp <- as(tmp, "BlockedInfinitySparseMatrix")
  tmp@groups <- e2@groups
  return(tmp)
})

# BISMs need to maintain grouping info when getting flipped
#' @export
t.BlockedInfinitySparseMatrix <- function(x) {
  tmp <- t.InfinitySparseMatrix(x)
  tmp <- as(tmp, "BlockedInfinitySparseMatrix")
  tmp@groups <- x@groups
  return(tmp)
}

##' @export
##' @rdname cbindrbind
cbind.BlockedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the blocked representation to a regular ISM and call the usual cbind method
  cbind(as.InfinitySparseMatrix(x), y, ...)
}

##' @export
##' @rdname cbindrbind
rbind.BlockedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the blocked representation to a regular ISM and call the usual cbind method
  rbind(as.InfinitySparseMatrix(x), y, ...)
}

#' Returns the dimension of each valid subproblem
#'
#' Returns a list containing the dimensions of all valid subproblems.
#' @param x A distance specification to get the sub-dimensions of.
#' @return A list of the dimensions of each valid subproblem. Any subproblems with 0 controls
#' or 0 treatments will be ignored. The names of the entries in the list will be the names of the
#' subproblems, if they exist.
#' @export
#' @docType methods
#' @rdname subdim-methods
#' @export
subdim <- function(x) {
  UseMethod("subdim")
}

#' @rdname subdim-methods
#' @export
subdim.InfinitySparseMatrix <- function(x) {
  list(dim(x))
}

#' @rdname subdim-methods
#' @export
subdim.matrix <- function(x) {
  list(dim(x))
}

#' @rdname subdim-methods
#' @export
subdim.BlockedInfinitySparseMatrix <- function(x) {
  out <- lapply(levels(x@groups), function(k) c(sum(row.names(x) %in% names(x@groups)[x@groups == k]),
                                                sum(colnames(x) %in% names(x@groups)[x@groups == k])))
  names(out) <- levels(x@groups)
  # drop off any subproblems lacking at least one treatment/control
  out[unlist(lapply(out, function(t) all(t > 0)))]
}

##' @rdname sort.ism
##' @export
sort.BlockedInfinitySparseMatrix <- function(x,
                                             decreasing=FALSE,
                                             ...,
                                             byCol=FALSE) {
  y <- sort.InfinitySparseMatrix(x, decreasing=decreasing, byCol=byCol)
  attr(y, "groups") <- attr(x, "groups")
  class(y) <- class(x)
  y
}


#' Returns the number of eligible matches for the distance.
#'
#' This will return a list of the number of finite entries in a distance
#' matrix. If the distance has no subgroups, it will be a list of length 1. If
#' the distance has subgroups (i.e. \code{x} is an
#' \code{BlockedInfinitySparseMatrix}, it will be a named list.)
#' @param x Any distance object.
#' @return A list counting the number of eligible matches in the distance.
#' @export
#' @docType methods
#' @rdname num_eligible_matches-methods
num_eligible_matches <- function(x) {
  UseMethod("num_eligible_matches")
}

#' @rdname num_eligible_matches-methods
#' @export
num_eligible_matches.optmatch.dlist <-function(x) {
    lapply(x, function(x) sum(is.finite(x)))
}

#' @rdname num_eligible_matches-methods
#' @export
num_eligible_matches.matrix <- function(x) {
  list(sum(is.finite(x)))
}

#' @rdname num_eligible_matches-methods
#' @export
num_eligible_matches.InfinitySparseMatrix <- function(x) {
  list(sum(is.finite(x@.Data)))
}

#' @usage \method{num_eligible_matches}{BlockedInfinitySparseMatrix}(x)
#' @rdname num_eligible_matches-methods
#' @export
num_eligible_matches.BlockedInfinitySparseMatrix <- function(x) {
  out <- lapply(levels(x@groups), function(k) sum(is.finite(x[x@groups == k]@.Data)))
  names(out) <- levels(x@groups)
  out
}

##' Displays a BlockedInfinitySparseMatrix
##'
##' Displays each block of the BlockedInfinitySparseMatrix separately.
##' @param object An BlockedInfinitySparseMatrix to print.
##' @return NULL
##' @export
setMethod("show", "BlockedInfinitySparseMatrix", function(object) { show(findSubproblems(object)) })
