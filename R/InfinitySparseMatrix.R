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
#' Usually, users will create distance specification using \code{\link{match_on}}, \code{\link{caliper}}, or
#' \code{\link{exactMatch}}. The ordering of units in an \code{InfinitySparseMatrix} is not guaranteed to be maintained after subsetting and/or other operations are performed.
#' @seealso \code{\link{match_on}}, \code{\link{caliper}}, \code{\link{exactMatch}}, \code{\link{fullmatch}},  \code{\link{pairmatch}}
#' @author Mark M. Fredrickson
#' @template ISMslotsTemplate
setClass("InfinitySparseMatrix",
         slots = c(cols = "integer",
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

  # #190 - dimension names must agree (agnostic to ordering)
  # If there are no dimension names, skip this check
  if ((!is.null(dimnames(e1)) & is.null(dimnames(e2))) |
        (is.null(dimnames(e1)) & !is.null(dimnames(e2)))) {
    warning("One matrix has dimnames and the other does not. Proper ordering is not guaranteed.")
  }
  else if (!is.null(dimnames(e1)) | !is.null(dimnames(e2))) {
    t1 <- dimnames(e1)$treated
    t2 <- dimnames(e2)$treated
    c1 <- dimnames(e1)$control
    c2 <- dimnames(e2)$control
    t1extra <- setdiff(t1, t2)
    t2extra <- setdiff(t2, t1)
    c1extra <- setdiff(c1, c2)
    c2extra <- setdiff(c2, c1)
    if (length(c(t1extra, t2extra, c1extra, c2extra)) > 0) {
      error <- ""
      if (length(t1extra > 0)) {
        error <- paste0(error, "    extra rows in first matrix: ",
                        paste0(t1extra, collapse = ", "), "\n")
      }
      if (length(t2extra > 0)) {
        error <- paste0(error, "    extra rows in second matrix: ",
                        paste0(t2extra, collapse = ", "), "\n")
      }
      if (length(c1extra > 0)) {
        error <- paste0(error, "    extra columns in first matrix: ",
                        paste0(c1extra, collapse = ", "), "\n")
      }
      if (length(c2extra > 0)) {
        error <- paste0(error, "    extra columns in second matrix: ",
                        paste0(c2extra, collapse = ", "), "\n")
      }
      stop(paste0("dimension names between matrices must be in agreement:\n",
                 error))
    }
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
##' @rdname ism.subset
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
                                  rownames = x@rownames[subset],
                                  dimension = c(sum(subset), sum(select))))
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

##' @param i Row indices.
##' @param j Col indices.
##' @param drop Ignored.
##' @rdname ism.subset
##' @export
setMethod("[", "InfinitySparseMatrix",
          function(x, i, j =NULL, ..., drop = TRUE) {
            if ( "drop" %in% names(match.call())) {
              warning("'drop' argument ignored for InfinitySparseMatrix object subsetting ")
            }
            # Handles [X] cases
            if (nargs() < 3) {
              if (missing(i)) {
                return(x)
              }
              return(x@.Data[i, ...])
            } else {
              # At this point we have two arguments, but one could be
              # null (e.g. [X,] or [,X], as opposed to [X])

              # when missing, replace with NULL to be handled below
              if (missing(i)) i <- NULL
              if (missing(j)) j <- NULL

              makelogical <- function(index, rowcol) {
                switch(class(index),
                       "numeric" = ,
                       "integer" = {
                         if (any(index < 0)) {
                           if (any(index >= 0)) {
                             stop("Cannot mix positive and negative subscripts")
                           }
                           !((1:dim(x)[rowcol]) %in% abs(index))
                         } else {
                           (1:dim(x)[rowcol]) %in% index
                         }
                       },
                       "character" = dimnames(x)[[rowcol]] %in% index,
                       "logical" = index,
                       "NULL" = rep(TRUE, dim(x)[rowcol]),
                       stop("Unrecognized class"))
              }

              subi <- makelogical(i, 1)
              subj <- makelogical(j, 2)

              subset(x, subset=subi, select=subj)
            }
          })


##' @param value replacement values
##' @rdname ism.subset
##' @export
setMethod("[<-", "InfinitySparseMatrix",
          function(x, i, j, value) {
#            if (is(x, "BlockedInfinitySparseMatrix")) x <- as.InfinitySparseMatrix(x)

            s <- sys.calls() # look at calling function to determine [X,X] vs [X]
            if (length(s[[length(s)-1]]) < 5) {
              # handles x[i] <- ...
              if (is.logical(i)) {
                i <- which(i)
              }
              x@.Data[i] <- value
              return(x)
            }

            if (missing(i)) i <- seq_len(nrow(x))
            if (missing(j)) j <- seq_len(ncol(x))

            makenumeric <- function(index, rowcol) {
              switch(class(index),
                     "numeric" = ,
                     "integer" = {
                       if (any(index < 0)) {
                         if (any(index >= 0)) {
                           stop("Cannot mix positive and negative subscripts")
                         }
                         which(!((1:dim(x)[rowcol]) %in% abs(index)))
                       } else {
                         index
                       }
                     },

                     "character" = which(dimnames(x)[[rowcol]] %in% index),
                     "logical" = which(index),
                     "NULL" = seq_len(dim(x)[rowcol]),
                     stop("Unrecognized class"))
            }

            numi <- makenumeric(i, 1)
            numj <- makenumeric(j, 2)


            # Create a list of all coordinates to be updated
            updateEntries <- expand.grid(numi, numj)
            if (nrow(updateEntries) %% length(value) != 0) {
              stop("number of items to replace is not a multiple of replacement length")
            }
            updateEntries <- cbind(updateEntries, as.vector(value))

            for (r in seq_len(nrow(updateEntries))) {
              if (is(x, "BlockedInfinitySparseMatrix")) {
                rowgroups <- x@groups[x@rownames]
                colgroups <- x@groups[x@colnames]
                # If we're trying to replace a cross-group term in BISM, conver to ISM
                if (rowgroups[updateEntries[r, 1]] != colgroups[updateEntries[r, 2]]) {
                  x <- as.InfinitySparseMatrix(x)
                }
              }
              # determine position in @.Data/@rows/@cols for replacement
              datalocation <- x@rows == updateEntries[r,1] &
                              x@cols == updateEntries[r,2]

              # Different steps depending on finite/infinite status of
              # replacement and current value

              # Replacement Inf, current Inf is ignored

              # Replacement Inf, current finite, delete it
              if (is.infinite(updateEntries[r, 3])) {
                x@rows <- x@rows[!datalocation]
                x@cols <- x@cols[!datalocation]
                x@.Data <- x@.Data[!datalocation]
              } else {
                if (any(datalocation)) {
                  # Replacement finite, current finite, replace it
                  x@.Data[datalocation] <- updateEntries[r,3]
                } else {
                  # Replacement finite, current Inf, add new entry
                  x@rows <- c(x@rows, as.integer(updateEntries[r,1]))
                  x@cols <- c(x@cols, as.integer(updateEntries[r,2]))
                  x@.Data <- c(x@.Data, as.integer(updateEntries[r,3]))
                }
              }
            }
            return(sort(x))
          })



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
      stop("Matrices must have matching colnames")
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
##' Internally, an \code{InfinitySparseMatrix} (Blocked or non) comprises of
##' vectors of values, row positions, and column positions. The ordering of
##' these vectors is not enforced. This function sorts the internal structure,
##' leaving the external structure unchanged (e.g. \code{as.matrix(ism)} and
##' \code{as.matrix(sort(ism))} will look identical despite sorting.)
##'
##' By default, the \code{InfinitySparseMatrix} is row-dominant, meaning the row
##' positions are sorted first, then column positions are sorted within each
##' row. Use argument \code{byCol} to change this.
##' @param x An \code{InfinitySparseMatrix} or
##'   \code{BlockedInfinitySparseMatrix}.
##' @param decreasing Logical. Should the sort be increasing or decreasing?
##'   Default \code{FALSE}.
##' @param ... Additional arguments ignored.
##' @param byCol Logical. Defaults to \code{FALSE}, so the returned ISM is
##'   row-dominant. \code{TRUE} returns a column-dominant ISM.
##' @return An object of the same class as \code{x} which is sorted according to
##'   \code{byCol}.
##' @rdname sort.ism
##' @export
sort.InfinitySparseMatrix <- function(x,
                                      decreasing=FALSE,
                                      ...,
                                      byCol = FALSE) {
  byCol <- as.logical(byCol)
  testByCol <- length(byCol) > 1 | any(is.na(byCol))
  if (testByCol) {
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

#' Blocked Infinity Sparse Matrix
#'
#' Blocked Infinity Sparse Matrices are similar to Infinity Sparse Matrices, but they also keep track of the groups of units via an additional slot, \code{groups}
#' @slot groups factor vector containing groups, with unit names as labels, when possible
#' @template ISMslotsTemplate
#'
#' @seealso \code{\link{match_on}}, \code{\link{exactMatch}}, \code{\link{fullmatch}},  \code{\link{InfinitySparseMatrix-class}}
#' @author Mark M. Fredrickson
setClass("BlockedInfinitySparseMatrix",
  slots = c(groups = "factor"),
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
  cbind.InfinitySparseMatrix(x, y, ...)
}

##' @export
##' @rdname cbindrbind
rbind.BlockedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the blocked representation to a regular ISM and call the usual cbind method
  rbind.InfinitySparseMatrix(x, y, ...)
}

#' Returns the dimension of each valid subproblem
#'
#' Returns a list containing the dimensions of all valid subproblems.
#' @param x A distance specification to get the sub-dimensions of.
#' @return A data frame listing the dimensions of each valid subproblem. Any subproblems with 0 controls
#' or 0 treatments will be ignored. The names of the entries in the list will be the names of the
#' subproblems, if they exist.  There will be two rows, named "treatments" and "controls".
#' @export
#' @docType methods
#' @rdname subdim-methods
#' @example inst/examples/subdim.R
#' @export
subdim <- function(x) {
  UseMethod("subdim")
}

#' @rdname subdim-methods
#' @export
subdim.InfinitySparseMatrix <- function(x) {
  data.frame("._"=dim(x), row.names=c("treatments", "controls"))
}

#' @rdname subdim-methods
#' @export
subdim.matrix <- function(x) {
  data.frame("._"=dim(x), row.names=c("treatments", "controls"))
}

#' @rdname subdim-methods
#' @export
subdim.BlockedInfinitySparseMatrix <- function(x) {
  out <- lapply(levels(x@groups), function(k) c(sum(row.names(x) %in% names(x@groups)[x@groups == k]),
                                                sum(colnames(x) %in% names(x@groups)[x@groups == k])))
  names(out) <- levels(x@groups)
  # drop off any subproblems lacking at least one possible treatment-control pairing
  filt <- vapply(levels(x@groups), function(l) {
      members <- names(x@groups[x@groups == l])
      row.members <- which(x@rownames %in% members)
      col.members <- which(x@colnames %in% members)
      ridx <- x@rows %in% row.members
      cidx <- x@cols %in% col.members
      any(ridx & cidx)
  },
  logical(1))
  out <- out[filt]
  out.cnms <- names(out)
  out <- as.data.frame(out)
  colnames(out) <- out.cnms
  row.names(out) <- c("treatments", "controls")
  out
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

##' This function generates a single block-diagonal distance matrix given
##' several distance matrices defined on subgroups.
##'
##' When you've generated several distances matrices on subgroups in your
##' analysis, you may wish to combine them into a single block-diagonal distance
##' matrix. The \code{dbind} function facilitates this.
##'
##' Any \code{BlockedInfinitySparseMatrix} include in \code{...} will be broken
##' into individual \code{InfinitySparseMatrix} before being joined back
##' together. For example, if \code{b} is a \code{BlockedInfinitySparseMatrix}
##' with 2 subgroups and \code{m} is a distance without subgroups, then
##' \code{dbind(b, m)} will be a \code{BlockedInfinitySparseMatrix} with 3
##' subgroups.
##'
##' If there are any shared names (either row or column) among all distances
##' passed in, by default all matrices will be renamed to ensure unique names by
##' appending "X." to each distance, where "X" is ascending lower case letters
##' ("a.", "b.", etc). Setting the \code{force_unique_names} argument to
##' \code{TRUE} errors on this instead.
##'
##' If the matrices need to be renamed and there are more than 26 separate
##' matrices, after the first 26 single "X." prefixs, they will continue as
##' "YX.", e.g "aa.", "ab.", "ac.". If more than 676 separate matrices, the
##' prefix wil continue to "ZYX.", e.g. "aaa.", "aab.", "aac.". This scheme
##' supports up to 18,278 unique matrices.
##'
##' Note that you do \strong{not} have to combine subgroup distances into a
##' single blocked distance using this function to ultimately obtain a single
##' matching set. Instead, take a look at the vignette
##' \code{vignette("matching-within-subgroups", package = "optmatch")} for
##' details on combining multiple matches.
##' @title Diagonally bind together subgroup-specific distances
##' @param ... Any number of distance objects which can be converted to
##'   \code{InfinitySparseMatrix}, such as class \code{matrix},
##'   \code{DenseMatrix}, \code{InfinitySparseMatrix}, or
##'   \code{BlockedInfinitySparseMatrix}, or \code{list}s containing distance
##'   objects.
##' @param force_unique_names Default \code{FALSE}. When row or column names are
##'   not unique among all distances, if \code{FALSE}, throw a warning and
##'   rename all rows and columns to ensure unique names. If \code{TRUE}, error
##'   on non-unique names.
##' @return A \code{BlockedInfinitySparseMatrix} containing a block-diagonal
##'   distance matrix. If only a single distance is passed to \code{dbind} and
##'   it is not already a \code{BlockedInfinitySparseMatrix}, the result will be
##'   an \code{InfinitySparseMatrix} instead.
##' @importFrom methods slot
##' @export
##' @examples
##' data(nuclearplants)
##' m1 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 0),
##'                caliper = 1)
##' m2 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 1),
##'                caliper = 1.3)
##' blocked <- dbind(m1, m2)
##'
##' dists <- list(m1, m2)
##'
##' blocked2 <- dbind(dists)
##' identical(blocked, blocked2)
dbind <- function(..., force_unique_names = FALSE) {
  mats <- list(...)

  if (length(mats) == 1 & !inherits(mats[[1]], "list")) {
    # If passed a single distance, return it as ISM or BISM
    return(.as.ism_or_bism(mats[[1]]))
  }

  # Below, if one of the entries is a BISM, we'll end up with a list of lists,
  # where the BISM is replaced with a list of ISMs. The following flattens this
  # into a single list of ISMs but seems overly complicated; basically for any
  # non-list inside `l`, it sticks them inside a sub-list, then the `unlist(...,
  # recursive = FALSE)` breaks the outer-most list away. Generally `unlist(...,
  # recursive = FALSE)` is suggested to work without additional concerns in most
  # cases, but that fails here because an ISM would be unlisted to return a
  # vector. I really don't like this solution but can't find anything better.
  #
  # Visualization of this
  # Starting:
  # List of 2
  #     $ InfinitySparseMatrix
  #     $ List of 2
  #         $ InfinitySparseMatrix
  #         $ InfinitySparseMatrix
  # Step 1:
  # List of 2
  #     $ List of 1
  #         $ InfinitySparseMatrix
  #     $ List of 2
  #         $ InfinitySparseMatrix
  #         $ InfinitySparseMatrix
  # Final result:
  # List of 3
  #     $ InfinitySparseMatrix
  #     $ InfinitySparseMatrix
  #     $ InfinitySparseMatrix
  flatten_list <- function(l) {
    unlist(lapply(l, function(x)
      if(!inherits(x, "list")) {
        # using `inherits` rather than `is.list` here because some objects
        # (specifically `data.frame`) return TRUE to `is.list`. This shouldn't
        # ever be an issue as input check above is restricted to types that
        # `is.list` returns FALSE for, but leaving with `inherits` for
        # future-proofing.
        list(x)
      }else {
        x
      }),
      recursive = FALSE)
  }


  # Convert all matrices to ISMs if they aren't already.
  mats <- lapply(mats, function(x) {
    if (is(x, "BlockedInfinitySparseMatrix")) {
      # Replace BISM with list of ISMs
      findSubproblems(x)
    } else if (inherits(x, "list")) {
      # If any entry in ... is a list,
      # 1) Convert all entries in that list to ISM while keeping BISM as BISM
      x <- lapply(x, .as.ism_or_bism)
      # 2) If we have any BISMs, split into list of ISMS
      x <- lapply(x, function(y) {
        if (is(y, "BlockedInfinitySparseMatrix")) {
          findSubproblems(y)
        } else {
          y
        }
      })
      # 3) pull list of lists into list
      flatten_list(x)
    } else {
      # This will error appropriately if some element in `mats` cannot be
      # converted to an ISM.
      .as.ism_or_bism(x)
    }
  })

  # If we were passed any BISMs, we have a list of lists of ISM, so flatten to a
  # single list.
  mats <- flatten_list(mats)

  # new row and column positions are based on current, incrementing by number of
  # rows/columns in all previous matrices.
  newcols <- lapply(mats, methods::slot, "cols")
  ncols <- vapply(lapply(mats, methods::slot, "dimension"), "[", 1, 2)
  for (i in 2:length(newcols)) {
    newcols[[i]] <- newcols[[i]] + sum(ncols[1:(i-1)])
  }
  newcols <- as.integer(do.call(c, newcols))

  newrows <- lapply(mats, methods::slot, "rows")
  nrows <- vapply(lapply(mats, methods::slot, "dimension"), "[", 1, 1)
  for (i in 2:length(newrows)) {
    newrows[[i]] <- newrows[[i]] + sum(nrows[1:(i-1)])
  }
  newrows <- as.integer(do.call(c, newrows))

  # names just get concatenated from all matrixes
  cnameslist <- lapply(mats, methods::slot, "colnames")
  newcolnames <- do.call(c, cnameslist)
  rnameslist <- lapply(mats, methods::slot, "rownames")
  newrownames <- do.call(c, rnameslist)
  if (any(duplicated(c(newcolnames)))) {
    if (force_unique_names == TRUE) {
      stop("Duplicated column or row names found.")
    }
    warning(paste("Duplicated column or row names found in matrices to be combined.\n",
                  "Renaming automatically to avoid issues; it is suggested to build",
                  "original matrices without this issue."))
    # If there are duplicated names, append a., b., c., etc to all names just to
    # ensure uniqueness.

    # To handle more than 26 entries, once we go through z., repeat with aa.,
    # ab., ac., ...., .zx, .zy, .zz, .aaa, .aab, .aac, etc. Supports up to 18278
    # separate entries
    doubleletters <- expand.grid(letters, letters)
    doubleletters <- paste0(doubleletters[, 2], doubleletters[, 1])
    tripleletters <- expand.grid(letters, doubleletters)
    tripleletters <- paste0(tripleletters[, 2], tripleletters[, 1])
    allletters <- c(letters, doubleletters, tripleletters)
    name_prefix <- paste0(allletters[seq_along(mats)], ".")
    cnameslist <- mapply(paste0, name_prefix, cnameslist, SIMPLIFY = FALSE)
    newcolnames <- do.call(c, cnameslist)
    rnameslist <- mapply(paste0, name_prefix, rnameslist, SIMPLIFY = FALSE)
    newrownames <- do.call(c, rnameslist)
  }

  # Adding all row dims and all column dims
  newdim <- as.integer(c(sum(vapply(lapply(mats, methods::slot, "dimension"), "[", 1, 1)),
                         sum(vapply(lapply(mats, methods::slot, "dimension"), "[", 1, 2))))

  # This needs to be much smarter, especially if any element is already a BISM
  groups <- as.factor(rep(0:(length(mats)-1), times =
                                      vapply(lapply(mats, slot, "colnames"), length, 1) +
                                        vapply(lapply(mats, slot, "rownames"), length, 1)))
  names(groups) <- do.call(c, Map(c, cnameslist, rnameslist))

  newdata <- do.call(c, mats)

  new("BlockedInfinitySparseMatrix",
      new("InfinitySparseMatrix",
          newdata,
          cols = newcols,
          rows = newrows,
          colnames = newcolnames,
          rownames = newrownames,
          dimension = newdim,
          call = call("match_on")),  # dummy call
      groups = groups)
}

##' @title Splits a BlockedInfinitySparseMatrix into a list of
##'   InfinitySparseMatrices
##' @param x a BlockedInfinitySparseMatrix
##' @param ... Ignored
##' @return A list of InfinitySparseMatrices
##' @export
as.list.BlockedInfinitySparseMatrix <- function(x, ...) {
  findSubproblems(x)
}

##' @export
as.list.InfinitySparseMatrix <- function(x, ...) {
  list(x)
}

##' @export
as.list.DenseMatrix <- function(x, ...) {
  list(as.InfinitySparseMatrix(x))
}

# (Internal) Converts item to ISM, but keeps as BISM if is BISM already.
.as.ism_or_bism <- function(x) {
  if (is(x, "BlockedInfinitySparseMatrix")) {
    return(x)
  }
  tryCatch(x <- as.InfinitySparseMatrix(x),
           error = function(e) {
             stop("Cannot convert object to InfinitySparseMatrices")
           })
  return(x)
}
