#' Stratified Infinity Sparse Matrix
#'
#' Stratified Infinity Sparse Matrices are similar to Infinity Sparse Matrices, but they also keep track of the groups of units via an additional slot, \code{groups}
#' @slot groups factor vector containing groups, with unit names as labels, when possible
#' @template ISMslotsTemplate
#'
#' @seealso \code{\link{match_on}}, \code{\link{exactMatch}}, \code{\link{fullmatch}},  \code{\link{InfinitySparseMatrix-class}}
#' @author Mark M. Fredrickson
setClass("StratifiedInfinitySparseMatrix",
  slots = c(rowgroups = "factor", colgroups = "factor"),
  contains = "InfinitySparseMatrix")

sismOpsHandleGroups <- function(groups1, groups2) {
  if (is.null(groups1) || is.null(groups2)) {
    newGroups <- NULL
  } else if (length(groups1) > length(groups2)) {
    newGroups <- groups1
  } else {
    newGroups <- groups2
  }
  return(newGroups)
}

# in both of the next two methods I use callGeneric(as(...), ...)
# I would have prefered callNextMethod, but I kept getting errors,
# so I manually made the call to the parent class.
setMethod("Ops", signature(e1 = "StratifiedInfinitySparseMatrix",
                           e2 = "StratifiedInfinitySparseMatrix"),
function(e1, e2) {
  tmp <- callGeneric(as(e1, "InfinitySparseMatrix"), as(e2, "InfinitySparseMatrix"))
  tmp <- as(tmp, "StratifiedInfinitySparseMatrix")
  tmp@rowgroups <- sismOpsHandleGroups(e1@rowgroups, e2@rowgroups)
  tmp@colgroups <- sismOpsHandleGroups(e1@colgroups, e2@colgroups)
  return(tmp)
})

# the case where SISM is first is covered above
setMethod("Ops", signature(e1 = "InfinitySparseMatrix",
                           e2 = "StratifiedInfinitySparseMatrix"),
function(e1, e2) {
  tmp <- callGeneric(e1, as(e2, "InfinitySparseMatrix"))
  tmp <- as(tmp, "StratifiedInfinitySparseMatrix")
  tmp@rowgroups <- e2@rowgroups
  tmp@colgroups <- e2@colgroups
  return(tmp)
})

# SISMs need to maintain grouping info when getting flipped
#' @export
t.StratifiedInfinitySparseMatrix <- function(x) {
  tmp <- t.InfinitySparseMatrix(x)
  tmp <- as(tmp, "StratifiedInfinitySparseMatrix")
  tmp@rowgroups <- x@colgroups
  tmp@colgroups <- x@rowgroups
  return(tmp)
}

##' @export
##' @rdname cbindrbind
cbind.StratifiedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the stratified representation to a regular ISM and call the usual cbind method
  cbind.InfinitySparseMatrix(x, y, ...)
}

##' @export
##' @rdname cbindrbind
rbind.StratifiedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the stratified representation to a regular ISM and call the usual cbind method
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
subdim.StratifiedInfinitySparseMatrix <- function(x) {
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
sort.StratifiedInfinitySparseMatrix <- function(x,
                                             decreasing=FALSE,
                                             ...,
                                             byCol=FALSE) {
  y <- sort.InfinitySparseMatrix(x, decreasing=decreasing, byCol=byCol)
  attr(y, "rowgroups") <- attr(x, "rowgroups")
  attr(y, "colgroups") <- attr(x, "colgroups")
  class(y) <- class(x)
  y
}


#' Returns the number of eligible matches for the distance.
#'
#' This will return a list of the number of finite entries in a distance
#' matrix. If the distance has no subgroups, it will be a list of length 1. If
#' the distance has subgroups (i.e. \code{x} is an
#' \code{StratifiedInfinitySparseMatrix}, it will be a named list.)
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

#' @usage \method{num_eligible_matches}{StratifiedInfinitySparseMatrix}(x)
#' @rdname num_eligible_matches-methods
#' @export
num_eligible_matches.StratifiedInfinitySparseMatrix <- function(x) {
  out <- lapply(levels(x@groups), function(k) sum(is.finite(x[x@groups == k]@.Data)))
  names(out) <- levels(x@groups)
  out
}

##' Displays a StratifiedInfinitySparseMatrix
##'
##' Displays each block of the StratifiedInfinitySparseMatrix separately.
##' @param object An StratifiedInfinitySparseMatrix to print.
##' @return NULL
##' @export
setMethod("show", "StratifiedInfinitySparseMatrix", function(object) { show(findSubproblems(object)) })

##' This function generates a single block-diagonal distance matrix given
##' several distance matrices defined on subgroups.
##'
##' When you've generated several distances matrices on subgroups in your
##' analysis, you may wish to combine them into a single block-diagonal distance
##' matrix. The \code{dbind} function facilitates this.
##'
##' Any \code{StratifiedInfinitySparseMatrix} include in \code{...} will be broken
##' into individual \code{InfinitySparseMatrix} before being joined back
##' together. For example, if \code{b} is a \code{StratifiedInfinitySparseMatrix}
##' with 2 subgroups and \code{m} is a distance without subgroups, then
##' \code{dbind(b, m)} will be a \code{StratifiedInfinitySparseMatrix} with 3
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
##' single stratified distance using this function to ultimately obtain a single
##' matching set. Instead, take a look at the vignette
##' \code{vignette("matching-within-subgroups", package = "optmatch")} for
##' details on combining multiple matches.
##' @title Diagonally bind together subgroup-specific distances
##' @param ... Any number of distance objects which can be converted to
##'   \code{InfinitySparseMatrix}, such as class \code{matrix},
##'   \code{DenseMatrix}, \code{InfinitySparseMatrix}, or
##'   \code{StratifiedInfinitySparseMatrix}, or \code{list}s containing distance
##'   objects.
##' @param force_unique_names Default \code{FALSE}. When row or column names are
##'   not unique among all distances, if \code{FALSE}, throw a warning and
##'   rename all rows and columns to ensure unique names. If \code{TRUE}, error
##'   on non-unique names.
##' @return A \code{StratifiedInfinitySparseMatrix} containing a block-diagonal
##'   distance matrix. If only a single distance is passed to \code{dbind} and
##'   it is not already a \code{StratifiedInfinitySparseMatrix}, the result will be
##'   an \code{InfinitySparseMatrix} instead.
##' @importFrom methods slot
##' @export
##' @examples
##' data(nuclearplants)
##' m1 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 0),
##'                caliper = 1)
##' m2 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 1),
##'                caliper = 1.3)
##' stratified <- dbind(m1, m2)
##'
##' dists <- list(m1, m2)
##'
##' stratified2 <- dbind(dists)
##' identical(stratified, stratified2)
dbind <- function(..., force_unique_names = FALSE) {
  mats <- list(...)

  if (length(mats) == 1 & !inherits(mats[[1]], "list")) {
    # If passed a single distance, return it as ISM or SISM
    return(.as.ism_or_sism(mats[[1]]))
  }

  # Below, if one of the entries is a SISM, we'll end up with a list of lists,
  # where the SISM is replaced with a list of ISMs. The following flattens this
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
    if (is(x, "StratifiedInfinitySparseMatrix")) {
      # Replace SISM with list of ISMs
      findSubproblems(x)
    } else if (inherits(x, "list")) {
      # If any entry in ... is a list,
      # 1) Convert all entries in that list to ISM while keeping SISM as SISM
      x <- lapply(x, .as.ism_or_sism)
      # 2) If we have any SISMs, split into list of ISMS
      x <- lapply(x, function(y) {
        if (is(y, "StratifiedInfinitySparseMatrix")) {
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
      .as.ism_or_sism(x)
    }
  })

  # If we were passed any SISMs, we have a list of lists of ISM, so flatten to a
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

  # This needs to be much smarter, especially if any element is already a SISM
  groups <- as.factor(rep(0:(length(mats)-1), times =
                                      vapply(lapply(mats, slot, "colnames"), length, 1) +
                                        vapply(lapply(mats, slot, "rownames"), length, 1)))
  names(groups) <- do.call(c, Map(c, cnameslist, rnameslist))

  newdata <- do.call(c, mats)

  new("StratifiedInfinitySparseMatrix",
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

##' @title Splits a StratifiedInfinitySparseMatrix into a list of
##'   InfinitySparseMatrices
##' @param x a StratifiedInfinitySparseMatrix
##' @param ... Ignored
##' @return A list of InfinitySparseMatrices
##' @export
as.list.StratifiedInfinitySparseMatrix <- function(x, ...) {
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

# (Internal) Converts item to ISM, but keeps as SISM if is SISM already.
.as.ism_or_sism <- function(x) {
  if (is(x, "StratifiedInfinitySparseMatrix")) {
    return(x)
  }
  tryCatch(x <- as.InfinitySparseMatrix(x),
           error = function(e) {
             stop("Cannot convert object to InfinitySparseMatrices")
           })
  return(x)
}

sismSubsetGroups <- function(groups, whichVec) {
  newGroups <- NULL
  if (!is.null(groups)) {
    newGroups <- groups[whichVec]
  }
  return(newGroups)
}

##' This matches the syntax and semantics of
##' subset for matrices.
##'
##' @title Subsetting for StratifiedInfinitySparseMatrices
##' @param x StratifiedInfinitySparseMatrix to be subset or bound.
##' @param subset Logical expression indicating rows to keep.
##' @param select Logical expression indicating columns to keep.
##' @param ... Other arguments are ignored.
##' @return An StratifiedInfinitySparseMatrix with only the selected elements.
##' @rdname sism.subset
##' @export
subset.StratifiedInfinitySparseMatrix <- function(x, subset, select, ...) {
  # TODO
  # remove duplication in favor of callNextMethod
  # https://www.rdocumentation.org/packages/methods/versions/3.6.2/topics/callNextMethod
  # subset groups also...
  return(new(
    "StratifiedInfinitySparseMatrix",
    callNextMethod(x, subet, select),
    rowgroups = sismSubsetGroups(x@rowgroups, subset),
    colgroups = sismSubsetGroups(x@colgroups, select)
  ))
}


#' Stratified Infinity Sparse Matrix
#'
#' Stratified Infinity Sparse Matrices are similar to Infinity Sparse Matrices, but they also keep track of the groups of units via an additional slot, \code{groups}
#' @slot groups factor vector containing groups, with unit names as labels, when possible
#' @template ISMslotsTemplate
#'
#' @seealso \code{\link{match_on}}, \code{\link{exactMatch}}, \code{\link{fullmatch}},  \code{\link{InfinitySparseMatrix-class}}
#' @author Mark M. Fredrickson
setClass("StratifiedInfinitySparseMatrix",
  slots = c(rowgroups = "factor", colgroups = "factor"),
  contains = "InfinitySparseMatrix")

sismOpsHandleGroups <- function(groups1, groups2) {
  if (is.null(groups1) || is.null(groups2)) {
    newGroups <- NULL
  } else if (length(groups1) > length(groups2)) {
    newGroups <- groups1
  } else {
    newGroups <- groups2
  }
  return(newGroups)
}

# in both of the next two methods I use callGeneric(as(...), ...)
# I would have prefered callNextMethod, but I kept getting errors,
# so I manually made the call to the parent class.
setMethod("Ops", signature(e1 = "StratifiedInfinitySparseMatrix",
                           e2 = "StratifiedInfinitySparseMatrix"),
function(e1, e2) {
  tmp <- callGeneric(as(e1, "InfinitySparseMatrix"), as(e2, "InfinitySparseMatrix"))
  tmp <- as(tmp, "StratifiedInfinitySparseMatrix")
  tmp@rowgroups <- sismOpsHandleGroups(e1@rowgroups, e2@rowgroups)
  tmp@colgroups <- sismOpsHandleGroups(e1@colgroups, e2@colgroups)
  return(tmp)
})

# the case where SISM is first is covered above
setMethod("Ops", signature(e1 = "InfinitySparseMatrix",
                           e2 = "StratifiedInfinitySparseMatrix"),
function(e1, e2) {
  tmp <- callGeneric(e1, as(e2, "InfinitySparseMatrix"))
  tmp <- as(tmp, "StratifiedInfinitySparseMatrix")
  tmp@rowgroups <- e2@rowgroups
  tmp@colgroups <- e2@colgroups
  return(tmp)
})

# SISMs need to maintain grouping info when getting flipped
#' @export
t.StratifiedInfinitySparseMatrix <- function(x) {
  tmp <- t.InfinitySparseMatrix(x)
  tmp <- as(tmp, "StratifiedInfinitySparseMatrix")
  tmp@rowgroups <- x@colgroups
  tmp@colgroups <- x@rowgroups
  return(tmp)
}

##' @export
##' @rdname cbindrbind
cbind.StratifiedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the stratified representation to a regular ISM and call the usual cbind method
  cbind.InfinitySparseMatrix(x, y, ...)
}

##' @export
##' @rdname cbindrbind
rbind.StratifiedInfinitySparseMatrix <- function(x, y, ...) {
  # demote the stratified representation to a regular ISM and call the usual cbind method
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
subdim.StratifiedInfinitySparseMatrix <- function(x) {
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
sort.StratifiedInfinitySparseMatrix <- function(x,
                                             decreasing=FALSE,
                                             ...,
                                             byCol=FALSE) {
  y <- sort.InfinitySparseMatrix(x, decreasing=decreasing, byCol=byCol)
  attr(y, "rowgroups") <- attr(x, "rowgroups")
  attr(y, "colgroups") <- attr(x, "colgroups")
  class(y) <- class(x)
  y
}


#' Returns the number of eligible matches for the distance.
#'
#' This will return a list of the number of finite entries in a distance
#' matrix. If the distance has no subgroups, it will be a list of length 1. If
#' the distance has subgroups (i.e. \code{x} is an
#' \code{StratifiedInfinitySparseMatrix}, it will be a named list.)
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

#' @usage \method{num_eligible_matches}{StratifiedInfinitySparseMatrix}(x)
#' @rdname num_eligible_matches-methods
#' @export
num_eligible_matches.StratifiedInfinitySparseMatrix <- function(x) {
  out <- lapply(levels(x@groups), function(k) sum(is.finite(x[x@groups == k]@.Data)))
  names(out) <- levels(x@groups)
  out
}

##' Displays a StratifiedInfinitySparseMatrix
##'
##' Displays each block of the StratifiedInfinitySparseMatrix separately.
##' @param object An StratifiedInfinitySparseMatrix to print.
##' @return NULL
##' @export
setMethod("show", "StratifiedInfinitySparseMatrix", function(object) { show(findSubproblems(object)) })

##' This function generates a single block-diagonal distance matrix given
##' several distance matrices defined on subgroups.
##'
##' When you've generated several distances matrices on subgroups in your
##' analysis, you may wish to combine them into a single block-diagonal distance
##' matrix. The \code{dbind} function facilitates this.
##'
##' Any \code{StratifiedInfinitySparseMatrix} include in \code{...} will be broken
##' into individual \code{InfinitySparseMatrix} before being joined back
##' together. For example, if \code{b} is a \code{StratifiedInfinitySparseMatrix}
##' with 2 subgroups and \code{m} is a distance without subgroups, then
##' \code{dbind(b, m)} will be a \code{StratifiedInfinitySparseMatrix} with 3
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
##' single stratified distance using this function to ultimately obtain a single
##' matching set. Instead, take a look at the vignette
##' \code{vignette("matching-within-subgroups", package = "optmatch")} for
##' details on combining multiple matches.
##' @title Diagonally bind together subgroup-specific distances
##' @param ... Any number of distance objects which can be converted to
##'   \code{InfinitySparseMatrix}, such as class \code{matrix},
##'   \code{DenseMatrix}, \code{InfinitySparseMatrix}, or
##'   \code{StratifiedInfinitySparseMatrix}, or \code{list}s containing distance
##'   objects.
##' @param force_unique_names Default \code{FALSE}. When row or column names are
##'   not unique among all distances, if \code{FALSE}, throw a warning and
##'   rename all rows and columns to ensure unique names. If \code{TRUE}, error
##'   on non-unique names.
##' @return A \code{StratifiedInfinitySparseMatrix} containing a block-diagonal
##'   distance matrix. If only a single distance is passed to \code{dbind} and
##'   it is not already a \code{StratifiedInfinitySparseMatrix}, the result will be
##'   an \code{InfinitySparseMatrix} instead.
##' @importFrom methods slot
##' @export
##' @examples
##' data(nuclearplants)
##' m1 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 0),
##'                caliper = 1)
##' m2 <- match_on(pr ~ cost, data = subset(nuclearplants, pt == 1),
##'                caliper = 1.3)
##' stratified <- dbind(m1, m2)
##'
##' dists <- list(m1, m2)
##'
##' stratified2 <- dbind(dists)
##' identical(stratified, stratified2)
dbind <- function(..., force_unique_names = FALSE) {
  mats <- list(...)

  if (length(mats) == 1 & !inherits(mats[[1]], "list")) {
    # If passed a single distance, return it as ISM or SISM
    return(.as.ism_or_sism(mats[[1]]))
  }

  # Below, if one of the entries is a SISM, we'll end up with a list of lists,
  # where the SISM is replaced with a list of ISMs. The following flattens this
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
    if (is(x, "StratifiedInfinitySparseMatrix")) {
      # Replace SISM with list of ISMs
      findSubproblems(x)
    } else if (inherits(x, "list")) {
      # If any entry in ... is a list,
      # 1) Convert all entries in that list to ISM while keeping SISM as SISM
      x <- lapply(x, .as.ism_or_sism)
      # 2) If we have any SISMs, split into list of ISMS
      x <- lapply(x, function(y) {
        if (is(y, "StratifiedInfinitySparseMatrix")) {
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
      .as.ism_or_sism(x)
    }
  })

  # If we were passed any SISMs, we have a list of lists of ISM, so flatten to a
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

  # This needs to be much smarter, especially if any element is already a SISM
  groups <- as.factor(rep(0:(length(mats)-1), times =
                                      vapply(lapply(mats, slot, "colnames"), length, 1) +
                                        vapply(lapply(mats, slot, "rownames"), length, 1)))
  names(groups) <- do.call(c, Map(c, cnameslist, rnameslist))

  newdata <- do.call(c, mats)

  new("StratifiedInfinitySparseMatrix",
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

##' @title Splits a StratifiedInfinitySparseMatrix into a list of
##'   InfinitySparseMatrices
##' @param x a StratifiedInfinitySparseMatrix
##' @param ... Ignored
##' @return A list of InfinitySparseMatrices
##' @export
as.list.StratifiedInfinitySparseMatrix <- function(x, ...) {
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

# (Internal) Converts item to ISM, but keeps as SISM if is SISM already.
.as.ism_or_sism <- function(x) {
  if (is(x, "StratifiedInfinitySparseMatrix")) {
    return(x)
  }
  tryCatch(x <- as.InfinitySparseMatrix(x),
           error = function(e) {
             stop("Cannot convert object to InfinitySparseMatrices")
           })
  return(x)
}

sismSubsetGroups <- function(groups, whichVec) {
  newGroups <- NULL
  if (!is.null(groups)) {
    newGroups <- groups[whichVec]
  }
  return(newGroups)
}

##' This matches the syntax and semantics of
##' subset for matrices.
##'
##' @title Subsetting for StratifiedInfinitySparseMatrices
##' @param x StratifiedInfinitySparseMatrix to be subset or bound.
##' @param subset Logical expression indicating rows to keep.
##' @param select Logical expression indicating columns to keep.
##' @param ... Other arguments are ignored.
##' @return An StratifiedInfinitySparseMatrix with only the selected elements.
##' @rdname sism.subset
##' @export
subset.StratifiedInfinitySparseMatrix <- function(x, subset, select, ...) {
  # TODO
  # remove duplication in favor of callNextMethod
  # https://www.rdocumentation.org/packages/methods/versions/3.6.2/topics/callNextMethod
  # subset groups also...
  return(new(
    "StratifiedInfinitySparseMatrix",
    callNextMethod(x, subet, select),
    rowgroups = sismSubsetGroups(x@rowgroups, subset),
    colgroups = sismSubsetGroups(x@colgroups, select)
  ))
}
