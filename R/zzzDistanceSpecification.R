################################################################################
### Objects for specifying distances, methods to manipulate them
################################################################################

# mdist.R provides some S3/S4 stuff for the optmatch.dlist class
#' @include mdist.R
NULL

### prepareMatching: DistanceSpecification -> arcs
### where arcs is a data.frame with 3 columns: control, treatment, distance

setGeneric("prepareMatching", function(distances)
  standardGeneric("prepareMatching"))

setMethod("prepareMatching", "matrix",
function(distances) {
  # note: treatment units are the rows, control units are the columns
  cs <- colnames(distances)
  rs <- rownames(distances)
  raw.distances <- as.vector(distances)

  control.reps <- length(raw.distances) / length(cs)
  treatment.reps <- length(raw.distances) / length(rs)

  controls <- rep(cs, each = control.reps)
  treatments <- rep(rs, treatment.reps)

  idx <- is.finite(raw.distances)

  # call factor() in order to drop unused levels
  tmp <- data.frame(control = factor(controls[idx]),
                    treated = factor(treatments[idx]),
                    distance = raw.distances[idx])

  return(tmp)
})

setMethod("prepareMatching", "InfinitySparseMatrix", function(distances) {
  if (is.null(distances@rownames) | is.null(distances@colnames)) {
    stop("Row and column names are required for matching.")
  }

  tmp <- data.frame(control = as.factor(distances@colnames[distances@cols]),
                    treated = as.factor(distances@rownames[distances@rows]),
                    distance = distances@.Data)

  return(tmp)

})

setMethod("prepareMatching", "optmatch.dlist",
function(distances) {
  Reduce(rbind, lapply(distances, function(grp) {
    prepareMatching(grp) # should be a matrix object
  }))
})


### subproblems: a DistanceSpecification may be split into smaller blocks
### returns false if there are no more subproblems, returns a list otherwise

setGeneric("subproblems", function(distances)
  standardGeneric("subproblems"))

setMethod("subproblems", "InfinitySparseMatrix", function(distances) FALSE)
setMethod("subproblems", "matrix", function(distances) FALSE)
setMethod("subproblems", "DenseMatrix", function(distances) FALSE)

# same method for matrices and ISMs
setMethod("subproblems", "BlockedInfinitySparseMatrix",
function(distances) {
  tmp <- lapply(levels(distances@groups), function(l) {
    members <- names(distances@groups[distances@groups == l])
    row.members <- which(distances@rownames %in% members)
    col.members <- which(distances@colnames %in% members)
    ridx <- distances@rows %in% row.members
    cidx <- distances@cols %in% col.members

    idx <- ridx & cidx

    makeInfinitySparseMatrix(distances[idx],
      rows = match(distances@rows[idx], row.members),
      cols = match(distances@cols[idx], col.members),
      rownames = distances@rownames[row.members],
      colnames = distances@colnames[col.members])
  })

  names(tmp) <- levels(distances@groups)

  ftmp <- Filter(function(x) { length(x) > 0 }, tmp)

  if (length(ftmp) == 0) {
    return(FALSE) # no subproblems found
  }

  return(ftmp)
})

setMethod("subproblems", "optmatch.dlist", function(distances) distances)

# a helper to get all the subproblems from the implied tree of a DS
findSubproblems <- function(d) {
  res <- subproblems(d)

  if (is.list(res)) {
    return(do.call(c, lapply(res, findSubproblems)))
  }

  return(list(d))

}

#' (Internal) Validate that objects are valid distance specifications.
#'
#' The functions \code{\link{fullmatch}} and \code{\link{pairmatch}}
#' create optimal matches of treated and control units given a
#' matrix (or similar representation) of distances between treated and
#' control units. These distance specifications must implement certain
#' generic functions. This function checks that all necessary methods
#' exist and the object can be used to specify distances in a way that
#' the matching functions can use.
#'
#' @param distance The object to test.
#' @param stopOnProblem If \code{TRUE} (default) the function will
#'   raise an error for invalid objects. Otherwise, it returns a
#'   logical.
#' @return logical
validDistanceSpecification <- function(distance, stopOnProblem = TRUE) {
  klass <- class(distance)[1] # in case there are multiple class names

  # we expect the following to be defined for the distance object
  # put any functions in this list that are called directly on distance
  methods <- c("dim", "dimnames", "prepareMatching", "subproblems")
  lapply(methods, function(m) {
    if (!hasMethod(m, klass)) {
      # skip the FUN = ... in the call stack
      if (stopOnProblem) {
        stop(paste("Invalid distance: object must have a", m, "method."), call. = F)
      } else {
        return(FALSE)
      }
    }
  })


  if (is.null(dimnames(distance))) {
    if (stopOnProblem) {
      stop("Invalid distance: object must have dimnames. Rows are treated, columns are control.")
    }
    return(FALSE)
  }

  sbprobs <- findSubproblems(distance)

  for(i in length(sbprobs)) {
    # for matrices we check with is.numeric

    if (!is.numeric(sbprobs[[i]])) {
      if (stopOnProblem) {
        stop("Invalid distance: object must be numeric.")
      }
      return(FALSE)
    }
  }

  if(!hasMethod(dim, class(distance)[1])) {
    if (stopOnProblem) {
      stop("Invalid distance: object must have dim method")
    }
    return(FALSE)
  }

  return(TRUE)
}
