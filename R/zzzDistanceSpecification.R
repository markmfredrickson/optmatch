################################################################################
### Objects for specifying distances, methods to manipulate them
################################################################################

# The following objects are valid DistanceSpecifications, that is, they layout
# how treatment and control units are related.

setClassUnion("DistanceSpecification", c("matrix", "InfinitySparseMatrix"))

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


### subproblems: a DistanceSpecification may be split into smaller blocks
### returns false if there are no more subproblems, returns a list otherwise

setGeneric("subproblems", function(distances)
  standardGeneric("subproblems"))

# same method for matrices and ISMs
setMethod("subproblems", "DistanceSpecification", function(distances) FALSE)

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

  ftmp <- Filter(function(x) { length(x) > 0 }, tmp)

  if (length(ftmp) == 0) {
    return(FALSE) # no subproblems found  
  }

  return(ftmp)
})

# a helper to get all the subproblems from the implied tree of a DS
findSubproblems <- function(d) {
  res <- subproblems(d)
  
  if (is.list(res)) {
    return(do.call(c, lapply(res, findSubproblems)))
  }

  return(list(d))

} 

### Validating and error checking for DistanceSpecification objects
setGeneric("validDistanceSpecifcation", function(distance, stopOnProblem = TRUE)
  standardGeneric("validDistanceSpecifcation"))

# for now, only one method for all DistSpec objects (matrix, ISM, BISM)
setMethod("validDistanceSpecifcation", "DistanceSpecification", 
function(distance, stopOnProblem = TRUE) {
  valid <- TRUE # innocent until proven guilty

  valid <- valid & !is.null(dimnames(distance))

  if (stopOnProblem & !valid) {
    stop("Distance must have dimnames. Rows are treated, columns are control.")  
  }

  # for matrices we check with is.numeric
  valid <- valid & is.numeric(distance)

  # for ISMs, non-numeric input can be turned into a zero length vector
  # luckily for us, this should be true of matrices as well.
  valid <- valid & length(distance) > 0

  if (stopOnProblem & !valid) {
    stop("Distance must be numeric.")  
  }

  return(valid)
})


