################################################################################
# Optmatch Class: the result of calling fullmatch()
################################################################################

#' Optmatch Class
#'
#' The \code{optmatch} class describes the results of an optimal full matching
#' (using either \code{\link{fullmatch}} or \code{\link{pairmatch}}). For the
#' most part, these objects can be treated as \code{factors}.
#'
#' \code{optmatch} objects descend from \code{factor}.
#' Elements of this vector correspond to members of the treatment and control
#' groups in reference to which the matching problem was posed, and are named
#' accordingly; the names are taken from the row and column names of
#' \code{distance}.  Each element of the vector is either \code{NA}, indicating
#' unavailability of any suitable matches for that element, or the
#' concatenation of: (i) a character abbreviation of the name of the subclass
#' (as encoded using \code{\link{exactMatch}}) (ii) the string \code{.}; and
#' (iii) a non-negative integer.  In this last place, positive whole numbers
#' indicate placement of the unit into a matched set and \code{NA} indicates
#' that all or part of the matching problem given to \code{fullmatch} was found
#' to be infeasible.  The functions \code{\link{matched}},
#' \code{\link{unmatched}}, and \code{\link{matchfailed}} distinguish these
#' scenarios.
#'
#' Secondarily, \code{fullmatch} returns various data about the matching
#' process and its result, stored as attributes of the named vector which is
#' its primary output.  In particular, the \code{exceedances} attribute gives
#' upper bounds, not necessarily sharp, for the amount by which the sum of
#' distances between matched units in the result of \code{fullmatch} exceeds
#' the least possible sum of distances between matched units in a feasible
#' solution to the matching problem given to \code{fullmatch}.  (Such a bound
#' is also printed by \code{print.optmatch} and \code{summary.optmatch}.)
#' @rdname optmatch
#' @name optmatch
#' @aliases optmatch-class
NA

####### Object Creation #########

#' (Internal) Create \code{optmatch} objects, the result of matching.
#'
#' This internal function is used to create the final output of the matching
#' functions (\code{\link{fullmatch}} and \code{\link{pairmatch}}). The
#' \code{optmatch} object descends from a \code{factor}, but contains additional
#' information relating to the quality of the match.
#'
#' @param distance A \code{DistanceSpecificaton} object used to create the
#'  match.
#' @param solutions A list of the results of the matching, one \code{list(cells,maxErr)} object per subproblem.
#' @param call The call to \code{fullmatch} or \code{pairmatch} to be displayed later.
#' @param data An object from which \code{names} or \code{row.names} will
#'  provide the order of the items in the match. If no names are attached to this object, the contents will be used as names.
#' @return \code{optmatch} object
#'
#' @seealso \code{\link{summary.optmatch}}
makeOptmatch <- function(distance,
                         solutions,
                         call,
                         data = NULL)
{
  # pull out just the matching vectors
  matching <- lapply(solutions, function(x) { x$cells })

  treated <- rownames(distance)

  grpnames <- names(matching)
  if (is.null(grpnames)) {
    grpnames <- 1:(length(matching))
  }

  optmatch.obj <- Reduce(mapply(function(label, groups) {
        tmp <- groups
        tmp[!is.na(groups)] <- paste(label, groups[!is.na(groups)],
          sep = ".")
        return(tmp)
        }, grpnames, matching), f = c)

  optmatch.obj <- as.factor(optmatch.obj)
  subproblems <- as.factor(unlist(mapply(function(label, group) { rep(label, length(group)) }, grpnames, matching)))
  names(optmatch.obj) <- names(subproblems) <- unlist(sapply(matching, names))


  # we try to get the order as row names, straight names, and finally from the
  # value of the data argument.
  optorder <- NULL
  if(!is.null(data)) {
    optorder <- row.names(data)

    if (is.null(optorder)) {
      optorder <- names(data)
    }

    if (is.null(optorder) & is.vector(data)) {
      optorder <- as.character(data)
    }

    if (is.null(optorder)) {
      # if we are here, the user tried to pass data, but we couldn't get names
      warning("Unable to find appropriate names in 'data' argument.")
    }
  }

  if (!is.null(optorder)) {
    optmatch.obj <- optmatch.obj[optorder]
    names(optmatch.obj) <- optorder
    subproblems <- subproblems[optorder]
  }

  class(optmatch.obj) <- c("optmatch", "factor")

  tmp <- sapply(solutions, function(x) { x$err })
  names(tmp) <- grpnames
  attr(optmatch.obj, "exceedances") <- tmp

  attr(optmatch.obj, "call") <- call

  attr(optmatch.obj, "contrast.group") <- names(optmatch.obj) %in% treated ### WHAT IS INROW?
  # TODO TURN ON WHEN MATCHED DISTANCES IS UPDATED
  attr(optmatch.obj, "matched.distances") <- matched.distances(optmatch.obj, distance)

  attr(optmatch.obj, "subproblem") <- subproblems

  return(optmatch.obj)
}


####### Subsetting and other manipulations #########

"[.optmatch" <-
  function(x, ..., drop=FALSE)
{
  y <- NextMethod("[")
  if  (!is.null(attr(x, "contrast.group"))) {
    cgs <- attr(x, "contrast.group")
    names(cgs) <- names(x)

    attr(y,"contrast.group") <- "["(cgs,...)
    names(attr(y, "contrast.group")) <-  NULL
  }

  ### The following is something of a kluge.  It would make more sense
  ### to remove matched distances that have been removed from the optmatch
  ### vector, but doing that is not straightforward, since the distances don't
  ### straighforwardly line up with the observations.  At present (version 0.6),
  ### the matched.distances attribute is only used in summary.optmatch;
  ### I have inserted code there to compensate for non-subsetting of the
  ### matched distances attribute in the case where matching has failed in some
  ### subclasses.
  if (!is.null(attr(x, "matched.distances"))) {
    attr(y, "matched.distances") <- attr(x, "matched.distances")
  }

  class(y) <- c("optmatch", "factor")

  return(y)
}

#' Returns the restrictions which were used to generate the match.
#'
#' If \code{mean.controls} was explicitly specified in the creation of the
#' optmatch object, it is returned; otherwise \code{omit.fraction} is given.
#'
#' Note that if the matching algorithm attempted to recover from initial
#' infeasible restrictions, the output from this function is likely different
#' from the original function call.
#'
#' @title optmatch_restrictions
#' @param obj An optmatch object
#' @return A list of \code{min.controls}, \code{max.controls} and either \code{omit.fraction} or \code{mean.controls}.
#' @author Josh Errickson
#' @export
optmatch_restrictions <- function(obj) {
  if (!is(obj, "optmatch")) {
    stop("Input must be an optmatch object")
  }
  if (is.null(attr(obj, "omit.fraction"))) {
    return(list("min.controls"=attr(obj, "min.controls"), "max.controls"=attr(obj, "max.controls"), "mean.controls"=attr(obj, "mean.controls")))
  } else {
    return(list("min.controls"=attr(obj, "min.controls"), "max.controls"=attr(obj, "max.controls"), "omit.fraction"=attr(obj, "omit.fraction")))
  }
}

#' Checks if the distance \code{newdist} is identical to the distance used to
#' generate the optmatch object \code{obj}.
#'
#' To save space, optmatch objects merely store a hash of the distance matrix
#' instead of the original object. This checks if the hash of \code{newdist} is
#' identical to the hash currently saved in \code{obj}.
#'
#' @param obj An optmatch object.
#' @param newdist A distance
#' @return Boolean whether the two distance specifications are identical.
#' @author Josh Errickson
#' @import digest
#' @export
optmatch_same_distance <- function(obj, newdist) {
  if (!is(obj, "optmatch")) {
    stop("obj must be an optmatch object")
  }
  if (!class(newdist) %in% c("BlockedInfinitySparseMatrix", "InfinitySparseMatrix", "DenseMatrix")) {
    stop("newdist must be a valid distance")
  }

  return(attr(obj, "hashed.distance") == digest(newdist))
}

#' Performs an update on an \code{optmatch} object.
#'
#' Note that passing \code{distance} again is strongly recommended. A warning
#' will be printed if the hash of the distance used to generate the
#' \code{optmatch} object differs from the hash of the new \code{distance}.
#'
#' @param optmatch \code{Optmatch} object to update.
#' @param ... Additional arguments to the call, or arguments with changed values.
#' @param evaluate If true evaluate the new call eslse return the call.
#' @return An updated \code{optmatch} object.
#' @author Josh Errickson
#' @import digest
#' @export
update.optmatch <- function(optmatch, ..., evaluate = TRUE) {
  if (is.null(call <- attr(optmatch, "call")))
    stop("optmatch must have a call attribute")
  extras <- match.call(expand.dots = FALSE)$...

  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  # check if distance in the new call is identical to distance in the original object
  newdigest <- digest(eval(parse(text=call[2]), envir=parent.frame()))
  if(newdigest != attr(optmatch, "hashed.distance")) {
    warning(paste("Distance given in update (", newdigest, ") is different than distance used to generate fullmatch (",
                  attr(optmatch,"hashed.distance"), ").", sep=''))
  }

  if (evaluate)
    eval(call, parent.frame())
  else call
}
