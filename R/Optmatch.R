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

# S4 class compatability
setOldClass(c("optmatch", "factor"))

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

  attr(optmatch.obj, "subproblem") <- subproblems

  return(optmatch.obj)
}


####### Subsetting and other manipulations #########

#' @export
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
  if  (!is.null(attr(x, "subproblem"))) {
    sps <- attr(x, "subproblem")
    # converting to character to avoid a bug with subsetting
    # names get dropped with as.numeric
    nms <- names(sps)
    sps <- as.character(sps)
    names(sps) <- nms

    attr(y, "subproblem") <- "["(sps,...)
    attr(y, "subproblem") <- as.factor(attr(y, "subproblem"))
  }

  # Per issue #107, `matched.distances` are dropped when subsetting if they exist.
  attr(x, "matched.distances") <- NULL

  class(y) <- c("optmatch", "factor")

  return(y)
}

#' Returns the restrictions which were used to generate the match.
#'
#' If \code{mean.controls} was explicitly specified in the creation of the
#' optmatch object, it is returned; otherwise \code{omit.fraction} is given.
#'
#' Note that if the matching algorithm attempted to recover from initial
#' infeasible restrictions, the output from this function may not be the same as
#' the original function call.
#'
#' @title optmatch_restrictions
#' @param obj An optmatch object
#' @return A list of \code{min.controls}, \code{max.controls} and either
#' \code{omit.fraction} or \code{mean.controls}.
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
#' Note that the distance is hashed with its \code{call} set to
#' \code{NULL}. (This avoids issues where, for example, \code{match_on(Z~X,
#' data=d, caliper=NULL)} and \code{match_on(Z~X, data=d)} produce identical
#' matches (since the default argument to \code{caliper} is \code{NULL}) but
#' distinct calls.)
#' @param obj An optmatch object.
#' @param newdist A distance
#' @return Boolean whether the two distance specifications are identical.
#' @export
optmatch_same_distance <- function(obj, newdist) {
  if (!is(obj, "optmatch")) {
    stop("obj must be an optmatch object")
  }
  if (!class(newdist) %in% c("BlockedInfinitySparseMatrix", "InfinitySparseMatrix", "DenseMatrix")) {
    stop("newdist must be a valid distance")
  }

  return(attr(obj, "hashed.distance") == dist_digest(newdist))
}

#' Performs an update on an \code{optmatch} object.
#'
#' NB: THIS CODE IS CURRENTLY VERY MUCH ALPHA AND SOMEWHAT UNTESTED, ESPECIALLY CALLING \code{update} ON AN
#' OPTMATCH OBJECT CREATED WITHOUT AN EXPLICIT DISTANCE ARGUMENT.
#'
#' Note that passing \code{data} again is strongly recommended. A warning will be printed if the hash of the data used to generate the
#' \code{optmatch} object differs from the hash of the new \code{data}.
#' @param optmatch \code{Optmatch} object to update.
#' @param ... Additional arguments to the call, or arguments with changed values.
#' @param evaluate If true evaluate the new call else return the call.
#' @return An updated \code{optmatch} object.
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

  if (evaluate) {
    newmatch <- eval(call, parent.frame())
    if (attr(newmatch, "hashed.distance") != attr(optmatch, "hashed.distance")) {
      warning(paste("Distance given in update (", attr(newmatch, "hashed.distance"),
                    ") is different than distance used to generate fullmatch (",
                    attr(optmatch,"hashed.distance"), ").", sep=''))
    }
    newmatch
  } else call

}

#' Compares the equality of optmatch objects, ignoring attributes and group names.
#'
#' This checks the equality of two optmatch objects. The only bits that matter are unit names
#' and the grouping. Other bits such as attributes, group names, order, etc are ignored.
#'
#' The names of the units can differ on any unmatched units, e.g., units whose value in the optmatch
#' object is \code{NA}. If matched objects have differing names, this is automatically \code{FALSE}.
#'
#' Note this ignores the names of the subgroups. So four members in subgroups either
#' \code{c("a", "a", "b", "b")} or \code{c("b", "b", "a", "a")} would be identical to this call.
#' @param o1 First optmatch object.
#' @param o2 Second optmatch object.
#' @return TRUE if the two matches have the same memberships.
compare_optmatch <- function(o1, o2) {
  if (length(setdiff(names(o1[!is.na(o1)]), names(o2[!is.na(o2)]))) > 0) {
    return(FALSE)
  }

  # Creates a list of the names of the members of each subgroup
  l1 <- lapply(levels(o1), function(x) sort(names(o1)[o1 == x]))
  l2 <- lapply(levels(o2), function(x) sort(names(o2)[o2 == x]))

  return(length(setdiff(l1,l2)) == 0)
}
