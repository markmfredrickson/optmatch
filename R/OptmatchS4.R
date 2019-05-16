################################################################################
# Optmatch Class: the result of calling fullmatch()
################################################################################

#' Optmatch Class
#'
#' The \code{Optmatch} class describes the results of an optimal full matching
#' (using either \code{\link{fullmatch}} or \code{\link{pairmatch}}). Often, these objects can be treated as \code{factors}.

#' Elements of the underlying factor vector correspond to members of the treatment and control
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
#' process and its result, stored as slots:
#' @slot node.data Data frame containing various information about units involved in matches, including name, price, whether or not it was a treatment or control unit, and
#' an identifier indicating to which subproblem each unit belonged. Column names are "name", "price", contrast.group (1 means treatment group, 0 means control group), and group (an id associated with each subproblem), respectively.
#' @slot prob.data Data frame containing information about subproblems. In particular, the \code{exceedances} column gives
#' upper bounds, not necessarily sharp, for the amount by which the sum of
#' distances between matched units in the result of \code{fullmatch} exceeds
#' the least possible sum of distances between matched units in a feasible
#' solution to the matching problem given to \code{fullmatch}. Other information details specified parameters about restrictions on control units.
#' @slot names Character vector indicating the names of treatment and control units
#' @slot call Formula that was used to generate the object.
#' @slot subproblem Vector indicating different subproblems
#' @slot hashed.distance See \code{hashed.distance} function
#' @slot matched.distance See \code{matched.distances} function
#' @rdname optmatch
#' @name optmatch
#' @aliases optmatch-class
NA


Optmatch <- setClass("Optmatch", representation(node.data = "data.frame", prob.data = "data.frame"
, names = "character"
, call = "call", subproblem = "factor", hashed.distance = "character", matched.distance = "array"), contains = "factor")

####### Object Creation #########

#' (Internal) Create \code{optmatch} objects, the result of matching.
#'
#' This internal function is used to create the final output of the matching
#' functions (\code{\link{fullmatch}} and \code{\link{pairmatch}}). The
#' \code{Optmatch} object is an S4 object that inherits from a \code{factor}, but contains additional
#' information relating to the quality of the match in various slots.
#'
#' @param distance A \code{DistanceSpecificaton} object used to create the
#'  match.
#' @param solutions A list of the results of the matching, one \code{list(cells,maxErr)} object per subproblem.
#' @param call The call to \code{fullmatch} or \code{pairmatch} to be displayed later.
#' @param data An object from which \code{names} or \code{row.names} will
#'  provide the order of the items in the match. If no names are attached to this object, the contents will be used as names.
#' @return \code{Optmatch} object with slots described in the class description.
#' @seealso \code{\link{summary.optmatch}}
#' @keywords internal
makeOptmatch <- function(distance,
                         solutions,
                         call,
                         data = NULL)
{
  # pull out just the matching vectors

  matching <- lapply(solutions, function(x) { x$cells })

  treated <- rownames(distance)
  control <- colnames(distance)
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
    subproblems <- subproblems[optorder]
    names(optmatch.obj) <- names(subproblems) <- optorder
  }
  #class(optmatch.obj) <- c("optmatch", "factor")
  s4.opt <- new("Optmatch", optmatch.obj)

  if(!is.null(call))
  {
    s4.opt@call <- call
  }

  cg <- rep(NA, length(names(optmatch.obj)))
  cg[names(optmatch.obj) %in% treated] <- 1
  cg[names(optmatch.obj) %in% colnames(distance)] <- 0

  s4.opt@node.data <- assemble_node.data(solutions, treated, control)
  s4.opt@prob.data <- assemble_prob.data(solutions)

  if( nrow(s4.opt@node.data) > 0 )
  {
    if(any(as.character(na.omit((s4.opt@node.data$name[s4.opt@node.data$contrast.group]))) %in% control))
    {
      s4.opt@node.data[match(treated, s4.opt@node.data$name), c('contrast.group')] <- TRUE
      s4.opt@node.data[match(control, s4.opt@node.data$name), c('contrast.group')] <- FALSE
      #not sure why we need this block of code-- fixes the reversal of t/c units that occurs in ~ 472 of fullmatch.R in some situations
    }

  }
  s4.opt@subproblem <- subproblems
  return(s4.opt)
}

####### Subsetting and other manipulations #########

# I'm essentially deprecating the "drop" argument -- I think it makes these objects more complicated to work with, and I don't think keeping unused levels as part of the object is necessarily a desired feature.
# Unused levels will always be dropped now

subsetOptmatchIndexDrop <- function(x, i, j, drop)
{
  sub.opt.data <- factor(x@.Data, labels = x@levels)
  names(sub.opt.data) <- x@names
  sub.opt.data <- sub.opt.data[i]
  #drop the unused levels, as dictated by DROP
  sub.opt.data <- droplevels(sub.opt.data)
  sub.opt <- new("Optmatch", sub.opt.data)
  sub.opt@prob.data <- x@prob.data
  sub.opt@node.data <- x@node.data
  sub.opt@hashed.distance <- x@hashed.distance
  sub.opt@call <- x@call
  if(length(x@subproblem)) {
    sub.opt@subproblem <- x@subproblem[i]
  }
  else
  {
    sub.opt@subproblem <- x@subproblem
  }
  return(sub.opt)
}

subsetOptmatchnoIndexDrop <- function(x, i, j, drop)
{
  sub.opt.data <- droplevels(factor(x@.Data, labels = x@levels))
  names(sub.opt.data) <- x@names
  sub.opt.data <- sub.opt.data
  #drop the unused levels, as dictated by DROP
  #sub.opt.data <- droplevels(sub.opt.data)
  sub.opt <- new("Optmatch", sub.opt.data)
  sub.opt@prob.data <- x@prob.data
  sub.opt@node.data <- x@node.data
  sub.opt@hashed.distance <- x@hashed.distance
  sub.opt@call <- x@call
  if(length(x@subproblem)) {
    sub.opt@subproblem <- x@subproblem
  }
  else
  {
    sub.opt@subproblem <- x@subproblem
  }
  return(sub.opt)
  #not quite sure what's happening here
}

subsetOptmatchnoIndexNoDrop <- function(x, i, j, drop)
{
  return(x)
}
#' @export
setMethod("[", signature(x = "Optmatch", i = "ANY", drop = "missing", j = "missing"), definition = subsetOptmatchIndexDrop)

#' @export
setMethod("[", signature(x = "Optmatch", i = "ANY", drop = "logical", j = "missing"), definition = subsetOptmatchIndexDrop)

#' @export
setMethod("[", signature(x = "Optmatch", i = "missing", drop = "missing", j = "missing"), definition = subsetOptmatchnoIndexNoDrop)

#' @export
setMethod("[", signature(x = "Optmatch", i = "missing", drop = "logical", j = "missing"), definition = subsetOptmatchnoIndexNoDrop)

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
  if (!is(obj, "Optmatch")) {
    stop("Input must be an Optmatch object")
  }

  if (all(is.na(slot(obj, "prob.data")$omit.fraction))) {
    a <- slot(obj, "prob.data")$min.control
    b <- slot(obj, "prob.data")$max.control
    c <- slot(obj, "prob.data")$mean.control
    names(a) <- if(!is.null(a)) slot(obj, "prob.data")$group
    names(b) <- if(!is.null(b)) slot(obj, "prob.data")$group
    names(c) <- if(!is.null(c)) slot(obj, "prob.data")$group
    d <- slot(obj, "prob.data")$omit.fraction
    names(d) <- if(!is.null(d)) slot(obj, "prob.data")$group
    return(list("min.controls"=a, "max.controls"=b, "mean.controls"=c, "omit.fraction" = d))
  } else {
    a <- slot(obj, "prob.data")$min.control
    b <- slot(obj, "prob.data")$max.control
    c <- slot(obj, "prob.data")$omit.fraction
    names(a) <- if(!is.null(a)) slot(obj, "prob.data")$group
    names(b) <- if(!is.null(b)) slot(obj, "prob.data")$group
    names(c) <- if(!is.null(c)) slot(obj, "prob.data")$group
    d <- slot(obj, "prob.data")$mean.controls
    names(d) <- if(!is.null(d)) slot(obj, "prob.data")$group
    return(list("min.controls"=a, "max.controls"=b, "omit.fraction"=c, "mean.controls" = d))
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
  if (!is(obj, "Optmatch")) {
    stop("obj must be an Optmatch object")
  }
  if (!class(newdist) %in% c("BlockedInfinitySparseMatrix", "InfinitySparseMatrix", "DenseMatrix")) {
    stop("newdist must be a valid distance")
  }

  return(slot(obj, "hashed.distance") == dist_digest(newdist))
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
#' @export
compare_optmatch <- function(o1, o2) {

  if (length(setdiff(names(o1[!is.na(o1)]), names(o2[!is.na(o2)]))) > 0) {
    return(FALSE)
  }

  # Creates a list of the names of the members of each subgroup
  l1 <- lapply(levels(o1), function(x) sort(names(o1)[o1 == x]))
  l2 <- lapply(levels(o2), function(x) sort(names(o2)[o2 == x]))

  return(length(setdiff(l1,l2)) == 0)
}

#' Function converting an S4-style Optmatch object into an equivalent s3 optmatch object, if possible
#' THIS NEEDS TO BE UPDATED IN ACCORDANCE WITH CHANGES/UPDATES TO FUTURE NODE.DATA CHANGES.
#' @export
as.optmatch <- function(Optmatch)
{
  if(!is(Optmatch, "Optmatch"))
  {
    stop("object is not an Optmatch object")
  }
  else
  {
    opt.obj <- droplevels(factor(Optmatch@.Data, labels = Optmatch@levels))
    class(opt.obj) <- c("optmatch", "factor")
    attr(opt.obj, "node.data") <- Optmatch@node.data
    attr(opt.obj, "prob.data") <- Optmatch@prob.data
    attr(opt.obj, "call") <- Optmatch@call
    attr(opt.obj, "subproblem") <- Optmatch@subproblem
    attr(opt.obj, "hashed.distance") <- Optmatch@hashed.distance
    names(opt.obj) <- Optmatch@names
    t <- Optmatch@node.data[match(Optmatch@names, Optmatch@node.data$name),c("contrast.group") ]
    names(t) <- Optmatch@node.data[match(Optmatch@names, Optmatch@node.data$name),c("name") ]
    attr(opt.obj, "contrast.group") <- t

    return(opt.obj)
  }
}




