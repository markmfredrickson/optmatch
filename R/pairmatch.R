#' Optimal 1:1 and 1:k matching
#'
#' Given a treatment group, a larger control reservoir, and a method for creating
#' discrepancies between each treatment and control unit (or optionally an
#' already created such discrepancy matrix), finds a pairing of treatment units
#' to controls that minimizes the sum of discrepancies.
#'
#' This is a wrapper to \code{\link{fullmatch}}; see its documentation for more
#' information, especially on additional arguments to pass, additional discussion
#' of valid input for parameter \code{x}, and feasibility recovery.
#'
#' If \code{remove.unmatchables} is \code{FALSE}, then if there are unmatchable
#' treated units then the matching as a whole will fail and no units will be
#' matched.  If \code{TRUE}, then this unit will be removed and the function will
#' attempt to match each of the other treatment units.  As of version 0.9-8,
#' if there are fewer matchable treated units than matchable controls then
#' \code{pairmatch} will attempt to place each into a matched pair each of the
#' matchable controls and a strict subset of the matchable treated units.
#' (Previously matching would have failed for subclasses of this structure.)
#'
#' Matching can still fail,
#' even with \code{remove.unmatchables} set to \code{TRUE},
#' if there is too much competition for certain controls; if you
#' find yourself in that situation you should consider full matching, which
#' necessarily finds a match for everyone with an eligible match somewhere.
#'
#' The units of the \code{optmatch} object returned correspond to members of the
#' treatment and control groups in reference to which the matching problem was
#' posed, and are named accordingly; the names are taken from the row and column
#' names of \code{distance} (with possible additions from the optional
#' \code{data} argument).  Each element of the vector is the concatenation of:
#' (i) a character abbreviation of \code{subclass.indices}, if that argument was
#' given, or the string '\code{m}' if it was not; (ii) the string \code{.}; and
#' (iii) a non-negative integer. Unmatched units have \code{NA} entries.
#' Secondarily, \code{fullmatch} returns various data about the matching process
#' and its result, stored as attributes of the named vector which is its primary
#' output.  In particular, the \code{exceedances} attribute gives upper bounds,
#' not necessarily sharp, for the amount by which the sum of distances between
#' matched units in the result of \code{fullmatch} exceeds the least possible sum
#' of distances between matched units in a feasible solution to the matching
#' problem given to \code{fullmatch}.  (Such a bound is also printed by
#' \code{print.optmatch} and by \code{summary.optmatch}.)
#'
#' @param x Any valid input to \code{match_on}. If \code{x} is a numeric vector,
#' there must also be passed a vector \code{z} indicating grouping. Both vectors
#' must be named.
#'
#' Alternatively, a precomputed distance may be entered.
#' @param controls The number of controls to be matched to each treatment
#' @param data Optional data set.
#' @param remove.unmatchables Should treatment group members for which there are
#' no eligible controls be removed prior to matching?
#' @param ... Additional arguments to pass to \code{\link{match_on}}
#' (e.g. \code{within})) or to \code{\link{fullmatch}} (e.g. \code{tol}).
#' It is an error to pass \code{min.controls},
#' \code{max.controls}, \code{mean.controls} or \code{omit.fraction} as
#' \code{pairmatch} must set these values.
#' @return A \code{\link{optmatch}} object (\code{factor}) indicating matched groups.
#' @references
#' Hansen, B.B. and Klopfer, S.O. (2006), \sQuote{Optimal full matching
#' and related designs via network flows}, \emph{Journal of Computational
#' and Graphical Statistics}, \bold{15}, 609--627.
#'
#' @seealso \code{\link{matched}}, \code{\link{caliper}}, \code{\link{fullmatch}}
#' @example inst/examples/pairmatch.R
#' @keywords nonparametric optimize
#' @export
pairmatch <- function(x,
                      controls = 1,
                      data = NULL,
                      remove.unmatchables = FALSE,
                      ...) {

  # if x does not exist then print helpful error msg
  x_str <- deparse(substitute(x))
  data_str <- deparse(substitute(data))
  tryCatch(x, error = function(e) {
    stop(missing_x_msg(x_str, data_str, ...))
  })

  # Check that max/min.controls and omit.fraction is not passed in ...
  dots <- names(match.call(expand.dots = TRUE))[-1] # first is always ""
  not.allowed <- c("min.controls", "max.controls", "mean.controls", "omit.fraction")
  found <- not.allowed %in% dots
  if (any(found)) {
    stop("Invalid argument(s) to pairmatch: ", paste(not.allowed[found],
      collapse = ", "))
  }

  if (is(x, "optmatch.dlist")) {
    warning("The use of 'optmatch.dlist' objects created by 'mdist()' is deprecated.\nPlease use 'match_on()' instead.")
  }

  UseMethod("pairmatch")
}

#' @export
pairmatch.default <- function(x,
                      controls = 1,
                      data = NULL,
                      remove.unmatchables = FALSE,
                      within = NULL,
                      ...) {

  if (!inherits(x, gsub("match_on.","",methods("match_on")))) {
    stop("Invalid input, must be a potential argument to match_on")
  }

  mfd <- if (!is.null(data)) {
    model.frame(data, na.action=na.pass)
  } else {
    if (inherits(x, "function")) {
      stop("A data argument must be given when passing a function")
    }
    model.frame(x, na.action=na.pass)
  }
  if (!is(mfd, "data.frame")) {
    stop("Please pass data argument")
  }
  m <- match_on(x, within=within, data=mfd, ...)
  out <- pairmatch(m,
                   controls=controls,
                   data=mfd,
                   remove.unmatchables=remove.unmatchables,
                   ...)
  attr(out, "call") <- match.call()
  out
}

#' @export
pairmatch.numeric <- function(x,
                      controls = 1,
                      data = NULL,
                      remove.unmatchables = FALSE,
                      z,
                      within = NULL,
                      ...) {

  m <- match_on(x, within=within, z=z, ...)
  out <- pairmatch(m,
                   controls=controls,
                   data=data,
                   remove.unmatchables=remove.unmatchables,
                   ...)
  attr(out, "call") <- match.call()
  out
}

#' @export
pairmatch.matrix <- function(x,
                             controls = 1,
                             data = NULL,
                             remove.unmatchables = FALSE,
                             within = NULL,
                             ...) {

  validDistanceSpecification(x) # will stop() on error

  if (!all(floor(controls) == controls) | !all(controls > 0)) {
    stop("Minimum controls must be greater than treated units")
  }

  subprobs <- findSubproblems(x)

  if (length(controls) > 1 & !(length(subprobs) == length(controls))) {
    stop(paste("Controls argument must have same length as the number of subproblems (",
      length(subprobs), ")", sep = ""))
  }

    if (!is.null(within)) warning("Ignoring non-null 'within' argument. When using 'pairmatch' with\n pre-formed distances, please combine them using '+'.")

  get_omf <- function(control, prob) {
    # hard coding type based trimming for now. this should probably
    # be a DistanceSpecification method, e.g. finiteRows()
    if (remove.unmatchables) {
      if (inherits(prob, "matrix")) {
          # drop any rows that are entirely NA
          prob <- prob[apply(prob, 1, function(row) {
              any(is.finite(row)) }),, drop=FALSE]
          # Now do the same for columns -- but only if
          # there are one or more rows left.
          # (Otherwise subsequent `apply()` quits.)
          if (nrow(prob)) {
              prob <- prob[,apply(prob, 2, function(col) {
                  any(is.finite(col)) }), drop=FALSE]
              }
      } else {
        # assuming an InfinitySparseMatrix here
        validrows <- which(1:(nrow(prob)) %in% prob@rows)
        prob@rownames <- prob@rownames[validrows]
          if (length(validrows)) {
              validcols <- which(1:(ncol(prob)) %in% prob@cols)
              prob@colnames <- prob@colnames[validcols]
              prob@dimension <- c(length(validrows), length(validcols))
              } else {
                  prob@dimension <- c(length(validrows), ncol(prob))
              }
      }
    }

    treatment_group_n <- nrow(prob)
    control_group_n <- ncol(prob)
    control_group_overage <- control_group_n - control * treatment_group_n
    treatment_group_overage <- treatment_group_n - control_group_n/control
      return(ifelse(control_group_overage>=0,
                    control_group_overage/control_group_n,
                    -1*treatment_group_overage/treatment_group_n))
  }

  omf <- mapply(controls, subprobs, FUN = get_omf)


  if(!remove.unmatchables) {
    saveopt <- options()$fullmatch_try_recovery
    options("fullmatch_try_recovery" = FALSE)
  }
  out <- fullmatch(x = x,
            min.controls = controls,
            max.controls = controls,
            omit.fraction = omf,
            data = data,
            ...)
  if(!remove.unmatchables) {
    options("fullmatch_try_recovery" = saveopt)
  }
  attr(out, "call") <- match.call()
  return(out)
}

#' @export
pairmatch.optmatch.dlist <- pairmatch.matrix
#' @export
pairmatch.InfinitySparseMatrix <- pairmatch.matrix
#' @export
pairmatch.BlockedInfinitySparseMatrix <- pairmatch.matrix

#' @aliases pairmatch
#' @rdname pairmatch
#' @export
pair <- pairmatch
