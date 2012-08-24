#' Optimal 1:1 and 1:k matching
#' 
#' Given a treatment group, a larger control reservoir, and discrepancies
#' between each treatment and control unit, finds a pairing of treatment
#' units to controls that minimizes the sum of discrepancies.
#' 
#' This is a wrapper to \code{\link{fullmatch}}; see its documentation for
#' more information, especially on addiitonal arguments to pass.
#' 
#' If \code{remove.unmatchables} is \code{FALSE}, then if there are
#' unmatchable treated units then the matching as a whole will fail and
#' no units will be matched.  If \code{TRUE}, then this unit will be
#' removed and the function will attempt to match each of the other
#' treatment units.  (In this case matching can still fail, if there is
#' too much competition for certain controls; if you find yourself in
#' that situation you should consider full matching, which necessarily
#' finds a match for everyone with an eligible match somewhere.)  }
#' \value{Primarily, a named vector of class \code{c('optmatch',
#' 'factor')}.  Elements of this vector correspond to members of the
#' treatment and control groups in reference to which the matching
#' problem was posed, and are named accordingly; the names are taken from
#' the row and column names of \code{distance}.  Each element of the
#' vector is the concatenation of: (i) a character abbreviation of
#' \code{subclass.indices}, if that argument was given, or the string
#' '\code{m}' if it was not; (ii) the string \code{.}; and (iii) a
#' nonnegative integer or the string \code{NA}.  In this last place,
#' positive whole numbers indicate placement of the unit into a matched
#' set, a number beginning with zero indicates a unit that was not
#' matched, and \code{NA} indicates that all or part of the matching
#' problem given to \code{fullmatch} was found to be infeasible.
#' Secondarily, \code{fullmatch} returns various data about the matching
#' process and its result, stored as attributes of the named vector which
#' is its primary output.  In particular, the \code{exceedances}
#' attribute gives upper bounds, not necessarily sharp, for the amount by
#' which the sum of distances between matched units in the result of
#' \code{fullmatch} exceeds the least possible sum of distances between
#' matched units in a feasible solution to the matching problem given to
#' \code{fullmatch}.  (Such a bound is also printed by
#' \code{print.optmatch} and by \code{summary.optmatch}.) 
#' 
#' @param distance A matrix of nonnegative discrepancies, 
#' each indicating the permissibility and desirability of matching the unit 
#' corresponding to its row (a 'treatment') to the unit
#' corresponding to its column (a 'control'); or a list of such matrices
#' made using \code{\link{mdist}}.  Finite
#' discrepancies indicate permissible matches, with smaller
#' discrepancies indicating more desirable matches. 
#' Matrix \code{distance}, or the matrix elements of \code{distance},
#' must have row and column names. 
#' @param controls The number of controls to be matched to each treatment
#' @param remove.unmatchables Should treatment group members for which there are no eligible controls be removed prior to matching?
#' @param ... Additional arguments to pass to \code{\link{fullmatch}}.
#' It is strongly suggested that you pass a \code{data} argument.
#' It is
#' an error to pass \code{min.controls}, \code{max.controls},
#' or \code{omit.fraction} as \code{pairmatch} must set these values.
#' @return \code{optmatch}, a factor like matching indicator
#' @references
#' Hansen, B.B. and Klopfer, S.O. (2006), \sQuote{Optimal full matching
#' and related designs via network flows}, \emph{Journal of Computational
#' and Graphical Statistics}, \bold{15}, 609--627.   
#' 
#' @seealso \code{\link{matched}}, \code{\link{caliper}}, \code{\link{fullmatch}}} 
#' @example inst/examples/pairmatch.R
#' @keywords nonparametric optimize
#' @export
pairmatch <- function(distance, controls = 1, remove.unmatchables = FALSE, ...) {

  # Should this checking be pushed to fullmatch to avoid duplication?
  if (!is(distance, "DistanceSpecification")) {
    stop("argument \'distance\' must be a DistanceSpecification")
  }

  if (!all(floor(controls) == controls) | !all(controls > 0)) {
    stop("Minimum controls must be greater than treated units")  
  }

  # Check that max/min.controls and omit.fraction is not passed in ...
  dots <- names(match.call(expand.dots = TRUE))[-1] # first is always ""
  not.allowed <- c("min.controls", "max.controls", "omit.fraction")
  found <- not.allowed %in% dots
  if (any(found)) {
    stop("Invalid argument(s) to pairmatch: ", paste(not.allowed[found],
      collapse = ", "))   
  }

  subprobs <- findSubproblems(distance)

  if (length(controls) > 1 & !(length(subprobs) == length(controls))) {
    stop(paste("Controls argument must have same length as the number of subproblems (", 
      length(subprobs), ")", sep = ""))
  }

  omf <- mapply(controls, subprobs, FUN = function(control, prob) {
    # hard coding type based trimming for now. this should probably
    # be a DistanceSpecification method, e.g. finiteRows()
    if (remove.unmatchables) {
      if (inherits(prob, "matrix")) {
      # drop any rows that are entirely NA
      prob <- prob[apply(prob, 1, function(row) {
        any(is.finite(row)) }),]
      } else { 
        # assuming an InfinitySparseMatrix here      
        validrows <- which(1:(nrow(prob)) %in% prob@rows)    
        prob@dimension <- c(length(validrows), ncol(prob))
        prob@rownames <- prob@rownames[validrows]
      }
    }

    # a similar procedure is used to remove all control rows that
    # are unreachable

    if (inherits(prob, "matrix")) {
      # drop any rows that are entirely NA
      prob <- prob[, apply(prob, 2, function(row) {
        any(is.finite(row)) })]
    } else { 
        # assuming an InfinitySparseMatrix here      
        validcols <- which(1:(ncol(prob)) %in% prob@cols)    
        prob@dimension <- c(nrow(prob), length(validcols))
        prob@colnames <- prob@colnames[validcols]
    }

    nt <- nrow(prob)  
    nc <- ncol(prob)
    return((nc - control * nt)/nc)
  })

  if (any(omf<0)) {
    stop('not enough controls in some subclasses')
  }
  
  fullmatch(distance = distance, 
            min.controls = controls,
            max.controls = controls, 
            omit.fraction = omf,
            ...)
}

