# Constant to control the maximum feasible (sub)problem
MAX_FEASIBLE <- 1e07 - 2

#' (Internal) Sets up the default values for maximum feasible problems
#
# @return NULL
setFeasibilityConstants <- function() {
  options("optmatch_max_problem_size" = MAX_FEASIBLE)    
}

#' (Internal) What is the maximum allowed problem size?
#'
#' To prevent users from starting excessively large matching problems, the
#' maximum problem size is limited by \code{options("optmatch_max_problem_size")}.
#' This function a quick helper to assist fetching this value as a scalar. If the
#' option isn't set, the function falls back to the default value, hard coded in
#' the \code{optmatch} package.
#'
#' @seealso \code{\link{options}}
#' @return logical
#' @examples getMaxProblemSize() > 1 & getMaxProblemSize() < 1e100
getMaxProblemSize <- function() {
  tmp <- options("optmatch_max_problem_size")[[1]]
  if (is.null(tmp[[1]])) {
    return(MAX_FEASIBLE)  
  }
  return(tmp[[1]])
}

#' Find the minimal exact match factors that will be feasible.
#'
#' The \code{\link{exactMatch}} function creates a smaller matching problem by
#' stratifying observations into smaller groups. For a problem that is larger
#' than #' maximum allowed size, #' \code{minExatMatch} provides a way to find the
#' smallest exact matching problem #' that will allow for matching.
#'
#' @param x The object for dispatching.
#' @param ... Additional arguments for methods.
#' @return DistanceSpecification As would be returned by \code{\link{exactMatch}}.
#' @export
setGeneric("minExactMatch", function(x, ...) standardGeneric("minExactMatch"))

#' The \code{formula} method takes a an argument of the form \code{Z ~ X1 + X2}, where
#' \code{Z} is indicates treatment or control status, and \code{X1} and \code{X2} are variables
#' can be converted to factors. Any additional arguments are passed to \code{\link{model.frame}}.
#' 
#' @param data A data.frame containing the variables in the formula \code{x}.
#’ @rdname minExactMatch-methods
setMethod("minExactMatch",
          signature = c("formula"),
function(x, ...) {
  parts <- as.character(x) # character vector of the form c("~", "Z", "X1 + X2")

  if (length(parts) != 3) {
    stop("Formula must be of the form Z ~ X1 + X2 + ...")  
  }

  lhs <- parts[2]
  rhs <- strsplit(parts[3], "\\+")[[1]] # ~ X1 + I(X2 == 2) =>  c("X1 ", " I(X2 == 2)")
  k <- length(rhs)

  for(i in 1:k) {
    fmla <- as.formula(paste(lhs, "~", paste(rhs[1:i], collapse = "+")),
                       env = attr(x, ".Environment")) 
  
    # generate a data.frame of treatment status and the blocking factor:
    z.b <- fmla2treatedblocking(fmla, ...)

    arcs <- aggregate(z.b$Z, by = list(z.b$B), function(grp) { sum(grp) * sum(1 - grp) })$x
    if (all(arcs < getMaxProblemSize())) {
      return(exactMatch(fmla, ...))  
    }
  }

  stop("Unable to create sufficiently small problem. Please provide more stratifying variables.")
})

#' (Internal) A helper function to turn formulas into treatment and blocking variables
#'
#' Given a function and any of the arguments noramally passed to model.frame,
#' this function will return a data.frame with two columns: a treatment indicator
#' and a blocking factor.
#'
#' @param x A formula
#' @param ... Arguments to be passed to model.frame (e.g. \code{data})
#' @return data.frame containing two columns: \code{Z} is a treatment indicator, 
#' \code{B} is a blocking factor
fmla2treatedblocking <- function(x, ...) {
  
  mf <- model.frame(x, ...)

  blocking <- interaction(mf[,-1])
  treatment <- toZ(mf[,1])
  
  df <- data.frame(Z = treatment, B = blocking)
  return(df)
}
