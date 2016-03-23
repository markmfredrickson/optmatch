################################################################################
# exactMatch: create InfinitySparseMatrices from factors
################################################################################

#' An exact match is one based on a factor. Within a level, all
#' observations are allowed to be matched. An exact match can be
#' combined with another distance matrix to create a set of matching
#' subproblems.
#'
#' \code{exactMatch} creates a block diagonal matrix of 0s and
#' \code{Inf}s. The pairs with 0 entries are within the same level of
#' the factor and legitimate matches.  \code{Inf} indicates units in
#' different levels. \code{exactMatch} replaces the
#' \code{structure.fmla} argument to several functions in previous
#' versions of optmatch.  For the \code{factor} method, the two
#' vectors \code{x} and \code{treatment} must be the same length. The
#' vector \code{x} is interpreted as indicating the grouping factors
#' for the data, and the vector \code{treatment} indicates whether a
#' unit is in the treatment or control groups.  At least one of these
#' two vectors must have names.  For the \code{formula} method, the
#' \code{data} argument may be omitted, in which case the method
#' attempts to find the variables in the environment from which the
#' function was called. This behavior, and the arguments \code{subset}
#' and \code{na.action}, mimics the behavior of \code{\link{lm}}.
#'
#' @title Generate an exact matching set of subproblems.
#' @author Mark M. Fredrickson
#'
#' @keywords nonparametric
#'
#' @param x A factor vector or a formula, used to select method.
#' @param treatment A vector the same length as \code{x} that can be
#'   coerced to a two level factor (e.g. a vector of 1s and 0s or a
#'   logical vector).
#' @param data A \code{data.frame} or \code{matrix} that contains the
#'   variables used in the formula \code{x}.
#' @param subset an optional vector specifying a subset of
#'   observations to be used
#' @param na.action A function which indicates what should happen when
#'   the data contain `NA's
#' @param ... Additional arguments for methods.
#'
#' @return A matrix like object, which is suitable to be given as
#'   \code{distance} argument to \code{\link{fullmatch}} or
#'   \code{\link{pairmatch}}. The exact match will be only zeros and
#'   \code{Inf} values, indicating a possible match or no possible
#'   match, respectively. It can be added to a another distance matrix
#'   to create a subclassed matching problem.
#'
#' @seealso \code{\link{caliper}}, \code{\link{antiExactMatch}},
#'   \code{\link{match_on}}, \code{\link{fullmatch}},
#'   \code{\link{pairmatch}}
#'
#' @export
#' @docType methods
#' @examples
#'
#' data(nuclearplants)
#'
#' ### First generate a standard propensity score
#' ppty <- glm(pr~.-(pr+cost), family = binomial(), data = nuclearplants)
#' ppty.distances <- match_on(ppty)
#'
#' ### Only allow matches within the partial turn key plants
#' pt.em <- exactMatch(pr ~ pt, data = nuclearplants)
#' as.matrix(pt.em)
#'
#' ### Blunt matches:
#' match.pt.em <- fullmatch(pt.em)
#' print(match.pt.em, grouped = TRUE)
#'
#' ### Combine the propensity scores with the subclasses:
#' match.ppty.em <- fullmatch(ppty.distances + pt.em)
#' print(match.ppty.em, grouped = TRUE)
#'
#' @rdname exactMatch
setGeneric("exactMatch",
  def = function(x, ...) {
    tmp <- standardGeneric("exactMatch")
    tmp@call <- match.call()
    return(tmp)
})

#' @export
#' @rdname exactMatch
setMethod(exactMatch, "vector", function(x, treatment) {
  if (length(x) != length(treatment)) {
    stop("Splitting vector and treatment vector must be same length")
  }

  # ham-handed way of saying, use x's names or use treatments's name
  # which ever is not null
  nms <- names(x)
  if (is.null(nms) & is.null(names(treatment))) {
    stop("Blocking or treatment factor must have names")
  } else {
    if(is.null(nms)) {
      nms <- names(treatment)
    }
  }

  # defensive programming
  x <- as.factor(x)
  treatment <- toZ(treatment)

  # the upper level is the treatment condition
  xT <- x[treatment]
  xC <- x[!treatment]

  csForTs <- lapply(xT, function(t) {
    which(t == xC)
  })

  cols <- unlist(csForTs)
  tmp <- sapply(csForTs, length)
  rows <- rep(1:(length(csForTs)), tmp)

  rns <- nms[treatment]
  cns <- nms[!treatment]

  tmp <- makeInfinitySparseMatrix(rep(0, length(rows)), cols = cols, rows =
    rows, rownames = rns, colnames = cns)

  tmp <- as(tmp, "BlockedInfinitySparseMatrix")
  tmp@groups <- x
  names(tmp@groups) <- nms
  return(tmp)
})

#' @export
#' @rdname exactMatch
setMethod(exactMatch, "formula", function(x, data = NULL, subset = NULL, na.action = NULL, ...) {
  # lifted pretty much verbatim from lm()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- x
  mf[[1L]] <- as.name("model.frame")

  mf <- eval(mf, parent.frame())

  blocking <- interaction(mf[,-1])
  treatment <- mf[,1]

  names(blocking) <- rownames(mf)
  names(treatment) <- rownames(mf)

  # formula is expected to be Z ~ B, where b is the blocking factor
  # and Z is treatment, Z ~ B1 + B2 ... is also allowed
  exactMatch(blocking, treatment) # use the factor based method
})

#' Specify a matching problem where units in a common factor cannot be matched.
#'
#' This function builds a distance specification where treated units
#' are infinitely far away from control units that share the same
#' level of a given factor variable. This can be useful for ensuring
#' that matched groups come from qualitatively different groups.
#'
#' The \code{\link{exactMatch}} function provides a way of specifying
#' a matching problem where only units within a factor level may be
#' matched. This function provides the reverse scenario: a matching
#' problem in which only units across factor levels are permitted to
#' match. Like \code{\link{exactMatch}}, the results of this function will
#' most often be used as a \code{within} argument to
#' \code{\link{match_on}} or another distance specification creation
#' function to limit the scope of the final distance specification
#' (i.e., disallowing any match between units with the same value on
#' the factor variable \code{x}).
#'
#' @param x A factor across which matches should be allowed.
#' @param z A treatment indicator factor (a numeric vector of 1 and 0,
#' a logical vector, or a 2 level factor).
#' @return A distance specification that encodes the across factor level constraint.
#' @seealso \code{\link{exactMatch}}, \code{\link{match_on}}, \code{\link{caliper}}, \code{\link{fullmatch}}, \code{\link{pairmatch}}
#' @export
#' @example inst/examples/antiExactMatch.R
antiExactMatch <- function(x, z) {
  z <- toZ(z)
  x <- as.factor(x)

  if (is.null(names(x)) && is.null(names(z))) {
    stop("Either 'x' or 'z' must have names")
  }

  nms <- names(x)
  if (is.null(nms)) {
    nms <- names(z)
  }

  controlnms <- nms[!z]
  treatednms <- nms[z]

  cid <- 1:length(controlnms)
  tid <- 1:length(treatednms)

  names(cid) <- controlnms
  names(tid) <- treatednms

  rowcols <- data.frame(rows = vector("numeric"), cols = vector("numeric"))

  for (l in levels(x)) {
    idx <- x == l
    in.treated      <- tid[nms[z & idx]]
    across.controls <- cid[nms[!z & !idx]]
    rowcols <- rbind(rowcols, expand.grid(rows = in.treated, cols = across.controls))
  }

  ret <- makeInfinitySparseMatrix(rep(0, dim(rowcols)[1]),
                                  rows = rowcols$rows,
                                  cols = rowcols$cols,
                                  rownames = treatednms,
                                  colnames = controlnms)

  return(ret)
}
