################################################################################
# exactMatch: create InfinitySparseMatrices from factors
################################################################################

setGeneric("exactMatch", 
  def = function(x, ...) {
    tmp <- standardGeneric("exactMatch")
    tmp@call <- match.call()
    return(tmp) 
})

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
#' match. Like \link{exactMatch}, the results of this function will
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
