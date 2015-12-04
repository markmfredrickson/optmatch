##' Functions to find the largest value of min.controls, or the
##' smallest value of max.controls, for which a full matching problem
##' is feasible.  These are determined by constraints embedded in the
##' matching problem's distance matrix.
##'
##' The function works by repeated application of full matching, so on
##' large problems it can be time-consuming.
##'
##' @note Essentially this is just a line search.  I've done several
##'   things to speed it up, but not everything that might be done.
##'   At present, not very thoroughly tested either: you might check
##'   the final results to make sure that \code{\link{fullmatch}}
##'   works with the values of \code{min.controls} (or
##'   \code{max.controls}) suggested by these functions, and that it
##'   ceases to work if you increase (decrease) those values.
##'   Comments appreciated.
##'
##' @title Set thinning and thickening caps for full matching
##' @param distance Either a matrix of non-negative, numeric
##'   discrepancies, or a list of such matrices. (See
##'   \code{\link{fullmatch}} for details.)
##' @param max.controls Optionally, set limits on the maximum number
##'   of controls per matched set.  (Only makes sense for
##'   \code{minControlsCap}.)
##' @param min.controls Optionally, set limits on the minimum number
##'   of controls per matched set.  (Only makes sense for
##'   \code{maxControlsCap}.)
##' @return For \code{minControlsCap},
##'   \code{strictest.feasible.min.controls} and
##'   \code{given.max.controls}. For \code{maxControlsCap},
##'   \code{given.min.controls} and
##'   \code{strictest.feasible.max.controls}.
##'
##'   \item{strictest.feasible.min.controls}{The largest values of the
##'   \code{\link{fullmatch}} argument \code{min.controls} that yield
##'   a full match;}
##'   \item{given.max.controls}{ The \code{max.controls} argument
##'   given to \code{minControlsCap} or, if none was given, a vector
##'   of \code{Inf}s. }
##'   \item{given.min.controls}{ The \code{min.controls} argument
##'   given to \code{maxControlsCap} or, if none was given, a vector
##'   of \code{0}s; }
##'   \item{strictest.feasible.max.controls}{The smallest values of
##'   the \code{\link{fullmatch}} argument \code{max.controls} that
##'   yield a full match.}
##'
##' @author Ben B. Hansen
##' @references Hansen, B.B. and S. Olsen Klopfer (2006),
##'   \sQuote{Optimal full matching and related designs via network
##'   flows}, \emph{Journal of Computational and Graphical Statistics}
##'   \bold{15}, 609--627.
##'
##' @seealso \code{\link{fullmatch}}
##' @keywords optimize
##' @export
##' @rdname minmaxctlcap
minControlsCap <- function(distance, max.controls=NULL)
{
  distance <- as.matrix(distance) # cast ISM to matrix, temporary
  if (!is.list(distance) & !is.matrix(distance))
    stop("argument \'distance\' must be a matrix or list")

  if (is.matrix(distance))
  {
    tmp <- maxControlsCap(t(distance),  min.controls=
                                          switch(1+is.null(max.controls),
                                                 ifelse(max.controls>=1, 1/ceiling(max.controls),
                                                        floor(1/max.controls) ), NULL))
  } else   {
    tmp <- maxControlsCap(lapply(distance, t), min.controls=
                                                 switch(1+is.null(max.controls),
               ifelse(max.controls>=1, 1/ceiling(max.controls),
                      floor(1/max.controls) ), NULL))
  }

  list(strictest.feasible.min.controls=
         (1/tmp$strictest.feasible.max.controls),
       given.max.controls=(1/tmp$given.min.controls))
}
