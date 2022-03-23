##' This is a special function used only in identifying the strata variables
##' when defining an \code{exactMatch} during a call to \code{fullmatch},
##' \code{pairmatch}, or \code{match_on}. It should not be called externally.
##'
##' @title Identify Stratafication Variables
##' @param ... any number of variables of the same length.
##' @return the variables with appropriate labels
##' @export
##' @examples
##' data(nuclearplants)
##' fullmatch(pr ~ cost + strata(pt), data = nuclearplants)
strata <- function (...)
{
  allf <- list(...)
  if (!all(vapply(allf, length, 1) == length(allf[[1]]))) {
    stop("Variables in `strata` must be same length")
  }
  d <- do.call(paste, allf)
  strata_NA <-  apply(vapply(allf, is.na, logical(length(d))), 1, any)
  d[strata_NA] <- NA
  d
}
