##' @title strata
##' @param ... any number of variables of the same length.
##' @return the variables with appropriate labels
##' @export
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
