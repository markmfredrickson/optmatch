##' @title strata
##' @param ... any number of variables of the same length.
##' @return the variables with appropriate labels
##' @export
strata <- function (...)
{
  allf <- list(...)
  d <- do.call(paste, allf)
  strata_NA <-  apply(sapply(allf, is.na), 1, any)
  d[strata_NA] <- NA
  d
}
