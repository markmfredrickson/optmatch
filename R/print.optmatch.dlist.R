#' @keywords internal
print.optmatch.dlist <- function(x, ...)
{
  if (is(x, "optmatch.dlist")) {
    warning("The use of 'optmatch.dlist' objects created by 'mdist()' is deprecated.\nPlease use 'match_on()' instead.")
  }

  cl <- match.call()
  m <- match('digits',names(cl),0)
  if (m) DIGITS <- cl[[m]] else DIGITS <- 2
  print(lapply(x, function(X) signif(X, DIGITS)))
}
