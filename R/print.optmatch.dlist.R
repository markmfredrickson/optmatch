print.optmatch.dlist <- function(x, ...)
{
  cl <- match.call()
  m <- match('digits',names(cl),0)
  if (m) DIGITS <- cl[[m]] else DIGITS <- 2
  print(lapply(x, function(X) signif(X, DIGITS)))
}
