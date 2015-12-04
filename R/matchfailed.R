#' @export
#' @rdname matched
matchfailed <- function(x) {
  failed <- !subproblemSuccess(x)
  grps <- attr(x, "subproblem")
  levels(grps) <- failed
  return(as.logical(grps))
}

#' (Internal) Report successful subproblems.
#'
#' \code{\link{fullmatch}} can break up a large matching problem into smaller
#' subproblems (for example, using strata defined by \code{\link{exactMatch}}).
#' This function lists the subproblems in a match and list whether at least on
#' treated unit was matched in subproblem. Subproblems that have no matched
#' treated units are said to have "failed."
#'
#' @param x The result of \code{\link{fullmatch}} or \code{\link{pairmatch}}.
#' @return A named logical vector indicating either success or failure for each subproblem.
subproblemSuccess <- function(x) {
  grps <- attr(x, "subproblem")
  failed <- sapply(split(x, grps), function(x) { all(is.na(x)) })
  return(!failed)
}
