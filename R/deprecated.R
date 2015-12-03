# Welcome to the function graveyard

#' Functions deprecated or removed from optmatch
#'
#' Over the course of time, several functions in optmatch have been removed in
#' favor of new interfaces and functions.
#' @name optmatch-defunct
#' @param ... All arguments ignored.
NULL

#' All functionality of the \code{pscore.dist} function has been moved into to
#' the \code{\link{mdist}} function. Additionally, this function will also act on
#' other objects, such as formulas.
#' The \code{\link{match_on}} function
#' also provides similar functionality, though with a different syntax.
#' @seealso \code{\link{mdist}}, \code{\link{match_on}}
#' @rdname optmatch-defunct
#' @export
pscore.dist <- function(...) {
  .Defunct(c("mdist", "match_on"), "optmatch")
}

#' All functionality of the \code{mahal.dist} function has been moved into to
#' the \code{\link{mdist}} function. Additionally, this function will also act on
#' other objects, such as \code{glm} objects. The \code{\link{match_on}} function
#' also provides similar functionality, though with a different syntax.
#' @rdname optmatch-defunct
#' @export
mahal.dist <- function(...) {
  .Defunct(c("mdist", "match_on"), "optmatch")
}
