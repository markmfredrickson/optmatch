#' @export
unmatched <- is.na

#' @export
matched <- function(x) !is.na(x)
