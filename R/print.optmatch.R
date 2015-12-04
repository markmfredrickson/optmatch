#' Printing \code{optmatch} objects.
#'
#' @param x The \code{optmatch} object, as returned by
#'   \code{\link{fullmatch}} or \code{\link{pairmatch}}.
#' @param grouped A logical indicating if the object should printed in
#'   the style of a named \code{factor} object (\code{grouped = TRUE})
#'   or as a table of group names and members.
#' @param quote A boolean indicating if the matched group names should
#'   be quoted or not (default is not to quote).
#' @param ... Arguments passed to \code{\link{print.default}}.
#' @seealso \code{\link{fullmatch}}, \code{\link{pairmatch}},
#'   \code{\link{print}}, \code{\link{summary.optmatch}}
#' @example inst/examples/print.optmatch.R
#' @method print optmatch
#' @export
print.optmatch <- function(x, quote = FALSE, grouped = FALSE, ...)
{
  if (length(x) <= 0) {
    cat("factor(0)\n")
  }
  else
  {
    if (grouped) {
      tmp <- aggregate(names(x), by = list(x), FUN = function(x) { paste(x,
        collapse = ", ")})
      colnames(tmp) <- c("Group", "Members")
      print(tmp, quote = FALSE, row.names = FALSE, ...)
    } else {
      rv <- as.character(x)
      names(rv) <- names(x)
      print(rv, quote = FALSE, ...)
    }
  }
  invisible(x)
}
