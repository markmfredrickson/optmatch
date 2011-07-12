print.optmatch <- function(x, quote = FALSE, grouped = FALSE, ...)
{
  if (length(x) <= 0) {
    cat("factor(0)\n")
  }
  else 
  {
    if (grouped) {
      tmp <- aggregate(names(x), by = list(x), FUN = as.character)
      colnames(tmp) <- c("Group", "Members")
      print(tmp, quote = quote, row.names = FALSE, ...)
    } else {
      rv <- as.character(x)
      names(rv) <- names(x)
      print(rv, quote = quote, ...)
    }
  }
  invisible(x)
}
