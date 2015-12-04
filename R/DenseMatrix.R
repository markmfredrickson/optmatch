# Simple wrapper around the stock R matrix class to keep a few extra booking
# details.

setClassUnion("OptionalCall", c("call", "NULL"))
setClass("DenseMatrix",
  representation(call = "OptionalCall"),
  contains = "matrix")

#' @export
as.matrix.DenseMatrix <- function(x, ...) { x@.Data }
