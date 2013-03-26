# Simple wrapper around the stock R matrix class to keep a few extra booking
# details.

setClassUnion("OptionalCall", c("call", "NULL"))
setClass("DenseMatrix", 
  representation(call = "OptionalCall"),
  contains = "matrix")

#' @S3method as.matrix DenseMatrix
as.matrix.DenseMatrix <- function(x, ...) { x@.Data }
