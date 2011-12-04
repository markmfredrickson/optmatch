setGeneric("caliper", function(x, width = 1, exclude = c(), compare = `<=`, ...)
  standardGeneric("caliper"))

setMethod("caliper", "InfinitySparseMatrix",
function(x, width = 1, exclude = c(), compare = `<=`, ...) {

  excluded.rows <- which(x@rownames %in% exclude)
  excluded.cols <- which(x@colnames %in% exclude)

  y <- discardOthers(x, compare(x, width) | 
                     x@rows %in% excluded.rows |
                     x@cols %in% excluded.cols)

  y@.Data <- rep(0, length(y@.Data))

  return(y)
})

setMethod("caliper", "matrix",
function(x, width = 1, exclude = c(), compare = `<=`, ...) {
  caliper(as.InfinitySparseMatrix(x), width = width, exclude = exclude, compare = compare, ...)  
})
