################################################################################
# Utility functions used throughout the package
################################################################################

### toZ: turn various objects into a logical vector indicating treatment
### accepts a variety of inputs and keeps names/rownames if present

setGeneric("toZ", function(x) {   
  
  if (any(is.na(x))) {
    stop("NAs not allowed in treatment indicator.")  
  }

  if (is.data.frame(x) | is.matrix(x)) {
    if (dim(x)[2] > 1) {
      stop("Treatment indicators must be vectors or single column
      matrices/data.frames")
    }

    nms <- rownames(x)
    x <- x[,1]
    names(x) <- nms
  }

  if (length(unique(x)) != 2) {
    stop(paste("Treatment indicator must have exactly 2 levels not",
      length(unique(x))))  
  }
  
  nms <- names(x)
  tmp <- standardGeneric("toZ")
  names(tmp) <- nms
  return(tmp)
})

# a noop
setMethod("toZ", "logical", identity)

# we already know it has two levels, so just call as.logical
setMethod("toZ", "numeric", function(x) as.logical(x))

setMethod("toZ", "character", function(x) toZ(as.factor(x)))

setMethod("toZ", "factor", function(x) {
  toZ(as.numeric(x) - 1)  
})



