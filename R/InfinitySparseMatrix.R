################################################################################
### InfinitySparseMatrix: A sparsem matrix class where infinity is the default
###                       instead of zero.
################################################################################

# instantiate with:
# a <- new("InfinitySparseMatrix", c(1,2,3), cols = c(1,2, 2), rows = c(1,2,4))
# To get a matrix like:
# 1   2
# Inf 3
setClass("InfinitySparseMatrix", 
  representation(cols = "numeric", rows = "numeric"),
  contains = "numeric")

### Basic Matrix-like Operations ###
dim.InfinitySparseMatrix <- function(x) { c(length(x@rows) - 1, max(x@cols)) }
