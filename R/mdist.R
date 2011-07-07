################################################################################
# Mdist: distance matrix creation functions
################################################################################

setGeneric("mdist", def = function(x, exclusions = NULL, ...)  standardGeneric("mdist"))

setMethod("mdist", "function", function(x, exclusions = NULL, z = NULL, data = NULL, ...) {

  if (is.null(data) | is.null(z)) {
    stop("Data and treatment indicator arguments are required.")
  }

  theFun <- match.fun(x)

  makedist(z, data, theFun, exclusions)
})


# mdist method: formula
setMethod("mdist", "formula", function(x, exclusions = NULL, data = NULL, subset = NULL, 
                                       s.matrix = NULL, ...) {
  if (length(x) != 3) {
    stop("Formula must have a left hand side.")  
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("x", "data", "subset"), # maybe later add "na.action"
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  names(mf)[names(mf) == "x"] <- "formula"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  if (dim(mf)[2] < 2) {
    stop("Formula must have a right hand side with at least one variable.")   
  }

  data <- subset(model.matrix(x, mf), T, -1) # drop the intercept

  z <- toZ(mf[,1])
  names(z) <- rownames(mf)
  
  if (!is.null(s.matrix)) {
    # TODO: error check the inv.cov matrix to make sure it is safe
    # should match dimension of mf
  } else {
    # default s.matrix is the inverse covariance matrix
    # the extra as.matrix() is that if there is only one variable, it will be
    # a matrix not a vector
    s.matrix <- solve(cov(as.matrix(data))) # don't need Z in the cov matrix
  }

  f <- function(treated, control) {
    n <- dim(treated)[1]
    tmp <- numeric(n) 
    for (i in 1:n) {
      tmp[i] <- sqrt(t(as.matrix(treated[i,] - control[i,])) %*% s.matrix %*% as.matrix(treated[i,] - control[i,]))
    }
    return(tmp)
  }



  makedist(z, data, f, exclusions)

})

# mdist method: glm
setMethod("mdist", "glm", function(x, exclusions = NULL, standardization.scale = mad, ...)
{
  stopifnot(all(c('y', 'linear.predictors','data') %in% names(x)))
  z <- x$y > 0
  pooled.sd <- if (is.null(standardization.scale)) {
    1 
  } else {
    szn.scale(x$linear.predictors, z ,standardization.scale)
  }

  lp.adj <- x$linear.predictors/pooled.sd

  f <- function(t, c) { abs(t - c) }
  
  makedist(z, lp.adj, f, exclusions)
})

szn.scale <- function(x, Tx, standardizer = mad, ...) {
  sqrt(((sum(!Tx) - 1) * standardizer(x[!Tx])^2 + 
        (sum(!!Tx) - 1) * standardizer(x[!!Tx])^2) / (length(x) - 2))
}

# mdist method: bigglm
setMethod("mdist", "bigglm", function(x, exclusions = NULL, data = NULL, standardization.scale = mad, ...)
{
  if (is.null(data)) {
    stop("data argument is required for computing mdists from bigglms")
  }

  if (!is.data.frame(data)) {
    stop("mdist doesn't understand data arguments that aren't data frames.")
  }

  theps <- predict(x, data, type = 'link', se.fit = FALSE)

  if (length(theps) != dim(data)[1]) {
    stop("predict.bigglm() returns a vector of the wrong length;
are there missing values in data?")
  }

  # this makes heavy use of the bigglm terms object, the original formula
  # if this implementation detail changes in later versions of bigglm,
  # look here for the problem.

  Data <-  model.frame(x$terms, data = data)
  z <- Data[, 1]
  pooled.sd <- if (is.null(standardization.scale)) {
    1
  } else { 
    szn.scale(theps, z, standardizer=standardization.scale,...)
  }
  
  psdiffs <- function(treatments, controls) {
    abs(treatments - controls) / pooled.sd
  }
  
  makedist(z, theps, psdiffs, exclusions)
      
})


### mdist method: numeric.
### (mdist can't work with numeric vectors at present,
### but it can return an informative error message).

setMethod("mdist", "numeric", function(x, exclusions = NULL, ...)
{

  stop("No mdist method for numerics.
  Consider using mdist(z ~ ps | strata, data = your.data)
  where ps is your numeric vector, z is your treatment assignment,
  and strata (optional) indicates a stratification variable, all
  columns in your.data")

})

# mdist methods for DistanceSpecifications
# apparently the class union is less important than the true
# type, so the numeric method above gets in the way
setMethod("mdist", "InfinitySparseMatrix", function(x, exclusions = NULL, ...) {
  return(x)
}) # just return the argument

setMethod("mdist", "matrix", function(x, exclusions = NULL, ...) {
  return(x)
}) # just return the argument

