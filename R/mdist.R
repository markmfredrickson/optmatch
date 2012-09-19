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
                                       inv.scale.matrix = NULL, COV = cov, ...) {
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
  
  if (!is.null(inv.scale.matrix)) {
    # TODO: error check the inv.cov matrix to make sure it is safe
    # should match dimension of mf
  } else {
    # default inv.scale.matrix is the inverse covariance matrix
    # the extra as.matrix() is that if there is only one variable, it will be
    # a matrix not a vector
    mt <- COV(data[z, ,drop=FALSE]) * (sum(z) - 1) / (length(z) - 2)
    mc <- COV(data[!z, ,drop=FALSE]) * (sum(!z) - 1) / (length(!z) - 2)

    inv.scale.matrix <- solve(mt + mc) # don't need Z in the cov matrix

    # the old mahal.dist wrapped the solve in a try() and used this if there
    # was failure. Is this a common issue? I'm waiting for a failure case
    # before turning this code on (and with adjustments to the different
    # variable names, etc.
    # 
    # if (inherits(icv,"try-error"))
    # {
    #    dnx <- dimnames(cv)
    #    s <- svd(cv)
    #    nz <- (s$d > sqrt(.Machine$double.eps) * s$d[1])
    #    if (!any(nz))
    #      stop("covariance has rank zero")

    #    icv <- s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
    #    dimnames(icv) <- dnx[2:1]
    # }

  }

  f <- function(treated, control) {
    n <- dim(treated)[1]
    tmp <- numeric(n) 
    for (i in 1:n) {
      tmp[i] <- t(as.matrix(treated[i,] - control[i,])) %*% inv.scale.matrix %*% as.matrix(treated[i,] - control[i,])
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

  mdist(lp.adj, z = z, exclusions = exclusions, ...)
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


# TODO: make this real roxygen
# Returns the absolute difference for treated and control units computed using
# the vector of scores \code{x}.
#
# @param z Vector of treatment assignments for each unit in \code{x}. Either
#   \code{x} or \code{z} must have names.
setMethod("mdist", "numeric", function(x, z, exclusions = NULL, ...)
{
  if(is.null(z) || missing(z)) {
    stop("You must supply a treatment indicator, 'z', when using the numeric mdist method.")
  }

  if(length(z) != length(x)) {
    stop("The scores, 'x', and the treatment vector, 'z', must be the same length.")
  }

  z <- toZ(z)

  f <- function(t, c) { abs(t - c) }
  makedist(z, x, f, exclusions)
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

