################################################################################
# Mdist: distance matrix creation functions
################################################################################

#' Create treated to control distances for matching problems
#' 
#' A generic function, with several supplied methods, for creating
#' matrices of distances between observations to be used in the match process.
#' Using these matrices, 
#' \code{pairmatch()} or \code{fullmatch()} can determine optimal matches.
#'
#' The \code{function} method takes as its \code{x} argument a function of two
#' arguments.  This function should in turn expect two arguments, a
#' \code{data.frame} carrying information about treatment group members and a
#' \code{data.frame} carrying information about control units.  The
#' \code{function} method also expects a \code{z} argument, a vector indicating
#' whether each unit is in the treated or control groups; on the basis of this
#' vector the \code{data} argument is what is chopped up and fed to the function
#' passed to \code{mdist} as argument \code{x}. While this might sound complicated
#' at first, it's quite flexible and, once you're used to it, easy to use. For
#' example, the simple 1-dimensional distance of treated and control units can be
#' implemented as: \code{d <- 1:20 ; names(d) <- letters[1:20] ; mdist(`-`, z =
#' rep(c(T,F),10), data = d)}.  See also additional examples below.  
#'
#' The formula method produces, by default, a Mahalanobis distance specification
#' based on the formula \code{treatment ~ var1 + var2 + ... }, where
#' \code{treatment} iexclusions = NULLs the treatment indicator. A Mahalanobis
#' distance scales the squared Euclidean distance by the inverse of the
#' covariance matrix. Other scale matrices can be supplied, for example, the
#' identity matrix will result in squared Euclidean distance. If you wish to use
#' an alternative function to compute the covariance, pass it in the \code{COV}
#' argument.
#'
#' An \code{glm} method takes an argument of class \code{glm} as
#' first argument.  It assumes that this object is a fitted propensity
#' model, extracting distances on the linear propensity score (logits of
#' the estimated conditional probabilities) and, by default, rescaling the distances
#' by the reciprocal of the pooled s.d. of treatment- and control-group propensity scores.
#' (The scaling uses \code{mad}, for resistance to outliers, by default; this can be
#' changed to the actual s.d., or rescaling can be skipped entirely, by
#' setting argument \code{standardization.scale} to \code{sd} or
#' \code{NULL}, respectively.)  The \code{bigglm}
#' method works analogously with \code{bigglm} objects, created by
#' the \code{bigglm} function from package \sQuote{biglm}, which can
#' handle bigger data sets than the ordinary glm function can. 
#'
#' The \code{numeric} method is simply a placeholder that returns an error
#' message and some helpful suggestions on how to create distance matrices, as
#' there is no obvious way to create distance matrix from just a numeric vector.
#'
#' The \code{matrix} and \code{InfinitySparseMatrix} just return their
#' arguments as these objects are already valid distance specifications.
#'
#' 
#' @param x An object defining how to create the distances
#' @param exclusions A \code{\link{DistanceSpecification}} such as the result
#' of \code{\link{exactMatch}} or \code{\link{caliper}}. Finite entries indicate
#' which distances to create. Including this argument can significantly speed up
#' computation for sparse matching problems.
#' @return Object of the class union \code{DistanceSpecification} (either a
#' \code{matrix} or a \code{InfinitySparseMatrix} or derived class), which is
#' suitable to be given as \code{distance} argument to \code{\link{fullmatch}}
#' or \code{\link{pairmatch}}. 
#' @seealso \code{\link{fullmatch}}, \code{\link{pairmatch}}, \code{\link{exactMatch}}, \code{\link{caliper}}
#' @references
#' P.~R. Rosenbaum and D.~B. Rubin (1985),
#' \sQuote{Constructing a control group using multivariate matched sampling
#'   methods that incorporate the propensity score},
#'  \emph{The American Statistician}, \bold{39} 33--38.
#' @export
#' @example inst/examples/match_on.R
#' @docType methods
#' @rdname match_on-methods
setGeneric("match_on", def = function(x, exclusions = NULL, ...)  standardGeneric("match_on"))

#' @param z A factor, logical, or binary vector indicating treatment (the higher level) and control (the lower level) for each unit in the study.
#' @param data A \code{data.frame} or \code{matrix} containing variables used by the method to construct the distance matrix.
#' @rdname match_on-methods
#' @aliases match_on,function,ANY-method
setMethod("match_on", "function", function(x, exclusions = NULL, z = NULL, data = NULL, ...) {

  if (is.null(data) | is.null(z)) {
    stop("Data and treatment indicator arguments are required.")
  }

  theFun <- match.fun(x)

  makedist(z, data, theFun, exclusions)
})


#' @param subset A subset of the data to use in creating the distance specification.
#' @param inv.scale.matrix A matrix that scales the distance computation. The
#' default Mahalanobis distance scales the squared Euclidean distance by
#' the inverse of the covariance matrix. Other scale matrices can be supplied, for example,
#' the identity matrix will result in squared Euclidean distance.
#' @param COV A covariance computing function. The default is \code{\link{cov}}
#' @rdname match_on-methods
setMethod("match_on", "formula", function(x, exclusions = NULL, data = NULL, subset = NULL, 
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

#' @param standardization.scale Standardizes the data based on the median absolute deviation (by default).
#' @rdname match_on-methods
setMethod("match_on", "glm", function(x, exclusions = NULL, standardization.scale = mad, ...)
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

#' @rdname match_on-methods
setMethod("match_on", "bigglm", function(x, exclusions = NULL, data = NULL, standardization.scale = mad, ...)
{
  if (is.null(data)) {
    stop("data argument is required for computing match_ons from bigglms")
  }

  if (!is.data.frame(data)) {
    stop("match_on doesn't understand data arguments that aren't data frames.")
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


### match_on method: numeric.
### (match_on can't work with numeric vectors at present,
### but it can return an informative error message).

#' @rdname match_on-methods
setMethod("match_on", "numeric", function(x, exclusions = NULL, ...)
{

  stop("No match_on method for numerics.
  Consider using match_on(z ~ ps | strata, data = your.data)
  where ps is your numeric vector, z is your treatment assignment,
  and strata (optional) indicates a stratification variable, all
  columns in your.data")

})

# match_on methods for DistanceSpecifications
# apparently the class union is less important than the true
# type, so the numeric method above gets in the way
#' @rdname match_on-methods
setMethod("match_on", "InfinitySparseMatrix", function(x, exclusions = NULL, ...) {
  return(x)
}) # just return the argument

#' @rdname match_on-methods
setMethod("match_on", "matrix", function(x, exclusions = NULL, ...) {
  return(x)
}) # just return the argument

