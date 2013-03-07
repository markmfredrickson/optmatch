################################################################################
# Mdist: distance matrix creation functions
################################################################################

#' Create treated to control distances for matching problems
#' 
#' A function with which to produce matching distances, for instance Mahalanobis
#' distances, propensity score discrepancies or calipers, or combinations thereof, for 
#' \code{\link{pairmatch}} or \code{\link{fullmatch}} to subsequently \dQuote{match on}.
#' Conceptually, the result of a call \code{match_on} is a treatment-by-control matrix of distances. 
#' Because these matrices can grow quite large, in practice \code{match_on} produces either an 
#' ordinary dense matrix or a special sparse matrix structure (that can make use of caliper and exact matching
#' constraints to reduce storage requirements).  Methods are supplied for these sparse structures, 
#'\code{InfinitySparseMatrix}es, so that they can be manipulated and modified in much the same way as dense matrices.
#'
#' \code{match_on} is generic. There are several supplied methods, all providing the same basic output: a matrix (or
#' similar) object with treated units on the rows and control units on the
#' columns. Each cell [i,j] then indicates the distance from a treated unit i to
#' control unit j. Entries that are \code{Inf} are said to be unmatchable. Such
#' units are guaranteed to never be in a matched set. For problems with many
#' \code{Inf} entries, so called sparse matching problems, \code{match_on} uses a
#' special data type that is more space efficient than a standard R \code{matrix}.
#' When problems are not sparse (i.e. dense), \code{match_on} uses the standard
#' \code{matrix} type.
#'
#' \code{match_on} methods differ on the types of arguments they take, making
#' the function a one-stop location of many different ways of specifying matches:
#' using functions, formulas, models, and even simple scores. Many of the methods
#' require additional arguments, detailed below. All methods take a \code{within}
#' argument, a distance specification made using \code{\link{exactMatch}} or
#' \code{\link{caliper}} (or some additive combination of these or other distance
#' creating functions). All \code{match_on} methods will
#' use the finite entries in the \code{within} argument as a guide for producing
#' the new distance. Any entry that is \code{Inf} in \code{within} will be
#' \code{Inf} in the distance matrix returned by \code{match_on}. This argument
#' can reduce the processing time needed to compute sparse distance matrices.
#'
#' The \code{match_on} function is similar to the older, but still supplied,
#' \code{\link{mdist}} function. Future development will concentrate on
#' \code{match_on}, but \code{mdist} is still supplied for users familiar with the
#' interface. For the most part, the two functions can be used interchangeably by
#' users.
#' 
#' @param x An object defining how to create the distances
#' @param within A valid distance specification, such as the result
#' of \code{\link{exactMatch}} or \code{\link{caliper}}. Finite entries indicate
#' which distances to create. Including this argument can significantly speed up
#' computation for sparse matching problems.
#' @return A distance specification (a matrix or similar object) which is
#' suitable to be given as the \code{distance} argument to \code{\link{fullmatch}}
#' or \code{\link{pairmatch}}. 
#' @param ... Other arguments for methods.
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
#' @aliases InfinitySparseMatrix-class
setGeneric("match_on", def = function(x, within = NULL, ...) {

  tmp <- standardGeneric("match_on")
  tmp@call <- match.call()
  return(tmp)

})

#' @details The \code{function} method takes as its \code{x} argument a
#' function of three arguments: \code{index}, \code{data}, and \code{z}. The \code{data} and \code{z}
#' arguments will be the same as those passed directly to
#' \code{match_on}. The \code{index} argument is a matrix of two columns,
#' representing the pairs of treated and control units that are valid comparisons
#' (given any \code{within} arguments). The first column is the row name or id of
#' the treated unit in the \code{data} object. The second column is the id for the
#' control unit, again in the \code{data} object. For each of these pairs, the
#' function should return the distance between the treated unit and control unit.
#' This may sound complicated, but is simple to use. For example, a function that
#' returned the absolute difference between to units using a vector of data would
#' be 
#' \code{f <- function(index, data, z) { abs(apply(index, 1, function(pair) { data[pair[1]] - data[pair[2]] })) }}.
#' (Note: This simple case is precisely handled by the \code{numeric} method.)
#' 
#' @param z A factor, logical, or binary vector indicating treatment (the higher level) and control (the lower level) for each unit in the study.
#' @param data A \code{data.frame} or \code{matrix} containing variables used by the method to construct the distance matrix.
#' @usage \S4method{match_on}{function}(x, within = NULL, z = NULL, data = NULL, ...)
#' @rdname match_on-methods
#' @aliases match_on,function-method
setMethod("match_on", "function", function(x, within = NULL, z = NULL, data = NULL, ...) {

  if (is.null(data) | is.null(z)) {
    stop("Data and treatment indicator arguments are required.")
  }

  theFun <- match.fun(x)

  makedist(z, data, theFun, within)
})


#' @details The formula method produces, by default, a Mahalanobis distance specification
#' based on the formula \code{Z ~ X1 + X2 + ... }, where
#' \code{Z}  the treatment indicator. A Mahalanobis
#' distance scales the Euclidean distance by the inverse of the
#' covariance matrix. Other options can be selected by the \code{method} argument. 
#' @param subset A subset of the data to use in creating the distance specification.
#' @param method A string indicating which method to use in computing the distances from the data. 
#' The current possibilities are \code{"mahalanobis", "euclidean"}. 
#' @usage \S4method{match_on}{formula}(x, within = NULL, data = NULL, subset = NULL, method = "mahalanobis", ...)
#' @rdname match_on-methods
#' @aliases match_on,formula-method
setMethod("match_on", "formula", function(x, within = NULL, data = NULL, subset = NULL, 
                                       method = "mahalanobis", ...) {
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
 
  f <- match.fun(paste("compute_", method, sep = ""))
  makedist(z, data, f, within)

})

compute_mahalanobis <- function(index, data, z) {
	
	mt <- cov(data[z, ,drop=FALSE]) * (sum(z) - 1) / (length(z) - 2)
	mc <- cov(data[!z, ,drop=FALSE]) * (sum(!z) - 1) / (length(!z) - 2)
	cv <- mt + mc
	
	inv.scale.matrix <- try(solve(cv))
	
	if (inherits(inv.scale.matrix,"try-error"))
	{
		dnx <- dimnames(cv)
		s <- svd(cv)
		nz <- (s$d > sqrt(.Machine$double.eps) * s$d[1])
		if (!any(nz))
		stop("covariance has rank zero")
		
		inv.scale.matrix <- s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
		dimnames(inv.scale.matrix) <- dnx[2:1]
	}
	
	nv <- nrow(index)
	
	result <- .C('mahalanobisHelper',
    as.integer(nv),
    as.integer(ncol(data)),
    t(data[index[, 1], ]),
    t(data[index[, 2], ]),
    t(inv.scale.matrix),
    result=numeric(nv), PACKAGE='optmatch')$result
	
	return(result)
}

# short alias if we need it
compute_mahal <- compute_mahalanobis

compute_euclidean <- function(index, data, z) {

  sqrt(apply(index, 1, function(pair) {

    pair.diff <- as.matrix(data[pair[1],] - data[pair[2],])
    t(pair.diff) %*% pair.diff
  
  }))
}

# short alias
compute_euclid <- compute_euclidean

#' @details The \code{glm} method accepts a fitted propensity
#' model, extracts distances on the linear propensity score (logits of
#' the estimated conditional probabilities), and rescales the distances
#' by the reciprocal of the pooled s.d. of treatment- and control-group propensity scores.
#' (The scaling uses \code{mad}, for resistance to outliers, by default; this can be
#' changed to the actual s.d., or rescaling can be skipped entirely, by
#' setting argument \code{standardization.scale} to \code{sd} or \code{NULL}, respectively.) 
#' The resulting distance matrix is the absolute difference between treated and
#' control units on the rescaled propensity scores. This method relies on the
#' \code{numeric} method, so you may pass a \code{caliper} argument.
#'
#' @param standardization.scale Standardizes the data based on the median absolute deviation (by default).
#' @usage \S4method{match_on}{glm}(x, within = NULL, standardization.scale = mad, ...)
#' @rdname match_on-methods
#' @aliases match_on,glm-method
setMethod("match_on", "glm", function(x, within = NULL, standardization.scale = mad, ...)
{
  stopifnot(all(c('y', 'linear.predictors','data') %in% names(x)))
  z <- x$y > 0
  pooled.sd <- if (is.null(standardization.scale)) {
    1 
  } else {
    szn.scale(x$linear.predictors, z ,standardization.scale)
  }

  lp.adj <- x$linear.predictors/pooled.sd

  match_on(lp.adj, within = within, z = z, ...)
})

szn.scale <- function(x, Tx, standardizer = mad, ...) {
  sqrt(((sum(!Tx) - 1) * standardizer(x[!Tx])^2 + 
        (sum(!!Tx) - 1) * standardizer(x[!!Tx])^2) / (length(x) - 2))
}

#' @details The \code{bigglm}
#' method works analogously to the \code{glm} method, but with \code{bigglm} objects, created by
#' the \code{bigglm} function from package \sQuote{biglm}, which can
#' handle bigger data sets than the ordinary glm function can.
#'
#' @usage \S4method{match_on}{bigglm}(x, within = NULL, data = NULL, standardization.scale = mad, ...)
#' @rdname match_on-methods
#' @aliases match_on,bigglm-method
setMethod("match_on", "bigglm", function(x, within = NULL, data = NULL, standardization.scale = mad, ...)
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

  theps <- as.vector(theps / pooled.sd)
  names(theps) <- rownames(data)

  match_on(theps, within = within, z = z, ... )    
})

#' @details The \code{numeric} method returns the absolute difference for treated and control units computed using
#' the vector of scores \code{x}. Either \code{x} or \code{z} must have names.
#' @param caliper The width of a caliper to fit on the difference of scores.
#'   This can improve efficiency versus first creating all the differences and
#'   then filtering out those entries that are larger than the caliper.
#' @usage \S4method{match_on}{numeric}(x, within = NULL, z, caliper = NULL, ...)
#' @rdname match_on-methods
#' @aliases match_on,numeric-method
setMethod("match_on", "numeric", function(x, within = NULL, z, caliper = NULL, ...) {

  if(missing(z) || is.null(z)) {
    stop("You must supply a treatment indicator, 'z', when using the numeric match_on method.")
  }

  if(length(z) != length(x)) {
    stop("The scores, 'x', and the treatment vector, 'z', must be the same length.")
  }

  z <- toZ(z)

  if(!is.null(caliper)) {
    allowed <- scoreCaliper(x, z, caliper)
    
    if (!is.null(within)) {
      within <- within + allowed
    } else {
      within <- allowed
    }
  }
  
  f <- function(index, data, z) { abs(apply(index, 1, function(pair) { data[pair[1]] - data[pair[2]] })) }

  makedist(z, x, f, within)
})

#' (Internal) Helper function to create an InfinitySparseMatrix from a set of scores, a treatment indicator, and a caliper width.
#'
#' @param x The scores, a vector indicating the 1-D location of each unit.
#' @param z The treatment assignment vector (same length as \code{x})
#' @param caliper The width of the caliper with respect to the scores \code{x}.
#' @return An \code{InfinitySparseMatrix} object, suitable to be passed to \code{\link{match_on}} as an \code{within} argument.
scoreCaliper <- function(x, z, caliper) {
  z <- toZ(z)

  treated <- x[z]
  k <- length(treated)
  control <- x[!z]

  # the following uses findInterval, which requires a sorted vector
  # there may be a speed increase in pulling out the guts of that function and calling them directly
  control <- sort(control)
  
  treatedids <- c()
  controlids <- c()
  
  # NB: for reasons unknown, you must add the double.eps in the function
  # call, saving it in a variable (e.g. width.eps <- width +
  # .Machine$double.eps) will not work.
  # The use of double.eps is to get a <= treated <= b intervals 

  stops <- findInterval(treated + caliper + .Machine$double.eps, control)
  starts <- length(control) - findInterval(-(treated - caliper -
                                             .Machine$double.eps), rev(-control))
  
  for (i in 1:k) {
    if (starts[i] < length(control) && stops[i] > 0 && starts[i] < stops[i]) {
      tmp <- seq(starts[i] + 1, stops[i])
      controlids <- c(controlids, tmp)
      treatedids <- c(treatedids, rep(i, length(tmp)))
    }
  }

  makeInfinitySparseMatrix(rep(0, length(treatedids)), controlids, treatedids, names(control), names(treated))
}

#' @details The \code{matrix} and \code{InfinitySparseMatrix} just return their
#' arguments as these objects are already valid distance specifications.
#' 
#' @rdname match_on-methods
#' @aliases match_on,InfinitySparseMatrix-method
setMethod("match_on", "InfinitySparseMatrix", function(x, within = NULL, ...) {
  return(x)
}) # just return the argument

#' @rdname match_on-methods
#' @aliases match_on,matrix-method
setMethod("match_on", "matrix", function(x, within = NULL, ...) {
  return(x)
}) # just return the argument

