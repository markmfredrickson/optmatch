################################################################################
# match_on: distance matrix creation functions
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
#' @param x An object defining how to create the distances. All methods require
#' some form of names (e.g. \code{names} for vectors or \code{rownames} for
#' matrix like objects)
#' @param within A valid distance specification, such as the result
#' of \code{\link{exactMatch}} or \code{\link{caliper}}. Finite entries indicate
#' which distances to create. Including this argument can significantly speed up
#' computation for sparse matching problems.
#' @param caliper The width of a caliper to use to exclude treated-control
#' pairs with values greater than the width. For some methods, there may be a
#' speed advantage to passing a width rather than using the
#' \code{\link{caliper}} function on an existing distance specification.
#' @param data An optional data frame.
#' @param ... Other arguments for methods.
#' @return A distance specification (a matrix or similar object) which is
#' suitable to be given as the \code{distance} argument to \code{\link{fullmatch}}
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
#' @aliases InfinitySparseMatrix-class
match_on <- function(x, within = NULL, caliper = NULL, data=NULL, ...) {
  # if x does not exist then print helpful error msg
  x_str <- deparse(substitute(x))
  data_str <- deparse(substitute(data))
  tryCatch(x, error = function(e) {
    stop(missing_x_msg(x_str, data_str, ...))})

  cl <- match.call()
  UseMethod("match_on")
}

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
#' returned the absolute difference between two units using a vector of data would
#' be
#' \code{f <- function(index, data, z) { abs(apply(index, 1, function(pair) { data[pair[1]] - data[pair[2]] })) }}.
#' (Note: This simple case is precisely handled by the \code{numeric} method.)
#'
#' @param z A factor, logical, or binary vector indicating treatment (the higher level) and control (the lower level) for each unit in the study.
#' @method match_on function
#' @rdname match_on-methods
match_on.function <- function(x, within = NULL, caliper = NULL, data = NULL, z = NULL, ...) {

  if (is.null(data) | is.null(z)) {
    stop("Data and treatment indicator arguments are required.")
  }
  if (!exists("cl")) cl <- match.call()

  theFun <- match.fun(x)

  tmp <- makedist(z, data, theFun, within)

  if (is.null(caliper)) {
    tmp@call <- cl
    return(tmp)
  }

  tmp <- tmp + optmatch::caliper(tmp, width = caliper)
  tmp@call <- cl
  return(tmp)
}


#' @details The formula method produces, by default, a Mahalanobis distance specification
#' based on the formula \code{Z ~ X1 + X2 + ... }, where
#' \code{Z} is the treatment indicator. The Mahalanobis distance is calculated as the
#' square root of d'Cd, where d is the vector of X-differences on a pair of observations and C
#' is an inverse (generalized inverse) of the pooled covariance of Xes. (The pooling is of the
#' covariance of X within the subset defined by \code{Z==0} and within the complement of that
#' subset. This is similar to a Euclidean distance calculated after reexpressing the Xes in
#' standard units, such that the reexpressed variables all have pooled SDs of 1; except that
#' it addresses redundancies among the variables by scaling down variables contributions in
#' proportion to their correlations with other included variables.)
#'
#' Euclidean distance is also available, via
#' \code{method="euclidean"}, and ranked, Mahalanobis distance, via
#' \code{method="rank_mahalanobis"}. Or, implement your own;
#' for hints as to how, refer to\cr
#' \url{https://github.com/markmfredrickson/optmatch/wiki/How-to-write-your-own-compute-method}
#'
#' @param subset A subset of the data to use in creating the distance specification.
#' @param method A string indicating which method to use in computing the distances from the data.
#' The current possibilities are \code{"mahalanobis", "euclidean", "rank_mahalanobis"}, or pass a user created distance function.
#' @method match_on formula
#' @rdname match_on-methods
match_on.formula <- function(x, within = NULL, caliper = NULL, data = NULL, subset = NULL, method = "mahalanobis", ...) {
  if (length(x) != 3) {
    stop("Formula must have a left hand side.")
  }
  if (!exists("cl")) cl <- match.call()

  mf <- match.call(expand.dots = FALSE)

  m <- match(c("x", "data", "subset"), # maybe later add "na.action"
             names(mf), 0L)
  mf <- mf[c(1L, m)]

  rm(m)

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


  if(is.character(method)){
    methodname <- method
  } else {
    methodname <- as.character(class(method))
  }

  which.method <- pmatch(methodname, c("mahalanobis", "euclidean", "rank_mahalanobis", "function"), 3)
  tmp <- switch(which.method,
		makedist(z, data, compute_mahalanobis, within),
		makedist(z, data, compute_euclidean, within),
    makedist(z, data, compute_rank_mahalanobis, within),
		makedist(z, data, match.fun(method), within)
		)
  rm(mf)

  if (is.null(caliper)) {
    tmp@call <- cl
    return(tmp)
  }

  tmp <- tmp + optmatch::caliper(tmp, width = caliper)
  tmp@call <- cl
  return(tmp)
}

# compute_mahalanobis computes mahalanobis distances between treatment and
# control pairs
#
# Arguments:
#   index: a 2 col array of rownames from the 'data' argument.
#     Col 1: treatment rownames
#     Col 2: control rownames
#   data: a matrix containing rows of treatment and control data. The
#     rownames are used in index to indicate which treatment and control pairs
#     get measured
#   z: a logical vector of length nrows(data); TRUE indicates treatment
#
# If called from the makedist function, index is most likely a cross product of
# treatment and control rownames.
#
# Value: a vector of distances a distance for each pair indicated in index
#
# This is the default method for calculating distances in the match_on methods.
# It calls a registered C routine mahalanobisHelper found in distances.c
# after computing the inverse of a covariate matrix

# compute_mahalanobis computes mahalanobis distances between treatment and
# control pairs
#
# Arguments:
#   index: a 2 col array of rownames from the 'data' argument.
#     Col 1: treatment rownames
#     Col 2: control rownames
#   data: a matrix containing rows of treatment and control data. The
#     rownames are used in index to indicate which treatment and control pairs
#     get measured
#   z: a logical vector of length nrows(data); TRUE indicates treatment
#
# If called from the makedist function, index is most likely a cross product of
# treatment and control rownames.
#
# Value: a vector of distances a distance for each pair indicated in index
#
# This is the default method for calculating distances in the match_on methods.
# It calls a registered C routine mahalanobisHelper found in distances.c
# after computing the inverse of a covariate matrix

compute_mahalanobis <- function(index, data, z) {
  if (!all(is.finite(data))) stop("Infinite or NA values detected in data for Mahalanobis computations.")

    mt <- cov(data[z, ,drop=FALSE]) * (sum(z) - 1) / (length(z) - 2)
    mc <- cov(data[!z, ,drop=FALSE]) * (sum(!z) - 1) / (length(!z) - 2)
    cv <- mt + mc
    rm(mt, mc)

    inv.scale.matrix <- try(solve(cv), silent = TRUE)

    if (inherits(inv.scale.matrix,"try-error")) {
      dnx <- dimnames(cv)
      s <- svd(cv)
      nz <- (s$d > sqrt(.Machine$double.eps) * s$d[1])
      if (!any(nz)) stop("covariance has rank zero")

      inv.scale.matrix <- s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
      dimnames(inv.scale.matrix) <- dnx[2:1]
      rm(dnx, s, nz)
    }

    rm(cv)

    return(.Call(mahalanobisHelper, data, index, inv.scale.matrix))
}


compute_euclidean <- function(index, data, z) {

  if (!all(is.finite(data))) stop("Infinite or NA values detected in data for distance computations.")

  return(.Call(mahalanobisHelper, data, index, diag(ncol(data))))
}


compute_rank_mahalanobis <- function(index, data, z) {
    if (!all(is.finite(data))) {
        stop("Infinite or NA values detected in data for Mahalanobis computations.")
    }

    return(
        .Call('r_smahal', index, data, z, PACKAGE='optmatch')
    )
}

#' @details The \code{glm} method assumes its first argument to be a fitted propensity
#' model. From this it extracts distances on the \emph{linear} propensity score: fitted values
#' of the linear predictor, the link function applied to the estimated conditional probabilities,
#' as opposed to the estimated conditional probabilities themselves (Rosenbaum \& Rubin, 1985).
#' For example, a logistic model (\code{glm} with \code{family=binomial()}) has the logit function
#' as its link, so from such models \code{match_on} computes distances in terms of logits of
#' the estimated conditional probabilities, i.e. the estimated log odds.
#'
#' Optionally these distances are also rescaled. The default is to rescale, by the reciprocal of
#' an outlier-resistant variant of the pooled s.d. of propensity scores.
#' (Outlier resistance is obtained by the application of \code{mad}, as opposed to \code{sd},
#' to linear propensity scores in the treatment; this can be
#' changed to the actual s.d., or rescaling can be skipped entirely, by
#' setting argument \code{standardization.scale} to \code{sd} or \code{NULL}, respectively.)
#' The overall result records absolute differences between treated and
#' control units on linear, possibly rescaled, propensity scores.
#'
#' In addition, one can impose a caliper in terms of these distances by providing a scalar as a
#' \code{caliper} argument, forbidding matches between treatment and control units differing in the
#' calculated propensity score by more than the specified caliper.  For example, Rosenbaum and Rubin's (1985)
#' caliper of one-fifth of a pooled propensity score s.d. would be imposed by specifying \code{caliper=.2},
#' in tandem either with the default rescaling or, to follow their example even more closely, with the
#' additional specification \code{standardization.scale=sd}. Propensity calipers are beneficial
#' computationally as well as statistically, for reasons indicated in the below discussion of
#' the \code{numeric} method.
#'
#' @param standardization.scale Standardizes the data based on the median absolute deviation (by default).
#' @method match_on glm
#' @rdname match_on-methods
match_on.glm <- function(x, within = NULL, caliper = NULL, data = NULL, standardization.scale = mad, ...) {

  stopifnot(all(c('y', 'linear.predictors','data') %in% names(x)))

  # If the data is given, using x$data intead of model.frame avoids issue #39
  if (is.data.frame(x$data)) {
    themf <- model.frame(x$data, na.action=na.pass)
    z <- themf[,all.vars(as.formula(x$formula))[[1]]] # the explicit cast is for characters
  } else {
    themf <- model.frame(x$formula, na.action=na.pass)
    z <- model.response(themf)
  }
  lp <- scores(x, newdata=themf)

  # If z has any missingness, drop it from both z and lp
  lp <- lp[!is.na(z)]
  z <- z[!is.na(z)]

  pooled.sd <- if (is.null(standardization.scale)) {
    1
  } else {
    match_on_szn_scale(lp, z, standardization.scale)
  }
  lp.adj <- lp/pooled.sd

  match_on(lp.adj, within = within, caliper = caliper, z = z, ...)
}

match_on_szn_scale <- function(x, Tx, standardizer = mad, ...) {
  if (is.function(standardizer)) {
    sqrt(((sum(!Tx) - 1) * standardizer(x[!Tx])^2 +
          (sum(!!Tx) - 1) * standardizer(x[!!Tx])^2) / (length(x) - 2))
  } else if (is.numeric(standardizer)) {
    sqrt(((sum(!Tx) - 1) * standardizer^2 +
          (sum(!!Tx) - 1) * standardizer^2) / (length(x) - 2))
  } else {
    stop("Invalid standardizer")
  }
}

#' @details The \code{bigglm}
#' method works analogously to the \code{glm} method, but with \code{bigglm} objects, created by
#' the \code{bigglm} function from package \sQuote{biglm}, which can
#' handle bigger data sets than the ordinary glm function can.
#'
#' @method match_on bigglm
#' @rdname match_on-methods
match_on.bigglm <- function(x, within = NULL, caliper = NULL, data = NULL, standardization.scale = mad, ...) {
  if (is.null(data)) {
    stop("data argument is required for computing match_ons from bigglms")
  }

  if (!is.data.frame(data)) {
    stop("match_on doesn't understand data arguments that aren't data frames.")
  }

  theps <- scores(x, data, type = 'link', se.fit = FALSE)

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
    match_on_szn_scale(theps, z, standardizer=standardization.scale,...)
  }

  theps <- as.vector(theps / pooled.sd)
  names(theps) <- rownames(data)

  match_on(theps, within = within, caliper = caliper, z = z, ... )
}

# Note that Details for the glm method, above, refers to the below for discussion of computational
# benefits of calipers -- if that's changed here, adjust there accordingly.
#' @details The \code{numeric} method returns absolute differences between treated and control units'
#' values of \code{x}. If a caliper is specified, pairings with \code{x}-differences greater than it
#' are forbidden.  Conceptually, those distances are set to \code{Inf}; computationally, if either of
#' \code{caliper} and \code{within} has been specified then only information about permissible pairings
#' will be stored, so the forbidden pairings are simply omitted. Providing a \code{caliper} argument here,
#' as opposed to omitting it and afterward applying the \code{\link{caliper}} function, reduces
#' storage requirements and may otherwise improve performance, particularly in larger problems.
#'
#' For the numeric method, \code{x} must have names.
#' @method match_on numeric
#' @rdname match_on-methods
match_on.numeric <- function(x, within = NULL, caliper = NULL, data = NULL, z, ...) {

  if(missing(z) || is.null(z)) {
    stop("You must supply a treatment indicator, 'z', when using the numeric match_on method.")
  }

  if(length(z) != length(x)) {
    stop("The scores, 'x', and the treatment vector, 'z', must be the same length.")
  }
  if (!exists("cl")) cl <- match.call()

  z <- toZ(z)

  if(!is.null(caliper)) {
    if (length(caliper) > 1) {
      stop("Argument `caliper` must be a scalar value, not a vector.")
    }

    allowed <- scoreCaliper(x, z, caliper)

    if (!is.null(within)) {
      within <- within + allowed
    } else {
      within <- allowed
    }
  }

  f <- function(index, data, z) { abs(data[index[,1]] - data[index[,2]]) }

  tmp <- makedist(z, x, f, within)
  tmp@call <- cl
  return(tmp)
}

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
#' @method match_on InfinitySparseMatrix
match_on.InfinitySparseMatrix <- function(x, within = NULL, caliper = NULL, data = NULL, ...) {

  if(is.null(caliper)) { return(x) }

  return(x + optmatch::caliper(x, width = caliper))
} # just return the argument

#' @rdname match_on-methods
#' @method match_on matrix
match_on.matrix <- function(x, within = NULL, caliper = NULL, data = NULL, ...) {
  if(is.null(caliper)) { return(x) }

  return(x + optmatch::caliper(x, width = caliper))
} # just return the argument
