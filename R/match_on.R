################################################################################
# match_on: distance matrix creation functions
################################################################################

#' Create treated to control distances for matching problems
#'
#' A function with which to produce matching distances, for instance Mahalanobis
#' distances, propensity score discrepancies or calipers, or combinations
#' thereof, for \code{\link{pairmatch}} or \code{\link{fullmatch}} to
#' subsequently \dQuote{match on}.  Conceptually, the result of a call
#' \code{match_on} is a treatment-by-control matrix of distances.  Because these
#' matrices can grow quite large, in practice \code{match_on} produces either an
#' ordinary dense matrix or a special sparse matrix structure (that can make use
#' of caliper and exact matching constraints to reduce storage requirements).
#' Methods are supplied for these sparse structures,
#' \code{InfinitySparseMatrix}es, so that they can be manipulated and modified
#' in much the same way as dense matrices.
#'
#' \code{match_on} is generic. There are several supplied methods, all providing
#' the same basic output: a matrix (or similar) object with treated units on the
#' rows and control units on the columns. Each cell [i,j] then indicates the
#' distance from a treated unit i to control unit j. Entries that are \code{Inf}
#' are said to be unmatchable. Such units are guaranteed to never be in a
#' matched set. For problems with many \code{Inf} entries, so called sparse
#' matching problems, \code{match_on} uses a special data type that is more
#' space efficient than a standard R \code{matrix}.  When problems are not
#' sparse (i.e. dense), \code{match_on} uses the standard \code{matrix} type.
#'
#' \code{match_on} methods differ on the types of arguments they take, making
#' the function a one-stop location of many different ways of specifying
#' matches: using functions, formulas, models, and even simple scores. Many of
#' the methods require additional arguments, detailed below. All methods take a
#' \code{within} argument, a distance specification made using
#' \code{\link{exactMatch}} or \code{\link{caliper}} (or some additive
#' combination of these or other distance creating functions). All
#' \code{match_on} methods will use the finite entries in the \code{within}
#' argument as a guide for producing the new distance. Any entry that is
#' \code{Inf} in \code{within} will be \code{Inf} in the distance matrix
#' returned by \code{match_on}. This argument can reduce the processing time
#' needed to compute sparse distance matrices.
#'
#' Details for each particular first type of argument follow:
#'
#' @param x A model formula, fitted glm or other object implicitly specifying a distance; see blurbs on specific methods in Details.
#' @param within A valid distance specification, such as the result of
#'   \code{\link{exactMatch}} or \code{\link{caliper}}. Finite entries indicate
#'   which distances to create. Including this argument can significantly speed
#'   up computation for sparse matching problems. Specify this filter either via
#'   \code{within} or via \code{strata} elements of a formula; mixing these
#'   methods will fail.
#' @param caliper The width of a caliper to use to exclude treated-control pairs
#'   with values greater than the width. For some methods, there may be a speed
#'   advantage to passing a width rather than using the \code{\link{caliper}}
#'   function on an existing distance specification.
#' @param exclude A list of units (treated or control) to exclude from the
#'     \code{caliper} argument (if supplied).
#' @param data An optional data frame.
#' @param ... Other arguments for methods.
#' @return A distance specification (a matrix or similar object) which is
#'   suitable to be given as the \code{distance} argument to
#'   \code{\link{fullmatch}} or \code{\link{pairmatch}}.
#' @seealso \code{\link{fullmatch}}, \code{\link{pairmatch}},
#'   \code{\link{exactMatch}}, \code{\link{caliper}}
#' @references P.~R. Rosenbaum and D.~B. Rubin (1985), \sQuote{Constructing a
#'   control group using multivariate matched sampling methods that incorporate
#'   the propensity score}, \emph{The American Statistician}, \bold{39} 33--38.
#' @export
#' @example inst/examples/match_on.R
#' @docType methods
#' @rdname match_on-methods
match_on <- function(x, within = NULL, caliper = NULL, exclude = NULL, data=NULL, ...) {
  # if x does not exist then print helpful error msg
  x_str <- deparse(substitute(x))
  data_str <- deparse(substitute(data))
  tryCatch(x, error = function(e) {
    stop(missing_x_msg(x_str, data_str, ...))})

  UseMethod("match_on")
}

#' @details \bold{First argument (\code{x}): \code{glm}.} The model is assumed to be
#'   a fitted propensity score model. From this it extracts distances on the
#'   \emph{linear} propensity score: fitted values of the linear predictor, the
#'   link function applied to the estimated conditional probabilities, as opposed
#'   to the estimated conditional probabilities themselves (Rosenbaum & Rubin,
#'   1985).  For example, a logistic model (\code{glm} with
#'   \code{family=binomial()}) has the logit function as its link, so from such
#'   models \code{match_on} computes distances in terms of logits of the
#'   estimated conditional probabilities, i.e. the estimated log odds.
#'
#'   Optionally these distances are also rescaled. The default is to rescale, by
#'   the reciprocal of an outlier-resistant variant of the pooled s.d. of
#'   propensity scores; see \code{\link{standardization_scale}}.  (The
#'   \code{standardization.scale} argument of this function can be used to
#'   change how this dispersion is calculated, e.g. to calculate an ordinary not
#'   an outlier-resistant s.d.; it will be passed down
#'   to \code{standardization_scale} as its \code{standardizer} argument.)
#'   To skip rescaling, set argument \code{standardization.scale} to 1.
#'   The overall result records
#'   absolute differences between treated and control units on linear, possibly
#'   rescaled, propensity scores.
#'
#'   In addition, one can impose a caliper in terms of these distances by
#'   providing a scalar as a \code{caliper} argument, forbidding matches between
#'   treatment and control units differing in the calculated propensity score by
#'   more than the specified caliper.  For example, Rosenbaum and Rubin's (1985)
#'   caliper of one-fifth of a pooled propensity score s.d. would be imposed by
#'   specifying \code{caliper=.2}, in tandem either with the default rescaling
#'   or, to follow their example even more closely, with the additional
#'   specification \code{standardization.scale=sd}. Propensity calipers are
#'   beneficial computationally as well as statistically, for reasons indicated
#'   in the below discussion of the \code{numeric} method.
#'
#'   One can also specify exactMatching criteria by using \code{strata(foo)} inside
#'   the formula to build the \code{glm}. For example, passing
#'   \code{glm(y ~ x + strata(s))} to \code{match_on} is equivalent to passing
#'   \code{within=exactMatch(y ~ strata(s))}. Note that when combining with
#'   the \code{caliper} argument, the standard deviation used for the caliper will be
#'   computed across all strata, not within each strata.
#'
#'   If data used to fit the glm have missing values in the left-hand side
#'   (dependent) variable, these observations are omitted from the output of
#'   match_on.  If there are observations with missing values in right hand
#'   side (independent) variables, then a re-fit of the model after imputing
#'   these variables using a simple scheme and adding indicator variables of
#'   missingness will be attempted, via the \code{\link{scores}} function.
#'
#' @param standardization.scale Function for rescaling of \code{scores(x)}, or
#'   \code{NULL}; defaults to \code{mad}. (See Details.)
#' @seealso \code{\link{scores}}
#' @method match_on glm
#' @rdname match_on-methods
#' @export
match_on.glm <- function(x,
                         within = NULL,
                         caliper = NULL,
                         exclude = NULL,
                         data = NULL,
                         standardization.scale = NULL,
                         ...) {

  stopifnot(all(c('y','data') %in% names(x)))

  # If the data is given, using x$data intead of model.frame avoids issue #39
  if (is.data.frame(x$data)) {
    themf <- model.frame(x$data, na.action=na.pass)
    z <- themf[,all.vars(as.formula(x$formula))[[1]]]
    # the explicit cast is for characters
  } else {
    themf <- model.frame(x$formula, na.action=na.pass)
    z <- model.response(themf)
  }
  lp <- scores(x, newdata=themf, ...)

  # If z has any missingness, drop it from both z and lp
  lp <- lp[!is.na(z)]
  z <- z[!is.na(z)]

  pooled.sd <-
    standardization_scale(lp, trtgrp=z, standardization.scale,
                       svydesign_ = x$'survey.design')

  lp.adj <- lp/pooled.sd

  if (!is.null(attr(terms(formula(x), special = "strata", data = data),
                    "specials")$strata)) {
    newwithin <- makeWithinFromStrata(formula(x), data)
    if (is.null(within)) {
      within <- newwithin$within
    } else {
      within <- newwithin$within + within
    }
  }


  out <-match_on(lp.adj,
                 within = within,
                 caliper = caliper,
                 exclude = exclude,
                 z = z,
                 ...)
  out@call <- match.call()
  return(out)
}

#' Pooled Dispersion for a Numeric Variable
#'
#' Dispersion as pooled across a treatment and a control group. By default,
#' the measure of dispersion calculated within each group is not the
#' ordinary standard deviation but rather the robust alternative
#' provided by \code{stats::mad}.
#'
#' @title pooled dispersion for a numeric variable
#' @param x numeric variable
#' @param trtgrp logical or numeric. If numeric, coerced to logical via \code{!}
#' @param standardizer function or numeric of length 1
#' @return numeric of length 1
#' @keywords internal
match_on_szn_scale <- function(x, trtgrp, standardizer = mad) {
  if (is.function(standardizer)) {
    sqrt(((sum(!trtgrp) - 1) * standardizer(x[!trtgrp])^2 +
          (sum(!!trtgrp) - 1) * standardizer(x[!!trtgrp])^2) / (length(x) - 2))
  } else if (is.numeric(standardizer)) {
    sqrt(((sum(!trtgrp) - 1) * standardizer^2 +
          (sum(!!trtgrp) - 1) * standardizer^2) / (length(x) - 2))
  } else {
    stop("Invalid standardizer")
  }
}

#' @details \bold{First argument (\code{x}): \code{bigglm}.} This method works
#'   analogously to the \code{glm} method, but with \code{bigglm} objects,
#'   created by the \code{bigglm} function from package \sQuote{biglm}, which
#'   can handle bigger data sets than the ordinary glm function can.
#'
#' @method match_on bigglm
#' @rdname match_on-methods
#' @export
match_on.bigglm <- function(x,
                            within = NULL,
                            caliper = NULL,
                            exclude = NULL,
                            data = NULL,
                            standardization.scale = NULL,
                            ...) {
  if (is.null(data)) {
    stop("data argument is required for computing match_ons from bigglms")
  }

  if (!is.data.frame(data)) {
    stop("match_on doesn't understand data arguments that aren't data frames.")
  }

  theps <- scores(x, data, type = 'link', se.fit = FALSE)

  if (length(theps) != dim(data)[1]) {
    stop("predict.bigglm() returns a vector of the wrong length; are there missing values in data?")
  }

  # this makes heavy use of the bigglm terms object, the original formula
  # if this implementation detail changes in later versions of bigglm,
  # look here for the problem.

  Data <-  model.frame(x$terms, data = data)
  z <- Data[, 1]
  pooled.sd <- standardization_scale(theps,
                                     trtgrp = z,
                                     standardizer = standardization.scale)

  theps <- as.vector(theps / pooled.sd)
  names(theps) <- rownames(data)

  out <- match_on(theps,
                  within = within,
                  caliper = caliper,
                  exclude = exclude,
                  z = z,
                  ... )
  out@call <- match.call()
  return(out)
}

#' @details \bold{First argument (\code{x}): \code{formula}.} The formula must have
#'   \code{Z}, the treatment indicator (\code{Z=0} indicates control group,
#'   \code{Z=1} indicates treatment group), on the left hand side, and any
#'   variables to compute a distance on on the right hand side. E.g. \code{Z ~ X1
#'   + X2}. The Mahalanobis distance is calculated as the square root of d'Cd,
#'   where d is the vector of X-differences on a pair of observations and C is an
#'   inverse (generalized inverse) of the pooled covariance of Xes. (The pooling
#'   is of the covariance of X within the subset defined by \code{Z==0} and
#'   within the complement of that subset. This is similar to a Euclidean
#'   distance calculated after reexpressing the Xes in standard units, such that
#'   the reexpressed variables all have pooled SDs of 1; except that it addresses
#'   redundancies among the variables by scaling down variables contributions in
#'   proportion to their correlations with other included variables.)
#'
#'   Euclidean distance is also available, via \code{method="euclidean"}, and
#'   ranked, Mahalanobis distance, via \code{method="rank_mahalanobis"}.
#'
#'   The treatment indicator \code{Z} as noted above must either be numeric
#'   (1 representing treated units and 0 control units) or logical
#'   (\code{TRUE} for treated, \code{FALSE} for controls). (Earlier versions of
#'   the software accepted factor variables and other types of numeric variable; you
#'   may have to update existing scripts to get them to run.)
#'
#'
#'   As an alternative to specifying a \code{within} argument, when \code{x} is
#'   a formula, the \code{strata} command can be used inside the formula to specify
#'   exact matching. For example, rather than using \code{within=exactMatch(y ~
#'   z, data=data)}, you may update your formula as \code{y ~ x + strata(z)}. Do
#'   not use both methods (\code{within} and \code{strata} simultaneously. Note
#'   that when combining with the \code{caliper} argument, the standard
#'   deviation used for the caliper will be computed across all strata, not
#'   separately by stratum.
#'
#'   A unit with NA treatment status (\code{Z}) is ignored and will not be included in the distance output.
#'  Missing values in variables on the right hand side of the formula are handled as follows. By default
#' \code{match_on} will (1) create a matrix of distances between observations which
#' have only valid values for **all** covariates and then (2) append matrices of Inf values
#' for distances between observations either of which has a missing values on any of the right-hand-side variables.
#' (I.e., observations with missing values are retained in the output, but
#' matches involving them are forbidden.)
#'
#' @param subset A subset of the data to use in creating the distance
#'   specification.
#' @param method A string indicating which method to use in computing the
#'   distances from the data.  The current possibilities are
#'   \code{"mahalanobis", "euclidean"} or \code{"rank_mahalanobis"}.
#' @method match_on formula
#' @rdname match_on-methods
#' @importFrom stats contrasts
#' @export
match_on.formula <- function(x,
                             within = NULL,
                             caliper = NULL,
                             exclude = NULL,
                             data = NULL,
                             subset = NULL,
                             method = "mahalanobis",
                             ...) {
  if (length(x) != 3) {
    stop("Formula must have a left hand side.")
  }

  if (grepl(".", as.character(x)[3], fixed = TRUE) &
      grepl("strata(", as.character(x)[3], fixed = TRUE)) {
    stop("Cannot use . expansion in formula when defining strata.")
  }

  mf <- match.call(expand.dots = FALSE)

  m <- match(c("x", "data", "subset"), # maybe later add "na.action"
             names(mf), 0L)
  mf <- mf[c(1L, m)]

  rm(m)

  # Get terms, with flags for any strata (to handle) or caliper (to break)
  t <- terms(x, specials = c("strata","caliper"), data = data)
  if (!is.null(attr(t, "specials")$strata)) {
    newwithin <- makeWithinFromStrata(x, data)
    x <- newwithin$x
    mf[[2]] <- x
    if (is.null(within)) {
      within <- newwithin$within
    } else {
      within <- newwithin$within + within
    }
  }

  # #114 - Catching user input of caliper in formula
  if (!is.null(attr(t, "specials")$caliper)) {
    calnames <- rownames(attr(t, 'factors'))[attr(t, "specials")$caliper]
    withinname <- as.list(match.call())$within
    if (is.null(withinname)) {
      withinname <- ""
    } else {
      withinname <- paste0(deparse(substitute(withinname)), " + ")
    }

    stop(paste0("Calipers should be applied via the `within` argument instead of in the formula.\n",
                "\tE.g. `within = ", ifelse(is.null(withinname), "", withinname),
                paste(calnames, collapse = " + "), "`"))
  }

  names(mf)[names(mf) == "x"] <- "formula"
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  if (dim(mf)[2] < 2) {
    stop("Formula must have a right hand side with at least one variable.")
  }

  tmpz <- toZ(mf[,1])
  tmpn <- rownames(mf)

  # If there are any NA treated members, throw them away first
  mf <- mf[!is.na(tmpz), ]
  tmpn <- tmpn[!is.na(tmpz)]
  tmpz <- tmpz[!is.na(tmpz)]

  mf <- na.omit(mf)

  dropped.t <- setdiff(tmpn[tmpz],  rownames(mf))
  dropped.c <- setdiff(tmpn[!tmpz], rownames(mf))

  # Switching contrts function, see #220
  data <- subset(model.matrix(x, mf,
                              contrasts.arg =
                                lapply(Filter(is.factor, mf),
                                       function(x) {
                                         stats::contrasts(x, contrasts = FALSE)/sqrt(2)
                                       })),
                 TRUE, -1) # drop the intercept

  z <- toZ(mf[,1])
  names(z) <- rownames(mf)

  if(is.character(method)){
    methodname <- method
  } else {
    methodname <- as.character(class(method))
  }

  which.method <- pmatch(methodname, c("mahalanobis", "euclidean", "rank_mahalanobis", "function"), 4)
  tmp <- switch(which.method,
		makedist(z, data, compute_mahalanobis, within),
		makedist(z, data, compute_euclidean, within),
    makedist(z, data, compute_rank_mahalanobis, within),
    {
      warning("Passing a user-defined `method` to `match_on.formula` is not supported and results are not guaranteed. User-defined distances should use `match_on.function` instead.")
      makedist(z, data, match.fun(method), within)
    }
		)
  rm(mf)

  if (length(dropped.t) > 0 || length(dropped.c)) {
    if(!is(tmp, "BlockedInfinitySparseMatrix"))
    {
      tmp <- as.InfinitySparseMatrix(tmp)
    }
    tmp@rownames <- c(tmp@rownames, dropped.t)
    tmp@colnames <- c(tmp@colnames, dropped.c)
    tmp@dimension <- c(length(tmp@rownames), length(tmp@colnames))
    if(is(tmp, "BlockedInfinitySparseMatrix"))
    {
      tmp@groups <- unlist(list(tmp@groups, within@groups[!names(within@groups) %in% names(tmp@groups)]))
    }

  }

  cl <- match.call()

  if (is.null(caliper)) {
    tmp@call <- cl
    return(tmp)
  }

  tmp <- tmp + optmatch::caliper(tmp, width = caliper, exclude = exclude)
  tmp@call <- cl

  return(tmp)
}

# Internal function. Given a formula (which includes strata elements)
# and the data frame, return an updated formula without the strata elements
# and an exactMatch using the strata.
#
# Note: Only use the \code{x} returned in match_on.formula; in match_on.glm we do
# not remove the strata variables from the model.
makeWithinFromStrata <- function(x, data)
{
  xs <- findStrata(x, data)

  em <- unlist(lapply(strsplit(xs$strata, "\\(|)|,"), "[", -1))
  lhs <- paste(xs$newx[[2]], collapse="")
  within <- exactMatch(as.formula(paste(lhs, "~", paste(em, collapse="+"))),
                             data=data)
  return(list(x= xs$newx, within=within))
}

findStrata <- function(x, data) {

  t <- terms(x, specials = "strata", data = data)

  strata <- rownames(attr(t, "factors"))[attr(t, "specials")$strata]

  if (length(strata) > 0) {
    x <- update(x, as.formula(paste("~ . - ", paste(strata, collapse="-"))))
    return(list(newx = x, strata = strata))
  }

  return(list(newx = x, strata = NULL))

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

compute_mahalanobis <- function(index, data, z) {
  if (!all(is.finite(data))) stop("Infinite or NA values detected in data for Mahalanobis computations.")

    mt <- cov(data[z, ,drop = FALSE]) * (sum(z) - 1) / (length(z) - 2)
    mc <- cov(data[!z, ,drop = FALSE]) * (sum(!z) - 1) / (length(!z) - 2)
    # Addressing #168
    if (sum(z) == 1) mt <- 0
    if (sum(!z) == 1) mc <- 0
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

    return(mahalanobisHelper(data, index, inv.scale.matrix))
}


compute_euclidean <- function(index, data, z) {

  if (!all(is.finite(data))) {
    stop("Infinite or NA values detected in data for distance computations.")
  }

  return(mahalanobisHelper(data, index, diag(ncol(data))))
}


compute_rank_mahalanobis <- function(index, data, z) {
    if (!all(is.finite(data))) {
        stop("Infinite or NA values detected in data for Mahalanobis computations.")
    }

    if (is.null(index)) return(sqrt(r_smahal(NULL, data, z)))

    if (is.null(rownames(data)) | !all(index %in% rownames(data)))
        stop("data must have row names matching index")

    # begin workaround solution to #128
    all_treated <- rownames(data)[as.logical(z)]
    all_control <- rownames(data)[!z]
    all_indices <- expand.grid(all_treated, all_control,
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    all_indices <- paste(all_indices[[1]], all_indices[[2]], sep="%@%")
    short_indices <- paste(index[,1], index[,2], sep="%@%")
    indices <- match(short_indices, all_indices)
    if (any(is.na(indices))) stop("Unanticipated problem. (Make sure row names of data don't use the string '%@%'.)")
    # Now, since `r_smahal` is ignoring its `index` argument anyway:
    rankdists <- sqrt(r_smahal(NULL, data, z))
    rankdists <- rankdists[indices]
    return(rankdists)
}

#' @details \bold{First argument (\code{x}): \code{function}.} The passed function
#'   must take arguments: \code{index}, \code{data}, and \code{z}. The
#'   \code{data} and \code{z} arguments will be the same as those passed directly
#'   to \code{match_on}. The \code{index} argument is a matrix of two columns,
#'   representing the pairs of treated and control units that are valid
#'   comparisons (given any \code{within} arguments). The first column is the row
#'   name or id of the treated unit in the \code{data} object. The second column
#'   is the id for the control unit, again in the \code{data} object. For each of
#'   these pairs, the function should return the distance between the treated
#'   unit and control unit.  This may sound complicated, but is simple to
#'   use. For example, a function that returned the absolute difference between
#'   two units using a vector of data would be \code{ f <- function(index, data,
#'   z) { abs(data[index[,1]] - data[index[,2]]) } }.  (Note: This simple case is
#'   precisely handled by the \code{numeric} method.)
#'
#' @param z A logical or binary vector indicating treatment and control for each
#'  unit in the study. TRUE or 1 represents a treatment unit, FALSE of 0 represents
#'  a control unit. Any unit with NA treatment status will be excluded from the
#'  distance matrix.
#' @method match_on function
#' @rdname match_on-methods
#' @export
match_on.function <- function(x,
                              within = NULL,
                              caliper = NULL,
                              exclude = NULL,
                              data = NULL,
                              z = NULL,
                              ...) {

  if (is.null(data) | is.null(z)) {
    stop("Data and treatment indicator arguments are required.")
  }

  theFun <- match.fun(x)

  data <- subset(data, !is.na(z))
  z <- z[!is.na(z)]

  tmp <- makedist(z, data, theFun, within)
  cl <- match.call()

  if (is.null(caliper)) {
    tmp@call <- cl
    return(tmp)
  }

  tmp <- tmp + optmatch::caliper(tmp, width = caliper, exclude = exclude)
  tmp@call <- cl
  return(tmp)
}

# Note that Details for the glm method, above, refers to the below for discussion of computational
# benefits of calipers -- if that's changed here, adjust there accordingly.
#' @details \bold{First argument (\code{x}): \code{numeric}.} This returns
#'   absolute differences between treated and control units' values of \code{x}.
#'   If a caliper is specified, pairings with \code{x}-differences greater than
#'   it are forbidden.  Conceptually, those distances are set to \code{Inf};
#'   computationally, if either of \code{caliper} and \code{within} has been
#'   specified then only information about permissible pairings will be stored,
#'   so the forbidden pairings are simply omitted. Providing a \code{caliper}
#'   argument here, as opposed to omitting it and afterward applying the
#'   \code{\link{caliper}} function, reduces storage requirements and may
#'   otherwise improve performance, particularly in larger problems.
#'
#'   For the numeric method, \code{x} must have names. If \code{z} is named
#'   it must have the same names as \code{x}, though it allows for a different
#'   ordering of names. \code{x}'s name ordering is considered canonical.
#' @method match_on numeric
#' @rdname match_on-methods
#' @export
match_on.numeric <- function(x, within = NULL, caliper = NULL, exclude = NULL, data = NULL, z, ...) {

  if(missing(z) || is.null(z)) {
    stop("You must supply a treatment indicator, 'z', when using the numeric match_on method.")
  }

  if(length(z) != length(x)) {
    stop("The scores, 'x', and the treatment vector, 'z', must be the same length.")
  }

  z <- toZ(z)

  # #189
  if (is.null(names(z))) {
    names(z) <- names(x)
  } else {
    # we already know length(x) == length(z)
    if (!all(names(x) %in% names(z))) {
      stop("names of x and z must be the same")
    }
    # Treat x's name ordering as canonical
    z <- z[names(x)]
  }

  x <- x[!is.na(z)]
  z <- z[!is.na(z)]

  missingX <- is.na(x)
  rnms <- names(z)
  dropped.t <- rnms[missingX & z]
  dropped.c <- rnms[missingX & !z]

  z <- z[!missingX]
  x <- x[!missingX]

  if(!is.null(caliper)) {
    if (length(caliper) > 1) {
      stop("Argument `caliper` must be a scalar value, not a vector.")
    }

    if (!is.null(exclude)) {
        x.tmp <- x[!names(x) %in% exclude]
        z.tmp <- z[!names(z) %in% exclude]
        allowed <- scoreCaliper(x.tmp, z.tmp, caliper)

        n1 <- sum(z)
        n0 <- sum(!z)

        z.exclude <- z[names(z) %in% exclude]

        n1_exclude <- sum(z.exclude)
        n0_exclude <- sum(!z.exclude)

        if (n1_exclude > 0) {
            n1_exclude_idx <- (n1 - n1_exclude + 1):n1
            n1_names <- names(z.exclude[z.exclude])
        } else {
            n1_exclude_idx <- NULL
            n1_names <- NULL
        }

        if (n0_exclude > 0) {
            n0_exclude_idx <- (n0 - n0_exclude + 1):n0
            n0_names <- names(z.exclude[!z.exclude])
        } else {
            n0_exclude_idx <- NULL
            n0_names <- NULL
        }

        allowed <- makeInfinitySparseMatrix(
            c(allowed@.Data,
              rep(0, n1_exclude * n0 + n0_exclude * n1 - (n1_exclude * n0_exclude))),
            rows = c(allowed@rows,
                     rep(n1_exclude_idx, n0),
                     rep(1:(n1 - n1_exclude), n0_exclude)),
            cols = c(allowed@cols,
                     rep(1:n0, each = n1_exclude),
                     rep(n0_exclude_idx, each = (n1 - n1_exclude))),
            rownames = c(allowed@rownames, n1_names),
            colnames = c(allowed@colnames, n0_names))

    } else {
        allowed <- scoreCaliper(x, z, caliper)
    }

    if (!is.null(within)) {
      within <- within + allowed
    } else {
      within <- allowed
    }
  }

  f <- function(index, data, z) { abs(data[index[,1]] - data[index[,2]]) }

  tmp <- makedist(z, x, f, within)
  # we dropped units with missing x values, now we reappply them

  if (length(dropped.t) > 0 || length(dropped.c)) {
    tmp <- as.InfinitySparseMatrix(tmp)
    tmp@rownames <- c(tmp@rownames, dropped.t)
    tmp@colnames <- c(tmp@colnames, dropped.c)
    tmp@dimension <- c(length(tmp@rownames), length(tmp@colnames))
  }

  tmp@call <- match.call()
  return(tmp)
}

#' (Internal) Helper function to create an InfinitySparseMatrix from a set of
#' scores, a treatment indicator, and a caliper width.
#'
#' @param x The scores, a vector indicating the 1-D location of each
#'   unit.
#' @param z The treatment assignment vector (same length as \code{x})
#' @param caliper The width of the caliper with respect to the scores
#'   \code{x}.
#' @return An \code{InfinitySparseMatrix} object, suitable to be
#'   passed to \code{\link{match_on}} as an \code{within} argument.
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

#' @details \bold{First argument (\code{x}): \code{matrix} or \code{InfinitySparseMatrix}.} These just return their
#'   arguments as these objects are already valid distance specifications.
#'
#' @rdname match_on-methods
#' @method match_on InfinitySparseMatrix
#' @export
match_on.InfinitySparseMatrix <- function(x,
                                          within = NULL,
                                          caliper = NULL,
                                          exclude = NULL,
                                          data = NULL, ...) {

  if (!is.null(data)) {
    x <- subset(x,
                subset = rownames(x) %in% rownames(data),
                select = colnames(x) %in% rownames(data))
  }

  cl <- match.call()

  if (is.null(caliper) & is.null(within)) {
    x@call <- cl
    return(x) # just return the argument
  }

  if (is.null(within)) { # If we're here, caliper arg is non-null
    out <- x + optmatch::caliper(x, width = caliper, exclude = exclude)
    out@call <- cl
    return(out)
  }

  ## If we're here, within is non-null, but caliper may or may not be null
  within  <-  within * 0
  if (!is.null(caliper)) {
    within  <- within +
      optmatch::caliper(x, width = caliper, exclude = exclude)
  }
  out <- x + within
  out@call <- cl
  return(out)
}

#' @rdname match_on-methods
#' @method match_on matrix
#' @export
match_on.matrix <- function(x,
                            within = NULL,
                            caliper = NULL,
                            exclude = NULL,
                            data = NULL,
                            ...) {

  if (!is.null(data)) {
    x <- subset(x,
                subset = rownames(x) %in% rownames(data),
                select = intersect(colnames(x), rownames(data)))
  }

  cl <- match.call()

  if (is.null(caliper) & is.null(within)) {
    if (is(x, "DenseMatrix")) {
      # If `x` is an actual matrix, don't store the call
      x@call <- cl
    }
    return(x) # just return the argument
  }

  if (is.null(within)) { # If we're here, caliper arg is non-null
    out <- x + optmatch::caliper(x, width = caliper, exclude = exclude)
    out@call <- cl
    return(out)
  }

  ## If we're here, within is non-null, but caliper may or may not be null
  within <- within * 0
  if (!is.null(caliper)) {
    within <- within + optmatch::caliper(x, width = caliper, exclude = exclude)
  }
  out <- x + within
  out@call <- cl
  return(out)
}

#' pooled dispersion for a numeric variable
#'
#' Dispersion as pooled across a treatment and a control group. By default,
#' the measure of dispersion calculated within each group is not the
#' ordinary standard deviation as in \code{stats::sd} but rather the robust alternative
#' encoded in \code{stats::mad}.  The dispersion measurements are combined
#' by squaring, averaging with weights proportional to one minus the sizes of
#' the groups and then taking square roots.  Used in \code{\link{match_on.glm}}.
#'
#' A non-NULL \code{svydesign_} parameter indicates that the dispersion
#' calculations are to be made respecting the weighting scheme implicit in
#' that \code{survey.design2} object. If \code{standardizer} is \code{NULL},
#' one gets a calculation in the style of \code{stats::mad} but with weights,
#' performed by \code{optmatch:::svy_sd}; for a pooling of weighted standard
#' deviations, one would pass a non-\code{NULL} \code{svydesign_} parameter along
#' with \code{standardizer=optmatch:::svy_sd}.
#' (More generally, the provided \code{standardizer}
#' function should accept as a sole argument a \code{survey.design2} object,
#' with \code{nrows(svydesign_$variables)} equal to the lengths of \code{x} and
#' \code{trtgrp}.  This object is expected to carry a numeric variable \sQuote{\code{x}},
#' and the \code{standardizer} function is to return the dispersion of this variable.)
#'
#' @param x numeric variable
#' @param trtgrp logical or numeric. If numeric, coerced to logical via \code{!}
#' @param standardizer function, \code{NULL} or numeric of length 1
#' @param svydesign_ ordinarily \code{NULL}, but may also be a
#'   \code{survey.design2}; see Details.
#' @return numeric of length 1
#' @export
#' @keywords internal
standardization_scale <- function(x, trtgrp, standardizer = NULL, svydesign_=NULL)
    {
    stopifnot(is.null(svydesign_) || is(svydesign_, "survey.design2"),
              is.null(standardizer) || is.function(standardizer) || is.numeric(standardizer)
              )
    if (is.numeric(standardizer))
    {
        if (length(standardizer)>1)
            warning("Multiple element standardizer, only the first is used")
        return(standardizer)
    }
    n_c <- sum(!trtgrp)
    n <- length(x)
    n_t <- n - n_c
    if (is.null(svydesign_))
        {
	if (is.null(standardizer)) standardizer <- stats::mad
	s_c <- standardizer(x[!trtgrp])
	s_t <- standardizer(x[as.logical(trtgrp)])
    } else {
        if (is.null(standardizer)) standardizer <- svy_mad
        des <- update(svydesign_, x=x, trtgrp=as.logical(trtgrp))
        des_t <- subset(des, trtgrp)
        des_c <- subset(des, !trtgrp)
        s_t <- standardizer(des_t)
        s_c <- standardizer(des_c)
    }
    s2_t <- s_t^2
    s2_c <- s_c^2
    sqrt(((n_t - 1) * s2_t +
          (n_c - 1) * s2_c) / (n - 2))
}

#' @keywords internal
svy_mad <- function(design)
{
  if (requireNamespace("survey", quietly = TRUE)) {
        med <- survey::svyquantile(~x, design, 0.5)[[1]][1]

        design <- update(design,
                        abs_dev=abs( design$variable$x - med )
                        )
        mad <- survey::svyquantile(~abs_dev, design, 0.5)[[1]][1]
        constant <- formals(stats::mad)$constant
        s2_t <- constant * mad
        return(s2_t)
  } else {
    stop("'survey' package must be installed")
  }
}
#' @keywords internal
svy_sd <- function(design)
{
  if (requireNamespace("survey", quietly = TRUE)) {
        var_ <- survey::svyvar(~x, design)
        return(sqrt(unname(var_)[1]))
  } else {
    stop("'survey' package must be installed")
  }
}

#' This method quells a warning when \code{optmatch::scores()}
#' is applied to a svyglm object.
#' @method model.frame svyglm
#' @keywords internal
model.frame.svyglm <- function (formula, ...)
{
    ans <- get_all_vars(formula, formula[["survey.design"]][["variables"]])
    attr(ans, "terms") <- terms(formula)
    ans
}
