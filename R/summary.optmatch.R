#' Display matching related statistics
#'
#' The summary function quantifies \code{optmatch} objects on the effective sample
#' size, the distribution of distances between matched units, and how well the
#' match reduces average differences.
#'
#' @param object The \code{optmatch} object to summarize.
#' @param propensity.model An optional propensity model (the result of
#'   a call to \code{glm}) to use when summarizing the match. Using
#'   the \code{RItools} package, an additional chi-squared test will
#'   be performed on the average differences between treated and
#'   control units on each variable used in the model. See the
#'   \code{xBalance} function in the \code{RItools} package for more
#'   details.
#' @param ... Additional arguments to pass to \code{xBalance} when
#'   also passing a propensity model.
#' @param min.controls To minimize the the display of a groups with
#'   many treated and few controls, all groups with more than 5
#'   treated units will be summarized as \dQuote{5+}. This is the
#'   reciprocal of the default value (1/5 = 0.2). Lower this value to
#'   see more groups.
#' @param max.controls Like \code{min.controls} sets maximum group
#'   sized displayed with respect to the number of controls. Raise
#'   this value to see more groups.
#' @param quantiles A points in the ECDF at which the distances
#'   between units will be displayed.
#' @return \code{optmatch.summary}
#' @seealso \code{\link{print.optmatch}}
#' @method summary optmatch
#' @rdname optmatch
#' @importFrom RItools xBalance
#' @export
summary.optmatch <- function(object,
                             propensity.model = NULL, ...,
                             min.controls=.2, max.controls=5,
                             quantiles=c(0,.5, .95, 1)
                             )
{
# Things to display in the summary method:
## effective sample size -- stratumStructure
## outlying propensity-score distances -- matched.distances()
## overall balance -- xBalance()
  so <- list()
  so$thematch <- object
  mfd <- is.na(object)
  if (all(mfd))
    {
      class(so) <- "summary.optmatch"
      so$matching.failed <- table(attr(object, "subproblem"), attr(object, "contrast.group"))
      dimnames(so$matching.failed)[[2]] <- c("z==0", "z==1")
      so$warnings <- c(so$warnings,
                       list("Matching failed.  (Restrictions impossible to meet?)\nEnter ?matchfailed for more info.")
                       )
      return(so)
    }
  subprobs <- attr(object, "subproblem")
  match.succeed <- tapply(mfd, subprobs, function(x) !all(x))
  so$matching.failed <- table(attr(object, "subproblem"), attr(object, "contrast.group"), exclude=names(match.succeed)[match.succeed])
  if (prod(dim(so$matching.failed)) == 0) {
    so$matching.failed <- NULL
  } else {
    dimnames(so$matching.failed)[[2]] <- c("z==0", "z==1")
  }
  so$matched.set.structures <- stratumStructure(object,min.controls=min.controls,max.controls=max.controls)
  so$effective.sample.size <- attr(so$matched.set.structures, "comparable.num.matched.pairs")

  # Per issue #106, allow the object to not carry a `matched.distances` attribute.
  if (!is.null(attr(object, "matched.distances"))) {
    matchdists <- attr(object, "matched.distances")[levels(object[!mfd, drop=TRUE])]
    matchdists <- unlist(matchdists)
    so$total.distance <- sum(matchdists)
    so$total.tolerances <- sum(unlist(attr(object, "exceedances")))
    so$matched.dist.quantiles <- quantile(matchdists, prob=quantiles)
  }

  if(!is.null(propensity.model) &&
     inherits(propensity.model, "glm")) {

    # users must save the model for reliable behavior.
    # we warn, instead of an error, but the user may get an error
    # from model.frame or later
    if(is.null(propensity.model$model)) {
      warning("This propensity seems to have been fit with 'model=FALSE'.\nI'm reconstructing the data set as best I can, but I might fail,\nor get a different data set than the one behind propensity.model.\nTo be sure, re-fit your propensity.model with 'model = TRUE'.")
    }

    # we need to handle the different ways of creating glm objects
    # 1: glm(Z ~ f(X), data = mydata)
    # 2: glm(Z ~ f(X)) # uses environment
    # 3: glm(fill.NAs(Z ~ f(x), data = mydata)) # This is no longer explicitly supported, see #103
    # each stores its data in different places.

    modelData <- NULL # be explicit for safety
    if (!is.null(propensity.model$data)) {
      # cases 1 and 2
      na.behavior <- FALSE
      modelData <- get_all_vars(propensity.model, data = propensity.model$data)

    } else {
      # case 3
      modelData <- propensity.model$model
      na.behavior <- TRUE
    }

    if (is.null(modelData)) {
      stop("summary.optmatch does not know how to process this type of model. Please file a bug report at https://github.com/markmfredrickson/optmatch/issues showing how you created your glm model.")
    }

    strata <- object[!mfd, drop=TRUE]
    data <- modelData[!mfd,]

    if (length(strata) != dim(data)[1]) {
      stop("'summary' method unable to recreate data. Consider passing 'data' argument to 'pairmatch' or 'fullmatch'.")
    }

    so$balance <- RItools::xBalance(fmla = formula(propensity.model),
                           strata = strata,
                           data = data,
                           report = c('adj.means', 'z.scores', 'chisquare.test'),
                           na.rm = na.behavior)

  } else if (!is.null(propensity.model)) so$warnings <-
    c(so$warnings,
      list("For covariate balance information, load the RItools package and\npass a (glm) propensity model to summary() as a second argument.")
      )

  class(so) <- "summary.optmatch"
  so
}

#' @export
print.summary.optmatch <- function(x,  digits= max(3, getOption("digits")-4),...) {
  if ('warnings' %in% names(x)) warns <- c(x$warnings, sep="\n")

  numsubprob <- length(levels(attr(x$thematch, "subproblem")))
  numsubprobfail <- if (is.null(x$matching.failed)) 0 else nrow(x$matching.failed)
  numobs <- length(x$thematch)
  numobsfail <- sum(x$matching.failed)

  if (numsubprob == numsubprobfail)
    {
      do.call(cat, warns)
      return(invisible(x))
    }

  if (numsubprobfail > 0)  {
    cat(paste("Matches made in", numsubprob - numsubprobfail,"of", numsubprob,
              "subgroups, accounting for", numobs - numobsfail, "of", numobs, "total observations.\n\n\n"))
  }

  attr(x$matched.set.structures, "comparable.num.matched.pairs") <- NULL
  cat("Structure of matched sets:\n")
  print(x$matched.set.structures)
  cat("Effective Sample Size: ", signif(x$effective.sample.size, digits), "\n")
  cat("(equivalent number of matched pairs).\n\n")
  if (!is.null(x$total.distance)) {
    cat("sum(matched.distances)=",
        signif(x$total.distance, digits),"\n",sep="")
    cat("(within",
        signif(x$total.tolerances,digits),
        "of optimum).\n")
    cat("Percentiles of matched distances:\n")
    print(signif(x$matched.dist.quantiles, digits))
  }

  if ('balance' %in% names(x))
    {
    cat("Balance test overall result:\n")
    b.o <- x$balance$overall
    row.names(b.o) <-  " "
    print(b.o, digits=3)
    }

  if ('warnings' %in% names(x))
    {
      cat("\n")
      do.call(cat, warns)
    }
  invisible(x)
}
