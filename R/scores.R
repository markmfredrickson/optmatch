#' Extract scores (propensity, prognostic,...) from a fitted model
#'
#' This is a wrapper for \code{predict}, adapted for use in matching.
#' Given a fitted model but no explicit \code{newdata} to \sQuote{predict}
#' from, it constructs its own \code{newdata}
#' in a manner that's generally better suited for matching.
#'
#' In contrast to \code{predict}, if \code{scores} isn't given an
#' explicit \code{newdata} argument then it attempts to reconstruct
#' one from the context in which it is called, rather than from its
#' first argument.  For example, if it's called within the
#' \code{formula} argument of a call to \code{glm}, its \code{newdata}
#' is the same data frame that \code{glm} evaluates that formula in,
#' as opposed to the model frame associated with \code{object}.
#' See Examples.
#'
#' The handling of missing independent variables also differs from
#' that of \code{predict}. If \code{newdata} (either the explicit
#' argument, or the implicit data generated from \code{object}) has
#' \code{NA} values, they're mean-imputed using
#' \code{\link{fill.NAs}}.  Also, missingness flags are added to the
#' formula of \code{object}, which is then re-fit, using \code{\link{fill.NAs}},
#' prior to calling \code{predict}.
#'
#' The mechanics of this re-fitting make it somewhat fragile, particularly for models
#' involving weights, offsets, or sample exclusions conveyed via a \code{subset} argument
#' to the model fitter. In such circumstances it's best to address missing observations before
#' passing \code{object} to \code{scores}, ensuring that \code{na.action(object)} is \code{NULL}.
#'
#' If \code{newdata} is specified and contains no missing data, \code{scores} returns the same value as
#' \code{predict}.
#'
#' @param object fitted model object determining scores to be generated.
#' @param newdata (optional) data frame containing variables with which scores are produced.
#' @param ... additional arguments passed to \code{predict}.
#' @return See individual \code{predict} functions.
#' @author Josh Errickson
#' @seealso \code{\link{predict}}
#' @export
#' @examples
#' data(nuclearplants)
#' pg <- lm(cost~., data=nuclearplants, subset=(pr==0))
#' # The following two lines produce identical results.
#' ps1 <- glm(pr~cap+date+t1+bw+predict(pg, newdata=nuclearplants), data=nuclearplants)
#' ps2 <- glm(pr~cap+date+t1+bw+scores(pg), data=nuclearplants)
scores <- function(object, newdata=NULL,...)
{
  # First, update object to use missingness
  olddata <- model.frame(object, na.action=na.pass, "weights")
  wts <- olddata$"(weights)"
  # If the formula is something like `y ~ . - x`, x is included in olddata.
  # This simplifies the formula and drops it
  olddata <- fill.NAs(model.frame(formula(terms(formula(object), simplify=TRUE)),
                                  na.action=na.pass, data=olddata))

  # rebuild the formula to handle expansion of factors and missing indicators
  newform <- reformulate(names(olddata)[-1], names(olddata)[1])
  # remove weights if it's hanging around.
  newform <- update(newform, . ~ . - `(weights)`)

  # Don't need subset anymore as the model.frame only pulled out that subset
  if (!is.null(wts)) {
    # For some reason, if wts is null, including weights=weights throws an error.
    # So this code is a bit duplicative.
    olddata$weights <- wts
    newobject <- update(object, formula.= newform, data=olddata, weights=weights,
                        subset=NULL)
  } else {
    newobject <- update(object, formula.= newform, data=olddata, subset=NULL)
  }

  # Now, let's get newdata if its missing
  if (is.null(newdata)) {
    newdata <- fill.NAs(get_all_vars(object, parent.frame()))
  } else {
    newdata <- fill.NAs(newdata)
  }

  eval(predict(newobject, newdata=newdata, ...))
}
