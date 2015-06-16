#' Extract scores (propensity, prognostic,...) from a fitted model
#'
#' This is a wrapper for \code{predict}, adapted for use in matching.  Given a
#' fitted model but no explicit \code{newdata} to \sQuote{predict} from, it
#' constructs its own \code{newdata} in a manner that's generally better suited
#' for matching.
#'
#' In contrast to \code{predict}, if \code{scores} isn't given an explicit
#' \code{newdata} argument then it attempts to reconstruct one from the context
#' in which it is called, rather than from its first argument.  For example, if
#' it's called within the \code{formula} argument of a call to \code{glm}, its
#' \code{newdata} is the same data frame that \code{glm} evaluates that formula
#' in, as opposed to the model frame associated with \code{object}.  See
#' Examples.
#'
#' The handling of missing independent variables also differs from that of
#' \code{predict} in two ways. First, if the data used to generate \code{object}
#' has \code{NA} values, they're mean-imputed using
#' \code{\link{fill.NAs}}. Secondly, if \code{newdata} (either the explicit
#' argument, or the implicit data generated from \code{object}) has \code{NA}
#' values, they're likewise mean-imputed using \code{\link{fill.NAs}}.  Also,
#' missingness flags are added to the formula of \code{object}, which is then
#' re-fit, using \code{\link{fill.NAs}}, prior to calling \code{predict}.
#'
#' If \code{newdata} is specified and contains no missing data, \code{scores}
#' returns the same value as \code{predict}.
#'
#' @param object fitted model object determining scores to be generated.
#' @param newdata (optional) data frame containing variables with which scores
#'   are produced.
#' @param ... additional arguments passed to \code{predict}.
#' @return See individual \code{predict} functions.
#' @author Josh Errickson
#' @seealso \code{\link{predict}}
#' @export
#' @examples
#' data(nuclearplants)
#' pg <- lm(cost~., data=nuclearplants, subset=(pr==0))
#' # The following two lines produce identical results.
#' ps1 <- glm(pr~cap+date+t1+bw+predict(pg, newdata=nuclearplants),
#'            data=nuclearplants)
#' ps2 <- glm(pr~cap+date+t1+bw+scores(pg), data=nuclearplants)
scores <- function(object, newdata=NULL, ...) {

  # if object does not exist then print helpful error msg
  object_str <- deparse(substitute(object))
  data_str <- deparse(substitute(newdata))
  tryCatch(object, error = function(e) {
    stop(missing_x_msg(object_str, data_str, ...))})

  UseMethod("scores")
}

scores.default <- function(object, newdata=NULL, ...) {
  # First, update object to use missingness
  olddata <- tryCatch(model.frame(object, na.action=na.pass, "weights"),
                      error = function(e) {
    warning(paste("Error gathering complete data.",
                  "If the data has missing cases, imputation will not be performed.",
                  "Either be explicit in including `data` arguments to objects, or perform",
                  "imputation beforehand."))
    model.frame(object, "weights")
  })
  wts <- olddata$"(weights)"

  # If the formula is something like `y ~ . - x`, x is included in olddata.
  # This simplifies the formula and drops it
  lhs <- deparse(terms(formula(object), simplify=TRUE)[[2]])
  rhs <- gsub("`", "", attr(terms(formula(object), simplify=TRUE), "term.labels"))
  olddata <- fill.NAs(olddata[, names(olddata) %in% c(rhs, lhs)])
  names(olddata) <- gsub("`", "", names(olddata))

  # rebuild the formula to handle expansion of factors and missing indicators
  vars <- paste0("`",names(olddata)[-1], "`")
  vars <- vars[!grepl(names(olddata)[1], vars, fixed=TRUE)]
  newform <- reformulate(vars, names(olddata[1]))
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
    newdata2 <- model.frame(formula(object), data=parent.frame(), na.action=na.pass)
  } else {
    newdata2 <- model.frame(formula(object), data=newdata, na.action=na.pass)
  }
  newdata2 <- fill.NAs(newdata2)
  names(newdata2) <- gsub("`", "", names(newdata2))

  # If we were given `newdata`, it may contain things we didn't capture yet
  # (Specifically, if fill.NAs was called on newdata, it'll contain some
  # xxx.NA columns)
  if (is.null(newdata)) {
    rhs <- gsub("`", "", attr(terms(newobject), "term.labels"))
    othervars <- rhs[!(rhs %in% names(newdata2))]
    tosearch <- if (length(othervars) > 0) {
      cbind(newdata2, model.frame(reformulate(othervars,), parent.frame()))
    } else {
      newdata2
    }
    newdata2 <- model.frame(formula(newobject),
                            data=tosearch,
                            na.action=na.pass)
  } else {
    newdata2 <- model.frame(formula(newobject), data=cbind(newdata2, newdata),
                            na.action=na.pass)
  }

  eval(predict(newobject, newdata=newdata2, ...))
}

#' Due to the nature of \code{bigglm} objects, the first stage imputation (if the data used
#' to generate \code{object} has \code{NA} values) does not occur.
scores.bigglm <- function(object, newdata=NULL, ...) {

  if (nrow(model.frame(object)) != object$n) {
    warning(paste("Model fit on data with missing values. Updating of model",
                  "with imputed data not supported for class bigglm.",
                  "Please perform imputation prior to calling scores."))
  }

  ## olddata <- model.frame(object, na.action=na.pass)
  ## wts <- weights(object)

  ## # If the formula is something like `y ~ . - x`, x is included in olddata.
  ## # This simplifies the formula and drops it
  ## lhs <- deparse(terms(formula(object), simplify=TRUE)[[2]])
  ## rhs <- attr(terms(formula(object), simplify=TRUE), "term.labels")
  ## olddata <- fill.NAs(olddata[, names(olddata) %in% c(rhs, lhs)])
  ## names(olddata) <- gsub("`", "", names(olddata))

  ## # rebuild the formula to handle expansion of factors and missing indicators
  ## vars <- paste0("`",names(olddata)[-1], "`")
  ## vars <- vars[!grepl(names(olddata)[1], vars, fixed=TRUE)]
  ## newform <- reformulate(vars, names(olddata[1]))

  ## # Don't need subset anymore as the model.frame only pulled out that subset
  ## if (!is.null(wts)) {
  ##   # For some reason, if wts is null, including weights=weights throws an error.
  ##   # So this code is a bit duplicative.
  ##   olddata$weights <- wts
  ##   newobject <- update(object, formula.= newform, data=olddata, weights=weights,
  ##                       subset=NULL)
  ## } else {
  ##   newobject <- update(object, formula.= newform, data=olddata, subset=NULL)
  ## }

  # Now, let's get newdata if its missing
  if (is.null(newdata)) {
    newdata2 <- model.frame(formula(object), data=parent.frame(), na.action=na.pass)
  } else {
    newdata2 <- model.frame(formula(object), data=newdata, na.action=na.pass)
  }
  newdata2 <- fill.NAs(newdata2)
  names(newdata2) <- gsub("`", "", names(newdata2))

  # If we were given `newdata`, it may contain things we didn't capture yet
  # (Specifically, if fill.NAs was called on newdata, it'll contain some
  # xxx.NA columns)
  if (is.null(newdata)) {
    rhs <- gsub("`", "", attr(terms(object), "term.labels"))
    othervars <- rhs[!(rhs %in% names(newdata2))]
    tosearch <- if (length(othervars) > 0) {
      cbind(newdata2, model.frame(reformulate(othervars,), parent.frame()))
    } else {
      newdata2
    }
    newdata2 <- model.frame(formula(object),
                            data=tosearch,
                            na.action=na.pass)
  } else {
    newdata2 <- model.frame(formula(object), data=cbind(newdata2, newdata),
                            na.action=na.pass)
  }

  eval(predict(object, newdata=newdata2, ...))

}
