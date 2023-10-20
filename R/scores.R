#' Extract scores (propensity, prognostic,...) from a fitted model
#'
#' This is a wrapper for \code{predict}, adapted for use in matching.  Given a
#' fitted model but no explicit \code{newdata} to \sQuote{predict} from, it
#' constructs its own \code{newdata} in a manner that's generally better suited
#' for matching.
#'
#' Like \code{predict}, its default predictions from a \code{glm} are on
#' the scale of the linear predictor, not the scale of the response; see
#' Rosenbaum \ Rubin (1985).  (This default can
#' be overridden by specifying \code{type="response"}.)
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
#' @references P.~R. Rosenbaum and D.~B. Rubin (1985), \sQuote{Constructing a
#'   control group using multivariate matched sampling methods that incorporate
#'   the propensity score}, \emph{The American Statistician}, \bold{39} 33--38.
#' @export
#' @examples
#' data(nuclearplants)
#' pg <- lm(cost~., data=nuclearplants, subset=(pr==0))
#' # The following two lines produce identical results.
#' ps1 <- glm(pr~cap+date+t1+bw+predict(pg, newdata=nuclearplants),
#'            data=nuclearplants)
#' ps2 <- glm(pr~cap+date+t1+bw+scores(pg), data=nuclearplants)
scores <- function(object, newdata=NULL, ...) {

  if (is(object, "bigglm")) {
    # bigglm has two problems:
    # 1) update.bigglm only allows adding new data, not replacing the old.
    # Solution: modify object$call$data and re-eval.
    # 2) model.frame(object, na.action=na.pass) fails without a data
    #    argument.
    # Don't have a good solution - we could require that `newdata` be given
    # but that's weird, because newdata and object$data don't have to be
    # the same.
    #
    # For now, bypass imputation into data if object is bigglm

    if (nrow(model.frame(object)) != object$n) {
      warning(paste("Imputation and refitting of bigglm objects",
                    "prior to prediction is not supported. Please",
                    "impute before calling scores."))
    }
    # Still using `fill.NAs` just to get expansion of factors the same
    olddata <- fill.NAs(model.matrix(formula(object),
                                     data=model.frame(object)))[,-1,drop=FALSE]
    colnames(olddata) <- gsub("`", "", colnames(olddata))
    olddata <- cbind(model.frame(object)[,1,drop=FALSE], olddata)

    call <- object$call
    call$data <- olddata
    call[[2]] <- formula(olddata)
    newobject <- eval(call)
  } else {
    mf <- tryCatch(model.frame(object, na.action=na.pass),
                   error = function(e) {
      fallback <- model.frame(object)
      warning(paste("Error gathering complete data.",
                    "If the data has missing cases, imputation will not be performed.",
                    "(Sometimes this can be fixed by supplying a `data` argument",
                    "when fitting the model that's to be passed to `scores`. Alternatively,",
                    "just take care of (impute) NAs before you fit that model.)")
              )
      fallback
    })
    wts <- mf$"(weights)"

    olddata <- fill.NAs(as.data.frame(model.matrix(formula(object),
                                                   data=mf)))[,-1,drop=FALSE]
    colnames(olddata) <- gsub("`", "", colnames(olddata))
    olddata <- tryCatch(cbind(model.frame(object, na.action=na.pass)[,1,drop=FALSE], olddata),
                        error = function(e) {
      cbind(model.frame(object)[,1,drop=FALSE], olddata)
    })

    if ( all('weights' != names(object$call)) ) {
      newobject <- update(object, formula.=formula(olddata),
                          data=olddata, subset=NULL)
    } else {
      newobject <- update(object, formula.=formula(olddata),
                          data=olddata, subset=NULL, weights=wts,
                          evaluate=FALSE)
      newobject$weights <- wts
      newobject <- eval(newobject)
      # This bit of silliness is because update evaluates its weight argument in
      # the wrong frame, and it was frustrating to find the correct one. This
      #  workaround replaces the weights with an actual vector in the call.
    }
  }

  # Now, let's get newdata if its missing
  if (is.null(newdata)) {
    newdata2 <- model.frame(formula(object), data=parent.frame(), na.action=na.pass)
  } else {
    newdata2 <- model.frame(formula(object), data=newdata, na.action=na.pass)
  }
  resp <- newdata2[,1,drop=FALSE]
  newdata2 <- fill.NAs(as.data.frame(model.matrix(formula(object),
                                                  data=newdata2)))[,-1,drop=FALSE]
  names(newdata2) <- gsub("`", "", names(newdata2))
  newdata2 <- cbind(resp, newdata2)

  eval(predict(newobject, newdata=newdata2, ...))
}


##' (Internal) Predict for CBPS objects 
##'
##' The CBPS package fits \sQuote{covariate balancing propensity score} for use in propensity score
##' weighting.  In the binary treatment case they can also be used for matching.  This method helps to 
##' implement that process in a manner consistent with use of propensity scores elsewhere in optmatch; see
##' \code{\link{scores}} documentation.
##'
##' @param object A CBPS object
##' @param newdata Unused.
##' @param type Return inverse logits of fitted values (the default) or fitted values themselves 
##' @param ... Unused.
##'
##' @return Inverse logit of the fitted values.
##' @importFrom stats plogis
predict.CBPS <- function(object, newdata=NULL, type=c("link", "response"), ...) {
    type <- match.arg(type)
    stopifnot(type %in% c("link", "response") )
    if (length(unique(object$y))>2) stop("Only binary treatments are supported")

    out <- if (type=="link") stats::plogis(object$fitted.values) else object$fitted.values
  names(out) <- rownames(object$x)
  return(out)
}
