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

  form <- formula(terms(formula(object), simplify=TRUE))
  nm <- all.vars(form)

  tmpform <- as.formula(paste("~", paste(nm, collapse="+")))

  fulldata <- eval(object$call$data, envir=attr(object$terms,".Environment"))

  mf <- model.frame(tmpform, data=fulldata,
                    subset = eval(object$call$subset, envir=fulldata),
                    na.action=na.pass)

  if (!all(complete.cases(mf))) {
    fill <- fill.NAs(mf)

    extras <- names(fill)[(ncol(mf)+1):ncol(fill)]

    newform <- formula(paste(row.names(attr(terms(form), "factors"))[1],
                             paste(attr(terms(form), "term.labels"),
                                   collapse="+"),
                             sep="~"))
    # Subset no longer needed before `model.frame` pulls out only the necessary
    # subset already.
    object <- update(object, newform, data=fill, subset = NULL)
    warning("Missing data found and imputed.")
  }

  if (is.null(newdata)) {
    newdata <- get_all_vars(object, data=parent.frame())

    fnd <- fill.NAs(newdata)

    p <- try(predict(object, newdata=fnd,...),
             silent=TRUE)
    if (!is(p, "try-error")) {
      return(p)
    }

    warning("Imputation into fulldata failed!")
  }

  predict(object, newdata=newdata, ...)
}
