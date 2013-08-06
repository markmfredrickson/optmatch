#' Wrapper for \code{predict} to cleanly look for new data to predict on.
#'
#' When called without a \code{newdata} argument, it will attempt to determine the correct
#' new data to predict on; e.g. in a \code{lm} or \code{glm} model, will use the data in that
#' model.
#'
#' If \code{newdata} (either the explicit argument, or the impicit data generated from
#' \code{object}) has \code{NA} values, imputation will be performed on the missing data via
#' the \code{\link{fill.NAs}} function and \code{object} will be refit using the imputed data
#' frame, before calling \code{predict}
#'
#' If \code{newdata} is specified and contains no missing data, this is identical to calling
#' \code{predict}.
#'
#' If the call to create \code{object} is involved, particularly if it includes optional
#' arguments such as \code{subset} or \code{weights} whose values reference the data, this
#' function may fail or otherwise have undesirable results if the \code{newdata} argument is
#' not given. It is therefore strongly recommended to include the \code{newdata} argument in
#' these sort of situations.
#'
#' @param object a model object from which prediction is desired.
#' @param newdata optionally, specifies a data frame in which to look for variables to predict
#' with. When omitted, attempts to intelligently use the correct data frame as opposed to
#' \code{predict} using the data originally used to create \code{object}.
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
  # If user didn't give newdata, extract from model call
  if (is.null(newdata)) {
    newdata <- get_all_vars(object, data=parent.frame())
  }
  if (!all(complete.cases(newdata))) {
    newdata2 <- eval(fill.NAs(formula(object), data=newdata))
    # so data will be used first from newdat2, then from newdata
    alldata <- cbind(newdata2,newdata)
    newobj <- eval(update(object, formula=formula(newdata2), data=alldata))
    thescores <- predict(newobj, newdata=alldata, ...)
    if (any(is.na(thescores))) warning("Couldn't figure out how get rid of NAs")
    return(thescores)
  }
  eval(predict(object, newdata=newdata,...))
}
