##' Wrapper for \code{predict} to cleanly look for new data to predict on.
##'
##' When called without a \code{newdata} argument, it will attempt to determine the correct
##' new data to predict on; e.g. in a \code{lm} or \code{glm} model, will use the data in that
##' model.
##'
##' Specifying \code{newdata} is identical to calling \code{predict}.
##'
##' It is not recommended to use \code{attach} and \code{detach} when using scores. The preferred methods
##' are using \code{with} or using the \code{newdata} parameter.
##'
##' @param object a model object from which prediction is desired.
##' @param newdata optionally, specifies a data frame in which to look for variables to predict
##' with. When omitted, attempts to intelligently use the correct data frame as opposed to \code{predict}
##' using the data originally used to credit \code{object}.
##' @param ... additional arguments passed to \code{predict}.
##' @return See individual \code{predict} functions.
##' @author Josh Errickson
##' @seealso \code{\link{predict}}
##' @examples
##' data(nuclearplants)
##' pg <- lm(cost~., data=nuclearplants, subset=(pr==0))
##' # The following two lines produce identical results.
##' ps1 <- glm(pr~cap+date+t1+bw+predict(pg, newdata=nuclearplants), data=nuclearplants)
##' ps2 <- glm(pr~cap+date+t1+bw+scores(pg), data=nuclearplants)
scores <- function(object, newdata=NULL,...)
{
  # If user didn't give newdata, extract from model call
  if (is.null(newdata)) {
    newdata <- model.frame(terms(object), data=parent.frame())
  }
  predict(object, newdata=newdata, na.action=na.omit,...)
}
