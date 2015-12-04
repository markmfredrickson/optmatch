#' @export
stratumStructure <- function(stratum, trtgrp=NULL, min.controls=0,max.controls=Inf) UseMethod("stratumStructure")

#' @export
stratumStructure.optmatch <- function(stratum,trtgrp, min.controls=0,max.controls=Inf) {
  trtgrp.arg.provided <- !missing(trtgrp) && !is.null(trtgrp)
  ZZ <- try(getZfromMatch(stratum), silent=TRUE)
  if (inherits(ZZ, "try-error") & !trtgrp.arg.provided)
    stop("stratum is of class optmatch but it has lost its contrast.group attribute; must specify trtgrp")

  if (inherits(ZZ, "try-error") & trtgrp.arg.provided)
    return(stratumStructure.default(stratum,trtgrp=trtgrp,min.controls=min.controls,max.controls=max.controls))

  if (trtgrp.arg.provided) # by implication, ZZ is not an error
  {
    warning("ignoring trtgrp argument to stratumStructure")
  }
  stratumStructure.default(stratum,trtgrp=ZZ,min.controls=min.controls,max.controls=max.controls)
}

#' @export
stratumStructure.default <- function(stratum,trtgrp,min.controls=0,max.controls=Inf) {
  if (is.null(trtgrp))
    stop("Unless stratum is of class \'optmatch\', stratumStructure() requires a trtgrp= argument.")

  if (!any(trtgrp<=0) | !any(trtgrp>0))
    warning("No variation in (trtgrp>0); was this intended?")

  stopifnot(is.numeric(min.controls), is.numeric(max.controls))

  if (length(min.controls)>1) warning("Only first element of min.controls will be used.")
  if (length(max.controls)>1) warning("Only first element of max.controls will be used.")


  stratum <- as.integer(as.factor(stratum))
  if (any(is.na(stratum)))
    stratum[is.na(stratum)] <- max(0, stratum, na.rm=TRUE) + 1:sum(is.na(stratum))

  ttab <- table(stratum,as.logical(trtgrp))
  comp.num.matched.pairs <- effectiveSampleSize(ttab)

  max.tx <- round(1/min.controls[1])
  max.controls <- round(max.controls[1])
  txsz <- pmin(ttab[,2], max.tx)
  ctlsz <- pmin(ttab[,1], max.controls)
  ans <- table(paste(txsz, ctlsz, sep=":"),
               dnn="stratum treatment:control ratios")

  tnn <- as.numeric(unlist(strsplit(names(ans), ":", fixed=FALSE)))
  i.ctl <- 2*(1:length(ans))
  i.tx <- 2*(1:length(ans))-1
  txnms <- as.character(tnn[i.tx])
  txnms[tnn[i.tx]==max.tx] <-
    paste(max.tx,"+", sep="")
  ctlnms <- as.character(tnn[i.ctl])
  ctlnms[tnn[i.ctl]==max.controls] <- paste(max.controls,"+",sep="")
  names(ans) <- paste(txnms, ctlnms, sep=":")

  onez <- tnn[i.tx]==1 & tnn[i.ctl]==0
  if (any(onez)) {
    tnn[i.tx][onez] <- Inf
    tnn[i.ctl][onez] <- 1
  }
  ans <- ans[order(-tnn[i.tx],tnn[i.ctl])]

  attr(ans, "comparable.num.matched.pairs") <- comp.num.matched.pairs
  class(ans) <- append(class(ans), "stratumStructure")
  ans
}

#' @export
print.stratumStructure <- function(x, ...) {
  attr(x, "comparable.num.matched.pairs") <- NULL
  print.table(x, ...)
}

getZfromMatch <- function(m) {
  if (!is.null(attr(m, "contrast.group"))) {
    # NB: originally, the next line called toZ(attr(...))
    # but this caused problems when there NAs induced into the match
    # by, for example, needing to make the match as long as a data.frame
    # that had missingness that was kicked out by glm() or other row-wise
    # deleting functions. For now, we ignore that problem in this function.
    return(attr(m, "contrast.group"))
  }

  stop("Unable to find 'contrast.group' attribute (treatment indicator)")
}

#' Compute the effective sample size of a match.
#'
#' The effective sample size is the sum of the harmonic means of the number
#' units in treatment and control for each matched group. For k matched pairs,
#' the effective sample size is k. As matched groups become more unbalanced, the
#' effective sample size decreases.
#'
#' @param x An \code{optmatch} object, the result of
#' \code{\link{fullmatch}} or \code{\link{pairmatch}}.
#' @param z A treatment indicator, a vector the same length as \code{match}.
#' This is only required if the \code{match} object does not contain the
#' contrast.group' attribute.
#' @seealso \code{\link{summary.optmatch}}, \code{\link{stratumStructure}}
#' @return The equivalent number of pairs in this match.
#' @export
effectiveSampleSize <- function(x, z = NULL) UseMethod("effectiveSampleSize")

#' @export
effectiveSampleSize.factor <- function(x, z = NULL) {
  if (is.null(z)) {
    z <- getZfromMatch(x)
  }

  effectiveSampleSize.default(x,z)
}

#' @export
effectiveSampleSize.default <- function(x, z = NULL) {

  if (missing(z) || is.null(z)) stop("default effectiveSampleSize method requires a treatment indicator, z")

  wasMatched <- !is.na(x)

  if (!any(wasMatched)) { return(0) }

  totals <- table(x[wasMatched], as.logical(z)[wasMatched])
  effectiveSampleSize.table(totals, z)
}

#' @export
effectiveSampleSize.table <- function(x, z = NULL) {

  stopifnot(length(dim(x))==2, ncol(x)==1 || ncol(x)==2,
            all(colnames(x) %in% c("FALSE","TRUE")))

  if (ncol(x)==1 | nrow(x)==0) return(0)

  sum(2/(1/x[,"FALSE"] + 1/x[,"TRUE"]))

}
