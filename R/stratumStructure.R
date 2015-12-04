##' Tabulate treatment:control ratios occurring in matched sets, and
##' the frequency of their occurrence.
##'
##' @title Return structure of matched sets
##'
##' @param stratum Matched strata, as returned by
##'   \code{\link{fullmatch}} or \code{\link{pairmatch}}
##' @param trtgrp Dummy variable for treatment group membership.  (Not
##'   required if \code{stratum} is an optmatch object, as returned by
##'   \code{\link{fullmatch}} or \code{\link{pairmatch}}.)
##' @param min.controls For display, the number of treatment group
##'   members per stratum will be truncated at the reciprocal of
##'   \code{min.controls}.
##' @param max.controls For display, the number of control group
##'   members will be truncated at \code{max.controls}.
##' @param x stratumStructure object to be printed.
##' @param ... Additional arguments to \code{print}.
##' @return A table showing frequency of occurrence of those
##'   treatment:control ratios that occur.
##'
##'   The \sQuote{effective sample size} of the stratification, in
##'   matched pairs.  Given as an attribute of the table, named
##'   \sQuote{\code{comparable.num.matched.pairs}}; see Note.
##'
##' @note For comparing treatment and control groups both of size 10,
##'   say, a stratification consisting of two strata, one with 9
##'   treatments and 1 control, has a smaller \sQuote{effective sample
##'   size}, intuitively, than a stratification into 10 matched pairs,
##'   despite the fact that both contain 20 subjects in
##'   total. \code{stratumStructure} first summarizes this aspect of
##'   the structure of the stratification it is given, then goes on to
##'   identify one number as the stratification's effective sample
##'   size.  The \sQuote{\code{comparable.num.matched.pairs}}
##'   attribute returned by \code{stratumStructure} is the sum of
##'   harmonic means of the sizes of the treatment and control
##'   subgroups of each stratum, a general way of calibrating such
##'   differences as well as differences in the number of subjects
##'   contained in a stratification.  For example, by this metric the
##'   9:1, 1:9 stratification is comparable to 3.6 matched pairs.
##'
##'   Why should effective sample size be calculated this way?  The
##'   phrase \sQuote{effective sample size} suggests the observations
##'   are taken to be similar in information content.  Modeling them
##'   as random variables, this suggests that they be assumed to have
##'   the same variance, \eqn{\sigma}{sigma}, conditional on what
##'   stratum they reside in.  If that is the case, and if also
##'   treatment and control observations differ in expectation by a
##'   constant that is the same for each stratum, then it can be shown
##'   that the optimum weights with which to combine treatment-control
##'   contrasts across strata, \eqn{s}{s}, are proportional to the
##'   stratum-wise harmonic means of treatment and control counts,
##'   \eqn{h_s = [(n_{ts}^{-1} + n_{cs}^{-1})/2]^{-1}}{h[s] =
##'   1/(0.5/n.t[s] + 0.5/n.c[s])} (Kalton, 1968).  The thus-weighted
##'   average of contrasts then has variance \eqn{2\sigma/\sum_s
##'   h_s}{2*sigma/sum(h)}. This motivates the use of \eqn{\sum_s
##'   h_s}{sum(h)} as a measure of effective sample size.  Since for a
##'   matched pair \eqn{s}{s}, \eqn{h_s=1}{h[s]=1}, \eqn{\sum_s
##'   h_s}{sum(h)} can be thought of as the number of matched pairs
##'   needed to attain comparable precision.  (Alternately, the
##'   stratification might be taken into account when comparing
##'   treatment and control groups using fixed effects in an ordinary
##'   least-squares regression, as in Hansen (2004). This leads to the
##'   same result.  A still different formulation, in which outcomes
##'   are not modeled as random variables but assignment to treatment
##'   or control is, again suggests the same weighting across strata,
##'   and a measure of precision featuring \eqn{\sum_s h_s}{sum(h)} in
##'   a similar role; see Hansen and Bowers (2008).
##' @author Ben B. Hansen
##' @references Kalton, G. (1968), \sQuote{Standardization: A
##'   technique to control for extraneous variables}, \emph{Applied
##'   Statistics}, \bold{17}, 118--136.
##'
##'   Hansen, B.B. (2004), \sQuote{Full Matching in an Observational
##'   Study of Coaching for the SAT}, \emph{Journal of the American
##'   Statistical Association}, \bold{99}, 609--618.
##'
##'   Hansen B.B. and Bowers, J. (2008), \sQuote{Covariate balance in
##'   simple, stratified and clustered comparative studies},
##'   \emph{Statistical Science}, \bold{23}, to appear.
##' @seealso \code{\link{matched}}, \code{\link{fullmatch}}
##' @examples
##' data(plantdist)
##' plantsfm <- fullmatch(plantdist) # A full match with unrestricted
##'                                  # treatment-control balance
##' plantsfm1 <- fullmatch(plantdist,min.controls=2, max.controls=3)
##' stratumStructure(plantsfm)
##' stratumStructure(plantsfm1)
##' stratumStructure(plantsfm, max.controls=4)
##'
##' @export
##' @rdname stratumStructure
stratumStructure <- function(stratum, trtgrp=NULL, min.controls=0,max.controls=Inf) UseMethod("stratumStructure")

##' @export
##' @rdname stratumStructure
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

##' @export
##' @rdname stratumStructure
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

##' @export
##' @rdname stratumStructure
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
#' @rdname effectiveSampleSize
effectiveSampleSize <- function(x, z = NULL) UseMethod("effectiveSampleSize")

#' @export
#' @rdname effectiveSampleSize
effectiveSampleSize.factor <- function(x, z = NULL) {
  if (is.null(z)) {
    z <- getZfromMatch(x)
  }

  effectiveSampleSize.default(x,z)
}

#' @export
#' @rdname effectiveSampleSize
effectiveSampleSize.default <- function(x, z = NULL) {

  if (missing(z) || is.null(z)) stop("default effectiveSampleSize method requires a treatment indicator, z")

  wasMatched <- !is.na(x)

  if (!any(wasMatched)) { return(0) }

  totals <- table(x[wasMatched], as.logical(z)[wasMatched])
  effectiveSampleSize.table(totals, z)
}

#' @export
#' @rdname effectiveSampleSize
effectiveSampleSize.table <- function(x, z = NULL) {

  stopifnot(length(dim(x))==2, ncol(x)==1 || ncol(x)==2,
            all(colnames(x) %in% c("FALSE","TRUE")))

  if (ncol(x)==1 | nrow(x)==0) return(0)

  sum(2/(1/x[,"FALSE"] + 1/x[,"TRUE"]))

}
