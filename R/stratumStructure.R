stratumStructure <- function(stratum,trtgrp=NULL,min.controls=0,
                             max.controls=Inf)
{
if (!inherits(stratum,"optmatch") & is.null(trtgrp))
  stop("stratum not of class \'optmatch\'; trtgrp must be specified")
if (!inherits(stratum,"optmatch"))
  warning("stratum not of class optmatch; was this intended?")
if (inherits(stratum, "optmatch") & is.null(attr(stratum, "contrast.group")) & is.null(trtgrp))
  stop("Argument 1 is of class optmatch but it has lost its contrast.group attribute; must specify trtgrp")
if (inherits(stratum, "optmatch") & !is.null(attr(stratum, "contrast.group")) & !is.null(trtgrp))
  warning("ignoring second argument to stratumStructure")

tgp <- getZfromMatch(stratum)

if (!any(tgp<=0) | !any(tgp>0))
   warning("No variation in (trtgrp>0); was this intended?")

stopifnot(is.numeric(min.controls), is.numeric(max.controls))

if (length(min.controls)>1) warning("Only first element of min.controls will be used.")
if (length(max.controls)>1) warning("Only first element of max.controls will be used.")

comp.num.matched.pairs <- effectiveSampleSize(stratum)

stratum <- as.integer(as.factor(stratum))
if (any(is.na(stratum)))
  stratum[is.na(stratum)] <- max(0, stratum, na.rm=TRUE) + 1:sum(is.na(stratum))

ttab <- table(stratum,as.logical(tgp))
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
if (any(onez))
  {
tnn[i.tx][onez] <- Inf
tnn[i.ctl][onez] <- 1
}
ans <- ans[order(-tnn[i.tx],tnn[i.ctl])]

attr(ans, "comparable.num.matched.pairs") <- comp.num.matched.pairs
ans
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

#' (Internal) Compute the effective sample size of a match.
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
#' @return The equivalent number of pairs in this match.
effectiveSampleSize <- function(x, z = NULL) {
  if (is.null(z)) {
    z <- getZfromMatch(x)
  }

  wasMatched <- !is.na(x)
  totals <- table(x[wasMatched], z[wasMatched])

  sum(2/(1/totals[,1] + 1/totals[,2]))
}
