##' .. content for \description{} (no empty lines) ..
##'  Given a vector or matrix of outcome measures, calculate the
##' treatment minus control difference of means within each matched set
##' and then average them, using weights to reflect differences in sizes or
##' compositions of matched sets, to produced matched estimates of the treatment
##' effect.  Also furnish standard errors and estimated covariances for the effects.
##' .. content for \details{} ..
##' The standard error estimate is essentially that of Abadie and Imbens
##' (2012, equation 5), generalizing that expression to accommodate matched
##' sets of varying structure and flexible weighting schemes.  It's determined with conditioning on the values of the treatment variable.  The same equation is adapted to calculation of covariances.
##'
##' In the presence of missing values, each matched set reports a treatment minus control mean difference for an outcome only if there is at least one non-missing treatment observationan and at least one non-missing control observation for that outcome.  Missing outcomes don't affect the weight associated with a matched set.  Effects are determined by a weighted average of non-missing matched differences.  Variances are also calculated from non-missing matched differences, with associated weights.  If multiple outcomes are presented, a correlation matrix is computed using the collection of matched sets for which a matched difference is available on each of the outcomes, and the covariance is determined from it and from the variances that were calculated for each outcome separately, with its distinctive missingness pattern.  An exception to this occurs if there are not two or more matched sets for which a matched difference is available on each of the outcomes; then correlations are calculated using all matched differences available for each of the two variables in question.  (This procedure risks producing a covariance that isn't nonnegative definite, which is why it's used sparingly.)
##' 
##' @title Means of matched differences, with standard errors
##' @param x a numeric vector or matrix of outcomes
##' @param matches optmatch object recording a match
##' @param weighting.scheme for averaging contributions from matched sets of differing sizes
##' @param keep.differences retain the matrix of matched differences?
##' @param ... additional arguments, not currently used.
##' @return
##' A list with entries \sQuote{effects}, a vector of length \code{ncol(x)} or, if \code{x} is a vector, 1; \sQuote{covariance}, a corresponding square matrix; and,
##' if \code{keep.differences=TRUE}, a matrix \sQuote{differences} of matched differences.
##' @references A. Abadie and G.W. Imbens (2012), A martingale representation for matching estimators, J. Amer. Statist. Assoc. 107/496 833-843. 
##' @author Ben B Hansen
matched_effects <- function(x, matches, weighting.scheme=c('ETT','harmonic')[1], keep.differences=FALSE,...)
{
stopifnot(inherits(matches, "optmatch"),
          'contrast.group' %in% names(attributes(matches)),
          is.character(weighting.scheme))
UseMethod("matched_effects")
}
matched_effects.default <- function(x, matches, weighting.scheme,keep.differences,...)
  {
x <- as.matrix(x)
colnames(x) <- if (ncol(x)>1) paste('y',1L:ncol(x), sep="") else 'y'
matched_effects_matrix(x, matches, weighting.scheme, keep.differences,...)
}
matched_effects.matrix <- function(x, matches, weighting.scheme, keep.differences,...)
  {
    stopifnot(is.numeric(x), nrow(x)==length(matches))
    Zz <- attr(matches, "contrast.group")
    wt.fct <- switch(weighting.scheme[1],
                     e=, et=, ett=, E=, ET=, ETT=function(z)sum(as.logical(z)),
                     h=,H=,ha=,Ha=,Harmonic=,
                     harmonic=function(z) (2*(length(z)-1)/length(z)))
    weights.unscaled <- tapply(Zz, matches, wt.fct) # speed me up!
    ms.diffs <- tapply(# please speed me up by using data.table instead!
                       cbind(Zz,x),
                       matches,
                       function(mat) { 
                         colMeans(mat[as.logical(mat[,1]),-1], na.rm=TRUE) -
                           colMeans(mat[!mat[,1],-1], na.rm=TRUE)
                       }, simplify=FALSE
                       )
    diffmat <- matrix(unlist(ms.diffs), length(ms.diffs), ncol(x), byrow=TRUE,
                      dimnames=list(names(ms.diffs), colnames(x)))
    weights.totals <- colSums( ifelse(is.na(diffmat), 0, weights.unscaled) )

    ans <- list()    
    ans$`effects` <- colSums(weights.unscaled*diffmat, na.rm=TRUE)/weights.totals 
    ans$`differences` <- if (keep.differences) diffmat else NULL
    diffmat <- diffmat -
      matrix(ans$`effects`, nrow(diffmat), ncol(diffmat), byrow=TRUE)
    
    vars <- colSums(weights.unscaled^2*diffmat^2, na.rm=TRUE)/
      weights.totals^2
    if (ncol(x)==1) { dim(vars) <- c(1,1) ; ans$`covariance` <- vars } else
    { # have to build a correlation matrix.  If there are enough complete cases,
      # base the correlation matrix on them --- ensures nonnegative-definiteness
      if (sum(acs <- complete.cases(diffmat))>1)
        {
          devprods <- apply(weights.unscaled[acs]*diffmat,1, function(x) x%o%x)
          csums <- rowSums(devprods)
          dim(csums) <- rep(ncol(diffmat), 2)
          corrmat <- cov2cor(csums)
        } else
      {
      corrmat <- diag(ncol(x)) #Separating b/c each entry can have different NA pattern.
      for (i in 2L:ncol(x)) for (j in 1L:(i-1))
        {
          acs <- complete.cases(diffmat[,c(i,j)])
          cij <- colSums(weights.unscaled[acs]^2*diffmat[acs,i]*
                         diffmat[acs,j])
          vs <- colSums(weights.unscaled[acs]^2*diffmat[acs,c(i,j)]^2)
          corrmat[i,j] <- corrmat[j,i] <- cij/prod(sqrt(vs))
        }
    }
      SEs <- sqrt(vars)
      ans$`covariance` <- corrmat * outer(SEs, SEs)
    }
    ans <- ans[c('effects','covariance','differences')]
    if (is.null(ans$'differences')) ans <- ans[-3L]
    class(ans) <- c("optmatch_matched_effects", class(ans))
    ans
  }

print.optmatch_matched_effects <- function(x,digits=getOption("digits"),...)
  {
cat("Effect estimates:\n")    
signif(cbind(`effects`=ans$effects,`se`=sqrt(diag(ans$covariance))), digits)
  }
