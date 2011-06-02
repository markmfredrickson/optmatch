\name{mdist}
\alias{mdist}
\title{Create matching distances}
\description{
  A generic function, with several supplied methods, for creating
  matrices of distances between observations to be used in the match process.
  Using these matrices, 
  \code{pairmatch()} or \code{fullmatch()} can determine optimal matches.
}
\usage{
  mdist(x, structure.fmla = NULL, ...)
}
\arguments{
  \item{x}{The object to use as the basis for forming the mdist.
  Methods exist for formulas, functions, and generalized linear models.}
  \item{structure.fmla}{A formula denoting the treatment variable on the 
  left hand side and an optional grouping expression on the right hand side. For
  example, \code{z ~ 1} indicates no grouping. \code{z ~ s} subsets the data only
  computing distances within the subsets formed by \code{s}. See method
  notes, below, for
  additional formula options.}
  \item{...}{Additional method arguments. Most methods require a 'data' argument.}
}
\details{
  The \code{mdist} method provides three ways to construct a matching distance
  (ie, a distance matrix or suitably organized list of such matrices):
  guided by a function, by a fitted model, or by a formula.  The class
  of the first argument given to \code{mdist} determines which of these
  methods is invoked.

  The \code{mdist.function} method takes a function of two arguments. When
  called, this function will receive the treatment observations as the first argument
  and the control observations as the second argument. As an example, the following
  computes the raw differences between values of \code{t1} for treatment
  units (here, nuclear plants with \code{pr==1}) and controls (here,
  plants with \code{pr==0}), returning the result as
  a distance matrix:
  \code{sdiffs <- function(treatments, controls) {
      abs(outer(treatments$t1, controls$t1, `-`))
    }
  }

  The \code{mdist.function} method does similar things as the
  earlier optmatch function \code{\link{makedist}}, although the
  interface is a bit different.
  
  The \code{mdist.formula} method computes the squared Mahalanobis
  distance between observations, with the right-hand side of the formula
  determining which variables contribute to the Mahalanobis distance. If
  matching is to be done within strata, the stratification can be
  communicated using either the \code{structure.fmla} argument (e.g. \code{~ grp}) or as part of
  the main formula (e.g. \code{z ~ x1 + x2 | grp}).

  An \code{mdist.glm} method takes an argument of class \code{glm} as
  first argument.  It assumes that this object is a fitted propensity
  model, extracting distances on the linear propensity score (logits of
  the estimated conditional probabilities) and, by default, rescaling the distances
  by the reciprocal of the pooled s.d. of treatment- and control-group propensity scores.
  (The scaling uses \code{mad}, for resistance to outliers, by default; this can be
  changed to the actual s.d., or rescaling can be skipped entirely, by
  setting argument \code{standardization.scale} to \code{sd} or
  \code{NULL}, respectively.)  A \code{mdist.bigglm}
  method works analogously with \code{bigglm} objects, created by
  the \code{bigglm} function from package \sQuote{biglm}, which can
  handle bigger data sets than the ordinary glm function can.  In
  contrast with \code{mdist.glm} it requires additional \code{data} and
  \code{structure.fmla} arguments.  (If you have enough data to have to
  use \code{bigglm}, then you'll probably have to subgroup before
  matching to avoid memory problems. So you'll have to use the
  \code{structure.fmla} argument anyway.)
}
\value{
  Object of class \code{optmatch.dlist}, which is suitable to be given
  as \code{distance} argument to \code{\link{fullmatch}} or
  \code{\link{pairmatch}}. For more information, see \code{\link{pscore.dist}}

}
\author{ Mark M. Fredrickson} %Ben B. Hansen

\references{
  P.~R. Rosenbaum and D.~B. Rubin (1985),
\sQuote{Constructing a control group using multivariate matched sampling
  methods that incorporate the propensity score},
 \emph{The American Statistician}, \bold{39} 33--38.
}

\seealso{ \code{\link{makedist}}, \code{\link{mahal.dist}}, \code{\link{fullmatch}}, \code{\link{pairmatch}},
          \code{\link{pscore.dist}}}
\examples{
data(nuclearplants)
mdist.examples <- list()
### Propensity score distances.
### Recommended approach:
(aGlm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants))
mdist.examples$ps1 <- mdist(aGlm)
### A second approach: first extract propensity scores, then separately
### create a distance from them.  (Useful when importing propensity
### scores from an external program.)
plantsPS <- predict(aGlm)
mdist.examples$ps2 <- mdist(pr~plantsPS, data=nuclearplants)^(1/2)
### Full matching on the propensity score.
fullmatch(mdist.examples$ps1)
fullmatch(mdist.examples$ps2)
### Because mdist.glm uses robust estimates of spread, 
### the results differ in detail -- but they are close enough
### to yield similar optimal matches.
all(fullmatch(mdist.examples$ps1)==fullmatch(mdist.examples$ps2)) # The same

### Mahalanobis distance:
mdist.examples$mh1 <- mdist(pr ~ t1 + t2, data = nuclearplants)

### Absolute differences on a scalar:

sdiffs <- function(treatments, controls) {
  abs(outer(treatments$t1, controls$t1, `-`))
}

(absdist <- mdist(sdiffs, structure.fmla = pr ~ 1, data = nuclearplants))

### Pair matching on the variable `t1`:
pairmatch(absdist)


### Propensity score matching within subgroups:
mdist.examples$ps3 <- mdist(aGlm, structure.fmla=~pt)
fullmatch(mdist.examples$ps3)

### Propensity score matching with a propensity score caliper:
mdist.examples$pscal <- mdist.examples$ps1 + caliper(1,mdist.examples$ps1)
fullmatch(mdist.examples$pscal) # Note that the caliper excludes some units

### A Mahalanobis distance for matching within subgroups:
mdist.examples$mh2 <- mdist(pr ~ t1 + t2 | pt, data = nuclearplants)
all.equal(mdist.examples$mh2,
          mdist(pr ~ t1 + t2, structure.fmla = ~ pt, data = nuclearplants))

### Mahalanobis matching within subgroups, with a propensity score
### caliper:
fullmatch(mdist.examples$mh2 + caliper(1, mdist.examples$ps3))

}

\keyword{nonparametric }% at least one, from doc/KEYWORDS
